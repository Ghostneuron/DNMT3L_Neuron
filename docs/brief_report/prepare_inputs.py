"""Extract fetal metadata and quality-filtered CpG matrices from public GEO files."""
import argparse
import gzip
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook


def series_header(path):
    cols = {}
    with gzip.open(path, "rt") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith("!series_matrix_table_begin"):
                return cols, line_no
            if line.startswith("!Sample_"):
                key, *values = line.rstrip("\n").replace('"', "").split("\t")
                if key == "!Sample_characteristics_ch1":
                    key = values[0].split(": ", 1)[0]
                    values = [value.split(": ", 1)[-1] for value in values]
                cols[key.removeprefix("!Sample_")] = values
    raise ValueError(f"No matrix in {path}")


def fetal_ages(path):
    workbook = load_workbook(path / "Mendioroz_2015_Tables.xlsx", read_only=True, data_only=True)
    sheet = workbook["Table S1 Cases and Controls"]
    ages = {}
    for sample, _, _, donor, _, tissue, diagnosis, age, sex in sheet.iter_rows(min_row=2, max_col=9, values_only=True):
        if tissue != "Fetal Cerebrum":
            continue
        weeks = re.fullmatch(r"\s*(\d+(?:\.\d+)?)(?:-(\d+(?:\.\d+)?))?\s*w\s*", age)
        assert weeks is not None
        endpoints = [float(value) for value in weeks.groups() if value is not None]
        ages[f"tb{sample}"] = (donor, int(diagnosis == "DS"), sex, float(np.mean(endpoints)), age.strip())
    workbook.close()
    assert len(ages) == 14
    return ages


def load_74486(path):
    wide = path / "GSE74486_precombat.txt.gz"
    metadata = []
    for record in (path / "GSE74486_brief.soft").read_text().split("^SAMPLE = ")[1:]:
        if "developmental stage: Fetus" not in record or "tissue: Cerebrum" not in record:
            continue
        metadata.append({
            "array": record.splitlines()[0].strip(),
            "donor": re.search(r"of patient ([^\n]+)", record)[1].strip(),
            "column": re.search(r"^!Sample_description = ((?:ds|tb)\d+)\s*$", record, re.M)[1],
            "region": "cerebrum",
            "sex": re.search(r"gender: (Male|Female)", record)[1].lower(),
            "dx": int("disease state: Down syndrome" in record),
        })
    meta = pd.DataFrame(metadata).set_index("array")
    assert meta.donor.nunique() == 14 and meta.dx.sum() == 8
    columns = meta.column.tolist()
    ages = fetal_ages(path)
    for array, sample in meta.iterrows():
        donor, dx, sex, age, source_age = ages[sample.column]
        assert (sample.donor, sample.dx, sample.sex[0].upper()) == (donor, dx, sex)
        meta.loc[array, "age"] = age
        meta.loc[array, "reported_age"] = source_age
    quality = [col + "_DetectionPval" for col in columns]
    dat = pd.read_csv(wide, sep="\t", usecols=["ID_REF", *columns, *quality], index_col="ID_REF", low_memory=False)
    beta = dat[columns].to_numpy(dtype="float32")
    beta[~(dat[quality].to_numpy(dtype=float) <= 0.01)] = np.nan
    return meta, pd.DataFrame(beta, index=dat.index, columns=meta.index)


def load_signal(path, order, methylated_first):
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    indices = [0, *(index for sample in order for index in (3 * sample + 1, 3 * sample + 2, 3 * sample + 3))]
    for sample in order:
        first, second, third = [header[3 * sample + k].lower() for k in (1, 2, 3)]
        assert ("unmethylated" not in first) == methylated_first
        assert ("unmethylated" in second) == methylated_first
        assert "detection pval" in third
    table = pd.read_csv(path, sep="\t", usecols=indices, index_col="ID_REF", low_memory=False)
    values = table.to_numpy(dtype=np.float32).reshape(len(table), len(order), 3)
    methylated = values[:, :, 0 if methylated_first else 1]
    unmethylated = values[:, :, 1 if methylated_first else 0]
    beta = methylated / (methylated + unmethylated + 100)
    beta[~(values[:, :, 2] <= 0.01)] = np.nan
    return pd.DataFrame(beta, index=table.index)


def infer_duplicates(beta, meta):
    fingerprints = beta.loc[beta.index.str.startswith("rs")].to_numpy()
    assert fingerprints.shape[0] == 65
    # Classify the three signal-ratio clusters, omitting their uncertain gaps.
    genotypes = np.full_like(fingerprints, np.nan)
    genotypes[fingerprints <= 0.20] = 0
    genotypes[(fingerprints >= 0.35) & (fingerprints <= 0.65)] = 1
    genotypes[fingerprints >= 0.80] = 2
    donors = list(meta.index)
    diagnostics = []
    for i, j in combinations(range(len(meta)), 2):
        xi, xj = genotypes[:, i], genotypes[:, j]
        joint = np.isfinite(xi) & np.isfinite(xj)
        overlap = int(joint.sum())
        rate = float(np.mean(xi[joint] == xj[joint])) if overlap else np.nan
        diagnostics.append((meta.index[i], meta.index[j], overlap, rate))
        if overlap < 45 or rate < 0.95:
            continue
        assert meta.iloc[i].dx == meta.iloc[j].dx
        assert meta.iloc[i].sex == meta.iloc[j].sex
        assert meta.iloc[i].age == meta.iloc[j].age
        previous = donors[j]
        donors = [donors[i] if donor == previous else donor for donor in donors]
    comparisons = pd.DataFrame(diagnostics, columns=["array_1", "array_2", "joint_snps", "genotype_agreement"])
    return donors, comparisons.sort_values(["genotype_agreement", "joint_snps"], ascending=False)


def load_73747(path, cache):
    header, _ = series_header(path / "GSE73747_series_matrix.txt.gz")
    meta = pd.DataFrame(header).set_index("geo_accession")
    meta = meta.loc[meta.tissue.str.contains("fetal"), ["title", "tissue", "disease state", "gender", "gestational.week"]].copy()
    meta.columns = ["title", "region", "dx", "sex", "age"]
    meta.region = meta.region.str.split().str[0]
    meta.dx = meta.dx.eq("trisomy").astype(int)
    meta.age = meta.age.astype(float)
    order = [int(re.match(r"Sample (\d+)", title)[1]) - 1 for title in meta.title]
    assert len(meta) == 59 and meta.dx.sum() == 24 and order == sorted(order)
    beta = load_signal(path / "GSE73747_signal.txt.gz", order, methylated_first=True)
    beta.columns = meta.index
    meta["donor"], comparisons = infer_duplicates(beta, meta)
    comparisons.to_csv(cache / "GSE73747_snp_duplicate_diagnostics.tsv", sep="\t", index=False)
    return meta, beta


def prepare_73747_series(inputs, cache):
    matrix = inputs / "GSE73747_series_matrix.txt.gz"
    _, first_row = series_header(matrix)
    quality = pd.read_pickle(cache / "GSE73747_beta.pkl")
    deposited = pd.read_csv(matrix, sep="\t", skiprows=first_row, usecols=["ID_REF", *quality.columns],
                            index_col="ID_REF", comment="!", low_memory=False)
    deposited = deposited.reindex(index=quality.index, columns=quality.columns)
    deposited = deposited.where(quality.notna()).astype("float32")
    deposited.to_pickle(cache / "GSE73747_series_beta.pkl")


def annotation(path):
    file = path / "GPL13534_annotation.soft.gz"
    with gzip.open(file, "rt") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith("!platform_table_begin"):
                break
        else:
            raise ValueError("No GPL13534 platform annotation")
        # Avoid reading the large SOFT sample records after the annotation table.
        selected = []
        for line in handle:
            if line.startswith("!platform_table_end"):
                break
            fields = line.rstrip("\n").split("\t")
            if not selected:
                header = {name: i for i, name in enumerate(fields)}
                keys = [header[name] for name in ("ID", "Relation_to_UCSC_CpG_Island")]
            selected.append([fields[i] for i in keys])
    gpl = pd.DataFrame(selected[1:], columns=["probe", "island"]).set_index("probe")
    zhou = pd.read_csv(path / "HM450_Zhou_manifest.tsv.gz", sep="\t",
                       usecols=["probeID", "CpG_chrm", "CpG_beg", "designType", "MASK_general"], index_col="probeID")
    keep = zhou.index.str.startswith("cg") & zhou.CpG_chrm.isin([f"chr{i}" for i in range(1, 23)]) & zhou.MASK_general.eq(False)
    zhou = zhou.loc[keep].drop(columns="MASK_general")
    assert zhou.index.is_unique and len(zhou) > 300_000
    gpl.island = gpl.island.str.replace(r"^[NS]_", "", regex=True).replace("", "OpenSea")
    return zhou.join(gpl, how="inner").sort_index()


def load_comparator(path):
    matrix = path / "GSE120558_EPIC_matrix.txt.gz"
    header, first_row = series_header(matrix)
    metadata = pd.DataFrame({"title": header["title"], "genotype": header["genotype"]}, index=header["geo_accession"])
    case_a = metadata.index[metadata.title.str.startswith("patient_1_fibroblast")]
    case_b = metadata.index[metadata.title.eq("patient_2_fibroblast")]
    controls = metadata.index[metadata.title.str.fullmatch(r"control_\d_fibroblast")]
    assert len(case_a) == 2 and len(case_b) == 1 and len(controls) == 3
    assert metadata.loc[case_a.union(case_b), "genotype"].str.contains("W330R").all()
    columns = [*case_a, *case_b, *controls]
    beta = pd.read_csv(matrix, sep="\t", skiprows=first_row, usecols=["ID_REF", *columns], index_col="ID_REF",
                       comment="!", low_memory=False)[columns].dropna()
    reference = beta[controls].mean(axis=1)
    a = beta[case_a].mean(axis=1) - reference
    b = beta[case_b].mean(axis=1) - reference
    effects = pd.DataFrame({"delta_a": a, "delta_b": b, "delta": (a + b) / 2})
    return metadata.loc[columns], effects


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path)
    parser.add_argument("cache", type=Path)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--series-only", action="store_true")
    modes.add_argument("--cerebrum-only", action="store_true")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.series_only:
        prepare_73747_series(args.inputs, args.cache)
        return
    if args.cerebrum_only:
        meta, _ = load_74486(args.inputs)
        meta.to_csv(args.cache / "GSE74486_sample_manifest.tsv", sep="\t")
        return
    ann = annotation(args.inputs)
    oe = load_signal(args.inputs / "GSE157078_signals.txt.gz", range(6), methylated_first=False)
    ann = ann.loc[ann.index.intersection(oe.dropna().index)]
    ann.to_pickle(args.cache / "eligible_annotations.pkl")
    for accession, (meta, beta) in [
        ("GSE74486", load_74486(args.inputs)), ("GSE73747", load_73747(args.inputs, args.cache)),
    ]:
        meta.to_csv(args.cache / f"{accession}_sample_manifest.tsv", sep="\t")
        beta.loc[ann.index.intersection(beta.index)].to_pickle(args.cache / f"{accession}_beta.pkl")
        print(accession, len(meta), meta.donor.nunique(), meta.groupby(["region", "dx"]).donor.nunique().to_dict(), flush=True)
    meta, effects = load_comparator(args.inputs)
    effects.loc[ann.index.intersection(effects.index)].to_pickle(args.cache / "dnmt3a_effects.pkl")
    meta.to_csv(args.cache / "dnmt3a_manifest.tsv", sep="\t")


if __name__ == "__main__":
    main()
