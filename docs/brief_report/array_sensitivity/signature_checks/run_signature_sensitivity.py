"""Reconstruct the frozen signature and run the dated post hoc checks."""
import argparse
import csv
import gzip
import importlib.util
from itertools import combinations
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd


COHORTS = ("GSE74486", "GSE73747")
ORIENTATION_SAMPLE = 3
OBSERVED = (3, 4, 5)  # Zero-based indices. GEO SAMPLE 4–6 are the overexpression arrays.


def read_module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def annotation_probes(path):
    with gzip.open(path, "rt", newline="", errors="replace") as stream:
        for line in stream:
            if line.startswith("IlmnID,"):
                header = next(csv.reader([line]))
                break
        else:
            raise ValueError("EPIC manifest header is missing")
        return {(row["Name"] or row["IlmnID"]) for row in csv.DictReader(stream, fieldnames=header)
                if (row["Name"] or row["IlmnID"]) and (row["CHR"] or "").strip() and (row["MAPINFO"] or "").strip().isdigit()}


def load_discovery(signals, manifest):
    columns = [field for i in range(1, 7)
               for field in (f"SAMPLE {i} Methylated signal", f"SAMPLE {i} Unmethylated Signal",
                             f"SAMPLE {i} Detection Pval")]
    all_signals = pd.read_csv(signals, sep="\t", index_col="ID_REF", low_memory=False)
    assert all_signals.index.is_unique and set(all_signals.columns) == set(columns)
    stats = all_signals[columns].to_numpy(dtype=float).reshape(-1, 6, 3)
    den = stats[:, :, 0] + stats[:, :, 1] + 100
    good = np.isfinite(stats).all(axis=(1, 2)) & (den > 0).all(axis=1) & (stats[:, :, 2] <= .01).all(axis=1)
    eligible = good & all_signals.index.isin(annotation_probes(manifest))
    beta = stats[eligible, :, 0] / den[eligible, :]
    print(f"{len(stats)} signal rows, {int(good.sum())} pass detection and {len(beta)} pass discovery QC", flush=True)
    return pd.DataFrame(beta, index=all_signals.index[eligible],
                        columns=[f"sample_{i}" for i in range(1, 7)])


def discovery_samples(soft, geo_inventory):
    chip_lookup = {gsm: (chip, pos) for gsm, chip, pos in re.findall(
        r"(GSM\d+)_(\d+)_(R\d+C\d+)_Grn\.idat\.gz", geo_inventory.read_text())}
    records = []
    for record in soft.read_text().split("^SAMPLE = ")[1:]:
        gsm = record.splitlines()[0].strip()
        sample = int(re.search(r"^!Sample_description = SAMPLE (\d+)$", record, re.M)[1])
        group = re.search(r"^!Sample_characteristics_ch1 = genotype: (.+)$", record, re.M)[1]
        chip, position = chip_lookup[gsm]
        records.append((sample, gsm, group, chip, position))
    table = pd.DataFrame(records, columns=["sample", "geo_accession", "group", "chip_id", "position"]).sort_values("sample")
    assert table["sample"].tolist() == list(range(1, 7))
    assert tuple(np.flatnonzero(table.group.eq("DNMT3L overexpression"))) == OBSERVED
    assert table.group.value_counts().to_dict() == {"CON": 3, "DNMT3L overexpression": 3}
    return table


def observed_signatures(beta, frozen, out):
    values = beta.to_numpy()
    rows = []
    signatures = {}
    for subset in combinations(range(6), 3):
        if ORIENTATION_SAMPLE not in subset:
            continue
        other = tuple(i for i in range(6) if i not in subset)
        delta = values[:, subset].mean(axis=1) - values[:, other].mean(axis=1)
        selected = np.abs(delta) >= .05
        signature = pd.Series(delta[selected], index=beta.index[selected], name="delta_beta")
        overlap = signature.index.intersection(frozen.index)
        rows.append({
            "observed": subset == OBSERVED,
            "putative_overexpression_samples": ",".join(str(i + 1) for i in subset),
            "putative_control_samples": ",".join(str(i + 1) for i in other),
            "original_overexpression_in_putative_group": len(set(subset).intersection(OBSERVED)),
            "threshold_count": len(signature),
            "gain_count": int((signature > 0).sum()),
            "loss_count": int((signature < 0).sum()),
            "overlap_with_observed": len(overlap),
            "same_sign_in_overlap": int((np.sign(signature.loc[overlap]) == np.sign(frozen.loc[overlap])).sum()),
        })
        if subset == OBSERVED:
            assert len(signature) == len(frozen) == 15602
            np.testing.assert_allclose(signature.reindex(frozen.index), frozen, atol=5e-7, rtol=5e-6)
            np.testing.assert_array_equal(np.sign(signature.reindex(frozen.index)), np.sign(frozen))
            case = beta.loc[frozen.index].to_numpy()[:, subset]
            control = beta.loc[frozen.index].to_numpy()[:, other]
            separated = np.where(frozen.to_numpy() > 0, case.min(axis=1) > control.max(axis=1),
                                 case.max(axis=1) < control.min(axis=1))
            signatures["effect_size"] = frozen
            signatures["complete_separation"] = frozen.loc[separated]
    assert len(rows) == 10
    pd.DataFrame(rows).sort_values("observed", ascending=False, kind="stable").to_csv(
        out / "label_split_counts.tsv", sep="\t", index=False)
    with tempfile.TemporaryDirectory(prefix="dnmt3l-limma-") as temp:
        matrix = Path(temp) / "beta.tsv"
        statistics = Path(temp) / "moderated.tsv"
        beta.to_csv(matrix, sep="\t", index_label="probe_id", float_format="%.12g")
        subprocess.run(["Rscript", str(Path(__file__).with_name("limma_signature.R")), str(matrix), str(statistics)], check=True)
        tests = pd.read_csv(statistics, sep="\t", index_col="probe_id")
        assert tests.index.equals(beta.index) and tests[["p", "q"]].notna().all().all()
    tests = tests.loc[frozen.index]
    for name, key in (("limma_nominal", "p"), ("limma_bh", "q")):
        signatures[name] = frozen.loc[tests[key] < .05]
    selected = pd.DataFrame({"delta_beta": frozen, "complete_separation": frozen.index.isin(signatures["complete_separation"].index)})
    selected.join(tests).to_csv(out / "discovery_filter_statistics.tsv.gz", sep="\t", index=True, index_label="probe_id",
                                 float_format="%.12g", compression={"method": "gzip", "mtime": 0})
    return signatures


def donor_checks(compendium, cache, signatures):
    sensitivity = read_module(compendium / "sensitivity/run_sensitivity.py")
    original = sensitivity.original_code(compendium / "analyze.py")
    annotation = pd.read_pickle(cache / "eligible_annotations.pkl")
    archived = pd.read_csv(compendium / "estimates.tsv", sep="\t").query('statistic == "DNMT3L_OR"').set_index("context")
    rows = []
    for cohort in COHORTS:
        meta, beta = original.load_cohort(cache, cohort, None, False)
        boot, _, _ = original.fit_weights(meta, np.random.default_rng(original.SEED), original.N_BOOT)
        boot = boot[:, :1 + original.N_BOOT]
        design, _ = original.designs(meta)
        permutations = sensitivity.residual_weights(design, np.delete(design, 1, axis=1), sensitivity.donor_blocks(meta))
        probes = sensitivity.probe_frame(original, beta, annotation, signatures["effect_size"])
        baseline = None
        for name, signature in signatures.items():
            probes["signature"] = signature.reindex(probes.index)
            n_up = int(probes.signature.gt(0).sum())
            n_down = int(probes.signature.lt(0).sum())
            record = {
                "cohort": cohort, "signature": name, "discovery_n": len(signature), "callable_n": n_up + n_down,
                "gain_n": n_up, "loss_n": n_down, "valid_bootstrap": 0, "valid_permutation": 0,
                "undefined_permutation": 0,
                **dict.fromkeys(("or", "ci_low", "ci_high", "ratio_to_original", "ratio_ci_low", "ratio_ci_high",
                                 "permutation_p_two_sided", "null_low", "null_high"), np.nan),
            }
            if n_up and n_down:
                signed, panel, n_strata = sensitivity.signature_panel(probes)
                y = np.ascontiguousarray(beta.loc[signed.index], dtype="float64")
                scores = sensitivity.or_from_weights(original, y, boot, panel, n_strata)
                if name == "effect_size":
                    baseline = scores
                    reference = archived.loc[cohort]
                    np.testing.assert_allclose([scores[0], *np.quantile(scores[1:], [.025, .975])],
                                               reference[["estimate", "ci_low", "ci_high"]].astype(float), rtol=2e-5)
                if np.isfinite(scores[0]) and scores[0] > 0:
                    valid_boot = scores[1:][np.isfinite(scores[1:]) & (scores[1:] >= 0)]
                    ratio = scores / baseline
                    valid_ratio = ratio[1:][np.isfinite(ratio[1:]) & (ratio[1:] >= 0)]
                    perm = sensitivity.or_from_weights(original, y, permutations, panel, n_strata)
                    valid_perm = perm[np.isfinite(perm) & (perm >= 0)]
                    undefined = len(perm) - len(valid_perm)
                    with np.errstate(divide="ignore"):
                        tail = np.count_nonzero(np.abs(np.log(valid_perm)) >= abs(np.log(scores[0])) - 1e-10)
                    record["or"] = scores[0]
                    record["ci_low"], record["ci_high"] = np.quantile(valid_boot, [.025, .975])
                    record["ratio_to_original"] = ratio[0]
                    record["ratio_ci_low"], record["ratio_ci_high"] = np.quantile(valid_ratio, [.025, .975])
                    record["valid_bootstrap"], record["valid_permutation"] = len(valid_boot), len(valid_perm)
                    record["undefined_permutation"] = undefined
                    record["permutation_p_two_sided"] = (1 + tail + undefined) / (1 + len(perm))
                    record["null_low"], record["null_high"] = np.quantile(valid_perm, [.025, .975])
            rows.append(record)
        print(f"Completed donor checks for {cohort}.", flush=True)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("signals", "manifest", "soft", "geo_inventory", "compendium", "cache", "out"):
        parser.add_argument(arg, type=Path)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    discovery_samples(a.soft, a.geo_inventory).to_csv(a.out / "array_positions.tsv", sep="\t", index=False)
    frozen = pd.read_csv(a.compendium / "GSE157078_candidate_dmcs.tsv", sep="\t", index_col="probe_id").delta_beta
    beta = load_discovery(a.signals, a.manifest)
    signatures = observed_signatures(beta, frozen, a.out)
    counts = [{"signature": name, "discovery_n": len(s), "gain_n": int(s.gt(0).sum()), "loss_n": int(s.lt(0).sum())}
              for name, s in signatures.items()]
    pd.DataFrame(counts).to_csv(a.out / "signature_counts.tsv", sep="\t", index=False)
    donor_checks(a.compendium, a.cache, signatures).to_csv(a.out / "donor_signature_sensitivity.tsv", sep="\t",
                                                           index=False, float_format="%.9g", na_rep="NA")


if __name__ == "__main__":
    main()
