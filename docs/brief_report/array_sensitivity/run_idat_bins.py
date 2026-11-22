"""Effect-size and context summaries on the fixed GSE157078 GEO support."""
from pathlib import Path
import argparse
import gzip

import numpy as np
import pandas as pd


WORKFLOWS = ("IDAT_raw", "Noob_background", "ssNoob_dye")
TYPES = ("Type I red", "Type I green", "Type II")
CONTEXTS = ("Island", "Shore", "Shelf", "Open sea")
BINS = ("0.05–<0.10", "0.10–<0.15", "≥0.15")


def load_effects(beta_path, signals, manifest):
    with gzip.open(manifest, "rt", encoding="utf-8-sig", errors="replace") as handle:
        header = next(i for i, line in enumerate(handle) if line.startswith("IlmnID,"))
    cols = ["IlmnID", "Name", "CHR", "MAPINFO", "Color_Channel",
            "Infinium_Design_Type", "Relation_to_UCSC_CpG_Island"]
    ann = pd.read_csv(manifest, skiprows=header, usecols=cols, dtype=str)
    ann = ann.loc[ann.CHR.notna() & ann.MAPINFO.str.fullmatch(r"\d+", na=False)]
    ids = ann.Name.fillna(ann.IlmnID)
    ann = ann.set_index(ids)
    assert ann.index.is_unique
    fields = ["ID_REF", *[f"SAMPLE {i} Detection Pval" for i in range(1, 7)]]
    geo = pd.read_csv(signals, sep="\t", usecols=fields, index_col="ID_REF")
    ids = geo.index[geo.le(.01).all(axis=1)].intersection(ann.index)
    beta = pd.read_csv(beta_path, sep="\t", index_col="probe_id").loc[ids]
    specs = ann.loc[ids]
    result = pd.DataFrame(index=ids)
    result["type"] = np.where(
        specs.Infinium_Design_Type.eq("II"), "Type II",
        specs.Color_Channel.map({"Red": "Type I red", "Grn": "Type I green", "Green": "Type I green"}))
    assert result.type.value_counts().to_dict() == dict(zip(TYPES, (92210, 49935, 723598)))
    ctx = specs.Relation_to_UCSC_CpG_Island.str.replace(r"^[NS]_", "", regex=True)
    result["context"] = ctx.fillna("Open sea").replace({"OpenSea": "Open sea"})
    assert set(result.context) == set(CONTEXTS)
    for mode in WORKFLOWS:
        values = beta[[f"{mode}_{i}" for i in range(1, 7)]].to_numpy()
        assert np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all()
        result[mode] = values[:, 3:].mean(axis=1) - values[:, :3].mean(axis=1)
    return result


def cells(effect):
    gain, loss = int((effect > 0).sum()), int((effect < 0).sum())
    return gain, loss, gain + loss, gain / (gain + loss) if gain + loss else np.nan


def summarize(data):
    bins, context = [], []
    for workflow in WORKFLOWS:
        for probe_type in ("All", *TYPES):
            scope = data if probe_type == "All" else data.loc[data.type.eq(probe_type)]
            e = scope[workflow]
            for label, low, high in zip(BINS, (.05, .10, .15), (.10, .15, np.inf)):
                values = e.loc[e.abs().ge(low) & e.abs().lt(high)]
                bins.append([workflow, probe_type, label, len(e), *cells(values)])
            for ctx in CONTEXTS:
                if workflow == "Noob_background":
                    continue
                values = e.loc[scope.context.eq(ctx)]
                context.append([workflow, probe_type, ctx, len(values), *cells(values[values.abs().ge(.05)])])
    columns = ["workflow", "probe_type", "category", "eligible", "gains", "losses", "selected", "gain_share"]
    return pd.DataFrame(bins, columns=columns), pd.DataFrame(context, columns=columns)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("idat_beta", "geo_signals", "epic_manifest", "out"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    data = load_effects(args.idat_beta, args.geo_signals, args.epic_manifest)
    bins, context = summarize(data)
    expected = [(13822, 1780), (18639, 2096), (4245, 4497)]
    for frame, modes in [(bins, WORKFLOWS), (context, ("IDAT_raw", "ssNoob_dye"))]:
        for mode in modes:
            rows = frame.loc[frame.workflow.eq(mode)]
            all_rows = rows.loc[rows.probe_type.eq("All"), ["gains", "losses"]].sum()
            subtypes = rows.loc[~rows.probe_type.eq("All"), ["gains", "losses"]].sum()
            assert tuple(all_rows) == tuple(subtypes) == expected[WORKFLOWS.index(mode)]
    args.out.mkdir(parents=True, exist_ok=True)
    bins.to_csv(args.out / "effect_size_by_probe_type.tsv", sep="\t", index=False, float_format="%.9g")
    context.to_csv(args.out / "gain_share_by_cpg_context.tsv", sep="\t", index=False, float_format="%.9g")
    selected = data.loc[data.ssNoob_dye.abs().ge(.05), "ssNoob_dye"]
    selected.to_csv(args.out / "ssnoob_signature.tsv", sep="\t", index_label="probe_id", header=["delta_beta"], float_format="%.15g")


if __name__ == "__main__":
    main()
