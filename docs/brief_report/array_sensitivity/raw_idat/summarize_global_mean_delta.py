#!/usr/bin/env python3
"""Summarize the global GSE157078 beta difference after preprocessing.

Usage:
    python summarize_global_mean_delta.py IDAT_BETA.tsv.gz GSE157078_signals.txt.gz OUTPUT.tsv

``IDAT_BETA.tsv.gz`` is the local intermediate written by ``process_idats.R``.
The result uses the original GEO quality universe: detection P <= 0.01 in all
six arrays. For every CpG, the difference is mean(D1, D2, D3) minus mean(Z1,
Z2, Z3). The global mean is the arithmetic mean of these CpG-level values.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOWS = {
    "Raw beta": "IDAT_raw",
    "Noob background only": "Noob_background",
    "ssNoob with dye correction": "ssNoob_dye",
}


def main(idat_beta: Path, geo_signals: Path, output: Path) -> None:
    p_columns = ["ID_REF", *[f"SAMPLE {i} Detection Pval" for i in range(1, 7)]]
    detection = pd.read_csv(geo_signals, sep="\t", compression="gzip", usecols=p_columns)
    passing = set(detection.loc[
        detection[[f"SAMPLE {i} Detection Pval" for i in range(1, 7)]].le(0.01).all(axis=1),
        "ID_REF",
    ])

    values = {
        label: {"n": 0, "control": 0.0, "overexpression": 0.0,
                "delta": 0.0, "absolute_delta": 0.0, "all_delta": []}
        for label in WORKFLOWS
    }
    for chunk in pd.read_csv(idat_beta, sep="\t", compression="gzip", chunksize=100_000):
        chunk = chunk.loc[chunk["probe_id"].isin(passing)]
        for label, prefix in WORKFLOWS.items():
            control = chunk[[f"{prefix}_{i}" for i in range(1, 4)]].to_numpy(float)
            overexpression = chunk[[f"{prefix}_{i}" for i in range(4, 7)]].to_numpy(float)
            delta = overexpression.mean(axis=1) - control.mean(axis=1)
            keep = np.isfinite(delta)
            delta = delta[keep]
            record = values[label]
            record["n"] += len(delta)
            record["control"] += np.nansum(control[keep])
            record["overexpression"] += np.nansum(overexpression[keep])
            record["delta"] += delta.sum()
            record["absolute_delta"] += np.abs(delta).sum()
            record["all_delta"].append(delta)

    rows = []
    for label, record in values.items():
        n = record["n"]
        delta = np.concatenate(record["all_delta"])
        rows.append({
            "workflow": label,
            "CpGs": n,
            "control_global_mean_beta": record["control"] / (3 * n),
            "overexpression_global_mean_beta": record["overexpression"] / (3 * n),
            "global_mean_delta_beta_OE_minus_control": record["delta"] / n,
            "median_CpG_delta_beta": np.median(delta),
            "mean_absolute_CpG_delta_beta": record["absolute_delta"] / n,
        })
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False, float_format="%.10f")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
