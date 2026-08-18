#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
import types
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    model_dir = Path(args.model_dir).resolve()
    sys.path.insert(0, str(model_dir))

    # The estimator does not use pandas, but scikit-learn may probe for it
    # while validating numpy inputs. Keep validation usable on systems where
    # pandas is absent or binary-incompatible with the sklearn extraction env.
    pandas = types.ModuleType("pandas")
    pandas.DataFrame = type("DataFrame", (), {})
    pandas.Series = type("Series", (), {})
    polars = types.ModuleType("polars")
    polars.DataFrame = type("DataFrame", (), {})
    polars.Series = type("Series", (), {})
    sys.modules.setdefault("pandas", pandas)
    sys.modules.setdefault("polars", polars)

    from estimator import TimingEventEstimator

    events = np.loadtxt(args.input, delimiter=",", dtype=np.float64)
    if events.ndim == 1:
      events = events[None, :]

    estimator = TimingEventEstimator(model_dir)
    out = estimator.predict_batch(events)

    keys = [
        "T0_native",
        "sigma0_native",
        "T1_native",
        "sigma1_native",
        "T_event_native",
        "sigma_event_native",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(keys + ["DeltaT_native"])
        for i in range(events.shape[0]):
            row = [float(out[key][i]) for key in keys]
            row.append(float(out["T1_native"][i] - out["T0_native"][i]))
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
