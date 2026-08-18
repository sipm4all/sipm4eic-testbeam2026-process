#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--events", type=int, default=1000)
    args = parser.parse_args()

    rng = np.random.default_rng(0x20260818)
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        for _ in range(args.events):
            base = rng.uniform(-200.0, 200.0)
            t1 = base + 0.1 * rng.normal(0.0, 0.35)
            event = np.empty(64, dtype=np.float64)
            for i in range(32):
                event[i] = base + 0.03 * i + rng.normal(0.0, 0.35) + 0.2 * rng.normal(0.0, 0.8)
            for i in range(32):
                event[32 + i] = t1 - 0.02 * i + rng.normal(0.0, 0.35) + 0.2 * rng.normal(0.0, 0.8)
            writer.writerow([f"{x:.17g}" for x in event])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
