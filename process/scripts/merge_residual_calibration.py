#!/usr/bin/env python3
"""Add exact [CHANNEL] residual offsets to a base calibration file."""

import argparse
from pathlib import Path


def channel_rows(path: Path):
    rows = {}
    section = None
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        text = line.split("#", 1)[0].strip()
        if not text:
            continue
        if text.startswith("["):
            section = text
            continue
        if section != "[CHANNEL]":
            continue
        fields = text.split()
        if len(fields) != 5 or "*" in fields[:4]:
            continue
        try:
            key = tuple(int(value) for value in fields[:4])
            value = float(fields[4])
        except ValueError as error:
            raise ValueError(f"{path}:{lineno}: malformed CHANNEL row") from error
        if key in rows:
            raise ValueError(f"{path}:{lineno}: duplicate CHANNEL row for {key}")
        rows[key] = value
    return rows


def merge(base_path: Path, residual_path: Path, output_path: Path):
    residual = channel_rows(residual_path)
    lines = base_path.read_text().splitlines()
    updated = set()
    section = None

    for index, line in enumerate(lines):
        text = line.split("#", 1)[0].strip()
        if text.startswith("["):
            section = text
            continue
        if section != "[CHANNEL]" or not text:
            continue
        fields = text.split()
        if len(fields) != 5 or "*" in fields[:4]:
            continue
        key = tuple(int(value) for value in fields[:4])
        if key in residual:
            lines[index] = " ".join(fields[:4] + [f"{float(fields[4]) + residual[key]:.9f}"])
            updated.add(key)

    missing = sorted(set(residual) - updated)
    lines.extend([
        "",
        "# Residual-only exact channel offsets added by merge_residual_calibration.py.",
        "[CHANNEL]",
        "# device fifo column pixel offset",
    ])
    lines.extend(
        f"{device} {fifo} {column} {pixel} {residual[(device, fifo, column, pixel)]:.9f}"
        for device, fifo, column, pixel in missing
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "# Calibration derived by adding residual CHANNEL offsets.\n"
        f"# Base: {base_path}\n"
        f"# Residuals: {residual_path}\n\n"
        + "\n".join(lines)
        + "\n"
    )
    print(f"residual rows: {len(residual)}")
    print(f"base rows updated: {len(updated)}")
    print(f"residual-only rows added: {len(missing)}")
    print(f"output: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path, help="base calibration.conf")
    parser.add_argument("--residual", required=True, type=Path, help="fit_calib residual CHANNEL file")
    parser.add_argument("--output", required=True, type=Path, help="derived calibration.conf")
    args = parser.parse_args()
    try:
        merge(args.base, args.residual, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
