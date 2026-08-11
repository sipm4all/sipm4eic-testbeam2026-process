#!/usr/bin/env python3
import argparse
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, order=True)
class Key:
    device: int
    fifo: int
    column: int
    pixel: int
    tdc: int

    @property
    def chip(self) -> int:
        return self.fifo // 4

    @property
    def fifo_in_chip(self) -> int:
        return self.fifo % 4

    @property
    def column_in_fifo(self) -> int:
        return self.column - 2 * self.fifo_in_chip


@dataclass
class Row:
    key: Key
    off: float
    iif: float
    source: str
    method: str = "measured"
    nsource: int = 0
    off_rms: float = 0.0
    iif_rms: float = 0.0


def parse_int_list(expr):
    out = []
    for part in expr.split(','):
        part = part.strip()
        if not part:
            continue
        if '..' in part:
            a, b = part.split('..', 1)
            a = int(a)
            b = int(b)
            if a > b:
                raise ValueError(f"invalid range {part}")
            out.extend(range(a, b + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def median(values):
    return statistics.median(values)


def rms(values):
    if not values:
        return 0.0
    m = sum(values) / len(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / len(values))


def load_tdc(filename):
    rows = {}
    duplicates = []
    current = None
    with open(filename) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('[') and line.endswith(']'):
                current = line[1:-1].strip()
                continue
            if current != 'TDC':
                continue
            fields = line.split()
            if len(fields) != 7:
                raise RuntimeError(f"malformed [TDC] row at {filename}:{lineno}: {raw.rstrip()}")
            try:
                device, fifo, column, pixel, tdc = map(int, fields[:5])
                off = float(fields[5])
                iif = float(fields[6])
            except ValueError as e:
                raise RuntimeError(f"bad [TDC] value at {filename}:{lineno}: {raw.rstrip()}") from e
            key = Key(device, fifo, column, pixel, tdc)
            row = Row(key, off, iif, f"{filename}:{lineno}")
            if key in rows:
                duplicates.append((key, rows[key].source, row.source))
            else:
                rows[key] = row
    if duplicates:
        msg = ["duplicate concrete [TDC] rows found:"]
        for key, first, second in duplicates[:20]:
            msg.append(f"  {key}: {first} and {second}")
        if len(duplicates) > 20:
            msg.append(f"  ... {len(duplicates) - 20} more")
        raise RuntimeError('\n'.join(msg))
    return rows


def valid_columns(fifo):
    first = 2 * (fifo % 4)
    return (first, first + 1)


def expected_keys(devices, fifos):
    keys = []
    for device in devices:
        for fifo in fifos:
            for column in valid_columns(fifo):
                for pixel in range(4):
                    for tdc in range(4):
                        keys.append(Key(device, fifo, column, pixel, tdc))
    return sorted(keys)


def equivalent_position_source(key, rows):
    # Same device and same repeated position across chips.
    return [r for r in rows.values()
            if r.key != key
            and r.key.device == key.device
            and r.key.fifo_in_chip == key.fifo_in_chip
            and r.key.column_in_fifo == key.column_in_fifo
            and r.key.pixel == key.pixel
            and r.key.tdc == key.tdc]


def global_equivalent_position_source(key, rows):
    return [r for r in rows.values()
            if r.key != key
            and r.key.fifo_in_chip == key.fifo_in_chip
            and r.key.column_in_fifo == key.column_in_fifo
            and r.key.pixel == key.pixel
            and r.key.tdc == key.tdc]


def same_device_tdc_source(key, rows):
    return [r for r in rows.values()
            if r.key != key
            and r.key.device == key.device
            and r.key.tdc == key.tdc]


def same_tdc_source(key, rows):
    return [r for r in rows.values()
            if r.key != key
            and r.key.tdc == key.tdc]


def fill_from_source(key, method, source):
    offs = [r.off for r in source]
    iifs = [r.iif for r in source]
    return Row(key=key,
               off=median(offs),
               iif=median(iifs),
               source="filled",
               method=method,
               nsource=len(source),
               off_rms=rms(offs),
               iif_rms=rms(iifs))


def complete(rows, expected):
    completed = dict(rows)
    filled = []
    unfilled = []
    methods = [
        ("same_device_equivalent_chip_position", equivalent_position_source),
        ("global_equivalent_chip_position", global_equivalent_position_source),
        ("same_device_tdc", same_device_tdc_source),
        ("global_tdc", same_tdc_source),
    ]
    for key in expected:
        if key in completed:
            continue
        row = None
        for method, selector in methods:
            source = selector(key, rows)
            if source:
                row = fill_from_source(key, method, source)
                break
        if row is None:
            unfilled.append(key)
        else:
            completed[key] = row
            filled.append(row)
    return completed, filled, unfilled


def global_tdc_rows(rows, expected):
    expected_set = set(expected)
    out = {}
    for tdc in range(4):
        selected = [r for k, r in rows.items() if k in expected_set and k.tdc == tdc]
        if not selected:
            continue
        offs = [r.off for r in selected]
        iifs = [r.iif for r in selected]
        out[tdc] = Row(key=Key(-1, -1, -1, -1, tdc),
                       off=median(offs),
                       iif=median(iifs),
                       source="global_fallback",
                       method="global_tdc_median",
                       nsource=len(selected),
                       off_rms=rms(offs),
                       iif_rms=rms(iifs))
    return out


def single_global_zero_off_row(rows, expected):
    expected_set = set(expected)
    selected = [r for k, r in rows.items() if k in expected_set]
    if not selected:
        return None
    iifs = [r.iif for r in selected]
    return Row(key=Key(-1, -1, -1, -1, -1),
               off=0.0,
               iif=median(iifs),
               source="global_fallback",
               method="single_global_zero_off",
               nsource=len(selected),
               off_rms=0.0,
               iif_rms=rms(iifs))


def phase_delta_summary(rows, expected, fallback_rows, single_global=None):
    expected_set = set(expected)
    measured = [r for k, r in rows.items() if k in expected_set]
    out = []

    def add(name, predictor):
        deltas = []
        for row in measured:
            for fine in (60, 70, 80, 90, 100, 110, 120):
                pred_off, pred_iif = predictor(row, fine)
                pred_phase = pred_off + pred_iif * fine
                true_phase = row.off + row.iif * fine
                deltas.append(pred_phase - true_phase)
        mean = sum(deltas) / len(deltas) if deltas else 0.0
        out.append({
            "name": name,
            "mean": mean,
            "rms": rms(deltas),
            "std": rms([x - mean for x in deltas]),
            "rms_ps": rms(deltas) * 3125.0,
        })

    if fallback_rows:
        add("per_tdc_global", lambda row, fine: (fallback_rows[row.key.tdc].off, fallback_rows[row.key.tdc].iif))
    if single_global is not None:
        add("single_global_off0", lambda row, fine: (0.0, single_global.iif))
    return out


def write_output(filename, rows, expected, single_global):
    expected_set = set(expected)
    with open(filename, 'w') as f:
        f.write("# TDC calibration file with single global fallback row\n")
        f.write("# generated by complete_tdc_calibration.py\n")
        f.write("# measured concrete rows are written first; wildcard fallback row is written last\n")
        f.write("[TDC]\n")
        f.write("# measured rows: device fifo column pixel tdc off iif\n")
        for key in sorted(k for k in rows if k in expected_set):
            row = rows[key]
            f.write(f"{key.device} {key.fifo} {key.column} {key.pixel} {key.tdc} {row.off:.10g} {row.iif:.10g}\n")

        f.write("\n[TDC]\n")
        f.write("# single global fallback row for missing concrete TDC calibrations\n")
        f.write("# device fifo column pixel tdc off iif\n")
        f.write("# '*' address fields have lower specificity than measured concrete rows\n")
        f.write("# single_global_off0 n={n} iif_rms={iif_rms:.8g}\n".format(
            n=single_global.nsource, iif_rms=single_global.iif_rms))
        f.write(f"* * * * * 0 {single_global.iif:.10g}\n")


def coverage_table(rows, expected):
    expected_by_device = Counter(k.device for k in expected)
    present_by_device = Counter(k.device for k in rows if k in set(expected))
    expected_by_device_chip = Counter((k.device, k.chip) for k in expected)
    present_by_device_chip = Counter((k.device, k.chip) for k in rows if k in set(expected))
    return expected_by_device, present_by_device, expected_by_device_chip, present_by_device_chip


def pattern_value(key, name):
    if name == "global_tdc":
        return (key.tdc,)
    if name == "device_tdc":
        return (key.device, key.tdc)
    if name == "fifo_mod_column_side_tdc":
        return (key.fifo_in_chip, key.column_in_fifo, key.tdc)
    if name == "fifo_mod_column_side_pixel_tdc":
        return (key.fifo_in_chip, key.column_in_fifo, key.pixel, key.tdc)
    if name == "chip_repeated_position_tdc":
        return (key.chip, key.fifo_in_chip, key.column_in_fifo, key.tdc)
    if name == "chip_repeated_position_pixel_tdc":
        return (key.chip, key.fifo_in_chip, key.column_in_fifo, key.pixel, key.tdc)
    if name == "device_repeated_position_tdc":
        return (key.device, key.fifo_in_chip, key.column_in_fifo, key.tdc)
    if name == "device_repeated_position_pixel_tdc":
        return (key.device, key.fifo_in_chip, key.column_in_fifo, key.pixel, key.tdc)
    if name == "absolute_fifo_column_tdc":
        return (key.fifo, key.column, key.tdc)
    if name == "absolute_fifo_column_pixel_tdc":
        return (key.fifo, key.column, key.pixel, key.tdc)
    raise RuntimeError(f"unknown pattern {name}")


def percentile_abs(values, fraction):
    if not values:
        return 0.0
    ordered = sorted(abs(v) for v in values)
    index = int(fraction * (len(ordered) - 1))
    return ordered[index]


def pattern_comparison(rows, expected):
    expected_set = set(expected)
    measured = [r for k, r in rows.items() if k in expected_set]
    patterns = [
        "global_tdc",
        "device_tdc",
        "fifo_mod_column_side_tdc",
        "fifo_mod_column_side_pixel_tdc",
        "chip_repeated_position_tdc",
        "chip_repeated_position_pixel_tdc",
        "device_repeated_position_tdc",
        "device_repeated_position_pixel_tdc",
        "absolute_fifo_column_tdc",
        "absolute_fifo_column_pixel_tdc",
    ]

    out = []
    for pattern in patterns:
        groups = defaultdict(list)
        for row in measured:
            groups[pattern_value(row.key, pattern)].append(row)

        off_res = []
        iif_res = []
        source_counts = []
        skipped = 0
        for row in measured:
            group = groups[pattern_value(row.key, pattern)]
            if len(group) <= 1:
                skipped += 1
                continue
            source_counts.append(len(group) - 1)
            offs = [r.off for r in group if r.key != row.key]
            iifs = [r.iif for r in group if r.key != row.key]
            off_res.append(row.off - median(offs))
            iif_res.append(row.iif - median(iifs))

        off_center = median(off_res) if off_res else 0.0
        iif_center = median(iif_res) if iif_res else 0.0
        out.append({
            "pattern": pattern,
            "groups": len(groups),
            "used": len(off_res),
            "skipped": skipped,
            "median_source_rows": median(source_counts) if source_counts else 0,
            "off_rms": rms(off_res),
            "off_mad": median([abs(v - off_center) for v in off_res]) if off_res else 0,
            "off_p90": percentile_abs(off_res, 0.90),
            "iif_rms": rms(iif_res),
            "iif_mad": median([abs(v - iif_center) for v in iif_res]) if iif_res else 0,
            "iif_p90": percentile_abs(iif_res, 0.90),
        })

    out.sort(key=lambda x: (x["iif_rms"], x["off_rms"]))
    return out


def write_report(filename, input_file, output_file, rows, expected, fallback, single_global):
    expected_set = set(expected)
    measured_expected = [k for k in rows if k in expected_set]
    missing = sorted(k for k in expected if k not in rows)
    extra = sorted(k for k in rows if k not in expected_set)
    exp_dev, pre_dev, exp_chip, pre_chip = coverage_table(rows, expected)

    with open(filename, 'w') as f:
        f.write("# TDC Calibration Fallback Report\n\n")
        f.write(f"input: `{input_file}`\n\n")
        f.write(f"output: `{output_file}`\n\n")
        f.write("## Summary\n\n")
        f.write(f"- measured TDC rows in expected geometry: {len(measured_expected)}\n")
        f.write(f"- expected concrete TDC rows: {len(expected)}\n")
        f.write(f"- missing concrete rows: {len(missing)}\n")
        f.write("- single global fallback rows written: 1\n")
        f.write(f"- rows outside expected geometry: {len(extra)}\n\n")
        f.write("Missing concrete rows are not filled explicitly. Instead, the output file appends one low-specificity wildcard row in a separate `[TDC]` section. Existing concrete rows remain more specific and therefore take precedence in `calibrator`.\n\n")

        f.write("## Single Global Fallback Row\n\n")
        f.write("| fallback | off | iif | source rows | iif rms |\n")
        f.write("|---|---:|---:|---:|---:|\n")
        f.write(f"| single_global_off0 | 0 | {single_global.iif:.8g} | {single_global.nsource} | {single_global.iif_rms:.8g} |\n")
        f.write("\n")

        f.write("## Phase Delta Comparison\n\n")
        f.write("This compares fallback phase against measured concrete calibration phase using deterministic fine values `60,70,...,120`. At 320 MHz, one clock is 3125 ps.\n\n")
        f.write("| fallback | mean clocks | RMS clocks | std clocks | RMS ps |\n")
        f.write("|---|---:|---:|---:|---:|\n")
        for item in phase_delta_summary(rows, expected, fallback, single_global):
            f.write("| {name} | {mean:.8g} | {rms:.8g} | {std:.8g} | {rms_ps:.3f} |\n".format(**item))
        f.write("\n")

        f.write("## Leave-One-Out Pattern Comparison\n\n")
        f.write("Each measured row is hidden and predicted from other measured rows in the same candidate pattern group. Lower RMS means the pattern is a better fallback model for measured data.\n\n")
        f.write("| pattern | groups | used | skipped | median source rows | off RMS | off p90 | iif RMS | iif p90 |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for item in pattern_comparison(rows, expected):
            f.write("| {pattern} | {groups} | {used} | {skipped} | {median_source_rows:.1f} | {off_rms:.8g} | {off_p90:.8g} | {iif_rms:.8g} | {iif_p90:.8g} |\n".format(**item))
        f.write("\n")

        f.write("## Coverage By Device\n\n")
        f.write("| device | measured | expected | missing | coverage |\n")
        f.write("|---:|---:|---:|---:|---:|\n")
        for device in sorted(exp_dev):
            measured = pre_dev[device]
            expected_n = exp_dev[device]
            missing_n = expected_n - measured
            cov = measured / expected_n if expected_n else 0.0
            f.write(f"| {device} | {measured} | {expected_n} | {missing_n} | {cov:.3f} |\n")
        f.write("\n")

        f.write("## Coverage By Device And Chip\n\n")
        f.write("| device | chip | measured | expected | missing | coverage |\n")
        f.write("|---:|---:|---:|---:|---:|---:|\n")
        for device, chip in sorted(exp_chip):
            measured = pre_chip[(device, chip)]
            expected_n = exp_chip[(device, chip)]
            missing_n = expected_n - measured
            cov = measured / expected_n if expected_n else 0.0
            f.write(f"| {device} | {chip} | {measured} | {expected_n} | {missing_n} | {cov:.3f} |\n")
        f.write("\n")

        f.write("## TDC Parameter Summary From Measured Rows\n\n")
        f.write("| tdc | n | off median | off rms | iif median | iif rms |\n")
        f.write("|---:|---:|---:|---:|---:|---:|\n")
        for tdc in range(4):
            selected = [r for k, r in rows.items() if k in expected_set and k.tdc == tdc]
            if not selected:
                continue
            offs = [r.off for r in selected]
            iifs = [r.iif for r in selected]
            f.write(f"| {tdc} | {len(selected)} | {median(offs):.8g} | {rms(offs):.8g} | {median(iifs):.8g} | {rms(iifs):.8g} |\n")
        f.write("\n")

        if missing:
            f.write("## Missing Concrete Rows Covered By Fallback\n\n")
            for k in missing:
                f.write(f"- {k.device} {k.fifo} {k.column} {k.pixel} {k.tdc}\n")
            f.write("\n")

        if extra:
            f.write("## Rows Outside Expected Geometry\n\n")
            for k in extra[:500]:
                f.write(f"- {k.device} {k.fifo} {k.column} {k.pixel} {k.tdc}\n")
            if len(extra) > 500:
                f.write(f"- ... {len(extra) - 500} more\n")


def main():
    ap = argparse.ArgumentParser(description="Add a single global [TDC] fallback row for missing concrete TDC calibrations.")
    ap.add_argument('--input', '-i', required=True, help='input calibration .conf file')
    ap.add_argument('--output', '-o', required=True, help='completed output calibration .conf file')
    ap.add_argument('--report', '-r', required=True, help='markdown report output')
    ap.add_argument('--devices', help='comma list/ranges of devices, e.g. 192..199. Default: devices observed in input')
    ap.add_argument('--fifos', help='comma list/ranges of FIFOs, e.g. 0..31. Default: FIFOs observed in input')
    args = ap.parse_args()

    rows = load_tdc(args.input)
    if not rows:
        raise SystemExit('ERROR: no [TDC] rows found')

    devices = parse_int_list(args.devices) if args.devices else sorted({k.device for k in rows})
    fifos = parse_int_list(args.fifos) if args.fifos else sorted({k.fifo for k in rows})
    expected = expected_keys(devices, fifos)
    fallback = global_tdc_rows(rows, expected)
    single_global = single_global_zero_off_row(rows, expected)
    missing = [k for k in expected if k not in rows]

    write_output(args.output, rows, expected, single_global)
    write_report(args.report, args.input, args.output, rows, expected, fallback, single_global)

    print(f"measured rows:        {len(rows)}")
    print(f"expected rows:        {len(expected)}")
    print(f"missing rows:         {len(missing)}")
    print("single global fallback rows: 1")
    print(f"output:               {args.output}")
    print(f"report:               {args.report}")


if __name__ == '__main__':
    main()
