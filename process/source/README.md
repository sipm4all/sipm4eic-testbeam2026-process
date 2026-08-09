# Compiled Processing Programs

This directory contains the C++ sources and `CMakeLists.txt` for the processing executables.

Build from the repository root:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

Executables are built in `process/build/` and installed to `process/bin/`.

## Programs

### dcalib

Runs TDC fine-time calibration for every unique ALCOR channel found in an input file and writes both ROOT diagnostics and a calibrator-compatible `[TDC]` text snippet:

```bash
process/bin/dcalib \
  --input decoded.root \
  --output dcalib.root \
  --calibration-output tdc_calibration.conf \
  --period 10000 \
  --min-pairs 100
```

For each channel identified by `device fifo column pixel`, it fits the four TDC parameters:

```text
off_0..off_3
iif_0..iif_3
```

The text output contains only a `[TDC]` section and can be concatenated with other `dcalib` outputs, provided the same concrete TDC row is not repeated. ROOT diagnostics are stored one directory per channel, each containing `hParam` and `hDelta`.

### calibrator

Creates or updates the calibrated `time` branch in the `alcor` tree:

```bash
process/bin/calibrator \
  --input decoded.root \
  --output calibrated.root \
  --config process/config/calibration/calibration_example.conf
```

For ALCOR hits:

```text
phase = off + iif * fine
time = coarse + 32768 * rollover - phase - channel_offset
```

For trigger tags:

```text
time = coarse + 32768 * rollover - trigger_offset
```

### sorter

Sorts one already decoded/calibrated single-lane ROOT file by `data.time` within spill boundaries:

```bash
process/bin/sorter --input calibrated.root --output sorted.root --window 32768
```

### after-pulse-suppressor

Suppresses hits that arrive too close to a previous hit from the same channel:

```bash
process/bin/after-pulse-suppressor --input sorted.root --output aps.sorted.root --window 50
```

Trigger tags and spill markers are propagated.

### merger

Merges multiple sorted streams and collapses duplicate spill markers:

```bash
process/bin/merger --input lane0.root lane1.root --output merged.root
```

Optional split-spill output:

```bash
process/bin/merger --input lane0.root lane1.root --output board.root --split-spills
```

This writes files such as `board.spill_0000.root`.

### trigger

Builds triggered frames using a declarative configuration:

```bash
process/bin/trigger \
  --input merged.root \
  --output frames.root \
  --config process/config/trigger/trigger_range.conf \
  --window 256
```

The persistent output tree is named `frames`. Each stored hit includes the original raw fields plus the calibrated `time` value that was used for triggering. The time is persisted in the category branches as `trigger_time`, `timing_time`, and `cherenkov_time`.

## Shared Headers

- `data_word.h`: common `data_t` representation for the `alcor` tree, optional calibrated time binding, and word-type helpers.
- `calibration.h`: parser and cached lookup API for timing calibration files.
