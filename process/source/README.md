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

The persistent output tree is named `frames`.

## Shared Headers

- `data_word.h`: common `data_t` representation for the `alcor` tree, optional calibrated time binding, and word-type helpers.
- `calibration.h`: parser and cached lookup API for timing calibration files.
