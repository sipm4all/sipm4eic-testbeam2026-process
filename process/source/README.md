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


### decoder

Decodes one raw per-FIFO `.dat` file into a ROOT file containing an `alcor` tree:

```bash
process/bin/decoder \
  --input alcdaq.fifo_13.dat \
  --output alcdaq.fifo_13.root \
  --allowed-spill-errors 0
```

The output tree branches are compatible with the downstream processing chain:

```text
device fifo type counter column pixel tdc rollover coarse fine
```

For control words such as START_SPILL, END_SPILL, and trigger tags, fields that do not belong to the control word payload are written with placeholder value `0` rather than `-1`. ALCOR hits likewise use `counter = 0` because the hit counter field is not meaningful for them. This keeps the decoded format compatible with a possible future move to unsigned field types. The word identity must be taken from `type`, not from placeholder hit fields. START_SPILL and END_SPILL use `fine = 0` for normal spills. If the decoder suppresses a spill payload because the spill exceeded `--allowed-spill-errors`, both the START_SPILL and END_SPILL words are written with `fine = 1`.

The decoder is stateful and conservative. For each ALCOR FIFO it searches for START_SPILL (`type == 7`), then accepts data only until the matching END_SPILL (`type == 15`) with the same counter is found. END_SPILL-looking words with the wrong counter are treated as decoding errors and skipped rather than accepted as real spill boundaries.

ALCOR hit words are accepted only when their decoded column is valid for the FIFO:

```text
column = 2 * (fifo % 4)
column = 2 * (fifo % 4) + 1
```

This rejects garbage words that decode into impossible FIFO/column combinations. Rollover markers `0x5c5c5c5c` are counted internally and are not written as hits.

`--allowed-spill-errors` controls spill-level payload quarantine. The default is strict:

```text
--allowed-spill-errors 0
```

If a completed spill has more errors than this threshold, the decoded output still writes the START_SPILL and END_SPILL markers for that spill, but omits all payload hits/tags from that spill. Both spill markers get `fine = 1` to tag the spill as suppressed. This preserves spill boundaries for downstream synchronization while preventing corrupted payload from entering the processing chain.

The decoder prints summary counters including spills found/written/emptied, wrong END_SPILL candidates, invalid-column hits, malformed words, and skipped words outside spills. The same summary is also written to a sidecar text file derived from the output name: `output.root` produces `output.summary`.

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

### cleaner

Clones the input `alcor` tree to a new ROOT file while dropping entries that fail the stream-cleaning predicate:

```bash
process/bin/cleaner --input decoded.root --output cleaned.root
```

The output tree preserves the input tree name and branch layout through ROOT `CloneTree(0)`. The current cleaning criterion removes malformed ALCOR hits whose `column` is not valid for the corresponding `fifo`. For a given FIFO, the allowed columns are `2 * (fifo % 4)` and `2 * (fifo % 4) + 1`. Spill markers, trigger tags, and other non-ALCOR words are preserved. The program checks ROOT `GetEntry()` return values and prints input/kept/dropped/output entry counts, including the number dropped by the column test.


### checker

Runs a read-only sanity check over one decoded per-FIFO `alcor` tree and writes an ASCII report:

```bash
process/bin/checker --input decoded.root --output decoded.check
```

The report contains the input entry count, counts of known word types, unknown-word count, and spill-boundary checks:

```text
start_spill_type7
end_spill_type15
alcor_hits_type1
trigger_tags_type9
unknown_words
spill_counter_consistent
open_spill_at_eof
spill_count_balance
errors
```

The current DAQ control-word convention is `type == 7` for START_SPILL and `type == 15` for END_SPILL. START and END spill counters are expected to pair with the same `counter` value, and counters are expected to increase by one spill-to-spill.

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
