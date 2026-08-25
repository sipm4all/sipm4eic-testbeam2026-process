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

For control words such as START_SPILL, END_SPILL, and trigger tags, fields that do not belong to the control word payload are written with placeholder value `0` rather than `-1`. ALCOR hits likewise use `counter = 0` because the hit counter field is not meaningful for them. This keeps the decoded format compatible with a possible future move to unsigned field types. The word identity must be taken from `type`, not from placeholder hit fields.

START_SPILL and END_SPILL markers written by the decoder use `fine = 0`. The decoder no longer writes artificial status-only spill markers for bad spills.

DAQ readout records marked by `0xdeadbeef` are never written as spill markers by the decoder. They are counted in the summary and skipped, because they do not provide reliable spill-boundary information.

The decoder is stateful and conservative. For each ALCOR FIFO it searches for START_SPILL (`type == 7`), then accepts data only until the matching END_SPILL (`type == 15`) with the same counter is found. END_SPILL-looking words with the wrong counter are treated as decoding errors and skipped rather than accepted as real spill boundaries. A spill is written to the output ROOT tree only after its matching END_SPILL has been found. If EOF is reached while a spill is open, that incomplete spill and its buffered payload are discarded; the decoded output remains spill-balanced and the summary records `incomplete_spills_discarded` and `incomplete_payload_words_discarded`.

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

If a completed spill has more errors than this threshold, the decoder suppresses the entire spill: it writes neither START_SPILL, payload, nor END_SPILL for that spill. The summary records `spills_suppressed_by_errors` and `payload_words_suppressed_by_errors`. This keeps decoded ROOT files internally spill-balanced and lets `checker.sh` identify FIFOs whose spill counts no longer match the rest of the run.

A separate DAQ-level suppression case is recognized when a buffer payload is exactly:

```text
START_SPILL 0xdeadbeef END_SPILL 0xdeadbeef
```

in the raw buffer. This record is counted as `daq_suppressed_records` and skipped. It is not written to the output tree.

The decoder prints summary counters including spills found/written, completed spills suppressed by errors, incomplete spills discarded at EOF, DAQ-suppressed records, wrong END_SPILL candidates, invalid-column hits, malformed words, and skipped words outside spills. The same summary is also written to a sidecar text file derived from the output name: `output.root` produces `output.summary`.

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

### coordinator

Creates or updates spatial coordinate branches in the `alcor` tree:

```bash
process/bin/coordinator \
  --input calibrated.root \
  --output coordinated.root
```

The program preserves the input tree entry order and writes exactly one output entry per input entry. It adds or recalculates:

```text
x
y
```

For ALCOR hits with `device == 200`, the coordinates are TIMING scintillator coordinates. The electronics channel is:

```cpp
eoch = pixel + 4 * column
```

and is converted to detector-output channel with the TIMING `eo2do` table. The TIMING pitch is:

```text
3.5 mm
```

so:

```cpp
x = 3.5 * (DO % 4)
y = 3.5 * (DO / 4)
```

For ALCOR hits with `device != 200`, the coordinates use the TESTBEAM2026 Cherenkov mapping from the legacy `mapping.h` logic:

```text
device, fifo -> chip -> PDU/matrix
column, pixel -> electronics-oriented channel
matrix mapping -> detector channel
detector channel -> global x/y
```

The Cherenkov electronics-oriented channel follows the same decoded convention as `data_word.h`:

```cpp
local_eoch = pixel + 4 * column
eoch = local_eoch + 32 * ((fifo / 4) % 2)
```

Control words and trigger tags are preserved with:

```text
x = 0
y = 0
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

The merger is strict about spill alignment. At every START_SPILL and END_SPILL boundary, all input streams must be positioned on the same marker type with the same `counter`. The END_SPILL counter must also match the current START_SPILL counter. If a FIFO is missing a spill, for example because the decoder skipped unrecoverable `0xdeadbeef` DAQ-suppressed records, the merger fails loudly and prints each stream state.

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

The persistent output tree is named `frames`. Each stored hit includes the original raw fields plus the calibrated `time` value used for triggering and the spatial coordinates `x` and `y` when they were added upstream by `coordinator`. These are persisted in category branches such as `trigger_time`, `timing_time`, `cherenkov_time`, and similarly `trigger_x`, `timing_x`, `cherenkov_x`, `trigger_y`, `timing_y`, `cherenkov_y`.

### timing

Adds trained TIMING event-estimator results to an existing triggered-frame ROOT file:

```bash
process/bin/timing \
  --input triggered.root \
  --output triggered.timing.root
```

The input must be the synchronized multi-tree output produced by `trigger`, with one entry per accepted frame in `frames`, `trigger`, `timing`, and `cherenkov`. The calibrated timing values are read from the `time[nhits]` branch of the `timing` tree. The program preserves all input trees and adds these scalar branches to the `timing` tree:

```text
timing_valid
T0
sigma0
T1
sigma1
T
sigmaT
```

`T0` and `sigma0` refer to TIMING0, `T1` and `sigma1` refer to TIMING1, and `T = (T0 + T1) / 2`. Values are in the estimator native time unit, `3.125 ns`. `timing_valid` is `1` only when the frame contains one usable time for every TIMING DO channel in both scintillators. If one or more channels are missing, the timing values are written as `NaN` and `timing_valid` is `0`.

### ring-finder

Finds zero or more Cherenkov ring candidates in each frame using calibrated
Cherenkov `x`, `y`, and trigger-relative `time` values:

```bash
process/bin/ring-finder \
  --input triggered.timing.root \
  --output triggered.rings.root \
  --time-window 5
```

The finder obtains spatial circle candidates with RANSAC. RANSAC is used only
to identify candidate inliers; the circle is then refit in the least-squares
sense using all accepted inliers. The candidate time is their arithmetic mean,
and only hits passing both the radial tolerance and ring time window are kept.
Accepted inliers are removed and the search repeats, allowing multiple
separated rings per frame.

The output preserves the synchronized input trees and adds a `ring` tree with
one entry per frame:

```text
nring
ring_x0[nring]
ring_y0[nring]
ring_r[nring]
ring_e[nring]
ring_phi[nring]
ring_time[nring]
ring_ninliers[nring]
```

The `ring` entry at index `i` corresponds to the `frames`, `trigger`,
`timing`, and `cherenkov` entries at index `i`. `ring_r` is the semi-major axis,
`ring_e` is the eccentricity, `ring_phi` is the rotation angle in radians,
and `ring_time` is the fitted ring time in native time units.
Setting `ring_e=0` gives a circle. Defaults are a radial tolerance of `5`,
a time window of `5` native units, and `8` rings maximum. Use
`--tolerance`, `--time-window`, `--min-inliers`, `--iterations`, and
`--max-rings` to change them. Candidate centers and radii are constrained by
`--min-x0`, `--max-x0`, `--min-y0`, `--max-y0`, `--min-radius`, and
`--max-radius` (defaults: tolerance `5` mm, center in `[-100,100]`, and
radius in `[1,200]`).
There is no eccentricity cut: valid ellipse candidates are retained, while
the original circle is used as a fallback if ellipse refinement loses inliers.

### offset_calibrator

Fits the peak position in each global-channel projection of the `hDeltaT` histogram produced by `macros/example/deltat.C` and writes channel timing offsets:

```bash
process/bin/offset_calibrator \
  --input deltat.root \
  --output timing_offsets.conf \
  --min-entries 1000
```

The calibration output contains a `[CHANNEL]` section. The convention is:

```text
calibrated_time = raw_time - offset
```

The tool is currently restricted to Cherenkov global channels `0..2047` (devices `192..199`). Channels with at least `--min-entries` entries in the `[-5,5]` delta-t interval receive a local Gaussian fit around their highest bin. The default is `1000` entries. A fit is accepted only when its peak-position uncertainty is no larger than `0.1` clock units. Channels failing either criterion use the median fitted position from the repeated hardware, in this order:

```text
channel position within PDU:  global_channel % 256
channel position within chip: global_channel % 32
channel position within FIFO: global_channel % 8
global median
```

The companion diagnostic file `<output>.root` contains `hOffset`, `hEntries` (the `[-5,5]` count), `hUncertainty`, `hSourceLevel`, the original `hDeltaT`, the shifted `hDeltaT_corrected`, and the `offsets` summary tree. The corrected histogram applies `delta_t_corrected = delta_t - offset` independently for every channel. Source levels are `256`, `32`, `8`, and `1` for PDU, chip, FIFO, and global fallback respectively; `0` means a direct fit. Inspect this diagnostic before using the generated calibration file in production.

## Shared Headers

- `data_word.h`: common `data_t` representation for the `alcor` tree, optional calibrated time binding, and word-type helpers.
- `calibration.h`: parser and cached lookup API for timing calibration files.
- `geometry.h`: TESTBEAM2026 spatial mapping helpers used by `coordinator`.

## TIMING Estimator

The self-contained trained TIMING event estimator lives under:

```text
timing/
```

It provides the native C++ `TimingEstimator` API and keeps the original Python package as the validation reference. See `timing/README.md` for the input order, calibration convention, model conversion details, and validation procedure.
