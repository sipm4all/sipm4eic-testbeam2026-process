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

### timing_ml.py

Experimental PyTorch study for the two TIMING scintillators. It learns:

- one residual offset per TIMING channel;
- event-by-event channel weights for TIMING0 and TIMING1;
- an estimator that minimizes the width of `timing0 - timing1`.

The script is intended for triggered-frame ROOT files, or for pre-extracted NPZ tensors:

```bash
process/bin/timing_ml.py \
  --input-root triggered.timing.root \
  --initial-calibration timing_offsets.conf \
  --output timing_ml_result.npz \
  --calibration-output timing_ml_offsets.conf
```

If `uproot` is not available, an NPZ file can be supplied with:

```bash
process/bin/timing_ml.py \
  --input-npz timing_tensor.npz \
  --initial-calibration timing_offsets.conf
```

The ROOT hit times are large absolute clock values, while useful timing differences are sub-clock. The script therefore shifts every event to a local time origin before converting to PyTorch `float32`; this is required to preserve fine timing information.

Diagnostics can be produced from a trained model:

```bash
process/bin/timing_ml.py \
  --input-npz timing_tensor.npz \
  --load-model timing_ml_model.pt \
  --epochs 0 \
  --diagnostics-output timing_ml_diagnostics.npz \
  --diagnostics-plot-output timing_ml_diagnostics.pdf
```

The diagnostics contain the learned TIMING0/TIMING1 weight pattern, first-hit-conditioned weight maps, average weight versus time from the first hit, effective channel counts, and the final `TIMING0 - TIMING1` histogram.

### timing_propagation_calib

Fits a simple isotropic light-propagation model independently for TIMING0 and TIMING1:

```bash
process/bin/timing_propagation_calib \
  --input triggered.timing.root \
  --calibration timing_offsets_from_cherenkov860.conf \
  --output timing_propagation_calib.root \
  --calibration-output timing_propagation_model.conf \
  --pitch 0.37 \
  --depth 1.0 \
  --slope 0.05 \
  --power 0.5
```

For each scintillator the fit parameters are:

- channel pitch in cm;
- effective light-emission depth in cm;
- light-propagation slope in clock/cm^power;
- light-propagation distance power.

The channel offsets are fixed to the supplied `[CHANNEL]` calibration file, typically the Cherenkov-860 timing calibration. For every TMinuit evaluation of the global propagation parameters, the emission position is refit independently in every event:

```text
per-event nuisance parameters:
    x_event
    y_event
    t0_event
```

For fixed global propagation parameters and fixed `(x_event, y_event)`, `t0_event` is solved analytically as the mean of `time_i - offset_i - propagation_i`. The program then performs a small grid-refinement search in `(x_event, y_event)` for each event. The residual minimized by TMinuit is:

```text
residual =
    time_i
  - offset_i
  - t0_event
  - slope * path_i(x_event, y_event)^power
```

with:

```text
path_i =
    sqrt((x_i - x_event)^2
       + (y_i - y_event)^2
       + depth^2)
```

where DO channel coordinates are `x_i = pitch * (DO % 4)`, `y_i = pitch * (DO / 4)`.

The ROOT output includes:

- `hResidualInitial_timing*`, residuals before the propagation-parameter fit;
- `hResidualFit_timing*`, residuals after the propagation-parameter fit;
- `hResidualVsPath_timing*`, residuals versus direct propagation path;
- `hYvsX_timing*`, fitted event-by-event emission positions;
- `hDeltaTimingEventT0`, the `TIMING0 - TIMING1` distribution using the fitted event emission time `t0_event` independently reconstructed in each scintillator;
- `hX0vsX1`, `hY0vsY1`, and `hPositionDistance`, correlations between the fitted TIMING0 and TIMING1 emission positions.

The text output is not a `[CHANNEL]` calibration file. It records the fitted propagation parameters:

```text
# timing pitch_cm depth_cm slope_clock_per_cm power
TIMING0 ...
TIMING1 ...
```

Useful controls:

```text
--max-events      limit statistics for fast studies
--max-calls       TMinuit MIGRAD call limit
--grid-steps      grid points per axis in the nested event-position fit
--min-channels    minimum channels required in an event fit
--require-trigger use only frames containing at least one trigger hit
```

Inspect the ROOT residual histograms and fitted `hYvsX_timing*` position maps to decide whether this propagation model is useful before applying it in a later estimator.

### timing_position_fit

Fits the light-emission position event by event with fixed channel offsets and fixed propagation geometry:

```bash
process/bin/timing_position_fit \
  --input triggered.timing.root \
  --calibration timing_offsets_from_cherenkov860.conf \
  --output timing_position_fit.root \
  --pitch 0.37 \
  --thickness 2.0 \
  --slope 0.05
```

For each frame and each TIMING scintillator, the program fits:

```text
x_event
y_event
t0_event
```

using:

```text
time_i =
    t0_event
  + offset_i
  + slope * 2 * sqrt((x_i - x_event)^2
                   + (y_i - y_event)^2
                   + thickness^2)
```

The channel offsets are fixed from the supplied calibration file. For fixed `(x_event, y_event)`, `t0_event` is solved analytically as the mean of `time_i - offset_i - propagation_i`; the program then performs a small grid-refinement search in `(x_event, y_event)`.

The ROOT output contains fitted position maps, per-scintillator residual RMS distributions, the `TIMING0 - TIMING1` distribution using fitted `t0_event`, and correlations between fitted TIMING0/TIMING1 positions.

## Shared Headers

- `data_word.h`: common `data_t` representation for the `alcor` tree, optional calibrated time binding, and word-type helpers.
- `calibration.h`: parser and cached lookup API for timing calibration files.
