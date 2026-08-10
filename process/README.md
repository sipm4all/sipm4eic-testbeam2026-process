# Processing Pipeline

The `process/` directory contains the compiled data-processing chain, workflow scripts, and configuration files.

Recommended order:

```text
decoded ROOT file
  -> calibrator
  -> sorter
  -> after-pulse-suppressor
  -> merger
  -> trigger
  -> triggered frame ROOT files
```

All compiled programs operate on ROOT files. The common raw input tree name is `alcor`. Triggered output files contain a `frames` tree with one entry per spill and flattened frame contents split into trigger, timing, and Cherenkov collections.

Directory contents:

```text
source/                 C++ programs and CMake build file
scripts/                shell scripts for larger batch workflows
config/calibration/     timing calibration files
config/trigger/         trigger-definition files
```

The key design point is that low-level timing calibration happens once, early, in `calibrator`. Downstream tools prefer the calibrated `time` branch when present and only use the old nominal timing expression for legacy uncalibrated files.

## Timing-Channel Calibration

`timing_calib` is a compiled version of the timing-channel calibration workflow. It reads triggered-frame ROOT files, uses only the `timing` hit collection, and derives 64 channel offsets for the two timing scintillators:

```text
TIMING0: device 200, FIFOs 0..3
TIMING1: device 200, FIFOs 4..7
```

The local timing channel is:

```cpp
channel = pixel + 4 * column
```

The program first performs iterative residual-to-scintillator-mean pre-calibration, keeping `fifo=0 column=0 pixel=0` fixed at zero after every iteration. It then runs a full 63-parameter `TMinuit` minimization. The objective includes both `timing0_mean - timing1_mean` and intra-scintillator channel residual terms, so channel residual maps are constrained during the global fit.

Example:

```bash
process/bin/timing_calib   --input triggered.timing.root   --output timing_calib.root   --calibration-output timing_channel_offsets.conf   --pre-iterations 3   --minimizer-calls 5000
```

For quick validation runs, `--max-frames N` limits how many frames are read. Use `--max-frames 0` or omit the option for the full file.

The text output is a `[CHANNEL]` calibration snippet using the same sign convention as `calibrator`:

```cpp
calibrated_time = raw_time - offset
```

It can be included or concatenated into a calibration configuration used by a later `calibrator` pass.


Timing calibration diagnostics now also include event-by-event spread histograms:

```text
hTiming0SpreadBefore / hTiming0SpreadAfter
hTiming1SpreadBefore / hTiming1SpreadAfter
hExpectedDeltaFromSpreadBefore / hExpectedDeltaFromSpreadAfter
```

For each accepted frame, the timing spread is the RMS of selected channel times inside one scintillator. The expected contribution to `RMS(timing0_mean - timing1_mean)` from independent channel jitter is estimated as:

```cpp
sqrt(spread0^2 / n0 + spread1^2 / n1)
```

Comparing this estimate with the observed `hDelta` RMS is a useful check of whether the two timing scintillators are limited by independent channel jitter or by an event-by-event effect common to many channels.

The timing diagnostic also writes first-order event-shape correlation histograms:

```text
hDeltaShapeCorrected
hDeltaVsSpread0 / hDeltaVsSpread1
hDeltaVsSlopeX0 / hDeltaVsSlopeX1 / hDeltaVsSlopeY0 / hDeltaVsSlopeY1
hDeltaVsLeftRight0 / hDeltaVsLeftRight1
```

`hDeltaShapeCorrected` subtracts a linear diagnostic model of `timing0_mean - timing1_mean` using per-event shape variables measured inside the two scintillators: channel spread, channel-time slope, left/right asymmetry, and even/odd asymmetry. This correction is diagnostic only; it is not written into the `[CHANNEL]` calibration snippet, because it represents event topology or light-propagation information rather than fixed channel offsets.

The event-shape diagnostics use the detector-channel geometry, not raw electronics order. The electronics channel is first mapped through:

```cpp
eo2do[32] = {22,20,18,16,24,26,28,30,25,27,29,31,23,21,19,17,
             9,11,13,15,7,5,3,1,6,4,2,0,8,10,12,14}
```

The detector channel is interpreted as a 4x8 matrix with `x = detector_channel % 4` and `y = detector_channel / 4`. The diagnostic variables include spread, x/y time gradients, and left/right asymmetry in this detector coordinate system.

The diagnostic ROOT output also includes first-hit timing comparisons:

```text
hDeltaFirstBefore
hDeltaFirstAfter
```

For these histograms, each scintillator time is the earliest selected channel time in that frame, rather than the robust average of selected channels. The same outlier-selected channel set is used, so this compares the timing estimator itself rather than changing event selection.

First-hit detector position diagnostics are also written:

```text
hDeltaVsFirstX0 / hDeltaVsFirstX1
hDeltaVsFirstY0 / hDeltaVsFirstY1
hFirstPosition0 / hFirstPosition1
hDeltaShapePositionCorrected
```

The first hit is interpreted as a beam-position proxy inside the 4x8 detector matrix. `hDeltaShapePositionCorrected` extends the diagnostic shape correction with the first-hit detector `(x,y)` coordinates for both scintillators.

First-hit channel correlations are stored as:

```text
hFirstDetectorChannel0Vs1
hFirstElectronicsChannel0Vs1
```

`hFirstDetectorChannel0Vs1` uses the `eo2do` detector-channel mapping and is the preferred plot for checking whether the earliest-light position proxy is correlated between TIMING0 and TIMING1. `hFirstElectronicsChannel0Vs1` is kept as a mapping cross-check.
