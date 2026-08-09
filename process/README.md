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
process/bin/timing_calib   --input triggered.timing.root   --output timing_calib.root   --calibration-output timing_channel_offsets.conf   --pre-iterations 5   --minimizer-calls 5000
```

For quick validation runs, `--max-frames N` limits how many frames are read. Use `--max-frames 0` or omit the option for the full file.

The text output is a `[CHANNEL]` calibration snippet using the same sign convention as `calibrator`:

```cpp
calibrated_time = raw_time - offset
```

It can be included or concatenated into a calibration configuration used by a later `calibrator` pass.
