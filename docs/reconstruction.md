# Physics Reconstruction Workflow

This document describes the standard workflow for reconstructing physics data with the repository scripts.

The intended sequence is:

```text
decoder.sh
  -> checker.sh
  -> process.sh
  -> trigger.sh
  -> timing.sh
```

Conceptually this expands to:

```text
raw per-FIFO data
  -> decoded ROOT files
  -> read-only sanity checks
  -> clean / calibrate / coordinate / sort / after-pulse suppress
  -> device-level split-spill merge
  -> run-level split-spill merge
  -> trigger frame production
  -> hadd per-trigger output
  -> optional trained TIMING estimator
```

The scripts use `/data/2026-testbeam/process/<run>/` as the working area. For physics data, `--run-type physics` is the default and can be omitted.

## Build First

From the repository root:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

This installs the executables in:

```text
process/bin/
```

## 1. Decode Raw Data

Decode the raw per-FIFO `.dat` files:

```bash
sipm4eic-testbeam2026-process/process/scripts/decoder.sh \
    --run <run>
```

If decoded ROOT files already exist, `decoder.sh` skips them by default. To regenerate existing decoded files, pass:

```bash
--overwrite
```

For physics data this reads:

```text
/data/2026-testbeam/actual/physics/<run>/<device>/raw/alcdaq.fifo_<N>.dat
```

and writes:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_<N>.root
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_<N>.summary
```

The decoder is intentionally strict. It writes only complete spills with matching START/END spill counters. Malformed completed spills are suppressed entirely when their error count exceeds `--allowed-spill-errors`, which defaults to zero. DAQ-readout-suppressed `0xdeadbeef` records are counted in the `.summary` file and skipped.

Use filters only when intentionally processing a subset:

```bash
sipm4eic-testbeam2026-process/process/scripts/decoder.sh \
    --run <run> \
    --devices kc705-200 rdo-{192..199} \
    --fifos {0..31}
```

## 2. Check Decoded Data

Run the checker on the decoded ROOT files:

```bash
sipm4eic-testbeam2026-process/process/scripts/checker.sh \
    --run <run>
```

This writes per-FIFO checks:

```text
/data/2026-testbeam/process/<run>/<device>/check/alcdaq.fifo_<N>.check
```

device-level checks:

```text
/data/2026-testbeam/process/<run>/<device>/check/<device>.check
```

and the run-level check:

```text
/data/2026-testbeam/process/<run>/check/<run>.check
```

Inspect the run-level file before continuing:

```bash
cat /data/2026-testbeam/process/<run>/check/<run>.check
```

Important fields are:

```text
consistent
errors
spill_count_uniform_start
spill_count_uniform_end
problem_check
```

A clean physics run should normally have:

```text
consistent: yes
errors: 0
spill_count_uniform_start: yes
spill_count_uniform_end: yes
```

If the run-level file lists `problem_check:` entries, inspect those per-device or per-FIFO reports before processing. The checker is read-only: it does not remove or modify decoded ROOT files.

## 3. Process And Merge

Run the calibrated processing workflow:

```bash
sipm4eic-testbeam2026-process/process/scripts/process.sh \
    --run <run> \
    --calibration <calibration.conf>
```

For each decoded FIFO file, `process.sh` runs:

```text
calibrator -> coordinator -> sorter -> after-pulse-suppressor
```

in parallel within each device.

`coordinator` is intentionally placed after `calibrator` and before `sorter`. The calibrated `time` branch and the spatial `x/y` branches are both per-hit enrichments of the original `alcor` stream. Running them before sorting keeps the later sorting, merging, and triggering stages focused on time ordering and event construction rather than detector geometry. When `trigger.sh` later builds triggered frames, `trigger.cc` copies `time`, `x`, and `y` into the split trigger/timing/Cherenkov frame collections.

It then merges each device with `merger --split-spills`, producing device-level split-spill files:

```text
/data/2026-testbeam/process/<run>/<device>/process/aps.sorted.spill_0000.root
/data/2026-testbeam/process/<run>/<device>/process/aps.sorted.spill_0001.root
...
```

Finally, it performs the run-level merge one spill at a time, producing:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0000.root
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0001.root
...
```

These run-level split-spill files are the input to `trigger.sh`.

By default `process.sh` does not overwrite existing outputs. Use:

```bash
--overwrite
```

only when intentionally regenerating existing products.

## 4. Build Triggered Frames

Run one or more trigger configurations on the merged run-level spill files:

```bash
sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
    --run <run> \
    --trigger process/config/trigger/<trigger.conf> <tag> \
    --window 256
```

Multiple trigger definitions can be run in one command:

```bash
sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
    --run <run> \
    --trigger process/config/trigger/trigger_range.conf range \
    --trigger process/config/trigger/trigger_set.conf set \
    --window 256
```

For each tag, `trigger.sh` runs `trigger` over every merged spill file:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0000.root
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0001.root
...
```

and writes temporary per-spill triggered files:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.spill_0000.root
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.spill_0001.root
...
```

It then combines them with `hadd` into:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.root
```

By default, triggered spill files are kept:

```text
CLEAN_TRIGGERED_SPILLS=0
```

Pass:

```bash
--clean-triggered-spills
```

to remove `triggered.<tag>.spill_*.root` after a successful `hadd`.

## 5. Add TIMING Estimator

After `trigger.sh`, the triggered file can be augmented with the trained TIMING estimator:

```bash
sipm4eic-testbeam2026-process/process/scripts/timing.sh \
    --run <run> \
    --trigger <tag>
```

By default, `timing.sh` reads:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.root
```

and writes:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.root
```

This stage adds per-frame `timing_valid`, `T0`, `sigma0`, `T1`, `sigma1`, `T`, and `sigmaT` branches. It preserves the existing triggered-frame tree and the `spill_participation` metadata.

For large triggered files, `timing.sh --parallel-spills` can instead run on the per-spill triggered files left by `trigger.sh`, then combine the timing-augmented spill files with `hadd`.

## 6. Check Delta-T Distributions

After triggered frames have been built, run the example `deltat.C` macro as a first reconstruction check.

For the finger-trigger workflow, use the trigger tag from device `200` as the reference:

```text
reference = {9, 200}
```

In the compact selector syntax used by `deltat.C`, this means:

```text
type   = 9
device = 200
fifo   = 32
column = *
pixel  = *
```

To compare all ALCOR hits against that trigger reference:

```bash
root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/<run>/trigger/triggered.fingers.root\", channel_selector_t(1, -1), channel_selector_t(9, 200), \"/data/2026-testbeam/process/<run>/trigger/deltat.fingers.trigger.root\")"
```

This fills:

```text
hDeltaT        delta_t = time_ALCOR - time_TRIGGER versus global channel index
hDeltaT_spill  delta_t versus spill id
hDeltaT_tdc0   delta_t versus fine for target TDC 0
hDeltaT_tdc1   delta_t versus fine for target TDC 1
hDeltaT_tdc2   delta_t versus fine for target TDC 2
hDeltaT_tdc3   delta_t versus fine for target TDC 3
```

The target selector `channel_selector_t(1, -1)` means all type-1 ALCOR hits, with `-1` used as the wildcard channel index. The reference selector `channel_selector_t(9, 200)` means trigger tags from device `200`.

If the triggered file has first been augmented with the trained TIMING estimator, `deltat.C` can also use the reconstructed TIMING as the reference:

```bash
sipm4eic-testbeam2026-process/process/scripts/timing.sh \
    --run <run> \
    --trigger fingers

root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/<run>/trigger/triggered.fingers.timing.root\", channel_selector_t(1, -1), timing_reference_t(\"T\"), \"/data/2026-testbeam/process/<run>/trigger/deltat.fingers.timing.root\")"
```

The timing-reference choices are `timing_reference_t("T")` for the combined estimator, `timing_reference_t("T0")` for TIMING0, and `timing_reference_t("T1")` for TIMING1. This mode uses only frames with `timing_valid == 1`.

## Normal Physics Example

Run `20260623-185238` is a useful physics example for this workflow. It was taken with an 11 GeV/c negative beam.

A full reconstruction of that run looks like:

```bash
run=20260623-185238
calib=sipm4eic-testbeam2026-process/process/config/calibration/baseline.conf

sipm4eic-testbeam2026-process/process/scripts/decoder.sh \
    --run "${run}"

sipm4eic-testbeam2026-process/process/scripts/checker.sh \
    --run "${run}"

cat "/data/2026-testbeam/process/${run}/check/${run}.check"

sipm4eic-testbeam2026-process/process/scripts/process.sh \
    --run "${run}" \
    --calibration "${calib}"

sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
    --run "${run}" \
    --trigger sipm4eic-testbeam2026-process/process/config/trigger/fingers.conf fingers \
    --window 256

root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/${run}/trigger/triggered.fingers.root\", channel_selector_t(1, -1), channel_selector_t(9, 200), \"/data/2026-testbeam/process/${run}/trigger/deltat.fingers.trigger.root\")"

sipm4eic-testbeam2026-process/process/scripts/timing.sh \
    --run "${run}" \
    --trigger fingers

root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/${run}/trigger/triggered.fingers.timing.root\", channel_selector_t(1, -1), timing_reference_t(\"T\"), \"/data/2026-testbeam/process/${run}/trigger/deltat.fingers.timing.root\")"
```

The final triggered output is:

```text
/data/2026-testbeam/process/20260623-185238/trigger/triggered.fingers.root
```

The delta-t check output is:

```text
/data/2026-testbeam/process/20260623-185238/trigger/deltat.fingers.trigger.root
/data/2026-testbeam/process/20260623-185238/trigger/deltat.fingers.timing.root
```

## Device Subsets

For normal physics reconstruction, run without `--devices` so all decoded devices in the run are merged together.

Use `--devices` only for intentional subset workflows, such as diagnostics:

```bash
sipm4eic-testbeam2026-process/process/scripts/process.sh \
    --run <run> \
    --devices rdo-192 \
    --calibration <calibration.conf>

sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
    --run <run> \
    --devices rdo-192 \
    --trigger <trigger.conf> diagnostic_rdo-192 \
    --window 256
```

In this mode, `process.sh` writes a subset-specific prefix:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.rdo-192.spill_0000.root
```

and `trigger.sh --devices rdo-192` consumes that same prefix.

## Output Content

The final triggered ROOT files contain a `frames` tree with one entry per spill. Each spill contains a common `nframes`, and each frame is stored in three flattened hit collections:

```text
trigger
timing
cherenkov
```

The output also propagates the `spill_participation` tree, which records the `(device,fifo)` sources that contributed to each merged spill. This is important when a decoded input stream has suppressed bad spills and therefore does not participate in every run-level merged spill.

For ROOT analysis, prefer the helper:

```cpp
#include "trigger_reader.h"
```

from:

```text
macros/lib/trigger_reader.h
```

See:

```text
macros/example/trigger_reader.C
macros/example/deltat.C
```

for minimal analysis examples.
