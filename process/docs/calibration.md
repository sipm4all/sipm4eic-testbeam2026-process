# Calibration Data Notes

This document collects working notes for calibration data taken during the 2026 beam test. It is based on logbook information plus practical command lines used to run the calibration scripts in this repository.

The main goal of the runs listed here is TDC calibration with `dcalib.sh`. The generated text snippets contain `[TDC]` rows that can later be combined with channel and trigger offsets to build a full calibrator configuration.

## Laser / Test-Pulse Runs

The relevant laser runs were recorded on `2026-06-18`. The logbook labels them as laser runs at 10 kHz with:

```text
opMode = 1
deltaThr = 10
```

However, from data analysis the observed rate appears to be about 50 kHz. At a 320 MHz clock frequency this corresponds to:

```text
clock period = 3.125 ns
50 kHz period = 20 us
20 us / 3.125 ns = 6400 clock cycles
```

Therefore the `dcalib.sh` commands below use:

```bash
--period 6400
```

The pulser was configured with:

```text
output 1 -> SPILL signals
output 2 -> laser triggers
```

The two pulser outputs were aligned. For each RDO, the electronics were configured so that one half of the chips received a test pulse while the other half saw laser light.

## ALCOR OpModes Used In These Runs

The run summary below uses ALCOR `opMode` values. For the calibration runs, the relevant meanings are:

```text
opMode = 1
  Normal SiPM signal-detection mode.
  Leading Edge (LE) detection and digitisation are enabled.

opMode = 4
  Normal SiPM signal-detection mode in ToT configuration.
  Both Leading Edge (LE) and Trailing Edge (TE) detection are enabled.

opMode = 2
  Test-pulse mode.
  A signal from the TESTPULSE input of the chip is sent directly into
  the digital part of the chip.
```

Here "normal" means that the chip detects real signals arriving from the SiPM inputs. In contrast, `opMode = 2` bypasses the analogue SiPM-input signal path and injects a controlled test pulse into the digital/TDC part of ALCOR.

For TDC calibration we use the chips in `opMode = 2`. Through the SPILL/pulser machinery, fixed-period pulses are sent to the ALCOR TDCs. Conceptually, one channel receives a periodic pulse train; in the real run this happens for many channels in parallel.

Each ALCOR channel has four TDCs. Consecutive pulses in one channel are digitised by these TDCs one after another. Because the injected pulse period is known, the measured time difference between consecutive hits should reproduce that period after applying the correct TDC fine-time calibration.

## What `dcalib` Fits

`dcalib` reads the `alcor` tree and groups ALCOR hits by electronics channel:

```text
device, fifo, column, pixel
```

Within each channel it builds pairs of consecutive hits, resetting the history at spill boundaries. For each pair it stores:

```text
current tdc, current fine
previous tdc, previous fine
coarse_delta = current coarse clock - previous coarse clock
```

Hits with invalid TDC indices are ignored. Pairs whose coarse difference is far from the expected period are not used in the chi-square.

For each channel, `dcalib` fits nine parameters:

```text
off_0, off_1, off_2, off_3
iif_0, iif_1, iif_2, iif_3
period
```

The fine-time phase for one TDC hit is:

```text
phase = off_tdc + iif_tdc * fine
```

For a pair of consecutive hits, the corrected time difference is:

```text
time_delta = coarse_delta - (phase_current - phase_previous)
```

The minimisation adjusts the TDC parameters so that:

```text
time_delta ~= fitted period
```

with the fitted period expected to be close to the command-line `--period` value. The first TDC offset, `off_0`, is fixed to remove the arbitrary common offset of the four TDCs in one channel. The output text file contains one `[TDC]` row for each fitted `(device, fifo, column, pixel, tdc)`.

## Run Summary

```text
20260618-183625
  chips 0,1,2,3: opMode = 2
  chips 4,5,6,7: opMode = 1

20260618-184325
  chips 0,1,2,3: opMode = 2
  chips 4,5,6,7: opMode = 4

20260618-185127
  chips 4,5,6,7: opMode = 2
  chips 0,1,2,3: opMode = 1

20260618-185857
  chips 4,5,6,7: opMode = 2
  chips 0,1,2,3: opMode = 4
```

For the processing scripts, chip number maps to FIFO groups as:

```text
chip = fifo / 4

chips 0,1     -> fifos 0..7
chips 0..3    -> fifos 0..15
chips 4..7    -> fifos 16..31
```

## Running TDC Calibration

The `dcalib.sh` workflow processes decoded test-pulse files from:

```text
/data/2026-testbeam/actual/testpulse/<run>/<device>/decoded/
```

and writes outputs under:

```text
/data/2026-testbeam/process/<run>/<device>/
```

For each FIFO input file, the workflow is:

```text
clean -> sort -> dcalib
```

The final calibration snippet is named, for example:

```text
dcalib.fifo_0.conf
```

Each snippet contains only a `[TDC]` section and can be used as an input fragment when assembling a full calibration file, provided duplicate concrete TDC rows are avoided.

## TIMING Readout Calibration

This command calibrates TDCs for the TIMING readout on `kc705-200`, using run `20260618-183625`.

The TIMING readout uses chips 0 and 1, corresponding to FIFOs `0..7`.

```bash
sipm4eic-testbeam2026-process/process/scripts/dcalib.sh \
    --run 20260618-183625 \
    --devices kc705-200 \
    --fifos {0..7} \
    --period 6400
```

## Cherenkov Readout Calibration: Chips 0..3

This command calibrates TDCs for chips `0,1,2,3` of the Cherenkov readout on all RDO devices, using run `20260618-183625`.

Chips `0..3` correspond to FIFOs `0..15`.

```bash
sipm4eic-testbeam2026-process/process/scripts/dcalib.sh \
    --run 20260618-183625 \
    --devices rdo-{192..199} \
    --fifos {0..15} \
    --period 6400
```

Run `20260618-184325` can also be used for the same chip group and should give compatible results.

## Cherenkov Readout Calibration: Chips 4..7

This command calibrates TDCs for chips `4,5,6,7` of the Cherenkov readout on all RDO devices, using run `20260618-185127`.

Chips `4..7` correspond to FIFOs `16..31`.

```bash
sipm4eic-testbeam2026-process/process/scripts/dcalib.sh \
    --run 20260618-185127 \
    --devices rdo-{192..199} \
    --fifos {16..31} \
    --period 6400
```

Run `20260618-185857` can also be used for the same chip group and should give compatible results.

## Merging The TDC Calibration Snippets

After the three `dcalib.sh` runs above, the output directory contains many per-FIFO files named like:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib.fifo_0.conf
```

Each file contains a `[TDC]` fragment. These fragments can be merged with:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh
```

The examples below first merge the output of each calibration campaign into one file, then merge those three files into one common TDC calibration file.

### Merge TIMING TDC Calibration

This merges the `kc705-200` TIMING calibration from run `20260618-183625`:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/20260618-183625/kc705-200/dcalib.fifo_{0..7}.conf \
    --output /data/2026-testbeam/process/tdc.timing.conf
```

### Merge Cherenkov TDC Calibration: Chips 0..3

This merges the RDO Cherenkov calibration for chips `0..3`, using run `20260618-183625`:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/20260618-183625/rdo-{192..199}/dcalib.fifo_{0..15}.conf \
    --output /data/2026-testbeam/process/tdc.cherenkov.chips_0_3.conf
```

### Merge Cherenkov TDC Calibration: Chips 4..7

This merges the RDO Cherenkov calibration for chips `4..7`, using run `20260618-185127`:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/20260618-185127/rdo-{192..199}/dcalib.fifo_{16..31}.conf \
    --output /data/2026-testbeam/process/tdc.cherenkov.chips_4_7.conf
```

### Merge All TDC Calibrations

This combines the three merged outputs into one common TDC calibration snippet:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/tdc.timing.conf \
    --input /data/2026-testbeam/process/tdc.cherenkov.chips_0_3.conf \
    --input /data/2026-testbeam/process/tdc.cherenkov.chips_4_7.conf \
    --output /data/2026-testbeam/process/tdc.20260618.conf
```

The merged `tdc.20260618.conf` contains only measured concrete `[TDC]` rows. Some channels may be missing because not all FIFO calibration snippets were available or not all channels had enough statistics. Complete the TDC file with one low-specificity global fallback row:

```bash
sipm4eic-testbeam2026-process/process/scripts/complete_tdc_calibration.py \
    --input /data/2026-testbeam/process/tdc.20260618.conf \
    --output sipm4eic-testbeam2026-process/process/config/calibration/tdc.20260618.conf \
    --report sipm4eic-testbeam2026-process/process/config/calibration/tdc.20260618.md
```

The completed output keeps all measured concrete rows and appends a separate `[TDC]` section containing a single fallback row:

```text
* * * * * 0 <global_iif>
```

The wildcard row has lower specificity than every measured concrete calibration row, so it is used only when a concrete TDC calibration is absent. The generated report documents missing rows, coverage by device/chip, leave-one-out pattern comparisons, and the phase error expected from the fallback.

The final `process/config/calibration/tdc.20260618.conf` still contains only the `[TDC]` calibration section. To run `calibrator`, this TDC file must be combined with `[CHANNEL]` and `[TRIGGER]` sections, or with explicit wildcard defaults for those sections.

## Notes And Checks

- The commands above assume the executables have already been built and installed with CMake.
- `dcalib.sh` currently hardcodes `MIN_PAIRS=1000`.
- `dcalib.sh` currently hardcodes `SORT_WINDOW=32768` for the sorting step before TDC calibration.
- The TDC calibration should be derived from otherwise uncalibrated data. The sorter stage may compute nominal hit time internally for ordering, but it does not create or persist a `time` branch when the input does not already have one.
- The output `.conf` fragments are per-FIFO `[TDC]` snippets. A complete calibration file still needs `[CHANNEL]` and `[TRIGGER]` sections, or suitable defaults, before it can be used by `calibrator`.
