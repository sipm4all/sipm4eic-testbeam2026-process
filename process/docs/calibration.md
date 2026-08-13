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

The full calibration preparation starts from raw test-pulse files under:

```text
/data/2026-testbeam/actual/testpulse/<run>/<device>/raw/
```

The first step is raw decoding with the strict decoder workflow. `decoder.sh` runs per-FIFO decode jobs in parallel for one device, waits for that device to finish, and then moves to the next device. Decoded ROOT files and decoder summaries are written under:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/
```

Decode the full test-pulse runs once:

```bash
sipm4eic-testbeam2026-process/process/scripts/decoder.sh \
    --run 20260618-183625 \
    --run-type testpulse

sipm4eic-testbeam2026-process/process/scripts/decoder.sh \
    --run 20260618-185127 \
    --run-type testpulse
```

The later `checker.sh` and `dcalib.sh` commands select the device/FIFO subsets relevant to each calibration campaign. Keeping decoding run-wide avoids having a partially decoded run directory that can confuse downstream workflows.

The default decoder policy is strict:

```text
--allowed-spill-errors 0
```

If a completed spill exceeds this threshold, the decoder suppresses the entire spill: it writes neither START_SPILL, payload, nor END_SPILL for that spill. Spill markers that are written use `fine = 0`. The decoder writes a spill only after the matching END_SPILL has been found. If EOF arrives while a spill is still open, the incomplete spill and its buffered payload are discarded, so the decoded ROOT file remains spill-balanced. If the DAQ readout already discarded a spill payload and wrote `0xdeadbeef` markers, the decoder counts those DAQ-suppressed records and skips them; it never writes fake spill markers for them. Each decoded file also has a sidecar `.summary` file with counts of spills written, completed spills suppressed by errors, incomplete spills discarded, DAQ-suppressed records, and decoding errors.

After decoding, run the independent checker workflow on the decoded files. This is a read-only preflight step; it does not modify the data and it is not part of `dcalib.sh`. It writes one ASCII `.check` report next to each decoded FIFO file and verifies basic stream consistency, including spill-marker counts and START/END spill-counter pairing.

Run the checker over the full decoded runs:

```bash
sipm4eic-testbeam2026-process/process/scripts/checker.sh \
    --run 20260618-183625 \
    --run-type testpulse

sipm4eic-testbeam2026-process/process/scripts/checker.sh \
    --run 20260618-185127 \
    --run-type testpulse
```

The calibration-specific device/FIFO selections are applied later by `dcalib.sh`. Running the checker over the complete decoded run gives one run-level consistency report for everything decoded, including devices and FIFOs that are not used for a particular TDC-calibration subset.

The checker also writes run-level FIFO selection lists:

```text
/data/2026-testbeam/process/<run>/<run>.good-fifos.list
/data/2026-testbeam/process/<run>/<run>.bad-fifos.list
```

The good/bad FIFO lists are diagnostics. They are useful to understand which FIFOs are compatible with the common spill structure, but the normal downstream scripts do not consume them automatically. If the checker identifies bad FIFOs, use the existing `--devices` and `--fifos` filters to select the hardware subset you want to process.

Each command follows the same hierarchy as the processing workflow: per-FIFO `.check` files are produced first, then one device-level `.check` file is written for each selected device, and finally one run-level `.check` file is written for the selected data set. The most important fields in each `.check` file are:

```text
start_spill_type7
end_spill_type15
unknown_words
spill_counter_consistent
open_spill_at_eof
spill_count_balance
errors
```

The calibration jobs should be launched only after the relevant `.check` files are understood. `dcalib.sh` reads decoded inputs from:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_*.root
```

and writes all per-FIFO calibration products and temporary files to:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib/
```

For each FIFO input file, the `dcalib.sh` workflow itself remains:

```text
clean -> sort -> dcalib
```

The final calibration snippet is named, for example:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib/dcalib.fifo_0.conf
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
/data/2026-testbeam/process/<run>/<device>/dcalib/dcalib.fifo_0.conf
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
    --input /data/2026-testbeam/process/20260618-183625/kc705-200/dcalib/dcalib.fifo_{0..7}.conf \
    --output /data/2026-testbeam/process/tdc.timing.conf
```

### Merge Cherenkov TDC Calibration: Chips 0..3

This merges the RDO Cherenkov calibration for chips `0..3`, using run `20260618-183625`:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/20260618-183625/rdo-{192..199}/dcalib/dcalib.fifo_{0..15}.conf \
    --output /data/2026-testbeam/process/tdc.cherenkov.chips_0_3.conf
```

### Merge Cherenkov TDC Calibration: Chips 4..7

This merges the RDO Cherenkov calibration for chips `4..7`, using run `20260618-185127`:

```bash
sipm4eic-testbeam2026-process/process/scripts/merge_calibration_file.sh \
    --input /data/2026-testbeam/process/20260618-185127/rdo-{192..199}/dcalib/dcalib.fifo_{16..31}.conf \
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
    --report sipm4eic-testbeam2026-process/process/config/calibration/tdc.20260618.md \
    --expect 200:0..7 \
    --expect 192..199:0..31
```

The completed output keeps all measured concrete rows and appends a separate `[TDC]` section containing a single fallback row:

```text
* * * * * 0 <global_iif>
```

The wildcard row has lower specificity than every measured concrete calibration row, so it is used only when a concrete TDC calibration is absent. The `--expect` options encode the actual calibration campaigns from the merge commands: TIMING device `200` FIFOs `0..7`, and Cherenkov RDO devices `192..199` FIFOs `0..31`. The generated report records those expectations explicitly and documents missing rows relative to them, coverage by device/chip, leave-one-out pattern comparisons, and the phase error expected from the fallback.

The final `process/config/calibration/tdc.20260618.conf` still contains only the `[TDC]` calibration section. To run `calibrator`, this TDC file must be combined with `[CHANNEL]` and `[TRIGGER]` sections, or with explicit wildcard defaults for those sections.

## Checking The Calibration With The Laser/Test-Pulse Runs

A useful closure check is to run the normal processing chain on the same two runs used to derive the TDC calibration. Before running these commands, the relevant raw files must already have been decoded with `decoder.sh`; `process.sh` reads from `/data/2026-testbeam/process/<run>/<device>/decoded/` and writes device-local products to `/data/2026-testbeam/process/<run>/<device>/process/`. It now stops after producing merged per-spill ROOT files. Run `trigger.sh` afterwards to produce triggered frames from those merged spill files.

These runs need one important qualification: they do not provide reliable timing synchronisation across different devices.

During these test-pulse runs, the hardware SPILL input of each device was used to send test pulses to ALCOR channels configured in `opMode = 2`. The data-taking spill was therefore generated in software. As a result, `kc705-200` and the individual RDO devices are not guaranteed to share a common synchronous time reference in these runs.

The calibration check is therefore a special per-RDO diagnostic workflow, not the normal run-level trigger workflow:

```text
process one rdo-N subset
trigger on one test-pulse channel in that same rdo-N
inspect only channels from the same rdo-N
ignore relative timing to other devices for this check
```

This is still a useful closure check. Within one RDO, the channels are read by the same electronics stream and are synchronous to the selected test-pulse trigger. Normal physics or run-level production should instead run `process.sh` without `--devices` and then run `trigger.sh` without `--devices`, producing run-level `triggered.<tag>.root` files.

The completed TDC file described above contains only `[TDC]` rows. For this check, build a complete calibrator configuration by adding explicit zero defaults for channel and trigger offsets:

```bash
cp sipm4eic-testbeam2026-process/process/config/calibration/tdc.20260618.conf \
   /data/2026-testbeam/process/calibration.20260618.check.conf

cat >> /data/2026-testbeam/process/calibration.20260618.check.conf <<'EOF'

[CHANNEL]
# device fifo column pixel offset
* * * * 0.0

[TRIGGER]
# device fifo offset
* * 0.0
EOF
```

This keeps the check focused on the TDC calibration only. If channel or trigger offsets are later measured independently, they can replace these wildcard zero defaults.

### Run `20260618-183625`

The ALCOR configuration for this run was:

```text
chips 0,1,2,3: opMode = 2
chips 4,5,6,7: opMode = 1
```

Therefore FIFOs `0..15` receive direct test pulses, while FIFOs `16..31` should record LASER light triggered at each test pulse. For each RDO, trigger on a test-pulse channel in that same RDO:

```text
device = RDO number
fifo   = 0
column = 0
pixel  = 0
```

Dedicated trigger files are provided for `rdo-192` through `rdo-199`:

```text
process/config/trigger/calib_check_20260618-183625_rdo-192.conf
...
process/config/trigger/calib_check_20260618-183625_rdo-199.conf
```

For this special same-RDO closure check, run the calibrated processing/merge step and then the trigger step once per RDO device:

```bash
for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/process.sh \
        --run 20260618-183625 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --calibration /data/2026-testbeam/process/calibration.20260618.check.conf

    sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
        --run 20260618-183625 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --trigger sipm4eic-testbeam2026-process/process/config/trigger/calib_check_20260618-183625_rdo-${dev}.conf calibcheck_rdo-${dev} \
        --window 32
done
```

The output triggered files are expected at:

```text
/data/2026-testbeam/process/20260618-183625/triggered.calibcheck_rdo-192.root
...
/data/2026-testbeam/process/20260618-183625/triggered.calibcheck_rdo-199.root
```

### Run `20260618-185127`

The ALCOR configuration for this run was the opposite:

```text
chips 4,5,6,7: opMode = 2
chips 0,1,2,3: opMode = 1
```

Therefore FIFOs `16..31` receive direct test pulses, while FIFOs `0..15` should record LASER light triggered at each test pulse. For each RDO, trigger on a test-pulse channel in that same RDO:

```text
device = RDO number
fifo   = 16
column = 0
pixel  = 0
```

Dedicated trigger files are provided for `rdo-192` through `rdo-199`:

```text
process/config/trigger/calib_check_20260618-185127_rdo-192.conf
...
process/config/trigger/calib_check_20260618-185127_rdo-199.conf
```

For this special same-RDO closure check, run the calibrated processing/merge step and then the trigger step once per RDO device:

```bash
for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/process.sh \
        --run 20260618-185127 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --calibration /data/2026-testbeam/process/calibration.20260618.check.conf

    sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
        --run 20260618-185127 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --trigger sipm4eic-testbeam2026-process/process/config/trigger/calib_check_20260618-185127_rdo-${dev}.conf calibcheck_rdo-${dev} \
        --window 32
done
```

The output triggered files are expected at:

```text
/data/2026-testbeam/process/20260618-185127/triggered.calibcheck_rdo-192.root
...
/data/2026-testbeam/process/20260618-185127/triggered.calibcheck_rdo-199.root
```

### Analysis With `deltat.C`

The triggered files can be checked with the example macro:

```text
macros/example/deltat.C
```

The macro reads the triggered-frame output, finds the requested trigger hit in each frame, and fills the 2D histogram `hDeltaT` with:

```text
delta_t = hit.time - trigger.time
```

Here `time` is the calibrated time stored by the processing chain in the triggered-frame output. The histogram x axis is the global channel index, and the y axis is `delta_t`. The histogram is normalised by the number of trigger hits found.

For run `20260618-183625`, use the same per-RDO test-pulse trigger channel as the processing step:

```bash
for dev in {192..199}; do
    root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/20260618-183625/triggered.calibcheck_rdo-${dev}.root\", 1, ${dev}, 0, 0, 0, \"/data/2026-testbeam/process/20260618-183625/deltat.calibcheck_rdo-${dev}.root\")"
done
```

For run `20260618-185127`, use the FIFO-16 test-pulse trigger channel:

```bash
for dev in {192..199}; do
    root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/20260618-185127/triggered.calibcheck_rdo-${dev}.root\", 1, ${dev}, 16, 0, 0, \"/data/2026-testbeam/process/20260618-185127/deltat.calibcheck_rdo-${dev}.root\")"
done
```

The useful closure check is whether channels in the same RDO produce narrow, stable `delta_t` structures relative to that RDO's direct test-pulse trigger. For `20260618-183625`, the LASER-side channels are expected on FIFOs `16..31`; for `20260618-185127`, they are expected on FIFOs `0..15`. Do not use these runs to judge timing offsets between different devices.

### What To Inspect

These validation runs are not used to derive new TDC constants. They are a closure check of the calibrated processing chain:

```text
decode -> calibrate -> sort -> merge -> trigger/frame -> deltat analysis
```

The trigger is a direct test-pulse channel in the same RDO being checked. The interesting check is the time distribution of same-RDO channels in the triggered frames. A narrow distribution around a stable relative time indicates that the calibrated `time` branch is being propagated correctly through calibration, sorting, merging, frame building, and analysis.

## Notes And Checks

- The commands above assume the executables have already been built and installed with CMake.
- `dcalib.sh` currently hardcodes `MIN_PAIRS=1000`.
- `dcalib.sh` currently hardcodes `SORT_WINDOW=32768` for the sorting step before TDC calibration.
- The TDC calibration should be derived from otherwise uncalibrated data. The sorter stage may compute nominal hit time internally for ordering, but it does not create or persist a `time` branch when the input does not already have one.
- The output `.conf` fragments are per-FIFO `[TDC]` snippets. A complete calibration file still needs `[CHANNEL]` and `[TRIGGER]` sections, or suitable defaults, before it can be used by `calibrator`.
