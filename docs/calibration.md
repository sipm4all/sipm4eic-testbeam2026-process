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
/data/2026-testbeam/process/<run>/check/<run>.good-fifos.list
/data/2026-testbeam/process/<run>/check/<run>.bad-fifos.list
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

For a normal run-level production workflow, do not use `--devices`: run `process.sh` once for the run, then run `trigger.sh` once for the run-level merged spills. The `--devices` mode below is only for the same-RDO calibration closure check, where the test-pulse trigger and the channels being inspected must belong to the same RDO.

Before `trigger.sh` was split out, `process.sh --devices rdo-X --trigger ...` did the full chain for one RDO: processing, same-RDO spill merge, trigger extraction, and final `hadd`. With the current split workflow, keep the same per-RDO diagnostic logic but run it in two stages: first process all RDOs, then trigger all RDOs.

```bash
for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/process.sh \
        --run 20260618-183625 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --calibration /data/2026-testbeam/process/calibration.20260618.check.conf
done

for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
        --run 20260618-183625 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --trigger sipm4eic-testbeam2026-process/process/config/trigger/calib_check_20260618-183625_rdo-${dev}.conf calibcheck_rdo-${dev} \
        --window 32
done
```

This produces one diagnostic file per RDO:

```text
/data/2026-testbeam/process/20260618-183625/trigger/triggered.calibcheck_rdo-192.root
...
/data/2026-testbeam/process/20260618-183625/trigger/triggered.calibcheck_rdo-199.root
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

Use the same two-stage per-RDO diagnostic workflow:

```bash
for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/process.sh \
        --run 20260618-185127 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --calibration /data/2026-testbeam/process/calibration.20260618.check.conf
done

for dev in {192..199}; do
    sipm4eic-testbeam2026-process/process/scripts/trigger.sh \
        --run 20260618-185127 \
        --run-type testpulse \
        --devices rdo-${dev} \
        --trigger sipm4eic-testbeam2026-process/process/config/trigger/calib_check_20260618-185127_rdo-${dev}.conf calibcheck_rdo-${dev} \
        --window 32
done
```

This produces one diagnostic file per RDO:

```text
/data/2026-testbeam/process/20260618-185127/trigger/triggered.calibcheck_rdo-192.root
...
/data/2026-testbeam/process/20260618-185127/trigger/triggered.calibcheck_rdo-199.root
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

For this calibration closure check, keep using the stored test-pulse trigger hit as the reference. This checks the calibrated hit times propagated by:

```text
decode -> calibrate -> sort -> merge -> trigger/frame
```

without adding the trained TIMING estimator as another processing stage.

For run `20260618-183625`, use the same test-pulse trigger channel as the diagnostic trigger. Example loop:

```bash
for dev in {192..199}; do
    root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/20260618-183625/trigger/triggered.calibcheck_rdo-${dev}.root\", 1, ${dev}, 0, 0, 0, \"/data/2026-testbeam/process/20260618-183625/trigger/deltat.calibcheck_rdo-${dev}.root\")"
done
```

For run `20260618-185127`, use the FIFO-16 test-pulse trigger channel:

```bash
for dev in {192..199}; do
    root -l -b -q "sipm4eic-testbeam2026-process/macros/example/deltat.C(\"/data/2026-testbeam/process/20260618-185127/trigger/triggered.calibcheck_rdo-${dev}.root\", 1, ${dev}, 16, 0, 0, \"/data/2026-testbeam/process/20260618-185127/trigger/deltat.calibcheck_rdo-${dev}.root\")"
done
```

The useful closure check is whether channels in the same RDO produce narrow, stable `delta_t` structures relative to that RDO's direct test-pulse trigger. For `20260618-183625`, the LASER-side channels are expected on FIFOs `16..31`; for `20260618-185127`, they are expected on FIFOs `0..15`. Do not use these runs to judge timing offsets between different devices.

### Deriving Channel Offsets With `fit_calib.C`

The `fit_calib.C` macro derives per-channel `[CHANNEL]` timing offsets from the
`hDeltaT` histogram produced by `deltat.C`. It is intentionally inclusive: the
statistics decision is made from the fraction

```text
integraler = histogram error / histogram integral
```

inside the `[-2,2]` interval. The current tiers are:

```text
integraler > 0.30       skip the channel
0.15 .. 0.30            use the restricted mean
0.08 .. 0.15            rebin by 4, then fit the peak
0.04 .. 0.08            rebin by 2, then fit the peak
<= 0.04                 fit at the original binning
```

The thresholds correspond approximately to increasing effective statistics.
When a fit fails, the restricted mean is used as a fallback. A skipped channel
is not shifted in the aligned diagnostic histogram and is omitted from the
generated calibration file.

Run it from the repository checkout with:

```bash
root -l -b -q \
  -e '.L macros/example/fit_calib.C' \
  -e 'fit_calib("deltat.ring.root", "fit_calib.ring.root", "timing_offsets.ring.conf")'
```

The third argument is optional. If it is omitted, the macro replaces the final
`.root` suffix of the output name with `.conf`. The ROOT output contains:

```text
hCalib             fitted/mean channel offsets
hMethod            method used for each channel (0 skipped, 1 mean,
                   2 rebin-4 fit, 3 rebin-2 fit, 4 full fit)
hDeltaT_aligned    hDeltaT after applying accepted channel offsets
```

The generated text file has this format:

```text
[CHANNEL]
# device fifo column pixel offset
```

Only accepted channels are written. The macro converts the global Cherenkov
channel index back to the electronics address using 2048 channels arranged as
eight devices, eight 32-channel chips per device, and four pixels per column.
The calibration convention is:

```text
calibrated_time = raw_time - offset
```

Inspect `hDeltaT_aligned` and the fit uncertainties before using the generated
`.conf` in a production `calibrator` run. The macro does not invent values for
channels with insufficient information; such channels remain available for a
later hardware-pattern or wildcard fallback.

### Suspected Lane Clock Instabilities

The aligned `hDeltaT` data were also inspected after grouping the global
Cherenkov channels in blocks of eight:

```text
lane_number = global_channel / 8
chip_number = lane_number / 4
lane_in_chip = lane_number % 4
```

Each lane corresponds to one electronics FIFO in the analysis channel mapping.
Several lanes show a second lane-wide population displaced by approximately one
native time unit from the main peak. The suspected lanes are:

```text
lane   chip   lane-in-chip   displacement
----   ----   ------------   ------------
90     22     2              +1
102    25     2              +1
117    29     1              -1
120    30     0              +1
137    34     1              +1
142    35     2              -1
173    43     1              -1
178    44     2              -1
191    47     3              -1
```

There are nine suspected lanes in total. They all belong to different chips:
no chip has more than one suspicious lane. These are suspected clock
instabilities or one-clock phase ambiguities, not final calibration corrections.
The effect is common to the eight channels in a lane and should therefore not
be treated as eight unrelated channel offsets without further validation.

The initial investigation used the following physics runs, but these runs are
not homogeneous. A lane compared against the timing-estimator reference
`T` does not necessarily produce the same `delta_t` distribution in all three
groups. The beam conditions were different:

```text
Group A: 11 GeV negative beam, collimators at 24, focus at 20 m

20260623-185238
20260623-191413
20260623-200510
20260623-202754
20260623-222714
20260623-224255

Group B: 11 GeV positive beam, collimators at 10, focus at 4 m

20260623-204903
20260623-210048
20260623-211302

Group C: 11 GeV negative beam, collimators at 24, focus at 4 m

20260623-212754
20260623-214144
20260623-221530
```

The reason for the observed run-to-run difference is currently unknown and is
under investigation. Until it is understood, lane stability studies and any
further calibration comparison should be performed separately within Groups
A, B, and C. A calibration or offset derived from one group must not be
assumed to describe the other groups.

For completeness, the complete run list used in the initial study was:

```text
20260623-185238
20260623-191413
20260623-200510
20260623-202754
20260623-204903
20260623-210048
20260623-211302
20260623-212754
20260623-214144
20260623-221530
20260623-222714
20260623-224255
```

For a focused check of lane 173, use the field selector for device 197,
FIFO 13, all columns and pixels. Global channels 1384 through 1391 map to
this FIFO. Use the trained timing estimator result `T` as the reference:

```bash
root -l -b -q \
  -e '.L macros/example/deltat.C' \
  -e 'deltat("ring.root", \
             field_selector_t(1,197,13,-1,-1), \
             timing_reference_t("T"), \
             "deltat_fifo173.root")'
```

This produces `hDeltaT` with:

```text
delta_t = time_fifo173 - T
```

The output can be inspected together with the FIFO-specific histograms
`hDeltaT_fifo090`, `hDeltaT_fifo102`, and so on, produced by `fit_calib.C`.

### Spill-Dependent Clock Jump

The spill-resolved histogram

```text
20260623-185238.deltat.fifo173.root:hDeltaT_spill
```

shows an additional effect in lane 173. The lane corresponds to device
197, FIFO 13. Restricting the per-spill calculation to the central
`delta_t` interval `[-2,2]`, its distribution relative to the timing-estimator
reference `T` changes from spill to spill: one state has a mean around
`+0.25` native time units, while the other is around `-0.73`. The separation is
close to one native clock unit.

This is therefore not simply a fixed FIFO offset. It is a spill-dependent clock
state or phase jump which occurs at spill transitions. In the inspected file,
the spill means alternate between these two broad levels; the exact mean varies
because of the finite statistics and the intrinsic distribution width. The
observation is important because a single channel or FIFO calibration constant
cannot remove this effect from all spills.

Until the mechanism is understood, calibration checks should retain the spill
dimension and should not collapse all spills into one `hDeltaT` distribution.
Potential future corrections will need either a spill-state indicator or a
spill-local timing adjustment.

The temporary diagnostic macro
`macros/example/clock_transition.C` detects the two spill states using only the
central `delta_t` interval `[-2,2]`. It performs a two-cluster fit of the
per-spill means, uses the midpoint of the two cluster means as the
classification threshold, keeps the majority state as reference, shifts only
the minority state by one native unit, and writes both a corrected spill
histogram and a run-specific `[CLOCK]` file. For lane 173 in the example run:

```bash
root -l -b -q \
  -e '.L macros/example/clock_transition.C' \
  -e 'clock_transition("20260623-185238.deltat.fifo173.root", \
                       "20260623-185238", 197, 13, \
                       "clock_transition.fifo173.root", \
                       "clock-corrections.fifo173.conf")'
```

The ROOT output contains `hDeltaT_spill` and
`hDeltaT_spill_aligned`. The generated correction file contains one compact
row of the form:

```text
[CLOCK]
# run device fifo correction spill...
20260623-185238 197 13 -1 4 6 9 ...
```

The correction sign is derived from the direction of the minority-state shift:
the calibrator applies `corrected_time = calibrated_time - correction`. Thus a
minority state one unit above the majority state receives `+1`, while a
minority state one unit below it receives `-1`. This macro is a diagnostic
starting point; the spill-state classification should be validated for each
run and lane before using the generated corrections in production processing.

For physics timing studies, after a triggered file has been augmented with:

```bash
sipm4eic-testbeam2026-process/process/scripts/timing.sh \
    --run <run> \
    --trigger fingers
```

which is equivalent to running:

```bash
process/bin/timing \
    --input triggered.fingers.root \
    --output triggered.fingers.timing.root
```

the same macro can use the trained TIMING estimator as the reference:

```cpp
deltat("triggered.fingers.timing.root",
       channel_selector_t(1, -1),
       timing_reference_t("T"),
       "deltat.fingers.timing.root");
```

Valid timing-reference names are `"T"`, `"T0"`, and `"T1"`. This mode requires the estimator branches `timing_valid`, `T0`, `T1`, and `T`, and uses only frames with `timing_valid == 1`.

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

### Run-Specific Clock Corrections

Clock instabilities are kept in a separate file from the stable TDC/channel
calibration. The optional file is passed to `calibrator` with `--clock`, and
the run identifier must be supplied with `--run`:

```bash
process/bin/calibrator \
    --input decoded.root \
    --output calibrated.root \
    --config process/config/calibration/20260819.calib.conf \
    --clock clock-corrections.conf \
    --run 20260623-185238
```

The format stores one fixed correction sign and an explicit list of affected
spills for each `(run, device, fifo)`:

```text
[CLOCK]
# run device fifo correction spill...
20260623-185238 197 13 +1 4 7 12 18
20260623-185238 194 26 +1 3 9 15
```

Clock-transition correction files use a fixed `+1` convention. The lower-time
state is the canonical state and only higher-time spills are listed for
correction. A missing row means no correction; unaffected spills are not
listed. The calibrator applies:

```text
corrected_time = calibrated_time - clock_correction
```

to ALCOR hits from the specified `(device,fifo)` while processing the listed
spill. A single clock file may contain rows for multiple runs; `--run` selects
which rows are active and rows for other runs are ignored. If the selected run
has no rows, no clock correction is applied. Duplicate `(run,device,fifo)` rows,
duplicate spills, and malformed spill numbers are rejected. Omitting `--clock`
preserves the normal calibration behavior.
