# Workflow Scripts

This directory contains shell scripts for running larger processing workflows. The scripts assume they are run from the repository checkout and that executables have been built in the `process/bin/` directory.

Build first:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

## dcalib.sh

`dcalib.sh` runs TDC calibration over all decoded FIFO files in a run. For each FIFO file the workflow is:

```text
clean -> sort -> dcalib
```

It uses the same input and output base directories as `process.sh`:

```bash
ipath="/data/2026-testbeam/actual/testpulse"
opath="/data/2026-testbeam/process"
```

Run:

```bash
process/scripts/dcalib.sh --run RUN_NAME --period 10000
process/scripts/dcalib.sh --run RUN_NAME --period 10000 --devices kc705-200 rdo-{192..195} --fifos {0..16}
```

Required command-line options:

```text
--run RUN          run name/directory
--period PERIOD    expected pulser period in coarse clock cycles
```

Optional filters:

```text
--devices all                    process all devices, default
--devices kc705-200 rdo-192      process selected device directories
--devices kc705-200 rdo-{192..195}
                                 process one device plus an inclusive RDO range
--devices kc705-200 "rdo-{192..195}"
                                 quoted prefixed ranges are also accepted

--fifos all                process all FIFOs, default
--fifos 0 4 8              process only selected FIFOs
--fifos {0..16}            process an inclusive FIFO range
--fifos "{0..16}"          quoted inclusive ranges are also accepted
```

The shell expands unquoted brace ranges before the script receives them. Quoted brace ranges are expanded by the script itself. Device filters match the actual device directory basename, such as `kc705-200` or `rdo-192`; FIFO filters remain numeric.

For each input file like:

```text
/data/2026-testbeam/actual/testpulse/<run>/<device>/decoded/alcdaq.fifo_0.root
```

the script writes outputs in the corresponding device directory:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/dcalib.fifo_0.conf
```

The final `.root` file contains compact TDC diagnostic histograms. The `.conf` file is the per-FIFO TDC calibration snippet; it contains only a `[TDC]` section in the format consumed by `calibrator`. These calibration outputs are not deleted by the script.

During each job, the script also creates intermediate files in the same device output directory:

```text
/data/2026-testbeam/process/<run>/<device>/cleaned.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/sorted.cleaned.fifo_0.root
```

`cleaned.*.root` is produced by `cleaner`; `sorted.cleaned.*.root` is produced by `sorter` and is the input to `dcalib`. These intermediate files are removed after `dcalib` succeeds.

`dcalib.sh` currently hardcodes:

```bash
SORT_WINDOW=32768
MIN_PAIRS=1000
```

`SORT_WINDOW` is passed to `sorter --window`. Channels with fewer than `MIN_PAIRS` adjacent-hit pairs are skipped by `dcalib`.


## merge_calibration_file.sh

`merge_calibration_file.sh` merges calibration text fragments into one calibration file:

```bash
process/scripts/merge_calibration_file.sh \
  --input dcalib.fifo_0.conf dcalib.fifo_1.conf \
  --output calibration.merged.conf
```

`--input` may be repeated, and each `--input` may be followed by more than one file. This supports shell-expanded patterns such as `--input rdo-{192..199}/dcalib.fifo_{0..15}.conf`. The script understands the standard calibration sections:

```text
[TDC]
[CHANNEL]
[TRIGGER]
```

Repeated section headers from input fragments are collapsed, so many per-FIFO `dcalib.fifo_*.conf` files, each containing a `[TDC]` section, become one file with a single `[TDC]` section. Comments and blank lines from input fragments are not copied; the merged file gets fresh section headers and column comments.

The script does not resolve duplicate calibration rows. If two fragments contain the same concrete calibration address, `calibrator` will later reject the merged file as ambiguous/duplicate.


## complete_tdc_calibration.py

`complete_tdc_calibration.py` inspects a merged `[TDC]` calibration file, reports missing concrete TDC rows, and writes a calibration file with global wildcard fallback rows:

```bash
process/scripts/complete_tdc_calibration.py \
  --input tdc.cherenkov.conf \
  --output tdc.cherenkov.with_global_fallback.conf \
  --report tdc.cherenkov.global_fallback.md
```

The expected geometry is derived from the devices and FIFOs found in the input unless `--devices` or `--fifos` are specified explicitly. Valid columns for each FIFO follow the DAQ mapping:

```text
column = 2 * (fifo % 4)
column = 2 * (fifo % 4) + 1
```

The script does not fill missing concrete rows. Instead it appends a separate `[TDC]` section containing one low-specificity wildcard row per TDC:

```text
* * * * tdc off iif
```

The fallback values are global medians computed from measured rows with the same TDC index. Concrete measured rows remain more specific and therefore take precedence in `calibrator`.

## process.sh

`process.sh` is the full processing workflow. It currently performs:

```text
calibrate each decoded FIFO file
sort each calibrated FIFO file
after-pulse suppress each sorted calibrated FIFO file
merge lanes per device with --split-spills
merge matching spills across devices
run trigger configurations per spill
hadd triggered spill files per trigger tag
```

Per-FIFO intermediate data files use stage-preserving prefixes:

```text
calibrated.fifo_0.root
sorted.calibrated.fifo_0.root
aps.sorted.calibrated.fifo_0.root
```

Run:

```bash
process/scripts/process.sh \
  --run RUN_NAME \
  --calibration process/config/calibration/calibration_example.conf \
  --trigger process/config/trigger/trigger_range.conf range \
  --trigger process/config/trigger/trigger_set.conf set \
  --window 256
```

Required command-line options:

```text
--run RUN                  run name/directory
--calibration FILE         timing calibration configuration
--trigger FILE TAG         trigger configuration and output tag; may be repeated
```

Optional command-line options:

```text
--run-type TYPE            input run type under /data/2026-testbeam/actual, default physics
                           supported values: physics, testpulse
--devices DEVICE ...       device directory names to process, default all
                           accepts names such as kc705-200, rdo-192, rdo-{192..199}
--overwrite                overwrite existing workflow outputs instead of skipping them
--window VALUE             trigger frame window, default 256
--help                     print usage
```

Short aliases are also accepted for the older single-argument options:

```text
-r, -c, -t, -w, -h
```

`--trigger` may be given multiple times. Tags must be unique because they are used in output filenames such as:

```text
triggered.<tag>.spill_0000.root
triggered.<tag>.root
```

Input files are read from:

```text
/data/2026-testbeam/actual/<run-type>/<run>/...
```

Use the default `--run-type physics` for normal physics data and `--run-type testpulse` for calibration/test-pulse runs. Use `--devices` when a workflow should process only selected device directories, for example a same-RDO calibration closure check. Final spill merging uses only the devices processed by the current invocation, so stale split-spill files from earlier per-device runs in the same output directory are ignored.

By default, `process.sh` does not overwrite existing workflow outputs. If an output from a previous stopped or failed processing attempt is found, the corresponding stage is skipped and the existing file is reused. Stages that do run successfully still perform the normal cleanup of their inputs/intermediate files. Pass `--overwrite` to force regeneration of existing outputs.

Important variables still configured near the top of the script:

```bash
actual_base            input decoded-data base directory before run type
opath                  output processing base directory
WRITE_LOGS             0 prints to terminal, 1 writes log files
CLEAN_DEVICE_SPILLS    remove intermediate device spill files
CLEAN_MERGED_SPILLS    remove merged spill files after triggering
CLEAN_TRIGGERED_SPILLS remove per-spill trigger output after hadd
```

When `WRITE_LOGS=0` and `CLEAN_DEVICE_SPILLS=1`, empty per-device output directories are removed after the device-level split-spill files have been consumed by the final merge. If logs are enabled, directories are kept so their log files remain available.

## trigger.sh

`trigger.sh` is a smaller helper that runs one trigger configuration over already merged split-spill files:

```bash
process/scripts/trigger.sh RUN_NAME
```

It uses `CONFIG` and `WINDOW` variables near the top of the script. For multi-configuration production processing, prefer `process.sh`.
