# Workflow Scripts

This directory contains shell scripts for running larger processing workflows. Raw decoding is handled by `decoder.sh` and `process/bin/decoder`. The scripts assume they are run from the repository checkout and that executables have been built in the `process/bin/` directory.

Build first:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```



## decoder.sh

`decoder.sh` decodes raw per-FIFO `.dat` files into ROOT `alcor` files. It follows the same hierarchy as the main processing workflow:

```text
FIFO decode jobs in one device
  -> wait for that device
  -> move to the next device
  -> finish the run
```

Run:

```bash
process/scripts/decoder.sh --run RUN_NAME
process/scripts/decoder.sh --run RUN_NAME --run-type testpulse --devices kc705-200 rdo-{192..195} --fifos {0..16}
process/scripts/decoder.sh --run RUN_NAME --allowed-spill-errors 1
```

Options:

```text
--run RUN                    run name/directory
--run-type TYPE              input run type, default physics; supported: physics, testpulse
--devices all                process all devices, default
--devices DEV ...            process selected device directories
--fifos all                  process all FIFOs, default
--fifos FIFO ...             process selected FIFO numbers or inclusive ranges
--allowed-spill-errors N     maximum errors before a spill payload is emptied, default 0
--overwrite                  overwrite existing decoded ROOT files
```

For an input file:

```text
/data/2026-testbeam/actual/<run-type>/<run>/<device>/raw/alcdaq.fifo_13.dat
```

it writes:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_13.root
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_13.summary
```

The `.summary` file is produced by `decoder` and records spills found/written/emptied, DAQ-readout-suppressed `0xdeadbeef` records, and raw decoding error counters. Normal spill markers use `fine = 0`; decoder-suppressed spills use `fine = 1`. `0xdeadbeef` records are counted and skipped, never written as fake spill markers.

## checker.sh

`checker.sh` runs the lightweight data sanity checker over decoded per-FIFO ROOT files produced by `decoder.sh`. It does not modify data. For each input file it writes one ASCII `.check` report next to the decoded ROOT file.

Run:

```bash
process/scripts/checker.sh --run RUN_NAME
process/scripts/checker.sh --run RUN_NAME --run-type testpulse --devices kc705-200 rdo-{192..195} --fifos {0..16}
```

The accepted filters mirror the `dcalib.sh` style:

```text
--run RUN            run name/directory
--run-type TYPE      accepted for compatibility with decoder.sh; decoded input is read from /data/2026-testbeam/process/<run>
--devices all        process all devices, default
--devices DEV ...    process selected device directories
--fifos all          process all FIFOs, default
--fifos FIFO ...     process selected FIFO numbers or inclusive ranges
```

For an input file:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_0.root
```

it writes per-FIFO reports such as:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_0.check
```

The workflow is hierarchical: per-FIFO checks are produced first, then one device-level check is written for that device, and finally one run-level check is written after all selected devices are complete. Aggregate reports are:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/<device>.check
/data/2026-testbeam/process/<run>/<run>.check
```

The device-level report sums entries and data-word counters from the selected FIFO reports. The logical spill counters remain per-stream quantities: `start_spill_type7` and `end_spill_type15` report the aggregate spill count, not the sum of identical spill markers across FIFOs. The run-level report follows the same convention when combining device reports. If inputs disagree, the aggregate spill-count fields report the most common input count while the min/max fields show the spread and `consistent: no` marks the aggregate as invalid.

Each `.check` file contains:

```text
entries
start_spill_type7
end_spill_type15
alcor_hits_type1
trigger_tags_type9
unknown_words
spill_counter_consistent
open_spill_at_eof
spill_count_balance
spill_count_uniform_start
spill_count_uniform_end
min_start_spill_type7
max_start_spill_type7
min_end_spill_type15
max_end_spill_type15
errors
```

Aggregate device/run `.check` files also include diagnostic pointers and explanations for failed inputs:

```text
problem_check: <path-to-fifo-or-device-check>
error: <path-to-fifo-or-device-check>: START_SPILL count range 0..0 differs from common count 100
```

Use `problem_check` to find the detailed per-FIFO or per-device report when an aggregate has `consistent: no`, nonzero `errors`, or a non-uniform spill count. The uniform spill-count checks compare the spill count per input stream. For per-FIFO reports this is the FIFO START/END count. For device aggregates this is the per-FIFO min/max stored in the device report, so a run-level check does not confuse an 8-FIFO device with a 32-FIFO device. `problem_check` entries are inputs that are internally inconsistent or deviate from the most common spill count; the following `error:` lines explain the aggregate-level reason.

Aggregation now treats missing or malformed numeric fields in input `.check` files as a workflow error. The script fails rather than silently converting such values to zero.

The DAQ end-of-spill word used by the current codebase is `type == 15`.

## dcalib.sh

`dcalib.sh` runs TDC calibration over all decoded FIFO files in a run. For each FIFO file the workflow is:

```text
clean -> sort -> dcalib
```

It reads the decoded ROOT files written by `decoder.sh` using the script-level `ipath` base directory, and writes calibration products under each device's `dcalib/` subdirectory using `opath`.

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
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_0.root
```

the script writes outputs in the corresponding device calibration directory:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib/dcalib.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/dcalib/dcalib.fifo_0.conf
```

The final `.root` file contains compact TDC diagnostic histograms. The `.conf` file is the per-FIFO TDC calibration snippet; it contains only a `[TDC]` section in the format consumed by `calibrator`. These calibration outputs are not deleted by the script.

During each job, the script also creates intermediate files in the same `dcalib/` directory:

```text
/data/2026-testbeam/process/<run>/<device>/dcalib/cleaned.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/dcalib/sorted.cleaned.fifo_0.root
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

`--input` may be repeated, and each `--input` may be followed by more than one file. This supports shell-expanded patterns such as `--input rdo-{192..199}/dcalib/dcalib.fifo_{0..15}.conf`. The script understands the standard calibration sections:

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
/data/2026-testbeam/process/<run>/<device>/process/calibrated.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/process/sorted.calibrated.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/process/aps.sorted.calibrated.fifo_0.root
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
--run-type TYPE            accepted for compatibility with decoder.sh; decoded input is read from /data/2026-testbeam/process/<run>
                           default: physics
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

Input files are read from decoded output produced by `decoder.sh`:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_*.root
```

Run `decoder.sh` with the appropriate `--run-type` first. `process.sh` keeps accepting `--run-type` so command lines can mirror the decoder command, but it no longer reads directly from `/data/2026-testbeam/actual`. Use `--devices` when a workflow should process only selected device directories, for example a same-RDO calibration closure check. Final spill merging uses only the devices processed by the current invocation, so stale split-spill files from earlier per-device runs in the same output directory are ignored.

By default, `process.sh` does not overwrite existing workflow outputs. If an output from a previous stopped or failed processing attempt is found, the corresponding stage is skipped and the existing file is reused. Stages that do run successfully still perform the normal cleanup of their inputs/intermediate files. Pass `--overwrite` to force regeneration of existing outputs.

Important variables still configured near the top of the script:

```bash
ipath                  decoded input base directory; input is read from ipath/RUN/DEVICE/decoded
opath                  processing output base directory
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
