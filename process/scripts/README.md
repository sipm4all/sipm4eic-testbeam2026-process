# Workflow Scripts

This directory contains shell scripts for running larger processing workflows. Raw decoding is handled by `decoder.sh` and `process/bin/decoder`. The scripts assume they are run from the repository checkout and that executables have been built in the `process/bin/` directory.

Build first:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

## merge_residual_calibration.py

`merge_residual_calibration.py` combines a base calibration with residual
`[CHANNEL]` offsets produced by `macros/example/fit_calib.C`. It adds residual
offsets to matching exact channels and appends residual-only channels as exact
overrides. The TDC and trigger sections come from the base file.

```bash
process/scripts/merge_residual_calibration.py \
  --base process/config/calibration/20260819.calib.conf \
  --residual fit_calib.deltat.20260821.conf \
  --output process/config/calibration/calibration.20260820.conf
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
--allowed-spill-errors N     maximum errors before a completed spill is suppressed, default 0
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

The `.summary` file is produced by `decoder` and records spills found/written, completed spills suppressed by errors, incomplete spills discarded at EOF, DAQ-readout-suppressed `0xdeadbeef` records, and raw decoding error counters. Spill markers written by the decoder use `fine = 0`. Bad completed spills are suppressed entirely; no artificial START/END marker pair is emitted for them. `0xdeadbeef` records are counted and skipped, never written as fake spill markers. If EOF occurs while a spill is open, that incomplete spill is not emitted to the decoded ROOT tree; this prevents downstream workflows from seeing unmatched START_SPILL markers.

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
/data/2026-testbeam/process/<run>/<device>/check/alcdaq.fifo_0.check
```

The workflow is hierarchical: per-FIFO checks are produced first, then one device-level check is written for that device, and finally one run-level check is written after all selected devices are complete. Aggregate reports are:

```text
/data/2026-testbeam/process/<run>/<device>/check/<device>.check
/data/2026-testbeam/process/<run>/check/<run>.check
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
discard_candidate_count: 2
would_pass_without_discard_candidates: yes
problem_check: <path-to-fifo-or-device-check>
discard_candidate: <path-to-fifo-or-device-check>
error: <path-to-fifo-or-device-check>: START_SPILL count range 0..0 differs from common count 100
```

Use `problem_check` to find the detailed per-FIFO or per-device report when an aggregate has `consistent: no`, nonzero `errors`, or a non-uniform spill count. The uniform spill-count checks compare the spill count per input stream. For per-FIFO reports this is the FIFO START/END count. For device aggregates this is the per-FIFO min/max stored in the device report, so a run-level check does not confuse an 8-FIFO device with a 32-FIFO device. `problem_check` entries are inputs that are internally inconsistent or deviate from the most common spill count; the following `error:` lines explain the aggregate-level reason.

`discard_candidate` entries are the subset of problematic inputs whose START/END spill counts differ from the common count. `would_pass_without_discard_candidates: yes` means that removing those candidates would leave a clean, spill-uniform aggregate. It does not delete or skip anything by itself; it only records the candidate list for a later workflow decision.

At run level, `checker.sh` also writes two machine-readable FIFO selection files:

```text
/data/2026-testbeam/process/<run>/check/<run>.good-fifos.list
/data/2026-testbeam/process/<run>/check/<run>.bad-fifos.list
```

The good list has columns:

```text
device fifo decoded_root check_file
```

The bad list has columns:

```text
device fifo decoded_root check_file reason
```

The selection is conservative. A FIFO is listed as good only if its device is clean as-is, or the device would pass after dropping its local discard candidates, and the resulting device spill count belongs to the run-level common spill count. The bad list records excluded FIFOs and a reason such as `discard_candidate`, `device_not_repairable`, or `run_spill_count_outlier`.

The good/bad FIFO lists are diagnostic outputs. Since the decoder now suppresses malformed or incomplete spills instead of emitting bad spill boundaries, downstream workflows normally process the selected decoded files directly and rely on checker summaries to decide whether a run/device should be trusted.

Aggregation now treats missing or malformed numeric fields in input `.check` files as a workflow error. The script fails rather than silently converting such values to zero.

The DAQ end-of-spill word used by the current codebase is `type == 15`.

## dcalib.sh

`dcalib.sh` runs TDC calibration over all decoded FIFO files in a run. For each FIFO file the workflow is:

```text
clean -> sort -> dcalib
```

It reads the decoded ROOT files written by `decoder.sh` using the script-level `ipath` base directory, and writes calibration products under each device's `dcalib/` subdirectory using `opath`. Use `--devices` and `--fifos` when a calibration campaign should process only selected hardware subsets.

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

`process.sh` is the calibrated processing and merge workflow. It performs:

```text
calibrate each decoded FIFO file
assign coordinates to each calibrated FIFO file
sort each coordinated/calibrated FIFO file
after-pulse suppress each sorted coordinated/calibrated FIFO file
merge lanes per device with --split-spills
merge matching spills across devices
```

It intentionally stops after producing merged per-spill ROOT files. Trigger/frame production is handled by `trigger.sh`, so trigger configurations can be rerun without recalibrating, sorting, AP suppressing, or merging again.

Per-FIFO intermediate data files use stage-preserving prefixes:

```text
/data/2026-testbeam/process/<run>/<device>/process/calibrated.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/process/coordinated.calibrated.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/process/sorted.coordinated.calibrated.fifo_0.root
/data/2026-testbeam/process/<run>/<device>/process/aps.sorted.coordinated.calibrated.fifo_0.root
```

Run:

```bash
process/scripts/process.sh \
  --run RUN_NAME \
  --calibration process/config/calibration/calibration_example.conf
```

Required command-line options:

```text
--run RUN                  run name/directory
--calibration FILE         timing calibration configuration
```

Optional command-line options:

```text
--run-type TYPE            accepted for compatibility with decoder.sh; decoded input is read from /data/2026-testbeam/process/<run>
                           default: physics
--devices DEVICE ...       device directory names to process, default all
                           accepts names such as kc705-200, rdo-192, rdo-{192..199}
--overwrite                overwrite existing workflow outputs instead of skipping them
--help                     print usage
```

Short aliases are also accepted for the remaining single-argument options:

```text
-r, -c, -h
```

Input files are read from decoded output produced by `decoder.sh`:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/alcdaq.fifo_*.root
```

Run `decoder.sh` with the appropriate `--run-type` first. `process.sh` keeps accepting `--run-type` so command lines can mirror the decoder command, but it no longer reads directly from `/data/2026-testbeam/actual`. Use `--devices` when a workflow should process only selected device directories, for example a same-RDO calibration closure check. Final spill merging uses only the devices processed by the current invocation, so stale split-spill files from earlier per-device runs in the same output directory are ignored.

The normal final merged output from `process.sh` is run-level data:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0000.root
/data/2026-testbeam/process/<run>/process/aps.sorted.spill_0001.root
...
```

This is the input normally consumed by `trigger.sh`. When `--devices` is used, `process.sh` writes a diagnostic or calibration-check subset with the device selection encoded in the prefix, for example:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.rdo-192.spill_0000.root
```

That subset mode is not the normal production input for run-level triggering.

By default, `process.sh` does not overwrite existing workflow outputs. If run-level merged spill files already exist for the selected device prefix, the whole process workflow is skipped because those files are the products consumed by `trigger.sh`. If an output from a previous stopped or failed processing attempt is found inside a lower-level stage, the corresponding stage is skipped and the existing file is reused. Stages that do run successfully still perform the normal cleanup of their inputs/intermediate files. Pass `--overwrite` to force regeneration of existing outputs.

Important variables still configured near the top of the script:

```bash
ipath               decoded input base directory; input is read from ipath/RUN/DEVICE/decoded
opath               processing output base directory
WRITE_LOGS          0 prints to terminal, 1 writes log files
CLEAN_DEVICE_SPILLS remove intermediate device spill files
CLEAN_MERGED_SPILLS remove final merged spill files; default 0 because trigger.sh consumes them
```

When `WRITE_LOGS=0` and `CLEAN_DEVICE_SPILLS=1`, empty per-device output directories are removed after the device-level split-spill files have been consumed by the final merge. If logs are enabled, directories are kept so their log files remain available.

### Merger Spill Alignment

`merger` aligns input streams by spill counter. At the start of each merge step it selects the lowest current START_SPILL counter among the input streams and merges only the streams that currently contain that counter. Streams whose next spill has a larger counter are left untouched until their counter becomes the selected counter. This supports decoded inputs where a bad or incomplete spill was suppressed entirely in one FIFO. For example:

```text
FIFO-1: spill 0, spill 1, spill 2
FIFO-2: spill 0,          spill 2
```

produces:

```text
merged spill 0: FIFO-1, FIFO-2
merged spill 1: FIFO-1
merged spill 2: FIFO-1, FIFO-2
```

The duplicate START_SPILL and END_SPILL words are still collapsed in the output `alcor` tree. In addition, every merged output file contains a `spill_participation` tree with one entry per merged spill:

```text
spill/I
counter/I
nsources/I
source_device[nsources]/I
source_fifo[nsources]/I
```

Only contributing sources are recorded. If an input file already contains `spill_participation`, the merger propagates those original sources into the new metadata tree. This lets the final cross-device merge preserve the FIFO-level participation recorded by the earlier per-device split merge.

## trigger.sh

`trigger.sh` is the triggered-frame workflow. It consumes the merged split-spill ROOT files produced by `process.sh`, runs one or more trigger configurations over each spill in parallel, and then combines each trigger tag with `hadd`.

Run:

```bash
process/scripts/trigger.sh \
  --run RUN_NAME \
  --trigger process/config/trigger/trigger_range.conf range \
  --trigger process/config/trigger/trigger_set.conf set \
  --window 256
```

Required command-line options:

```text
--run RUN                  run name/directory
--trigger FILE TAG         trigger configuration and output tag; may be repeated
```

Optional command-line options:

```text
--run-type TYPE            accepted for symmetry with process.sh, default physics
--devices DEVICE ...       same device subset used with process.sh, default all
--input-prefix PREFIX      explicit merged-spill prefix, default derived from --devices
--overwrite                overwrite existing trigger outputs instead of skipping them
--window VALUE             trigger frame window, default 256
--clean-triggered-spills   remove triggered.<tag>.spill_*.root after hadd
--help                     print usage
```

By default:

```bash
CLEAN_TRIGGERED_SPILLS=0
```

so the per-spill triggered files are kept after the final `hadd`. This is intentional: they are useful for debugging individual spills and checking `spill_participation` propagation. Pass `--clean-triggered-spills` only when those intermediate triggered files should be removed.

The input prefix must match the output prefix from `process.sh`. For normal run-level triggering, do not pass `--devices`; both workflows use:

```text
aps.sorted
```

and the final run-level triggered file is:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.root
```

Use `trigger.sh --devices ...` only for special diagnostic workflows where `process.sh` was also intentionally run with the same `--devices` subset, for example same-RDO calibration closure checks. In that case, pass the same `--devices` list to `trigger.sh`:

```bash
process/scripts/process.sh \
  --run RUN_NAME \
  --devices rdo-192 \
  --calibration calibration.conf

process/scripts/trigger.sh \
  --run RUN_NAME \
  --devices rdo-192 \
  --trigger trigger.conf calibcheck_rdo-192 \
  --window 32
```

This reads:

```text
/data/2026-testbeam/process/<run>/process/aps.sorted.rdo-192.spill_*.root
```

and writes:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.calibcheck_rdo-192.spill_0000.root
/data/2026-testbeam/process/<run>/trigger/triggered.calibcheck_rdo-192.root
```

## timing.sh

`timing.sh` is the post-trigger TIMING-estimator workflow. It consumes triggered-frame ROOT files produced by `trigger.sh` and adds one set of timing-estimator branches per frame.

Run on the final `hadd` output:

```bash
process/scripts/timing.sh \
  --run RUN_NAME \
  --trigger fingers
```

Required command-line options:

```text
--run RUN                  run name/directory
--trigger TAG              trigger output tag; may be repeated
```

Optional command-line options:

```text
--run-type TYPE            accepted for symmetry with the other workflows, default physics
--parallel-spills          process triggered.<tag>.spill_*.root files in parallel, then hadd
--jobs N                   maximum parallel spill jobs with --parallel-spills, default 8
--overwrite                overwrite existing timing outputs instead of skipping them
--clean-timing-spills      remove triggered.<tag>.timing.spill_*.root after hadd
--help                     print usage
```

Default mode reads:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.root
```

and writes:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.root
```

This is the preferred mode unless the final triggered file is large enough that per-spill parallelism is useful.

With `--parallel-spills`, the script reads:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.spill_*.root
```

writes:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.spill_*.root
```

and combines them into:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.root
```

By default:

```bash
CLEAN_TIMING_SPILLS=0
```

so the per-spill timing files are kept. Pass `--clean-timing-spills` only when those intermediate files should be removed after the final `hadd`.
