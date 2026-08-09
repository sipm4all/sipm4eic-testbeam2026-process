# Workflow Scripts

This directory contains shell scripts for running larger processing workflows. The scripts assume they are run from the repository checkout and that executables have been built in the `process/bin/` directory.

Build first:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

## process.sh

`process.sh` is the full processing workflow. It currently performs:

```text
calibrate each decoded FIFO file
sort each calibrated FIFO file
after-pulse suppress each sorted FIFO file
merge lanes per device with --split-spills
merge matching spills across devices
run trigger configurations per spill
hadd triggered spill files per trigger tag
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
--window VALUE             trigger frame window, default 256
--help                     print usage
```

Short aliases are also accepted:

```text
-r, -c, -t, -w, -h
```

`--trigger` may be given multiple times. Tags must be unique because they are used in output filenames such as:

```text
triggered.<tag>.spill_0000.root
triggered.<tag>.root
```

Important variables still configured near the top of the script:

```bash
ipath                  input decoded-data base directory
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
