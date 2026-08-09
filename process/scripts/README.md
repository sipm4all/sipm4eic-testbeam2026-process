# Workflow Scripts

This directory contains shell scripts for running larger processing workflows. The scripts assume they are run from the repository checkout and that executables have been built in the top-level `build/` directory.

Build first:

```bash
cmake -S process/source -B build
cmake --build build -j
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
process/scripts/process.sh RUN_NAME
```

Important variables near the top of the script:

```bash
ipath                  input decoded-data base directory
opath                  output processing base directory
CALIBRATION_CONFIG     calibration configuration file
TRIGGER_CONFIGS        list of trigger configuration files
TRIGGER_TAGS           output tag for each trigger configuration
TRIGGER_WINDOW         frame half-window for trigger
WRITE_LOGS             0 prints to terminal, 1 writes log files
CLEAN_DEVICE_SPILLS    remove intermediate device spill files
CLEAN_MERGED_SPILLS    remove merged spill files after triggering
CLEAN_TRIGGERED_SPILLS remove per-spill trigger output after hadd
```

## trigger.sh

`trigger.sh` is a smaller helper that runs one trigger configuration over already merged split-spill files:

```bash
process/scripts/trigger.sh RUN_NAME
```

It uses `CONFIG` and `WINDOW` variables near the top of the script. For multi-configuration production processing, prefer `process.sh`.
