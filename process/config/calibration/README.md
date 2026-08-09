# Timing Calibration Configuration

Calibration files are whitespace-separated text files with three sections:

```text
[TDC]
# device fifo column pixel tdc off iif

[CHANNEL]
# device fifo column pixel offset

[TRIGGER]
# device fifo offset
```

Blank lines are allowed. Comments begin with `#`.

Address fields may be exact integers or `*` wildcards. Calibration values must be numeric; wildcards are not allowed in value columns.

## Matching

If multiple rules match the same concrete TDC/channel/trigger source, the most-specific match wins. Specificity is the number of non-wildcard address fields. Equal-specificity matches are ambiguous and cause `calibrator` to fail.

## ALCOR Timing Formula

For ALCOR hits (`type == 1`):

```text
phase = off + iif * fine
time = coarse + 32768 * rollover - phase - channel_offset
```

`off` and `iif` are TDC-level parameters identified by:

```text
device fifo column pixel tdc
```

`channel_offset` is identified by:

```text
device fifo column pixel
```

## Trigger Timing Formula

For trigger tags (`type == 9`):

```text
time = coarse + 32768 * rollover - trigger_offset
```

Trigger offsets are identified by:

```text
device fifo
```

A positive offset means the source arrives late in raw time and is shifted earlier by calibration.

See `calibration_example.conf` for a complete template with global defaults and overrides.


## Producing TDC Rows with dcalib

`process/bin/dcalib` can generate the `[TDC]` rows from a decoded ROOT file:

```bash
process/bin/dcalib \
  --input decoded.root \
  --output dcalib.root \
  --calibration-output tdc_calibration.conf
```

The text output contains only:

```text
[TDC]
# device fifo column pixel tdc off iif
...
```

Outputs from multiple `dcalib` runs can be concatenated with `cat` to create a larger TDC calibration snippet. If two files contain the same concrete TDC row, `calibrator` will reject the final configuration as a duplicate rule.
