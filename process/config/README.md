# Configuration Files

Configuration files are split by purpose:

```text
calibration/     low-level timing calibration for calibrator
trigger/         declarative trigger definitions for trigger
```

Calibration belongs upstream in `calibrator`. Trigger selector `time_offset` values should normally be zero when running on calibrated data, unless an analysis intentionally applies an additional high-level relative timing selection.
