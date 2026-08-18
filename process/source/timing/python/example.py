
import numpy as np
from estimator import TimingEventEstimator

est = TimingEventEstimator()

# One event: exactly 64 calibrated values.
event = np.loadtxt("one_event.txt")
result = est.predict(event)
print(result.as_dict())

# Full timing.txt / timing_events_860.txt:
data = np.loadtxt("timing_events_860.txt")
out = est.predict_batch(data)
np.savetxt(
    "event_times.txt",
    np.column_stack([
        out["T0_native"],
        out["sigma0_ps"],
        out["T1_native"],
        out["sigma1_ps"],
        out["T_event_native"],
        out["sigma_event_ps"],
    ]),
    header="T0_native sigma0_ps T1_native sigma1_ps T_event_native sigma_event_ps",
)
