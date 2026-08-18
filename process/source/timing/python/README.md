# TIMING event-time estimator

Black-box estimator trained on `timing_events_860.txt`. The packaged time estimator is the exact final two-model boosted-tree blend used in the analysis; the uncertainty estimator is the calibrated per-scintillator heteroscedastic network.

## Input

Exactly 64 calibrated floating-point times per event:

- columns 0..31: TIMING0 DO0..DO31
- columns 32..63: TIMING1 DO0..DO31

The input must already use the Cherenkov-860 convention:

`calibrated_time = hit.time - offset_860`

Native time unit: **3.125 ns**.

## Output

For each event:

- `T0`: best estimated time from TIMING0
- `sigma0`: predicted event-by-event uncertainty of TIMING0
- `T1`: best estimated time from TIMING1
- `sigma1`: predicted event-by-event uncertainty of TIMING1
- `T_event = (T0 + T1)/2`
- `sigma_event = sqrt(sigma0^2 + sigma1^2)/2`

The inference API provides times in native units and ns, and sigmas in native units and ps.

## Python

```python
from estimator import TimingEventEstimator

est = TimingEventEstimator("/path/to/timing_event_estimator")
r = est.predict(event64)

print(r.T0_native)
print(r.sigma0_ps)
print(r.T1_native)
print(r.sigma1_ps)
print(r.T_event_native)
print(r.sigma_event_ps)
```

For many events, use `predict_batch(array_N_by_64)`.

## Important uncertainty note

The quadrature combination `sqrt(sigma0^2 + sigma1^2)` is calibrated against the measured `DeltaT = T1 - T0`.

Without an external timing reference, the *global* decomposition of that variance into TIMING0 and TIMING1 is not uniquely measurable. The training uses a neutral gauge convention. Event-by-event changes in each predicted sigma remain useful, while `sigma_event` is the quantity directly calibrated by the TIMING0/TIMING1 comparison.

## Dependencies

- numpy
- scikit-learn
- joblib
- torch
