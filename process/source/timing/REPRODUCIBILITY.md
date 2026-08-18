# TIMING Event-Time Estimator Package
## Provenance, contents, training procedure, validation, and reproduction guide

**Package:** `timing_event_estimator.zip`  
**Reference dataset:** `timing_events_860.txt`  
**Date of package creation:** 2026-08-18  
**Purpose:** extract the best available per-event time estimate and uncertainty from each of the two 32-channel TIMING scintillators, then form a combined event time.

---

## 1. Executive summary

The package contains a trained black-box estimator for two scintillators, `TIMING0` and `TIMING1`. Each scintillator is read out by 32 SiPM channels. For one event, the input is therefore 64 calibrated timing values:

```text
0..31   TIMING0 DO0..DO31
32..63  TIMING1 DO0..DO31
```

The input values must already use the Cherenkov-860 calibration convention:

```text
calibrated_time = hit.time - offset_860
```

The native time unit is:

```text
1 native unit = 3.125 ns = 3125 ps
```

For each event the estimator returns:

```text
T0       best time estimate from TIMING0
sigma0   predicted event-by-event uncertainty of TIMING0

T1       best time estimate from TIMING1
sigma1   predicted event-by-event uncertainty of TIMING1

T_event  = (T0 + T1) / 2
sigma_event = sqrt(sigma0^2 + sigma1^2) / 2
```

The package was developed by optimizing the width of

```text
DeltaT = T1 - T0
```

on a large sample of 288,837 events and validating on a held-out test sample that was not used for training or model selection.

The final packaged timing estimator gives, on the baseline selected test sample:

```text
sigma68(DeltaT)       = 187.999 ps
sigma68(T_event)      = 93.999 ps

RMS(DeltaT)           = 227.494 ps
RMS(T_event)          = 113.747 ps
```

The predicted per-event uncertainty is useful for selecting the best-timed events. On the same held-out test sample:

```text
fraction retained      sigma68(T_event)
----------------------------------------
90%                     88.73 ps
80%                     84.83 ps
72%                     82.75 ps
70%                     82.13 ps
60%                     79.74 ps
50%                     77.35 ps
30%                     72.92 ps
20%                     70.58 ps
10%                     68.37 ps
 5%                     66.58 ps
```

These efficiencies are relative to the baseline quality-selected sample used to train the final timing models. The baseline quality selection retains about 72% of all events.

The final package is intended primarily as a **frozen inference artifact**. It is suitable for conversion to native C++/ROOT inference, provided that the C++ implementation is checked event-by-event against the Python reference.

---

# 2. Detector and input data

## 2.1 Physical system

The TIMING system consists of two scintillators, `TIMING0` and `TIMING1`.

Each scintillator is approximately:

```text
20 x 20 x 30 mm
```

One 20 x 30 mm face is coupled to a matrix of 32 SiPM channels arranged as a 4 x 8 grid.

For detector output channel `DO`:

```text
x = DO % 4
y = DO / 4
```

so the detector channel coordinates are:

```text
x = 0..3
y = 0..7
```

The two scintillators are close together along the beam direction and are crossed by the same particle in an event.

The intended final observable is not tracking. The purpose of the system is a precise event-time measurement:

```text
T_event = (T0 + T1) / 2
```

with the resolution inferred from the difference:

```text
DeltaT = T1 - T0
```

Under the assumption that the measurement errors of the two scintillators are independent,

```text
sigma_event = sigma(DeltaT) / 2
```

without requiring the two scintillators to have equal individual resolutions.

---

## 2.2 Reference input file

The full reference sample used during development was:

```text
timing_events_860.txt
```

It contains:

```text
288,837 events
64 floating-point values per event
one event per line
no NaN values
```

The SHA-256 checksum of the exact reference file available during package creation was:

```text
d9591c775bf0b7bae9d0ce42d5fe7e0b41a42d655dabe5a5c4fadeecf785a7e5
```

The mounted filename during package creation was `timing_events_860(1).txt`; its contents correspond to the full uploaded `timing_events_860.txt` sample.

The column order is fixed:

```text
columns  0..31  = TIMING0 DO0..DO31
columns 32..63  = TIMING1 DO0..DO31
```

---

## 2.3 Input calibration

The timing values are already channel calibrated.

The calibration source was:

```text
/home/preghenella/CODEX/sorting/timing_offsets_from_cherenkov860.conf
```

The applied convention before writing the text file was:

```text
calibrated_time = hit.time - offset_860
```

The packaged estimator **must not apply these offsets again**.

The original electronics-output to detector-output mapping used before writing the file was:

```cpp
eo2do[32] = {
  22,20,18,16, 24,26,28,30,
  25,27,29,31, 23,21,19,17,
   9,11,13,15,  7, 5, 3, 1,
   6, 4, 2, 0,  8,10,12,14
};
```

The package assumes the file is already in detector-output order, i.e. DO0..DO31 for each scintillator.

---

# 3. Why the estimator is translation-equivariant

The absolute timing values contain large clock values. The physically meaningful information used by the estimator is the relative timing pattern inside each scintillator.

For each scintillator, the package first sorts the 32 channel times and defines an anchor:

```text
B = mean of the 10 earliest channel times
```

All model features are then constructed from differences relative to `B`, for example:

```text
t_i - B
sorted_t_k - B
sorted_t_k - sorted_t_1
```

The trained model predicts only a correction to `B`.

Therefore, if all 32 channels of one scintillator are shifted by a common amount `delta`:

```text
t_i -> t_i + delta
```

then all relative features are unchanged and:

```text
T_i -> T_i + delta
sigma_i -> unchanged
```

This property was explicitly verified after packaging.

For a 5,000-event numerical check, the packaged estimator reproduced the reference implementation with:

```text
maximum |DeltaT_package - DeltaT_reference| < 0.001 ps
```

and arbitrary independent time shifts applied to TIMING0 and TIMING1 were reproduced to floating-point precision.

---

# 4. Data split

A fixed pseudo-random split was used throughout the final analysis.

Random generator:

```python
rng = np.random.default_rng(20260818)
```

The 288,837 events were shuffled once and split as:

```text
training raw:     173,302 events   60%
validation raw:    57,767 events   20%
test raw:          57,768 events   20%
```

The test sample was kept separate from training and model selection.

---

# 5. Baseline quality selection used for final model training

Before the final time models were trained, a detector-consistency selection was applied. This selection was developed during the exploratory analysis and was used to define a relatively clean training domain.

The package itself does **not** apply these cuts automatically. They are documented here because they define the population on which the final timing models and uncertainty model were calibrated.

## 5.1 Position estimate

For each scintillator, a position proxy was reconstructed from the earliest 10 channels.

The 10 earliest channels were ranked by time, and their detector coordinates were combined using exponentially decreasing rank weights:

```python
w_k = exp(-k / 1.6)
```

with `k = 0..9`.

The resulting position estimate was:

```text
p = weighted mean of DO coordinates of the first 10 hits
```

TIMING1 was affine-aligned to TIMING0 using only the raw training split.

The difference variables were:

```text
DeltaX = x1_aligned - x0
DeltaY = y1_aligned - y0
```

The robust central core obtained from the training sample was approximately:

```text
mu_X    = -0.02177
mu_Y    = -0.01067
sigma_X =  0.30243 grid units
sigma_Y =  0.21904 grid units
```

The position-consistency variable was:

```text
Rpos = sqrt(
    ((DeltaX - mu_X) / sigma_X)^2 +
    ((DeltaY - mu_Y) / sigma_Y)^2
)
```

and the cut was:

```text
Rpos < 2
```

---

## 5.2 First-hit consistency

For each detector, the earliest channel was found:

```text
first_DO = argmin(channel times)
```

The distance in detector-grid coordinates between the first-hit DO and the reconstructed position was calculated.

The event variable was:

```text
Qfirst = max(
    distance(first_DO_TIMING0, position_TIMING0),
    distance(first_DO_TIMING1, position_TIMING1)
)
```

The central distribution was characterized by:

```text
mu_Q    = 0.4834567871
sigma_Q = 0.1353945306
```

and the final cut was:

```text
Qfirst < mu_Q + 2*sigma_Q
       < 0.7542458482
```

---

## 5.3 Selected sample sizes

After both cuts:

```text
training:    124,755 events
validation:   41,748 events
test:         41,553 events
```

The total efficiency of these cuts over the full 288,837-event sample was:

```text
72.03%
```

The efficiency in the held-out raw test split was:

```text
71.93%
```

---

# 6. Exploratory studies that led to the final black-box approach

The final black-box design was not assumed at the beginning. Several physically motivated estimators were investigated first.

These studies are not required for inference, but they explain why the final package uses the full 32-channel timing pattern.

## 6.1 Earliest-N averaging

For each scintillator, the channel times were sorted:

```text
t_(1) <= t_(2) <= ... <= t_(32)
```

The simple estimator:

```text
T_N = mean of the first N times
```

showed a broad optimum around the first 10–12 hits.

A full two-dimensional scan of averaging an ordered window `M..N` showed that discarding the earliest hit did not improve the resolution. The optimum remained on the `M = 0` boundary.

This established that the earliest hits contain useful information, but simple averaging was not optimal.

## 6.2 Position/propagation corrections

Per-channel position-dependent empirical timing corrections improved the simple first-N estimator, demonstrating that detector geometry and optical propagation mattered.

Position-dependent weighting and local covariance models gave additional, but smaller, improvements.

## 6.3 Full timing pattern

The strongest improvement came when the estimator was allowed to use the entire relative timing pattern across the 32 SiPMs.

This indicated that later hits and channel-to-channel correlations contain information about event-specific optical fluctuations even if those late hits are not directly useful to average as time measurements.

At that point the optimization goal was simplified to:

```text
Use whatever information exists in the 32-channel pattern to minimize DeltaT,
subject to exact time-translation equivariance.
```

That is the philosophy of the final package.

---

# 7. Final time estimator

The final per-scintillator timing estimate is an ensemble of two independently trained `HistGradientBoostingRegressor` models.

For each scintillator:

```text
T_i = anchor B_i + learned correction
```

The model uses only that scintillator's own 32 channel measurements.

No information from the opposite scintillator is needed at inference time to calculate `T0` or `T1`.

---

## 7.1 Anchor

For either detector:

```python
s = np.sort(t, axis=1)
B = s[:, :10].mean(axis=1)
```

where `t` contains the 32 calibrated channel times.

---

## 7.2 Model A features

Model A uses 62 translation-invariant features.

### 32 channel-relative times

```text
t_i - B, i = 0..31
```

### 24 sorted relative times

```text
s_k - B, k = 1..24
```

### 6 order-statistic summary variables

Using zero-based Python indexing as in the implementation:

```text
s[0]  - B
s[1]  - s[0]
s[3]  - s[0]
s[7]  - s[0]
s[15] - s[0]
s[23] - s[0]
```

Total:

```text
32 + 24 + 6 = 62 features
```

These features are stored and evaluated as `float32` after the differences have been formed from `float64` input times.

---

## 7.3 Model A hyperparameters

The stored `HistGradientBoostingRegressor` parameters are:

```text
loss                 = squared_error
learning_rate        = 0.07
max_iter             = 180
max_leaf_nodes       = 31
min_samples_leaf     = 80
l2_regularization    = 2.0
max_bins             = 255
max_depth            = None
random_state         = 42
early_stopping       = auto
validation_fraction  = 0.1
n_iter_no_change     = 10
tol                   = 1e-7
```

Separate Model-A regressors are stored for TIMING0 and TIMING1.

---

## 7.4 Model B features

Model B uses 52 features:

### 32 channel-relative times

```text
t_i - B, i = 0..31
```

### 20 sorted relative times

```text
s_k - B, k = 1..20
```

Total:

```text
32 + 20 = 52 features
```

These 52 features are standardized separately for TIMING0 and TIMING1 using the stored `StandardScaler` objects.

---

## 7.5 Model B hyperparameters

The stored Model-B regressors use:

```text
loss                 = absolute_error
learning_rate        = 0.06
max_iter             = 140
max_leaf_nodes       = 63
min_samples_leaf     = 120
l2_regularization    = 6.0
max_bins             = 255
max_depth            = None
random_state         = 73
early_stopping       = auto
validation_fraction  = 0.1
n_iter_no_change     = 10
tol                   = 1e-7
```

Again, separate regressors are stored for TIMING0 and TIMING1.

---

# 8. How the time regressors were trained

The only directly measurable timing target available from the two scintillators is their difference.

Define:

```text
dB = B1 - B0
```

The goal is to learn detector-specific corrections such that:

```text
DeltaT = T1 - T0
```

is as narrow as possible.

For an additive pair of detector models:

```text
T0 = B0 - p0(X0)
T1 = B1 + p1(X1)
```

therefore:

```text
DeltaT = dB + p1(X1) + p0(X0)
```

The regression target is:

```text
y = -dB
```

and the two detector functions are fitted by backfitting: one detector is fitted to the residual left after the current estimate from the other detector.

Conceptually:

```python
p0 = 0
p1 = 0

repeat:
    fit p1(X1) to y - p0
    fit p0(X0) to y - p1
```

Model A used the 62-feature representation and **3 backfitting cycles**. Model B used the 52-feature standardized representation and **2 backfitting cycles**.

For Model A the stored regressors use `random_state = 42`; for Model B they use `random_state = 73`. The final stored objects are the regressors from the last backfitting cycle for each detector.

---

# 9. Final time-model ensemble

The final model is a validation-selected blend:

```text
a = 0.65
```

with:

```text
65% Model A
35% Model B
```

The exact inference equations stored in the package are:

```text
T0 = B0 - [0.65 * p0_A + 0.35 * p0_B]

T1 = B1 + [0.65 * p1_A + 0.35 * p1_B]
```

and therefore:

```text
DeltaT = T1 - T0
```

The final event time is:

```text
T_event = (T0 + T1) / 2
```

---

# 10. Per-event uncertainty estimator

In addition to `T0` and `T1`, the package predicts an uncertainty for each scintillator:

```text
sigma0
sigma1
```

Each uncertainty model uses only the 52 standardized features of its own scintillator.

---

## 10.1 Sigma network architecture

For both TIMING0 and TIMING1:

```text
input: 52

Dense(52 -> 64)
SiLU

Dense(64 -> 32)
SiLU

Dense(32 -> 1)
```

The final network output is transformed using:

```python
variance = softplus(network_output) + 1e-6
sigma = sqrt(variance)
```

The trained weights are stored in:

```text
sigmanet_timing0.pt
sigmanet_timing1.pt
```

---

## 10.2 Sigma-network training sample

The uncertainty networks were trained **after** the final time estimator had been fixed.

The out-of-sample validation residuals of the final time estimator were used:

```text
r = DeltaT_predicted
```

The 41,748-event validation subset was randomly divided again for uncertainty training/calibration.

Random generator:

```python
rng = np.random.default_rng(20260818 + 99)
```

i.e.

```text
seed = 20260917
```

Approximately 70% of the validation events were used to train the uncertainty networks and 30% were used for early stopping and variance calibration.

---

## 10.3 Heteroscedastic loss

For one event the networks produce detector variances:

```text
v0
v1
```

The predicted DeltaT variance is:

```text
v = v0 + v1
```

The training loss is a Gaussian heteroscedastic negative log-likelihood:

```text
0.5 * [ log(v) + r^2 / v ]
```

A small gauge regularization term encourages the global average variance assigned to TIMING0 and TIMING1 to remain comparable:

```text
0.02 * [ mean(log(v0)) - mean(log(v1)) ]^2
```

This gauge is necessary because the two-detector comparison directly constrains only:

```text
v0 + v1
```

not the unique global partition between the two detectors.

---

## 10.4 Sigma training optimizer

The uncertainty networks were trained with:

```text
optimizer          = AdamW
learning rate      = 2e-3
weight decay       = 1e-4
batch size         = 4096
maximum epochs     = 120
early-stop patience= 10
gradient clipping  = 5.0
```

The saved weights are the best state according to the held-out uncertainty-calibration subset.

---

## 10.5 Absolute variance calibration

After training, a single multiplicative scale factor was applied to both predicted detector sigmas so that the quadrature prediction matched the observed validation residual scale.

Stored scale:

```text
1.0084866285324097
```

The final package therefore uses:

```text
sigma_i = raw_sigma_i * 1.0084866285324097
```

The event-level predicted uncertainty is:

```text
sigma_event = sqrt(sigma0^2 + sigma1^2) / 2
```

---

# 11. Interpretation and identifiability of sigma0 and sigma1

This point is important.

From the two-counter data alone, the directly measured quantity is:

```text
Var(DeltaT) = Var(T1 - T0)
```

If the individual detector errors are independent:

```text
Var(DeltaT) = sigma0^2 + sigma1^2
```

Therefore the quadrature sum is experimentally constrained.

However, without an external timing reference, the global division:

```text
sigma0 versus sigma1
```

is not uniquely identifiable.

The uncertainty training therefore uses a neutral gauge convention to divide the variance between the two detectors.

Consequences:

1. `sigma_event = sqrt(sigma0^2 + sigma1^2)/2` is the quantity that is directly calibrated by the two-scintillator comparison.
2. Event-to-event variations of `sigma0` and `sigma1` are useful for identifying which individual scintillator is predicted to be well or poorly measured.
3. A statement such as “TIMING0 has an absolute global resolution of exactly X ps and TIMING1 exactly Y ps” cannot be established uniquely from this dataset alone.

An independent external time reference would be required to measure the absolute resolution of the two detectors separately.

---

# 12. Uncertainty calibration performance

The uncertainty prediction was checked on the untouched 41,553-event test sample.

When test events were sorted into deciles of predicted `sigma_event`, the measured RMS changed monotonically with the prediction.

Representative values were:

```text
mean predicted sigma_event   measured event RMS
------------------------------------------------
 66.5 ps                      73.7 ps
 77.8 ps                      79.2 ps
 85.0 ps                      83.1 ps
 91.3 ps                      89.1 ps
 97.8 ps                      93.5 ps
104.9 ps                     100.7 ps
113.3 ps                     103.8 ps
124.5 ps                     118.3 ps
142.7 ps                     140.9 ps
218.9 ps                     197.7 ps
```

The predicted uncertainty should primarily be interpreted as an RMS/variance prediction.

Because the DeltaT distributions retain non-Gaussian tails, the central-68% width is generally smaller than the RMS-based uncertainty.

---

# 13. Baseline packaged performance

The following numbers were recomputed directly with the final packaged estimator on 2026-08-18, using the exact held-out test indices and baseline quality selection documented above.

Test events:

```text
41,553
```

Results:

```text
sigma68(DeltaT)   = 187.9987 ps
sigma68(T_event)  =  93.9994 ps

RMS(DeltaT)       = 227.4937 ps
RMS(T_event)      = 113.7469 ps
```

The package does not automatically reject events.

If the user selects events by the package's predicted `sigma_event`, the performance improves continuously as lower-quality events are rejected.

For example:

```text
retained   sigma_event cut   measured sigma68(T_event)
-------------------------------------------------------
90%        < 157.27 ps        88.73 ps
80%        < 131.78 ps        84.83 ps
72%        < 120.47 ps        82.75 ps
70%        < 118.21 ps        82.13 ps
60%        < 108.91 ps        79.74 ps
50%        < 101.13 ps        77.35 ps
30%        <  88.21 ps        72.92 ps
20%        <  81.56 ps        70.58 ps
10%        <  73.60 ps        68.37 ps
 5%        <  67.75 ps        66.58 ps
```

These retained fractions are relative to the 41,553-event baseline selected test sample.

---

# 14. Residual channel-offset test

The input data already use the Cherenkov-860 per-channel calibration.

A later explicit test added one free residual offset per DO channel after the black-box timing estimator.

The fitted residual channel offsets were small:

```text
TIMING0 offset RMS ~ 5.9 ps
TIMING1 offset RMS ~ 5.0 ps

largest |residual offset|:
TIMING0 ~ 14.3 ps
TIMING1 ~ 11.9 ps
```

Adding these residual offsets changed the event-time resolution by less than 1 ps.

Therefore the final package does **not** include an extra residual channel-offset correction.

The existing Cherenkov-860 calibration is considered sufficiently accurate for the trained model.

---

# 15. Exact package contents

The ZIP contains a top-level directory:

```text
timing_event_estimator/
```

with 13 files.

## 15.1 Files

```text
README.md
estimator.py
example.py
metadata.json
requirements.txt

scaler_timing0.joblib
scaler_timing1.joblib

timeA_timing0.joblib
timeA_timing1.joblib

timeB_timing0.joblib
timeB_timing1.joblib

sigmanet_timing0.pt
sigmanet_timing1.pt
```

---

## 15.2 Meaning of each file

### `estimator.py`

Canonical Python inference implementation.

It defines:

```text
TimingEventEstimator
TimingEstimate
SigmaNet
```

and contains the exact feature construction and final ensemble formulas.

For any future port to C++, this file should be treated as the numerical reference.

### `metadata.json`

Machine-readable metadata containing:

```text
input order
native time unit
calibration convention
feature definitions
ensemble blend coefficients
sigma calibration factor
```

### `README.md`

Short user-facing package instructions.

### `example.py`

Minimal Python example showing single-event and batch inference.

### `requirements.txt`

Minimal Python package dependencies.

### `scaler_timing0.joblib`, `scaler_timing1.joblib`

The two `StandardScaler` instances used for the 52-feature Model-B and sigma-network inputs.

### `timeA_timing0.joblib`, `timeA_timing1.joblib`

The Model-A `HistGradientBoostingRegressor` objects using the 62-feature representation.

### `timeB_timing0.joblib`, `timeB_timing1.joblib`

The Model-B `HistGradientBoostingRegressor` objects using the 52-feature standardized representation.

### `sigmanet_timing0.pt`, `sigmanet_timing1.pt`

PyTorch `state_dict` files for the two 52 -> 64 -> 32 -> 1 uncertainty networks.

---

# 16. Package checksums

The final ZIP SHA-256 is:

```text
3290c488ad6b73d628d2bdec22018e6767380621a1c2420a9ae960da49bfa9de
```

Individual package-file SHA-256 checksums at creation time:

```text
86176736bca2766828c7fbe09f3e808ef3b18b6f7a9624f51ba20bbf5f172349  README.md
d18280bb7ea37c67175aa66cfb20f44cfd9a73e715217c0adaf3dd8899576a99  estimator.py
7ce348ff3db0afe266a960fbb749c01c74fb447a96b64b97825358399d286a6e  example.py
d4f1402c8e6960d30c3b9e6ae14f24f8601d1858b2cca496a19f82dbfc9d8b59  metadata.json
31f97c6b1a1c9030f7c007b689a53f2a98be1412655665a395b158f7c02f6c1b  requirements.txt
e7d114c37b13f0ecc8852767f79b2519a468e6e2403dce83760dc1282b2f8423  scaler_timing0.joblib
2852678e5d9091aba187031164409add0eaa6a7a40e172f1022e52bda1841af4  scaler_timing1.joblib
c2e53b96c0ac0f0d7231f03cf00e70386d0016ff112d13f1761a9a54e24b8fd4  timeA_timing0.joblib
60d59f9743ae0f476a837348c55ec09f01991aa0b755d2fb8ec289c0289a8136  timeA_timing1.joblib
001b9feb2eba220dcc714abac9cc387efd55e347b43a489e94e00870e3067804  timeB_timing0.joblib
b334539592300bb04a3ed3bed9b9a224f196c502afd03d36b81e427060290878  timeB_timing1.joblib
75b9b2d4e0e985a85f69cdc7eadeb933a7a1bf648410ab5c251f9f1c5531853c  sigmanet_timing0.pt
a92042f2dc619d422e9c2b263d13753a92f9aebaca800313f136ec895ec892f2  sigmanet_timing1.pt
```

**Note:** if this document is added to the package and the ZIP is recreated, the ZIP checksum will necessarily change. Individual files should be re-hashed after any repackaging.

---

# 17. Software environment used when the package was inspected/finalized

The runtime environment during final package inspection was:

```text
Python       3.13.5
NumPy        2.3.5
scikit-learn 1.8.0
joblib       1.5.3
PyTorch      2.10.0+cpu
```

For long-term reproducibility it is recommended to preserve these versions, for example with a dedicated environment/lock file or container image.

The current `requirements.txt` intentionally contains only package names and does not pin exact versions.

---

# 18. Exact inference algorithm

For one scintillator with input vector:

```text
t[0..31]
```

the reference implementation performs:

```python
# keep float64 for absolute times
s = sort(t)
B = mean(s[0:10])

# Model A features
XA = concatenate(
    t - B,
    s[0:24] - B,
    [
        s[0]  - B,
        s[1]  - s[0],
        s[3]  - s[0],
        s[7]  - s[0],
        s[15] - s[0],
        s[23] - s[0],
    ]
).astype(float32)

# Model B / sigma features
XB = concatenate(
    t - B,
    s[0:20] - B
).astype(float32)

XB_scaled = detector_specific_scaler.transform(XB)

pA = detector_specific_ModelA.predict(XA)
pB = detector_specific_ModelB.predict(XB_scaled)
```

For TIMING0:

```text
T0 = B0 - (0.65*p0_A + 0.35*p0_B)
```

For TIMING1:

```text
T1 = B1 + (0.65*p1_A + 0.35*p1_B)
```

For the uncertainty:

```python
z = SigmaNet(XB_scaled)
raw_variance = softplus(z) + 1e-6
raw_sigma = sqrt(raw_variance)
sigma = raw_sigma * 1.0084866285324097
```

Finally:

```text
T_event = (T0 + T1) / 2

sigma_event = sqrt(sigma0^2 + sigma1^2) / 2
```

---

# 19. Minimal Python use

```python
from estimator import TimingEventEstimator

est = TimingEventEstimator("/path/to/timing_event_estimator")

result = est.predict(event64)

print(result.T0_native)
print(result.sigma0_ps)

print(result.T1_native)
print(result.sigma1_ps)

print(result.T_event_native)
print(result.sigma_event_ps)
```

Batch operation:

```python
import numpy as np

events = np.loadtxt("timing_events_860.txt")

out = est.predict_batch(events)

T0 = out["T0_native"]
T1 = out["T1_native"]

DeltaT = T1 - T0

Tevent = out["T_event_native"]
sigma_event_ps = out["sigma_event_ps"]
```

---

# 20. Reproducing the baseline validation

A future user wishing to verify the frozen package should perform the following steps.

## Step 1 — Verify the input file

Check:

```text
288837 rows
64 columns
all finite
correct DO channel order
already Cherenkov-860 calibrated
```

If the original reference file is available, check:

```text
SHA256 = d9591c775bf0b7bae9d0ce42d5fe7e0b41a42d655dabe5a5c4fadeecf785a7e5
```

## Step 2 — Recreate the split

```python
rng = np.random.default_rng(20260818)

idx = np.arange(288837)
rng.shuffle(idx)

ntrain = int(0.60 * len(idx))
nval   = int(0.20 * len(idx))

train_raw = idx[:ntrain]
val_raw   = idx[ntrain:ntrain+nval]
test_raw  = idx[ntrain+nval:]
```

This gives:

```text
173302 / 57767 / 57768
```

## Step 3 — Recreate the baseline quality selection

Reconstruct the rank-weighted positions, align TIMING1 to TIMING0 using the raw training set, calculate `Rpos` and `Qfirst`, and apply:

```text
Rpos < 2
Qfirst < 0.7542458482
```

Expected selected sizes:

```text
train = 124755
val   = 41748
test  = 41553
```

Small floating-point/version differences in the robust position fit should be checked if a few events near a cut boundary differ.

## Step 4 — Run frozen inference on selected test events

Expected:

```text
sigma68(DeltaT)  ~ 187.999 ps
RMS(DeltaT)      ~ 227.494 ps
```

and therefore:

```text
sigma68(T_event) ~ 93.999 ps
RMS(T_event)     ~ 113.747 ps
```

## Step 5 — Check uncertainty ranking

Calculate:

```text
sigma_event = sqrt(sigma0^2 + sigma1^2) / 2
```

and retain the lowest predicted-sigma events.

For example, retaining the lowest 50% of predicted `sigma_event` should give approximately:

```text
sigma_event cut ~ 101.13 ps
sigma68(T_event) ~ 77.35 ps
```

on the reference held-out test sample.

---

# 21. Reproducing the training from scratch

There are two different meanings of “reproduce” and they should not be confused.

## 21.1 Exact inference reproduction

This is fully reproducible from the package.

The package contains:

```text
all trained tree parameters
all scaler means/scales
all neural-network weights
all inference formulas
all feature definitions
```

A future Python or C++ implementation should be able to reproduce the frozen estimator numerically event by event.

This is the recommended meaning of reproducibility for production use.

## 21.2 Exact retraining reproduction

The training procedure is documented above, including the principal NumPy split seeds and the stored model hyperparameters.

However, the original exploratory training session did **not** explicitly freeze every source of low-level numerical randomness for the PyTorch sigma-network training.

In particular, no explicit `torch.manual_seed(...)` was recorded in the final exploratory uncertainty-training step.

Therefore:

```text
retraining the sigma networks from scratch is expected to reproduce the method and performance,
but is not guaranteed to reproduce the saved .pt weights bit-for-bit.
```

This is why the trained `.pt` state dictionaries in the ZIP should be treated as the canonical uncertainty model.

For any future retraining campaign, the training script should explicitly set:

```python
np.random.seed(...)
torch.manual_seed(...)
```

and, if CUDA is used, the appropriate deterministic CUDA settings.

The retraining script, exact environment, training-data checksum, and resulting model hashes should then be archived together.

---

# 22. Important methodological caveat: common-mode time bias

The final time estimator is optimized through:

```text
DeltaT = T1 - T0
```

This constrains differential timing very strongly.

However, with only the two scintillators there is no external measurement of the true event time.

A correction that is identical in both detectors:

```text
T0 -> T0 + h(event)
T1 -> T1 + h(event)
```

would cancel in `DeltaT` but shift:

```text
T_event = (T0 + T1)/2
```

The translation-equivariant anchor and independent per-scintillator inference strongly constrain pathological solutions, but a truly external reference is required to measure any common-mode event-time bias directly.

Therefore:

```text
sigma(DeltaT)/2
```

should be interpreted as the event-time resolution under the standard assumption that the residual errors of TIMING0 and TIMING1 are independent and do not contain a significant shared common-mode component.

This is the same physical assumption normally used when extracting the resolution of a two-counter average from their time difference.

---

# 23. Recommended repository integration

For the repository:

```text
sipm4all/sipm4eic-testbeam2026-process
```

the recommended location is:

```text
process/source/timing/
```

Suggested structure:

```text
process/source/timing/
├── README.md
├── REPRODUCIBILITY.md
├── python/
│   ├── estimator.py
│   ├── example.py
│   ├── metadata.json
│   ├── requirements.txt
│   ├── scaler_timing0.joblib
│   ├── scaler_timing1.joblib
│   ├── timeA_timing0.joblib
│   ├── timeA_timing1.joblib
│   ├── timeB_timing0.joblib
│   ├── timeB_timing1.joblib
│   ├── sigmanet_timing0.pt
│   └── sigmanet_timing1.pt
├── include/
│   └── TimingEstimator.h
├── src/
│   └── TimingEstimator.cxx
└── test/
    └── validate_timing_estimator.C
```

The ZIP itself does not need to be committed if its contents are unpacked and versioned.

The Python estimator should remain in the repository even after a C++ implementation exists, because it is the authoritative numerical reference for regression testing.

---

# 24. Requirements for a C++/ROOT port

A future native C++ implementation should preserve exactly:

1. input DO order;
2. float64 subtraction of the large absolute times before conversion to float32 feature values;
3. the first-10 anchor;
4. exact Model-A feature order;
5. exact Model-B feature order;
6. exact scaler parameters;
7. exact tree thresholds and leaf values;
8. exact 0.65 / 0.35 ensemble blend;
9. exact sigma-network dense weights and biases;
10. `SiLU` activation;
11. `softplus` output transformation;
12. sigma calibration factor;
13. independent time-translation equivariance.

Validation against Python should compare:

```text
T0
sigma0
T1
sigma1
T_event
sigma_event
DeltaT
```

on at least 1,000 common events.

The conversion should not be considered complete until event-by-event numerical agreement has been demonstrated.

---

# 25. Suggested validation tolerances for a C++ port

The frozen Python package reproduced the final analysis reference at effectively floating-point precision.

For a native C++ port, practical targets are:

```text
DeltaT absolute difference:
preferably < 1e-3 ps

T0 / T1:
as close to floating-point precision as the tree/network conversion permits

sigma0 / sigma1:
preferably < 1e-3 ps
```

If a generated inference backend produces slightly different floating-point accumulation ordering, a somewhat looser tolerance can be accepted only after showing that the statistical DeltaT distribution is unchanged.

---

# 26. What must not be changed silently

The following changes define a new estimator and should require a new model/version:

```text
changing the input calibration
changing DO ordering
changing the first-10 anchor definition
changing feature ordering
changing scaler parameters
changing tree parameters
changing neural-network parameters
changing the ensemble coefficient
changing the sigma calibration factor
training on a different dataset
changing the event selection used for training
```

If any of these are changed, the model should receive a new version identifier and a new provenance document.

---

# 27. Recommended future reproducibility improvements

The frozen package is sufficient for exact inference reproduction.

For future training campaigns, the following should additionally be archived:

```text
train.py
exact Python lock file / environment
training dataset checksum
all random seeds
Git commit hash
training log
validation plots
final model checksums
```

A future package should ideally include a machine-executable script:

```text
reproduce_training.py
```

which rebuilds the model from the original timing file and writes a manifest containing all checksums.

---

# 28. Final reference equations

Per detector:

```text
T0 = B0 - [0.65*p0_A(X0_A) + 0.35*p0_B(X0_B)]

T1 = B1 + [0.65*p1_A(X1_A) + 0.35*p1_B(X1_B)]
```

Per-event uncertainties:

```text
sigma0 = SigmaNet0(X0_B) * 1.0084866285324097

sigma1 = SigmaNet1(X1_B) * 1.0084866285324097
```

Final event quantities:

```text
DeltaT = T1 - T0

T_event = (T0 + T1) / 2

sigma_event = sqrt(sigma0^2 + sigma1^2) / 2
```

Native time unit:

```text
1 unit = 3.125 ns = 3125 ps
```

---

# 29. Canonical source of truth

For production and future conversion work, the priority order should be:

1. **trained files in this package**;
2. **`estimator.py`**, which defines the exact inference;
3. **`metadata.json`**, which records the essential model configuration;
4. this provenance/reproduction document;
5. historical exploratory analysis notes.

If any future implementation disagrees with `estimator.py` on the same 64 calibrated input values, the future implementation should be considered incorrect until the discrepancy is understood.

