# TIMING Event Estimator

This directory contains a native C++ port of the trained Python timing estimator in `python/`.

The estimator is a black-box model trained for two 32-channel TIMING scintillators. The trained package was imported from `timing_event_estimator.zip`; the ZIP file itself is not stored in the repository.

Do not retrain, edit, or silently regenerate the model parameters. The Python package under `python/` is the numerical reference implementation.

## Input

Each event is represented by exactly 64 calibrated floating-point timing values:

```text
0..31   TIMING0 DO0..DO31
32..63  TIMING1 DO0..DO31
```

Input times must already use the Cherenkov-860 timing convention:

```text
calibrated_time = hit.time - offset_860
```

No additional calibration is applied by `TimingEstimator`.

The native time unit is:

```text
3.125 ns
```

## C++ API

The public API is:

```cpp
#include "TimingEstimator.h"

std::array<double, 64> times;

TimingEstimator estimator;
TimingResult result = estimator.estimate(times);
```

`TimingResult` contains:

```cpp
double T0;      // TIMING0 estimate, native units
double sigma0;  // TIMING0 uncertainty, native units
double T1;      // TIMING1 estimate, native units
double sigma1;  // TIMING1 uncertainty, native units
double T;       // (T0 + T1) / 2, native units
double sigmaT;  // sqrt(sigma0^2 + sigma1^2) / 2, native units
```

Convenience conversion helpers are also available:

```cpp
result.T0_ns();
result.T1_ns();
result.T_ns();
result.sigma0_ps();
result.sigma1_ps();
result.sigmaT_ps();
```

## Model Conversion

The Python estimator has two model families.

The time estimators are four scikit-learn `HistGradientBoostingRegressor` models:

```text
timeA_timing0.joblib
timeB_timing0.joblib
timeA_timing1.joblib
timeB_timing1.joblib
```

They are converted to native C++ decision-tree arrays in:

```text
src/TimingEstimatorData.inc
```

The C++ prediction is exactly the scikit-learn raw tree traversal:

```text
prediction = baseline + sum(reached_leaf_value)
```

The sigma estimators are two small PyTorch networks:

```text
52 -> 64 SiLU -> 32 SiLU -> 1
```

Their weights are also stored in `src/TimingEstimatorData.inc` and evaluated directly in C++ with float arithmetic.

There is no Python, scikit-learn, joblib, or PyTorch dependency during normal C++ event processing.

## Feature Construction

For each scintillator independently:

1. Sort the 32 channel times.
2. Define anchor `B` as the mean of the 10 earliest channel times.
3. Build Model A features:

```text
32 channel times relative to B
24 earliest sorted times relative to B
6 order-statistic summaries
```

4. Build Model B / sigma features:

```text
32 channel times relative to B
20 earliest sorted times relative to B
```

5. Standardize the 52 Model B / sigma features with the trained per-scintillator scaler.

This preserves the translation-equivariant structure of the Python estimator. If all TIMING0 inputs are shifted by `delta`, then `T0` shifts by `delta` and `sigma0` is unchanged. TIMING1 behaves independently in the same way.

## Output Convention

The estimator returns times and sigmas in native units. Since one native unit is `3.125 ns`, sigma values can be converted to ps by multiplying by `3125`.

The final event estimate is:

```text
T      = (T0 + T1) / 2
sigmaT = sqrt(sigma0^2 + sigma1^2) / 2
```

## Validation

Build the project:

```bash
cmake -S process/source -B process/build
cmake --build process/build --target validate_timing_estimator -j
```

This builds:

```text
process/build/libtiming_estimator.so
process/build/validate_timing_estimator
```

For interactive ROOT use, load the shared library and include the API header:

```cpp
gSystem->Load("process/build/libtiming_estimator.so");
#include "process/source/timing/include/TimingEstimator.h"
```

Run the validation wrapper:

```bash
process/source/timing/test/validate_timing_estimator.sh
```

The wrapper:

1. Generates deterministic 64-channel test events.
2. Runs the Python reference estimator from `python/`.
3. Runs the native C++ estimator on the same events.
4. Prints max absolute and RMS differences for:

```text
T0
sigma0
T1
sigma1
T
sigmaT
DeltaT = T1 - T0
```

On the development machine, the 1000-event validation gave:

```text
tested events: 1000
            T0 max_abs=5.68434188608e-14 rms=1.35143524006e-14
        sigma0 max_abs=4.71565167137e-07 rms=1.10062483677e-07
            T1 max_abs=5.68434188608e-14 rms=1.15623259972e-14
        sigma1 max_abs=2.21455370708e-07 rms=6.35453483922e-08
             T max_abs=5.68434188608e-14 rms=1.09895194901e-14
        sigmaT max_abs=2.16150634613e-07 rms=5.19970482372e-08
        DeltaT max_abs=8.52651282912e-14 rms=1.79761081434e-14
```

The time estimators therefore reproduce Python at numerical precision. The sigma differences are from direct C++ float network evaluation versus PyTorch and are far below practical analysis precision.

The validation also checks independent time-translation equivariance for TIMING0 and TIMING1.

## Python Environment Note

The bundled Python reference package requires:

```text
numpy
scikit-learn
joblib
torch
```

The model files were serialized with scikit-learn 1.8. If the system Python has an older scikit-learn, use a temporary Python environment for validation. This affects validation only; normal C++ inference does not use Python.
