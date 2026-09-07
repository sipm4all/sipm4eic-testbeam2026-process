# IRT and Synthetic Ellipse Data

`irt` performs inverse ray tracing on the current triggered data format. It
reads the `cherenkov` tree and, when fitting is requested, uses the accepted
entries in the `ring` tree to select the fit hits.

## Synthetic ellipse generator

`irt_ellipse` is a validation utility for the IRT geometry fit. It reads the
`hResults` histogram produced by `macros/example/hitmap_fit_ellipse.C` and
extracts:

```text
X0, Y0, A, B, theta, sigmaRho
```

It generates a configurable number of synthetic frames. For every hit it
draws an independent random azimuth uniformly around the ellipse, then draws
an independent Gaussian normalized radial fluctuation:

```text
rho = 1 + Gaussian(0, sigmaRho)
```

The generated file contains frame-aligned `frames`, `cherenkov`, and `ring`
trees. The ring entry is synthetic metadata describing the generated ellipse,
so the ordinary `irt --fit` program can consume the file without a special fit
path. `Nbkg` is deliberately not used.

Example:

```bash
process/build/irt_ellipse \
  --input hitmap.root \
  --output synthetic.root \
  --events 10000 \
  --hits 32 \
  --seed 12345
```

The ellipse parameters can be overridden individually. An override takes
precedence over `hResults`:

```bash
--x0 VALUE --y0 VALUE --a VALUE --b VALUE \
--theta VALUE --sigma-rho VALUE
```

For example, an ideal ellipse without radial smearing is generated with:

```bash
--sigma-rho 0
```

## Running the normal IRT fit

The generated file is then processed by the standard executable:

```bash
process/build/irt \
  --input synthetic.root \
  --output reconstructed.root \
  --config process/config/geometry/nominal.conf \
  --fit-config process/config/geometry/fit.sps1.conf \
  --fit \
  --max-frames 10000 \
  --fit-max-hits 2000
```

This two-step procedure is useful because it tests the real IRT input and fit
path. It does not use the ellipse generator as an alternative fitting
implementation.
