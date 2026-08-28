# IRT Geometry

The track is fixed in the laboratory frame:

```ini
[track]
origin = (0.0, 0.0, 0.0)
direction = (0.0, 0.0, 1.0)
emission_z = 2534.0
```

The mirror and detector are described in assembly coordinates. `irt` applies
the inverse rigid assembly rotation to the track before ray tracing. This is
equivalent to rotating the complete assembly when the same transformation is
applied to all assembly points and directions.

## Mirror

```ini
[mirror]
center = (1150.08, 0.0, 939.0)
radius = 2203.01
pivot_center = (327.9, 0.0, 3110.6)
rotation_vector = (0.0, 0.0, 0.0)
```

## Detector

```ini
[detector]
sphere_center = (1838.0, 0.0, 1414.0)
sphere_radius = 1100.0
theta = -0.66041320
phi = 0.0
rotation_vector = (0.0, 0.0, 0.0)
```

The angles locate the tangent detector plane on the sphere. Cherenkov hit
`x` and `y` values are local coordinates relative to the plane center.

## Assembly rotation

```ini
[assembly]
rotation_x_pivot = (0.0, 0.0, 2500.0)
rotation_x = 0.0
rotation_y_pivot = (0.0, 0.0, 0.0)
rotation_y = -0.09
```

The x rotation is followed by the y rotation in the assembly convention. The
inverse track transformation is applied in reverse order. The negative y sign
is the convention selected by the 10,000-frame sign comparison.

## Fit configuration

Fit rows have the form:

```text
parameter start min max status
```

where `status` is `fixed` or `free`. Run a configured fit with:

```bash
process/build/irt \
    --input input.root \
    --output output.root \
    --config process/config/geometry/nominal.conf \
    --fit-config process/config/geometry/fit.conf
```

Use bounded deltas or explicitly constrained parameters when possible, because
track and assembly parameters can be strongly correlated.
