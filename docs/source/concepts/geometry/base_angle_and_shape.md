---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Base Angle & Shape

On the sinusoidal plane every ITACaRT cell is the same parallelogram, and its
acute angle is exactly forty-five degrees. On the ellipsoid it is not. The base
angle is what that forty-five degrees becomes once the cell is carried to the
curved surface, and it is measured rather than predicted.

```{code-cell} python
from itacart import metrics
[(f"NE({x:04d}/0000)", round(metrics.cell_base_angle(f"NE({x:04d}/0000)"), 6))
 for x in (1, 500, 1000, 2000)]
```

Along the equator the angle stays close to forty-five, drifting by about seven
hundredths of a degree across two thousand columns. Away from it the drift is
not small:

```{code-cell} python
[(f"NE(0100/{y:04d})", round(metrics.cell_base_angle(f"NE(0100/{y:04d})"), 6))
 for y in (0, 100, 300, 600, 900)]
```

## What is being measured, and between what

The angle is taken at the **anchor**, between two lines, and neither choice is
arbitrary.

One arm is the **parallel** — its due-east-and-due-west tangent at the anchor.
A parallel is not a geodesic, and that matters: reading the base as a geodesic
chord instead biases the result by half the convergence of the meridians across
the base. The bias is zero at the equator, where a parallel *is* a geodesic,
and grows with latitude — about 0.023 degree at row 300 and 0.284 at row 900,
where it is the entire difference between 89.426479 and 89.710055.

The other arm is the **geodesic** from the anchor to the corner the shear
displaces. That corner is identifiable in a parallelogram, a triangle and a
trapezoid alike, so the metric is defined for every cell rather than for the
interior only.

## It is a function of lattice position

The base angle is not determined by latitude alone. ITACaRT is
constructed on the ellipsoid, so the angle follows the cell's
position in the ellipsoidal lattice.

Two cells at the same latitude in different columns do not
share an angle:

```{code-cell} python
[(c, round(metrics.cell_base_angle(c), 6))
 for c in ("NE(0100/0374)", "NE(1400/0374)", "SE(1400/0374)")]
```

The two northeastern cells sit on the same row and differ by more
than sixty degrees. Their longitudinal positions in the lattice
are different, so equal latitude does not imply equal shape.

The southeastern cell mirrors its northern counterpart and has the
same angle. This is a consequence of the symmetry of the ITACaRT
construction across the equator.

 `cell_base_angle` therefore evaluates the geometry of the cell
itself rather than treating latitude as a lookup variable.

```{figure} ../../../_static/f4/f4_18_base_crossover.png
:alt: Where the base and height relationship crosses over
:width: 100%

`docs/_static/f4/f4_18_base_crossover.png`
```

```{figure} ../../../_static/fpaper/fpaper_08_angular_distortion.png
:alt: Angular distortion of the parallelogram grid within one quadrant
:width: 100%

Figure 8 of the paper: the distortion near the prime meridian and away
from it, at same quadrant observed (a) in proximity to the prime
meridian within the Italian Peninsula, (b) at mid-latitudes and
mid-longitudes in the East China Sea, and (c) at high latitudes
and longitudes in the Kamchatka Peninsula.
`docs/_static/fpaper/fpaper_08_angular_distortion.png`
```

## Why the shape is allowed to vary

A grid could insist on a constant angle and pay for it in area. ITACaRT insists
on area and pays for it in angle. Shape distortion changes what a cell looks
like; area distortion would change what a cell is worth, and in cadastral work
the second is the one with consequences.

So the lean is a measured property of each cell, the equal-area guarantee is
the invariant, and [Equal Area](../grid_fundamentals/equal_area.md) is where
that ordering is argued.
