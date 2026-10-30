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

# Area, Perimeter & Compactness

Two shape metrics are defined over ITACaRT cells, and both exist for one
purpose: comparing this grid against other discrete global grid systems on
their own terms. Each is a transcription of a published convention, and each
states where the convention stops being trustworthy. A third measure, at the
end of the page, is of a different kind: how far a position written as a cell
can lie from the point the cell stands for.

```{code-cell} python
from itacart import metrics
cell = "SE(1400/0374)"
metrics.compactness(cell), metrics.normalized_cell_area(cell)
```

## Compactness

The isoperimetric quotient, $4 \pi A / p^2$. It is unitless, a circle scores
one, and less compact shapes fall towards zero. This is equation (1) of section
2.2 of Kmoch, Vasilyev, Virro and Uuemaa (2022) in "Area and shape
distortions in open-source discrete global grid systems".

The point worth knowing is where $A$ and $p$ are measured. They are **planar**,
taken in the oblique ellipsoidal Lambert Azimuthal Equal Area projection
centred on the cell's own centroid — step 3 of that paper's workflow. They are
not geodesic quantities, and they are not measured on the sinusoidal plane
either.

That is a deliberate choice to be comparable rather than merely accurate. A
geodesic perimeter would be a better description of the cell; it would also be
incommensurable with every number the comparison paper publishes, which would
defeat the purpose of computing the metric at all. For a ten-kilometre cell the
difference against the geodesic quantities is second-order in the ratio of cell
size to Earth radius.

The quotient is defined for every cell whose ring does not wrap the whole
parallel. The cells that do are all inside the polar caps: the one that holds
the pole, at every resolution, whose ring runs pole, -180, +180 and so visits
the same point twice, and, from resolution 4 down, many cells around it whose
rings wrap the parallel once or many times over. Everything else is measured,
down to the centimetre and down to the thinnest trapezoid the border cuts:

```{code-cell} python
from itacart.exceptions import GeometryError
sliver = "NE(1984/0090(3(E3)))"
try:
    metrics.compactness("NE(0000/1000)")
except GeometryError as exc:
    refused = type(exc).__name__
refused, round(metrics.compactness(sliver), 4)
```

The trapezoid has a base of 6.1 centimetres at resolution 3 and still has an
area and a perimeter, so it has a quotient. Near the pole the projection runs
on the authalic latitude, and its cosine is formed from the area between the
latitude and the pole rather than from an arcsine, which lost a digit for every
decade of latitude towards the pole and most of a centimetre cell with them.

## Normalized area

A cell's area divided by the nominal area of its resolution — the second metric
of the same section, with one documented substitution.

The published convention divides each cell's area by the **mean area of all
cells at the same resolution**. That mean is not computable here: resolution 1
alone holds roughly four million cells, and every finer level multiplies that.
The denominator used instead is the area the resolution table assigns to the
level.

The substitution is safe, and it is unsafe to leave unstated. It is safe
because the sinusoidal plane is equal-area, so the nominal value is what the
mean converges to. Stating it matters because a reader comparing these figures
against the paper's is entitled to know that one denominator has been swapped
for another, even when the two agree.

```{code-cell} python
[(c, round(metrics.normalized_cell_area(c), 6))
 for c in ("SE(1400/0374)", "NE(0000/0500)", "SE(1930/0196)")]
```

An interior parallelogram and a meridian triangle both score one — they carry
the nominal area exactly. The absorbing cell scores about a half, and that is
the metric reporting a real fact rather than a rounding artefact.
[Effective vs Nominal Geometry](effective_vs_nominal_geometry.md) is what to do
about it.

## Positional uncertainty of a cell

A position recorded as a cell is known to within that cell, and it is read
back at the cell's representative position, its centroid.
{func}`itacart.cell_max_uncertainty` measures the worst case of that
substitution on the ground: the geodesic distance on the WGS84 ellipsoid from
the centroid to the farthest vertex of the cell. It answers in three forms —
the radius in metres, the largest distance between two vertices in metres, and
the radius over the side of the resolution:

```{code-cell} python
import itacart
cell = "NE(0100/0900)"
[round(itacart.cell_max_uncertainty(cell, output), 3)
 for output in ("radius", "diameter", "ratio")]
```

On the sinusoidal plane every regular parallelogram is the same figure, so a
measure taken there would be a constant of the resolution. On the ground it is
not. The cell's shape there is set by $u = |\lambda \sin \varphi|$, with
$\lambda$ in radians from Greenwich, and to first order in the side $l$ the
radius is

$$
l \cdot \max\left(\tfrac{1}{2}\sqrt{1 + u^2},\
\sqrt{(1 - u/2)^2 + 1/4}\right)
$$

which is $1.118\,l$ on the prime meridian, $l/\sqrt{2}$ where $u = 1$ — the
cell is a square there — and $1.65\,l$ as $u$ reaches $\pi$ near the
antimeridian; in Chukotka, past it, $u$ and the radius grow a little more.
Three cells of the same resolution show the range:

```{code-cell} python
[(c, round(itacart.cell_max_uncertainty(c, "ratio"), 4))
 for c in ("NE(0001/0000)", "NE(0100/0900)", "NE(0176/0943)")]
```

A finite cell departs from that law in proportion to its side, and most near
the pole; the function's documentation gives the departures measured, by band
of latitude, down to the rows just before the polar row and past $u = \pi$ in
Chukotka.

The maxima are taken at the vertices, not over the continuous boundary. Away
from the poles the two agree. The cell that holds a pole is where they part
at every resolution: its base is a whole parallel, and its two base vertices,
at −180 and +180 degrees, are one point, so its diameter over the vertices is
half the diameter of its boundary. That is why its diameter comes out shorter
than its radius:

```{code-cell} python
cap = "NE(0000/1000)"
[round(itacart.cell_max_uncertainty(cap, output), 1)
 for output in ("radius", "diameter")]
```

Under a cap, some spellings outside the cap's physical refinement reach
further along their boundary than at their vertices as well, in radius as in
diameter; the function's documentation gives the figures measured.

The number is a property of the tessellation, not the accuracy of a survey: it
says how far the cell's centre can be from a position the cell holds, and
nothing about how well that position was measured.

## Reading these numbers

Both metrics are descriptive. Neither is a quality score, and neither should be
used to decide whether a cell is fit for a computation — that question is
answered by asking the construction directly:

```{code-cell} python
import itacart
[(c, itacart.is_equal_area_cell(c))
 for c in ("SE(1400/0374)", "NE(0000/0500)", "SE(1930/0196)")]
```

Each metric's page in {py:mod}`itacart.metrics` states which convention it
transcribes and the conditions under which the transcription holds. Read it
before publishing a comparison.
