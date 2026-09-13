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
states where the convention stops being trustworthy.

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
