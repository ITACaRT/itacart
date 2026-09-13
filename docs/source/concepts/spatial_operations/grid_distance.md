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

# Grid Distance

Distance measured in cell steps, not in metres.

```{code-cell} python
import itacart
from itacart import topology
a = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
b = itacart.geo_to_cell(-46.6300, -23.5500, resolution=9)
topology.grid_distance(a, b)
```

Three hundred and sixty-two steps at one-metre resolution between two positions
a few hundred metres apart. The number counts lattice moves under a metric, and
the metric is a choice:

```{code-cell} python
(topology.grid_distance(a, b),
 topology.grid_distance(a, b, metric="manhattan"))
```

`chebyshev` allows diagonal moves and so needs fewer of them; `manhattan`
allows only axis moves. Both are exact counts of their own kind of step.

## What it is good for

Anything where the question is really about adjacency at a chosen scale:
whether two parcels are within a few cells of each other, how far a traversal
has to run, how wide a buffer in cells should be. It is integer, it is computed
from the indices alone, and it never touches geometry.

It is also the natural companion to rings and disks — a cell is in
`grid_disk(origin, k)` exactly when its grid distance from the origin is at
most `k` under the same metric.

## What it is not

**It is not a geodesic distance.** Multiplying steps by a cell edge length gives
a figure that looks like metres and is not one: the path is a staircase on the
lattice, not a great-circle arc, and a chebyshev diagonal covers more ground
than an axis step of the same count. If you need metres between two positions,
take the centroids and use the geodesic solutions in {py:mod}`itacart.geodesy`.

**It is not comparable across resolutions.** A distance of ten at resolution 9
and a distance of ten at resolution 11 describe very different separations,
because the step is the cell. Always carry the resolution alongside the number.

**It is not reliable across the domain border.** Where cells are missing,
absorbed or deflected, the lattice is not a clean square grid and step counts
stop meaning what they mean in the interior.
[Absorbing Boundary Cells](../boundaries/absorbing_boundary_cells.md) and
[Meridian and Equator](../boundaries/meridian_and_equator.md) are where that
irregularity is described.

## A working rule

Use grid distance for questions about the grid. Use geodesy for questions about
the Earth. The two agree in spirit and never in units.
