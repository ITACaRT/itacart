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

# Geographic vs Grid Geometry

A cell has an exported geometry and it has an identity, and they are not
interchangeable. This page states the scope of what round-trips, and what a
plausible-looking shortcut costs. It is a declared limitation of this version,
stated rather than implied.

## Export and recovery are an inverse pair. Export and refill are not.

`cells_to_geojson` writes the index into the Feature `id` member and the cell's
vertices into the geometry. `recover_from_geojson` reads the index back out of
the Feature. The pair returns the cell it was handed, and it does so **because
the index travelled with the Feature** rather than being inferred from the
coordinates.

That is the scope of the claim. A recovery handed a Feature whose geometry
belongs to another cell still answers the original index, and a named test does
exactly that.

The geometry is exact in its own right, and separately: the exported vertices,
joined again as straight segments in the sinusoidal plane, rebuild the cell to
within a unit in the last place.

```{figure} ../../../_static/f7/f7_07_figure_7a_round_trip.png
:alt: The exported vertices rebuilding the cell
:width: 100%

`docs/_static/f7/f7_07_figure_7a_round_trip.png`
```

**Taking the exported polygon for a new geographic geometry** — joining its
positions along the ellipsoid and filling the result — is neither half of that
pair. It is reingestion of a rendering, and nothing here promises it returns
the cell.

```{figure} ../../../_static/f9/f9_15_polygon_departs_from_cell.png
:alt: Where a reingested polygon departs from the cell
:width: 100%

`docs/_static/f9/f9_15_polygon_departs_from_cell.png`
```

## What reingestion costs, measured

Row 999 at resolution 1 holds exactly six cells: four in the last column, one
per quadrant, each absorbing the border, and two prime-meridian triangles,
which exist in the eastern quadrants only because the western quadrants have no
meridian column. Six is the enumerated total of that row, not a sample of it.

```{code-cell} python
from itacart import boundary
import itacart
census = []
for quadrant in ("NE", "NW", "SE", "SW"):
    last = boundary.last_lattice_column(quadrant, 999, 10000.0)
    for column in (0, last):
        candidate = f"{quadrant}({column:04d}/0999)"
        try:
            census.append((candidate, boundary.cell_shape(candidate)))
        except itacart.exceptions.ITACaRTError:
            pass
len(census)
```

With orthodromic densification at one kilometre, the reingested outline covers
**33.1111 per cent** of an absorber and **51.6779 per cent** of a triangle.

Both figures are readings of that step rather than properties of the grid.
Refined, they settle near 33.03 and 51.57 per cent, and the declared readings
and the settled ones are pinned together, so neither can be read as the other.

```{figure} ../../../_static/f5/f5_11_round_trip_border.png
:alt: The round trip at the border
:width: 100%

`docs/_static/f5/f5_11_round_trip_border.png`
```

## What each route answers

Handed such an outline:

- {func}`~itacart.polyfill` under `center` **raises** rather than answering an
  empty cover, because the figure that comes back holds no cell centre.
- Under `intersects` the cell is in the answer, along with the others the
  outline touches.
- `recover_from_geojson` returns the cell, because it reads the `id`.

An ordinary cell loses less than a thousandth of its figure on the same route.
The refusal belongs to **this row**, not to the route.

```{figure} ../../../_static/f7/f7_08_border_families_refused.png
:alt: Which border families are refused on reingestion
:width: 100%

`docs/_static/f7/f7_08_border_families_refused.png`
```

## The rule to carry away

If you need the cell back, carry the index — that is what the `id` member is
for, and it is exact everywhere including the border. If you need the shape,
the exported geometry is exact as a shape. What you must not do is treat the
shape as a route back to the identity: it works almost everywhere, and where it
fails it fails by a third of a cell.
