<p align="center">
  <img src="https://raw.githubusercontent.com/itacart/itacart/main/docs/_static/logo.png"
       alt="ITACaRT" width="520">
</p>

**ITA Cadastral Ellipsoidal Reference Tessellation** — an equal-area parallelogram
Discrete Global Grid System (DGGS) for terrestrial cadastral mapping, tessellated
directly on the WGS84 ellipsoid.

[![CI](https://github.com/itacart/itacart/actions/workflows/ci.yml/badge.svg)](https://github.com/itacart/itacart/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/itacart/itacart/branch/main/graph/badge.svg)](https://codecov.io/gh/itacart/itacart)
[![PyPI](https://img.shields.io/pypi/v/itacart)](https://pypi.org/project/itacart/)
[![Python](https://img.shields.io/pypi/pyversions/itacart)](https://pypi.org/project/itacart/)
[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-green.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://img.shields.io/badge/DOI-10.14393%2Frbcv77n0a--79281-blue)](https://doi.org/10.14393/rbcv77n0a-79281)

> **Status:** alpha. The public API is defined and documented; implementations are
> landing module by module.

## What makes it different

Most DGGS project a polyhedron onto a sphere. ITACaRT tessellates the WGS84
ellipsoid directly, trading computational convenience for geodetic fidelity —
the trade that cadastral work actually needs, since parcel area carries legal
and fiscal weight.

- **Equal-area cells** everywhere except at controlled antemeridian exceptions
- **Decimal hierarchy** — 14 levels from 10 km to 1 cm, with whole-number metric
  areas (1 km², 100 m², 1 cm²) that make one token equal one unit of area
- **Cartesian-like addressing** so surveyors can work with it without
  relearning their intuitions
- **Compositional index** that represents a whole vector feature in one string
  instead of an unstructured list of atomic identifiers

## Install

```bash
pip install itacart
```

Optional extras:

```bash
pip install "itacart[geo]"       # GeoDataFrame export
pip install "itacart[parallel]"  # parallel cell filling
```

`shapely` is the only hard runtime dependency. The ellipsoidal sinusoidal
projection is computed directly from the paper's equations, so PROJ is never
called at runtime.

## Quick start

```python
import itacart

# Address a position at 1 cm resolution
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)

# Inverse geometry
lon, lat = itacart.cell_to_centroid(cell)
ring = itacart.cell_to_boundary(cell, close=True)

# Hierarchy
parent = itacart.get_parent(cell)
children = list(itacart.get_children(parent, flatten=True))   # 25 at odd resolutions

# Topology
neighbours = itacart.grid_disk(cell, k_distance=1)

# Vector features
from shapely.geometry import Polygon
parcel = Polygon([...])
index = itacart.polyfill(parcel, resolution=13, compact=True)
area_m2 = itacart.count_internal_cells(parcel, 13) * itacart.nominal_cell_area(13)
```

## Compositional indices

A single string addresses one cell or a whole region:

```
SE(1400/0374(3(C2(3))))              one cell
...4(C1,C2)                          two siblings under one parent
```

Parentheses descend a level; commas separate siblings. Every function accepts
this form. Functions returning one value per cell return a list aligned with
`decompose()` order:

```python
region = "NE(0001/0002(1(A1,A2,A3)))"
itacart.count_cells(region)          # 3
itacart.get_parent(region)           # 3 entries, positionally aligned
itacart.grid_disk(region)            # 3 lists, one per input cell
```

Functions that answer one value per cell collapse that alignment for a single
cell: `get_parent` of one cell is a string, not a list of one. `get_children`
is the exception and keeps its grouping either way, so
`list(get_children(cell))` has length one and the children are inside it — pass
`flatten=True` to read them directly.

## Boundary behaviour

Three departures from the uniform grid are part of the specification, not
edge cases to be handled later:

| Condition | Cell shape | Equal area |
|---|---|---|
| Interior | Parallelogram | Yes |
| Prime meridian | Isosceles triangle | Yes |
| Antemeridian / extension edge | Trapezoid | **No** |

Because trapezoidal cells break the area guarantee, area queries distinguish
nominal from effective:

```python
itacart.nominal_cell_area(13)        # 0.0001 m² — the Table 1 value
itacart.effective_cell_area(cell)    # actual area, clipping accounted for
itacart.is_equal_area_cell(cell)     # gate area-sensitive computations on this
```

Two extension zones carry the eastern quadrants across 180° over inhabited
land — Fiji and the Chukotka/Wrangel group. Antarctica is excluded.

```python
itacart.extension_zone(cell)         # "FIJI" | "CHUKOTKA" | None
```

## Resolutions

| Res | Edge | Area | Refinement | Codes | Visualization | Analysis |
|----:|------|------|-----------:|-------|---------------|----------|
| 1 | 10 km | 100 km² | base | `0000/0000` | 1:100 000 000 | 1:20 000 000 |
| 3 | 1 km | 1 km² | 1→25 | `A1`–`E5` | 1:10 000 000 | 1:2 000 000 |
| 9 | 1 m | 1 m² | 1→25 | `A1`–`E5` | 1:10 000 | 1:2 000 |
| 13 | 1 cm | 1 cm² | 1→25 | `A1`–`E5` | 1:100 | 1:20 |

Even resolutions refine 1→4, odd resolutions 1→25. Full table via
`itacart.resolution_table()`.

## OGC conformance

ITACaRT meets DGGS Core (Topic 21 / ISO 19170-1) in full and EAERS partially.
The divergences are deliberate: direct surface tessellation rules out the
polyhedral interface, and the antemeridian trapezoids are accepted so the grid
covers inhabited land. `itacart.conformance()` reports the full table, and the
conformance test suite verifies it in CI.

## Citing

```bibtex
@article{silva2025itacart,
  author  = {Silva, Israel Nunes and Dietzsch, Gabriel
             and Shiguemori, Elcio Hideiti},
  title   = {{ITACaRT}: An Equal-Area Parallelogram Discrete Global Grid System
             for Terrestrial Cadastral Mapping---Designed for Usability and
             Blockchain Integration},
  journal = {Revista Brasileira de Cartografia},
  volume  = {77},
  year    = {2025},
  doi     = {10.14393/rbcv77n0a-79281}
}
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).

---

Release history in [CHANGELOG.md](CHANGELOG.md). How to build, test and contribute in [CONTRIBUTING.md](CONTRIBUTING.md).
