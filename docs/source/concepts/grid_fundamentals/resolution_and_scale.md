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

# Resolution & Scale

The resolution ladder is decimal. That is the reason a cell count can be
read directly as an area.

```{code-cell} python
import itacart
table = itacart.resolution_table()
len(table)
```

Fourteen entries. Resolution 0 addresses the quadrant and carries no metric
size; the sized levels run from resolution 1 down to resolution 13.

```{code-cell} python
[(row["resolution"], row["cell_size_m"], row["cell_area_m2"])
 for row in table if row["resolution"] in (1, 3, 5, 9, 13)]
```

Ten kilometres, one kilometre, one hundred metres, one metre, one centimetre —
and the areas are the squares of those: one hundred square kilometres down to
one square centimetre. Every one is a round metric quantity.

## Why the refinement alternates

Getting from ten to one in metric takes two steps, not one: ten halves to five,
and five divides by five to one. A grid refines by the same two steps, applied
to each axis.

```{code-cell} python
itacart.constants.REFINEMENT_RATIO[2:10]
```

```{figure} ../../../_static/fpaper/fpaper_02_hierarchical_refinement.png
:alt: The 1-to-25 and 1-to-4 hierarchical refinement
:width: 70%

Figure 2 of the paper: the two refinements, side by side.
`docs/_static/fpaper/fpaper_02_hierarchical_refinement.png`
```

Halving each axis gives four children; dividing each axis by five gives
twenty-five. Alternating them walks the decimal ladder exactly, which a single
constant ratio cannot do — quadtrees give powers of two and never land on a
round metric size, and a uniform 1-to-25 overshoots.

```{figure} ../../../_static/f2/f2_01_resolution_alternation.png
:alt: Alternation between the two refinement ratios across resolutions
:width: 100%

The alternation across the ladder.
`docs/_static/f2/f2_01_resolution_alternation.png`
```

```{figure} ../../../_static/f2/f2_02_refinement_grids.png
:alt: The 1-to-4 and 1-to-25 refinement grids
:width: 100%

The two refinement grids: four children on even levels, twenty-five on odd.
`docs/_static/f2/f2_02_refinement_grids.png`
```

The two ratios are also why the index alphabet changes between levels. A level
that refines into four needs four codes; a level that refines into twenty-five
needs twenty-five, which is where the `A1`–`E5` alphabet comes from.
[Index Structure](../index_and_hierarchy/index_structure.md) reads an address
against these alphabets.

## Choosing a resolution

Cadastral work spans a wide range of scales — the paper takes roughly 1:500 for
urban plots out to 1:10 000 for large rural regions as the range a cadastral
hierarchy has to cover. The table carries a suggested scale for each
resolution, separately for visualisation and for analysis, because the
resolution at which a grid is legible is not the resolution at which it is
trustworthy to measure:

```{code-cell} python
[(row["resolution"], row["visualization_scale"], row["analysis_scale"])
 for row in table if row["resolution"] in (1, 9, 13)]
```

The analysis scale is consistently finer than the visualisation scale — a cell
you can see is coarser than a cell you should compute with.
[Choosing Resolution from Mapping Scale](../cadastral_workflows/choosing_resolution_from_mapping_scale.md)
turns this into a procedure for a real parcel.

## The ceiling

Resolution 13 is the finest level. Asking for children beyond it is refused
rather than approximated:

```{code-cell} python
from itacart.exceptions import MaxResolutionError
finest = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
try:
    list(itacart.get_children(finest, flatten=True))
except MaxResolutionError as exc:
    print(exc)
```

One centimetre is below the precision any GNSS-referenced cadastral survey
delivers, so the ceiling is not a limitation in practice — it is where the grid
stops promising more than the measurement can support.
