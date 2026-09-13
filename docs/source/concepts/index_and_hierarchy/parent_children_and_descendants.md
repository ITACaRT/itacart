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

# Parent, Children & Descendants

Ascent is truncation. Descent is enumeration. The asymmetry between those two
words is the whole content of this page.

```{code-cell} python
import itacart
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
parent = itacart.get_parent(cell)
parent
```

The parent is the index with its last component removed — nothing is computed
from coordinates, and that is what the standard means when it asks that
topological relations be derivable from the identifier.

```{code-cell} python
itacart.get_ancestors(cell)[:3]
```

Coarse to fine, ending just above the cell itself.

## Descent counts rather than multiplies

Going down, the refinement ratio tells you what to expect: four children at an
even level, twenty-five at an odd one.

```{code-cell} python
children = list(itacart.get_children(parent, flatten=True))
len(children)
```

But the ratio is a prediction, not the answer. At the domain border a cell can
hold fewer children than the rule implies, because the border truncates the
parent's own footprint:

```{code-cell} python
len(list(itacart.get_children("SE(1930/0196)", flatten=True)))
```

Three, where the ratio says four. The package obtains the count by enumerating
what actually exists rather than by multiplying, and this is why. A cell whose
outer side has been carried onto the border has been measured to hold between
two and six children. [Absorbing Boundary Cells](../boundaries/absorbing_boundary_cells.md)
is where those cells come from.

## The shape that does not collapse

`get_children` keeps a grouping that its neighbours drop. `get_parent` answers
a scalar for one cell; `get_ancestors` answers a flat chain. `get_children`
answers one list *per input cell*, and it does so even when the input is one
cell:

```{code-cell} python
groups = list(itacart.get_children(parent))
len(groups), len(groups[0])
```

Length one, with the children inside it. That is not a bug and the docstring
states it: the grouping is what preserves positional alignment when the input
is compositional, and it does not disappear just because the input happens to
hold a single cell. Pass `flatten=True` to read the children directly, as every
example on this page does.

## Descendants at depth

`get_children` steps one level by default and accepts a target resolution;
`get_descendants` streams a whole level at once:

```{code-cell} python
itacart.count_cells(itacart.get_parent(parent)), sum(
    1 for _ in itacart.get_descendants(itacart.get_parent(parent), target_res=9))
```

Streaming matters here. Descending several odd levels multiplies by
twenty-five each time, so a list would be the wrong shape long before the
computation is.

## Where the cell sits among its siblings

```{code-cell} python
itacart.child_position(cell)
```

The ordinal is positional within the parent, which is what the neighbour rules
at even and odd resolutions arithmetic on —
[Neighbors & Grid Traversal](../spatial_operations/neighbors_and_grid_traversal.md)
uses it.
