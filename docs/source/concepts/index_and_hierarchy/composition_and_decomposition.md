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

# Composition & Decomposition

A compositional index folds many cells into one string. Decomposing unfolds it;
composing folds it back. The two are inverses, and the order is stable — which
is what everything else in the package leans on.

```{code-cell} python
from itacart import index as ix
region = "NE(0001/0002(1(A1,A2,A3)))"
ix.decompose(region)
```

Three atomic indices, each a complete address in its own right. Folding them
back returns what we started with:

```{code-cell} python
ix.compose(ix.decompose(region)) == region
```

## The order is a contract, not an accident

`decompose` emits cells in a fixed order, and every function that answers one
value per cell answers in that same order.

```{figure} ../../../_static/f2/f2_06_decompose_order.png
:alt: The ordering decompose emits
:width: 100%

`docs/_static/f2/f2_06_decompose_order.png`
```

That alignment is what makes the notation usable in a pipeline. Ask for the
parents of a region and you get one parent per cell, positionally matched:

```{code-cell} python
import itacart
list(zip(ix.decompose(region), itacart.get_parent(region)))
```

Without the guarantee you would have to carry a separate list of which answer
belongs to which cell — which is exactly the unstructured-list problem the
compositional index exists to avoid.

## Counting without unfolding

A region can address a very large number of cells, and you often want the count
rather than the cells:

```{code-cell} python
ix.count_cells(region), ix.count_cells("SE(1400/0374)")
```

`count_cells` reads the structure rather than materialising the members. When
you do want them one at a time, `iter_cells` streams in `decompose` order
instead of building the whole list:

```{code-cell} python
next(iter(ix.iter_cells(region)))
```

## Taking one address apart and putting it back

`split_components` and `join_components` are the atomic-level pair, working on
a single cell's path rather than on a set:

```{code-cell} python
parts = ix.split_components("SE(1400/0374(3(C2(3))))")
parts, ix.join_components(parts)
```

This is the level at which you inspect a resolution, swap a refinement code, or
build an address from parts you already hold. The set-level pair and the
component-level pair do not mix: `decompose` is about how many cells a string
names, `split_components` is about how deep one cell's address goes.

## The parse tree underneath

Both pairs are views of one structure, which `parse` exposes directly:

```{code-cell} python
ix.parse("NE(0001/0002(1(A1,A2)))")
```

```{figure} ../../../_static/f2/f2_05_figure7_parse_tree.png
:alt: The parse tree of a compositional index
:width: 100%

`docs/_static/f2/f2_05_figure7_parse_tree.png`
```

Reaching for `parse` is rarely necessary — the functions above cover the
ordinary cases — but it is there when you need to walk the structure rather
than consume it.
