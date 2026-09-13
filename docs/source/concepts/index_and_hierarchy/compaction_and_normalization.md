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

# Compaction & Normalization

Two different tidyings, often confused. Compaction changes which cells a string
names — coarser ones, covering the same ground. Normalization changes only how
the same cells are written.

```{code-cell} python
import itacart
from itacart import index as ix

parent = "NE(0500/0100(1))"
siblings = ix.compose(list(itacart.get_children(parent, flatten=True)))
ix.count_cells(siblings)
```

Twenty-five siblings that exhaust their parent. Compaction replaces them with
it:

```{code-cell} python
itacart.compact_cells(siblings)
```

One cell where there were twenty-five, covering identical ground. The
replacement happens only when the sibling set is exhaustive — a partial set
stays as it is, because folding it would claim ground that was not in the
input.

Compaction recurses: a parent produced this way can itself join a complete set
one level up. Expanding is the inverse, to a uniform resolution:

```{code-cell} python
sum(1 for _ in itacart.uncompact_cells(parent, target_res=3))
```

## Why compaction is worth doing

A parcel covered at a fine resolution can run to a great many cells, and most
of its interior is solid. Compaction leaves the fine cells along the boundary
and folds the interior into as few coarse cells as it can, which is usually an
order-of-magnitude reduction in what has to be stored, hashed or transmitted.
[Parcel Representation](../cadastral_workflows/parcel_representation.md) covers
the trade-off in cadastral terms.

The saving is safe precisely because compaction is lexical: it manipulates the
index tree and does not consult geometry, so it cannot introduce a geometric
error. What it can do at the domain border is inherit one —
[Lexical Ancestry vs Physical Refinement](lexical_ancestry_vs_physical_refinement.md)
is where that is stated.

## Normalization writes the same cells one way

```{code-cell} python
ix.normalize("NE(0001/0002(1(A2,A1,A1)))")
```

Duplicates collapse, siblings are ordered, and the result is canonical: two
strings naming the same set normalize to the same text.

```{figure} ../../../_static/f2/f2_07_normalisation_rewrites.png
:alt: The rewrites normalization performs
:width: 100%

`docs/_static/f2/f2_07_normalisation_rewrites.png`
```

That property is what makes an index comparable and hashable. Two records of
the same parcel, written by different hands in different orders, agree only if
there is a canonical form to agree in — which is the whole basis of using an
index as a key. [Canonical Geometry](../cadastral_workflows/canonical_geometry.md)
carries the idea through to the geometry itself.

## Which one you want

Use **normalization** before comparing, hashing or storing an index. Use
**compaction** when you want fewer cells for the same ground. They compose:
compact first, then normalize the result, and you have the smallest canonical
spelling of a region.
