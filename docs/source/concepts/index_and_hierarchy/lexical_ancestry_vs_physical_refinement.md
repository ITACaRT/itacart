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

# Lexical Ancestry vs Physical Refinement

Ancestry in ITACaRT is a property of the index string. It is *usually* also a
statement about geometry, and the cases where it is not are declared here
rather than discovered later.

```{code-cell} python
import itacart
parent = "SE(1930/0196)"
child = list(itacart.get_children(parent, flatten=True))[0]
itacart.is_ancestor(parent, child)
```

`is_ancestor` answers by reading the strings: one index is an ancestor of
another when it is a prefix of it. `contains` and `common_ancestor` work the
same way. No geometry is consulted, which is what makes them cheap and what
makes them exact about the thing they actually decide.

## Where the two part company

For an interior cell, lexical ancestry and geometric containment agree: a child
lies inside its parent, and the twenty-five or four children tile it exactly.

At the domain border they do not, and the reason is geometric rather than
incidental. The border is a smooth convex curve on the plane. An absorbing cell
replaces its outer side by the chord of that curve across its own height. A
parent spans one cell height and takes one chord; each child spans half that
height and takes two — and a chord of a convex curve over a shorter span lies
further out than a chord over a longer one.

So a child of an absorbing cell reaches past its parent's effective ring. By
construction, not by rounding, and no enumeration order changes it.

The package declares this as a limitation of this version rather than treating
it as a defect awaiting repair, and the declaration is pinned by measurement.
Across 572 lateral parents — the last column of every seventh row in four
quadrants — every parent is affected and exactly two children of each leave,
which is 1 144 children. The share of a child's own area falling outside its
parent has a median of 0.038142 per cent, a mean of 0.056931, a ninetieth
percentile of 0.104096, and a worst case of 2.411015 per cent at
`SE(1930/0196)` — the cell used above. The excess has a closed form: the
sagitta of the border across the parent's height, times half that height, to
four decimal places.

Every figure in that paragraph is fixed by a named test, and those tests fail
if the absorbing side stops being a per-cell chord.

## What this means in practice

**Hierarchical queries remain correct.** `is_ancestor`, `contains`,
`compact_cells` and `uncompact_cells` are lexical operations with lexical
semantics. They answer questions about the index tree, and the index tree is
exactly what they describe.

**Areal reasoning across levels needs care at the border.** If you are summing
child areas to obtain a parent area, the sum can exceed the parent's effective
figure at absorbing cells. The guard is the same one that governs area
everywhere in this grid:

```{code-cell} python
itacart.is_equal_area_cell(parent), round(itacart.effective_cell_area(parent), 1)
```

`False`, and the effective figure is roughly half the nominal hundred square
kilometres — this cell is absorbing the border, so neither the equal-area
guarantee nor tidy parent-child area arithmetic applies to it.
[Effective vs Nominal Geometry](../geometry/effective_vs_nominal_geometry.md)
is the general treatment.

**The interior is unaffected.** An ordinary cell's children tile it, and the
distinction on this page never surfaces. That is why it is worth stating
explicitly: a property that holds almost everywhere is the kind that gets
assumed everywhere.
