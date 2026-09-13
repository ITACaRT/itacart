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

# Cell Validity & Existence

Three things can be true of an index string, and they are not the same thing.
It can **parse**. It can be **valid**. And the cell it names can **exist**. A
string can pass the first two and fail the third.

```{code-cell} python
from itacart import index as ix
ix.is_valid_index("XX(abc)")
```

Not a quadrant code, so it does not parse.

```{code-cell} python
ix.is_valid_index("NE(0001/0002(A1))")
```

This one parses. It fails anyway, because resolution 2 is an even level and
wants a digit from `1` to `4`; `A1` belongs to the odd alphabet. Validity is
about whether every code is legal *for its own level*.

```{figure} ../../../_static/f2/f2_03_alphabet_acceptance.png
:alt: Which alphabet each level accepts
:width: 100%

`docs/_static/f2/f2_03_alphabet_acceptance.png`
```

## Valid and still nonexistent

```{code-cell} python
ix.is_valid_index("NW(0000/0500)")
```

Valid. Every component is legal, the quadrant is real, the coordinates are in
range. But the cell is not there:

```{code-cell} python
import itacart
from itacart.exceptions import NonExistentCellError
try:
    itacart.cell_to_centroid("NW(0000/0500)")
except NonExistentCellError as exc:
    print(exc)
```

The reason is in the specification rather than in the implementation. Cells
straddling the prime meridian are isosceles triangles, mirrored across the
meridian, and they are indexed in the eastern quadrants only — the western
quadrants have no meridian column, because it would name the same ground twice.
So at resolution 1 a western cell with an X index of zero is structurally
nonexistent. [Meridian and Equator](../boundaries/meridian_and_equator.md) is
where that geometry is worked out.

## Why the package refuses rather than repairs

A nonexistent cell could be silently mapped to its eastern twin. The package
raises instead, and the choice is deliberate: an address that names no cell is
almost always a bug upstream — an arithmetic slip in a traversal, a coordinate
that escaped its domain, a record built by hand. Answering it plausibly hides
the fault and propagates it into stored data.

The same reasoning governs the rest of the error surface. Every exception the
package raises descends from one root, so a pipeline can guard itself without
catching stray failures from elsewhere:

```{code-cell} python
from itacart.exceptions import ITACaRTError
issubclass(NonExistentCellError, ITACaRTError)
```

## The distinctions that matter most

`InvalidIndexError` — the string is malformed. `InvalidRefinementCodeError` — a
code is from the wrong alphabet for its level. `NonExistentCellError` — the
address is well formed and names nothing. `DomainError` — coordinates or a cell
fall outside the addressable domain altogether. `MaxResolutionError` and
`MinResolutionError` — an attempt to refine past 13 or ascend above 0.

Catching the root is right for a pipeline boundary. Catching a specific one is
right where you can actually do something different — retrying a traversal at
a boundary, say, rather than aborting it. The full hierarchy is in
{py:mod}`itacart.exceptions`.
