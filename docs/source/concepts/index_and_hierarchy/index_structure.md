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

# Index Structure

An index is a string, and the string is the address at every resolution at
once. Reading one left to right is descending the hierarchy.

```{code-cell} python
from itacart import index as ix
ix.split_components("SE(1400/0374(3(C2(3))))")
```

Five components, one per level.

```{figure} ../../../_static/f2/f2_04_index_anatomy.png
:alt: The anatomy of a compositional index
:width: 100%

`docs/_static/f2/f2_04_index_anatomy.png`
```

```{figure} ../../../_static/fpaper/fpaper_03_indexing_by_resolution.png
:alt: Indexing at resolution 0, resolution 1, even and odd resolutions
:width: 100%

Figure 3 of the paper: the four kinds of component, in order —
quadrant, Cartesian pair, the digits 1 to 4, and the A1–E5 alphabet.
`docs/_static/fpaper/fpaper_03_indexing_by_resolution.png`
```

## The four kinds of component

**The quadrant**, a two-letter code: `NE`, `NW`, `SE`, `SW`. The globe is
divided by the prime meridian and the equator, and this is resolution 0.

**The base cell**, a pair of integers `XXXX/YYYY` relative to the quadrant's
origin. This is resolution 1, the ten-kilometre level, and it is addressed by
Cartesian coordinates in the projection — which is the point. A surveyor
reading `1400/0374` is reading a column and a row.

On the WGS84 ellipsoid a quadrant holds roughly two thousand base cells along
the equator and a thousand along the central meridian. Because the quadrants
mirror one another, the coordinates are never negative.

**Even-resolution codes**, a single digit `1` to `4`. An even level refines its
parent into four, in a two-by-two arrangement.

**Odd-resolution codes**, a letter-digit pair from `A1` to `E5`. An odd level
refines its parent into twenty-five, in a five-by-five arrangement.

Which alphabet is legal at a given position is fixed by that position, not
chosen. Offering a digit where the level wants a pair is a malformed index, not
an alternative spelling:

```{code-cell} python
ix.is_valid_index("NE(0001/0002(1))"), ix.is_valid_index("NE(0001/0002(A1))")
```

Resolution 2 is even, so `1` is legal there and `A1` is not.

## The punctuation

A parenthesis descends one level. A comma separates siblings at the same level.
Those two marks are what make the notation compositional — one string can hold
one cell or a whole feature:

```{code-cell} python
ix.is_atomic("SE(1400/0374(3(C2(3))))"), ix.is_atomic("NE(0001/0002(1(A1,A2,A3)))")
```

The second addresses three cells. Every function in the package accepts either
form; [Composition & Decomposition](composition_and_decomposition.md) covers
what comes back.

## Reading the quadrant and the base cell out

Two shortcuts save parsing the string yourself:

```{code-cell} python
ix.quadrant_of("SE(1400/0374(3(C2(3))))"), ix.base_cell_of("SE(1400/0374(3(C2(3))))")
```

## Why the string and not a number

Most grids assign a cell an atomic identifier — an integer from a space-filling
curve, or a bit-packed path. Those are compact and opaque: the relationship
between two cells is recoverable only by decoding.

ITACaRT keeps the hierarchy visible in the text. A parent is a prefix. A
sibling differs in one component. A feature is a set of shared prefixes with
branches. That readability is the usability criterion the grid was designed
against, and it is also what lets a cadastral record carry its own geometry in
a form a person can audit.
