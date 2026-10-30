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

# The Integer Index

One 64-bit integer per cell, and back. A compositional index is a string, and a
string is the wrong key for a database column, a sort, a join or a hash table
of a hundred million cells. The integer index gives every cell a number that
carries its resolution, sorts in the grid's canonical order, and reaches its
parent by arithmetic.

```{code-cell} python
import itacart

cell = itacart.geo_to_cell(-45.8722, -23.2124, 9)
number = itacart.cell_to_int(cell)
cell, number, itacart.int_to_cell(number) == cell
```

## How the number is built

The construction follows the NUNIQ scheme of the IVOA's Multi-Order Coverage
maps. Each resolution has a *box* of positions, and the positions of one
resolution follow those of the one above:

```
position(quadrant) = q                                  NE, NW, SE, SW
position(level 1)  = (q × 2005 + column) × 1001 + row
position(level k)  = position(level k − 1) × r_k + code_k
integer            = OFFSET[k] + position
```

Here `r_k` is 4 at an even resolution and 25 at an odd one, and the code is
the refinement code's place in the alphabet of its level. `OFFSET[k]` is the
size of every box above resolution `k`, so the fourteen resolutions, the
quadrants included, occupy fourteen consecutive ranges:

```{code-cell} python
[itacart.int_resolution(itacart.cell_to_int(c)) for c in ("NE", "NE(0000/0000)", cell)]
```

The largest integer is 8.4335·10¹⁸, below 2⁶³ ≈ 9.2234·10¹⁸. Every integer
fits a signed 64-bit column, the type that databases and columnar formats
offer everywhere, where an unsigned one is not.

```{code-cell} python
from itacart.integer_index import INTEGER_LIMIT

INTEGER_LIMIT, INTEGER_LIMIT < 2**63, round(100 * (1 - INTEGER_LIMIT / 2**63), 2)
```

## Order and ancestry

Within a resolution, the integers sort in the canonical order of
{func}`itacart.compose`. Dividing a position by the ratios walks up its
spelling, so the parent is arithmetic:

```{code-cell} python
itacart.int_to_cell(itacart.int_parent(number)) == itacart.get_parent(cell)
```

**Ancestry is lexical, as it is for strings.** The parent of a border cell's
surplus child is the integer of the prefix it is spelled under, and that prefix
names no cell. The physical parent is a question about the border, and lives
with the strings in {func}`itacart.get_children`:

```{code-cell} python
surplus = "NE(2004/0000(3))"
prefix = itacart.int_parent(itacart.cell_to_int(surplus))
try:
    itacart.int_to_cell(prefix)
except itacart.NonExistentCellError as error:
    print(error)
```

## The box is wider than the grid

The box holds positions that are not cells: the western zero column, which no
quadrant west of the meridian has; the columns a row lacks past the border;
column 2004, which is only ever a prefix; and the codes the polar row
collapses. Those positions have integers and name no cell, so the
[existence contract](cell_validity_and_existence.md) runs both ways. A
spelling that names no cell gets no integer, and an integer that names no cell
gets no spelling:

```{code-cell} python
try:
    itacart.cell_to_int("NW(0000/0110)")
except itacart.NonExistentCellError as error:
    print(error)
```

## What is accepted

Anything that declares itself an integer is accepted, such as a NumPy `int64`
read out of an array. A float is refused, even an integral one, because a float
that reaches an index has usually lost digits on the way. A bool is refused
too: `True` is not a cell. An iterable of integers answers with a list aligned
with it.

```{code-cell} python
import numpy as np

itacart.int_to_cell(np.int64(number)) == cell, itacart.int_resolution([number, 0])
```

`str_to_int` and `int_to_str` are aliases, under the names H3 gives the same
pair of conversions. For whole regions rather than single cells, the integers
become ranges: [Intervals](intervals.md).
