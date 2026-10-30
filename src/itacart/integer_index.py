"""The integer index: one 64-bit integer per cell, and back.

A compositional index is a string, and a string is the wrong key for a
database column, a sort, a join or a hash table of a hundred million
cells. This module gives every cell an integer, in the manner of the
NUNIQ scheme of the IVOA's Multi-Order Coverage maps for HEALPix: one
number that carries the resolution, sorts in the grid's canonical order
within a resolution, and reaches its parent by arithmetic.

The integer of a cell at resolution ``k`` is ``OFFSET[k] + position``.
The *position* counts the box of every address the grammar can write at
that resolution, in canonical order: the quadrant first, in the order
``NE, NW, SE, SW``; then the resolution-1 pair, column before row as the
siblings of :func:`itacart.compose` are ordered, over 2 005 columns and
1 001 rows; then one refinement code per level, by its place in the
level's alphabet, four at an even level and twenty-five at an odd one::

    position(quadrant)  = q
    position(level 1)   = (q * 2005 + column) * 1001 + row
    position(level k)   = position(level k-1) * r_k + code_k

and ``OFFSET[k]`` is the size of every box above resolution ``k``, so the
fourteen resolutions occupy fourteen consecutive ranges, the four
quadrants of resolution 0 first. The largest integer is 8.4335e18, below
2**63 = 9.2234e18, so every integer fits a signed 64-bit column: the
type databases and columnar formats offer everywhere, where the unsigned
one is not.

The box is wider than the grid. It holds the western zero column, which
no quadrant west of the meridian has; the columns a row lacks past the
border; column 2004, which is only ever a prefix; and the codes the polar
row collapses away. Those positions have integers and name no cell, so
the existence contract runs both ways here: a spelling that names no cell
gets no integer, and an integer that names no cell gets no spelling. Both
refusals are :class:`~itacart.exceptions.NonExistentCellError`.

**Ancestry is lexical, as it is for strings.** :func:`int_parent` divides
the position, which is :func:`itacart.get_parent` in arithmetic: the
parent of ``NE(2004/0000(3))`` is the integer of ``NE(2004/0000)``, a
prefix that names no cell, although the cell's physical parent is
``NE(2003/0000)``. The physical relation is :func:`itacart.get_children`
and lives with the strings.

The names :func:`str_to_int` and :func:`int_to_str` are the ones H3 gives
the same pair of conversions, kept as aliases for the same reason the
other H3 names are.
"""

from __future__ import annotations

import operator
from bisect import bisect_right
from typing import Iterable

from ._existence import require_existing_cells
from .boundary import is_valid_cell
from .constants import MAX_RESOLUTION, QUADRANTS, refinement_alphabet
from .exceptions import (
    InvalidIndexError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)
from .index import is_atomic, iter_cells, split_components

__all__ = [
    "cell_to_int",
    "int_to_cell",
    "int_resolution",
    "int_parent",
    "str_to_int",
    "int_to_str",
]

#: Columns of the resolution-1 box: 0 to 2004, the last one a prefix only.
BOX_COLUMNS = 2005
#: Rows of the resolution-1 box: 0 to 1000, the last one the polar caps.
BOX_ROWS = 1001

#: Positions under one position of the level above, from resolution 0 down.
#: Resolution 1 refines a quadrant into its whole column-by-row box.
RATIOS: tuple[int, ...] = (len(QUADRANTS), BOX_COLUMNS * BOX_ROWS) + tuple(
    len(refinement_alphabet(resolution)) for resolution in range(2, MAX_RESOLUTION + 1)
)


def _sizes() -> tuple[int, ...]:
    sizes = [RATIOS[0]]
    for ratio in RATIOS[1:]:
        sizes.append(sizes[-1] * ratio)
    return tuple(sizes)


#: Positions in the box of each resolution, 0 to 13.
BOX_SIZES: tuple[int, ...] = _sizes()

#: First integer of each resolution, and one past the last of resolution 13.
OFFSETS: tuple[int, ...] = tuple(sum(BOX_SIZES[:k]) for k in range(len(BOX_SIZES) + 1))

#: One past the largest integer. Below ``2**63``, which is the point.
INTEGER_LIMIT: int = OFFSETS[-1]

_QUADRANT_POSITION = {quadrant: position for position, quadrant in enumerate(QUADRANTS)}


def _position(components: list[str]) -> tuple[int, int]:
    """Resolution and box position of one spelling, cell or prefix."""
    position = _QUADRANT_POSITION[components[0]]
    if len(components) > 1:
        column, row = (int(part) for part in components[1].split("/"))
        position = (position * BOX_COLUMNS + column) * BOX_ROWS + row
    for level, code in enumerate(components[2:], start=2):
        position = position * RATIOS[level] + refinement_alphabet(level).index(code)
    return len(components) - 1, position


def _integer(cell: str) -> int:
    resolution, position = _position(split_components(cell))
    return OFFSETS[resolution] + position


def _as_integer(value: object) -> int | None:
    """``value`` as a Python int if it is an integer, else ``None``.

    Anything that declares itself an integer through ``__index__`` is
    taken -- a NumPy ``int64`` read out of an array is the common case --
    and a float is not, even an integral one, because a float that
    reaches an index has usually lost digits on the way. A bool is an int
    in Python and is refused anyway: ``True`` is not a cell, and accepting
    it would decode a filter mask as a quadrant.
    """
    if isinstance(value, bool):
        return None
    try:
        return operator.index(value)  # type: ignore[arg-type]
    except TypeError:
        return None


def _checked(value: object) -> int:
    """An integer of the box, or the refusal of something that is not one."""
    integer = _as_integer(value)
    if integer is None:
        raise InvalidIndexError(
            f"an integer index must be an int, not {type(value).__name__}"
        )
    if not 0 <= integer < INTEGER_LIMIT:
        raise InvalidIndexError(
            f"{integer} is outside the integer index, which runs from 0 to "
            f"{INTEGER_LIMIT - 1}"
        )
    return integer


def _resolution(value: int) -> int:
    return bisect_right(OFFSETS, value) - 1


def _spelling(value: int) -> str:
    """The address an integer of the box writes, whether or not a cell."""
    resolution = _resolution(value)
    position = value - OFFSETS[resolution]
    codes: list[str] = []
    for level in range(resolution, 1, -1):
        position, code = divmod(position, RATIOS[level])
        codes.append(refinement_alphabet(level)[code])
    if resolution == 0:
        return QUADRANTS[position]
    quadrant, pair = divmod(position, BOX_COLUMNS * BOX_ROWS)
    column, row = divmod(pair, BOX_ROWS)
    inner = "".join(f"({code}" for code in reversed(codes))
    return f"{QUADRANTS[quadrant]}({column:04d}/{row:04d}{inner}" + ")" * (
        len(codes) + 1
    )


def _many(value: int | Iterable[int]) -> tuple[list[int], bool]:
    """The integers a call was handed, and whether it was handed one."""
    if _as_integer(value) is not None or isinstance(value, bool):
        return [_checked(value)], True
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise InvalidIndexError(
            f"an integer index must be an int or ints, not {type(value).__name__}"
        )
    return [_checked(item) for item in value], False


def cell_to_int(cell: str) -> int | list[int]:
    """The integer of a cell.

    Args:
        cell: Compositional index string.

    Returns:
        The integer for a single cell, or a list aligned with
        :func:`itacart.decompose` when the index holds several.

    Raises:
        InvalidIndexError: If the index is malformed.
        NonExistentCellError: If any cell of the index names no cell.
            The box has integers for addresses that are not cells, and
            none of them is handed out for a spelling the predicate
            denies.

    Example:
        >>> cell_to_int("NE")
        0
        >>> cell_to_int("NE(0000/0000)")
        4
        >>> cell_to_int("NE(0625/0451(1,2))") == [
        ...     cell_to_int("NE(0625/0451(1))"), cell_to_int("NE(0625/0451(2))")]
        True
    """
    require_existing_cells(cell)
    values = [_integer(atom) for atom in iter_cells(cell)]
    return values[0] if is_atomic(cell) else values


def int_to_cell(value: int | Iterable[int]) -> str | list[str]:
    """The cell an integer names.

    Args:
        value: An integer of the index, or several.

    Returns:
        The atomic index for one integer, or a list of them, aligned
        with the input, for several.

    Raises:
        InvalidIndexError: If a value is not an int, or lies outside the
            integer index.
        NonExistentCellError: If a value lies in the box but its address
            names no cell, as the western zero column, a column past the
            border, the prefix column 2004 and the codes the polar row
            collapses away all do.

    Example:
        >>> int_to_cell(cell_to_int("SE(1400/0374(3(C2)))"))
        'SE(1400/0374(3(C2)))'
    """
    values, single = _many(value)
    cells = []
    for item in values:
        spelling = _spelling(item)
        if not is_valid_cell(spelling):
            raise NonExistentCellError(
                f"{item} is the integer of {spelling!r}, which names no cell"
            )
        cells.append(spelling)
    return cells[0] if single else cells


def int_resolution(value: int | Iterable[int]) -> int | list[int]:
    """The resolution an integer carries, 0 to 13.

    Every integer of the index has one, cell or not: the resolutions
    occupy consecutive ranges, so this is a search over fourteen offsets.

    Raises:
        InvalidIndexError: If a value is not an int, or lies outside the
            integer index.

    Example:
        >>> int_resolution(cell_to_int("NE(0625/0451(3(C2)))"))
        3
    """
    values, single = _many(value)
    resolutions = [_resolution(item) for item in values]
    return resolutions[0] if single else resolutions


def int_parent(
    value: int | Iterable[int], target_res: int | None = None
) -> int | list[int]:
    """The lexical parent of an integer, by division of its position.

    The arithmetic counterpart of :func:`itacart.get_parent`, and lexical
    in the same way: the result is the integer of the prefix one level up
    (or at ``target_res``), and that prefix may name no cell. At the
    domain border it does not, for a child spelled under the next
    column's prefix; :func:`int_to_cell` refuses it there exactly as
    :func:`itacart.is_valid_cell` denies the string.

    Args:
        value: An integer of the index, or several.
        target_res: Resolution of the ancestor wanted. Defaults to one
            level up.

    Returns:
        The parent integer, or a list aligned with the input.

    Raises:
        InvalidIndexError: If a value is not an int, or lies outside the
            integer index.
        MinResolutionError: If a value is a quadrant, which has no parent.
        ResolutionError: If ``target_res`` is not coarser than a value's
            resolution, or is negative.

    Example:
        >>> int_parent(cell_to_int("NE(0625/0451(3))")) == cell_to_int("NE(0625/0451)")
        True
    """
    values, single = _many(value)
    parents = []
    for item in values:
        resolution = _resolution(item)
        if resolution == 0:
            raise MinResolutionError(f"{item} is a quadrant and has no parent")
        wanted = resolution - 1 if target_res is None else target_res
        if not 0 <= wanted < resolution:
            raise ResolutionError(
                f"resolution {wanted} is not coarser than {item}, which is at "
                f"resolution {resolution}"
            )
        position = item - OFFSETS[resolution]
        for level in range(resolution, wanted, -1):
            position //= RATIOS[level]
        parents.append(OFFSETS[wanted] + position)
    return parents[0] if single else parents


str_to_int = cell_to_int
int_to_str = int_to_cell
