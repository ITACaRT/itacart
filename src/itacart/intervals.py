"""Intervals: a region as sorted ranges of integers, with set operations.

The Multi-Order Coverage maps of the IVOA store a region of the HEALPix
sky as a list of half-open ranges of cell numbers at one depth, and do
union, intersection and difference on those lists by merging them. This
module does the same for ITACaRT, on the integers of
:mod:`itacart.integer_index`.

An interval list is a list of ``(start, stop)`` pairs of integers, sorted,
disjoint and with no two touching, every integer of one resolution: the
*depth* of the list, which the integers carry with them. A cell coarser
than the depth stands for its descendants at the depth, and away from the
domain border those are a single range, because the integers of a cell's
descendants are consecutive. The overlap of two parcels, or of an obstacle
and a flight corridor, is then one intersection of two lists, in time
linear in their lengths.

**The semantics is the region's, and a region in this package is
physical.** A cell whose outer side has been carried onto the domain
border has children spelled under the next column's prefix, outside the
range its own integers span, and the range would also hold addresses past
the border that name no cell. Such a cell is expanded through its physical
children, as :func:`itacart.get_children` gives them, down to the depth.
Everywhere else the range is exact: a cell that does not absorb the border
has exactly its lexical descendants as children, at every depth.

**The depth is chosen, as a coverage map's is, and it is not fixed at 13.**
The border family thins with depth but does not end, so a cell that
absorbs the border has one absorbing descendant per row at every
resolution below it -- a thousand at resolution 7 under a resolution-1
cell, a million at resolution 13. Pinning every list to resolution 13 would
turn one border cell into a million ranges. The default depth is the
finest resolution the index holds, which keeps the cost proportional to
the region's own boundary at its own resolution. Two lists of different
depths are brought to the finer one before they are combined, which is
exact, because a coarser list names whole cells.

Back from intervals to an index, the answer is what
:func:`itacart.compact_cells` writes, so its TreeBlob is the one
:func:`itacart.serialize_to_blob` writes for the region. That holds for
every index whose cells do not overlap; an index that writes a cell and,
under another column's prefix, a child of it, compacts with the child kept
and comes back from its intervals without it.
"""

from __future__ import annotations

from typing import Iterable, Sequence

from ._existence import require_existing_cells
from .boundary import absorbs_border, is_valid_cell
from .constants import MAX_RESOLUTION
from .exceptions import InvalidIndexError, NonExistentCellError, ResolutionError
from .hierarchy import compact_cells, get_children
from .index import compose, iter_cells, split_components
from .integer_index import (
    OFFSETS,
    RATIOS,
    _as_integer,
    _position,
    _resolution,
    _spelling,
)

__all__ = [
    "cells_to_intervals",
    "intervals_to_cells",
    "intervals_union",
    "intervals_intersection",
    "intervals_difference",
    "intervals_contain",
    "intervals_overlap",
    "intervals_count",
]

Interval = tuple[int, int]


def _span(coarse: int, fine: int) -> int:
    """Positions at ``fine`` under one position at ``coarse``."""
    span = 1
    for level in range(coarse + 1, fine + 1):
        span *= RATIOS[level]
    return span


def _merged(intervals: Iterable[Interval]) -> list[Interval]:
    """Sorted, and every overlapping or touching pair joined."""
    out: list[Interval] = []
    for start, stop in sorted(intervals):
        if out and start <= out[-1][1]:
            if stop > out[-1][1]:
                out[-1] = (out[-1][0], stop)
        else:
            out.append((start, stop))
    return out


def _canonical(intervals: Iterable[Sequence[int]]) -> tuple[list[Interval], int | None]:
    """An interval list in canonical form, and its depth.

    Any iterable of pairs is taken and sorted and merged, so a list that
    was concatenated or written by hand is accepted; what is refused is
    anything that is not a list of ranges of one depth.
    """
    if isinstance(intervals, (str, bytes)):
        raise InvalidIndexError("an interval list is a list of (start, stop) pairs")
    pairs: list[Interval] = []
    depth: int | None = None
    for pair in intervals:
        try:
            start, stop = pair
        except (TypeError, ValueError):
            raise InvalidIndexError(
                f"an interval is a (start, stop) pair, not {pair!r}"
            ) from None
        ends = (_as_integer(start), _as_integer(stop))
        if ends[0] is None or ends[1] is None:
            raise InvalidIndexError(f"interval ends must be ints: {pair!r}")
        start, stop = ends[0], ends[1]
        if not start < stop:
            raise InvalidIndexError(f"an interval must be non-empty: {pair!r}")
        if start < 0 or stop > OFFSETS[-1]:
            raise InvalidIndexError(f"{pair!r} lies outside the integer index")
        here = _resolution(start)
        if stop > OFFSETS[here + 1]:
            raise InvalidIndexError(f"{pair!r} crosses from one resolution to the next")
        if depth is None:
            depth = here
        elif here != depth:
            raise InvalidIndexError(
                f"an interval list holds one resolution, and this one holds "
                f"{depth} and {here}"
            )
        pairs.append((start, stop))
    return _merged(pairs), depth


def _expand(cell: str, depth: int, out: list[Interval]) -> None:
    resolution, position = _position(split_components(cell))
    if resolution == depth:
        value = OFFSETS[depth] + position
        out.append((value, value + 1))
        return
    if resolution == 0:
        raise ResolutionError(
            f"{cell!r} is a quadrant; its region at a finer resolution is the "
            "resolution-1 lattice, which is not enumerated"
        )
    if not absorbs_border(cell):
        span = _span(resolution, depth)
        start = OFFSETS[depth] + position * span
        out.append((start, start + span))
        return
    for child in get_children(cell, flatten=True):
        _expand(str(child), depth, out)


def cells_to_intervals(index: str, resolution: int | None = None) -> list[Interval]:
    """The region an index names, as an interval list.

    Args:
        index: Compositional index string.
        resolution: Depth of the list. Defaults to the finest resolution
            the index holds. A cell coarser than it stands for its
            physical descendants at it.

    Returns:
        The canonical interval list: sorted, disjoint, no two touching.

    Raises:
        InvalidIndexError: If the index is malformed.
        NonExistentCellError: If any cell of the index names no cell.
        ResolutionError: If ``resolution`` is outside 0 to 13, or coarser
            than a cell of the index, or if a quadrant is to be expanded:
            its descendants are the resolution-1 lattice, which the
            hierarchy does not enumerate either.

    Example:
        >>> from itacart import cell_to_int
        >>> cells_to_intervals("NE(0625/0451(1))", resolution=3) == [
        ...     (cell_to_int("NE(0625/0451(1(A1)))"),
        ...      cell_to_int("NE(0625/0451(1(E5)))") + 1)]
        True
    """
    require_existing_cells(index)
    cells = list(iter_cells(index))
    finest = max(len(split_components(cell)) - 1 for cell in cells)
    depth = finest if resolution is None else resolution
    if not 0 <= depth <= MAX_RESOLUTION:
        raise ResolutionError(f"resolution {depth} is outside 0 to {MAX_RESOLUTION}")
    if depth < finest:
        raise ResolutionError(
            f"resolution {depth} is coarser than a cell of the index, which "
            f"reaches resolution {finest}"
        )
    out: list[Interval] = []
    for cell in cells:
        _expand(cell, depth, out)
    return _merged(out)


def _blocks(start: int, stop: int, level: int, depth: int, out: list[Interval]) -> None:
    """The largest aligned blocks of ``[start, stop)``, coarse to fine.

    Positions are at ``depth``; a block at ``level`` spans the positions of
    one cell of that level. Each block found is appended as ``(level,
    position)``. The two remainders on either side of the blocks are each
    shorter than one block, so the recursion visits at most a few hundred
    blocks per interval however long it is.
    """
    if start >= stop:
        return
    span = _span(level, depth)
    first = -(-start // span)
    last = stop // span
    if first >= last:
        _blocks(start, stop, level + 1, depth, out)
        return
    _blocks(start, first * span, level + 1, depth, out)
    out.extend((level, position) for position in range(first, last))
    _blocks(last * span, stop, level + 1, depth, out)


def _cells_of_block(level: int, position: int, depth: int, out: list[str]) -> None:
    """The cells a wholly covered block stands for, or a refusal.

    A block whose root is a cell that does not absorb the border is that
    cell. A root that absorbs it is replaced by its lexical children,
    every one of which must be a cell: a range covering a child that names
    no cell covers an address past the border. A root that names no cell
    is refused for the same reason, which holds because a prefix that
    names no cell always holds such an address below it.
    """
    spelling = _spelling(OFFSETS[level] + position)
    if not is_valid_cell(spelling):
        raise NonExistentCellError(
            f"the intervals cover {spelling!r}, which names no cell"
        )
    if level == 0 and depth > 0:
        raise NonExistentCellError(
            f"the intervals cover the whole box of quadrant {spelling!r}, which "
            "holds addresses that name no cell"
        )
    if level == depth or not absorbs_border(spelling):
        out.append(spelling)
        return
    ratio = RATIOS[level + 1]
    for code in range(ratio):
        _cells_of_block(level + 1, position * ratio + code, depth, out)


def intervals_to_cells(intervals: Iterable[Sequence[int]]) -> str:
    """The index of the region an interval list covers.

    Returns:
        The compositional index :func:`itacart.compact_cells` writes for
        the region, so its TreeBlob is the region's TreeBlob.

    Raises:
        InvalidIndexError: If the list is empty, or is not a list of
            ranges of one resolution.
        NonExistentCellError: If a range covers an address that names no
            cell: the western zero column, a column past the border, the
            prefix column 2004, or the codes the polar row collapses away.

    Example:
        >>> intervals_to_cells(cells_to_intervals("NE(0625/0451(1,2,3,4))"))
        'NE(0625/0451)'
    """
    canonical, depth = _canonical(intervals)
    if depth is None:
        raise InvalidIndexError("an empty interval list names no cell")
    blocks: list[Interval] = []
    for start, stop in canonical:
        _blocks(start - OFFSETS[depth], stop - OFFSETS[depth], 0, depth, blocks)
    cells: list[str] = []
    for level, position in blocks:
        _cells_of_block(level, position, depth, cells)
    return compact_cells(compose(cells))


def _aligned(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> tuple[list[Interval], list[Interval]]:
    """Both lists in canonical form, at the finer of their two depths."""
    one, one_depth = _canonical(first)
    two, two_depth = _canonical(second)
    if one_depth is None or two_depth is None or one_depth == two_depth:
        return one, two
    if one_depth < two_depth:
        return cells_to_intervals(intervals_to_cells(one), two_depth), two
    return one, cells_to_intervals(intervals_to_cells(two), one_depth)


def intervals_union(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> list[Interval]:
    """Every position in either list.

    Lists of different depths are brought to the finer depth first.

    Raises:
        InvalidIndexError: If either is not a list of ranges of one
            resolution.
        NonExistentCellError: If the coarser list covers an address that
            names no cell, which its refinement finds.

    Example:
        >>> a = cells_to_intervals("NE(0625/0451(1))")
        >>> b = cells_to_intervals("NE(0625/0451(2))")
        >>> intervals_to_cells(intervals_union(a, b))
        'NE(0625/0451(1,2))'
    """
    one, two = _aligned(first, second)
    return _merged(one + two)


def intervals_intersection(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> list[Interval]:
    """Every position in both lists, by one merge of the two.

    Raises:
        InvalidIndexError: If either is not a list of ranges of one
            resolution.
        NonExistentCellError: If the coarser list covers an address that
            names no cell, which its refinement finds.

    Example:
        >>> a = cells_to_intervals("NE(0625/0451(1,2))")
        >>> b = cells_to_intervals("NE(0625/0451(2,3))")
        >>> intervals_to_cells(intervals_intersection(a, b))
        'NE(0625/0451(2))'
    """
    one, two = _aligned(first, second)
    out: list[Interval] = []
    i = j = 0
    while i < len(one) and j < len(two):
        start = max(one[i][0], two[j][0])
        stop = min(one[i][1], two[j][1])
        if start < stop:
            out.append((start, stop))
        if one[i][1] < two[j][1]:
            i += 1
        else:
            j += 1
    return out


def intervals_difference(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> list[Interval]:
    """Every position of the first list that is not in the second.

    Raises:
        InvalidIndexError: If either is not a list of ranges of one
            resolution.
        NonExistentCellError: If the coarser list covers an address that
            names no cell, which its refinement finds.

    Example:
        >>> a = cells_to_intervals("NE(0625/0451)")
        >>> b = cells_to_intervals("NE(0625/0451(4))")
        >>> intervals_to_cells(intervals_difference(a, b))
        'NE(0625/0451(1,2,3))'
    """
    one, two = _aligned(first, second)
    out: list[Interval] = []
    j = 0
    for start, stop in one:
        while j < len(two) and two[j][1] <= start:
            j += 1
        cursor = start
        k = j
        while k < len(two) and two[k][0] < stop:
            if two[k][0] > cursor:
                out.append((cursor, two[k][0]))
            cursor = max(cursor, two[k][1])
            k += 1
        if cursor < stop:
            out.append((cursor, stop))
    return out


def intervals_contain(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> bool:
    """Whether every position of the second list is in the first.

    Example:
        >>> intervals_contain(cells_to_intervals("NE(0625/0451)"),
        ...                   cells_to_intervals("NE(0625/0451(3(C2)))"))
        True
    """
    one, two = _aligned(first, second)
    return not intervals_difference(two, one)


def intervals_overlap(
    first: Iterable[Sequence[int]], second: Iterable[Sequence[int]]
) -> bool:
    """Whether the two lists share a position.

    Example:
        >>> intervals_overlap(cells_to_intervals("NE(0625/0451(1))"),
        ...                   cells_to_intervals("NE(0625/0451(2))"))
        False
    """
    return bool(intervals_intersection(first, second))


def intervals_count(intervals: Iterable[Sequence[int]]) -> int:
    """How many cells of the list's depth the list covers.

    Example:
        >>> intervals_count(cells_to_intervals("NE(0625/0451(3))", resolution=4))
        100
    """
    canonical, _ = _canonical(intervals)
    return sum(stop - start for start, stop in canonical)
