"""Intervals, measured against the hierarchy that already knows the answer.

The oracle is the package itself. A cell's interval list at a depth is the
set of integers of :func:`itacart.uncompact_cells` at that depth; the way
back is :func:`itacart.compact_cells`; and a set operation on two lists is
the same operation on those two sets of integers. Each claim is checked
over the whole border family, the caps and their descendants, and regions
the package produces, because the border and the caps are where a range
of integers stops being the region a cell names.
"""

from __future__ import annotations

import pytest
from shapely.geometry import Polygon

import itacart
from itacart import integer_index as ii
from itacart import intervals as iv
from itacart.constants import CELL_SIZE_M
from itacart.exceptions import InvalidIndexError, NonExistentCellError, ResolutionError
from tests.unit.test_totality import CORPUS as TOTALITY_CORPUS

QUADRANTS = ("NE", "NW", "SE", "SW")
SIDE_1 = CELL_SIZE_M[1]
assert SIDE_1 is not None


def _positions(intervals: list[tuple[int, int]]) -> list[int]:
    return [value for start, stop in intervals for value in range(start, stop)]


def _oracle(index: str, depth: int) -> list[int]:
    """The integers of a region's physical descendants at ``depth``."""
    integers = [
        ii.cell_to_int(cell)
        for atom in itacart.decompose(index)
        for cell in itacart.uncompact_cells(atom, depth)
    ]
    assert all(isinstance(value, int) for value in integers)
    return sorted(set(integers))  # type: ignore[arg-type]


def _border_family() -> list[str]:
    """The last column of every row below the caps, in four quadrants."""
    cells = []
    for quadrant in QUADRANTS:
        for row in range(1000):
            column = itacart.last_lattice_column(quadrant, row, SIDE_1)
            cells.append(f"{quadrant}({column:04d}/{row:04d})")
    return cells


BORDER = _border_family()


def _square(lon: float, lat: float, side: float) -> Polygon:
    return Polygon(
        [(lon, lat), (lon + side, lat), (lon + side, lat + side), (lon, lat + side)]
    )


# --------------------------------------------------------------------------
# One cell, one range, away from the border
# --------------------------------------------------------------------------


def test_a_cell_away_from_the_border_is_one_range_of_its_descendants() -> None:
    """Consecutive integers, exactly as many as the descendants.

    Interior cells, a meridian triangle and the cells of both sides of the
    meridian, at depths one to three levels below them.
    """
    for cell in (
        "NE(0625/0451)",
        "SW(0007/0008(4(C3)))",
        "NE(0000/0300)",
        "NW(0001/0300(2))",
        "SE(1400/0374(3(C2(3))))",
    ):
        assert itacart.absorbs_border(cell) is False
        resolution = itacart.get_resolution(cell)
        for depth in range(resolution, min(resolution + 3, 13) + 1):
            intervals = iv.cells_to_intervals(cell, depth)
            assert len(intervals) == 1
            assert _positions(intervals) == _oracle(cell, depth)


def test_the_whole_border_family_agrees_with_the_hierarchy() -> None:
    """Four thousand cells that absorb the border, at the next depth.

    Each is expanded through its physical children, some of them spelled
    under the next column's prefix, and the integers are those of
    :func:`itacart.uncompact_cells`. Coming back gives the cell itself,
    which is what :func:`itacart.compact_cells` writes for it.
    """
    assert len(BORDER) == 4000
    split = 0
    for cell in BORDER:
        intervals = iv.cells_to_intervals(cell, 2)
        assert _positions(intervals) == _oracle(cell, 2), cell
        assert iv.intervals_to_cells(intervals) == cell
        split += len(intervals) > 1
    assert split > 0


def test_the_border_family_agrees_two_levels_down() -> None:
    """Every tenth row, at depth 3, where the family is twenty-five wide."""
    for cell in BORDER[::10]:
        intervals = iv.cells_to_intervals(cell, 3)
        assert _positions(intervals) == _oracle(cell, 3), cell
        assert iv.intervals_to_cells(intervals) == cell


def test_the_caps_and_their_descendants_agree_with_the_hierarchy() -> None:
    """The caps absorb, and their descendants have collapsed codes."""
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        for depth in (1, 2, 3, 4, 5):
            intervals = iv.cells_to_intervals(cap, depth)
            assert _positions(intervals) == _oracle(cap, depth)
            assert iv.intervals_to_cells(intervals) == cap
        for cell in itacart.get_descendants(cap, 4):
            intervals = iv.cells_to_intervals(cell, 5)
            assert _positions(intervals) == _oracle(cell, 5)
            assert iv.intervals_to_cells(intervals) == cell


# --------------------------------------------------------------------------
# The way back
# --------------------------------------------------------------------------


def _round_trip_corpus() -> list[str]:
    entries = [
        entry for family in sorted(TOTALITY_CORPUS) for entry in TOTALITY_CORPUS[family]
    ]
    entries.append(
        "NW(0625/0451(1(E1(3(B2(4(A2,B2,B3,B4,C2,C3,C4,C5,D1,D2,D3,D4,D5,E1,E2,E3,"
        "E4,E5)),B3(3(C1,D1,D2,E1,E2,E3,E4)),C2(1(A5,B5,C5,D4,D5,E4,E5),2,3(A4,A5,"
        "B3,B4,B5,C3,C4,C5,D2,D3,D4,D5,E2,E3,E4,E5),4),C3(1,2(A1,B1,B2,B3,C1,C2,C3,"
        "C4,D1,D2,D3,D4,E1,E2,E3,E4,E5),3,4),D2(1(A1,A2,A3,A4,A5,B2,B3,B4,B5,C3,C4,"
        "C5,D5),2,4(A3,A4,A5,B5)),D3(1,2,3(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,C2,C3,C4,"
        "C5,D4,D5),4),D4(1(D1,E1),3(A1,B1,C1,D1,E1)))))))"
    )
    entries.append(
        itacart.compose(list(itacart.get_children("NE(2003/0000)", flatten=True)))
    )
    entries.append(itacart.polyfill(_square(-0.004, 10.0, 0.008), 6, compact=True))
    return entries


def test_the_way_back_is_compact_cells_and_the_regions_treeblob() -> None:
    """The index and the bytes, for regions the package produces.

    The totality corpus, the paper's Figure 7, the five children of a
    border cell and a compacted fill across the meridian. Each comes back
    as :func:`itacart.compact_cells` writes it, and its TreeBlob is the
    one :func:`itacart.serialize_to_blob` writes for the original.
    """
    corpus = _round_trip_corpus()
    assert len(corpus) == 35 + 3
    assert all(
        itacart.is_valid_cell(cell) for x in corpus for cell in itacart.decompose(x)
    )
    for index in corpus:
        back = iv.intervals_to_cells(iv.cells_to_intervals(index))
        assert back == itacart.compact_cells(index), index
        assert itacart.encode_tree(back) == itacart.serialize_to_blob(index)


def test_an_index_with_a_child_under_another_prefix_comes_back_without_it() -> None:
    """The one family the way back does not reproduce, by construction.

    The index writes a border cell and, under the next column's prefix,
    one of its children. Compaction keeps both, since nothing folds; the
    region is the cell alone, and that is what comes back.
    """
    redundant = "NE(2003/0000,2004/0000(3))"
    assert itacart.compact_cells(redundant) == redundant
    assert iv.intervals_to_cells(iv.cells_to_intervals(redundant)) == "NE(2003/0000)"


# --------------------------------------------------------------------------
# Set operations, against sets of integers
# --------------------------------------------------------------------------


def _pairs() -> list[tuple[str, str]]:
    west = itacart.polyfill(_square(10.0, 10.0, 0.012), 5, compact=True)
    east = itacart.polyfill(_square(10.006, 10.004, 0.012), 5, compact=True)
    border = list(itacart.get_children("NE(2003/0000)", flatten=True))
    own = itacart.compose([cell for cell in border if cell.startswith("NE(2003/0000(")])
    cap_part = itacart.compose(list(itacart.get_descendants("NE(0000/1000)", 3))[:4])
    return [
        (west, east),
        ("NE(2003/0000)", own),
        ("NE(2003/0000)", "NE(2004/0000(3))"),
        ("NE(0000/1000)", cap_part),
        (west, "NE(2003/0000)"),
    ]


def test_set_operations_agree_with_sets_of_integers() -> None:
    """Union, intersection, difference, containment, overlap and count."""
    for first, second in _pairs():
        depth = max(
            itacart.get_resolution(cell)
            for index in (first, second)
            for cell in itacart.decompose(index)
        )
        a = iv.cells_to_intervals(first, depth)
        b = iv.cells_to_intervals(second, depth)
        sa, sb = set(_oracle(first, depth)), set(_oracle(second, depth))
        assert _positions(iv.intervals_union(a, b)) == sorted(sa | sb)
        assert _positions(iv.intervals_intersection(a, b)) == sorted(sa & sb)
        assert _positions(iv.intervals_difference(a, b)) == sorted(sa - sb)
        assert _positions(iv.intervals_difference(b, a)) == sorted(sb - sa)
        assert iv.intervals_contain(a, b) == (sb <= sa)
        assert iv.intervals_overlap(a, b) == bool(sa & sb)
        assert iv.intervals_count(a) == len(sa)


def test_the_overlap_of_two_parcels_is_one_intersection() -> None:
    """Two fills that overlap, and the overlap as an index."""
    west, east = _pairs()[0]
    overlap = iv.intervals_intersection(
        iv.cells_to_intervals(west), iv.cells_to_intervals(east)
    )
    assert overlap
    index = iv.intervals_to_cells(overlap)
    assert set(itacart.uncompact_cells(index, 5)) == set(
        itacart.uncompact_cells(west, 5)
    ) & set(itacart.uncompact_cells(east, 5))


def test_lists_of_different_depths_meet_at_the_finer() -> None:
    """A coarser list stands for whole cells, so refining it is exact."""
    coarse = iv.cells_to_intervals("NE(0625/0451(1))", 2)
    fine = iv.cells_to_intervals("NE(0625/0451(2(A1)))", 3)
    union = iv.intervals_union(coarse, fine)
    assert union == iv.intervals_union(
        iv.cells_to_intervals("NE(0625/0451(1))", 3), fine
    )
    assert iv.intervals_to_cells(union) == "NE(0625/0451(1,2(A1)))"
    assert iv.intervals_union(fine, coarse) == union
    assert iv.intervals_difference(fine, coarse) == fine
    assert iv.intervals_union([], fine) == fine
    assert iv.intervals_intersection(coarse, []) == []
    assert iv.intervals_contain(coarse, []) is True


# --------------------------------------------------------------------------
# The existence contract, and the form of a list
# --------------------------------------------------------------------------


def test_a_range_over_an_address_that_names_no_cell_is_refused() -> None:
    """Past the border, the western zero column, a collapsed polar code."""
    for denied, depth in (
        ("NE(2004/0000)", 1),
        ("NW(0000/0451)", 1),
        ("NE(1520/0451)", 1),
        ("NE(0000/1000(2))", 2),
    ):
        value = ii._integer(denied)
        with pytest.raises(NonExistentCellError):
            iv.intervals_to_cells([(value, value + 1)])
        assert ii.int_resolution(value) == depth
    quadrant = (ii.OFFSETS[1], ii.OFFSETS[1] + ii.BOX_SIZES[1] // 4)
    with pytest.raises(NonExistentCellError, match="whole box"):
        iv.intervals_to_cells([quadrant])
    # The fifth child of NE(2003/0000) is NE(2004/0000(3)); the code after it
    # under the same prefix lies past the border.
    border = iv.cells_to_intervals("NE(2003/0000)", 2)
    assert iv.intervals_to_cells(border) == "NE(2003/0000)"
    widened = border[:1] + [(border[1][0], border[1][1] + 1)]
    with pytest.raises(NonExistentCellError, match="2004/0000"):
        iv.intervals_to_cells(widened)


def test_a_quadrant_is_a_range_at_its_own_depth_only() -> None:
    """Its descendants are the resolution-1 lattice, not enumerated."""
    assert iv.cells_to_intervals("NE") == [(0, 1)]
    assert iv.intervals_to_cells([(0, 1)]) == "NE"
    with pytest.raises(ResolutionError):
        iv.cells_to_intervals("NE", 1)


def test_the_depth_must_reach_every_cell() -> None:
    with pytest.raises(ResolutionError):
        iv.cells_to_intervals("NE(0625/0451(3))", 1)
    with pytest.raises(ResolutionError):
        iv.cells_to_intervals("NE(0625/0451)", 14)
    with pytest.raises(NonExistentCellError):
        iv.cells_to_intervals("NW(0000/0451)")


def test_a_list_is_accepted_in_any_order_and_refused_when_it_is_no_list() -> None:
    """Unsorted, overlapping and touching ranges are merged; nonsense is not."""
    a, b = iv.cells_to_intervals("NE(0625/0451(1))", 3), iv.cells_to_intervals(
        "NE(0625/0451(2))", 3
    )
    shuffled = [b[0], a[0], (a[0][0], a[0][0] + 1)]
    assert iv.intervals_count(shuffled) == 50
    assert iv.intervals_union(shuffled, []) == [(a[0][0], b[0][1])]
    for bad in (
        "NE(0625/0451)",
        [(1, 2, 3)],
        [(1.0, 2)],
        [(5, 5)],
        [(-1, 2)],
        [(ii.OFFSETS[1] - 1, ii.OFFSETS[1] + 1)],
        [(0, 1), (ii.OFFSETS[1], ii.OFFSETS[1] + 1)],
        [(0, ii.OFFSETS[-1] + 1)],
    ):
        with pytest.raises(InvalidIndexError):
            iv.intervals_count(bad)  # type: ignore[arg-type]
    with pytest.raises(InvalidIndexError):
        iv.intervals_to_cells([])
