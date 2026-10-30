"""The integer index, measured as a bijection and as an order.

Four claims, each by enumeration and each in four quadrants, because the
enumerator of this grid has been validated in one quadrant and been wrong
in the other three twice. The integer and the cell go back and forth
without loss over every cell of whole resolution-1 rows, the border family,
the caps and their descendants, and chains to resolution 13. The integers
of one resolution sort in the canonical order the TreeBlob writes. The
parent is lexical, by division, and agrees with the string's. And the
existence contract holds in both directions: what names no cell gets no
integer, and what is no cell's integer gets no spelling.
"""

from __future__ import annotations

import pytest

import itacart
from itacart import integer_index as ii
from itacart.constants import CELL_SIZE_M
from itacart.exceptions import (
    InvalidIndexError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)
from itacart.hierarchy import _sort_key
from itacart.serialization.tree_blob import recompose_to_prefix_form

QUADRANTS = ("NE", "NW", "SE", "SW")
SIDE_1 = CELL_SIZE_M[1]
assert SIDE_1 is not None

#: Whole resolution-1 rows: the equator, a middle row, the last row of
#: lattice cells and the row of the caps.
WHOLE_ROWS = (0, 451, 999, 1000)


def _row_positions(quadrant: str, row: int) -> list[str]:
    """Every column the box holds in one row, cell or not."""
    return [f"{quadrant}({column:04d}/{row:04d})" for column in range(ii.BOX_COLUMNS)]


def _canonical(cells: list[str]) -> list[str]:
    """The order the TreeBlob writes, which does not compact.

    The codec does not serialize the polar caps in this version, so the
    cells of row 1000 are ordered by the key compaction uses instead, which
    is the second implementation of the same order in the package.
    """
    if any(cell[8:12] == "1000" for cell in cells):
        return sorted(cells, key=lambda c: _sort_key(itacart.split_components(c)))
    return itacart.decompose(recompose_to_prefix_form(itacart.compose(cells)))


def _chains() -> list[str]:
    """Chains through physical children from four quadrants to resolution 13."""
    starts = []
    for quadrant in QUADRANTS:
        for row in (0, 300, 900, 999):
            top = itacart.last_lattice_column(quadrant, row, SIDE_1)
            for column in (0, 1, 1000, top):
                cell = f"{quadrant}({column:04d}/{row:04d})"
                if itacart.is_valid_cell(cell):
                    starts.append(cell)
    starts += ["NE(0000/1000)", "SE(0000/1000)"]
    cells = []
    for position, start in enumerate(starts):
        cell = start
        cells.append(cell)
        for depth in range(2, 14):
            children = list(itacart.get_children(cell, flatten=True))
            cell = children[(7 * position + depth) % len(children)]
            cells.append(cell)
    return cells


CHAINS = _chains()


# --------------------------------------------------------------------------
# The box and its type
# --------------------------------------------------------------------------


def test_the_index_fits_a_signed_64_bit_integer() -> None:
    """By arithmetic, with the box sizes it rests on declared.

    Four quadrants at resolution 0; four quadrants of 2 005 columns by
    1 001 rows at resolution 1; and a hundred times more for every two
    levels below, so resolution 13 holds 8 028 020 times ten to the
    twelfth. The sum of the fourteen boxes is the limit, and the largest
    integer is one below it.
    """
    assert ii.BOX_SIZES[0] == 4
    assert ii.BOX_SIZES[1] == 4 * 2005 * 1001 == 8028020
    assert ii.BOX_SIZES[13] == 8028020 * 10**12
    assert ii.OFFSETS[0] == 0 and ii.OFFSETS[1] == 4
    assert ii.INTEGER_LIMIT == sum(ii.BOX_SIZES) == 8433475555555150104
    assert ii.INTEGER_LIMIT - 1 < 2**63 - 1
    assert round(1 - ii.INTEGER_LIMIT / 2**63, 4) == 0.0856


def test_resolution_zero_is_the_four_quadrants() -> None:
    """The quadrants are cells to the predicate, so they have integers."""
    assert [ii.int_to_cell(value) for value in range(4)] == list(QUADRANTS)
    assert [ii.cell_to_int(quadrant) for quadrant in QUADRANTS] == [0, 1, 2, 3]
    assert ii.int_parent(ii.cell_to_int("SW(0001/0451)")) == 3
    with pytest.raises(MinResolutionError):
        ii.int_parent(0)


# --------------------------------------------------------------------------
# The bijection, by enumeration
# --------------------------------------------------------------------------


def test_every_position_of_whole_rows_is_a_cell_and_back_or_refused() -> None:
    """Whole rows in four quadrants: cells round-trip, the rest is refused.

    The column zero west of the meridian, the columns a row lacks past the
    border, the prefix column 2004 and every column of row 1000 but the
    caps' have integers in the box and name no cell. Each is refused by
    :func:`~itacart.integer_index.int_to_cell`, and each cell comes back.
    """
    cells = 0
    refused = 0
    for quadrant in QUADRANTS:
        for row in WHOLE_ROWS:
            values = []
            for position in _row_positions(quadrant, row):
                value = ii._integer(position)
                values.append(value)
                if itacart.is_valid_cell(position):
                    assert ii.cell_to_int(position) == value
                    assert ii.int_to_cell(value) == position
                    cells += 1
                else:
                    with pytest.raises(NonExistentCellError):
                        ii.int_to_cell(value)
                    refused += 1
            assert values == sorted(values)
            assert len(set(values)) == ii.BOX_COLUMNS
    assert cells + refused == 4 * len(WHOLE_ROWS) * 2005
    # East of the meridian a row holds column zero and the west does not;
    # row 999 holds two cells in the east and one in the west, and row
    # 1000 only the cap, in the east.
    assert cells == 2 * (2004 + 1520 + 2 + 1) + 2 * (2003 + 1519 + 1 + 0) == 14100
    assert refused == 4 * len(WHOLE_ROWS) * 2005 - cells


def test_the_border_family_and_its_children_round_trip() -> None:
    """The last column of every row in four quadrants, and its children.

    The children include the ones spelled under the next column's prefix,
    so the parent check here meets the case where the lexical parent is
    an address and not a cell.
    """
    lexical_parent_is_no_cell = 0
    for quadrant in QUADRANTS:
        for row in range(0, 1000, 7):
            column = itacart.last_lattice_column(quadrant, row, SIDE_1)
            border = f"{quadrant}({column:04d}/{row:04d})"
            value = ii.cell_to_int(border)
            assert ii.int_to_cell(value) == border
            for child in itacart.get_children(border, flatten=True):
                child_value = ii.cell_to_int(child)
                assert ii.int_to_cell(child_value) == child
                prefix = str(itacart.get_parent(child))
                assert ii.int_parent(child_value) == ii._integer(prefix)
                if not itacart.is_valid_cell(prefix):
                    lexical_parent_is_no_cell += 1
                    with pytest.raises(NonExistentCellError):
                        ii.int_to_cell(ii.int_parent(child_value))
    assert lexical_parent_is_no_cell > 0


def test_the_caps_and_their_descendants_round_trip() -> None:
    """The caps exist in two quadrants, and their row collapses codes.

    Every descendant to resolution 5 comes back, and the codes the polar
    row collapses away -- the off-meridian codes of the sub-row that holds
    the pole -- have integers that name no cell.
    """
    descendants = [
        cell
        for cap in ("NE(0000/1000)", "SE(0000/1000)")
        for resolution in (2, 3, 4, 5)
        for cell in itacart.get_descendants(cap, resolution)
    ]
    assert len(descendants) == 2498
    for cell in descendants + ["NE(0000/1000)", "SE(0000/1000)"]:
        assert ii.int_to_cell(ii.cell_to_int(cell)) == cell
    for collapsed in ("NE(0000/1000(2))", "NE(0000/1000(3))", "SE(0000/1000(4))"):
        assert itacart.is_valid_cell(collapsed) is False
        with pytest.raises(NonExistentCellError):
            ii.int_to_cell(ii._integer(collapsed))
    for missing in ("NW(0000/1000)", "SW(0000/1000)", "NE(0001/1000)"):
        with pytest.raises(NonExistentCellError):
            ii.int_to_cell(ii._integer(missing))


def test_chains_to_resolution_13_round_trip_with_their_parents() -> None:
    """Every level of every chain: the integer, its resolution, its parent."""
    assert len(CHAINS) == 13 * 50
    for cell in CHAINS:
        value = ii.cell_to_int(cell)
        resolution = itacart.get_resolution(cell)
        assert ii.int_to_cell(value) == cell
        assert ii.int_resolution(value) == resolution
        assert ii.OFFSETS[resolution] <= value < ii.OFFSETS[resolution + 1]
        assert ii.int_parent(value) == ii._integer(str(itacart.get_parent(cell)))
        assert ii.int_parent(value, target_res=0) == ii.cell_to_int(cell[:2])


# --------------------------------------------------------------------------
# The order
# --------------------------------------------------------------------------


def test_the_integers_of_one_resolution_sort_in_canonical_order() -> None:
    """Sorting by integer is sorting the way the TreeBlob writes.

    At each resolution of the chains, the cells sorted by integer are the
    cells in the order the canonical prefix form lists them, which orders
    quadrants, then column before row, then codes by their alphabet.
    """
    for resolution in range(1, 14):
        cells = sorted({c for c in CHAINS if itacart.get_resolution(c) == resolution})
        assert len(cells) > 30
        by_integer = sorted(cells, key=lambda c: ii.cell_to_int(c))
        assert by_integer == _canonical(cells)
        lattice = [cell for cell in cells if cell[8:12] != "1000"]
        assert [c for c in by_integer if c in lattice] == _canonical(lattice)


def test_an_interior_family_is_a_run_of_consecutive_integers() -> None:
    """Twenty-five siblings and four, in alphabet order, with no gaps."""
    for parent in ("NE(0625/0451(1))", "NE(0625/0451)", "SW(0007/0008(4(C3)))"):
        children = list(itacart.get_children(parent, flatten=True))
        values = [ii.cell_to_int(child) for child in children]
        assert values == list(range(values[0], values[0] + len(values)))


def test_a_border_family_is_not_one_run() -> None:
    """The physical children of a border cell lie in more than one range.

    The integer is lexical, so a child spelled under another prefix sits
    with that prefix's family. ``SW(0001/0999(4))`` has thirty children,
    twenty-one of them under its own spelling.
    """
    parent = "SW(0001/0999(4))"
    children = list(itacart.get_children(parent, flatten=True))
    values = sorted(ii.cell_to_int(child) for child in children)
    runs = 1 + sum(1 for a, b in zip(values, values[1:]) if b != a + 1)
    assert len(children) == 30
    assert len([c for c in children if c.startswith("SW(0001/0999(4(")]) == 21
    assert runs > 1


# --------------------------------------------------------------------------
# The existence contract, both ways
# --------------------------------------------------------------------------


def test_a_spelling_that_names_no_cell_gets_no_integer() -> None:
    """Refused with the package's exception, the existence one."""
    for denied in ("NW(0000/0451)", "NE(2004/0000)", "NE(1520/0451)", "NE(0001/1000)"):
        with pytest.raises(NonExistentCellError):
            ii.cell_to_int(denied)


def test_something_that_is_not_an_integer_of_the_index_is_refused() -> None:
    """Never ``ValueError``: an integer index is an index."""
    for bad in (True, False, 4.0, "4", None, -1, ii.INTEGER_LIMIT):
        with pytest.raises(InvalidIndexError):
            ii.int_to_cell(bad)  # type: ignore[arg-type]
    with pytest.raises(InvalidIndexError):
        ii.int_resolution([4, 5.0])  # type: ignore[list-item]
    with pytest.raises(ResolutionError):
        ii.int_parent(ii.cell_to_int("NE(0625/0451(3))"), target_res=2)
    with pytest.raises(ResolutionError):
        ii.int_parent(ii.cell_to_int("NE(0625/0451(3))"), target_res=-1)


def test_an_integer_declared_by_numpy_is_an_integer() -> None:
    """A value read out of an integer array is taken; a float is not."""
    numpy = pytest.importorskip("numpy")
    assert ii.int_to_cell(numpy.int64(4)) == "NE(0000/0000)"
    assert ii.int_to_cell(numpy.array([4, 5])) == ["NE(0000/0000)", "NE(0000/0001)"]
    with pytest.raises(InvalidIndexError):
        ii.int_to_cell(numpy.float64(4.0))


# --------------------------------------------------------------------------
# Vectorised semantics and names
# --------------------------------------------------------------------------


def test_a_compositional_index_answers_in_decompose_order() -> None:
    """One integer per cell, aligned with :func:`itacart.decompose`."""
    index = "NE(0625/0451(2,1),0626/0451),SW(0007/0008)"
    values = ii.cell_to_int(index)
    assert isinstance(values, list)
    assert values == [ii.cell_to_int(cell) for cell in itacart.decompose(index)]
    assert ii.int_to_cell(values) == itacart.decompose(index)
    assert ii.int_resolution(values) == [2, 2, 1, 1]
    parents = ii.int_parent(values)
    assert isinstance(parents, list) and parents[0] == parents[1]


def test_the_h3_names_are_the_same_functions() -> None:
    assert ii.str_to_int is ii.cell_to_int
    assert ii.int_to_str is ii.int_to_cell
    assert itacart.str_to_int("NE") == 0
