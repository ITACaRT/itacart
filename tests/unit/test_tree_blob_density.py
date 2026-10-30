"""Every figure the TreeBlob specification cites, measured.

The specification's density and information figures were measured once,
written down, and then left to the reader to trust. Three of them were
true only at one parity of the leaf resolution, which the text did not
say, so two sentences that were each true contradicted each other when
read together. This module fixes each figure by a named test, so the
next reading of the page can be checked against the package rather than
against memory.
"""

from __future__ import annotations

import math

import pytest
from shapely.geometry import Polygon

import itacart
from itacart.serialization import tree_blob


def _bits_per_leaf(index: str, leaves: int) -> float:
    return len(tree_blob.encode_tree(index)) * 8 / leaves


def _square(side: float) -> Polygon:
    return Polygon(
        [
            (10.0, 10.0),
            (10.0 + side, 10.0),
            (10.0 + side, 10.0 + side),
            (10.0, 10.0 + side),
        ]
    )


def _leaf_bits(resolution: int) -> int:
    """Bits a leaf at ``resolution`` writes: its component and its count."""
    component = 21 if resolution == 1 else (2 if resolution % 2 == 0 else 5)
    return component + (5 if (resolution + 1) % 2 == 1 and resolution < 13 else 3)


def _node_bits(resolution: int) -> int:
    """Bits an inner node at ``resolution`` writes, by the same rule."""
    return _leaf_bits(resolution)


def _asymptote(leaf_resolution: int) -> float:
    """Bits per leaf of an unbounded complete refinement ending at a parity."""
    total, share, resolution = 0.0, 1.0, leaf_resolution
    while resolution >= 2:
        total += _node_bits(resolution) / share
        share *= 4 if resolution % 2 == 0 else 25
        resolution -= 1
    return total


# --------------------------------------------------------------------------
# The floor and the two asymptotes
# --------------------------------------------------------------------------


def test_the_floor_is_seven_bits_per_leaf() -> None:
    """Two bits of component and five of count, at an even resolution.

    A leaf writes its component and a count field of zero at the width
    its children would have used, so no leaf of any tree costs less.
    """
    assert min(_leaf_bits(resolution) for resolution in range(1, 14)) == 7
    assert [_leaf_bits(r) for r in (1, 2, 3, 12, 13)] == [24, 7, 8, 7, 8]


def test_complete_refinement_costs_8_4_bits_at_odd_and_9_1_at_even() -> None:
    """The densest trees, and the figure depends on the leaves' parity.

    Ten thousand leaves under one cell: at resolution 7 under a resolution-3
    cell, and at resolution 6 under a resolution-2 cell. The asymptotes are
    the series of the field widths, 8.3636 and 9.0909.
    """
    odd = list(itacart.get_descendants("NE(0625/0451(3(C2)))", 7))
    even = list(itacart.get_descendants("NE(0625/0451(3))", 6))
    assert len(odd) == len(even) == 10000
    assert round(_bits_per_leaf(itacart.compose(odd), 10000), 2) == 8.37
    assert round(_bits_per_leaf(itacart.compose(even), 10000), 2) == 9.10
    assert round(_asymptote(13), 4) == 8.3636
    assert round(_asymptote(12), 4) == 9.0909


@pytest.mark.parametrize(
    "resolution, side, cells, bits",
    [(7, 0.019, 43953, 8.40), (8, 0.0095, 43741, 9.15)],
)
def test_a_contiguous_fill_converges_by_the_parity_of_its_leaves(
    resolution: int, side: float, cells: int, bits: float
) -> None:
    """A real cover approaches the asymptote of its parity from above.

    A square of about two kilometres near 10 degrees north and east,
    filled without compaction: some forty-four thousand cells at an odd
    resolution and at an even one.
    """
    index = itacart.polyfill(_square(side), resolution, compact=False)
    assert itacart.count_cells(index) == cells
    assert round(_bits_per_leaf(index, cells), 2) == bits


def test_twenty_five_complete_siblings_cost_11_2_bits_per_leaf() -> None:
    """The frame is paid once, and twenty-five leaves already dilute it."""
    siblings = "NE(0625/0451(1(" + ",".join(itacart.refinement_alphabet(3)) + ")))"
    assert len(tree_blob.encode_tree(siblings)) == 35
    assert _bits_per_leaf(siblings, 25) == 11.2


# --------------------------------------------------------------------------
# Scattered leaves and single leaves
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "resolution, bits",
    [(5, 54.2), (9, 84.2), (11, 99.1), (12, 106.2), (13, 114.2)],
)
def test_a_scattered_leaf_pays_for_its_own_path(resolution: int, bits: float) -> None:
    """One leaf per resolution-1 cell, three hundred of them in one quadrant.

    Each leaf writes every node from its resolution-1 cell down, so the
    cost is the path: twenty-four bits at resolution 1, seven per even
    level and eight per odd one, plus the frame spread over the three
    hundred. It passes a hundred only from resolution 12 on.
    """
    leaves = []
    for column in range(100, 400):
        cell = f"NE({column:04d}/0451)"
        while itacart.get_resolution(cell) < resolution:
            cell = next(iter(itacart.get_children(cell, flatten=True)))
        leaves.append(cell)
    path = sum(_node_bits(level) for level in range(1, resolution + 1))
    expected = 8 * math.ceil((300 * path + 24 + 18) / 8) / 300
    measured = _bits_per_leaf(itacart.compose(leaves), 300)
    assert measured == expected
    assert round(measured, 1) == bits


def test_one_leaf_is_twenty_bytes_at_resolution_13_and_nine_at_1() -> None:
    """The frame dominates a single leaf, and the node is the smaller object."""
    deep = "NE(0625/0451(3(C2(1(A1(2(B3(4(D5(1(E2(3(A4)))))))))))))"
    assert len(tree_blob.encode_tree(deep)) == 20
    assert len(tree_blob.encode_node(deep)) == 10
    assert len(tree_blob.encode_tree("NE(0625/0451)")) == 9
    assert len(tree_blob.encode_node("NE(0625/0451)")) == 4


# --------------------------------------------------------------------------
# The node against its information content
# --------------------------------------------------------------------------


def test_a_resolution_13_address_carries_62_8_bits_of_box_and_62_1_of_cells() -> None:
    """Two figures, and the specification once gave only the first.

    The box is every position the fields can name: four quadrants, 2 005
    columns by 1 001 rows, and the alphabets of twelve levels. The cells
    that exist are fewer, 5 100 221 at resolution 1 and a hundred more for
    every two levels below it.
    """
    box = math.log2(4 * 2005 * 1001 * 4**6 * 25**6)
    cells = math.log2(5100221 * 10**12)
    assert round(box, 2) == 62.80
    assert round(cells, 2) == 62.15
    assert round(80 / box, 2) == 1.27


def test_a_resolution_13_index_is_fifty_five_characters_seven_times_the_box() -> None:
    """The text the specification measures the format against.

    Every resolution-13 spelling has the same length: two letters of
    quadrant, the opening bracket and nine characters of pair, then a
    bracket and a code per level -- one character at the six even levels,
    two at the six odd ones -- and thirteen closing brackets. As eight-bit
    text that is seven times the information the address carries.
    """
    deep = "NE(0625/0451(3(C2(1(A1(2(B3(4(D5(1(E2(3(A4)))))))))))))"
    assert len(deep) == 2 + 1 + 9 + 6 * 2 + 6 * 3 + 13 == 55
    assert len(str(itacart.geo_to_cell(-46.6328862, -23.5508962, 13))) == 55
    box = math.log2(4 * 2005 * 1001 * 4**6 * 25**6)
    assert round(55 * 8 / box) == 7


def test_a_resolution_13_node_is_73_bits_of_fields_and_7_of_padding() -> None:
    """The node writes no count field.

    A marker byte with the resolution, two bits of quadrant, twenty-one of
    resolution-1 pair, two for each of six even levels and five for each
    of six odd ones, then zeros to the byte.
    """
    deep = "NE(0625/0451(3(C2(1(A1(2(B3(4(D5(1(E2(3(A4)))))))))))))"
    node = tree_blob.encode_node(deep)
    bits = "".join(f"{byte:08b}" for byte in node)
    fields = 8 + 2 + 21 + 6 * 2 + 6 * 5
    assert fields == 73
    assert len(bits) == 80
    assert node[0] == 0xA0 | 13
    assert bits[fields:] == "0000000"
    assert 21 - math.log2(2005 * 1001) + 30 - 6 * math.log2(25) == pytest.approx(
        2.20, abs=0.01
    )
