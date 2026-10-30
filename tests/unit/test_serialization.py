"""Tests for :mod:`itacart.serialization`.

Ported from the reference implementation's own suite, which covers both
formats: 63 test functions over the tree codec and 71 over the geometry
codec. The note that used to stand here claimed nothing was portable and
cited a search for ``neighbor``, ``grid_disk`` and ``adjacen`` — a
topology measurement pasted into a serialization file. Searching for the
names this module actually defines finds the two files.

Ported tests are rewritten rather than copied: the reference implements a
partial grid, without the meridian triangle and without the trapezoidal
cell's surplus child, so its index space runs from one where this one
runs from zero and stops one column short of what the grammar admits.
"""

from __future__ import annotations

import pytest

import itacart
from itacart.exceptions import (
    InvalidIndexError,
    MalformedBlobError,
    ResolutionError,
)
from itacart.serialization import tree_blob

PAPER_CELL = "NE(0625/0451)"
DEEP_CELL = "NE(0625/0451(3(C2(1(A1(2(B3(4(D5(1(E2(3(A4)))))))))))))"


def _children_of(cell: str) -> list[str]:
    """Flatten the sibling groups :func:`itacart.get_children` yields."""
    flat: list[str] = []
    for group in itacart.get_children(cell):
        flat.extend(group if isinstance(group, list) else [group])
    return flat


def _last_column(quadrant: str, row: int) -> int:
    return itacart.last_lattice_column(quadrant, row, itacart.cell_size(1))


# --------------------------------------------------------------------------
# The eight formal properties
# --------------------------------------------------------------------------


def test_encoding_is_deterministic() -> None:
    """The same index encodes to the same bytes every time."""
    index = "NE(0625/0451(2,1),0626/0451)"
    assert tree_blob.encode_tree(index) == tree_blob.encode_tree(index)


def test_decode_of_encode_is_the_recomposed_form() -> None:
    """Decoding recovers the canonical form, not the input spelling."""
    index = "NE(0625/0451(2,1))"
    decoded = tree_blob.decode_tree(tree_blob.encode_tree(index))
    assert decoded == tree_blob.recompose_to_prefix_form(index)
    assert decoded == "NE(0625/0451(1,2))"


def test_encoding_is_invariant_under_recomposition() -> None:
    """Recomposing before encoding changes nothing."""
    index = "NE(0625/0451(1),0625/0451(2))"
    recomposed = tree_blob.recompose_to_prefix_form(index)
    assert tree_blob.encode_tree(index) == tree_blob.encode_tree(recomposed)


def test_recomposition_is_idempotent() -> None:
    """Canonical form is a fixed point."""
    once = tree_blob.recompose_to_prefix_form("NE(0625/0451(4,2,1))")
    assert tree_blob.recompose_to_prefix_form(once) == once


def test_round_trip_is_bit_exact() -> None:
    """Encoding a decoded blob reproduces the blob."""
    blob = tree_blob.encode_tree("NE(0625/0451(2,1),0626/0451)")
    assert tree_blob.encode_tree(tree_blob.decode_tree(blob)) == blob


def test_the_blob_is_addressed_by_its_leaf_set() -> None:
    """Walk-order and prefix-index spellings of one leaf set agree."""
    walk = "NE(0625/0451(1),0625/0451(2))"
    prefix = "NE(0625/0451(2,1))"
    assert tree_blob.encode_tree(walk) == tree_blob.encode_tree(prefix)


def test_a_resolution_13_leaf_occupies_ten_bytes() -> None:
    """The size property the format is named for."""
    node = tree_blob.encode_node(DEEP_CELL)
    assert tree_blob.resolution_of_binary(node) == 13
    assert len(node) == 10


def test_binary_prefix_truncation_is_monotone() -> None:
    """Truncating twice equals truncating once to the coarser level."""
    node = tree_blob.encode_node("NE(0625/0451(3(C2)))")
    twice = tree_blob.prefix_at_resolution_binary(
        tree_blob.prefix_at_resolution_binary(node, 2), 1
    )
    assert twice == tree_blob.prefix_at_resolution_binary(node, 1)
    assert twice == tree_blob.encode_node(PAPER_CELL)


# --------------------------------------------------------------------------
# The index space, which the reference implementation models partially
# --------------------------------------------------------------------------


def test_the_meridian_triangle_encodes() -> None:
    """Column zero is a cell, and it is the triangle on the meridian.

    The reference implementation has no triangle: its grid is
    parallelograms only, so it drops the meridian column and numbers the
    first parallelogram as column one. Porting that offset would make
    every column and row in this package encode as its neighbour.
    """
    for quadrant in ("NE", "SE"):
        cell = f"{quadrant}(0000/0000)"
        assert itacart.is_triangular_cell(cell)
        assert itacart.is_valid_cell(cell)
        assert tree_blob.decode_tree(tree_blob.encode_tree(cell)) == cell


def test_the_equatorial_row_encodes() -> None:
    """Row zero is a row. The reference numbers rows from one."""
    cell = "NE(0625/0000)"
    assert tree_blob.decode_tree(tree_blob.encode_tree(cell)) == cell


def test_the_meridian_triangle_is_absent_west_of_the_meridian() -> None:
    """The triangle is one cell addressed from the east, not two.

    The reference implementation refuses column zero everywhere, which
    is right in the west for the wrong reason and wrong in the east.
    """
    for quadrant in ("NW", "SW"):
        for probe in ("0000/0451", "0000/0451(1)", "0000/0451(1(A1))"):
            assert not itacart.is_valid_cell(f"{quadrant}({probe})")
        with pytest.raises(InvalidIndexError, match=f"not addressable in {quadrant}"):
            tree_blob.encode_tree(f"{quadrant}(0000/0451)")


def test_the_codec_does_not_re_derive_the_boundary() -> None:
    """Format validity and cell existence are different questions.

    Column 2004 is well formed wherever the grammar admits it, but it
    only *exists* where the last column is 2003, and the last column
    retreats with the cosine of the latitude. Deriving that per cell
    inside the codec would put a second copy of the boundary model in
    the serializer. The existence oracle stays where it already is.
    """
    well_formed = "NE(2004/0451(1))"
    tree_blob.encode_tree(well_formed)
    assert not itacart.is_valid_cell(well_formed)
    assert itacart.last_lattice_column("NE", 451, itacart.cell_size(1)) == 1519


def test_the_polar_caps_are_refused_as_a_restriction_of_the_codecs() -> None:
    """Row 1000 is in the grid, and this version does not serialize it.

    The two caps exist, so a message calling their row outside the grid
    would state a property of the grid that is false. Both codecs refuse
    them, and both say why in the same words, because the GeometryBlob
    reads its vertices through the tree parser. The cells used are the
    real caps: a column of row 1000 that names no cell would be refused
    either way and could not tell the two messages apart.
    """
    from itacart.exceptions import GeometryError

    message = (
        "row 1000 holds the polar caps, which this version of the binary "
        "codecs does not serialize"
    )
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        assert itacart.is_valid_cell(cap) is True
        for encoder in (
            tree_blob.encode_tree,
            tree_blob.encode_node,
            tree_blob.recompose_to_prefix_form,
            itacart.serialize_to_blob,
        ):
            with pytest.raises(InvalidIndexError, match=message):
                encoder(cap)
        with pytest.raises(GeometryError, match=message):
            itacart.encode_geometry([[cap]], "POINT", resolution=1)


def test_a_row_past_the_caps_is_outside_the_grid() -> None:
    """From row 1001 on the old message stays, because it is true there."""
    for row in (1001, 9999):
        with pytest.raises(InvalidIndexError, match=f"row {row} is outside the grid"):
            tree_blob.encode_tree(f"NE(0000/{row:04d})")


def _with_field(blob: bytes, offset: int, width: int, value: int) -> bytes:
    """``blob`` with the ``width``-bit field at bit ``offset`` set to ``value``."""
    bits = "".join(f"{byte:08b}" for byte in blob)
    bits = bits[:offset] + f"{value:0{width}b}" + bits[offset + width :]
    return bytes(int(bits[i : i + 8], 2) for i in range(0, len(bits), 8))


def test_a_decoded_cap_row_is_refused_with_the_same_reason() -> None:
    """The decoders share the reason with the parser, row for row.

    A blob is never written with row 1000, so these are blobs made to lie.
    The row field sits after the three header bytes, the two quadrant
    bits, the sixteen-bit group count and the eleven-bit column in a tree
    blob, and after the resolution byte, the quadrant and the column in a
    node blob. Row 1023, the largest the field holds, is the control: it
    keeps the message that says it is outside the grid.
    """
    tree = tree_blob.encode_tree("NE(0625/0000)")
    node = tree_blob.encode_node("NE(0625/0000)")
    tree_row, node_row = 24 + 2 + 16 + 11, 8 + 2 + 11
    assert tree_blob.validate_tree(_with_field(tree, tree_row, 10, 999)) is None
    assert tree_blob.decode_node(_with_field(node, node_row, 10, 999)) == (
        "NE(0625/0999)"
    )
    for row, reason in ((1000, "row 1000 holds the polar caps"), (1023, "outside")):
        with pytest.raises(MalformedBlobError, match=reason):
            tree_blob.validate_tree(_with_field(tree, tree_row, 10, row))
        with pytest.raises(MalformedBlobError, match=reason):
            tree_blob.decode_node(_with_field(node, node_row, 10, row))


# --------------------------------------------------------------------------
# The trapezoid exception, enumerated rather than sampled
# --------------------------------------------------------------------------


def test_the_trapezoids_surplus_child_encodes() -> None:
    """A trapezoidal parent yields one child past the refinement ratio.

    The surplus child lands in column 2004, one past the last column
    that is a cell. Both facts break a codec that reads the refinement
    ratio as a child count and the last column as an upper bound.
    """
    parent = f"NE({_last_column('NE', 0):04d}/0000)"
    assert itacart.is_trapezoidal_cell(parent)
    children = _children_of(parent)
    assert len(children) == itacart.refinement_ratio(2) + 1
    assert any(child.startswith("NE(2004/") for child in children)
    index = itacart.compose(children)
    blob = tree_blob.encode_tree(index)
    assert tree_blob.count_vertices(blob) == len(children)
    assert tree_blob.decode_tree(blob) == tree_blob.recompose_to_prefix_form(index)


@pytest.mark.slow
def test_every_last_column_of_every_row_encodes() -> None:
    """Enumerated over all 1000 rows, not sampled.

    The surplus child is one cell out of thousands and lives only in the
    last column of a row, so sampling rows at random does not reach the
    family reliably. Each row's last column is computed and every child
    of it is encoded and read back.
    """
    checked = 0
    for row in range(1000):
        parent = f"NE({_last_column('NE', row):04d}/{row:04d})"
        for child in _children_of(parent):
            assert tree_blob.decode_node(tree_blob.encode_node(child)) == child
            checked += 1
    assert checked > 1000


def test_column_2004_is_a_prefix_and_never_a_cell() -> None:
    """It may carry refinements; it may not stand alone."""
    assert itacart.is_valid_cell("NE(2004/0000(3))")
    assert not itacart.is_valid_cell("NE(2004/0000)")
    tree_blob.encode_tree("NE(2004/0000(3))")
    with pytest.raises(InvalidIndexError, match="column 2004 is not addressable"):
        tree_blob.encode_tree("NE(2004/0000)")


def test_the_column_past_the_prefix_column_is_refused() -> None:
    with pytest.raises(InvalidIndexError, match="outside the grid"):
        tree_blob.encode_tree("NE(2005/0000(3))")


def test_the_child_count_field_is_bounded_by_its_width() -> None:
    """Five children fit where the refinement ratio says four.

    The count field is three bits wide at an even child resolution, so
    the surplus child costs nothing. Reading the ratio as the bound is
    what refuses the whole eastern border.
    """
    assert tree_blob._max_children_at(2) == 7
    assert tree_blob._max_children_at(3) == 31
    assert tree_blob._max_children_at(2) > itacart.refinement_ratio(2)


# --------------------------------------------------------------------------
# Canonical forms: recompose is not normalize
# --------------------------------------------------------------------------


def test_recompose_orders_siblings_without_compacting_them() -> None:
    """A complete sibling set survives recomposition and not normalisation.

    ``normalize`` collapses four children into their parent, which
    preserves coverage and changes the leaf set. A blob's identity is
    its leaf set, so the codec cannot use that form.
    """
    complete = "NE(0625/0451(2,1,4,3))"
    assert tree_blob.recompose_to_prefix_form(complete) == "NE(0625/0451(1,2,3,4))"
    assert itacart.normalize(complete) == "NE(0625/0451)"


def test_a_cell_written_whole_absorbs_its_descendants_on_every_port() -> None:
    """The codec keeps the region a nested spelling names.

    ``NE(0500/0300,0500/0300(1))`` is a valid index of two cells, one
    inside the other. Merging the two occurrences of the parent used to
    hand the whole one the other's children, so every port stored the
    daughter and lost the parent: a quarter of the region. The index layer
    absorbs the daughter instead, and the codec now does the same, so the
    four ports agree with the plain parent byte for byte, written in
    either order and nested more than one level deep.
    """
    parent = "NE(0500/0300)"
    nested = (
        "NE(0500/0300,0500/0300(1))",
        "NE(0500/0300(1),0500/0300)",
        "NE(0500/0300(1(A1)),0500/0300)",
        "NE(0500/0300,0500/0300(1(A1),2))",
    )
    for index in nested:
        assert itacart.count_cells(index) > 1
        assert itacart.compose(itacart.decompose(index)) == parent
        assert tree_blob.encode_tree(index) == tree_blob.encode_tree(parent)
        assert tree_blob.recompose_to_prefix_form(index) == parent
        assert tree_blob.serialize_to_blob(index, compact=False) == (
            tree_blob.encode_tree(parent)
        )
        assert tree_blob.encode_node(index) == tree_blob.encode_node(parent)

    # One level down, the whole occurrence is the daughter, and it absorbs
    # the granddaughter while its sibling stays.
    deeper = "NE(0500/0300(1,1(A1),2(B2)))"
    assert tree_blob.recompose_to_prefix_form(deeper) == "NE(0500/0300(1,2(B2)))"
    assert itacart.compose(itacart.decompose(deeper)) == "NE(0500/0300(1,2(B2)))"


def test_absorption_is_lexical_like_the_index_layer() -> None:
    """A daughter spelled under another column's prefix is not absorbed.

    ``NE(2004/0000(3))`` lies inside ``NE(2003/0000)`` and is spelled under
    the prefix of the next column, so neither the index layer nor the
    codec sees it as written under the parent. Both keep the two leaves,
    which loses nothing: the parent already covers the daughter.
    """
    border = "NE(2003/0000,2004/0000(3))"
    assert itacart.normalize(border) == border
    assert itacart.compact_cells(border) == border
    assert tree_blob.recompose_to_prefix_form(border) == border
    assert tree_blob.decode_tree(tree_blob.encode_tree(border)) == border


def test_the_vertex_parser_keeps_its_merge() -> None:
    """The GeometryBlob is not this repair's to change.

    It parses vertices through the tree parser's merge, and that merge
    keeps its old rule when it is not asked to absorb. A two-cell nested
    spelling handed in as one vertex still reads as the daughter, as it
    did in 1.0.0; asking for absorption reads it as the parent.
    """
    from itacart.serialization.tree_blob import _merge, _parse, _render

    nested = "NE(0500/0300,0500/0300(1))"
    blob = itacart.encode_geometry([[nested]], "POINT", resolution=2)
    assert itacart.decode_geometry(blob)["rings"] == [["NE(0500/0300(1))"]]
    assert _render(_merge(_parse(nested))) == "NE(0500/0300(1))"
    assert _render(_merge(_parse(nested), absorb=True)) == "NE(0500/0300)"


def test_compacting_changes_the_blob_on_purpose() -> None:
    """``serialize_to_blob`` compacts, and so encodes a different set."""
    index = "NE(0625/0451(1,2,3,4))"
    assert tree_blob.serialize_to_blob(index, compact=True) != tree_blob.encode_tree(
        index
    )
    assert tree_blob.serialize_to_blob(index, compact=False) == tree_blob.encode_tree(
        index
    )


def test_a_blob_stores_a_region_and_not_a_filling() -> None:
    """The fill resolution is not in the bytes because it never was.

    Two fills run at different resolutions can produce the same string,
    so no decoder can recover which one ran. This is a declared property
    of the format, not a defect to be repaired later: ``max(resolution)``
    is a lower bound on the fill resolution and never the value itself.
    """
    index = "NE(0625/0451)"
    blob = tree_blob.encode_tree(index)
    assert tree_blob.deserialize_from_blob(blob) == index
    at_two = list(itacart.uncompact_cells(index, 2))
    at_three = list(itacart.uncompact_cells(index, 3))
    assert len(at_two) != len(at_three)
    expanded = itacart.compose(at_two)
    assert tree_blob.decode_tree(tree_blob.encode_tree(expanded)) != index


# --------------------------------------------------------------------------
# Node blobs and binary operations
# --------------------------------------------------------------------------


def test_node_round_trip() -> None:
    for cell in (PAPER_CELL, "NE(0625/0451(3))", "SW(0001/0999(4(E5)))"):
        assert tree_blob.decode_node(tree_blob.encode_node(cell)) == cell


def test_encode_node_refuses_a_branching_index() -> None:
    with pytest.raises(InvalidIndexError, match="exactly one cell"):
        tree_blob.encode_node("NE(0625/0451(1,2))")


def test_ancestry_holds_on_encoded_nodes() -> None:
    parent = tree_blob.encode_node(PAPER_CELL)
    child = tree_blob.encode_node("NE(0625/0451(3))")
    assert tree_blob.is_ancestor_binary(parent, child)
    assert not tree_blob.is_ancestor_binary(child, parent)
    assert not tree_blob.is_ancestor_binary(
        tree_blob.encode_node("NW(0625/0451)"), child
    )


def test_ancestry_accepts_a_tree_as_the_descendant() -> None:
    parent = tree_blob.encode_node(PAPER_CELL)
    tree = tree_blob.encode_tree("NE(0625/0451(1,2))")
    assert tree_blob.is_ancestor_binary(parent, tree)
    assert not tree_blob.is_ancestor_binary(
        tree_blob.encode_node("NE(0626/0451)"), tree
    )


def test_an_ancestor_of_a_tree_must_dominate_every_leaf() -> None:
    """Covering one branch is not covering the tree.

    Found by running the reference suite against this port: an ``any``
    quantifier here would call a resolution-5 node the ancestor of a
    tree that merely touches it.
    """
    parent = tree_blob.encode_node("NW(0625/0451(1(E1(3(B2)))))")
    both = tree_blob.encode_tree("NW(0625/0451(1(E1(3(B2(4(A2)),B3(4(C5)))))))")
    only_b2 = tree_blob.encode_tree("NW(0625/0451(1(E1(3(B2(4(A2)))))))")
    assert not tree_blob.is_ancestor_binary(parent, both)
    assert tree_blob.is_ancestor_binary(parent, only_b2)


def test_an_empty_descendant_blob_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="non-empty bytes"):
        tree_blob.is_ancestor_binary(tree_blob.encode_node(PAPER_CELL), b"")


def test_prefix_at_a_finer_resolution_is_refused() -> None:
    node = tree_blob.encode_node(PAPER_CELL)
    with pytest.raises(ResolutionError, match="not a prefix"):
        tree_blob.prefix_at_resolution_binary(node, 2)


def test_leaves_stream_without_text() -> None:
    blob = tree_blob.encode_tree("NE(0625/0451(1,2),0626/0451)")
    leaves = list(tree_blob.iter_leaves(blob))
    assert len(leaves) == 3 == tree_blob.count_vertices(blob)
    assert tree_blob.decode_node(leaves[0]) == "NE(0625/0451(1))"


# --------------------------------------------------------------------------
# Malformed input
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "index",
    [
        "",
        "XX(0625/0451)",
        "NE(0625/0451",
        "NE(0625/0451))",
        "NE(0625/0451(9))",
        "NE(0625/0451(3(ZZ)))",
        "NE()",
    ],
)
def test_malformed_indices_are_refused(index: str) -> None:
    with pytest.raises(InvalidIndexError):
        tree_blob.encode_tree(index)


def test_a_non_string_index_is_refused() -> None:
    with pytest.raises(InvalidIndexError, match="must be a string"):
        tree_blob.encode_tree(b"NE(0625/0451)")  # type: ignore[arg-type]


def test_descending_past_the_finest_resolution_is_refused() -> None:
    too_deep = DEEP_CELL.replace("A4)", "A4(1))")
    with pytest.raises(InvalidIndexError, match="cannot descend"):
        tree_blob.encode_tree(too_deep)


@pytest.mark.parametrize(
    ("blob", "message"),
    [
        (b"", "shorter than"),
        (b"\x00\x10\x00", "bad magic"),
        (b"\xc7\x00\x00", "unsupported format version"),
        (b"\xc7\x11\x00", "reserved flag"),
        (b"\xc7\x10\x01", "reserved header"),
        (b"\xc7\x10\x00\x00\x00\x00", "carries no cell"),
        (b"\xc7\x10\x00\x00\x40", "ends inside a field"),
    ],
)
def test_malformed_blobs_are_refused(blob: bytes, message: str) -> None:
    with pytest.raises(MalformedBlobError, match=message):
        tree_blob.validate_tree(blob)


def test_a_non_bytes_blob_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="must be bytes"):
        tree_blob.validate_tree("not bytes")  # type: ignore[arg-type]


def test_trailing_padding_must_be_zero() -> None:
    blob = bytearray(tree_blob.encode_tree("NE(0625/0451(3))"))
    blob[-1] |= 0b0000_0001
    with pytest.raises(MalformedBlobError, match="padding is not zero"):
        tree_blob.validate_tree(bytes(blob))


def test_unread_bytes_after_the_body_are_refused() -> None:
    blob = tree_blob.encode_tree(PAPER_CELL) + b"\x00\x00"
    with pytest.raises(MalformedBlobError, match="unread bits"):
        tree_blob.validate_tree(blob)


def test_an_empty_node_blob_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="empty node blob"):
        tree_blob.resolution_of_binary(b"")


def test_a_node_blob_with_the_wrong_magic_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="bad node magic"):
        tree_blob.resolution_of_binary(b"\x00\x00")


def test_a_node_resolution_outside_the_grid_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="resolution 15 out of range"):
        tree_blob.resolution_of_binary(b"\xaf\x00")


def test_a_resolution_13_leaf_inside_a_tree_still_emits_its_count() -> None:
    """The count field of the deepest leaf has no children to describe."""
    blob = tree_blob.encode_tree(DEEP_CELL)
    assert tree_blob.count_vertices(blob) == 1
    assert tree_blob.decode_tree(blob) == DEEP_CELL


def test_a_corrupt_row_in_a_tree_blob_is_refused() -> None:
    """Ten bits hold 1023, the grid holds 1000, so a blob can lie."""
    blob = bytearray(tree_blob.encode_tree("NE(0625/0451)"))
    blob[6] |= 0b0001_1111  # row bits 0-4, leaving the column alone
    blob[7] |= 0b1111_1000  # row bits 5-9, leaving the count at zero
    with pytest.raises(MalformedBlobError, match="row .* outside the grid"):
        tree_blob.validate_tree(bytes(blob))


def test_a_corrupt_column_in_a_node_blob_is_refused() -> None:
    """Eleven bits hold 2047, the grammar admits 2004."""
    node = bytearray(tree_blob.encode_node("NE(0625/0451)"))
    node[1] |= 0b0011_1111
    node[2] |= 0b1111_1111
    with pytest.raises(MalformedBlobError, match="not addressable"):
        tree_blob.decode_node(bytes(node))


def test_a_corrupt_row_in_a_node_blob_is_refused() -> None:
    node = bytearray(tree_blob.encode_node("NE(0625/0451)"))
    node[2] |= 0b0000_0111  # row bits 0-2
    node[3] |= 0b1111_1110  # row bits 3-9, leaving the padding bit at zero
    with pytest.raises(MalformedBlobError, match="row .* outside the grid"):
        tree_blob.decode_node(bytes(node))


def test_a_reserved_refinement_code_in_a_node_blob_is_refused() -> None:
    node = bytearray(tree_blob.encode_node("NE(0625/0451(3(A1)))"))
    node[-1] |= 0b1111_1111
    with pytest.raises(MalformedBlobError, match="reserved refinement code"):
        tree_blob.decode_node(bytes(node))


def test_a_reserved_refinement_code_is_refused() -> None:
    """Codes 25 to 31 of a five-bit field address nothing."""
    blob = bytearray(tree_blob.encode_tree("NE(0625/0451(3(A1)))"))
    blob[9] |= 0b0111_1100  # the five-bit resolution-3 field, at bit 73
    with pytest.raises(MalformedBlobError, match="reserved"):
        tree_blob.validate_tree(bytes(blob))


def test_a_blob_claiming_children_below_the_finest_resolution_is_refused() -> None:
    """The count field of the deepest leaf may only say zero."""
    blob = bytearray(tree_blob.encode_tree(DEEP_CELL))
    blob[18] |= 0b0000_0001  # first bit of the final three-bit count field
    blob[19] |= 0b1100_0000  # its remaining two bits
    with pytest.raises(MalformedBlobError, match="admits no children"):
        tree_blob.validate_tree(bytes(blob))


def test_a_corrupt_column_in_a_tree_blob_is_refused() -> None:
    """Eleven bits hold 2047; the grammar admits 2004, and only refined."""
    blob = bytearray(tree_blob.encode_tree("NE(0625/0451)"))
    blob[5] |= 0b1111_1111  # column bits 0-7
    blob[6] |= 0b1110_0000  # column bits 8-10
    with pytest.raises(MalformedBlobError, match="not addressable"):
        tree_blob.validate_tree(bytes(blob))


# --------------------------------------------------------------------------
# More than one quadrant in a single blob
# --------------------------------------------------------------------------


def test_an_index_may_name_more_than_one_quadrant() -> None:
    """`compose` writes this form, so the codec has to read it.

    A blob carries one group per quadrant. Refusing the multi-quadrant
    spelling would make the canonical output of a core function
    unencodable, which is a defect in the codec and not in the index.
    """
    index = "NE(0001/0001),NE(0002/0001),NW(0003/0004)"
    blob = tree_blob.encode_tree(index)
    assert tree_blob.count_vertices(blob) == 3
    assert tree_blob.decode_tree(blob) == tree_blob.recompose_to_prefix_form(index)
    assert tree_blob.decode_tree(blob) == "NE(0001/0001,0002/0001),NW(0003/0004)"


def test_a_region_straddling_the_equator_encodes() -> None:
    """Crossing a quadrant boundary is ordinary, not exotic.

    Any region over the equator lands in two quadrants, and `polyfill`
    produces the two-group index directly.
    """
    from shapely.geometry import Polygon

    region = Polygon([(10.0, -0.02), (10.04, -0.02), (10.04, 0.02), (10.0, 0.02)])
    index = itacart.polyfill(region, 4, compact=False)
    assert sorted({cell[:2] for cell in itacart.decompose(index)}) == ["NE", "SE"]
    blob = tree_blob.encode_tree(index)
    assert tree_blob.count_vertices(blob) == itacart.count_cells(index)
    assert tree_blob.decode_tree(blob) == tree_blob.recompose_to_prefix_form(index)


def test_compose_output_is_always_encodable() -> None:
    """All four quadrants at once, which is the widest an index gets."""
    cells = [
        "NE(0001/0001)",
        "NE(0002/0001)",
        "NW(0003/0004)",
        "SE(0005/0006)",
        "SW(0007/0008)",
    ]
    index = itacart.compose(cells)
    blob = tree_blob.encode_tree(index)
    assert tree_blob.decode_tree(blob) == index
    assert tree_blob.count_vertices(blob) == len(cells)


def test_quadrant_groups_are_ordered_so_the_blob_stays_addressable() -> None:
    """Two spellings of the same cross-quadrant set give one blob."""
    forward = "NE(0001/0001),NW(0003/0004)"
    backward = "NW(0003/0004),NE(0001/0001)"
    assert tree_blob.encode_tree(forward) == tree_blob.encode_tree(backward)


def test_a_blob_repeating_a_quadrant_is_refused() -> None:
    """Groups are merged before encoding, so a repeat can only be damage."""
    good = bytearray(tree_blob.encode_tree("NE(0001/0001),NW(0003/0004)"))
    good[8] &= 0b1100_1111  # the second group's quadrant, at bit 66, back to NE
    with pytest.raises(MalformedBlobError, match="appears twice"):
        tree_blob.validate_tree(bytes(good))


def test_a_node_blob_still_encodes_exactly_one_cell() -> None:
    with pytest.raises(InvalidIndexError, match="exactly one cell"):
        tree_blob.encode_node("NE(0001/0001),NW(0003/0004)")
