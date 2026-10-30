"""The TreeBlob that 1.0.0 published does not move.

A TreeBlob is addressed by its bytes: a hash over them is a key, and a key
written into storage under 1.0.0 must still name the same region under
every later version. So the contract is the bytes. This module fixes them
by hash over an enumerated corpus, for the three ways a TreeBlob is made:
:func:`itacart.encode_tree`, :func:`itacart.serialize_to_blob` with its
default compaction, and :func:`itacart.geometry_to_tree` from a
GeometryBlob. The GeometryBlob's own bytes are fixed alongside, because
it reads its vertices through the same parser.

The digests were computed with the code of release 1.0.0 and are compared
here with the code under test. One family of input is excluded, and the
exclusion is asserted rather than assumed: an index that writes a cell and
one of its descendants, whose bytes 1.0.0 got wrong -- it kept the
descendant and lost the rest of the cell -- and which now encodes the cell.
"""

from __future__ import annotations

import hashlib
from typing import Iterable

import itacart
from tests.unit.test_totality import CORPUS as TOTALITY_CORPUS

QUADRANTS = ("NE", "NW", "SE", "SW")

#: Index literals the serialization suites pin, less the malformed ones,
#: the polar caps and the nested ones, which have their own tests.
SERIALIZATION_VECTORS = (
    "NE(0001/0001)",
    "NE(0001/0001,0002/0001),NW(0003/0004)",
    "NE(0002/0001)",
    "NE(0500/0300(1))",
    "NE(0500/0300(1,2(B2)))",
    "NE(0500/0300)",
    "NE(0625/0000)",
    "NE(0625/0451(1))",
    "NE(0625/0451(1),0625/0451(2))",
    "NE(0625/0451(1,2))",
    "NE(0625/0451(1,2),0626/0451)",
    "NE(0625/0451(1,2,3,4))",
    "NE(0625/0451(2,1))",
    "NE(0625/0451(2,1),0626/0451)",
    "NE(0625/0451(2,1,4,3))",
    "NE(0625/0451(3(A1)))",
    "NE(0625/0451(3(C2(1(A1(2(B3(4(D5(1(E2(3(A4)))))))))))))",
    "NE(0625/0451(3(C2)))",
    "NE(0625/0451(3))",
    "NE(0625/0451(4,2,1))",
    "NE(0625/0451)",
    "NE(0625/0999)",
    "NE(0626/0451)",
    "NE(2003/0000,2004/0000(3))",
    "NE(2004/0000(3))",
    "NE(2004/0451(1))",
    "NW(0003/0004),NE(0001/0001)",
    "NW(0625/0451(1(E1(3(B2(4(A2)))))))",
    "NW(0625/0451(1(E1(3(B2(4(A2)),B3(4(C5)))))))",
    "NW(0625/0451(1(E1(3(B2)))))",
    "NW(0625/0451)",
    "SE(0005/0006)",
    "SW(0001/0999(4(E5)))",
    "SW(0007/0008)",
    "SE(1400/0374(3(C2(3(C2(4(C1)))))))",
    "SE(1400/0374(3(C2(3))))",
    "NW(0625/0451(1(E1(3(B2(4(A2,B2,B3,B4,C2,C3,C4,C5,D1,D2,D3,D4,D5,E1,E2,E3,"
    "E4,E5)),B3(3(C1,D1,D2,E1,E2,E3,E4)),C2(1(A5,B5,C5,D4,D5,E4,E5),2,3(A4,A5,"
    "B3,B4,B5,C3,C4,C5,D2,D3,D4,D5,E2,E3,E4,E5),4),C3(1,2(A1,B1,B2,B3,C1,C2,C3,"
    "C4,D1,D2,D3,D4,E1,E2,E3,E4,E5),3,4),D2(1(A1,A2,A3,A4,A5,B2,B3,B4,B5,C3,C4,"
    "C5,D5),2,4(A3,A4,A5,B5)),D3(1,2,3(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,C2,C3,C4,"
    "C5,D4,D5),4),D4(1(D1,E1),3(A1,B1,C1,D1,E1)))))))",
)


def _last_columns() -> list[str]:
    """The last column of every row below the caps, in four quadrants."""
    side = itacart.cell_size(1)
    cells = []
    for quadrant in QUADRANTS:
        for row in range(1000):
            column = itacart.last_lattice_column(quadrant, row, side)
            cells.append(f"{quadrant}({column:04d}/{row:04d})")
    return cells


def _families() -> list[str]:
    """Meridian column, first columns, children of a border cell, all levels."""
    entries = []
    for quadrant in QUADRANTS:
        for row in (0, 1, 451, 999):
            first = 0 if quadrant[1] == "E" else 1
            entries.append(f"{quadrant}({first:04d}/{row:04d})")
            entries.append(f"{quadrant}({first + 1:04d}/{row:04d})")
    for parent in ("NE(2003/0000)", "SW(1784/0300)", "NE(0000/0300)", "SE(0001/0999)"):
        children = list(itacart.get_children(parent, flatten=True))
        entries.append(itacart.compose(children))
        entries.extend(children)
    codes: list[str] = []
    for code in ("3", "C2", "1", "A1", "2", "B3", "4", "D5", "1", "E2", "3", "A4"):
        codes.append(code)
        path = "".join(f"({token}" for token in codes)
        entries.append(f"NE(0625/0451{path}" + ")" * (len(codes) + 1))
    return entries


def corpus() -> list[str]:
    """Every entry, in a fixed order."""
    entries = [
        entry for family in sorted(TOTALITY_CORPUS) for entry in TOTALITY_CORPUS[family]
    ]
    return entries + list(SERIALIZATION_VECTORS) + _families() + _last_columns()


def _nested(index: str) -> bool:
    """Whether ``index`` writes a cell and also something spelled under it."""
    cells = itacart.decompose(index)
    return any(
        other.startswith(cell.rstrip(")") + "(") for cell in cells for other in cells
    )


def _exists(index: str) -> bool:
    valid = itacart.is_valid_cell(index)
    return valid if isinstance(valid, bool) else all(valid)


def _rings() -> list[tuple[list[list[str]], str, int]]:
    """GeometryBlob inputs: vertex rings at two resolutions, and points."""
    from shapely.geometry import LineString, Polygon

    shapes = []
    for resolution in (7, 13):
        square = Polygon(
            [(10.0, 10.0), (10.001, 10.0), (10.001, 10.001), (10.0, 10.001)]
        )
        shapes.append(
            ([itacart.vertex_to_cell(square, resolution)], "POLYGON", resolution)
        )
        line = LineString([(-0.001, -0.001), (0.0015, 0.002), (0.003, -0.0005)])
        shapes.append(
            ([itacart.vertex_to_cell(line, resolution)], "LINESTRING", resolution)
        )
    for cell in ("NE(0625/0451)", "SW(0007/0008)", "NE(2003/0000)", "SE(0000/0999)"):
        shapes.append(([[cell]], "POINT", 1))
    return shapes


def _digest(blobs: Iterable[bytes]) -> str:
    hasher = hashlib.sha256()
    for blob in blobs:
        hasher.update(len(blob).to_bytes(4, "big"))
        hasher.update(blob)
    return hasher.hexdigest()


def _digests() -> dict[str, str]:
    entries = corpus()
    existing = [entry for entry in entries if _exists(entry)]
    geometries = [
        itacart.encode_geometry(rings, kind, resolution=resolution)
        for rings, kind, resolution in _rings()
    ]
    return {
        "encode_tree": _digest(itacart.encode_tree(entry) for entry in entries),
        "serialize_to_blob": _digest(
            itacart.serialize_to_blob(entry) for entry in existing
        ),
        "encode_geometry": _digest(geometries),
        "geometry_to_tree": _digest(
            itacart.geometry_to_tree(blob) for blob in geometries
        ),
    }


#: Computed with release 1.0.0, at ``a56c5d2``.
RELEASE_1_0_0 = {
    "encode_tree": (
        "67a8d85f2920c6ca1e638c1f925e323c" "fd0d4811c2cd7fc63890371fca85c5b1"
    ),
    "serialize_to_blob": (
        "afff2bc19f1cbf119d6a895329a15ee7" "ce7b79f32a76d7c2068a9686b1a7b486"
    ),
    "encode_geometry": (
        "db0ce6e822e5b92011bc655d7d7434a1" "a18d5fde0f128517cec0321ee5af29cb"
    ),
    "geometry_to_tree": (
        "1e3fec4a167151fb6a9cd761398995e0" "7794f2c27c05c10fc3f7be4a65aa924b"
    ),
}


def test_the_corpus_is_what_it_says_and_holds_no_nested_entry() -> None:
    """Sizes declared, and the one excluded family is absent by assertion."""
    entries = corpus()
    assert len(entries) == 35 + len(SERIALIZATION_VECTORS) + len(_families()) + 4000
    assert len(entries) == 4138
    assert sum(_exists(entry) for entry in entries) == 4133
    assert not any(_nested(entry) for entry in entries)
    assert _nested("NE(0500/0300,0500/0300(1))")
    assert not _nested("NE(2003/0000,2004/0000(3))")


def test_the_bytes_of_release_1_0_0_do_not_move() -> None:
    """The four digests are the ones 1.0.0 produced on the same corpus."""
    assert _digests() == RELEASE_1_0_0


def test_the_digest_sees_a_single_changed_byte() -> None:
    """The control: one flipped bit anywhere changes the digest."""
    blobs = [itacart.encode_tree(entry) for entry in SERIALIZATION_VECTORS]
    flipped = list(blobs)
    flipped[-1] = flipped[-1][:-1] + bytes([flipped[-1][-1] ^ 1])
    assert _digest(flipped) != _digest(blobs)
