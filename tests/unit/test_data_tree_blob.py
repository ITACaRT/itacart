"""The DataTreeBlob: values on a lexical tree of cells, read in place.

The two cases the format was designed against are built here in full and
enumerated, not sampled: a dense field of ten thousand leaves under one
resolution-3 cell, averaged up the tree; and a drone obstacle map of ten
thousand one-metre leaves under one resolution-5 cell, with two buildings,
a tree canopy over one of them and a strip nobody surveyed. Every
aggregate is checked against the leaves by an independent count, and
every number the format promises -- an exact mean, a block's share, a
width, a refusal -- is pinned here.
"""

from __future__ import annotations

import functools
import math
import random
import struct
from fractions import Fraction
from typing import Any

import pytest

import itacart
from itacart import (
    DataVariable,
    data_tree_schema,
    decode_data_tree,
    encode_data_tree,
    read_data_cell,
    read_data_path,
    slice_data_tree,
    update_data_tree,
)
from itacart.exceptions import (
    IncompatibleProfileError,
    InvalidIndexError,
    MalformedBlobError,
    NonAtomicIndexError,
    NonExistentCellError,
    ResolutionError,
    SerializationError,
)
from itacart.integer_index import OFFSETS, RATIOS
from itacart.serialization import data_tree_blob as dtb

ITA = (-45.8722, -23.2124)
DENIED_WEST = "NW(0000/0110)"


def _f32(value: float) -> float:
    return float(struct.unpack("<f", struct.pack("<f", value))[0])


def _lattice(cell: str, root_level: int) -> tuple[int, int]:
    """Column and row of a leaf inside its root, read from its codes.

    A two-by-two family at an even level, five-by-five at an odd one, so
    four levels below a resolution-3 or resolution-5 root make a hundred by
    a hundred. Deterministic and cheap, which is all a synthetic raster
    needs.
    """
    level = itacart.get_resolution(cell)
    position = itacart.cell_to_int(cell) - OFFSETS[level]
    digits = []
    for finer in range(level, root_level, -1):
        position, digit = divmod(position, RATIOS[finer])
        digits.append((finer, digit))
    u = v = 0
    for finer, digit in reversed(digits):
        side = 2 if RATIOS[finer] == 4 else 5
        u, v = u * side + digit % side, v * side + digit // side
    return u, v


# --------------------------------------------------------------------------
# Case A: a dense field, averaged up the tree
# --------------------------------------------------------------------------

DENSE = (
    DataVariable("n", "u32", op="sum", absent="zero", leaf_max=1),
    DataVariable("z", "f32", op="mean", weight="n"),
    DataVariable("cls", "u8", op="mode", weight="n"),
    DataVariable("qual", "u8", species="level"),
)


@functools.lru_cache(maxsize=None)
def dense_case() -> tuple[str, dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    root = itacart.geo_to_cell(ITA[0], ITA[1], 3)
    data = {}
    for leaf in itacart.get_descendants(root, 7):
        u, v = _lattice(leaf, 3)
        z = _f32(
            612.0
            + 38.0
            * math.sin(2.1 * math.pi * u / 100)
            * math.cos(1.7 * math.pi * v / 100)
            + 55.0 * math.exp(-((u - 35) ** 2 + (v - 45) ** 2) / 200.0)
        )
        data[leaf] = {"n": 1, "z": z, "cls": 0 if z < 590 else 1 if z < 620 else 2}
    level5 = sorted({itacart.get_parent(leaf, 5) for leaf in data})
    observed: dict[str, dict[str, Any]] = {
        cell: {"qual": (k * 7) % 4} for k, cell in enumerate(level5)
    }
    # A parent and a daughter with values of their own, at other levels.
    observed[itacart.get_parent(level5[0], 4)] = {"qual": 9}
    observed[next(iter(itacart.get_children(level5[0], flatten=True)))] = {"qual": 3}
    return root, data, observed


@functools.lru_cache(maxsize=None)
def dense_blob() -> bytes:
    _, data, observed = dense_case()
    return encode_data_tree(data, DENSE, observed=observed)


# --------------------------------------------------------------------------
# Case B: a drone obstacle map
# --------------------------------------------------------------------------

DRONE = (
    DataVariable("obs", "u32", op="sum", absent="zero", leaf_max=1),
    DataVariable("occ", "u32", op="sum", leaf_max=1),
    DataVariable("bld", "u32", op="sum", leaf_max=1),
    DataVariable("veg", "u32", op="sum", leaf_max=1),
    DataVariable("hmax", "u16", op="max"),
    DataVariable("fence", "u8", species="level"),
)


def _obstacles(u: int, v: int) -> list[tuple[str, int]]:
    """Two buildings and a canopy over the second, heights in decimetres."""
    hits = []
    if 12 <= u < 46 and 10 <= v < 42:
        hits.append(("building", 240))
    if (55 <= u < 90 and 50 <= v < 62) or (55 <= u < 68 and 62 <= v < 88):
        hits.append(("building", 95))
    if (u - 70) ** 2 + (v - 70) ** 2 < 121:
        hits.append(("vegetation", 140))
    return hits


def _surveyed(u: int, v: int) -> bool:
    return not (u >= 80 and v < 32)


@functools.lru_cache(maxsize=None)
def drone_case() -> tuple[str, dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    root = itacart.geo_to_cell(ITA[0], ITA[1], 5)
    data: dict[str, dict[str, Any]] = {}
    for leaf in itacart.get_descendants(root, 9):
        u, v = _lattice(leaf, 5)
        if not _surveyed(u, v):
            data[leaf] = {"obs": 0, "occ": None, "bld": None, "veg": None, "hmax": None}
            continue
        hits = _obstacles(u, v)
        data[leaf] = {
            "obs": 1,
            "occ": int(bool(hits)),
            "bld": int(any(kind == "building" for kind, _ in hits)),
            "veg": int(any(kind == "vegetation" for kind, _ in hits)),
            "hmax": max((height for _, height in hits), default=0),
        }
    level7 = sorted({itacart.get_parent(leaf, 7) for leaf in data})
    observed = {
        cell: {"fence": 2 if k in (0, 1, 10, 11) else 1 if k % 9 == 4 else 0}
        for k, cell in enumerate(level7)
    }
    return root, data, observed


@functools.lru_cache(maxsize=None)
def drone_blob(prune: bool = False) -> bytes:
    _, data, observed = drone_case()
    return encode_data_tree(data, DRONE, observed=observed, prune=prune)


@functools.lru_cache(maxsize=None)
def _groups(case: str) -> dict[int, list[dict[str, Any]]]:
    """Every node of a case's tree, by integer, with the leaf values beneath it."""
    _, data, _ = dense_case() if case == "dense" else drone_case()
    groups: dict[int, list[dict[str, Any]]] = {}
    for leaf, values in data.items():
        integer = itacart.cell_to_int(leaf)
        level = itacart.int_resolution(integer)
        groups.setdefault(integer, []).append(values)
        for ancestor in range(level - 1, -1, -1):
            groups.setdefault(itacart.int_parent(integer, ancestor), []).append(values)
    return groups


def _under(case: str, cell: str) -> list[dict[str, Any]]:
    return _groups(case)[itacart.cell_to_int(cell)]


# --------------------------------------------------------------------------
# The dense field
# --------------------------------------------------------------------------


def test_the_dense_field_comes_back_exactly_and_from_any_order() -> None:
    """Ida e volta exata, and the bytes do not depend on the order of the pairs."""
    root, data, observed = dense_case()
    assert len(data) == 10_000
    blob = dense_blob()
    decoded = decode_data_tree(blob)
    for leaf, values in data.items():
        assert decoded[leaf] == values
    assert decoded[root]["n"] == 10_000

    shuffled = list(data.items())
    random.Random(4).shuffle(shuffled)
    assert encode_data_tree(shuffled, DENSE, observed=list(observed.items())) == blob

    schema = data_tree_schema(blob)
    assert schema["resolution"] == 7 and schema["deepest"] == 7
    assert schema["nodes"] == (1, 1, 1, 1, 4, 100, 400, 10_000)
    assert schema["variables"] == DENSE and schema["pruned"] is False


def test_the_dense_field_aggregates_exactly() -> None:
    """Counts are counts, the mean is the exact weighted mean rounded once.

    The mean is rounded at every level, so the root's drifts from the exact
    mean of all ten thousand leaves; the drift is pinned at one f32 ulp.
    """
    decoded = decode_data_tree(dense_blob())
    nodes = [cell for cell in decoded if itacart.get_resolution(cell) in (3, 4, 5, 6)]
    assert len(nodes) == 505
    for cell in nodes:
        leaves = _under("dense", cell)
        assert decoded[cell]["n"] == len(leaves)
        level = itacart.get_resolution(cell)
        daughters = [
            decoded[d] for d in itacart.get_children(cell, flatten=True) if d in decoded
        ]
        exact = sum(Fraction(d["z"]) * d["n"] for d in daughters) / sum(
            d["n"] for d in daughters
        )
        assert decoded[cell]["z"] == dtb._to_f32(exact.numerator, exact.denominator)
        tally: dict[int, int] = {}
        for d in daughters:
            tally[d["cls"]] = tally.get(d["cls"], 0) + d["n"]
        best = max(tally.values())
        assert decoded[cell]["cls"] == min(k for k, c in tally.items() if c == best)
        if level == 3:
            truth = sum(Fraction(leaf["z"]) for leaf in leaves) / len(leaves)
            ulp = math.ulp(_f32(float(truth))) * 2**29  # an f32 ulp at this size
            assert abs(decoded[cell]["z"] - float(truth)) <= ulp


def test_a_parent_and_its_daughter_keep_their_own_values() -> None:
    """A level variable is observed, not derived, at every level it lives on."""
    _, _, observed = dense_case()
    blob = dense_blob()
    decoded = decode_data_tree(blob)
    for cell, values in observed.items():
        assert decoded[cell]["qual"] == values["qual"]
        assert read_data_cell(blob, cell)["qual"] == values["qual"]
    parent = next(c for c in observed if itacart.get_resolution(c) == 4)
    daughter = next(c for c in observed if itacart.get_resolution(c) == 6)
    assert (decoded[parent]["qual"], decoded[daughter]["qual"]) == (9, 3)
    # Where it was not observed, a level variable is unknown; where it has no
    # values at all, it is not there.
    other = next(
        c for c in decoded if itacart.get_resolution(c) == 4 and c not in observed
    )
    assert decoded[other]["qual"] is None
    assert "qual" not in decoded[next(iter(dense_case()[1]))]


def test_a_read_lands_on_one_ancestor_per_level() -> None:
    """The path of a cell at resolution r has r + 1 nodes, the quadrant first."""
    root, data, _ = dense_case()
    blob = dense_blob()
    leaf = sorted(data)[1234]
    path = read_data_path(blob, leaf)
    assert [itacart.get_resolution(cell) for cell, _ in path] == list(range(8))
    assert path[0][0] == "SW" and path[3][0] == root and path[-1][0] == leaf
    decoded = decode_data_tree(blob)
    assert all(values == decoded[cell] for cell, values in path)
    assert read_data_path(blob, itacart.cell_to_int(leaf)) == path

    # A deep chain, down to resolution 13, read at every depth.
    chain = itacart.geo_to_cell(ITA[0], ITA[1], 13)
    deep = encode_data_tree(
        {chain: {"n": 1}}, [DataVariable("n", "u32", op="sum", leaf_max=1)]
    )
    for level in range(14):
        ancestor = chain if level == 13 else itacart.get_parent(chain, level)
        found = read_data_path(deep, ancestor)
        assert len(found) == level + 1 and found[-1] == (ancestor, {"n": 1})


def test_reads_take_spellings_integers_and_many_cells() -> None:
    root, data, _ = dense_case()
    blob = dense_blob()
    leaves = sorted(data)[:3]
    composed = itacart.compose(leaves)
    expected = [data[leaf] for leaf in itacart.decompose(composed)]
    assert read_data_cell(blob, composed) == expected
    assert read_data_cell(blob, [itacart.cell_to_int(c) for c in leaves]) == [
        data[c] for c in leaves
    ]
    assert read_data_cell(blob, leaves) == [data[c] for c in leaves]
    paths = read_data_path(blob, composed)
    assert len(paths) == 3 and all(
        path[-1][1] == data[c] for path, c in zip(paths, itacart.decompose(composed))
    )
    for bad, error in (
        (DENIED_WEST, NonExistentCellError),
        ("garbage", InvalidIndexError),
        (1.5, InvalidIndexError),
        (True, InvalidIndexError),
        (b"NE", InvalidIndexError),
        (-1, InvalidIndexError),
        (OFFSETS[1] + 2004 * 1001, NonExistentCellError),  # NE(2004/0000)
    ):
        with pytest.raises(error):
            read_data_cell(blob, bad)
    with pytest.raises(NonExistentCellError):
        read_data_path(blob, DENIED_WEST)


def test_a_cell_the_blob_does_not_hold_answers_with_absence() -> None:
    blob = dense_blob()
    elsewhere = "NE(0625/0451(1))"
    assert read_data_cell(blob, elsewhere) == {"n": 0, "z": None, "cls": None}
    assert read_data_path(blob, elsewhere) == []
    sibling = itacart.get_parent(dense_case()[0], 2)
    near = next(
        c for c in itacart.get_children(sibling, flatten=True) if c != dense_case()[0]
    )
    assert [cell for cell, _ in read_data_path(blob, near)] == [
        "SW",
        itacart.get_parent(sibling),
        sibling,
    ]
    # Finer than the observation resolution is refused, not guessed.
    with pytest.raises(
        ResolutionError, match="finer than the observation resolution 7"
    ):
        read_data_cell(blob, itacart.geo_to_cell(ITA[0], ITA[1], 8))


# --------------------------------------------------------------------------
# The drone map
# --------------------------------------------------------------------------


def test_the_drone_map_counts_union_classes_and_what_was_observed() -> None:
    """Binary leaves, integer counts above them, and the percentage derived."""
    root, data, observed = drone_case()
    blob = drone_blob()
    decoded = decode_data_tree(blob)
    for leaf, values in data.items():
        assert decoded[leaf] == values
    for cell, values in decoded.items():
        if itacart.get_resolution(cell) == 9:
            continue
        leaves = _under("drone", cell)
        known = [leaf for leaf in leaves if leaf["occ"] is not None]
        assert values["obs"] == len(known)
        expected_occ = sum(leaf["occ"] for leaf in known) if known else None
        assert values["occ"] == expected_occ
        heights = [leaf["hmax"] for leaf in known]
        assert values["hmax"] == (max(heights) if heights else None)
    top = decoded[root]
    union = sum(1 for v in data.values() if v["occ"])
    both = sum(1 for v in data.values() if v["bld"] and v["veg"])
    assert top["occ"] == union and top["bld"] + top["veg"] - top["occ"] == both > 0
    assert top["obs"] == 10_000 - sum(1 for v in data.values() if v["obs"] == 0)
    # Coverage over what was observed, against coverage over the cell: the
    # unsurveyed strip is why the two differ.
    assert (top["occ"], top["obs"], top["bld"], top["veg"]) == (2069, 9360, 1846, 373)
    assert round(100 * top["occ"] / top["obs"], 2) == 22.1
    assert round(100 * top["occ"] / 10_000, 2) == 20.69
    for cell, values in observed.items():
        assert decoded[cell]["fence"] == values["fence"]
    assert all(
        "fence" not in v for c, v in decoded.items() if itacart.get_resolution(c) != 7
    )


def test_unknown_is_not_free() -> None:
    """A node nobody surveyed is unknown, not empty, and counts no observation."""
    root, data, observed = drone_case()
    blob = drone_blob()
    unseen = [leaf for leaf, values in data.items() if values["obs"] == 0]
    corner = itacart.get_parent(unseen[0], 7)
    assert all(leaf["occ"] is None for leaf in _under("drone", corner))
    assert read_data_cell(blob, corner) == {
        "obs": 0,
        "occ": None,
        "bld": None,
        "veg": None,
        "hmax": None,
        "fence": observed[corner]["fence"],
    }


def test_pruning_folds_only_families_that_agree_on_every_variable() -> None:
    """The bridge's first error, closed: a full family with two heights stays.

    Every leaf reads the same from the pruned blob as from the unpruned one:
    an intensive value is inherited, a count is the block's even share.
    """
    root, data, _ = drone_case()
    full, pruned = drone_blob(), drone_blob(True)
    assert (len(full), len(pruned)) == (73_632, 11_168)
    assert data_tree_schema(pruned)["pruned"] is True
    assert data_tree_schema(pruned)["nodes"] == (1, 1, 1, 1, 1, 1, 4, 100, 96, 1075)
    leaves = sorted(data)
    assert read_data_cell(pruned, leaves) == read_data_cell(full, leaves)
    # A family whose leaves are all occupied but disagree on height.
    mixed = next(
        cell
        for cell in decode_data_tree(full)
        if itacart.get_resolution(cell) == 8
        and all(leaf["occ"] == 1 for leaf in _under("drone", cell))
        and len({leaf["hmax"] for leaf in _under("drone", cell)}) > 1
    )
    beneath = next(iter(itacart.get_children(mixed, flatten=True)))
    assert read_data_path(pruned, beneath)[-1][0] == beneath
    # Nothing is folded at a level where a level variable lives.
    assert data_tree_schema(pruned)["nodes"][7] == 100


def test_a_fractional_union_cannot_be_counted_back() -> None:
    """The bridge's second error, as the arithmetic that makes it one.

    Per-class counts added together exceed the union by exactly the cells
    where classes overlap, so only the binary union leaf, read by OR, gives
    the covered count; and an update by OR changes every ancestor by one.
    """
    root, data, _ = drone_case()
    blob = bytearray(drone_blob())
    top = read_data_cell(blob, root)
    assert top["bld"] + top["veg"] > top["occ"]
    free = next(leaf for leaf, values in sorted(data.items()) if values["occ"] == 0)
    before = [values["occ"] for _, values in read_data_path(blob, free)]
    update_data_tree(blob, free, {"occ": 1, "veg": 1, "hmax": 50})
    after = [values["occ"] for _, values in read_data_path(blob, free)]
    assert [b - a for a, b in zip(before, after)] == [1] * len(before)


def test_the_denominator_of_a_mean_is_the_weight() -> None:
    """The bridge's third error: the mean is over what the weight counts.

    Five of twenty-five daughters observed at 100 and twenty at 0 give 20;
    the same five observed and twenty unknown give 100, over the five.
    """
    n = DataVariable("n", "u32", op="sum", absent="zero", leaf_max=1)
    cover = DataVariable("cover", "f64", op="mean", weight="n")
    parent = "NE(0625/0451(1))"
    daughters = list(itacart.get_children(parent, flatten=True))
    assert len(daughters) == 25
    observed_all = {
        d: {"n": 1, "cover": 100.0 if i < 5 else 0.0} for i, d in enumerate(daughters)
    }
    observed_five = {
        d: ({"n": 1, "cover": 100.0} if i < 5 else {"n": 0, "cover": None})
        for i, d in enumerate(daughters)
    }
    assert read_data_cell(encode_data_tree(observed_all, [n, cover]), parent) == {
        "n": 25,
        "cover": 20.0,
    }
    assert read_data_cell(encode_data_tree(observed_five, [n, cover]), parent) == {
        "n": 5,
        "cover": 100.0,
    }


# --------------------------------------------------------------------------
# Blocks, border and cap
# --------------------------------------------------------------------------

COUNT = DataVariable("n", "u32", op="sum", absent="zero", leaf_max=1)
TALL = DataVariable("h", "f32", op="max")


def test_a_leaf_above_the_resolution_stands_for_its_block() -> None:
    block = "NE(0625/0451(1(C3(2(A1)))))"
    blob = encode_data_tree({block: {"n": 100, "h": 12.5}}, [COUNT, TALL], resolution=7)
    assert itacart.get_resolution(block) == 5
    finest = next(iter(itacart.get_descendants(block, 7)))
    middle = itacart.get_parent(finest)
    assert read_data_cell(blob, finest) == {"n": 1, "h": 12.5}
    assert read_data_cell(blob, middle) == {"n": 25, "h": 12.5}
    sums = encode_data_tree(
        {block: {"s": 1.0}}, [DataVariable("s", "f32", op="sum")], resolution=7
    )
    assert read_data_cell(sums, finest) == {"s": _f32(0.01)}


@pytest.mark.parametrize(
    "data, resolution, message",
    [
        ({"NE(2003/0000)": {"n": 25}}, 2, "absorbs the border"),
        ({"NE(0625/0451(1))": {"n": 24}}, 3, "divided evenly"),
        ({"NE(0625/0451(1))": {"n": 50}}, 3, "outside 0 to 1 per cell over 25"),
    ],
)
def test_a_block_is_refused_where_it_cannot_be_one(
    data: dict[str, dict[str, Any]], resolution: int, message: str
) -> None:
    with pytest.raises(SerializationError, match=message):
        encode_data_tree(data, [COUNT], resolution=resolution)


def test_the_border_is_lexical_and_says_so() -> None:
    """A border cell aggregates its own stem; the surplus hangs from a prefix.

    The prefix ``NE(2004/0000)`` names no cell. It is a node of the blob and
    holds the aggregate of its subtree, and it is never emitted: not by the
    decoder, not on a path.
    """
    parent = "NE(2003/0000)"
    leaves = list(itacart.get_descendants(parent, 3))
    area = DataVariable("area", "f64", op="sum", absent="zero")
    blob = encode_data_tree(
        {
            leaf: {"n": 1, "area": float(itacart.effective_cell_area(leaf))}
            for leaf in leaves
        },
        [COUNT, area],
    )
    decoded = decode_data_tree(blob)
    assert len(leaves) == 125 and "NE(2004/0000)" not in decoded
    assert all(itacart.is_valid_cell(cell) for cell in decoded)
    assert decoded[parent]["n"] == 97 and decoded["NE"]["n"] == 125
    share = decoded[parent]["area"] / itacart.effective_cell_area(parent)
    assert round(share, 4) == 0.7764
    surplus_area = sum(
        itacart.effective_cell_area(leaf)
        for leaf in leaves
        if leaf.startswith("NE(2004/0000")
    )
    assert round(surplus_area / itacart.effective_cell_area(parent), 4) == 0.2239
    surplus = next(leaf for leaf in leaves if leaf.startswith("NE(2004/0000"))
    chain = [itacart.get_parent(surplus, level) for level in range(3)] + [surplus]
    assert chain[1] == "NE(2004/0000)"
    assert [cell for cell, _ in read_data_path(blob, surplus)] == [
        cell for cell in chain if itacart.is_valid_cell(cell)
    ]
    assert data_tree_schema(blob)["nodes"] == (1, 2, 7, 125)


def test_the_cap_enters_and_aggregates_all_of_itself() -> None:
    cap = "NE(0000/1000)"
    level = [cap]
    for _ in range(2, 6):
        level = [c for p in level for c in itacart.get_children(p, flatten=True)]
    blob = encode_data_tree({leaf: {"n": 1} for leaf in level}, [COUNT])
    decoded = decode_data_tree(blob)
    assert len(level) == 1194 and decoded[cap]["n"] == 1194
    assert all(itacart.is_valid_cell(cell) for cell in decoded)
    assert data_tree_schema(blob)["nodes"] == (1, 1, 3, 16, 56, 1194)
    # Twenty of those nodes are prefixes inside the cap that name no cell.
    assert len(decoded) == 1 + 1 + 3 + 16 + 56 + 1194 - 20


# --------------------------------------------------------------------------
# In place
# --------------------------------------------------------------------------


def test_an_update_in_place_is_the_encoding_of_the_new_data() -> None:
    """Byte for byte: the widths depend on the structure, not on the values."""
    root, data, observed = drone_case()
    blob = bytearray(drone_blob())
    leaf = sorted(data)[4321]
    changed = {**data[leaf], "occ": 1, "bld": 1, "hmax": 333}
    update_data_tree(blob, leaf, {"occ": 1, "bld": 1, "hmax": 333})
    assert bytes(blob) == encode_data_tree(
        {**data, leaf: changed}, DRONE, observed=observed
    )
    fence = next(iter(observed))
    update_data_tree(memoryview(blob), fence, {"fence": 7})
    assert read_data_cell(bytes(blob), fence)["fence"] == 7

    _, dense, notes = dense_case()
    field = bytearray(dense_blob())
    target = sorted(dense)[17]
    update_data_tree(field, itacart.cell_to_int(target), {"z": 1000.0, "cls": 2})
    assert bytes(field) == encode_data_tree(
        {**dense, target: {**dense[target], "z": 1000.0, "cls": 2}},
        DENSE,
        observed=notes,
    )
    update_data_tree(field, target, {})
    assert bytes(field) == encode_data_tree(
        {**dense, target: {**dense[target], "z": 1000.0, "cls": 2}},
        DENSE,
        observed=notes,
    )


def test_an_unknown_written_in_place_clears_its_bit() -> None:
    root, data, observed = drone_case()
    blob = bytearray(drone_blob())
    leaf = next(leaf for leaf, values in sorted(data.items()) if values["occ"] == 1)
    update_data_tree(blob, leaf, {"occ": None, "obs": 0, "hmax": None})
    assert read_data_cell(blob, leaf)["occ"] is None
    assert bytes(blob) == encode_data_tree(
        {**data, leaf: {**data[leaf], "occ": None, "obs": 0, "hmax": None}},
        DRONE,
        observed=observed,
    )


@pytest.mark.parametrize(
    "cell, values, error, message",
    [
        ("leaf", {"nope": 1}, SerializationError, "no variable 'nope'"),
        ("parent", {"occ": 1}, SerializationError, "aggregate of its daughters"),
        ("leaf", {"fence": 1}, SerializationError, "no values at resolution 9"),
        ("leaf", {"occ": 2}, SerializationError, "outside 0 to 1 per cell"),
        ("leaf", {"hmax": 70_000}, SerializationError, "outside u16"),
        ("leaf", {"hmax": "tall"}, SerializationError, "not an integer"),
        ("elsewhere", {"occ": 1}, SerializationError, "holds no node"),
        ("leaf", [("occ", 1)], SerializationError, "mapping of variable names"),
        ("many", {"occ": 1}, NonAtomicIndexError, "more than one cell"),
        (DENIED_WEST, {"occ": 1}, NonExistentCellError, "names no cell"),
    ],
)
def test_an_update_refuses_and_writes_nothing(
    cell: str, values: Any, error: type[Exception], message: str
) -> None:
    root, data, _ = drone_case()
    leaf = sorted(data)[10]
    target = {
        "leaf": leaf,
        "parent": itacart.get_parent(leaf),
        "elsewhere": "NE(0625/0451(1))",
        "many": itacart.compose(sorted(data)[:2]),
    }.get(cell, cell)
    blob = bytearray(drone_blob())
    before = bytes(blob)
    with pytest.raises(error, match=message):
        update_data_tree(blob, target, values)
    assert bytes(blob) == before


def test_an_update_needs_a_writable_buffer() -> None:
    with pytest.raises(SerializationError, match="writable buffer"):
        update_data_tree(drone_blob(), sorted(drone_case()[1])[0], {"occ": 1})


def test_a_sum_that_would_overflow_its_type_is_refused_before_writing() -> None:
    small = DataVariable("s", "u8", op="sum")
    cells = list(itacart.get_children("NE(0625/0451)", flatten=True))
    blob = bytearray(encode_data_tree({c: {"s": 60} for c in cells}, [small]))
    before = bytes(blob)
    with pytest.raises(SerializationError, match="overflows u8"):
        update_data_tree(blob, cells[0], {"s": 120})
    assert bytes(blob) == before
    with pytest.raises(SerializationError, match="overflows u8; declare leaf_max"):
        encode_data_tree({c: {"s": 70} for c in cells}, [small])


# --------------------------------------------------------------------------
# Slices
# --------------------------------------------------------------------------


def test_a_slice_is_the_sub_tree_with_its_path() -> None:
    root, data, observed = drone_case()
    full = drone_blob()
    node = sorted(c for c in observed)[37]
    sliced = slice_data_tree(full, node)
    decoded, original = decode_data_tree(sliced), decode_data_tree(full)
    stem = node.rstrip(")")
    beneath = [c for c in original if c == node or c.startswith(stem + "(")]
    assert len(beneath) == 1 + 4 + 100
    for cell in beneath:
        assert decoded[cell] == original[cell]
    path = [c for c in decoded if c not in beneath]
    assert [itacart.get_resolution(c) for c in path] == [0, 1, 2, 3, 4, 5, 6]
    for cell in path:
        assert {k: v for k, v in decoded[cell].items() if k != "fence"} == {
            k: v for k, v in original[node].items() if k != "fence"
        }
    assert data_tree_schema(sliced)["nodes"] == (1, 1, 1, 1, 1, 1, 1, 1, 4, 100)
    assert slice_data_tree(full, node) == sliced
    assert slice_data_tree(full, itacart.cell_to_int(node)) == sliced
    leaf = beneath[-1]
    assert decode_data_tree(slice_data_tree(full, leaf))[leaf] == original[leaf]


def test_a_slice_of_a_node_the_blob_does_not_hold_is_refused() -> None:
    with pytest.raises(SerializationError, match="holds no node"):
        slice_data_tree(drone_blob(), "NE(0625/0451(1))")
    with pytest.raises(NonAtomicIndexError):
        slice_data_tree(drone_blob(), itacart.compose(sorted(drone_case()[1])[:2]))


# --------------------------------------------------------------------------
# Absence, widths and the schema
# --------------------------------------------------------------------------


def test_the_validity_bitmap_keeps_the_whole_domain() -> None:
    """255 is a u8 like any other: unknown is a bit, not a value."""
    code = DataVariable("code", "u8", op="max")
    zero = DataVariable("z", "i16", op="min", absent="zero")
    cells = list(itacart.get_children("NE(0625/0451)", flatten=True))
    blob = encode_data_tree(
        {cells[0]: {"code": 255, "z": -3}, cells[1]: {"code": None}, cells[2]: {}},
        [code, zero],
    )
    assert read_data_cell(blob, cells[0]) == {"code": 255, "z": -3}
    assert read_data_cell(blob, cells[1]) == {"code": None, "z": 0}
    assert read_data_cell(blob, cells[2]) == {"code": None, "z": 0}
    assert read_data_cell(blob, "NE(0625/0451)") == {"code": 255, "z": -3}
    assert read_data_cell(blob, cells[3]) == {"code": None, "z": 0}


def test_a_count_widens_with_the_structure_and_never_overflows() -> None:
    """Declared u8, a count of ten thousand at the root: the width is per level."""
    blob = dense_blob()
    assert read_data_cell(blob, dense_case()[0])["n"] == 10_000
    reader = dtb._Blob(blob)
    assert reader.widths["n"] == ["u16", "u16", "u16", "u16", "u16", "u8", "u8", "u8"]
    assert reader.widths["z"] == ["f32"] * 8
    assert reader.widths["qual"] == [None, None, None, None, "u8", "u8", "u8", None]
    wide = DataVariable("w", "u64", op="sum", leaf_max=2**64 - 1)
    with pytest.raises(SerializationError, match="does not fit 64 bits"):
        encode_data_tree(
            {c: {"w": 1} for c in itacart.get_children("NE(0625/0451)", flatten=True)},
            [wide],
        )


def test_the_operators_aggregate_as_declared() -> None:
    cells = list(itacart.get_children("NE(0625/0451)", flatten=True))
    n = DataVariable("n", "u8", op="sum", absent="zero", leaf_max=1)
    variables = [
        n,
        DataVariable("s", "f64", op="sum"),
        DataVariable("lo", "i32", op="min"),
        DataVariable("mean", "i32", op="mean", weight="n"),
        DataVariable("mode", "u8", op="mode", weight="n"),
        DataVariable("raw", "f64"),
    ]
    blob = encode_data_tree(
        {
            cells[0]: {"n": 1, "s": 0.1, "lo": 5, "mean": 1, "mode": 4, "raw": 1.5},
            cells[1]: {"n": 1, "s": 0.2, "lo": -5, "mean": 2, "mode": 3, "raw": 2.5},
            cells[2]: {"n": 0, "s": 0.3, "lo": 7, "mean": 100, "mode": 9},
            cells[3]: {"n": 1, "s": 1e300, "lo": None, "mean": None, "mode": 3},
        },
        variables,
    )
    parent = read_data_cell(blob, "NE(0625/0451)")
    exact_sum = Fraction(0.1) + Fraction(0.2) + Fraction(0.3) + Fraction(1e300)
    assert parent == {
        "n": 3,
        "s": float(exact_sum),
        "lo": -5,
        "mean": 2,  # (1 + 2) / 2 = 1.5, half to even; the zero weight is skipped
        "mode": 3,
        "raw": None,
    }
    assert read_data_cell(blob, cells[0])["raw"] == 1.5


@pytest.mark.parametrize(
    "variables, message",
    [
        ([], "1 to 255 variables, not 0"),
        ([DataVariable(f"v{i}", "u8") for i in range(256)], "not 256"),
        (["n"], "variables are DataVariable"),
        ([DataVariable("", "u8")], "not 1 to 255 bytes"),
        ([DataVariable("x" * 256, "u8")], "not 1 to 255 bytes"),
        ([DataVariable("a", "u8"), DataVariable("a", "u8")], "declared twice"),
        ([DataVariable("a", "u128")], "type 'u128'"),
        ([DataVariable("a", "u8", species="tree")], "species 'tree'"),
        ([DataVariable("a", "u8", op="median")], "operator 'median'"),
        ([DataVariable("a", "u8", species="level", op="sum")], "never aggregated"),
        ([DataVariable("a", "u8", absent="maybe")], "absence 'maybe'"),
        ([DataVariable("a", "u8", op="max", weight="b")], "only a mean or a mode"),
        ([DataVariable("a", "f32", op="mean")], "needs a weight"),
        ([DataVariable("a", "f32", op="mean", weight="a")], "needs a weight"),
        (
            [
                DataVariable("w", "f32", op="sum"),
                DataVariable("a", "f32", op="mean", weight="w"),
            ],
            "not an integer pyramid sum",
        ),
        ([DataVariable("a", "i8", op="sum", leaf_max=1)], "unsigned pyramid sum"),
        ([DataVariable("a", "u8", op="sum", leaf_max=0)], "not a positive u8"),
        ([DataVariable("a", "u8", op="sum", leaf_max=256)], "not a positive u8"),
        ([DataVariable("a", "u8", op="sum", leaf_max=True)], "not a positive u8"),
    ],
)
def test_the_schema_is_checked(variables: list[Any], message: str) -> None:
    with pytest.raises(SerializationError, match=message):
        encode_data_tree({"NE(0625/0451)": {}}, variables)


def test_one_variable_may_be_given_bare() -> None:
    blob = encode_data_tree({"NE(0625/0451)": {"n": 1}}, COUNT, resolution=1)
    assert data_tree_schema(blob)["variables"] == (COUNT,)


@pytest.mark.parametrize(
    "data, message",
    [
        ({"NE(0625/0451)": {"z": float("nan")}}, "not a finite number"),
        ({"NE(0625/0451)": {"z": "high"}}, "not a number"),
        ({"NE(0625/0451)": {"z": True}}, "not a number"),
        ({"NE(0625/0451)": {"z": 1e39}}, "overflows f32"),
        ({"NE(0625/0451)": {"n": 1.0}}, "not an integer"),
        ({"NE(0625/0451)": {"k": -1}}, "outside u8"),
        ({"NE(0625/0451)": {"q": 1}}, "is a level variable and goes in observed"),
        ({"NE(0625/0451)": {"nope": 1}}, "no variable 'nope'"),
        ({"NE(0625/0451)": [1]}, "mapping of variable names"),
        ("NE(0625/0451)", "mapping of cells to values"),
        ([("NE(0625/0451)",)], "pairs of a cell and its values"),
        ([("NE(0625/0451)", {}), (itacart.cell_to_int("NE(0625/0451)"), {})], "twice"),
        (
            {"NE(0625/0451)": {"n": 1}, "NE(0625/0451(1))": {"n": 1}},
            "has daughters whose aggregate it holds",
        ),
        ({}, "needs at least one cell"),
    ],
)
def test_the_data_is_checked(data: Any, message: str) -> None:
    variables = [
        DataVariable("n", "u32", op="sum"),
        DataVariable("z", "f32", op="max"),
        DataVariable("k", "u8", op="max"),
        DataVariable("q", "u8", species="level"),
    ]
    with pytest.raises(SerializationError, match=message):
        encode_data_tree(data, variables)


def test_observations_are_checked_like_data() -> None:
    q = DataVariable("q", "u8", species="level")
    n = DataVariable("n", "u32", op="sum")
    with pytest.raises(
        SerializationError, match="is a pyramid variable and goes in data"
    ):
        encode_data_tree({}, [n, q], observed={"NE(0625/0451)": {"n": 1}})
    only = encode_data_tree({}, [n, q], observed={"NE(0625/0451(1))": {"q": 4}})
    assert decode_data_tree(only)["NE(0625/0451(1))"] == {"n": None, "q": 4}
    assert data_tree_schema(only)["resolution"] == 2


@pytest.mark.parametrize(
    "cell, error",
    [
        (DENIED_WEST, NonExistentCellError),
        ("NE(0625/0451(9))", InvalidIndexError),
        ("NE(0625/0451(1,2))", NonAtomicIndexError),
        (2.0, InvalidIndexError),
        (None, InvalidIndexError),
        (OFFSETS[1] + 2004 * 1001, NonExistentCellError),
    ],
)
def test_the_cells_are_checked(cell: Any, error: type[Exception]) -> None:
    with pytest.raises(error):
        encode_data_tree({cell: {"n": 1}}, [COUNT])


@pytest.mark.parametrize("resolution", [-1, 14, True, 2.0, "7"])
def test_the_resolution_is_checked(resolution: Any) -> None:
    with pytest.raises(ResolutionError, match="observation resolution is 0 to 13"):
        encode_data_tree({"NE(0625/0451)": {"n": 1}}, [COUNT], resolution=resolution)


def test_a_cell_finer_than_the_resolution_is_refused() -> None:
    with pytest.raises(
        ResolutionError, match="finer than the observation resolution 1"
    ):
        encode_data_tree({"NE(0625/0451(1))": {"n": 1}}, [COUNT], resolution=1)


def test_a_whole_quadrant_can_be_the_data() -> None:
    blob = encode_data_tree({"SE": {"n": 7}}, [DataVariable("n", "u8", op="sum")])
    assert decode_data_tree(blob) == {"SE": {"n": 7}}
    assert read_data_path(blob, "SE") == [("SE", {"n": 7})]
    assert data_tree_schema(blob)["nodes"] == (1,)
    assert read_data_cell(blob, "NE") == {"n": None}


# --------------------------------------------------------------------------
# Malformed blobs
# --------------------------------------------------------------------------

SMALL = (
    DataVariable("n", "u32", op="sum", absent="zero", leaf_max=1),
    DataVariable("h", "u16", op="max"),
    DataVariable("q", "u8", species="level"),
)


@functools.lru_cache(maxsize=None)
def small_blob() -> bytes:
    data = {}
    for base in ("NE(0625/0451)", "NE(0626/0451)", "SW(0100/0200)"):
        for k, child in enumerate(itacart.get_children(base, flatten=True)):
            data[child] = {"n": 1, "h": 10 * k}
    return encode_data_tree(data, SMALL, observed={"NE(0625/0451)": {"q": 5}})


def _parts(blob: bytes) -> dict[str, Any]:
    reader = dtb._Blob(blob)
    sections = []
    for offset, fmt, count in sorted(reader._specs.values()):
        sections.append(bytes(blob[offset : offset + struct.calcsize(fmt) * count]))
    return {
        "variables": reader.variables,
        "resolution": reader.resolution,
        "deepest": reader.deepest,
        "counts": list(reader.counts),
        "widths": {k: list(v) for k, v in reader.widths.items()},
        "flags": int(reader.pruned),
        "sections": sections,
    }


def _rebuilt(**changes: Any) -> bytes:
    parts = {**_parts(small_blob()), **changes}
    return dtb._frame(
        parts["variables"],
        parts["resolution"],
        parts["deepest"],
        parts["counts"],
        parts["widths"],
        parts["flags"],
        parts["sections"],
    )


def _section(key: Any) -> tuple[int, int]:
    offset, fmt, count = dtb._Blob(small_blob())._specs[key]
    return offset, struct.calcsize(fmt) * count


def _patched(key: Any, index: int, value: int, fmt: str = "B") -> bytes:
    blob = bytearray(small_blob())
    offset, _ = _section(key)
    struct.pack_into("<" + fmt, blob, offset + index * struct.calcsize(fmt), value)
    return bytes(blob)


def test_the_small_blob_is_what_the_corruptions_start_from() -> None:
    blob = small_blob()
    assert _rebuilt() == blob
    assert data_tree_schema(blob)["nodes"] == (2, 3, 12)
    assert len(decode_data_tree(blob)) == 17


@pytest.mark.parametrize(
    "blob, error, message",
    [
        (b"", MalformedBlobError, "not a DataTreeBlob"),
        (b"\xc7" + bytes(40), MalformedBlobError, "not a DataTreeBlob"),
        (b"\xc9\x20" + bytes(40), IncompatibleProfileError, "version 2"),
        (b"\xc9\x12" + bytes(40), IncompatibleProfileError, "flags 0x2"),
        (b"\xc9\x10\x01\x02\x01" + bytes(40), MalformedBlobError, "deepest level 2"),
        (b"\xc9\x10\x02\x02\x01\x05", MalformedBlobError, "truncated or malformed"),
        (b"\xc9\x10\x02\x02\x01\x01\xff" + bytes(40), MalformedBlobError, "truncated"),
    ],
)
def test_a_malformed_frame_is_refused(
    blob: bytes, error: type[Exception], message: str
) -> None:
    with pytest.raises(error, match=message):
        data_tree_schema(blob)


def _frame_with(variable_bytes: bytes) -> bytes:
    """A frame whose single variable is written by hand."""
    return b"\xc9\x10\x01\x01\x01" + variable_bytes + bytes(64)


@pytest.mark.parametrize(
    "variable, message",
    [
        (
            b"\x01n" + bytes([99, 0, 1, 1, 0xFF]) + bytes(8) + b"\x03\x00\x00\x00",
            "codes",
        ),
        (
            b"\x01n" + bytes([0, 1, 1, 1, 0xFF]) + bytes(8) + b"\x01\x00\x00",
            "invalid schema",
        ),
        (b"\x01n" + bytes([0, 1, 0, 1, 0xFF]) + bytes(8) + b"\x04\x00\x00", "levels"),
        (
            b"\x01n" + bytes([0, 0, 1, 1, 0xFF]) + bytes(8) + b"\x03\x00\x63\x00",
            "a width",
        ),
        (
            b"\x01n" + bytes([0, 0, 1, 1, 0xFF]) + bytes(8) + b"\x01\x00\x00",
            "misses a level",
        ),
        (
            b"\x01n" + bytes([0, 0, 1, 1, 0xFF]) + bytes(8) + b"\x03\x00\x00\x01",
            "widths of",
        ),
        (
            b"\x01n"
            + bytes([0, 0, 1, 1, 0xFF])
            + struct.pack("<Q", 1)
            + b"\x03\x00\x00\x04",
            "widths of",
        ),
    ],
)
def test_a_malformed_schema_is_refused(variable: bytes, message: str) -> None:
    with pytest.raises(MalformedBlobError, match=message):
        data_tree_schema(_frame_with(variable))


def test_the_sections_are_checked_against_the_frame() -> None:
    with pytest.raises(MalformedBlobError, match="a level without nodes"):
        data_tree_schema(_rebuilt(counts=[2, 0, 12]))
    with pytest.raises(MalformedBlobError, match="sections"):
        data_tree_schema(_rebuilt(sections=_parts(small_blob())["sections"][:-1]))
    with pytest.raises(MalformedBlobError, match="a section at"):
        blob = bytearray(small_blob())
        table = blob.index(struct.pack("<I", _section(("q", 0))[0]))
        struct.pack_into("<I", blob, table, _section(("q", 0))[0] + 8)
        data_tree_schema(bytes(blob))
    with pytest.raises(MalformedBlobError, match="bytes, expected"):
        data_tree_schema(small_blob() + bytes(8))


@pytest.mark.parametrize(
    "blob, message",
    [
        (lambda: _patched(("q", 0), 1, 0), "the quadrant codes"),
        (lambda: _patched(("first", 0), 1, 9, "I"), "the directory of level 0"),
        (lambda: _patched(("p1", 0), 1, 0, "I"), "the resolution-1 positions"),
        (lambda: _patched(("first", 0), 1, 3, "I"), "under the wrong quadrant"),
        (lambda: _patched(("mask", 1), 0, 0b0111), "the mask of node 0"),
        (lambda: _patched(("cell", 1), 0, 0), "a cell flag"),
        (lambda: _patched(("value", "h", 2), 3, 3, "H"), "is not the max"),
        (lambda: _patched(("value", "n", 1), 0, 5), "is not the sum"),
        (lambda: _patched(("value", "n", 2), 0, 2), "outside 0 to 1 per cell"),
        (lambda: _patched(("value", "q", 1), 1, 9), "an unknown slot"),
    ],
)
def test_the_decoder_checks_every_invariant(blob: Any, message: str) -> None:
    with pytest.raises(MalformedBlobError, match=message):
        decode_data_tree(blob())


def test_a_read_through_a_broken_directory_is_refused() -> None:
    broken = _patched(("first", 1), 0, 200, "I")
    cell = sorted(
        c for c in decode_data_tree(small_blob()) if itacart.get_resolution(c) == 2
    )[0]
    with pytest.raises(MalformedBlobError):
        read_data_cell(broken, cell)
    with pytest.raises(MalformedBlobError):
        read_data_path(broken, cell)
    with pytest.raises(MalformedBlobError):
        slice_data_tree(broken, cell)
    with pytest.raises(MalformedBlobError):
        update_data_tree(bytearray(broken), cell, {"h": 1})


def test_the_blob_is_any_contiguous_buffer() -> None:
    import array

    blob = small_blob()
    words = array.array("I")
    words.frombytes(blob)
    assert decode_data_tree(words) == decode_data_tree(blob)
    assert decode_data_tree(memoryview(blob)) == decode_data_tree(blob)
    with pytest.raises(SerializationError, match="bytes-like"):
        decode_data_tree("blob")  # type: ignore[arg-type]
    with pytest.raises(SerializationError, match="contiguous"):
        decode_data_tree(memoryview(blob + blob)[::2])


def test_a_big_endian_host_reads_the_same(monkeypatch: pytest.MonkeyPatch) -> None:
    """The format is little-endian; a host that is not reads through struct."""
    blob = small_blob()
    cell = sorted(c for c in decode_data_tree(blob) if itacart.get_resolution(c) == 2)[
        5
    ]
    expected = (
        decode_data_tree(blob),
        read_data_path(blob, cell),
        slice_data_tree(blob, cell),
    )
    monkeypatch.setattr(dtb, "_LITTLE_ENDIAN_HOST", False)
    copy = bytes(bytearray(blob))  # a new object, parsed afresh
    assert (
        decode_data_tree(copy),
        read_data_path(copy, cell),
        slice_data_tree(copy, cell),
    ) == expected
    target = bytearray(blob)
    update_data_tree(target, cell, {"h": 77})
    monkeypatch.setattr(dtb, "_LITTLE_ENDIAN_HOST", True)
    little = bytearray(blob)
    update_data_tree(little, cell, {"h": 77})
    assert target == little
    with pytest.raises(IndexError):
        dtb._Swapped(memoryview(blob), 0, "B", 2)[2]


# --------------------------------------------------------------------------
# Exact arithmetic
# --------------------------------------------------------------------------


def test_the_f32_rounding_is_correct_and_ties_to_even() -> None:
    rnd = random.Random(19)
    for _ in range(2000):
        mantissa = rnd.getrandbits(53) | 1 << 52
        exponent = rnd.randint(-160, 120)
        value = math.ldexp(mantissa, exponent - 52)
        fraction = Fraction(value)
        try:
            expected = _f32(value)
        except OverflowError:
            with pytest.raises(SerializationError, match="overflows f32"):
                dtb._to_f32(fraction.numerator, fraction.denominator)
            continue
        assert dtb._to_f32(fraction.numerator, fraction.denominator) == expected
        assert dtb._to_f32(-fraction.numerator, fraction.denominator) == -expected
    tie = Fraction(1) + Fraction(1, 2**24)  # halfway between 1 and the next f32
    assert dtb._to_f32(tie.numerator, tie.denominator) == 1.0
    above = Fraction(1) + Fraction(3, 2**24)
    assert dtb._to_f32(above.numerator, above.denominator) == 1.0 + 2**-22
    assert dtb._to_f32(0, 7) == 0.0
    assert dtb._to_f32(1, 2**151) == 0.0
    assert dtb._to_f32(3, 2**150) == 2.0**-148  # a tie in the subnormals, to even


def test_the_other_types_round_once_to_even() -> None:
    assert [dtb._rounded("i32", Fraction(k, 2)) for k in (5, 7, -5, -7, 3)] == [
        2,
        4,
        -2,
        -4,
        2,
    ]
    assert dtb._rounded("f64", Fraction(1, 3)) == 1 / 3
    with pytest.raises(SerializationError, match="overflows f64"):
        dtb._rounded("f64", Fraction(10**400))
    none = DataVariable("x", "u8", op="none", absent="zero")
    assert dtb._aggregate(none, [1, 2], [1, 1]) == 0
    mean = DataVariable("m", "f64", op="mean", weight="w")
    assert dtb._aggregate(mean, [1.0, 2.0], [0, None]) is None


def test_a_slice_of_a_quadrant_or_of_a_shallow_leaf() -> None:
    block = "NE(0625/0451(1))"
    deep = list(itacart.get_descendants("NE(0625/0451(2))", 4))
    blob = encode_data_tree(
        {block: {"h": 2.5}, **{leaf: {"h": float(k)} for k, leaf in enumerate(deep)}},
        [TALL],
    )
    assert decode_data_tree(slice_data_tree(blob, "NE")) == decode_data_tree(blob)
    shallow = decode_data_tree(slice_data_tree(blob, block))
    assert list(shallow) == ["NE", "NE(0625/0451)", block]
    assert shallow[block] == {"h": 2.5} and shallow["NE"] == {"h": 2.5}
    assert read_data_cell(blob, deep[0].replace("(2(", "(1(", 1)) == {"h": 2.5}


def test_reads_beside_the_data() -> None:
    blob = dense_blob()
    root = dense_case()[0]
    assert read_data_path(blob, "SW(0471/0256)") == [
        ("SW", decode_data_tree(blob)["SW"])
    ]
    near = next(
        c
        for c in itacart.get_children(itacart.get_parent(root), flatten=True)
        if c != root
    )
    beside = next(iter(itacart.get_descendants(near, 5)))
    assert read_data_cell(blob, beside) == {
        "n": 0,
        "z": None,
        "cls": None,
        "qual": None,
    }


def test_an_f32_that_rounds_past_the_largest_overflows() -> None:
    with pytest.raises(SerializationError, match="overflows f32"):
        dtb._to_f32(2**129, 1)
    almost = Fraction(2**128) - Fraction(2**103)  # rounds up to 2**128
    with pytest.raises(SerializationError, match="overflows f32"):
        dtb._to_f32(almost.numerator, almost.denominator)


def test_the_sizes_the_specification_quotes() -> None:
    """Every size on the specification page, in bits per leaf of ten thousand."""
    import zlib

    def split(blob: bytes) -> tuple[float, ...]:
        sizes: dict[str, int] = {}
        for key, (_, fmt, count) in dtb._Blob(blob)._specs.items():
            size = struct.calcsize(fmt) * count
            sizes[key[0]] = sizes.get(key[0], 0) + size + -size % 8
        structure = sum(sizes.get(k, 0) for k in ("q", "p1", "first", "mask", "cell"))
        return tuple(
            round(8 * part / 10_000, 1)
            for part in (structure, sizes["value"], sizes["valid"], len(blob))
        )

    assert len(dense_blob()) == 70_808
    assert split(dense_blob()) == (3.1, 50.9, 2.3, 56.6)
    assert split(drone_blob()) == (3.2, 50.7, 4.4, 58.9)
    assert split(drone_blob(True)) == (1.2, 6.4, 0.7, 8.9)
    leaves = itacart.compose(dense_case()[1])
    assert round(8 * len(itacart.encode_tree(leaves)) / 10_000, 1) == 8.4
    assert 30 < len(drone_blob()) / len(zlib.compress(drone_blob(), 9)) < 35
    assert 1.9 < len(dense_blob()) / len(zlib.compress(dense_blob(), 9)) < 2.1
