"""DataTreeBlob: values on a tree of cells, read and updated in place.

A TreeBlob holds a set of cells and answers *where*. A DataTreeBlob holds a
tree of cells and, at every node of it, the values of a fixed set of
variables, so it answers *what is there*: at a cell and at every one of its
ancestors, in one walk down the tree, on the bytes, without decoding the
rest. The two formats live side by side. The TreeBlob is the canonical,
content-addressed form of a region; this one is the fast, hierarchical form
of data about cells, and it has no canonical form to hash.

**The tree is lexical**, like the ancestry of the index. Its nodes are the
cells written into it and every prefix of their spellings. Away from the
domain border that is the physical tree. On the border, the surplus
children of a cell are spelled under the next column's prefix, which names
no cell, so they hang under that prefix and not under their physical
parent; inside a polar cap the same happens at intermediate levels. Such a
prefix is a node of the blob, because the walk needs it, and it holds the
aggregate of what hangs under it. No function here emits it, and a path
through it skips it. The aggregate of a border cell therefore covers its
own stem only, a share of the cell that depends on how deep the data goes.

**Variables.** A blob carries up to 255 fixed-width columns, declared with
:class:`DataVariable`. A *pyramid* variable is observed at the leaves and
held at every internal node as an operator over the daughters -- their
sum, mean, maximum, minimum or mode -- computed exactly and rounded once to
the variable's type. A *level* variable is observed wherever the caller
gives it a value, at a parent and at its daughter alike, and is never
aggregated. A mean or a mode is weighted by an integer sum declared before
it, normally the count of the cells where the value was observed, so a
node's mean is over what was observed rather than over what could have
been.

**Absence.** Each variable declares what a missing value means: zero, or
unknown. An unknown value is a cleared bit in a validity bitmap beside the
values, so the type keeps its whole domain. Operators skip the unknown
daughters, and a node none of whose daughters is known is unknown itself.
Operators run over the daughters the tree holds: a daughter it does not
hold is not a zero.

**Blocks.** The blob declares an observation resolution, the finest of its
data. A leaf above it stands for every cell beneath it: reading inside it
returns its values, with a sum divided evenly over the cells. Such a leaf
must be a cell that does not absorb the border, so that the cells beneath
it are exactly the product of the refinement ratios.
``encode_data_tree(..., prune=True)`` makes blocks out of complete
families of leaves that agree on every variable, and out of nothing else.

**Widths.** An unsigned sum declared with ``leaf_max`` is written at each
level in the narrowest unsigned type its structure can fill, so a count
costs one byte at the leaves and two or four higher up, and an update in
place can never overflow it.

**Reading.** Every read works on the bytes. The frame of a ``bytes`` blob,
which cannot change, is parsed on the first read and kept for the next
ones -- the last sixteen blobs read are kept -- while a ``bytearray`` is
parsed on every call, because it can change between them. A read handed a
compositional index, or several integers, answers for all of them on one
parse.

The wire format is specified on the DataTreeBlob encoding page of the
documentation.
"""

from __future__ import annotations

import math
import struct
import sys
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
from typing import Any, Iterable, Iterator, Mapping, NamedTuple, Sequence, Union

from ..boundary import absorbs_border, is_valid_cell
from ..constants import MAX_RESOLUTION, QUADRANTS
from ..exceptions import (
    IncompatibleProfileError,
    InvalidIndexError,
    MalformedBlobError,
    NonAtomicIndexError,
    NonExistentCellError,
    ResolutionError,
    SerializationError,
)
from ..index import iter_cells, split_components
from ..integer_index import (
    BOX_COLUMNS,
    BOX_ROWS,
    OFFSETS,
    RATIOS,
    _as_integer,
    _checked,
    _position,
    _resolution,
    _spelling,
)

__all__ = [
    "DataVariable",
    "encode_data_tree",
    "decode_data_tree",
    "read_data_path",
    "read_data_cell",
    "slice_data_tree",
    "update_data_tree",
    "data_tree_schema",
]

#: First byte of every DataTreeBlob: the ``0xC_`` family, after the
#: TreeBlob's ``0xC7`` and the GeometryBlob's ``0xC8``.
MAGIC = 0xC9

#: Version written in the high nibble of the second byte.
VERSION = 1

#: Flag bit of the low nibble: the writer pruned complete families into
#: blocks. Informational; the other three bits are reserved and zero.
FLAG_PRUNED = 0x1

#: The fixed-width types, in the order of their codes on the wire.
DTYPES: dict[str, str] = {
    "u8": "B",
    "u16": "H",
    "u32": "I",
    "u64": "Q",
    "i8": "b",
    "i16": "h",
    "i32": "i",
    "i64": "q",
    "f32": "f",
    "f64": "d",
}
_DTYPE_CODES = tuple(DTYPES)
_UNSIGNED = ("u8", "u16", "u32", "u64")

SPECIES = ("pyramid", "level")
OPERATORS = ("none", "sum", "mean", "max", "min", "mode")
ABSENCE = ("unknown", "zero")

#: Every finite double is a whole multiple of ``2**-1074``, so scaling by
#: ``2**1074`` turns a sum of doubles into a sum of integers, exactly.
_DYADIC = 1074

_FLT_MAX = float.fromhex("0x1.fffffep+127")

#: Positions of the resolution-1 box, the four quadrants included.
_LEVEL1_BOX = len(QUADRANTS) * BOX_COLUMNS * BOX_ROWS

#: A memoryview cast reads in the host's byte order, and the format is
#: little-endian. A big-endian host reads through ``_Swapped`` instead.
_LITTLE_ENDIAN_HOST = sys.byteorder == "little"

Value = Union[int, float, None]


class DataVariable(NamedTuple):
    """One column of a DataTreeBlob.

    Args:
        name: Unique name, 1 to 255 bytes of UTF-8. It lives in the schema;
            no text is stored per cell.
        dtype: One of ``u8``, ``u16``, ``u32``, ``u64``, ``i8``, ``i16``,
            ``i32``, ``i64``, ``f32``, ``f64``. A categorical variable is an
            integer code.
        species: ``"pyramid"``, observed at the leaves and aggregated at
            every internal node, or ``"level"``, observed at any node and
            never aggregated.
        op: For a pyramid variable, ``"sum"``, ``"mean"``, ``"max"``,
            ``"min"``, ``"mode"`` (ties go to the smallest value) or
            ``"none"``, which leaves internal nodes without a value. A level
            variable takes ``"none"``.
        absent: What a missing value means: ``"unknown"`` or ``"zero"``.
        weight: Required for ``"mean"`` and ``"mode"``, refused otherwise:
            the name of an integer pyramid sum declared earlier, whose
            value at each daughter weights it. Normally the count of the
            cells where the averaged value was observed.
        leaf_max: Optional, for an unsigned pyramid sum: the largest value
            one observation cell holds. With it, each level is written in
            the narrowest unsigned type its structure can fill.

    Example:
        >>> DataVariable("height", "u16", op="max").op
        'max'
    """

    name: str
    dtype: str
    species: str = "pyramid"
    op: str = "none"
    absent: str = "unknown"
    weight: str | None = None
    leaf_max: int | None = None


# --------------------------------------------------------------------------
# Exact arithmetic
# --------------------------------------------------------------------------


def _range(dtype: str) -> tuple[int, int]:
    bits = 8 * struct.calcsize(DTYPES[dtype])
    if dtype in _UNSIGNED:
        return 0, (1 << bits) - 1
    return -(1 << (bits - 1)), (1 << (bits - 1)) - 1


def _is_float(dtype: str) -> bool:
    return dtype[0] == "f"


def _zero(dtype: str) -> int | float:
    return 0.0 if _is_float(dtype) else 0


def _absent(variable: DataVariable) -> Value:
    return None if variable.absent == "unknown" else _zero(variable.dtype)


def _dyadic(value: float) -> int:
    """``value * 2**1074`` as an exact integer."""
    numerator, denominator = value.as_integer_ratio()
    return numerator << (_DYADIC - (denominator.bit_length() - 1))


def _exact(value: int | float) -> Fraction:
    if isinstance(value, float):
        return Fraction(_dyadic(value), 1 << _DYADIC)
    return Fraction(value)


def _to_f32(numerator: int, denominator: int) -> float:
    """The binary32 value nearest to ``numerator / denominator``, ties to even."""
    if numerator == 0:
        return 0.0
    sign = -1.0 if numerator < 0 else 1.0
    numerator = abs(numerator)

    def scaled(exponent: int) -> tuple[int, int, int]:
        if exponent >= 0:
            divisor = denominator << exponent
            return (*divmod(numerator, divisor), divisor)
        return (*divmod(numerator << -exponent, denominator), denominator)

    # The first guess puts the quotient in [2**23, 2**25); one step fixes it.
    exponent = numerator.bit_length() - denominator.bit_length() - 24
    if scaled(exponent)[0] >= 1 << 24:
        exponent += 1
    exponent = max(exponent, -149)
    quotient, rest, divisor = scaled(exponent)
    if 2 * rest > divisor or (2 * rest == divisor and quotient & 1):
        quotient += 1
    value = math.ldexp(quotient, exponent)
    if value > _FLT_MAX:
        raise SerializationError(f"{sign * value:g} overflows f32")
    return sign * value


def _rounded(dtype: str, exact: Fraction) -> int | float:
    """An exact rational, rounded once to a variable's type, ties to even."""
    if dtype == "f32":
        return _to_f32(exact.numerator, exact.denominator)
    if dtype == "f64":
        try:
            # Integer true division is correctly rounded.
            return exact.numerator / exact.denominator
        except OverflowError as exc:
            raise SerializationError("the result overflows f64") from exc
    quotient, rest = divmod(exact.numerator, exact.denominator)
    if 2 * rest > exact.denominator or (2 * rest == exact.denominator and quotient & 1):
        quotient += 1
    return quotient


def _aggregate(
    variable: DataVariable, values: Sequence[Value], weights: Sequence[Value]
) -> Value:
    """A pyramid variable's operator over its daughters' values.

    Unknown daughters are skipped, and so are the daughters a mean or a mode
    gives no weight; with none left, the node takes the variable's absence.
    """
    if variable.op == "none":
        return _absent(variable)
    known: list[tuple[int | float, int]] = []
    for value, weight in zip(values, weights):
        if value is None:
            continue
        units = int(weight or 0)
        if variable.op in ("mean", "mode") and units <= 0:
            continue
        known.append((value, units))
    if not known:
        return _absent(variable)
    if variable.op == "max":
        return max(value for value, _ in known)
    if variable.op == "min":
        return min(value for value, _ in known)
    if variable.op == "mode":
        tally: dict[int | float, int] = {}
        for value, units in known:
            tally[value] = tally.get(value, 0) + units
        best = max(tally.values())
        return min(value for value, count in tally.items() if count == best)
    if variable.op == "sum":
        if not _is_float(variable.dtype):
            return sum(int(value) for value, _ in known)
        total = sum(_dyadic(float(value)) for value, _ in known)
        return _rounded(variable.dtype, Fraction(total, 1 << _DYADIC))
    weight_total = sum(units for _, units in known)
    exact = sum((_exact(value) * units for value, units in known), Fraction(0))
    return _rounded(variable.dtype, exact / weight_total)


# --------------------------------------------------------------------------
# Existence, measured once per node of the grid
# --------------------------------------------------------------------------


@lru_cache(maxsize=1 << 16)
def _status(level: int, position: int) -> int:
    """0 names no cell; 1 a cell that absorbs the border; 2 a regular cell.

    Every lexical descendant of a regular cell is a cell, which is what lets
    one measurement stand for a whole subtree. The facts are the grid's and
    never change, so they are cached.
    """
    spelling = _spelling(OFFSETS[level] + position)
    if not is_valid_cell(spelling):
        return 0
    return 1 if absorbs_border(spelling) else 2


def _exists(level: int, position: int) -> bool:
    """Whether the node at ``(level, position)`` names a cell."""
    if level == 0:
        return True
    chain = _chain(level, position)
    for ancestor in range(1, level):
        if _status(ancestor, chain[ancestor]) == 2:
            return True
    return _status(level, position) != 0


def _regular(level: int, position: int) -> bool:
    """A cell whose descendants are the full product of the ratios."""
    return level > 0 and _status(level, position) == 2


def _chain(level: int, position: int) -> list[int]:
    """Box positions of a node and its lexical ancestors, resolution 0 first."""
    out = [0] * (level + 1)
    out[level] = position
    for finer in range(level, 0, -1):
        position //= RATIOS[finer]
        out[finer - 1] = position
    return out


def _cells_between(coarse: int, fine: int) -> int:
    """Cells of resolution ``fine`` under one regular cell of ``coarse``."""
    count = 1
    for finer in range(coarse + 1, fine + 1):
        count *= RATIOS[finer]
    return count


# --------------------------------------------------------------------------
# Cells and values as the caller hands them
# --------------------------------------------------------------------------


def _cell_key(cell: Any) -> tuple[int, int, str]:
    """Resolution, box position and spelling of one existing cell."""
    if isinstance(cell, str):
        resolution, position = _position(split_components(cell))
        spelling = cell
    elif _as_integer(cell) is not None or isinstance(cell, bool):
        value = _checked(cell)
        resolution = _resolution(value)
        position = value - OFFSETS[resolution]
        spelling = _spelling(value)
    else:
        raise InvalidIndexError(
            f"a cell is a spelling or an integer, not {type(cell).__name__}"
        )
    if not _exists(resolution, position):
        raise NonExistentCellError(f"index names no cell: {spelling}")
    return resolution, position, spelling


def _cells(cell: Any) -> tuple[list[tuple[int, int, str]], bool]:
    """The cells a read was handed, and whether it was handed exactly one."""
    if isinstance(cell, str):
        if "," not in cell:
            return [_cell_key(cell)], True
        return [_cell_key(atom) for atom in iter_cells(cell)], False
    if _as_integer(cell) is not None or isinstance(cell, bool):
        return [_cell_key(cell)], True
    if isinstance(cell, (bytes, bytearray)) or not isinstance(cell, Iterable):
        raise InvalidIndexError(
            f"a cell is a spelling or an integer, not {type(cell).__name__}"
        )
    return [_cell_key(item) for item in cell], False


def _one_cell(cell: Any) -> tuple[int, int, str]:
    if isinstance(cell, str) and "," in cell:
        raise NonAtomicIndexError(f"{cell!r} addresses more than one cell")
    return _cell_key(cell)


def _leaf_value(variable: DataVariable, value: Any, where: str) -> Value:
    """A value as the caller gave it, checked and brought to its type."""
    if value is None:
        return _absent(variable)
    if _is_float(variable.dtype):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise SerializationError(
                f"{variable.name} at {where}: {value!r} is not a number"
            )
        number = float(value)
        if not math.isfinite(number):
            raise SerializationError(
                f"{variable.name} at {where}: {value!r} is not a finite number"
            )
        if variable.dtype == "f64":
            return number
        try:
            return float(struct.unpack("<f", struct.pack("<f", number))[0])
        except OverflowError as exc:
            raise SerializationError(
                f"{variable.name} at {where}: {value!r} overflows f32"
            ) from exc
    integer = _as_integer(value)
    if integer is None:
        raise SerializationError(
            f"{variable.name} at {where}: {value!r} is not an integer"
        )
    low, high = _range(variable.dtype)
    if variable.leaf_max is None and not low <= integer <= high:
        raise SerializationError(
            f"{variable.name} at {where}: {integer} is outside {variable.dtype}"
        )
    return integer


def _pairs(source: Any, label: str) -> Iterator[tuple[Any, Any]]:
    if isinstance(source, Mapping):
        yield from source.items()
        return
    if isinstance(source, (str, bytes, bytearray)) or not isinstance(source, Iterable):
        raise SerializationError(
            f"{label} is a mapping of cells to values, or pairs of them, "
            f"not {type(source).__name__}"
        )
    for entry in source:
        if not (isinstance(entry, (tuple, list)) and len(entry) == 2):
            raise SerializationError(
                f"{label} holds pairs of a cell and its values, not {entry!r}"
            )
        yield entry[0], entry[1]


# --------------------------------------------------------------------------
# The schema
# --------------------------------------------------------------------------


def _check_variables(variables: Iterable[DataVariable]) -> tuple[DataVariable, ...]:
    if isinstance(variables, DataVariable):
        variables = (variables,)
    declared = tuple(variables)
    if not 1 <= len(declared) <= 255:
        raise SerializationError(
            f"a DataTreeBlob carries 1 to 255 variables, not {len(declared)}"
        )
    seen: dict[str, DataVariable] = {}
    for variable in declared:
        if not isinstance(variable, DataVariable):
            raise SerializationError(
                f"variables are DataVariable, not {type(variable).__name__}"
            )
        name = variable.name
        if not isinstance(name, str) or not 1 <= len(name.encode()) <= 255:
            raise SerializationError(f"variable name {name!r} is not 1 to 255 bytes")
        if name in seen:
            raise SerializationError(f"variable {name!r} is declared twice")
        if variable.dtype not in DTYPES:
            raise SerializationError(
                f"{name}: type {variable.dtype!r} is not one of {', '.join(DTYPES)}"
            )
        if variable.species not in SPECIES:
            raise SerializationError(
                f"{name}: species {variable.species!r} is not pyramid or level"
            )
        if variable.op not in OPERATORS:
            raise SerializationError(
                f"{name}: operator {variable.op!r} is not one of "
                f"{', '.join(OPERATORS)}"
            )
        if variable.species == "level" and variable.op != "none":
            raise SerializationError(
                f"{name}: a level variable is never aggregated, so its operator "
                f"is 'none', not {variable.op!r}"
            )
        if variable.absent not in ABSENCE:
            raise SerializationError(
                f"{name}: absence {variable.absent!r} is not unknown or zero"
            )
        _check_weight(variable, seen)
        _check_leaf_max(variable)
        seen[name] = variable
    return declared


def _check_weight(variable: DataVariable, seen: Mapping[str, DataVariable]) -> None:
    name = variable.name
    if variable.op not in ("mean", "mode"):
        if variable.weight is not None:
            raise SerializationError(f"{name}: only a mean or a mode takes a weight")
        return
    weight = seen.get(variable.weight) if isinstance(variable.weight, str) else None
    if weight is None:
        raise SerializationError(
            f"{name}: a {variable.op} needs a weight, the name of an integer "
            f"pyramid sum declared before it"
        )
    if weight.species != "pyramid" or weight.op != "sum" or _is_float(weight.dtype):
        raise SerializationError(
            f"{name}: weight {weight.name!r} is not an integer pyramid sum"
        )


def _check_leaf_max(variable: DataVariable) -> None:
    leaf_max = variable.leaf_max
    if leaf_max is None:
        return
    name = variable.name
    if not (
        variable.species == "pyramid"
        and variable.op == "sum"
        and variable.dtype in _UNSIGNED
    ):
        raise SerializationError(f"{name}: leaf_max is for an unsigned pyramid sum")
    if (
        isinstance(leaf_max, bool)
        or not isinstance(leaf_max, int)
        or not 1 <= leaf_max <= _range(variable.dtype)[1]
    ):
        raise SerializationError(
            f"{name}: leaf_max {leaf_max!r} is not a positive {variable.dtype}"
        )


# --------------------------------------------------------------------------
# The tree a writer builds
# --------------------------------------------------------------------------


class _Tree:
    """Nodes per level, daughters as runs, and every variable's values."""

    def __init__(
        self,
        variables: tuple[DataVariable, ...],
        resolution: int,
        levels: list[list[int]],
        values: dict[str, list[list[Value] | None]],
        pruned: bool,
    ) -> None:
        self.variables = variables
        self.resolution = resolution
        self.levels = levels
        self.values = values
        self.pruned = pruned
        self.deepest = len(levels) - 1
        self.runs = _runs(levels)

    def slot(self, level: int, position: int) -> int:
        return bisect_left(self.levels[level], position)

    def is_leaf(self, level: int, slot: int) -> bool:
        if level == self.deepest:
            return True
        first, last = self.runs[level][slot]
        return first == last

    def row(self, name: str, level: int) -> list[Value]:
        row = self.values[name][level]
        assert row is not None
        return row


def _runs(levels: Sequence[Sequence[int]]) -> list[list[tuple[int, int]]]:
    """``runs[k][i]``: the daughters of node ``i`` of level ``k``, in ``k + 1``.

    Every node has its parent in the level above, so the daughters of
    consecutive nodes are consecutive runs below.
    """
    runs = []
    for level in range(len(levels) - 1):
        ratio = RATIOS[level + 1]
        below = levels[level + 1]
        out = []
        cursor = 0
        for position in levels[level]:
            first = cursor
            while cursor < len(below) and below[cursor] < (position + 1) * ratio:
                cursor += 1
            out.append((first, cursor))
        runs.append(out)
    return runs


def _ingest(
    variables: tuple[DataVariable, ...], data: Any, observed: Any
) -> tuple[dict[tuple[int, int], dict[str, Value]], dict[tuple[int, int], str]]:
    """The caller's pairs, keyed by node, every cell and value checked."""
    by_name = {variable.name: variable for variable in variables}
    given: dict[tuple[int, int], dict[str, Value]] = {}
    spellings: dict[tuple[int, int], str] = {}
    for source, species, label in (
        (data, "pyramid", "data"),
        (observed, "level", "observed"),
    ):
        if source is None:
            continue
        seen: set[tuple[int, int]] = set()
        for cell, values in _pairs(source, label):
            resolution, position, spelling = _cell_key(cell)
            key = (resolution, position)
            if key in seen:
                raise SerializationError(f"{label} names {spelling} twice")
            seen.add(key)
            if not isinstance(values, Mapping):
                raise SerializationError(
                    f"{label} at {spelling}: the values are a mapping of variable "
                    f"names, not {type(values).__name__}"
                )
            entry = given.setdefault(key, {})
            spellings[key] = spelling
            for name, value in values.items():
                variable = by_name.get(name)
                if variable is None:
                    raise SerializationError(
                        f"{label} at {spelling}: there is no variable {name!r}"
                    )
                if variable.species != species:
                    other = "observed" if species == "pyramid" else "data"
                    raise SerializationError(
                        f"{label} at {spelling}: {name!r} is a {variable.species} "
                        f"variable and goes in {other}"
                    )
                entry[name] = _leaf_value(variable, value, spelling)
    return given, spellings


def _build(
    variables: tuple[DataVariable, ...],
    given: dict[tuple[int, int], dict[str, Value]],
    spellings: dict[tuple[int, int], str],
    resolution: Any,
    pruned: bool = False,
) -> _Tree:
    if not given:
        raise SerializationError("a DataTreeBlob needs at least one cell")
    pyramids = {v.name for v in variables if v.species == "pyramid"}
    deepest = max(level for level, _ in given)
    if resolution is None:
        observed = [
            level for (level, _), entry in given.items() if pyramids & set(entry)
        ]
        resolution = max(observed) if observed else deepest
    if (
        isinstance(resolution, bool)
        or not isinstance(resolution, int)
        or not 0 <= resolution <= MAX_RESOLUTION
    ):
        raise ResolutionError(
            f"the observation resolution is 0 to {MAX_RESOLUTION}, not {resolution!r}"
        )
    if deepest > resolution:
        finest = max(given, key=lambda key: key[0])
        raise ResolutionError(
            f"{spellings[finest]} is finer than the observation resolution "
            f"{resolution}"
        )
    nodes: list[set[int]] = [set() for _ in range(deepest + 1)]
    for level, position in given:
        for ancestor, ancestor_position in enumerate(_chain(level, position)):
            nodes[ancestor].add(ancestor_position)
    levels = [sorted(row) for row in nodes]
    values: dict[str, list[list[Value] | None]] = {}
    for variable in variables:
        if variable.species == "pyramid":
            values[variable.name] = [[None] * len(row) for row in levels]
            continue
        present = {
            level for (level, _), entry in given.items() if variable.name in entry
        }
        values[variable.name] = [
            [_absent(variable)] * len(row) if level in present else None
            for level, row in enumerate(levels)
        ]
    tree = _Tree(variables, resolution, levels, values, pruned)
    for (level, position), entry in given.items():
        slot = tree.slot(level, position)
        leaf = tree.is_leaf(level, slot)
        for name, value in entry.items():
            if name in pyramids and not leaf:
                raise SerializationError(
                    f"data at {spellings[(level, position)]}: {name!r} is a pyramid "
                    f"variable, and this node has daughters whose aggregate it holds"
                )
            tree.row(name, level)[slot] = value
    for level in range(deepest, -1, -1):
        for slot, position in enumerate(levels[level]):
            if tree.is_leaf(level, slot):
                _settle_leaf(tree, level, slot, given.get((level, position)))
                continue
            first, last = tree.runs[level][slot]
            for variable in variables:
                if variable.species == "pyramid":
                    aggregate = _aggregate(
                        variable,
                        tree.row(variable.name, level + 1)[first:last],
                        _weights(tree, variable, level + 1, first, last),
                    )
                    _check_sum(variable, aggregate, level, position)
                    tree.row(variable.name, level)[slot] = aggregate
    return tree


def _weights(
    tree: _Tree, variable: DataVariable, level: int, first: int, last: int
) -> list[Value]:
    if variable.weight is None:
        return [1] * (last - first)
    return tree.row(variable.weight, level)[first:last]


def _settle_leaf(
    tree: _Tree, level: int, slot: int, entry: dict[str, Value] | None
) -> None:
    """Absent pyramid values filled in, and a leaf's sums checked."""
    position = tree.levels[level][slot]
    for variable in tree.variables:
        if variable.species != "pyramid":
            continue
        row = tree.row(variable.name, level)
        if entry is None or variable.name not in entry:
            row[slot] = _absent(variable)
            continue
        _check_leaf(variable, row[slot], level, position, tree.resolution)


def _check_leaf(
    variable: DataVariable, value: Value, level: int, position: int, resolution: int
) -> None:
    """A leaf above the resolution is a block, and its sum spreads evenly."""
    if value is None:
        return
    where = _spelling(OFFSETS[level] + position)
    cells = _cells_between(level, resolution)
    if level < resolution and not _regular(level, position):
        raise SerializationError(
            f"{where} is a leaf above the observation resolution {resolution}, "
            f"so it stands for every cell beneath it; it absorbs the border, "
            f"where those cells are not a regular block"
        )
    if variable.op != "sum" or _is_float(variable.dtype):
        return
    assert isinstance(value, int)
    if value % cells:
        raise SerializationError(
            f"{variable.name} at {where}: a block of {cells} cells holds a sum "
            f"divided evenly over them, and {value} is not"
        )
    if variable.leaf_max is not None and not 0 <= value <= variable.leaf_max * cells:
        raise SerializationError(
            f"{variable.name} at {where}: {value} is outside 0 to {variable.leaf_max} "
            f"per cell over {cells} cell(s)"
        )


def _check_sum(variable: DataVariable, value: Value, level: int, position: int) -> None:
    """A sum of a declared type has to fit it; one with ``leaf_max`` always does."""
    if (
        value is None
        or variable.op != "sum"
        or variable.leaf_max is not None
        or _is_float(variable.dtype)
    ):
        return
    low, high = _range(variable.dtype)
    if not low <= int(value) <= high:
        raise SerializationError(
            f"{variable.name}: the sum {value} at "
            f"{_spelling(OFFSETS[level] + position)} overflows {variable.dtype}; "
            f"declare leaf_max, or a wider type"
        )


def _prune(tree: _Tree) -> dict[tuple[int, int], dict[str, Value]]:
    """Complete families of equal leaves, folded into blocks, as surviving pairs.

    Bottom-up, so a block can join its sisters into a larger one. A family
    folds only under a parent that is a regular cell, only when every
    daughter the ratio allows is present and is a leaf agreeing with its
    sisters on every variable, and never at a level where a level variable
    lives, whose values there would be lost.
    """
    deepest = tree.deepest
    leaf = [
        [tree.is_leaf(level, slot) for slot in range(len(row))]
        for level, row in enumerate(tree.levels)
    ]
    alive = [[True] * len(row) for row in tree.levels]
    levels_with_observations = {
        level
        for variable in tree.variables
        if variable.species == "level"
        for level, row in enumerate(tree.values[variable.name])
        if row is not None
    }
    for level in range(deepest - 1, -1, -1):
        if level + 1 in levels_with_observations:
            continue
        for slot, (first, last) in enumerate(tree.runs[level]):
            if (
                last - first == RATIOS[level + 1]
                and all(leaf[level + 1][first:last])
                and _regular(level, tree.levels[level][slot])
                and all(
                    _family_agrees(tree.row(v.name, level + 1)[first:last])
                    for v in tree.variables
                    if v.species == "pyramid"
                )
            ):
                for daughter in range(first, last):
                    alive[level + 1][daughter] = False
                leaf[level][slot] = True
    survivors: dict[tuple[int, int], dict[str, Value]] = {}
    for level, row in enumerate(tree.levels):
        for slot, position in enumerate(row):
            if alive[level][slot] and leaf[level][slot]:
                survivors[(level, position)] = {
                    v.name: tree.row(v.name, level)[slot]
                    for v in tree.variables
                    if v.species == "pyramid"
                }
    return survivors


def _family_agrees(values: Sequence[Value]) -> bool:
    head = values[0]
    return all(_same(value, head) for value in values[1:])


def _same(one: Value, other: Value) -> bool:
    if one is None or other is None:
        return one is None and other is None
    if isinstance(one, float) or isinstance(other, float):
        return struct.pack("<d", float(one)) == struct.pack("<d", float(other))
    return one == other


# --------------------------------------------------------------------------
# Writing the bytes
# --------------------------------------------------------------------------


def _mask_format(level: int) -> str:
    """Mask of a node at ``level``: one bit per daughter its ratio allows."""
    return "B" if RATIOS[level + 1] == 4 else "I"


def _pack(fmt: str, items: Sequence[Any]) -> bytes:
    return struct.pack(f"<{len(items)}{fmt}", *items)


def _bits(flags: Sequence[bool]) -> bytes:
    out = bytearray((len(flags) + 7) // 8)
    for index, flag in enumerate(flags):
        if flag:
            out[index >> 3] |= 1 << (index & 7)
    return bytes(out)


def _masks(tree: _Tree, level: int) -> list[int]:
    ratio = RATIOS[level + 1]
    below = tree.levels[level + 1]
    masks = []
    for position, (first, last) in zip(tree.levels[level], tree.runs[level]):
        mask = 0
        for daughter in range(first, last):
            mask |= 1 << (below[daughter] - position * ratio)
        masks.append(mask)
    return masks


def _widths(tree: _Tree) -> dict[str, list[str | None]]:
    """The type each variable is written in at each level, or ``None``."""
    cells: list[list[int]] = [[0] * len(row) for row in tree.levels]
    for level in range(tree.deepest, -1, -1):
        for slot in range(len(tree.levels[level])):
            if tree.is_leaf(level, slot):
                cells[level][slot] = _cells_between(level, tree.resolution)
            else:
                first, last = tree.runs[level][slot]
                cells[level][slot] = sum(cells[level + 1][first:last])
    widths: dict[str, list[str | None]] = {}
    for variable in tree.variables:
        row: list[str | None] = []
        for level in range(tree.deepest + 1):
            if tree.values[variable.name][level] is None:
                row.append(None)
            elif variable.leaf_max is None:
                row.append(variable.dtype)
            else:
                bound = variable.leaf_max * max(cells[level])
                fitting = [d for d in _UNSIGNED if bound <= _range(d)[1]]
                if not fitting:
                    raise SerializationError(
                        f"{variable.name}: {bound} at resolution {level} does not "
                        f"fit 64 bits"
                    )
                row.append(fitting[0])
        widths[variable.name] = row
    return widths


def _section_count(deepest: int, presence: Sequence[tuple[int, bool]]) -> int:
    """Sections a blob of this depth and this variable layout holds."""
    count = 1 + (deepest >= 1)
    count += sum(1 + 2 * (level >= 1) for level in range(deepest))
    return count + sum(1 + unknown for _, unknown in presence)


def _frame(
    variables: tuple[DataVariable, ...],
    resolution: int,
    deepest: int,
    counts: Sequence[int],
    widths: Mapping[str, Sequence[str | None]],
    flags: int,
    sections: Sequence[bytes],
) -> bytes:
    names = [variable.name for variable in variables]
    head = bytearray(
        [MAGIC, (VERSION << 4) | flags, resolution, deepest, len(variables)]
    )
    for variable in variables:
        name = variable.name.encode()
        head += bytes([len(name)]) + name
        head += bytes(
            [
                _DTYPE_CODES.index(variable.dtype),
                SPECIES.index(variable.species),
                OPERATORS.index(variable.op),
                ABSENCE.index(variable.absent),
                names.index(variable.weight) if variable.weight else 0xFF,
            ]
        )
        head += struct.pack("<Q", variable.leaf_max or 0)
        present = [level for level, width in enumerate(widths[variable.name]) if width]
        head += struct.pack("<H", sum(1 << level for level in present))
        head += bytes(
            _DTYPE_CODES.index(str(widths[variable.name][level])) for level in present
        )
    head += _pack("I", counts)
    head += struct.pack("<H", len(sections))
    table = len(head)
    head += bytes(4 * len(sections))
    body = bytearray()
    start = -len(head) % 8 + len(head)
    offsets = []
    for section in sections:
        offsets.append(start + len(body))
        body += section
        body += bytes(-len(body) % 8)
    struct.pack_into(f"<{len(sections)}I", head, table, *offsets)
    head += bytes(start - len(head))
    return bytes(head + body)


def _serialize(tree: _Tree) -> bytes:
    widths = _widths(tree)
    deepest = tree.deepest
    sections = [bytes(tree.levels[0])]
    if deepest >= 1:
        sections.append(_pack("I", tree.levels[1]))
    for level in range(deepest):
        runs = tree.runs[level]
        sections.append(_pack("I", [first for first, _ in runs] + [runs[-1][1]]))
        if level >= 1:
            sections.append(_pack(_mask_format(level), _masks(tree, level)))
            sections.append(
                _bits(
                    [
                        tree.is_leaf(level, slot) or _exists(level, position)
                        for slot, position in enumerate(tree.levels[level])
                    ]
                )
            )
    for level in range(deepest + 1):
        for variable in tree.variables:
            width = widths[variable.name][level]
            if width is None:
                continue
            row = tree.row(variable.name, level)
            zero = _zero(width)
            sections.append(
                _pack(
                    DTYPES[width], [zero if value is None else value for value in row]
                )
            )
            if variable.absent == "unknown":
                sections.append(_bits([value is not None for value in row]))
    return _frame(
        tree.variables,
        tree.resolution,
        deepest,
        [len(row) for row in tree.levels],
        widths,
        FLAG_PRUNED if tree.pruned else 0,
        sections,
    )


# --------------------------------------------------------------------------
# Reading the bytes
# --------------------------------------------------------------------------


class _Swapped:
    """A little-endian array read through ``struct``, for a big-endian host."""

    def __init__(self, view: memoryview, offset: int, fmt: str, count: int) -> None:
        self._view = view
        self._offset = offset
        self._fmt = "<" + fmt
        self._size = struct.calcsize(fmt)
        self._count = count

    def __len__(self) -> int:
        return self._count

    def __getitem__(self, index: int) -> Any:
        if not 0 <= index < self._count:
            raise IndexError(index)
        return struct.unpack_from(
            self._fmt, self._view, self._offset + index * self._size
        )[0]

    def __setitem__(self, index: int, value: Any) -> None:
        struct.pack_into(
            self._fmt, self._view, self._offset + index * self._size, value
        )


class _Sections:
    """One kind of section, indexed by level or by ``(name, level)``."""

    def __init__(self, blob: _Blob, kind: str) -> None:
        self._blob = blob
        self._kind = kind

    def _key(self, index: Any) -> tuple[Any, ...]:
        return (self._kind, *index) if isinstance(index, tuple) else (self._kind, index)

    def __getitem__(self, index: Any) -> Any:
        return self._blob.section(self._key(index))

    def get(self, index: Any) -> Any:
        key = self._key(index)
        return self._blob.section(key) if key in self._blob._specs else None


class _Blob:
    """The frame of a DataTreeBlob, parsed, and views on its sections."""

    def __init__(self, buffer: Any, writable: bool = False) -> None:
        try:
            view = memoryview(buffer)
        except TypeError as exc:
            raise SerializationError(
                f"a DataTreeBlob is bytes-like, not {type(buffer).__name__}"
            ) from exc
        if not view.c_contiguous:
            raise SerializationError("a DataTreeBlob is a contiguous buffer")
        if view.format != "B":
            view = view.cast("B")
        if writable and view.readonly:
            raise SerializationError(
                "update_data_tree writes in place and needs a writable buffer, such "
                "as a bytearray"
            )
        self.view = view
        try:
            self._parse()
        except (IndexError, ValueError, struct.error) as exc:
            raise MalformedBlobError(f"truncated or malformed frame: {exc}") from exc

    # -- the frame
    def _parse(self) -> None:
        view = self.view
        if len(view) < 5 or view[0] != MAGIC:
            raise MalformedBlobError("not a DataTreeBlob: the first byte is not 0xC9")
        if view[1] >> 4 != VERSION or view[1] & 0x0E:
            raise IncompatibleProfileError(
                f"DataTreeBlob version {view[1] >> 4} with flags {view[1] & 0x0F:#x}; "
                f"this build reads version {VERSION}"
            )
        self.pruned = bool(view[1] & FLAG_PRUNED)
        self.resolution, self.deepest, count = view[2], view[3], view[4]
        _require(
            self.deepest <= self.resolution <= MAX_RESOLUTION,
            f"deepest level {self.deepest} and resolution {self.resolution}",
        )
        cursor = 5
        raw = []
        for _ in range(count):
            length = view[cursor]
            name = bytes(view[cursor + 1 : cursor + 1 + length]).decode()
            cursor += 1 + length
            codes = list(view[cursor : cursor + 5])
            (leaf_max,) = struct.unpack_from("<Q", view, cursor + 5)
            (presence,) = struct.unpack_from("<H", view, cursor + 13)
            cursor += 15
            present = [level for level in range(16) if presence >> level & 1]
            widths = list(view[cursor : cursor + len(present)])
            cursor += len(present)
            raw.append((name, codes, leaf_max, present, widths))
        variables = []
        for name, codes, leaf_max, _, _ in raw:
            dtype, species, op, absence, weight = codes
            _require(
                dtype < len(_DTYPE_CODES)
                and species < len(SPECIES)
                and op < len(OPERATORS)
                and absence < len(ABSENCE)
                and (weight == 0xFF or weight < len(raw)),
                f"the codes of variable {name!r}",
            )
            variables.append(
                DataVariable(
                    name,
                    _DTYPE_CODES[dtype],
                    SPECIES[species],
                    OPERATORS[op],
                    ABSENCE[absence],
                    None if weight == 0xFF else raw[weight][0],
                    leaf_max or None,
                )
            )
        try:
            self.variables = _check_variables(variables)
        except SerializationError as exc:
            raise MalformedBlobError(f"invalid schema: {exc}") from exc
        self.widths: dict[str, list[str | None]] = {}
        for variable, (_, _, _, present, widths) in zip(self.variables, raw):
            row: list[str | None] = [None] * (self.deepest + 1)
            _require(
                not present or max(present) <= self.deepest,
                f"the levels of variable {variable.name!r}",
            )
            for level, code in zip(present, widths):
                _require(code < len(_DTYPE_CODES), f"a width of {variable.name!r}")
                row[level] = _DTYPE_CODES[code]
            _check_widths(variable, row)
            self.widths[variable.name] = row
        self.counts = list(struct.unpack_from(f"<{self.deepest + 1}I", view, cursor))
        cursor += 4 * (self.deepest + 1)
        _require(all(self.counts), "a level without nodes")
        (sections,) = struct.unpack_from("<H", view, cursor)
        presence = [
            (level, variable.absent == "unknown")
            for level in range(self.deepest + 1)
            for variable in self.variables
            if self.widths[variable.name][level]
        ]
        _require(
            sections == _section_count(self.deepest, presence),
            f"{sections} sections",
        )
        self.offsets = list(struct.unpack_from(f"<{sections}I", view, cursor + 2))
        self.frame_end = cursor + 2 + 4 * sections
        self._open()

    def _open(self) -> None:
        """Every section checked against the next offset; views made on demand."""
        layout: list[tuple[Any, str, int]] = [(("q", 0), "B", self.counts[0])]
        if self.deepest >= 1:
            layout.append((("p1", 0), "I", self.counts[1]))
        for level in range(self.deepest):
            layout.append((("first", level), "I", self.counts[level] + 1))
            if level >= 1:
                layout.append(
                    (("mask", level), _mask_format(level), self.counts[level])
                )
                layout.append((("cell", level), "B", (self.counts[level] + 7) // 8))
        for level in range(self.deepest + 1):
            for variable in self.variables:
                width = self.widths[variable.name][level]
                if width:
                    count = self.counts[level]
                    layout.append(
                        (("value", variable.name, level), DTYPES[width], count)
                    )
                    if variable.absent == "unknown":
                        bitmap = ("valid", variable.name, level)
                        layout.append((bitmap, "B", (count + 7) // 8))
        bound = -self.frame_end % 8 + self.frame_end
        specs: dict[Any, tuple[int, str, int]] = {}
        for offset, (key, fmt, count) in zip(self.offsets, layout):
            _require(offset == bound, f"a section at {offset}, expected {bound}")
            specs[key] = (offset, fmt, count)
            bound = offset + struct.calcsize(fmt) * count
            bound += -bound % 8
        _require(len(self.view) == bound, f"{len(self.view)} bytes, expected {bound}")
        self._specs = specs
        self._made: dict[Any, Any] = {}
        self.first = _Sections(self, "first")
        self.masks = _Sections(self, "mask")
        self.is_cell_bits = _Sections(self, "cell")
        self.values = _Sections(self, "value")
        self.valid = _Sections(self, "valid")

    def section(self, key: Any) -> Any:
        view = self._made.get(key)
        if view is None:
            view = self._made[key] = self._view(*self._specs[key])
        return view

    @property
    def quadrants(self) -> Any:
        return self.section(("q", 0))

    @property
    def level1(self) -> Any:
        return self.section(("p1", 0))

    def _view(self, offset: int, fmt: str, count: int) -> Any:
        if _LITTLE_ENDIAN_HOST:
            size = struct.calcsize(fmt)
            window: Any = self.view[offset : offset + size * count]
            return window.cast(fmt)
        return _Swapped(self.view, offset, fmt, count)

    # -- the tree
    def walk(self, level: int, chain: Sequence[int]) -> list[int]:
        """Slots of the cell's lexical ancestors the blob holds, level 0 first."""
        quadrants = self.quadrants
        for slot in range(len(quadrants)):
            if quadrants[slot] == chain[0]:
                break
        else:
            return []
        out = [slot]
        if level == 0 or self.deepest == 0:
            return out
        level1 = self.level1
        low, high = self.first[0][slot], self.first[0][slot + 1]
        slot = bisect_left(level1, chain[1], low, high)
        if slot == high or level1[slot] != chain[1]:
            return out
        out.append(slot)
        for finer in range(2, min(level, self.deepest) + 1):
            digit = chain[finer] - chain[finer - 1] * RATIOS[finer]
            mask = self.masks[finer - 1][slot]
            if not mask >> digit & 1:
                break
            slot = self.first[finer - 1][slot] + (mask & ((1 << digit) - 1)).bit_count()
            out.append(slot)
        return out

    def daughters(self, level: int, slot: int) -> tuple[int, int]:
        if level == self.deepest:
            return 0, 0
        return self.first[level][slot], self.first[level][slot + 1]

    def is_leaf(self, level: int, slot: int) -> bool:
        first, last = self.daughters(level, slot)
        return first == last

    def is_cell(self, level: int, slot: int) -> bool:
        bits = self.is_cell_bits.get(level)
        return bits is None or bool(bits[slot >> 3] >> (slot & 7) & 1)

    def value(self, name: str, level: int, slot: int) -> Value:
        valid = self.valid.get((name, level))
        if valid is not None and not valid[slot >> 3] >> (slot & 7) & 1:
            return None
        return self.values[(name, level)][slot]  # type: ignore[no-any-return]

    def node(self, level: int, slot: int) -> dict[str, Value]:
        return {
            variable.name: self.value(variable.name, level, slot)
            for variable in self.variables
            if self.widths[variable.name][level]
        }

    def write(self, name: str, level: int, slot: int, value: Value) -> None:
        width = self.widths[name][level]
        assert width is not None
        valid = self.valid.get((name, level))
        if valid is not None:
            bit = 1 << (slot & 7)
            if value is None:
                valid[slot >> 3] &= ~bit & 0xFF
            else:
                valid[slot >> 3] |= bit
        self.values[(name, level)][slot] = _zero(width) if value is None else value


def _require(condition: bool, what: str) -> None:
    if not condition:
        raise MalformedBlobError(f"malformed DataTreeBlob: {what}")


def _check_widths(variable: DataVariable, row: Sequence[str | None]) -> None:
    present = [width for width in row if width]
    if variable.species == "pyramid":
        _require(len(present) == len(row), f"pyramid {variable.name!r} misses a level")
    if variable.leaf_max is None:
        _require(
            all(width == variable.dtype for width in present),
            f"the widths of {variable.name!r}",
        )
    else:
        _require(
            all(width in _UNSIGNED for width in present),
            f"the widths of {variable.name!r}",
        )


@lru_cache(maxsize=16)
def _parsed(blob: bytes) -> _Blob:
    return _Blob(blob)


def _blob(buffer: Any, writable: bool = False) -> _Blob:
    """The parsed blob. A ``bytes`` blob cannot change, so its frame is parsed
    once and kept for the next reads (the last sixteen); a writable buffer is
    parsed on every call."""
    if not writable and type(buffer) is bytes:
        return _parsed(buffer)
    return _Blob(buffer, writable)


def _read_guard(exc: Exception) -> MalformedBlobError:
    return MalformedBlobError(f"malformed DataTreeBlob: {exc}")


# --------------------------------------------------------------------------
# Full decoding, with every invariant of the format checked
# --------------------------------------------------------------------------


def _levels(blob: _Blob) -> list[list[int]]:
    """Box positions of every node, rebuilt from the masks and checked."""
    deepest = blob.deepest
    quadrants = list(blob.quadrants)
    _require(
        all(a < b for a, b in zip(quadrants, quadrants[1:])) and quadrants[-1] < 4,
        "the quadrant codes",
    )
    levels = [quadrants]
    for level in range(deepest):
        first = list(blob.first[level])
        _require(
            first[0] == 0
            and first[-1] == blob.counts[level + 1]
            and all(a <= b for a, b in zip(first, first[1:])),
            f"the directory of level {level}",
        )
    if deepest >= 1:
        level1 = list(blob.level1)
        _require(
            all(a < b for a, b in zip(level1, level1[1:])) and level1[-1] < _LEVEL1_BOX,
            "the resolution-1 positions",
        )
        box = BOX_COLUMNS * BOX_ROWS
        first = blob.first[0]
        for slot, quadrant in enumerate(quadrants):
            _require(
                all(
                    level1[i] // box == quadrant
                    for i in range(first[slot], first[slot + 1])
                ),
                "a resolution-1 node under the wrong quadrant",
            )
        levels.append(level1)
    for level in range(1, deepest):
        ratio = RATIOS[level + 1]
        below: list[int] = []
        first = blob.first[level]
        for slot, (position, mask) in enumerate(zip(levels[level], blob.masks[level])):
            _require(
                mask < 1 << ratio and mask.bit_count() == first[slot + 1] - first[slot],
                f"the mask of node {slot} at level {level}",
            )
            while mask:
                low = mask & -mask
                below.append(position * ratio + low.bit_length() - 1)
                mask ^= low
        levels.append(below)
    return levels


def _decode(blob: _Blob) -> dict[str, dict[str, Value]]:
    """Leaves first, so a leaf's own fault is reported before its parent's."""
    levels = _levels(blob)
    rows: list[list[tuple[str, dict[str, Value]]]] = []
    for level in range(blob.deepest, -1, -1):
        row = []
        for slot, position in enumerate(levels[level]):
            leaf = blob.is_leaf(level, slot)
            cell = _exists(level, position)
            _require(
                blob.is_cell(level, slot) == cell and (cell or not leaf),
                f"node {slot} at level {level}: a cell flag, or a leaf naming no cell",
            )
            _check_node(blob, level, slot, position, leaf)
            if cell:
                row.append(
                    (_spelling(OFFSETS[level] + position), blob.node(level, slot))
                )
        rows.append(row)
    return {spelling: values for row in reversed(rows) for spelling, values in row}


def _check_node(blob: _Blob, level: int, slot: int, position: int, leaf: bool) -> None:
    """Unknown slots hold zero, leaves hold sums that fit, internals agree."""
    for variable in blob.variables:
        width = blob.widths[variable.name][level]
        if width is None:
            continue
        value = blob.value(variable.name, level, slot)
        if value is None:
            _require(
                blob.values[(variable.name, level)][slot] == 0,
                f"an unknown slot of {variable.name!r} holding a value",
            )
        if variable.species != "pyramid":
            continue
        if leaf:
            try:
                _check_leaf(variable, value, level, position, blob.resolution)
            except SerializationError as exc:
                raise MalformedBlobError(f"malformed DataTreeBlob: {exc}") from exc
            continue
        first, last = blob.daughters(level, slot)
        expected = _aggregate(
            variable,
            [blob.value(variable.name, level + 1, d) for d in range(first, last)],
            [
                blob.value(variable.weight, level + 1, d) if variable.weight else 1
                for d in range(first, last)
            ],
        )
        _require(
            _same(_as_stored(width, expected), value),
            f"{variable.name!r} at node {slot} of level {level} is not the "
            f"{variable.op} of its daughters",
        )


def _as_stored(width: str, value: Value) -> Value:
    """A value as the bytes hold it, for comparing with what was read."""
    if value is None or width != "f32":
        return value
    return float(struct.unpack("<f", struct.pack("<f", value))[0])


# --------------------------------------------------------------------------
# The public functions
# --------------------------------------------------------------------------

Buffer = Union[bytes, bytearray, memoryview]
Pairs = Union[Mapping[Any, Mapping[str, Any]], Iterable[tuple[Any, Mapping[str, Any]]]]


def encode_data_tree(
    data: Pairs,
    variables: Iterable[DataVariable],
    *,
    observed: Pairs | None = None,
    resolution: int | None = None,
    prune: bool = False,
) -> bytes:
    """Encode values on cells as a DataTreeBlob.

    The input is pairs, never an index: an index folds a daughter into its
    parent, and the daughter's own value would go with it.

    Args:
        data: The pyramid values, as a mapping of cells to ``{name: value}``
            or as pairs of them. A cell is a spelling or an integer of the
            integer index. Every cell given is a leaf of the tree: the values
            of its ancestors are aggregated, never given. ``None`` is the
            variable's absence.
        variables: The columns, as :class:`DataVariable`.
        observed: The level values, in the same form, at any node.
        resolution: The observation resolution: the unit of every count and
            the floor of every read. Defaults to the finest cell in
            ``data``. A leaf above it stands for the block of cells beneath
            it.
        prune: Fold every complete family of leaves that agree on every
            variable into its parent, as a block, bottom-up. Off by default.

    Returns:
        The blob. The same input gives the same bytes, in whatever order the
        pairs arrive.

    Raises:
        InvalidIndexError: If a cell is malformed.
        NonExistentCellError: If a cell names no cell.
        NonAtomicIndexError: If an entry holds more than one cell.
        ResolutionError: If ``resolution`` is outside 0 to 13, or a cell is
            finer than it.
        SerializationError: If the schema or a value is invalid: a name
            unknown or given twice, a value outside its type, a pyramid value
            at a node that has daughters, or a block that absorbs the border
            or whose sum does not divide evenly over its cells.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = encode_data_tree(
        ...     {"NE(0625/0451(1))": {"n": 1}, "NE(0625/0451(2))": {"n": 1}}, [n]
        ... )
        >>> read_data_cell(blob, "NE(0625/0451)")
        {'n': 2}
    """
    declared = _check_variables(variables)
    given, spellings = _ingest(declared, data, observed)
    tree = _build(declared, given, spellings, resolution)
    if prune:
        levels = {v.name for v in declared if v.species == "level"}
        rebuilt = {key: dict(entry) for key, entry in _prune(tree).items()}
        for key, entry in given.items():
            observations = {name: v for name, v in entry.items() if name in levels}
            if observations:
                rebuilt.setdefault(key, {}).update(observations)
        names = {key: _spelling(OFFSETS[key[0]] + key[1]) for key in rebuilt}
        tree = _build(declared, rebuilt, names, tree.resolution, pruned=True)
    return _serialize(tree)


def decode_data_tree(blob: Buffer) -> dict[str, dict[str, Value]]:
    """Every cell of a DataTreeBlob with its values, every invariant checked.

    Internal nodes are there with their own values, so a parent and its
    daughter both come back. A prefix that names no cell -- the stem a
    border cell's surplus children hang from -- is a node of the blob and is
    not returned.

    Args:
        blob: The bytes.

    Returns:
        A mapping of spellings to ``{name: value}``, coarsest level first and
        in canonical order within a level. Each node carries the variables
        that have values at its resolution; an unknown value is ``None``.

    Raises:
        MalformedBlobError: If the bytes break the format anywhere, including
            an aggregate that is not its daughters' or a leaf naming no cell.
        IncompatibleProfileError: If the blob is of a version this build
            does not read.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = encode_data_tree({"NE(0625/0451(3))": {"n": 1}}, [n])
        >>> decode_data_tree(blob)
        {'NE': {'n': 1}, 'NE(0625/0451)': {'n': 1}, 'NE(0625/0451(3))': {'n': 1}}
    """
    return _decode(_blob(blob))


def data_tree_schema(blob: Buffer) -> dict[str, Any]:
    """What a DataTreeBlob declares, read from its frame.

    Args:
        blob: The bytes.

    Returns:
        ``resolution`` (the observation resolution), ``deepest`` (the
        finest level holding nodes), ``pruned`` (whether the writer folded
        families into blocks), ``variables`` (a tuple of
        :class:`DataVariable`) and ``nodes`` (the node count per level).

    Raises:
        MalformedBlobError: If the frame is malformed.
        IncompatibleProfileError: If the version is not one this build reads.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = encode_data_tree({"NE(0625/0451(3))": {"n": 1}}, [n])
        >>> data_tree_schema(blob)["nodes"]
        (1, 1, 1)
    """
    reader = _blob(blob)
    return {
        "resolution": reader.resolution,
        "deepest": reader.deepest,
        "pruned": reader.pruned,
        "variables": reader.variables,
        "nodes": tuple(reader.counts),
    }


def read_data_path(
    blob: Buffer, cell: str | int | Iterable[int]
) -> list[tuple[str, dict[str, Value]]] | list[list[tuple[str, dict[str, Value]]]]:
    """A cell's ancestors that the blob holds, with their values, read in place.

    One walk down the tree: a binary search at resolution 1, then one step
    of constant time per level, each landing on an ancestor. The path is
    lexical; a prefix that names no cell is walked through and left out.

    Args:
        blob: The bytes.
        cell: A spelling or an integer of the integer index. A compositional
            index, or an iterable of integers, reads each cell.

    Returns:
        ``(spelling, {name: value})`` from the quadrant down to the deepest
        node the blob holds on the cell's path; the cell itself comes last
        only when the blob holds it. A list of such paths, aligned with the
        cells, when more than one was given.

    Raises:
        InvalidIndexError: If a cell is malformed.
        NonExistentCellError: If a cell names no cell.
        MalformedBlobError: If the frame is malformed.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = encode_data_tree({"NE(0625/0451(3))": {"n": 1}}, [n])
        >>> [spelling for spelling, _ in read_data_path(blob, "NE(0625/0451(3))")]
        ['NE', 'NE(0625/0451)', 'NE(0625/0451(3))']
    """
    reader = _blob(blob)
    keys, single = _cells(cell)
    try:
        paths = [_path(reader, level, position) for level, position, _ in keys]
    except (IndexError, struct.error) as exc:
        raise _read_guard(exc) from exc
    return paths[0] if single else paths


def read_data_cell(
    blob: Buffer, cell: str | int | Iterable[int]
) -> dict[str, Value] | list[dict[str, Value]]:
    """The values of a cell, read in place.

    A cell the blob holds answers with its own values. A cell inside a block
    -- beneath a leaf above the observation resolution -- answers with the
    block's: the same value, or for a sum its even share. Any other cell
    answers with each variable's absence.

    Args:
        blob: The bytes.
        cell: A spelling or an integer of the integer index. A compositional
            index, or an iterable of integers, reads each cell.

    Returns:
        ``{name: value}`` with every pyramid variable and the level variables
        that have values at the cell's resolution. A list of them, aligned
        with the cells, when more than one was given.

    Raises:
        InvalidIndexError: If a cell is malformed.
        NonExistentCellError: If a cell names no cell.
        ResolutionError: If a cell is finer than the observation resolution.
        MalformedBlobError: If the frame is malformed.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", absent="zero", leaf_max=1)
        >>> blob = encode_data_tree({"NE(0625/0451(3))": {"n": 1}}, [n])
        >>> read_data_cell(blob, "NE(0625/0451(4))")
        {'n': 0}
    """
    reader = _blob(blob)
    keys, single = _cells(cell)
    try:
        answers = [_cell_values(reader, *key) for key in keys]
    except (IndexError, struct.error) as exc:
        raise _read_guard(exc) from exc
    return answers[0] if single else answers


def slice_data_tree(blob: Buffer, cell: str | int) -> bytes:
    """The sub-tree of a node, with the path that reaches it, as a new blob.

    The descendants of a node are one contiguous run of every level, so the
    slice copies runs. The ancestors on the path keep their level values,
    and their pyramid values become the aggregate of the slice alone.

    Args:
        blob: The bytes.
        cell: A node the blob holds, as a spelling or an integer.

    Returns:
        A DataTreeBlob with the same variables and observation resolution.

    Raises:
        InvalidIndexError: If the cell is malformed.
        NonExistentCellError: If it names no cell.
        NonAtomicIndexError: If it holds more than one cell.
        SerializationError: If the blob holds no node for it.
        MalformedBlobError: If the frame is malformed.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = encode_data_tree(
        ...     {"NE(0625/0451(1))": {"n": 1}, "NE(0625/0451(2))": {"n": 1}}, [n]
        ... )
        >>> list(decode_data_tree(slice_data_tree(blob, "NE(0625/0451(2))")))
        ['NE', 'NE(0625/0451)', 'NE(0625/0451(2))']
    """
    reader = _blob(blob)
    level, position, spelling = _one_cell(cell)
    chain = _chain(level, position)
    try:
        slots = reader.walk(level, chain)
        if len(slots) != level + 1:
            raise SerializationError(f"the blob holds no node for {spelling}")
        tree = _sliced(reader, level, chain, slots)
    except (IndexError, struct.error) as exc:
        raise _read_guard(exc) from exc
    return _serialize(tree)


def update_data_tree(
    buffer: bytearray | memoryview, cell: str | int, values: Mapping[str, Any]
) -> None:
    """Write values at one node in place, and refresh the aggregates above it.

    Values have fixed widths and fixed places, so an update rewrites a slot
    per level of the path and nothing else. A pyramid value is written at a
    leaf; a level value at any node that has values at its resolution.
    Adding or removing a node changes the structure, and is a new encoding.

    Args:
        buffer: The blob, writable: a ``bytearray`` or a writable
            ``memoryview``.
        cell: A node the blob holds, as a spelling or an integer.
        values: ``{name: value}``; ``None`` is the variable's absence.

    Raises:
        InvalidIndexError: If the cell is malformed.
        NonExistentCellError: If it names no cell.
        NonAtomicIndexError: If it holds more than one cell.
        SerializationError: If the buffer is read-only, the blob holds no node
            for the cell, a name is unknown or has no values at the cell's
            resolution, a pyramid value targets a node with daughters, or a
            value is invalid. Nothing is written when anything is refused.
        MalformedBlobError: If the frame is malformed.

    Example:
        >>> n = DataVariable("n", "u32", op="sum", leaf_max=1)
        >>> blob = bytearray(encode_data_tree({"NE(0625/0451(3))": {"n": 0}}, [n]))
        >>> update_data_tree(blob, "NE(0625/0451(3))", {"n": 1})
        >>> read_data_cell(blob, "NE(0625/0451)")
        {'n': 1}
    """
    reader = _blob(buffer, writable=True)
    level, position, spelling = _one_cell(cell)
    if not isinstance(values, Mapping):
        raise SerializationError(
            f"the values are a mapping of variable names, not {type(values).__name__}"
        )
    chain = _chain(level, position)
    try:
        slots = reader.walk(level, chain)
        if len(slots) != level + 1:
            raise SerializationError(f"the blob holds no node for {spelling}")
        pending = _pending(reader, level, chain, slots, spelling, values)
    except (IndexError, struct.error) as exc:
        raise _read_guard(exc) from exc
    for (name, depth, slot), value in pending.items():
        reader.write(name, depth, slot, value)


# --------------------------------------------------------------------------
# What the public functions are made of
# --------------------------------------------------------------------------


def _path(
    reader: _Blob, level: int, position: int
) -> list[tuple[str, dict[str, Value]]]:
    chain = _chain(level, position)
    return [
        (_spelling(OFFSETS[depth] + chain[depth]), reader.node(depth, slot))
        for depth, slot in enumerate(reader.walk(level, chain))
        if reader.is_cell(depth, slot)
    ]


def _cell_values(
    reader: _Blob, level: int, position: int, spelling: str
) -> dict[str, Value]:
    if level > reader.resolution:
        raise ResolutionError(
            f"{spelling} is finer than the observation resolution {reader.resolution}"
        )
    slots = reader.walk(level, _chain(level, position))
    if len(slots) == level + 1:
        return reader.node(level, slots[-1])
    reached = len(slots) - 1
    block = reached >= 0 and reader.is_leaf(reached, slots[-1])
    out: dict[str, Value] = {}
    for variable in reader.variables:
        if variable.species == "pyramid":
            if block:
                stored = reader.value(variable.name, reached, slots[-1])
                out[variable.name] = _inherited(variable, stored, reached, level)
            else:
                out[variable.name] = _absent(variable)
        elif level <= reader.deepest and reader.widths[variable.name][level]:
            out[variable.name] = _absent(variable)
    return out


def _inherited(variable: DataVariable, value: Value, block: int, level: int) -> Value:
    """A block's value as one cell of it reads it: a sum is shared evenly."""
    if value is None or variable.op != "sum":
        return value
    cells = _cells_between(block, level)
    if isinstance(value, int):
        return value // cells
    return _rounded(variable.dtype, _exact(value) / cells)


def _sliced(reader: _Blob, level: int, chain: Sequence[int], slots: list[int]) -> _Tree:
    """The tree of a slice: the path, then every run beneath the node."""
    runs = [(slot, slot + 1) for slot in slots]
    levels = [[chain[depth]] for depth in range(level + 1)]
    first, last = runs[-1]
    for depth in range(level, reader.deepest):
        low, high = reader.first[depth][first], reader.first[depth][last]
        if low == high:
            break
        ratio = RATIOS[depth + 1]
        below: list[int] = []
        if depth == 0:
            below = [reader.level1[slot] for slot in range(low, high)]
        else:
            for slot, position in zip(range(first, last), levels[depth]):
                mask = reader.masks[depth][slot]
                while mask:
                    bit = mask & -mask
                    below.append(position * ratio + bit.bit_length() - 1)
                    mask ^= bit
        levels.append(below)
        runs.append((low, high))
        first, last = low, high
    values: dict[str, list[list[Value] | None]] = {}
    for variable in reader.variables:
        values[variable.name] = [
            (
                [reader.value(variable.name, depth, slot) for slot in range(low, high)]
                if reader.widths[variable.name][depth]
                else None
            )
            for depth, (low, high) in enumerate(runs)
        ]
    tree = _Tree(reader.variables, reader.resolution, levels, values, reader.pruned)
    for depth in range(level - 1, -1, -1):
        for variable in reader.variables:
            if variable.species == "pyramid":
                tree.row(variable.name, depth)[0] = _aggregate(
                    variable,
                    tree.row(variable.name, depth + 1)[0:1],
                    _weights(tree, variable, depth + 1, 0, 1),
                )
    return tree


def _pending(
    reader: _Blob,
    level: int,
    chain: Sequence[int],
    slots: list[int],
    spelling: str,
    values: Mapping[str, Any],
) -> dict[tuple[str, int, int], Value]:
    """Every slot an update writes, computed before anything is written."""
    slot = slots[-1]
    leaf = reader.is_leaf(level, slot)
    by_name = {variable.name: variable for variable in reader.variables}
    pending: dict[tuple[str, int, int], Value] = {}
    changed: set[str] = set()
    for name, raw in values.items():
        variable = by_name.get(name)
        if variable is None:
            raise SerializationError(f"there is no variable {name!r}")
        if reader.widths[name][level] is None:
            raise SerializationError(
                f"{name!r} has no values at resolution {level} in this blob"
            )
        if variable.species == "pyramid" and not leaf:
            raise SerializationError(
                f"{name!r} at {spelling} is the aggregate of its daughters; write "
                f"at a leaf"
            )
        value = _leaf_value(variable, raw, spelling)
        if variable.species == "pyramid":
            _check_leaf(variable, value, level, chain[level], reader.resolution)
            changed.add(name)
        pending[(name, level, slot)] = value

    def current(name: str, depth: int, at: int) -> Value:
        key = (name, depth, at)
        return pending[key] if key in pending else reader.value(name, depth, at)

    affected = [
        variable
        for variable in reader.variables
        if variable.species == "pyramid"
        and (variable.name in changed or variable.weight in changed)
    ]
    for depth in range(level - 1, -1, -1):
        first, last = reader.daughters(depth, slots[depth])
        for variable in affected:
            weight = variable.weight
            aggregate = _aggregate(
                variable,
                [current(variable.name, depth + 1, d) for d in range(first, last)],
                [
                    current(weight, depth + 1, d) if weight else 1
                    for d in range(first, last)
                ],
            )
            _check_sum(variable, aggregate, depth, chain[depth])
            pending[(variable.name, depth, slots[depth])] = aggregate
    for name, depth, at in pending:
        reader.values[(name, depth)][at]  # every slot is there before any is written
    return pending
