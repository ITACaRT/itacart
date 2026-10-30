"""TreeBlob: dense binary encoding of a compositional index as a cell set.

Canonical prefix-index form, content-addressable. A resolution-13 leaf
occupies exactly 10 bytes.

Formal properties: determinism, ``decode(encode(x)) == recompose(x)``,
``encode(x) == encode(recompose(x))``, idempotent recomposition,
bit-exact round trip, content-addressability, 10 bytes at resolution 13,
and monotone :func:`prefix_at_resolution_binary`.

The encoding depends only on the index being prefix-structured with a
known per-level alphabet, so it ports to other DGGS with little change.

Wire format
-----------

A *tree blob* encodes a whole compositional tree as one :class:`bytes`::

    [magic = 0xC7]                        1 byte
    [version (4 bits) | flags (4 bits)]   1 byte
    [quadrant (2 bits) | reserved (6)]    1 byte
    [resolution-1 child count, 16 bits]   2 bytes
    [bit-packed children, pre-order]      variable

Packing is big-endian, MSB-first within every byte. Each node carries a
child-count field whose width follows the resolution of its children; a
count of zero marks a leaf.

A *node blob* encodes one root-to-cell path and is the form returned by
:func:`iter_leaves`::

    [0xA0 | resolution]                   1 byte
    [quadrant + column/row + refinements] variable

Index space
-----------

Columns and rows are written as they appear in the index, without
offset. Column ``0`` is the meridian triangle and exists only in the
eastern quadrants; column ``2004`` is a prefix that spells the fifth
child of a trapezoidal cell in the last column, and is never a
resolution-1 cell of its own. The refinement ratio of a level is not the
number of children a node has: a trapezoidal parent yields one more than
the ratio, at every depth, so the count field is bounded by its own
width and not by the ratio.

Row ``1000`` holds the two polar caps, ``NE(0000/1000)`` and
``SE(0000/1000)``. The grid has them; this version of the binary codecs
does not serialize them. That is a restriction of the codecs and not a
property of the grid, so a row-1000 index is refused with a message that
says so, while a row past it is refused as outside the grid, which it is.
The GeometryBlob reads its vertices through the same parser and refuses
the caps with the same message.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator, Sequence, Union

from ..constants import QUADRANTS, refinement_alphabet
from ..exceptions import InvalidIndexError, MalformedBlobError, ResolutionError

__all__ = [
    "encode_tree",
    "decode_tree",
    "encode_node",
    "decode_node",
    "serialize_to_blob",
    "deserialize_from_blob",
    "is_ancestor_binary",
    "prefix_at_resolution_binary",
    "resolution_of_binary",
    "validate_tree",
    "iter_leaves",
    "count_vertices",
    "recompose_to_prefix_form",
]

# ---------------------------------------------------------------------------
# Format constants
# ---------------------------------------------------------------------------

MAGIC_TREE = 0xC7
MAGIC_NODE_HIGH_NIBBLE = 0xA
FORMAT_VERSION = 1

MIN_RESOLUTION = 1
MAX_RESOLUTION = 13

#: Greatest column that is a resolution-1 cell in its own right.
MAX_CELL_COLUMN = 2003
#: Greatest column the grammar admits. It spells the fifth child of a
#: trapezoidal cell in column 2003 and only ever appears refined.
MAX_PREFIX_COLUMN = 2004
MIN_COLUMN = 0
MIN_ROW = 0
#: Greatest row this version of the codecs writes. It is not the last row
#: of the grid: see :data:`POLAR_CAP_ROW`.
MAX_ROW = 999
#: The row of the two polar caps, one past :data:`MAX_ROW`. The grid has
#: it and this version of the codecs does not serialize it.
POLAR_CAP_ROW = 1000
_POLAR_CAP_REFUSAL = (
    "row 1000 holds the polar caps, which this version of the binary codecs "
    "does not serialize"
)

X_BITS = 11  # 2**11 = 2048 > 2004
Y_BITS = 10  # 2**10 = 1024 > 999

COUNT_WIDTH_EVEN_CHILDREN = 3
COUNT_WIDTH_ODD_CHILDREN = 5
COUNT_WIDTH_ROOT_CHILDREN = 16
#: Magic, version and the quadrant-group count. Group headers follow, packed.
HEADER_SIZE_BYTES = 3
#: Resolution of the synthetic root that holds one node per quadrant.
GLOBE_RESOLUTION = -1
ROOT_MAX_CHILDREN = (1 << COUNT_WIDTH_ROOT_CHILDREN) - 1

_QUADRANT_TO_BITS = {name: bits for bits, name in enumerate(QUADRANTS)}
_BITS_TO_QUADRANT = {bits: name for name, bits in _QUADRANT_TO_BITS.items()}

_QUADRANT_RE = re.compile(r"(" + "|".join(QUADRANTS) + r")\(")
_RES1_RE = re.compile(r"(\d{1,4})/(\d{1,4})")

Component = Union[str, int, "tuple[int, int]"]


def _pair(component: Component) -> "tuple[int, int]":
    """Narrow a resolution-1 component to its column and row."""
    if not isinstance(component, tuple):  # pragma: no cover - internal invariant
        raise InvalidIndexError(f"resolution-1 component must be a pair: {component!r}")
    return component


# ---------------------------------------------------------------------------
# Level alphabets, derived rather than transcribed
# ---------------------------------------------------------------------------


def _alphabet(resolution: int) -> tuple[str, ...]:
    """Refinement tokens addressing a node at ``resolution``."""
    return refinement_alphabet(resolution)


def _component_width(resolution: int) -> int:
    """Bit width of a refinement field. Resolution 1 uses ``X_BITS + Y_BITS``."""
    size = len(_alphabet(resolution))
    width = 1
    while (1 << width) < size:
        width += 1
    return width


def _count_width_for_children_at(child_res: int) -> int:
    """Width of the count field for children living at ``child_res``."""
    if child_res > MAX_RESOLUTION:
        return COUNT_WIDTH_EVEN_CHILDREN
    if len(_alphabet(child_res)) <= 4:
        return COUNT_WIDTH_EVEN_CHILDREN
    return COUNT_WIDTH_ODD_CHILDREN


def _count_width_for_a_leaf_at(resolution: int) -> int:
    """Width of the zero-valued count field a leaf must still emit."""
    return _count_width_for_children_at(resolution + 1)


def _max_children_at(child_res: int) -> int:
    """Greatest child count the format can express at ``child_res``.

    Bounded by the width of the count field, not by
    :func:`~itacart.resolutions.refinement_ratio`. The ratio describes a
    node in the grid interior; a trapezoidal parent yields one child
    more than the ratio, and that surplus does not attenuate with depth.
    Encoding the ratio as a hard limit here would refuse the whole
    eastern border of the grid.
    """
    return (1 << _count_width_for_children_at(child_res)) - 1


#: Quadrants holding the meridian triangle. It is one cell, not two, so
#: it is addressed from the east and column zero is empty to the west.
EASTERN_QUADRANTS = ("NE", "SE")


def _row_refusal(row: int) -> "str | None":
    """Why ``row`` cannot be written, or ``None`` when it can.

    Two reasons, and they are different statements. Row 1000 is part of
    the grid -- it holds the polar caps -- and the refusal is this
    version's restriction. Any other row past :data:`MAX_ROW` is outside
    the grid. The parser and both decoders ask here, so the three places
    that read a row give the same answer.
    """
    if MIN_ROW <= row <= MAX_ROW:
        return None
    if row == POLAR_CAP_ROW:
        return _POLAR_CAP_REFUSAL
    return f"row {row} is outside the grid"


def _column_is_addressable(column: int, refined: bool, quadrant: str) -> bool:
    """Whether ``column`` may stand in ``quadrant``.

    Two structural exceptions, both constant and neither dependent on
    latitude. Column 0 is the meridian triangle, a single cell addressed
    from the east, so it is absent west of the meridian at every
    resolution. Column 2004 is the trapezoid exception: the grammar
    admits it so the fifth child of a trapezoidal cell in column 2003
    can be spelled, and it is a prefix only, never a cell.

    Whether a column exists at a *given row* is a separate question,
    since the last column retreats with the cosine of the latitude.
    That is what :func:`itacart.is_valid_cell` answers; this codec
    validates the wire format and does not re-derive the boundary.
    """
    if column == MIN_COLUMN:
        return quadrant in EASTERN_QUADRANTS
    if column <= MAX_CELL_COLUMN:
        return True
    return column == MAX_PREFIX_COLUMN and refined


# ---------------------------------------------------------------------------
# Bit-level primitives
# ---------------------------------------------------------------------------


class _BitWriter:
    """MSB-first bit accumulator."""

    def __init__(self) -> None:
        self._bytes = bytearray()
        self._current = 0
        self._filled = 0

    def write(self, value: int, width: int) -> None:
        if value < 0 or value >= (1 << width):  # pragma: no cover - internal guard
            raise MalformedBlobError(f"value {value} does not fit in {width} bits")
        for shift in range(width - 1, -1, -1):
            self._current = (self._current << 1) | ((value >> shift) & 1)
            self._filled += 1
            if self._filled == 8:
                self._bytes.append(self._current)
                self._current = 0
                self._filled = 0

    def write_byte_aligned(self, value: int) -> None:
        self._bytes.append(value)

    def to_bytes(self) -> bytes:
        if self._filled:
            return bytes(self._bytes) + bytes([self._current << (8 - self._filled)])
        return bytes(self._bytes)


class _BitReader:
    """MSB-first bit reader over a blob."""

    def __init__(self, data: bytes, start_bit: int = 0) -> None:
        self._data = data
        self._pos = start_bit

    def read(self, width: int) -> int:
        if self._pos + width > len(self._data) * 8:
            raise MalformedBlobError("blob ends inside a field")
        value = 0
        for _ in range(width):
            byte = self._data[self._pos >> 3]
            value = (value << 1) | ((byte >> (7 - (self._pos & 7))) & 1)
            self._pos += 1
        return value

    def remaining_bits(self) -> int:
        return len(self._data) * 8 - self._pos


# ---------------------------------------------------------------------------
# Internal tree
# ---------------------------------------------------------------------------


@dataclass
class _Node:
    """One node of a compositional tree."""

    resolution: int
    component: Component
    children: "list[_Node]" = field(default_factory=list)


class _Parser:
    """Recursive-descent parser for the compositional index grammar."""

    def __init__(self, source: str) -> None:
        self.source = re.sub(r"\s+", "", source)
        self.pos = 0
        self.quadrant = ""

    def _fail(self, message: str) -> "InvalidIndexError":
        snippet = self.source[self.pos : self.pos + 16]
        return InvalidIndexError(f"{message} at position {self.pos}: {snippet!r}")

    def _consume(self, text: str) -> None:
        if not self.source.startswith(text, self.pos):
            raise self._fail(f"expected {text!r}")
        self.pos += len(text)

    def parse(self) -> _Node:
        """Parse the whole index into a globe root holding quadrant groups.

        An index may name more than one quadrant, comma-separated at the
        top level, and :func:`itacart.compose` produces exactly that for
        any cell set that crosses the equator or the meridian. A region
        straddling the equator is ordinary, not exotic, so the grammar
        the codec accepts has to be the grammar the package writes.
        """
        globe = _Node(resolution=GLOBE_RESOLUTION, component="")
        while True:
            globe.children.append(self._parse_quadrant())
            if self.pos < len(self.source) and self.source[self.pos] == ",":
                self.pos += 1
                continue
            break
        if self.pos != len(self.source):
            raise self._fail("trailing characters after the closing parenthesis")
        return globe

    def _parse_quadrant(self) -> _Node:
        match = _QUADRANT_RE.match(self.source, self.pos)
        if not match:
            raise self._fail("expected a quadrant code followed by '('")
        quadrant = match.group(1)
        self.quadrant = quadrant
        self.pos += len(quadrant) + 1
        root = _Node(resolution=0, component=quadrant)
        self._parse_children(root, child_res=1)
        self._consume(")")
        return root

    def _parse_children(self, parent: _Node, child_res: int) -> None:
        parent.children.append(self._parse_node(child_res))
        while self.pos < len(self.source) and self.source[self.pos] == ",":
            self.pos += 1
            parent.children.append(self._parse_node(child_res))

    def _parse_node(self, resolution: int) -> _Node:
        node = self._parse_res1() if resolution == 1 else self._parse_code(resolution)
        if self.pos < len(self.source) and self.source[self.pos] == "(":
            if resolution >= MAX_RESOLUTION:
                raise self._fail(f"cannot descend below resolution {MAX_RESOLUTION}")
            self.pos += 1
            self._parse_children(node, child_res=resolution + 1)
            self._consume(")")
        if resolution == 1:
            column = _pair(node.component)[0]
            if not _column_is_addressable(column, bool(node.children), self.quadrant):
                raise self._fail(
                    f"column {column} is not addressable in {self.quadrant}"
                )
        return node

    def _parse_res1(self) -> _Node:
        match = _RES1_RE.match(self.source, self.pos)
        if not match:
            raise self._fail("expected a resolution-1 pair XXXX/YYYY")
        column = int(match.group(1))
        row = int(match.group(2))
        if not MIN_COLUMN <= column <= MAX_PREFIX_COLUMN:
            raise self._fail(f"column {column} is outside the grid")
        refusal = _row_refusal(row)
        if refusal is not None:
            raise self._fail(refusal)
        self.pos += len(match.group(0))
        return _Node(resolution=1, component=(column, row))

    def _parse_code(self, resolution: int) -> _Node:
        alphabet = _alphabet(resolution)
        width = max(len(token) for token in alphabet)
        token = self.source[self.pos : self.pos + width]
        if token not in alphabet:
            raise self._fail(f"expected a refinement code of resolution {resolution}")
        self.pos += len(token)
        return _Node(resolution=resolution, component=token)


def _parse(index: str) -> _Node:
    """Parse an index string into the internal tree."""
    if not isinstance(index, str):
        raise InvalidIndexError(f"index must be a string, got {type(index).__name__}")
    if not index:
        raise InvalidIndexError("empty index")
    return _Parser(index).parse()


def _render(globe: _Node) -> str:
    """Render the internal tree back to an index string."""
    parts: list[str] = []

    def emit(node: _Node) -> None:
        if node.resolution == 1:
            column, row = _pair(node.component)
            parts.append(f"{column:04d}/{row:04d}")
        else:
            parts.append(str(node.component))
        if node.children:
            parts.append("(")
            for position, child in enumerate(node.children):
                if position:
                    parts.append(",")
                emit(child)
            parts.append(")")

    for position, root in enumerate(globe.children):
        if position:
            parts.append(",")
        parts.append(str(root.component))
        parts.append("(")
        for offset, child in enumerate(root.children):
            if offset:
                parts.append(",")
            emit(child)
        parts.append(")")
    return "".join(parts)


def _sort_key(node: _Node) -> "tuple[int, ...]":
    if node.resolution == GLOBE_RESOLUTION + 1:
        return (_QUADRANT_TO_BITS[str(node.component)],)
    if node.resolution == 1:
        return _pair(node.component)
    return (_alphabet(node.resolution).index(str(node.component)),)


def _merge(node: _Node, absorb: bool = False) -> _Node:
    """Merge duplicate siblings and order them, depth first.

    With ``absorb``, a sibling written whole -- with no children -- takes
    in every other occurrence of the same component, so a cell absorbs any
    descendant also written under it. That is the rule
    :func:`itacart.compose`, :func:`itacart.normalize` and
    :func:`itacart.compact_cells` already follow, and it is what the
    tree-blob entry points ask for. Without it, merging the two
    occurrences in ``NE(0500/0300,0500/0300(1))`` hands the whole-cell
    occurrence the other one's children and the parent is lost: the
    output covers a quarter of the input.

    The default stays as it was because the GeometryBlob parses its
    vertices through here, and that format's behaviour is not this
    function's to change.
    """
    merged: "dict[str, _Node]" = {}
    order: list[str] = []
    whole: "set[str]" = set()
    for child in node.children:
        key = f"{child.resolution}:{child.component!r}"
        if key not in merged:
            merged[key] = _Node(child.resolution, child.component)
            order.append(key)
        if not child.children:
            whole.add(key)
        merged[key].children.extend(child.children)
    if absorb:
        for key in whole:
            merged[key].children = []
    node.children = [_merge(merged[key], absorb) for key in order]
    node.children.sort(key=_sort_key)
    return node


def recompose_to_prefix_form(index: str) -> str:
    """Rewrite an index in canonical prefix form.

    The textual counterpart of what :func:`encode_tree` does internally:
    every distinct prefix appears once, siblings are ordered, and the
    leaf set is preserved exactly. Idempotent.

    One rewrite touches the leaf set, and it preserves the region. A cell
    written whole absorbs any descendant also written under it, as
    :func:`itacart.compose` and :func:`itacart.normalize` do, so
    ``NE(0500/0300,0500/0300(1))`` recomposes to ``NE(0500/0300)``. The
    absorption is lexical, like theirs: a descendant spelled under another
    column's prefix -- the children a border cell fathers to the east --
    is not under the cell's own spelling and stays as written, which loses
    nothing because both leaves are kept.

    This is not :func:`itacart.normalize`, which collapses a complete
    sibling set into its parent. Compaction preserves coverage while
    changing the leaf set, and a blob's identity is its leaf set, so the
    two canonical forms answer different questions and have different
    fixed points.

    Args:
        index: Compositional index string.

    Returns:
        The index in canonical prefix form.

    Raises:
        InvalidIndexError: If the index is malformed, or if it names row
            1000, which holds the polar caps: this version of the binary
            codecs does not serialize them.

    Example:
        >>> recompose_to_prefix_form("NE(0625/0451(2,1))")
        'NE(0625/0451(1,2))'
    """
    return _render(_merge(_parse(index), absorb=True))


# ---------------------------------------------------------------------------
# Tree encoding
# ---------------------------------------------------------------------------


def encode_tree(index: str) -> bytes:
    """Encode a compositional index as a TreeBlob.

    Canonicalisation is the encoder's job, not the caller's. The index is
    recomposed into canonical prefix form here, so equivalent spellings
    produce identical bytes and calling
    :func:`recompose_to_prefix_form` beforehand is redundant rather than
    required. It stays public because two of the format's properties are
    stated in terms of it and because a reader verifying a hash needs to
    name the form the bytes are in. A cell written whole absorbs any
    descendant also written under it, so an index naming a cell and one
    of its descendants encodes the cell -- the region it names.

    Note that :func:`itacart.normalize` is *not* that form: it collapses
    a complete sibling set into its parent, which preserves the region
    while changing the leaf set. Passing a normalised index encodes a
    different, smaller cell set, on purpose.

    Args:
        index: Compositional index string.

    Returns:
        The binary blob.

    Raises:
        InvalidIndexError: If the index is malformed, or if it names row
            1000, which holds the polar caps: this version of the binary
            codecs does not serialize them.

    Example:
        >>> encode_tree("NE(0625/0451(2,1))") == encode_tree("NE(0625/0451(1,2))")
        True
    """
    return _encode_globe(_merge(_parse(index), absorb=True))


def _encode_globe(globe: _Node) -> bytes:
    """Emit the header and one group per quadrant the index names."""
    groups = globe.children
    writer = _BitWriter()
    writer.write_byte_aligned(MAGIC_TREE)
    writer.write_byte_aligned((FORMAT_VERSION << 4) & 0xFF)
    writer.write_byte_aligned(((len(groups) - 1) << 6) & 0xFF)
    for root in groups:
        quadrant = str(root.component)
        count = len(root.children)
        if count > ROOT_MAX_CHILDREN:  # pragma: no cover - 65 536 in one group
            raise InvalidIndexError(f"too many resolution-1 cells: {count}")
        writer.write(_QUADRANT_TO_BITS[quadrant], 2)
        writer.write(count, COUNT_WIDTH_ROOT_CHILDREN)
        for child in root.children:
            _encode_node_bits(writer, child, quadrant)
    return writer.to_bytes()


def _write_component(writer: _BitWriter, node: _Node, quadrant: str) -> None:
    """Write one component. Its validity is already settled.

    Every path here runs through :class:`_Parser` or through the decoder,
    and both check the column, the row and the refinement alphabet before
    a node exists. Repeating those checks would put three copies of one
    rule in the module, and a copy that cannot fire is a copy that cannot
    be kept honest.
    """
    if node.resolution == 1:
        column, row = _pair(node.component)
        writer.write(column, X_BITS)
        writer.write(row, Y_BITS)
        return
    alphabet = _alphabet(node.resolution)
    code = alphabet.index(str(node.component))
    writer.write(code, _component_width(node.resolution))


def _encode_node_bits(writer: _BitWriter, node: _Node, quadrant: str) -> None:
    _write_component(writer, node, quadrant)
    count = len(node.children)
    if count == 0:
        writer.write(0, _count_width_for_a_leaf_at(node.resolution))
        return
    child_res = node.resolution + 1
    limit = _max_children_at(child_res)
    if count > limit:  # pragma: no cover - see below, it cannot fire
        raise InvalidIndexError(
            f"too many children at resolution {child_res}: {count} > {limit}"
        )
    # Siblings are merged by component before this point, so a node has at
    # most one child per token of its level's alphabet: four or twenty-five.
    # The count field holds seven or thirty-one. The guard above is therefore
    # a statement that the field is wide enough, not a defence against input.
    writer.write(count, _count_width_for_children_at(child_res))
    for child in node.children:
        _encode_node_bits(writer, child, quadrant)


# ---------------------------------------------------------------------------
# Tree decoding
# ---------------------------------------------------------------------------


def decode_tree(blob: bytes) -> str:
    """Decode a TreeBlob into a compositional index.

    Round-trips to the canonical form of the original, not necessarily
    to its input spelling.

    Args:
        blob: The binary blob.

    Returns:
        The compositional index string in canonical prefix form.

    Raises:
        MalformedBlobError: If the blob fails structural validation.

    Example:
        >>> decode_tree(encode_tree("NE(0625/0451(2,1))"))
        'NE(0625/0451(1,2))'
    """
    return _render(_decode_to_tree(blob))


def _decode_to_tree(blob: bytes) -> _Node:
    if not isinstance(blob, (bytes, bytearray)):
        raise MalformedBlobError(f"blob must be bytes, got {type(blob).__name__}")
    if len(blob) < HEADER_SIZE_BYTES:
        raise MalformedBlobError("blob is shorter than the tree header")
    if blob[0] != MAGIC_TREE:
        raise MalformedBlobError(f"bad magic byte 0x{blob[0]:02X}")
    version = (blob[1] >> 4) & 0xF
    if version != FORMAT_VERSION:
        raise MalformedBlobError(f"unsupported format version {version}")
    if blob[1] & 0xF:
        raise MalformedBlobError("reserved flag bits are set")
    if blob[2] & 0b0011_1111:
        raise MalformedBlobError("reserved header bits are set")
    globe = _Node(resolution=GLOBE_RESOLUTION, component="")
    reader = _BitReader(bytes(blob), start_bit=HEADER_SIZE_BYTES * 8)
    seen: "set[str]" = set()
    for _ in range(((blob[2] >> 6) & 0b11) + 1):
        quadrant = _BITS_TO_QUADRANT[reader.read(2)]
        if quadrant in seen:
            raise MalformedBlobError(f"quadrant {quadrant} appears twice")
        seen.add(quadrant)
        count = reader.read(COUNT_WIDTH_ROOT_CHILDREN)
        if count == 0:
            raise MalformedBlobError(f"quadrant {quadrant} carries no cell")
        root = _Node(resolution=0, component=quadrant)
        for _ in range(count):
            root.children.append(_decode_node_bits(reader, 1, quadrant))
        globe.children.append(root)
    _check_trailing(reader)
    return globe


def _check_trailing(reader: _BitReader) -> None:
    left = reader.remaining_bits()
    if left >= 8:
        raise MalformedBlobError(f"{left} unread bits after the body")
    if left and reader.read(left) != 0:
        raise MalformedBlobError("trailing padding is not zero")


def _decode_node_bits(reader: _BitReader, expected_res: int, quadrant: str) -> _Node:
    if expected_res == 1:
        column = reader.read(X_BITS)
        row = reader.read(Y_BITS)
        refusal = _row_refusal(row)
        if refusal is not None:
            raise MalformedBlobError(refusal)
        node = _Node(resolution=1, component=(column, row))
    else:
        alphabet = _alphabet(expected_res)
        code = reader.read(_component_width(expected_res))
        if code >= len(alphabet):
            raise MalformedBlobError(
                f"refinement code {code} is reserved at resolution {expected_res}"
            )
        node = _Node(resolution=expected_res, component=alphabet[code])
    count = reader.read(_count_width_for_a_leaf_at(expected_res))
    if count:
        if expected_res >= MAX_RESOLUTION:
            raise MalformedBlobError(f"resolution {MAX_RESOLUTION} admits no children")
        for _ in range(count):
            node.children.append(_decode_node_bits(reader, expected_res + 1, quadrant))
    if expected_res == 1:
        column = _pair(node.component)[0]
        if not _column_is_addressable(column, bool(node.children), quadrant):
            raise MalformedBlobError(
                f"column {column} is not addressable in {quadrant}"
            )
    return node


def validate_tree(blob: bytes) -> None:
    """Validate blob structure, raising on the first problem found.

    Args:
        blob: The binary blob.

    Raises:
        MalformedBlobError: With a message identifying the failure.

    Example:
        >>> validate_tree(encode_tree("NE(0625/0451)")) is None
        True
    """
    _decode_to_tree(blob)


# ---------------------------------------------------------------------------
# Node blobs
# ---------------------------------------------------------------------------


def _single_path(index: str) -> "tuple[str, list[_Node]]":
    globe = _merge(_parse(index), absorb=True)
    if len(globe.children) > 1:
        raise InvalidIndexError("a node blob encodes exactly one cell")
    root = globe.children[0]
    path: list[_Node] = []
    node = root
    while node.children:
        if len(node.children) > 1:
            raise InvalidIndexError("a node blob encodes exactly one cell")
        node = node.children[0]
        path.append(node)
    return str(root.component), path


def encode_node(cell: str) -> bytes:
    """Encode one atomic index as a standalone node.

    Args:
        cell: Atomic index string.

    Returns:
        The node bytes, 10 for a resolution-13 cell.

    Raises:
        InvalidIndexError: If the index is malformed or not atomic, or if
            it names row 1000, which holds the polar caps: this version of
            the binary codecs does not serialize them.

    Example:
        >>> len(encode_node("NE(0625/0451)"))
        4
    """
    quadrant, path = _single_path(cell)
    writer = _BitWriter()
    writer.write_byte_aligned((MAGIC_NODE_HIGH_NIBBLE << 4) | path[-1].resolution)
    writer.write(_QUADRANT_TO_BITS[quadrant], 2)
    for node in path:
        _write_component(writer, node, quadrant)
    return writer.to_bytes()


def decode_node(blob: bytes) -> str:
    """Decode a standalone node.

    Args:
        blob: The node bytes.

    Returns:
        The atomic index string.

    Raises:
        MalformedBlobError: If the bytes are not a valid node.

    Example:
        >>> decode_node(encode_node("NE(0625/0451(3))"))
        'NE(0625/0451(3))'
    """
    quadrant, components = _decode_node_components(blob)
    parts = [quadrant, "("]
    for position, (resolution, component) in enumerate(components):
        if position:
            parts.append("(")
        if resolution == 1:
            column, row = _pair(component)
            parts.append(f"{column:04d}/{row:04d}")
        else:
            parts.append(str(component))
    parts.append(")" * len(components))
    return "".join(parts)


def _decode_node_components(blob: bytes) -> "tuple[str, list[tuple[int, Component]]]":
    resolution = resolution_of_binary(blob)
    reader = _BitReader(bytes(blob), start_bit=8)
    quadrant = _BITS_TO_QUADRANT[reader.read(2)]
    components: "list[tuple[int, Component]]" = []
    for level in range(1, resolution + 1):
        if level == 1:
            column = reader.read(X_BITS)
            row = reader.read(Y_BITS)
            refusal = _row_refusal(row)
            if refusal is not None:
                raise MalformedBlobError(refusal)
            if not _column_is_addressable(column, resolution > 1, quadrant):
                raise MalformedBlobError(
                    f"column {column} is not addressable in {quadrant}"
                )
            components.append((1, (column, row)))
        else:
            alphabet = _alphabet(level)
            code = reader.read(_component_width(level))
            if code >= len(alphabet):
                raise MalformedBlobError(f"reserved refinement code {code}")
            components.append((level, alphabet[code]))
    _check_trailing(reader)
    return quadrant, components


def resolution_of_binary(node: bytes) -> int:
    """Resolution of an encoded node.

    Args:
        node: Encoded node.

    Returns:
        Resolution level, 1 to 13.

    Raises:
        MalformedBlobError: If the bytes are not a valid node.

    Example:
        >>> resolution_of_binary(encode_node("NE(0625/0451(3))"))
        2
    """
    if not isinstance(node, (bytes, bytearray)) or not node:
        raise MalformedBlobError("empty node blob")
    if (node[0] >> 4) != MAGIC_NODE_HIGH_NIBBLE:
        raise MalformedBlobError(f"bad node magic 0x{node[0]:02X}")
    resolution = node[0] & 0xF
    if not MIN_RESOLUTION <= resolution <= MAX_RESOLUTION:
        raise MalformedBlobError(f"node resolution {resolution} out of range")
    return resolution


def prefix_at_resolution_binary(node: bytes, resolution: int) -> bytes:
    """Truncate an encoded node to a coarser resolution.

    Monotone: truncating twice equals truncating once to the coarser of
    the two levels.

    Args:
        node: Encoded node.
        resolution: Target resolution level.

    Returns:
        The truncated node bytes.

    Raises:
        MalformedBlobError: If ``node`` is not a valid node.
        ResolutionError: If ``resolution`` is finer than the node's.

    Example:
        >>> deep = encode_node("NE(0625/0451(3))")
        >>> prefix_at_resolution_binary(deep, 1) == encode_node("NE(0625/0451)")
        True
    """
    quadrant, components = _decode_node_components(node)
    if not MIN_RESOLUTION <= resolution <= len(components):
        raise ResolutionError(f"resolution {resolution} is not a prefix of this node")
    writer = _BitWriter()
    writer.write_byte_aligned((MAGIC_NODE_HIGH_NIBBLE << 4) | resolution)
    writer.write(_QUADRANT_TO_BITS[quadrant], 2)
    for level, component in components[:resolution]:
        _write_component(writer, _Node(level, component), quadrant)
    return writer.to_bytes()


def is_ancestor_binary(parent_node: bytes, child: bytes) -> bool:
    """Test ancestry directly on encoded nodes.

    Avoids decoding, which matters when the check runs inside a smart
    contract or a tight indexing loop.

    When ``child`` is a tree blob the relation is quantified over *every*
    leaf: an ancestor of a tree prefixes the whole tree, not merely some
    branch of it. A parent covering one leaf of a two-leaf tree is not an
    ancestor of that tree.

    Args:
        parent_node: Encoded candidate ancestor.
        child: Encoded candidate descendant, node blob or tree blob.

    Returns:
        ``True`` if the ancestry relation holds for every leaf given.

    Raises:
        MalformedBlobError: If either operand is malformed.

    Example:
        >>> is_ancestor_binary(
        ...     encode_node("NE(0625/0451)"), encode_node("NE(0625/0451(3))")
        ... )
        True
    """
    quadrant, path = _decode_node_components(parent_node)
    if not isinstance(child, (bytes, bytearray)) or not child:
        raise MalformedBlobError("descendant blob must be non-empty bytes")
    if child[0] == MAGIC_TREE:
        return all(_covers(quadrant, path, leaf) for leaf in iter_leaves(child))
    return _covers(quadrant, path, child)


def _covers(
    quadrant: str,
    path: "Sequence[tuple[int, Component]]",
    candidate: bytes,
) -> bool:
    other_quadrant, other_path = _decode_node_components(candidate)
    if quadrant != other_quadrant or len(path) > len(other_path):
        return False
    return list(path) == list(other_path[: len(path)])


# ---------------------------------------------------------------------------
# Streaming
# ---------------------------------------------------------------------------


def iter_leaves(blob: bytes) -> Iterator[bytes]:
    """Stream encoded leaf nodes without decoding to text.

    Args:
        blob: The binary blob.

    Yields:
        Encoded leaf node bytes.

    Raises:
        MalformedBlobError: If the blob fails structural validation.

    Example:
        >>> len(list(iter_leaves(encode_tree("NE(0625/0451(1,2))"))))
        2
    """
    globe = _decode_to_tree(blob)

    def walk(node: _Node, prefix: "list[_Node]", quadrant: str) -> Iterator[bytes]:
        here = [*prefix, node]
        if node.children:
            for child in node.children:
                yield from walk(child, here, quadrant)
            return
        writer = _BitWriter()
        writer.write_byte_aligned((MAGIC_NODE_HIGH_NIBBLE << 4) | here[-1].resolution)
        writer.write(_QUADRANT_TO_BITS[quadrant], 2)
        for step in here:
            _write_component(writer, step, quadrant)
        yield writer.to_bytes()

    for root in globe.children:
        for child in root.children:
            yield from walk(child, [], str(root.component))


def count_vertices(blob: bytes) -> int:
    """Count leaves in a blob without materialising them.

    Args:
        blob: The binary blob.

    Returns:
        Leaf count.

    Raises:
        MalformedBlobError: If the blob fails structural validation.

    Example:
        >>> count_vertices(encode_tree("NE(0625/0451(1,2))"))
        2
    """

    def count(node: _Node) -> int:
        if not node.children:
            return 1
        return sum(count(child) for child in node.children)

    globe = _decode_to_tree(blob)
    return sum(count(child) for root in globe.children for child in root.children)


# ---------------------------------------------------------------------------
# Storage-facing wrappers
# ---------------------------------------------------------------------------


def serialize_to_blob(index: str, compact: bool = True) -> bytes:
    """Encode an index, compacting first by default.

    Convenience wrapper over :func:`encode_tree`; the name callers
    coming from a storage or blockchain context tend to reach for.

    Compaction changes the leaf set, and a blob's identity is its leaf
    set, so two calls differing only in ``compact`` produce different
    blobs on purpose.

    Args:
        index: Compositional index string.
        compact: Apply :func:`itacart.compact_cells` first.

    Returns:
        The binary blob.

    Raises:
        InvalidIndexError: If the index is malformed, or if it names row
            1000, which holds the polar caps: this version of the binary
            codecs does not serialize them. They exist, so this is the one
            refusal here of a cell the grid has.
        NonExistentCellError: If any cell of the index names no cell. A
            blob is a record of cells, so writing one for an address with
            no cell under it would put the fold into storage, where the
            next reader has no way left to tell.
    """
    from .._existence import require_existing_cells
    from ..hierarchy import compact_cells

    require_existing_cells(index)
    return encode_tree(compact_cells(index) if compact else index)


def deserialize_from_blob(blob: bytes) -> str:
    """Decode a blob to the index it stores.

    A TreeBlob stores a region, not a filling: the resolution a fill ran
    at is not recoverable from the bytes because it was never in them.
    Callers needing a uniform resolution hold that number themselves and
    pass it to :func:`itacart.uncompact_cells`.

    Args:
        blob: The binary blob.

    Returns:
        The compositional index string.

    Raises:
        MalformedBlobError: If the blob fails structural validation.
    """
    return decode_tree(blob)
