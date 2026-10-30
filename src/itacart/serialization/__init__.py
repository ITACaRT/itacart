"""Binary encodings of ITACaRT indices.

Three formats, with different invariants:

**TreeBlob** (:mod:`itacart.serialization.tree_blob`) encodes the
compositional tree as a *set* of cells in canonical prefix-index form. It
is content-addressable and reaches 10 bytes for a resolution-13 leaf.
Order, ring topology and edge model are all discarded, which makes it the
right structure for coverage, ancestry and indexing queries, and the wrong
one for identity. The design generalises beyond ITACaRT: any DGGS with a
prefix-structured hierarchical index can adopt it.

**GeometryBlob** (:mod:`itacart.serialization.geometry_blob`) encodes an
*ordered sequence* of vertices with edge model, ring topology,
densification parameters and OGC SFA type. It preserves everything
TreeBlob drops, and is the form to hash when a geometry has to be
identified rather than merely covered.

The two are bridged by :func:`~itacart.serialization.geometry_blob.geometry_to_tree`,
in one direction only: a GeometryBlob determines a TreeBlob, but the same
TreeBlob arises from many geometries.

**DataTreeBlob** (:mod:`itacart.serialization.data_tree_blob`) encodes
*values on cells*: a tree of cells with a fixed set of variables at every
node, aggregated from the leaves or observed at any level. It is read and
updated in place, a cell and its ancestors in one walk, and it has no
canonical form: it is the format for data with a value at every level,
where the TreeBlob is the format for a region.

Provenance: ``itacart_core/binary_index.py`` and ``geometry_blob.py``.
"""

from __future__ import annotations

from .data_tree_blob import (
    DataVariable,
    data_tree_schema,
    decode_data_tree,
    encode_data_tree,
    read_data_cell,
    read_data_path,
    slice_data_tree,
    update_data_tree,
)
from .geometry_blob import (
    decode_geometry,
    encode_geometry,
    geometry_hash,
    geometry_to_tree,
    read_geometry_type,
    validate_geometry,
)
from .tree_blob import (
    count_vertices,
    decode_node,
    decode_tree,
    deserialize_from_blob,
    encode_node,
    encode_tree,
    is_ancestor_binary,
    iter_leaves,
    prefix_at_resolution_binary,
    recompose_to_prefix_form,
    resolution_of_binary,
    serialize_to_blob,
    validate_tree,
)

__all__ = [
    "encode_tree",
    "decode_tree",
    "encode_node",
    "decode_node",
    "recompose_to_prefix_form",
    "serialize_to_blob",
    "deserialize_from_blob",
    "is_ancestor_binary",
    "prefix_at_resolution_binary",
    "resolution_of_binary",
    "validate_tree",
    "iter_leaves",
    "count_vertices",
    "encode_geometry",
    "decode_geometry",
    "validate_geometry",
    "geometry_hash",
    "geometry_to_tree",
    "read_geometry_type",
    "DataVariable",
    "encode_data_tree",
    "decode_data_tree",
    "read_data_path",
    "read_data_cell",
    "slice_data_tree",
    "update_data_tree",
    "data_tree_schema",
]
