---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Serialization

Two binary formats, for two different questions. **TreeBlob** stores which
cells. **GeometryBlob** stores a geometry's identity — its vertices, ring
topology, edge model and type.

```{code-cell} python
import itacart
from itacart import serialization, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
cover = geometry.polyfill(parcel, resolution=7)
blob = serialization.serialize_to_blob(cover)
itacart.count_cells(cover), len(blob)
```

Four cells in seventeen bytes, and the round trip is exact:

```{code-cell} python
serialization.deserialize_from_blob(blob) == itacart.normalize(itacart.compact_cells(cover))
```

Note what it compares against. `serialize_to_blob` compacts before encoding by
default, so what comes back is the compacted, normalised form of what went in —
the same ground, spelled as small as it goes.
[Compaction & Normalization](../index_and_hierarchy/compaction_and_normalization.md)
is why that is the right default.

## Why not just store the string

A compositional index is already compact and already readable, and for most
purposes it is the right thing to store. The binary forms earn their place in
three situations.

**When the answer is needed without decoding.** Structural questions can be
answered on the bytes:

```{code-cell} python
node = serialization.encode_node(itacart.decompose(cover)[0])
(len(node),
 serialization.resolution_of_binary(node),
 serialization.is_ancestor_binary(
     serialization.prefix_at_resolution_binary(node, 5), node))
```

Seven bytes, and ancestry tested without ever building a string.

**When the set is large.** `count_vertices` and `iter_leaves` work on the blob
without materialising the cells:

```{code-cell} python
serialization.count_vertices(blob)
```

**When it has to be hashed.** A hash over text is a hash over a spelling. A
hash over a canonical blob is a hash over the thing.

## Validation is explicit

```{code-cell} python
serialization.validate_tree(blob) is None
```

`validate_tree` and `validate_geometry` raise on the first structural problem
they find rather than returning a verdict. A malformed blob is a
`MalformedBlobError`; an impossible combination of header fields is an
`IncompatibleProfileError`. Both descend from `SerializationError`, which
descends from the package root.

## The specifications

The wire formats are specified in full, byte by byte, and those specifications
are the authority rather than this page:

- [TreeBlob binary encoding](binary_encoding_spec.md)
- [GeometryBlob binary encoding](geometry_encoding_spec.md)

Read them before writing an encoder in another language, before changing a
header field, or when a blob from elsewhere does not decode.
