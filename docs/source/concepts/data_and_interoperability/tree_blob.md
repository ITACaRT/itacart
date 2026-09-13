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

# TreeBlob

A cell set, encoded as a prefix tree. It answers *where*, compactly, and it
answers some questions without being decoded at all.

```{code-cell} python
import itacart
from itacart import serialization, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
cover = geometry.polyfill(parcel, resolution=7)
blob = serialization.encode_tree(cover)
itacart.count_cells(cover), len(blob)
```

## Why a tree and not a list

A compositional index is already a tree — shared prefixes with branches — and
the encoding keeps that structure instead of flattening it. Cells under one
parent share their parent's bytes once, so a dense region costs far less than
the same count of scattered cells.

This is the same saving compaction gives at the text level, made structural.
Compacting before encoding, which `serialize_to_blob` does by default,
compounds the two.

## Reading without decoding

```{code-cell} python
serialization.count_vertices(blob)
```

`iter_leaves` streams encoded leaf nodes in the same way, so a large set can be
walked without ever building its strings. For a single node, the structural
questions are answered on the bytes:

```{code-cell} python
node = serialization.encode_node(itacart.decompose(cover)[0])
(serialization.resolution_of_binary(node),
 serialization.is_ancestor_binary(
     serialization.prefix_at_resolution_binary(node, 4), node))
```

## Child counts are not refinement ratios

The format does not assume four or twenty-five children. It stores the count it
found, because at the domain border a cell holds fewer — between two and six —
and an encoder that multiplied instead of counting would produce a blob that
decodes to cells which do not exist.

That is [Absorbing Boundary Cells](../boundaries/absorbing_boundary_cells.md)
reaching into the wire format, and the specification gives it a section of its
own.

## One quadrant per blob

A TreeBlob holds one quadrant. A geometry whose vertices span more than one
needs more than one blob — which is a consequence of the quadrant being the
root of the index, and of the grid not wrapping:
[Antimeridian](../boundaries/antimeridian.md).

## Validation

```{code-cell} python
serialization.validate_tree(blob) is None
```

Raises on the first structural problem rather than returning a verdict. A blob
that survives validation decodes to a set of cells that exist.

The wire format, the index space, the density figures and the relationship to
GeometryBlob are specified in
[TreeBlob binary encoding](binary_encoding_spec.md).
