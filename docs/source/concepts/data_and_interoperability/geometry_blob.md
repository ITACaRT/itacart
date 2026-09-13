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

# GeometryBlob

The form to use when a geometry must be **identified** rather than merely
covered. It carries ordered vertices, ring topology, the edge model, the
densification parameters and the OGC type — everything that makes two encodings
of the same parcel the same thing.

```{code-cell} python
from itacart import serialization, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
ring = geometry.vertex_to_cell(parcel, resolution=13)
blob = serialization.encode_geometry([ring], geometry_type="POLYGON")
len(blob), serialization.read_geometry_type(blob)
```

`read_geometry_type` answers from the nine-byte header without decoding the
body — useful when you are routing blobs rather than reading them.

## What the header declares

```{code-cell} python
profile = serialization.decode_geometry(blob)
{k: profile[k] for k in profile if k != "rings"}
```

Every one of those fields is a decision that a consumer cannot recover from
coordinates alone. [Edge Semantics](edge_semantics.md) is the page for the
first two.

## Orientation is meaning; rotation is spelling

Ring orientation is significant and is preserved exactly. Canonicalisation may
rotate a closed ring to its least starting vertex — it never reverses one,
because rotating is a change of spelling and reversing is a change of meaning.

```{figure} ../../../_static/f7/f7_10_canonical_rotation.png
:alt: Canonical rotation of a closed ring
:width: 100%

`docs/_static/f7/f7_10_canonical_rotation.png`
```

Which ring is the exterior and which are holes is carried by **part order**:
the first ring of a component is its exterior. That is a structural convention
rather than something inferred from winding, and it is why the format can hold
multi-part geometries without ambiguity.

## Identity, and hashing

```{code-cell} python
serialization.geometry_hash(blob)[:12]
```

The hash is taken over the canonical form, which is what makes it an identity
rather than a checksum of one particular spelling. Two records of the same
parcel, encoded independently, hash the same — and that is the property a
registry or a ledger needs.

[Canonical Geometry](../cadastral_workflows/canonical_geometry.md) is the
cadastral treatment of the same idea.

## The bridge to TreeBlob, and why it is one-way

```{code-cell} python
tree = serialization.geometry_to_tree(blob)
len(blob), len(tree)
```

`geometry_to_tree` derives the TreeBlob of a GeometryBlob's vertex set. The
derivation is **one-way by construction**: order, ring topology, edge model and
geometry type all vanish, so distinct geometries over the same vertices yield
one TreeBlob.

Coverage survives; identity does not. This is a property the format is built
on, and the package states it as a test that asserts the loss rather than one
that regrets it.

Note that the derived blob is not necessarily smaller — above it is larger. The
two forms differ in what they can answer, not in how much they cost.

So: keep the GeometryBlob if you will ever need to say *which* geometry. Derive
the TreeBlob when you only need to say *where*.

The wire format is specified byte by byte in
[GeometryBlob binary encoding](geometry_encoding_spec.md).
