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

# Canonical Geometry

Two records of the same parcel, written by different hands, agree only if there
is a canonical form to agree in. This page is about getting geometry into one.

```{code-cell} python
import itacart
from itacart import serialization, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
ring = geometry.vertex_to_cell(parcel, resolution=13)
blob = serialization.encode_geometry([ring], geometry_type="POLYGON")
serialization.geometry_hash(blob).hex()[:12]
```

The hash is taken over the canonical form. Encode the same parcel again,
starting the ring at a different vertex, and the hash is the same.

## What canonicalisation changes, and what it must not

A closed ring can be written starting from any of its vertices. Rotating it to
the least starting vertex is a change of **spelling**. Reversing it swaps
interior for exterior — a change of **meaning**.

So canonicalisation rotates and never reverses:

```{code-cell} python
profile = serialization.decode_geometry(blob)
profile["canonical_profile"]
```

```{figure} ../../../_static/f7/f7_10_canonical_rotation.png
:alt: Canonical rotation of a closed ring
:width: 100%

`docs/_static/f7/f7_10_canonical_rotation.png`
```

Ring orientation is preserved exactly, and which ring is the exterior is
carried by part order — the first ring of a component is its exterior — rather
than inferred from winding.

## Three levels of canonical form

They are separate operations and a cadastral record usually needs all three.

**The index.** `normalize` collapses duplicates and orders siblings, so two
spellings of the same cell set become the same text —
[Compaction & Normalization](../index_and_hierarchy/compaction_and_normalization.md).

```{code-cell} python
from itacart import index as ix
ix.normalize("NE(0001/0002(1(A2,A1,A1)))")
```

**The rings.** `canonicalize_rings` normalises a set of rings to a single
spelling per geometry, which is what the encoder applies before hashing.

**The blob.** The GeometryBlob's own canonical profile, declared in its header,
so a consumer knows which normalisation the bytes went through.

## Why a hash over text is not enough

An index string is already readable and already compact, and it is tempting to
hash it directly. Do not, for a stored parcel identity.

A hash over text is a hash over a spelling: whitespace, sibling order,
compaction state and resolution all change the string without changing the
parcel. `normalize` fixes the first two. The blob fixes the rest, because it
encodes the structure rather than its rendering.

And a cover is the wrong thing to hash anyway. It is lossy about the boundary,
so two genuinely different parcels can share a cover —
[Parcel Representation](parcel_representation.md).

## What the hash is for

A canonical geometry hash provides a compact identifier for the exact
geometry represented by a `GeometryBlob`. A registry can publish that
identifier so that a geometry can later be checked against the committed
record, and two authorities that already hold copies of the same parcel
can compare their hashes without exchanging the geometries themselves.

When the hash is recorded in an append-only or otherwise tamper-evident
ledger, a later change to the geometry produces a different hash and can
therefore be detected.

Canonicalization is what makes the hash useful as a stable identifier of
the geometry rather than merely of one particular serialization. Without
a canonical form, geometrically equivalent representations — for example,
the same ring starting at different vertices — may produce different byte
sequences and therefore different hashes. The hash would still verify
those exact bytes, but it could no longer be used directly to establish
that two equivalent encodings represent the same canonical geometry.
