# GeometryBlob binary encoding

GeometryBlob is the binary form of an **ordered** ITACaRT vertex
sequence. It preserves what TreeBlob discards — vertex order, ring
topology, edge model, densification parameters and OGC type — which
makes it the form to hash when a geometry must be *identified* rather
than merely covered.

## 1. Header

Nine bytes, fixed.

| Bytes | Bits | Field | Values |
| --- | --- | --- | --- |
| `0x00` | 8 | Magic | `0xC8` |
| `0x01` | 8 | Format version | `0x01` |
| `0x02` hi | 4 | Canonical profile | `0x0` MIN_LEX_CYCLIC_ROTATION |
| `0x02` lo | 4 | Edge model | `0x0` PLANAR_SINUSOIDAL_STRAIGHT, `0x1` WGS84_GEODESIC |
| `0x03` hi | 4 | Densification model | `0x0` NONE, `0x1` ORTHODROMIC_VINCENTY, `0x2` ASSUMED_BY_PRODUCER |
| `0x03` lo | 4 | Resolution mode | `0x0` UNIFORM, `0x1` MULTI |
| `0x04`–`0x05` | 16 | `max_segment_m` | 0–65 535 |
| `0x06` hi | 4 | Uniform resolution | 1–13 |
| `0x06` lo | 4 | Geometry type | OGC SFA, see below |
| `0x07`–`0x08` | 16 | Component count | 1 for single types, N for MULTI\* |

### Geometry type codes

Codes follow OGC Simple Feature Access, where **code 0 is reserved**.

| Code | Type |
| --- | --- |
| `0x0` | reserved |
| `0x1` | POINT |
| `0x2` | LINESTRING |
| `0x3` | POLYGON |
| `0x4` | MULTIPOINT |
| `0x5` | MULTILINESTRING |
| `0x6` | MULTIPOLYGON |
| `0x7` | GEOMETRYCOLLECTION, reserved |

Numbering this enum from zero shifts every type by one. The shift is
silent: the blob still decodes, into the wrong type, and its identity
hash still verifies against itself. The codes are pinned by a named
test for that reason.

## 2. Body

After the header, for each component: a 16-bit part count, then for each
part a 4-bit ring role, a reserved 4-bit nibble that must be zero, a
16-bit vertex count, and the vertex stream.

Ring roles are `0x0` EXTERIOR, `0x1` INTERIOR, `0x2` LINE, `0x3` POINT.
The role is recoverable from the geometry type and the part index, so
the decoder reads it and does not depend on it.

Structure per type:

| Type | Components | Parts per component |
| --- | --- | --- |
| POINT | 1 | 1 part, 1 vertex |
| LINESTRING | 1 | 1 part |
| POLYGON | 1 | exterior first, then holes |
| MULTIPOINT | N | 1 part, 1 vertex each |
| MULTILINESTRING | N | 1 part each |
| MULTIPOLYGON | N | exterior first, then that polygon's holes |

MULTIPOLYGON is the one type whose argument nests one level deeper — a
list of polygons, each a list of rings — because a flat list of rings
cannot say which hole belongs to which polygon. Flattening it would
attach every hole to the first polygon and produce a blob that decodes
cleanly into the wrong geometry, so the flat spelling is refused.

## 3. Vertex stream

Vertices are delta-encoded against their predecessor. Vertex 0 is
written in full:

- 2 bits quadrant
- 11 bits column, 10 bits row, written without offset
- one refinement field per level: 2 bits at even resolutions, 5 at odd

Every later vertex writes a 4-bit `shared_levels` field first — the
number of hierarchy levels it has in common with its predecessor — and
then only the levels below that. A `shared_levels` of zero means no
sharing and the vertex is written in full.

At resolution 13 a fully-written vertex costs 65 bits; a vertex sharing
all thirteen levels costs the 4-bit field alone.

## 4. Structural preconditions

The encoder refuses these rather than encoding them:

| Type | Requirement |
| --- | --- |
| POINT | exactly one vertex |
| LINESTRING | at least two vertices, no consecutive duplicates |
| POLYGON, each ring | at least three vertices, no closing repetition, no consecutive duplicates |
| MULTI\* | each component satisfies its single-geometry rule |

**Ring closure is implicit.** The last vertex does not repeat the first,
and an explicitly closed ring is refused. Accepting one and then
rotating it to canonical form moves the repeated vertex into the middle
of the sequence, which turns a closed ring into one that is neither
closed nor valid.

The encoder does **not** validate cadastral admissibility. Whether a
given record type may carry a POLYGON is a question for the processing
layer; GeometryBlob is a data structure without domain knowledge.

## 5. Compatibility matrix

Two combinations are contradictions and are refused:

| Edge model | Densification | Verdict |
| --- | --- | --- |
| WGS84_GEODESIC | NONE | refused — claims a curve preserved by a producer that ran nothing |
| PLANAR_SINUSOIDAL_STRAIGHT | ORTHODROMIC_VINCENTY | refused — the edges are straight |
| WGS84_GEODESIC | ORTHODROMIC_VINCENTY or ASSUMED_BY_PRODUCER | accepted |
| PLANAR_SINUSOIDAL_STRAIGHT | NONE or ASSUMED_BY_PRODUCER | accepted |

## 6. Canonicalisation and identity

The declared profile is minimum lexicographic cyclic rotation: closed
rings are rotated to their least starting vertex, so the same ring
written from a different vertex canonicalises to the same sequence.
{func}`itacart.canonicalize_rings` performs it.

**Winding is not touched.** Rotating is a change of spelling; reversing
would be a change of meaning. The header declares one transformation and
the encoder applies exactly that one, so a ring and its reverse are two
geometries and hash apart.

Which ring is the exterior and which are holes is carried by **part
order** -- the first ring of a component is its exterior -- and not by
reading the direction. Direction is preserved to the byte all the same,
because many producers and consumers use it to carry that same
distinction, and a ring handed on reversed may be read downstream as the
opposite kind of ring. The common convention is counter-clockwise for an
exterior boundary and clockwise for a hole:

```
exterior = [A, B, C, D]   # counter-clockwise
hole     = [E, H, G, F]   # clockwise: E, F, G, H reversed
```

Nothing here reverses a ring to make it agree with that convention, and
nothing here checks that it does.

{func}`itacart.geometry_hash` returns `keccak256` of the blob. Not
SHA3-256: the two differ in their padding byte, and substituting one for
the other would silently change every identity. The digest is meaningful
only for canonicalised blobs — without canonicalisation the same
geometry written from a different starting vertex hashes differently,
which is the whole point of the profile.

## 7. Bridge to TreeBlob

{func}`itacart.geometry_to_tree` decodes the blob, takes the set of
vertex indices in order of first appearance with duplicates removed, and
encodes that set as a TreeBlob.

The derivation is one-way by construction. Order, ring topology, edge
model and geometry type all vanish, so distinct geometries over the same
vertices yield one TreeBlob. Coverage survives; identity does not. This
is a property the format is built on, and it is stated as a test that
asserts the loss rather than one that regrets it.

A TreeBlob holds one quadrant, so a geometry whose vertices span more
than one is refused rather than truncated.

## 8. Density

Density is a curve in the depth of the shared prefix. On a synthetic
resolution-13 polygon whose vertices share levels 1 to 7, the amortised
cost is about 25 bits per vertex for a uniform blob and about 29 for a
multi-resolution one, easing towards those figures as vertex count
rises. The extremes are structural: 4 bits per vertex when all thirteen
levels are shared, 65 when none are.

For cadastral polygons in the 500 to 1000 vertex range, that puts blobs
between roughly 1.6 and 3.6 kilobytes.

## 9. Provenance

Ported from `itacart_core/geometry_blob.py`. Four things changed.

- The **index space** is zero-based on both axes, for the reason given
  in the TreeBlob specification.
- The reference carries its own bit reader and writer, its own
  densification and its own geodesics, duplicating code this package
  already has. The port imports {func}`itacart.densify_segment`,
  {func}`itacart.canonicalize_rings` and the geodesy module instead.
- The reference's canonicalisation also **enforces ring orientation** by
  signed area, which its header does not declare. That step is not
  ported: winding is meaning, not spelling.
- MULTIPOLYGON nesting is expressed rather than refused.
