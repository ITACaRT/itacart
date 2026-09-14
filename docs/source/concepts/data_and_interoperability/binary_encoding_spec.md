# TreeBlob binary encoding

TreeBlob is the dense binary form of a compositional index. It stores a
**set of cells** — a region — and nothing else.

It is compact because it stores that region as a **tree** rather than as a
list, which is what a discrete global grid makes possible: cells sharing a
prefix share its bytes, and a prefix is written once however many leaves
hang from it. A contiguous cover costs around **8.4 bits per cell** —
against the fifty-five characters one resolution-13 index takes as text —
and the floor of seven bits per leaf is a property of the format, not of a
lucky region. Section 7 is that curve, and the cost of the frame that
carries it.

The format exists so that a region can be identified by its bytes: two
inputs describing the same leaf set produce byte-identical blobs, and a
content hash over those bytes is therefore a stable key.

## 1. What the blob holds, and what it does not

A TreeBlob holds the quadrant, the set of leaf cells, and the tree
structure that reaches them. It does not hold vertex order, ring
topology, or the resolution some fill was run at.

That last omission is deliberate and worth stating plainly, because it
is easy to mistake for a gap. A blob is mixed-resolution by
construction: a compacted region carries leaves at many levels. "The
fill resolution" is a property of the *procedure* that produced a set,
not of the set. The greatest resolution among the leaves is a lower
bound on any fill resolution and never the value itself, so no decoder
can recover it. Callers who need a uniform expansion hold that number
themselves and pass it to {func}`itacart.uncompact_cells`.

## 2. Canonical form

Two spellings of one leaf set must produce one blob, so the encoder
rewrites its input before emitting bytes.

*Walk order* is the spelling a contour traversal produces: a prefix is
shared only with the previous leaf along the boundary, so the same
component can appear more than once among siblings. *Prefix-index form*
is canonical: every distinct prefix appears once, and siblings are
ordered.

{func}`itacart.recompose_to_prefix_form` performs that rewrite and is
public, because two of the format's properties are stated in terms of
it.

**The encoder canonicalises; the caller need not.** {func}`itacart.encode_tree`
recomposes its input before emitting bytes, so calling the rewrite first
is redundant rather than required. Nothing is lost by doing it, and
nothing is gained.

It is **not** {func}`itacart.normalize`. Normalisation collapses a
complete sibling set into its parent, which preserves the region while
changing the leaf set. A blob's identity is its leaf set, so the codec
needs the form that reorders without compacting. The two canonical forms
answer different questions and have different fixed points:

```
NE(0625/0451(2,1,4,3))  ->  recompose  ->  NE(0625/0451(1,2,3,4))
                        ->  normalize  ->  NE(0625/0451)
```

## 3. Wire format

A **tree blob** encodes a whole compositional tree as one `bytes` value.

| Bytes | Bits | Field | Values |
| --- | --- | --- | --- |
| `0x00` | 8 | Magic | `0xC7` |
| `0x01` hi | 4 | Format version | `0x1` |
| `0x01` lo | 4 | Flags | reserved, must be zero |
| `0x02` hi | 2 | Quadrant group count minus one | 0–3, so one to four groups |
| `0x02` lo | 6 | Reserved | must be zero |

Then one **quadrant group** per quadrant the index names, bit-packed
from byte 3 onwards:

| Bits | Field | Values |
| --- | --- | --- |
| 2 | Quadrant | bit 0 west, bit 1 south |
| 16 | Resolution-1 child count | at least 1 |
| … | Bit-packed children, pre-order | variable |

Groups are ordered by quadrant code and no quadrant may appear twice, so
a cell set has one spelling in bytes however it was written.

An index naming more than one quadrant is ordinary rather than exotic:
{func}`itacart.compose` produces that form for any cell set crossing the
equator or the meridian, and {func}`itacart.polyfill` produces it for
any region over the equator. A codec that held one quadrant could not
encode the canonical output of either.

Packing is MSB-first within every byte, and the body is padded to a byte
boundary with zeros. Trailing padding that is not zero is a malformed
blob, as is any whole unread byte after the body.

Each node writes its component, then a **child-count field** whose width
follows the resolution of its children. A count of zero marks a leaf; a
leaf still writes the field, at the width its children would have used.

| Level | Component | Count field for children |
| --- | --- | --- |
| Resolution 1 | 11 bits column, 10 bits row | 3 bits |
| Even resolution | 2 bits | 5 bits |
| Odd resolution | 5 bits | 3 bits |

A **node blob** encodes one root-to-cell path and is the form
{func}`itacart.iter_leaves` returns:

| Bytes | Field |
| --- | --- |
| `0x00` | `0xA0` in the high nibble, resolution in the low nibble |
| … | 2 bits quadrant, then the components of each level |

## 4. Index space

Columns and rows are written as they appear in the index, **without
offset**. Both axes are zero-based.

| Axis | Range | Field |
| --- | --- | --- |
| Column | 0–2003 as cells, 2004 as a prefix | 11 bits |
| Row | 0–999 | 10 bits |

Two structural exceptions, both constant, neither dependent on latitude:

**Column 0 is the meridian triangle.** It is one cell, addressed from
the east, so it exists in the NE and SE quadrants and is absent from NW
and SW at every resolution. A codec that refuses column 0 outright is
right in the west for the wrong reason and wrong in the east.

**Column 2004 is the trapezoid exception.** The grammar admits it so
that the fifth child of a trapezoidal cell in column 2003 can be
spelled. It is a prefix only, never a resolution-1 cell, so it is
accepted refined and refused bare.

Whether a column exists *at a given row* is a different question, since
the last column of a row retreats with the cosine of the latitude. That
is what {func}`itacart.is_valid_cell` answers. This codec validates the
wire format and does not re-derive the boundary, because a second copy
of the boundary model inside the serializer is a second copy that can
disagree with the first.

## 5. Child counts are not refinement ratios

{func}`itacart.refinement_ratio` describes a node in the interior of the
grid. A trapezoidal parent yields one child more than the ratio, and the
surplus does not attenuate with depth. Reading the ratio as a hard bound
refuses the whole eastern border.

The count field is therefore bounded by its own width — seven at an even
child resolution, thirty-one at an odd one — and not by the ratio. Since
siblings are merged by component before encoding, a node has at most one
child per token of its level's alphabet, four or twenty-five, so the
field is always wide enough.

## 6. Properties

Eight properties hold, and each is pinned by a named test in
`tests/unit/test_serialization.py`.

1. **Determinism.** The same index encodes to the same bytes.
2. `decode(encode(x)) == recompose(x)`.
3. `encode(x) == encode(recompose(x))`.
4. **Recomposition is idempotent.**
5. **Round trip is bit-exact**: re-encoding a decoded blob reproduces it.
6. **Content addressing.** Any two spellings of one leaf set agree.
7. **A resolution-13 node is ten bytes.** The property is about
   {func}`~itacart.serialization.tree_blob.encode_node`, the standalone node, and
   not about a blob wrapping it.
8. **Prefix truncation is monotone**: truncating twice equals truncating
   once to the coarser of the two levels.

## 7. Density

Density is a curve in how much prefix consecutive leaves share, not a
constant.

Every leaf writes its own component and its own terminating count field,
so the floor is **seven bits per leaf** — two bits of component plus
five of count, at an even resolution. No tree of any shape goes below
it. The densest trees the format admits, complete refinement with no
gaps, measure about 8.4 bits per leaf. A contiguous fill converges from
above towards roughly nine as the region grows. Leaves scattered across
many resolution-1 cells cost upwards of a hundred, because each pays for
its own path from the root.

Any figure below seven bits per leaf is unreachable and describes
something other than this format.

### What the frame costs

Those figures are the cost of the **body**. A blob also carries a fixed
frame: three header bytes — magic, version and flags, group count — and,
per quadrant group, two bits of quadrant plus a sixteen-bit resolution-1
child count. The frame is paid once, whatever it wraps.

For a single leaf the frame dominates. Measured against this
implementation, a blob holding one resolution-13 leaf is **twenty bytes**,
of which ten are the node; at resolution 1 it is nine bytes around a
four-byte node. The crossover is quick — twenty-five complete siblings
cost 11.2 bits per leaf, and a contiguous fill of forty-five thousand
cells 8.39, which is the complete-refinement figure above reached from a
real cover.

So a blob is the wrong container for one cell and the right one from a
handful upwards. Where a single cell has to travel alone, the node is the
smaller object and carries the same address.

Unlike the seven-bit floor, which follows from the format, the twenty
bytes is a reading of this implementation and no test pins it.

### The node against its information content

Ten bytes is close to the floor rather than far from it. A resolution-13 address
carries about **62.8 bits** of information: two of quadrant, twenty-one of
the resolution-1 pair, twelve across the six even levels and twenty-eight
across the six odd ones. Eighty bits is **1.27 times** that, and the
difference is byte alignment plus the count field every node writes. The
index string, at fifty-five characters, is seven times it.

## 8. Relationship to GeometryBlob

{func}`itacart.geometry_to_tree` derives a TreeBlob from a GeometryBlob.
The derivation is one-way by construction: the result is a function of
the vertex **set** alone, so distinct vertex orders, ring topologies,
edge models and geometry types over the same vertices all yield the same
TreeBlob. Coverage survives; identity does not.

## 9. Provenance

The encoding is a port of the reference implementation in
`itacart_core/binary_index.py`. Three of its behaviours were corrected
rather than carried over, because the reference models a partial grid —
parallelograms only, without the meridian triangle or the trapezoidal
cell's surplus child, a limitation its own exception module states.

- Its index space runs from one on both axes, where this one runs from
  zero; carrying the offset would make every column and row encode as
  its neighbour.
- It refuses column 2004, which makes the trapezoid's fifth child
  inexpressible.
- It bounds child counts by the refinement ratio, which refuses the
  eastern border.
- Its grammar admits one quadrant per blob, which refuses every region
  that crosses the equator.

The format version was raised to `0x1` to mark the break. Blobs written
under version `0x0` are refused rather than read as if they agreed.
