# DataTreeBlob binary encoding

A DataTreeBlob stores **values on cells**: a tree of cells and, at every node
of it, the values of a fixed set of variables. A cell and all of its
ancestors are read in one walk down the tree, on the bytes, and a value is
rewritten in place. Where a [TreeBlob](binary_encoding_spec.md) answers
*where*, this format answers *what is there* at every level at once.

The two formats live side by side and keep separate contracts. The TreeBlob
is canonical and content-addressed: two spellings of one region give one
blob, and its hash is a key. The DataTreeBlob has no canonical form and is not
meant to be hashed: its job is fast access to data on the hierarchical grid,
including the data of the parent cells.

## 1. What the blob holds, and what it does not

A DataTreeBlob holds:

- the nodes of a tree of cells, level by level, from the quadrants down;
- a schema of up to 255 fixed-width variables;
- the value of each variable at each node where the variable lives;
- a validity bitmap for every variable whose absence means "unknown";
- the **observation resolution**, the finest of the data and the unit of
  every count.

It does not hold geometry, vertex order, text per cell, or lists of features
per cell. A categorical attribute is an integer code; the names of the
variables are the only text, and they live in the schema.

## 2. The tree is lexical

The nodes are the cells written into the blob and every prefix of their
spellings, which is the ancestry of the index (see
[Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md)).
Away from the domain border, that is the physical tree, and a node's
aggregate covers the node.

**On the border it is not.** The surplus children of a cell that absorbs the
border are spelled under the next column's prefix, which names no cell, so
they hang under that prefix and not under their physical parent. The prefix
is a node of the blob, because the walk needs it, and it holds the aggregate
of what hangs under it. No function emits it: the decoder leaves it out, and
a path through it skips it.

So the aggregate of a border cell covers its own stem only, and the share
depends on how deep the data goes. With its five physical children at
resolution 2, the trapezoid `NE(2003/0000)` aggregates 80.05% of its area.
With its descendants at resolution 3, it aggregates 77.64%; the remaining
22.39% sits under the prefix `NE(2004/0000)`. The two add up to 100.03%,
because the effective areas of a border cell's descendants add up to a little
more than its own. Both figures at resolution 3 are pinned by
`test_the_border_is_lexical_and_says_so`.

**The polar caps enter.** Row 1000 is in the box of the integer index, and a
cap refines like any other cell. Its own aggregate covers all of it: in the
same test module, the cap `NE(0000/1000)` refined to resolution 5 aggregates
all 1 194 of its descendants. Inside a cap, the same lexical shift happens
at intermediate levels: twenty of the nodes of that tree are prefixes that
name no cell.

## 3. Variables

Each variable is declared with {class}`itacart.DataVariable`.

| Field | Values |
| --- | --- |
| `dtype` | `u8`, `u16`, `u32`, `u64`, `i8`, `i16`, `i32`, `i64`, `f32`, `f64`, coded 0 to 9 in that order |
| `species` | `pyramid` (0): observed at the leaves, aggregated at every internal node. `level` (1): observed at any node, never aggregated |
| `op` | `none` (0), `sum` (1), `mean` (2), `max` (3), `min` (4), `mode` (5). A level variable takes `none` |
| `absent` | `unknown` (0) or `zero` (1): what a value nobody gave means |
| `weight` | required by `mean` and `mode`, refused by the others: an integer pyramid sum declared earlier |
| `leaf_max` | optional, for an unsigned pyramid sum: the most one observation cell can hold |

**Aggregation is local and exact.** An internal node holds its operator over
its daughters' stored values, computed in exact arithmetic and rounded once to
the variable's type, ties to even. A float sum is the exact sum of the
daughters, rounded once; a mean is the exact weighted mean. Rounded at every
level, the mean at the root of the dense case in
`test_the_dense_field_aggregates_exactly` stays within one f32 ulp of the
exact mean of its ten thousand leaves.

**Operators skip what is unknown.** An unknown daughter is left out, a mean or
a mode also leaves out the daughters whose weight is zero, and a node with no
daughter left takes the variable's absence. Operators run over the daughters
the tree holds: a daughter it does not hold is not a zero.

**The denominator of a mean is its weight.** A mean is weighted by a count,
normally the count of the cells where the averaged value was observed. Five
of twenty-five daughters observed at 100 and twenty at 0 average 20; the
same five observed and twenty unknown average 100, over the five. Both are
pinned by `test_the_denominator_of_a_mean_is_the_weight`.

**A pyramid value is given at a leaf only.** At a node with daughters it would
contradict the aggregate, and the encoder refuses it. A parent and its
daughter can both carry values of their own through a level variable, which
is observed wherever it is given.

## 4. Blocks and pruning

**A leaf above the observation resolution stands for its block.** Reading any
cell beneath such a leaf returns the leaf's values, except that a sum is
divided evenly over the cells: an integer sum must divide exactly by the
number of observation cells in the block, and a float sum is divided and
rounded once.

**A block is a regular cell.** It must be a cell that does not absorb the
border, so that the cells beneath it are exactly the product of the
refinement ratios. A leaf that absorbs the border, or an integer sum that does
not divide evenly, is refused.

**Pruning is optional and off by default.** With `prune=True`, the encoder
folds a complete family of leaves into its parent, bottom-up, so that a
block can join its sisters into a larger one. It folds only when:

- every daughter the ratio allows is present and is a leaf;
- the daughters agree on every variable, value and validity alike;
- the parent is a regular cell;
- no level variable has values at the daughters' level, because those
  values would be lost.

Folding by one variable alone would destroy the others: a family fully covered
by obstacles of two heights stays a family.

## 5. Wire format

Every multi-byte field is little-endian.

### 5.1 The frame

| Bytes | Field | Values |
| --- | --- | --- |
| `0x00` | Magic | `0xC9` |
| `0x01` hi | Version | `1` |
| `0x01` lo | Flags | bit 0: the writer asked for pruning; bits 1 to 3 reserved and zero |
| `0x02` | Observation resolution | 0 to 13 |
| `0x03` | Deepest level `K` | 0 to the observation resolution |
| `0x04` | Variable count `V` | 1 to 255 |

Then one record per variable, in schema order:

| Bytes | Field |
| --- | --- |
| 1 | Name length `L`, then `L` bytes of UTF-8 |
| 1 | Type code |
| 1 | Species code |
| 1 | Operator code |
| 1 | Absence code |
| 1 | Index of the weight variable, or `0xFF` |
| 8 | `leaf_max`, or 0 for none |
| 2 | Presence mask: bit `k` set when the variable has values at level `k` |
| 1 per level present | Type code of the values at that level |

A pyramid variable has values at every level from 0 to `K`. A level variable
has them at the levels where it was observed. An unsigned sum with `leaf_max`
is written at each level in the narrowest unsigned type that holds
`leaf_max` times the most observation cells a node of that level covers, so
an update in place never overflows it. Every other variable is written in its
declared type.

The frame continues with `K + 1` node counts (`u32`), the number of sections
(`u16`) and their absolute offsets (`u32` each), and is padded with zeros to
a multiple of eight bytes.

### 5.2 The sections

Each section starts on a multiple of eight bytes and is padded with zeros. In
order:

| Section | Content |
| --- | --- |
| Quadrants | `n₀` bytes: the quadrant codes, ascending (NE 0, NW 1, SE 2, SW 3) |
| Resolution 1 | `n₁` × `u32`: the box positions of the resolution-1 nodes, ascending, the quadrant included |
| For `k` = 0 to `K − 1`: directory | `(n_k + 1)` × `u32`: the index of each node's first daughter at level `k + 1`, and the count of that level |
| For `k` ≥ 1: masks | `n_k` × `u8` when the daughters are at an even resolution (4 possible), `u32` at an odd one (25): bit `d` set when daughter `d` of the alphabet is present |
| For `k` ≥ 1: cell flags | `⌈n_k / 8⌉` bytes: bit set when the node names a cell |
| For `k` = 0 to `K`, each variable present | `n_k` values in the level's type, then, for an `unknown` absence, a validity bitmap of `⌈n_k / 8⌉` bytes |

Bits are numbered from the least significant bit of byte 0: node `i` is bit
`i mod 8` of byte `⌊i / 8⌋`. An unknown value has its bit cleared and holds
zero in its slot, so every type keeps its whole domain. A level holds its
nodes in canonical order, which is the order of the integer index.

## 6. Reading and writing

**The walk.** To read a cell at resolution `r`, find its quadrant, binary
search its resolution-1 position, and then take one step of constant time per
level:

```
iₖ = first[k − 1][iₖ₋₁] + popcount(mask[k − 1][iₖ₋₁] & (2ᵈ − 1))
```

where `d` is the cell's code at level `k`. Every step lands on an ancestor, so
the values of the parents come at no extra cost: one slot per level in each
vector.

**On the bytes.** The reader casts the sections to typed views and decodes
nothing else. A `bytes` blob cannot change, so its frame is parsed on the
first read and kept for the next ones; the last sixteen are kept. A
`bytearray` is parsed on every call. A host that is not little-endian reads
through `struct` instead of casting.

**In place.** {func}`itacart.update_data_tree` writes a value at a leaf, or a
level value at any node, and recomputes the aggregate of every ancestor from
its daughters. Every slot is computed before any is written, and a refusal
writes nothing. Because the widths follow the structure and not the values,
an update in place on a blob encoded without pruning produces exactly the
bytes a new encoding of the new data would, as
`test_an_update_in_place_is_the_encoding_of_the_new_data` pins. A pruned blob
stays valid and holds the right values, but an update can leave a family that a
new pruned encoding would fold: the format has no canonical form, and pruning
is the writer's choice. Adding or removing a node changes the structure, and is
a new encoding.

**Slices.** The descendants of a node form one contiguous run of every level,
so {func}`itacart.slice_data_tree` copies runs. The ancestors on the path keep
their level values, and their pyramid values become the aggregate of the
slice alone.

## 7. Size and speed

The two cases below are built in `tests/unit/test_data_tree_blob.py`. Each has
ten thousand leaves, and their sizes are pinned there.

| Case | Variables | Bytes | Structure | Values | Validity | Total |
| --- | --- | --- | --- | --- | --- | --- |
| Dense field, resolution 3 to 7 | count, f32 mean, mode, level | 70 808 | 3.1 | 50.9 | 2.3 | 56.6 bits per leaf |
| Drone map, resolution 5 to 9 | four counts, max height, level | 73 632 | 3.2 | 50.7 | 4.4 | 58.9 bits per leaf |
| The same map, pruned | | 11 168 | 1.2 | 6.4 | 0.7 | 8.9 bits per leaf |

- **The values dominate.** The structure costs about three bits per leaf: a
  mask and a directory entry per internal node, against the 8.4 bits per leaf
  a TreeBlob spends on structure alone.
- **Pruning is where large maps shrink.** Free space and the inside of a
  building fold into blocks.
- **Compression belongs to the transport.** With zlib, the drone map
  compresses about thirty to one and the dense field about two to one.

**Speed is a reading, not a pinned figure.** Measured in pure Python on the
development machine, and to be read as an order of magnitude:

- the walk alone costs about half a microsecond per level;
- reading one cell by its integer costs about 15 µs, the existence check of
  the input included, and by its spelling about 30 µs;
- an update in place costs a few hundred microseconds, most of it in parsing
  the frame of a `bytearray`.
