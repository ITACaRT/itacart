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

# DataTreeBlob

Values on cells, held at every level of the tree. A TreeBlob answers *where*;
a DataTreeBlob answers *what is there*, at a cell and at each of its
ancestors at once.

The example on this page is an obstacle map for a drone. A ten-metre cell at
ITA is refined to its hundred one-metre cells. A building 24 m tall covers
thirty-five of them, heights are kept in decimetres, and ten cells were never
surveyed.

```{code-cell} python
import itacart
from itacart import DataVariable

block = itacart.geo_to_cell(-45.8722, -23.2124, 7)
leaves = list(itacart.get_descendants(block, 9))

variables = [
    DataVariable("obs", "u16", op="sum", absent="zero", leaf_max=1),
    DataVariable("occ", "u16", op="sum", leaf_max=1),
    DataVariable("hmax", "u16", op="max"),
]
data = {}
for k, leaf in enumerate(leaves):
    if k >= 90:  # nobody surveyed the last ten
        data[leaf] = {"obs": 0, "occ": None, "hmax": None}
    else:
        building = k < 25 or 50 <= k < 60
        data[leaf] = {"obs": 1, "occ": int(building), "hmax": 240 if building else 0}

blob = itacart.encode_data_tree(data, variables)
len(leaves), len(blob)
```

## A cell and its ancestors, in one walk

Each node of the tree holds its own values: at a leaf, what was observed; above
it, the aggregate of its daughters. Reading a cell walks down the tree from the
quadrant, one step per level, and every step lands on an ancestor.

```{code-cell} python
for cell, values in itacart.read_data_path(blob, leaves[3]):
    print(f"{cell:40s} {values}")
```

Every node above the ten-metre cell holds the same counts, because the map
covers nothing else: 90 cells observed, 35 of them occupied, and a highest
obstacle of 240 dm. The one-metre cell's family of twenty-five lies entirely
inside the building. The percentage is derived when it is read, from the two
integers, so it says what it is a share of:

```{code-cell} python
top = itacart.read_data_cell(blob, block)
round(100 * top["occ"] / top["obs"], 1), round(100 * top["occ"] / len(leaves), 1)
```

Of what was observed, 38.9% is occupied. Of the whole cell, 35% is *known* to
be, and the ten cells nobody surveyed are not free: they are unknown.

## Unknown is not zero

Each variable declares what a missing value means. `obs` declares `"zero"`:
a cell nobody surveyed counts no observation. `occ` and `hmax` declare
`"unknown"`, and an unknown value is a cleared bit beside the values, so every
`u16` stays available as a value:

```{code-cell} python
itacart.read_data_cell(blob, leaves[95])
```

Operators skip what is unknown, and a mean is weighted by a count declared
with it, so it is a mean over what was observed. The page on
[the wire format](data_tree_encoding_spec.md) gives the arithmetic and the
example of the denominator that this choice settles.

## A parent and its daughter, each with its own value

A *pyramid* variable is observed at the leaves and aggregated above them. A
*level* variable is observed wherever it is given, at any node, and never
aggregated. A geofence class observed on the ten-metre cell is information the
leaves do not carry:

```{code-cell} python
fence = DataVariable("fence", "u8", species="level")
fenced = itacart.encode_data_tree(
    data, variables + [fence], observed={block: {"fence": 2}, leaves[0]: {"fence": 1}}
)
itacart.read_data_cell(fenced, block)["fence"], itacart.read_data_cell(fenced, leaves[0])["fence"]
```

## Blocks

Free space and the inside of a building repeat the same values over many
cells. With `prune=True`, the encoder folds every complete family of leaves
that agree on **every** variable into its parent, which then stands for the
whole block. A family with two heights stays a family, however fully covered.

```{code-cell} python
pruned = itacart.encode_data_tree(data, variables, prune=True)
schema = itacart.data_tree_schema(pruned)
len(pruned), schema["nodes"], itacart.read_data_cell(pruned, leaves[3]) == data[leaves[3]]
```

Two families fold: the one inside the building and one that is entirely free,
leaving 50 one-metre leaves of the hundred. A cell inside a block reads the
block's values, with a count shared evenly among its cells, so every cell
reads exactly as it did. Folding by occupancy alone would also have merged
families whose cells are all occupied at different heights.

## Updating in place

Values have fixed widths at fixed places, so an update rewrites one slot per
level of the path. A new obstacle on a free cell changes every count above it
by one:

```{code-cell} python
writable = bytearray(blob)
itacart.update_data_tree(writable, leaves[30], {"occ": 1, "hmax": 120})
itacart.read_data_cell(writable, block)
```

Adding or removing a cell changes the structure, and is a new encoding.

## The tree is lexical

The nodes are the cells and every prefix of their spellings, as in the
ancestry of the index. On the domain border, the surplus children of a cell
hang under a prefix that names no cell, so a border cell's aggregate covers its
own stem only: 80.05% of the trapezoid `NE(2003/0000)` with the data at
resolution 2. Such a prefix is never returned, and the polar caps enter like
any other cell. Both are measured in the
[specification](data_tree_encoding_spec.md), with the byte layout, the sizes
of the two cases the format was designed against, and the cost of a read.
