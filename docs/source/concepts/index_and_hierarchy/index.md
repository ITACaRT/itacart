# Index & Hierarchy

How a cell is named, and how names relate to one another across resolutions.

The index is the part of ITACaRT a user handles most, and it is doing more work
than an identifier usually does. It is human-readable, it carries the cell's
whole ancestry in its own text, it addresses a set of cells as easily as one,
and parent and neighbour relations are computed from the string rather than
from coordinates.

That density is deliberate — a discrete global grid is characterised by its
indexing method more than by anything else — but it means the notation has to
be read carefully. These pages read it, in order, and end with the two
distinctions that catch people: what ancestry in the string does and does not
claim about geometry, and the difference between an index that is legal and a
cell that exists.

```{toctree}
:maxdepth: 1

index_structure
composition_and_decomposition
parent_children_and_descendants
lexical_ancestry_vs_physical_refinement
compaction_and_normalization
cell_validity_and_existence
```
