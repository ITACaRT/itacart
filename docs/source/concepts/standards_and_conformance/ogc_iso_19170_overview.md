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

# OGC / ISO 19170 Overview

What the standard asks of a discrete global grid reference system, in the terms
the standard itself uses.

```{code-cell} python
import itacart
report = itacart.conformance()
report["standard"], list(report["classes"])
```

OGC 20-040r3, published as ISO 19170-1 — Topic 21, *Discrete Global Grid
Systems*. Two conformance classes concern a system like this one: **Core**, the
baseline every DGGS must meet, and **EAERS**, the Equal Area Earth Reference
System class.

## What Core asks for

Core is about being a reference system at all. Its requirements cover the
things a client has to be able to rely on: a harmonized data model, a defined
coordinate reference system, simple cell geometry, a designated direct position
for each cell, globally unique addresses, quantization functions that take
positions to cells, and topological query functions that give parent, child and
neighbour relationships.

The last of those is worth drawing out, because it shapes the whole index
design. The standard expects these relations to be **computed from the cell
identifiers**, not from geometry — which is why ancestry in ITACaRT is a
property of the index string and neighbours are arithmetic on codes.
[Index Structure](../index_and_hierarchy/index_structure.md) is that design;
[Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md)
is where it is qualified.

## What EAERS asks for

EAERS is stricter and narrower: it asks that the system be built from a base
unit polyhedron mapped to the surface, that it present the polyhedral
interface of the harmonized model, and that cell areas be equal within a stated
budget.

Two of those three are where ITACaRT diverges, and it diverges by
construction. The grid tessellates the ellipsoid directly, so there is no
polyhedron and no polyhedral interface to describe.
[Design Decisions](design_decisions.md) is the argument for that choice.

## How the system describes itself

The standard expects a reference system to be able to announce what it is, and
the engine does:

```{code-cell} python
description = itacart.describe()
{k: description[k] for k in ("identifier", "title", "standard")}
```

```{code-cell} python
description["tessellation"], description["constraints"]
```

Note the `constraints` entry: `cellEqualSized` is asserted **with its
exceptions named in the same object**. A client reading this description learns
both that the grid is equal-area and exactly which families are not, without
having to read a paper.

```{code-cell} python
description["crs"]
```

## Reading a conformance claim

Three words in the standard do different work, and the package's report keeps
them apart: `conformant` is the formal yes-or-no against a *shall*; `coverage`
is `full`, `partial` or `none` and describes how much of the requirement is
substantively met; `justification` says why.

A requirement can be substantively almost entirely met and still formally
non-conformant, because a *shall* does not accept "nearly".
[Conformance Status](conformance_status.md) is where that distinction does the
most work.
