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

# Design Decisions

Four decisions a user needs to know about, because each one costs something
visible and each one was taken deliberately.

## Direct tessellation, and the EAERS divergence

ITACaRT maps no base unit polyhedron to the surface. The initial tessellation
is direct on the ellipsoid, for absolute geodetic fidelity.

That single choice is the whole of the EAERS divergence. Requirement 22 is not
met because there is no polyhedron; requirement 20 is only partial because the
harmonized model's polyhedral interface has nothing to describe. The package's
own justification puts it plainly: this is **a property of the reference system
rather than a gap**.

The cost is real — a grid that cannot present the polyhedral interface cannot
be dropped into tooling that expects one. The gain is that a parcel's area is
computed on the figure the survey was referred to, without a spherical
approximation in between. For cadastral work, where area carries legal weight,
that is the trade worth making.

## Ancestry is lexical

A parent is a prefix; ancestry is decided by reading strings, not geometry.
That is what the standard asks for — topological relations computed from the
identifiers — and it is what makes the operations cheap and exact about what
they decide.

It also means that at the domain border, lexical ancestry and geometric
containment part company: a child of an absorbing cell reaches past its
parent's effective ring, by construction. This is declared as a limitation of
this version, with its measurements, at
[Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md).

The alternative would be to make ancestry geometric, which would make every
hierarchical query a spatial computation and would still not repair the
underlying geometry. The decision is to keep the operation honest about being
lexical and to state where that is not enough.

## The antemeridian is not crossed

Each quadrant's lattice runs outward from the prime meridian and stops. Keeping
uniform areas all the way around is impractical, so ITACaRT accepts a
discontinuity at one meridian rather than distorting every cell to force a
closure.

The meridian it accepts the discontinuity on is almost entirely ocean. The
inhabited land it crosses is handled by two named extension zones — Fiji and
Chukotka — and Antarctica is excluded on the grounds of limited cadastral use.

```{code-cell} python
import itacart
itacart.describe()["boundary_treatments"]
```

The cost lands on ocean and Antarctic cells, which have unequal areas, and on
GeoJSON strictness inside the two zones, where longitudes run past 180 degrees.
[Antimeridian](../boundaries/antimeridian.md) and
[Extension Zones](../boundaries/extension_zones.md).

## Coordinates are never rounded

Exported coordinates carry full precision, and the reason is arithmetic rather
than fastidiousness.

Six decimal places is roughly ten centimetres. That is smaller than a cell only
down to resolution 10; from resolution 11 the **whole cell is narrower than one
rounding step**, and any rounding collapses it into a degenerate polygon that is
still syntactically valid GeoJSON — valid, and meaningless.

Emitting full precision costs bytes and never lies. Given that the grid's
finest resolution is one centimetre and its cadastral purpose is measurement,
the bytes are the cheaper of the two.

## Areal invalidity is refused, not repaired

The thread running through all of the above. Where the package meets something
it cannot answer correctly, it raises instead of returning a plausible value:

```{code-cell} python
from itacart import geometry
from itacart.exceptions import ITACaRTError
from shapely.geometry import Polygon

tiny = Polygon([(-46.63290, -23.55090), (-46.63289, -23.55090),
                (-46.63289, -23.55089), (-46.63290, -23.55089)])
try:
    geometry.polyfill(tiny, resolution=5)
except ITACaRTError as exc:
    print(type(exc).__name__, "--", exc)
```

An empty cover silently returned would be read as "no overlap" and would reach
a land registry as a parcel with no extent. A nominal area returned for a
clipped cell would misstate a holding by half a cell. A nonexistent cell mapped
quietly to its eastern twin would bury an upstream arithmetic fault in stored
data.

In each case the plausible answer is the dangerous one, because it looks right.
The refusal is louder, and it fails at the moment the mistake is cheapest to
fix.
