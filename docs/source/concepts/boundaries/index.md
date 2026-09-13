# Boundaries & Global Continuity

What happens where the grid meets the seams of the globe.

A lattice is regular and a planet is not. ITACaRT resolves the mismatch in four
places, and in all four the resolution is written into the specification rather
than left to the implementation: the prime meridian and the equator, which
divide the quadrants; the antemeridian, which the grid does not cross except in
two named zones; the outer edge of each quadrant, where cells absorb what is
left over; and the poles, where the lattice thins to almost nothing.

These are the pages where the guarantees elsewhere in this documentation are
qualified. Read them before trusting an area computation near a seam.

```{toctree}
:maxdepth: 1

meridian_and_equator
antimeridian
extension_zones
absorbing_boundary_cells
polar_caps
```
