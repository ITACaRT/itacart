# API reference

Everything below is generated from the source at build time. The page for a
module is its module docstring plus one entry per public name, taken from the
docstring of that name — there is no second copy of the API to keep in sync.
The facade re-exports every public name, so `itacart.geo_to_cell` and
`itacart.cells.geo_to_cell` are the same function and either import works.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: autosummary/module.rst
   :recursive:

   itacart
```

## Reading order

The package is layered. Each layer below depends on the ones above it and on
nothing beneath, so reading downward never requires a forward reference. Every
module is paired with the conceptual page that argues what the module
implements: the docstring states the contract, the concept page states why the
contract is what it is.

### Foundations

{py:mod}`itacart.constants`
: Immutable values transcribed from the paper — the WGS84 ellipsoid, the
  resolution table, the refinement alphabets, the antemeridian extension zones.
  Data, not behaviour, with a single exception noted on the page.
  See [The ITACaRT Grid](../concepts/grid_fundamentals/the_itacart_grid.md).

{py:mod}`itacart.exceptions`
: Every error the package raises, rooted at `ITACaRTError`, so a caller can
  guard a whole pipeline with one `except` and never catch a stray
  `ValueError`. Which failures are refusals rather than bugs is the subject of
  [Cell Validity & Existence](../concepts/index_and_hierarchy/cell_validity_and_existence.md).

{py:mod}`itacart.geodesy`
: The ellipsoidal sinusoidal projection and the geodesic solutions, computed
  from the paper's equations. PROJ is not called at run time.
  See [The Sinusoidal Construction](../concepts/grid_fundamentals/the_sinusoidal_construction.md).

{py:mod}`itacart.resolutions`
: The resolution table and the refinement ratios, as helpers rather than as raw
  data. See [Resolution & Scale](../concepts/grid_fundamentals/resolution_and_scale.md).

### Addressing

{py:mod}`itacart.index`
: Parsing, composition and canonical form of the compositional index — the
  bridge between one string and the cells it addresses.
  See [Index Structure](../concepts/index_and_hierarchy/index_structure.md) and
  [Composition & Decomposition](../concepts/index_and_hierarchy/composition_and_decomposition.md).

{py:mod}`itacart.hierarchy`
: Ascent, descent and compaction across resolutions. Note that the functions
  here divide into those that answer one value per cell and those that answer
  many, and that the two families collapse differently for a single cell; each
  docstring states which it is. See
  [Parent, Children & Descendants](../concepts/index_and_hierarchy/parent_children_and_descendants.md),
  [Compaction & Normalization](../concepts/index_and_hierarchy/compaction_and_normalization.md)
  and [Lexical Ancestry vs Physical Refinement](../concepts/index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md).

{py:mod}`itacart.cells`
: Quantization and inverse geometry — the core mapping between an address and a
  position. See [Point-to-Cell Indexing](../concepts/spatial_operations/point_to_cell_indexing.md)
  and [Cell Geometry](../concepts/geometry/cell_geometry.md).

### Working against geography

{py:mod}`itacart.geometry`
: Vector geometry against the grid: filling a polygon, mapping vertices,
  canonical form. See [Cell Filling / Polyfill](../concepts/spatial_operations/cell_filling_polyfill.md),
  [Containment Modes](../concepts/spatial_operations/containment_modes.md) and
  [Canonical Geometry](../concepts/cadastral_workflows/canonical_geometry.md).

{py:mod}`itacart.topology`
: Neighbourhood and adjacency on the parallelogram lattice.
  See [Neighbors & Grid Traversal](../concepts/spatial_operations/neighbors_and_grid_traversal.md)
  and [Grid Distance](../concepts/spatial_operations/grid_distance.md).

{py:mod}`itacart.boundary`
: Prime meridian, antemeridian and extension zones — where the uniform grid
  stops being uniform. The whole of
  [Boundaries & Global Continuity](../concepts/boundaries/index.md) is the
  argument behind this module.

### Measuring

{py:mod}`itacart.metrics`
: Cell shape metrics, for comparison against other DGGS implementations.
  Compactness and normalized area follow the conventions of Kmoch, Vasilyev,
  Virro and Uuemaa (2022) so that the numbers are comparable with that paper's
  figures rather than merely accurate; the cell base angle is the ellipsoidal
  image of the 45 degree angle the ITACaRT paper draws on the sinusoidal plane.
  Each function's page states which convention it transcribes and where the
  convention stops being trustworthy. See
  [Area, Perimeter & Compactness](../concepts/geometry/area_perimeter_and_compactness.md)
  and [Base Angle & Shape](../concepts/geometry/base_angle_and_shape.md).

### Moving data

{py:mod}`itacart.serialization`
: The binary encodings of indices and geometry.
  See [Serialization](../concepts/data_and_interoperability/serialization.md),
  [GeometryBlob](../concepts/data_and_interoperability/geometry_blob.md) and
  [TreeBlob](../concepts/data_and_interoperability/tree_blob.md).

{py:mod}`itacart.interop`
: Export to standard geospatial formats. What an export carries and what it
  necessarily drops is the subject of
  [Geographic vs Grid Geometry](../concepts/data_and_interoperability/geographic_vs_grid_geometry.md);
  the format itself is covered in [GeoJSON](../concepts/data_and_interoperability/geojson.md).

### Standards

{py:mod}`itacart.engine`
: The OGC-facing engine object, which is the shape the standard expects a
  reference system to present. What it conforms to, and where it does not, is
  in [Conformance Status](../concepts/standards_and_conformance/conformance_status.md).
