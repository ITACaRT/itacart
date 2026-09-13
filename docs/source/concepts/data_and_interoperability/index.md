# Data & Interoperability

Getting data in and out, and the semantics that survive the trip.

Two kinds of export live here and they answer different questions. **Standard
formats** — GeoJSON, WKT, GeoDataFrame — hand the grid to tools that know
nothing about it. **Binary encodings** — TreeBlob and GeometryBlob — keep a
cell set or a geometry in a form that is compact, hashable and verifiable.

Running through both is one distinction worth reading before either: what an
export carries, and what it does not. A rendering of a cell is not the cell,
and the page that says so is the second one below.

```{toctree}
:maxdepth: 1

geojson
geographic_vs_grid_geometry
edge_semantics
serialization
geometry_blob
tree_blob
binary_encoding_spec
geometry_encoding_spec
```
