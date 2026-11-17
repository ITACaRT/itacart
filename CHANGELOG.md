# Changelog

Format based on [Keep a Changelog](https://keepachangelog.com/1.1.0/);
versioning follows [SemVer](https://semver.org/).

## [1.2.0] - 2026-11-17

The count asks a whole level of its descent at a time. The counts are the ones
release 1.1.0 returns; what changes is how long they take, how workers start,
and a dependency the count now imports. With them come corrections to the
documentation and the declaration of a limitation at the poles.

### Changed

- `count_internal_cells` evaluates its descent a level at a time, in arrays.
  The nodes, the squares and the predicates are those of release 1.1.0: every
  square is built with the same floating-point arithmetic in the same order
  and asked the same prepared predicates, with one Shapely call for up to
  50 000 squares instead of one call per square. A leaf is decided by its
  centre, which is the answer the three questions of 1.1.0 give on a valid
  geometry; where the polygon's piece on the lattice is not a valid geometry,
  the leaf is asked those three questions. The counts are identical. They are
  pinned against twenty counts taken with 1.1.0 and against the recursive
  descent of 1.1.0, which stays in the package as the reference: on the public
  cases, and on 600 random polygons in the plane of the lattice, a number the
  release check raised to 15 000 with no divergence. Serially, on one Linux
  machine with two cores, Python 3.12, Shapely 2.2 and GEOS 3.14, a square
  kilometre at resolution 13 went from 34 s to 0.74 s, and a one-degree box at
  resolution 9 from 43 s to 1.25 s. The gain is the ordinary cells'. The walks
  of the last lattice column of a row, the polar row, the caps and the
  prime-meridian column still descend one node at a time, and a patch of the
  polar row that took 34 s at resolution 11 in 1.1.0 still does.
- **`count_internal_cells(..., n_jobs > 1)` spreads a large descent over
  processes instead of threads, started by `spawn` on every system, Linux
  included. A script whose count reaches the processes now needs an
  `if __name__ == "__main__":` guard around its entry point, on every system:
  without it each worker runs the script again as it starts, and the count
  fails with `BrokenProcessPool`.** A notebook needs no guard, and a count too
  small to repay starting the workers stays in the calling process and needs
  none either. `spawn`, and not the platform's default, because a process
  forked from a parent with another thread alive -- a notebook kernel, a
  server -- can deadlock, and Python warns that it may. In 1.1.0 two threads
  made the count 4.1 times slower than one, for the one-degree box at
  resolution 7; on two processes the box at resolution 9 takes 0.73 s against
  1.25 s on one. The count does not depend on `n_jobs`.
- `polyfill` runs in the calling thread whatever `n_jobs` says. In 1.1.0 it
  spread its descent over threads, and two threads made it 3.9 times slower
  than one, for the one-degree box at resolution 5. The parameter is accepted
  and validated as before, and the fill is the same.
- `numpy` is a declared runtime dependency, from 1.14, the floor Shapely 2.0
  already requires, so installing the package installs nothing new. The
  batched count imports it.
- The package summary reads "Mapping — Designed" where it read
  "Mapping-Designed".

### Removed

- The `parallel` extra. It installed `joblib`, which no module of the package
  ever imported. Asked for `itacart[parallel]`, pip still installs the
  package, and warns that it provides no such extra.

### Fixed

- The documentation promised what the package did not do. The installation
  page and the README said the `parallel` extra filled cells across several
  processes, and the page on cell filling that spreading the fill over workers
  paid off for large parcels. In 1.1.0 the fill used threads with or without
  the extra, and was slower on two of them. The pages now say what `n_jobs`
  does.
- The page on the TreeBlob said a blob holds one quadrant. It holds one group
  per quadrant the index names, and the page now encodes an index over three
  quadrants and reads it back whole.
- The docstring of `cell_to_boundary` said a trapezoid could keep three or
  five vertices. A ring has three or four: the trapezoid's outer side is
  carried onto the border rather than cut by it.

### Known limitations

Stated in the modules that own it, `hierarchy` and `cells`, and pinned by named
tests that fail when it is lifted.

- **The cell that holds a pole at resolution 10 is not refined.**
  `get_children` raises `GeometryError` for it as its result is iterated, and
  so do `cell_to_children`, which is the same function, and `get_descendants`,
  `uncompact_cells` and `cells_to_intervals` that go from resolution 10 to 11
  there; `count_internal_cells`, `polyfill`, `from_geojson` and
  `from_geodataframe` raise it for a geometry that reaches the cell, from
  resolution 11. Inside the cell, at resolutions 11 and 12, `geo_to_cell` can
  answer with a spelling that names no cell, and so can `latlng_to_cell`,
  `sinusoidal_to_cell`, `vertex_to_cell`, and `from_geojson` and
  `from_geodataframe` for a point: `is_valid_cell` denies the spelling, and
  the functions that take a cell then refuse it. That is an exception to the
  rule that no function answers with a cell that does not exist. Four
  positions are known to be answered so, at longitude ±155.0390625 and
  latitude ±89.99999962902984.

  Releases 1.0.0 and 1.1.0 behaved the same, in every name they had.

## [1.1.0] - 2026-10-30

Three features: an integer for every cell, a region as ranges of those
integers, and a format for data at every level of the tree. A measure joins
them: how far a position written as a cell can lie from the cell's centre. With
them come the corrections release 1.0.0 needed.

### Added

Numbering

- `integer_index`: one 64-bit integer per cell, in the manner of the NUNIQ
  scheme of the IVOA's coverage maps. `cell_to_int`, `int_to_cell`,
  `int_resolution` and `int_parent`, with `str_to_int` and `int_to_str` as
  aliases under the names H3 gives them. The integer carries the resolution,
  the four quadrants of resolution 0 included, sorts in canonical order within
  a resolution and reaches the lexical parent by arithmetic. The largest is
  8.4335·10¹⁸, below 2⁶³, so it fits a signed 64-bit column. The existence
  contract holds both ways: a spelling that names no cell gets no integer, and
  an integer that names no cell gets no spelling. A NumPy integer is accepted;
  a float and a bool are refused.
- `intervals`: a region as sorted, half-open ranges of those integers at a
  chosen depth, with `cells_to_intervals`, `intervals_to_cells`,
  `intervals_union`, `intervals_intersection`, `intervals_difference`,
  `intervals_contain`, `intervals_overlap` and `intervals_count`. The semantics
  is the region's, so a cell that absorbs the border is expanded through its
  physical children. For an index whose cells do not overlap, the way back is
  what `compact_cells` writes, and its TreeBlob is the region's, byte for byte.

Data on cells

- `serialization.data_tree_blob`: the DataTreeBlob, values on a lexical tree
  of cells, read and updated in place. `encode_data_tree`, `decode_data_tree`,
  `read_data_path`, `read_data_cell`, `slice_data_tree`, `update_data_tree`,
  `data_tree_schema` and the `DataVariable` declaration.
  - Up to 255 fixed-width variables per blob. A *pyramid* variable is observed
    at the leaves and, unless its operator is none, held at every internal
    node as the sum, mean, maximum, minimum or mode of the daughters, computed
    exactly and rounded once. A *level* variable is observed at any node, so a
    parent and its daughter keep values of their own.
  - A mean or a mode is weighted by a declared count, so it is over what was
    observed and not over what could have been.
  - Each variable declares whether a missing value is zero or unknown. An
    unknown value is a cleared bit in a validity bitmap, so every type keeps
    its whole domain.
  - A leaf above the observation resolution stands for its block, and
    `prune=True` folds complete families of leaves that agree on every
    variable, whose parent is a cell that does not absorb the border, at
    levels that carry no level variable.
  - A cell and its ancestors are read in one walk down the tree, on the bytes.
    An update rewrites one slot per level of the path; on a blob encoded
    without pruning, it leaves the bytes that a new encoding of the new data
    would write.
  - The polar caps are in. On the border the tree is lexical, as the index's
    ancestry is, and the aggregate of a border cell covers its own stem; a
    prefix that names no cell is a node of the blob and is never returned.

Metrics

- `cell_max_uncertainty`: how far a position written as a cell can lie from
  the cell's centroid, measured on the WGS84 ellipsoid with `inverse_geodesic`.
  `output="radius"`, the default, is the distance from the centroid to the
  farthest vertex, which is the reach `cells_to_geometry` states for a rebuilt
  vertex; `"diameter"` is the largest distance between two vertices; `"ratio"`
  is the radius over the side of the resolution. Every cell that exists is
  answered, at a cost of six geodesics at most. For a regular parallelogram the
  radius follows `l · max(√(1 + u²)/2, √((1 − u/2)² + 1/4))`, with
  `u = |λ sin φ|`: from `l/√2` where the cell is a square on the ground to
  `1.65 l` at the antimeridian, and a little more past it in Chukotka, so cells
  of one resolution differ by a factor of more than two. The departures from
  that law are measured and stated, by band of latitude. The maxima are taken
  at the vertices, not over the boundary: in the cell that holds a pole, whose
  base is a whole parallel, the diameter is half the boundary's, and under a
  cap a spelling outside the cap's physical refinement can reach further along
  its boundary than at its vertices.

Contracts and instruments

- The half of the existence contract that says "a cell that exists is
  accepted" is measured: every public name that takes a cell, except
  `to_geodataframe`, which needs the optional `geo` extra, is handed every
  valid cell of the corpus, and each refusal is declared, with its exception
  and its reason, in a set pinned by equality, with a control that sees an
  undeclared refusal reported. The readers of a blob are not handed the polar
  caps, which no encoder writes.
- The bytes of release 1.0.0 are pinned by hash over an enumerated corpus, for
  `encode_tree`, `serialize_to_blob`, `encode_geometry` and
  `geometry_to_tree`.
- The installed metadata and `pyproject.toml` are tested against
  `__version__`.
- `cell_base_angle` and `compactness` are tested from resolution 1 to 13, on
  chains in the four quadrants, the meridian triangles and the polar caps, on
  an enumeration of the border slivers, on every cell of the polar caps down
  to resolution 5 and on a sample around the poles below it: the angle
  against the mirrored twin, the closed-form limit and an independent
  integration of the azimuth, and compactness for the cells it refuses and
  for its authalic latitude against a 60-digit reference.

Documentation

- Concept pages for the integer index, the intervals and the DataTreeBlob, and
  the DataTreeBlob's byte-level specification; a section on the positional
  uncertainty of a cell in the page on area, perimeter and compactness.

### Changed

- `inverse_geodesic`: the convergence of the Vincenty iteration is judged
  relative to the line's angular length, up to one radian, with a floor of
  four units in the last place, and the `tol` argument is read the same way.
  A short line now converges to the right azimuth. Over 2 000 random lines per
  length, against GeographicLib, the largest azimuth error falls from
  3.5·10⁻² degree to 1.1·10⁻⁵ at 1 cm, and from 3.0·10⁻⁴ to 1.5·10⁻⁷ at 1 m.
  Results change for lines shorter than a radian, about 6 400 km, and less
  the longer the line: the azimuth by the errors above at a centimetre or a
  metre, by up to 3.5·10⁻⁷ degree at 1 km and by less than a nanodegree from
  1 000 km, the distance by a few micrometres. Lines longer than a radian, the
  near-antipodal ones included, are solved exactly as before, and so is ITA to
  Curitiba.
- `densify_segment` and `densify_orthodromic` inherit that change, and so do
  `polyfill` and `count_internal_cells`, which densify a polygon before
  filling it. Interior points move by micrometres, and an edge whose length
  falls within micrometres of a multiple of the segment length can be cut into
  one piece more or fewer, so a filled index or count, and the TreeBlob of the
  index, can differ from 1.0.0's at fine resolutions.
- `cell_base_angle` and `compactness` hold from resolution 1 to 13. Release
  1.0.0 compared vertices against a micro-degree, about 11 cm, more than a
  cell measures from resolution 11 down: the angle came out 0 for every such
  cell, compactness refused most of them, and thin border trapezoids were
  refused from resolution 3. The angle's frame read the quadrant's letter,
  which answered 0 for half the cells under every meridian triangle and polar
  cap; it now follows the cell's chirality. The angle is 0 for a triangle and
  for nothing else. Compactness refuses the cells whose ring wraps the whole
  parallel, and only those, all of them inside the polar caps: the one that
  holds the pole, at every resolution, and from resolution 4 down many cells
  around it whose rings wrap the parallel once or many times over. Down to
  resolution 10, release 1.0.0 returned a number for those that was not their
  quotient. The authalic latitude is computed without cancellation near the
  pole.
- Their values move as well. Away from the cells that used to answer 0, the
  angle, whose short geodesic now converges, moves by up to 3·10⁻⁸ degree at
  resolution 1 and by up to about 6·10⁻⁴ at resolution 10: the São Paulo cell
  of the getting-started, at resolution 9, measured 55.988883 degrees and
  measures 55.988672. Compactness moves in its last digits for almost every
  cell, by up to about 10⁻¹¹ relatively at resolution 1 and a few 10⁻⁷ at
  resolution 10.
- TreeBlob: in `encode_tree`, `recompose_to_prefix_form`, `encode_node` and
  `serialize_to_blob(compact=False)`, an index that writes a cell and also a
  descendant under the cell's own spelling now encodes the cell, absorbing the
  descendant as `compose` and `normalize` do. Release 1.0.0 kept the
  descendant and lost the rest of the cell. The bytes change for that input
  only; every other input encodes as it did.
- The binary codecs still do not serialize the polar caps in row 1000, and
  their refusal now says that this is a restriction of the codecs, not of the
  grid, in the TreeBlob and the GeometryBlob alike. Rows 1001 and beyond are
  still refused as outside the grid.
- Continuous integration runs the current major versions of its actions, and
  no longer triggers on a `develop` branch that does not exist. The pre-commit
  hooks build their environments with Python 3.12 and run at the commit stage
  unless a hook names another.

### Fixed

- The TreeBlob specification states its density by the parity of the leaf
  resolution, the cost of scattered leaves by resolution, row 1000 as a
  restriction of the codecs, the information a node carries and the fields a
  node blob writes, and the number of behaviours 1.0.0 corrected; every figure
  it quotes is pinned by a named test.
- Four pages said compaction was lexical; it is physical, and they now show
  the border cell `NE(2003/0000)` it folds.
- The documentation of `cells_to_geometry` measured the reach's departure from
  its law, for parallelograms and for the meridian triangles, on a lattice
  that skipped the last rows before the polar row, where the departure is
  largest, and quoted it as a bound. It now quotes departures measured down to
  those rows, by band of latitude.

## [1.0.0] - 2026-10-22

The first functional release. Where `0.1.0a4` defined the public surface in
stubs, this release implements it, holds every published figure to a named
test, and documents it.

### Added

Grid and geodesy

- Ellipsoidal sinusoidal projection and geodesic solutions implemented from
  the paper. PROJ is absent from the runtime and from every extra; the
  solutions are held to a Runge-Kutta integration of the geodesic equations
  and to Clairaut's relation rather than to a second library.
- `resolutions` and `cells`: fourteen levels from 10 km to 1 cm, alternating
  one-to-four and one-to-twenty-five refinement, with whole-number metric
  areas.
- `boundary`: the prime meridian column, the antemeridian trapezoids, the
  polar row and the extension zones, each enumerated rather than sampled.
- `nominal_cell_area` and `effective_cell_area` kept distinct, and
  `cell_to_anchor` (requirement 12) and `cell_to_centroid` kept distinct.

Navigation and geometry

- `hierarchy`: ascent, descent and compaction, with lexical ancestry and
  physical refinement answered by separate functions and the difference
  between them documented in the module.
- `topology`: neighbourhood, with border deflection.
- `geometry`: `polyfill`, densification and canonicalization. Filling
  descends the absorbing and polar families instead of refusing them, so
  the last column and the polar row are covered.
- `metrics`: `compactness`, `cell_base_angle` and `normalized_cell_area`,
  validated against regular hexagon, square and triangle before being
  applied to a cell.

Encoding and interoperability

- `serialization`: TreeBlob and GeometryBlob. Two inputs describing the same
  leaf set produce byte-identical blobs, so a content hash over those bytes
  is a stable key for the region.
- `interop`: GeoJSON and WKT export, GeoDataFrame conversion, and recovery
  of an index from a Feature. Coordinates are never rounded: from resolution
  11 a cell is narrower than one six-decimal rounding step.
- H3 v4 compatible aliases.

Engine, conformance and contracts

- `engine`: `describe`, `crs` and `conformance`.
- OGC DGGS Core / EAERS conformance suite, written against OGC 20-040r3 —
  the edition ISO 19170-1 follows — with every previously expected failure
  closed.
- An existence-and-spelling contract: a public entry point handed an index
  that is well formed but names no cell refuses it rather than computing on
  it.
- A cell answers as a point at its centroid; `from_geojson` quantizes Point
  across the three containment modes; `recover_from_geojson` refuses rather
  than returning an empty list; `LineString` is refused for having no area.

Documentation

- A published documentation site covering getting started, concepts and
  guides, and the API reference. Its examples are versioned MyST and execute
  during the build, so a page that claims a result has produced it.

### Changed

- Licence migrated from MIT to BSD 3-Clause. Both are permissive and neither
  restricts use; the third clause of BSD 3-Clause is the difference, and it
  withholds the right to use the copyright holder's name to endorse derived
  products without written permission.
- `normalize` no longer rewrites a spelling in the one way that would change
  the region it names.
- Cell validity now rejects a ring that is not simple in effective geometry,
  not only one that is malformed as a spelling.
- A polygon or multipolygon that is not a valid geometry is refused before
  anything reads it, with the reason the geometry engine gives, and is never
  repaired: a zero-width buffer or `make_valid` would decide which area was
  meant.

### Fixed

- The parent of a border cell is resolved on two axes. Child enumeration no
  longer assumes a cell keeps its own stem plus one to the east, which had
  made the hierarchy of absorbing cells disagree with the geometry.
- The geodesic walk across the pole no longer loses the polar row.
- `count_internal_cells` counts the last column, the polar row and the caps
  in full. It walked them with the fill's ceiling as its budget and returned
  the running total at the moment of passing it as though it were the count.
- `build-system` required `setuptools>=61`, incompatible with the SPDX form
  of `license` (PEP 639). Raised to `>=77`.
- `LICENSE` was missing, breaking `license-files` in `pyproject.toml`.

### Known limitations

Both are stated in the module that owns them, and the figures in each are
pinned by named tests.

- **An exported GeoJSON polygon is not the cell.** RFC 7946 defines the line
  between two positions as a straight segment in longitude and latitude,
  while a cell is defined by straight edges in the sinusoidal plane. Export
  and recovery are an inverse pair because the index travels with the
  Feature; export and refill are not, and nothing here promises that
  reingesting a rendering returns the cell.
- **Lexical ancestry and physical refinement part company at the domain
  border.** A cell whose outer side has been carried onto the border reaches
  east past its own nominal column, so some of its children are spelled
  under a prefix that names no cell. Those children exist and are addressed
  uniquely; their physical parent is found by refinement, not by slicing the
  string.

## [0.1.0a4] - 2026-08

### Changed
- `requires-python` raised from 3.8 to 3.10.
- Package layout defined under `src/itacart/`.

## [0.1.0a3]

- Placeholder release on PyPI.
