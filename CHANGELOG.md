# Changelog

Format based on [Keep a Changelog](https://keepachangelog.com/1.1.0/);
versioning follows [SemVer](https://semver.org/).

## [1.0.0]

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

### Fixed

- The parent of a border cell is resolved on two axes. Child enumeration no
  longer assumes a cell keeps its own stem plus one to the east, which had
  made the hierarchy of absorbing cells disagree with the geometry.
- The geodesic walk across the pole no longer loses the polar row.
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
