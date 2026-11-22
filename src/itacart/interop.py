"""Export to standard geospatial formats.

Satisfies OGC DGGS Core requirements 18 and 19 (interoperability
functions). The paper marks these as met by design; concrete exporters are
what turn that into a demonstrable claim.

Provenance: new; the notebook did this ad hoc with GeoPandas.

Three edge conventions meet in this module and they are not the same. A cell
is defined by straight edges in the *sinusoidal plane*, and that is the
convention its area is computed under. Geodesic paths are a second
convention, used wherever a distance on the ellipsoid is meant. GeoJSON
imposes a third: RFC 7946 section 3.1.1 defines the line between two
positions as a straight segment in longitude and latitude, and warns that it
may differ markedly from the path along the curved surface. **The exported
polygon is therefore not the cell.** The difference is a property of the
format, not a defect to repair, and it is measured in the test suite rather
than assumed.

Coordinates are never rounded. Six decimal places are roughly ten
centimetres, which is smaller than a cell only down to resolution 10; from
resolution 11 the whole cell is narrower than one rounding step and any
rounding collapses it into a degenerate polygon that is still syntactically
valid GeoJSON. Emitting full precision costs bytes and never lies.

**Export and recovery are an inverse pair. Export and refill are not.**
:func:`cells_to_geojson` writes the index into the Feature ``id`` member and
the cell's vertices into the geometry; :func:`recover_from_geojson` reads the
index back out of the Feature. The pair returns the cell it was handed, and
it does so because the index travelled with the Feature rather than being
inferred from the coordinates. That is the scope of the claim and it is
stated rather than implied: a recovery handed a Feature whose geometry
belongs to another cell still answers the original index, and a named test
does exactly that. The geometry is exact in its own right, and separately --
the exported vertices, joined again as straight segments in the sinusoidal
plane, rebuild the cell to within a unit in the last place.

Taking the exported polygon for a new geographic geometry -- joining its
positions along the ellipsoid and filling the result -- is neither half of
that pair. It is reingestion of a rendering, and nothing here promises it
returns the cell.

**What reingestion costs, measured.** Row 999 at resolution one holds exactly
six cells: four in the last column, one per quadrant, each absorbing the
border, and two prime-meridian triangles, which exist in the eastern
quadrants only because the western quadrants have no meridian column. Six is
the enumerated total of that row, not a sample of it. With orthodromic
densification at one kilometre the reingested outline covers 33.1111 per cent
of an absorber and 51.6779 per cent of a triangle. **Both figures are
readings of that step rather than properties of the grid**: refined, they
settle near 33.03 and 51.57 per cent, and the declared readings and the
settled ones are pinned together, so neither can be read as the other. Handed
such an outline, :func:`~itacart.polyfill` under ``center`` never answers the
cell. For four of the six it raises rather than answering an empty cover,
because the figure that comes back holds no cell centre. For the two absorbers
of the eastern quadrants it answers the cap of their own hemisphere: every
absorber of the row comes back with an edge along the meridian, the cap's
centroid lies exactly on that edge, and the centre rule moves a centre on an
outline a micrometre east before testing it, which puts that centroid inside
the two eastern figures and outside the two western ones. Under
``intersects`` the cell is in the answer along with others the outline
touches; :func:`recover_from_geojson` returns the cell. An ordinary cell loses
less than a thousandth of its figure on the same route, so the loss belongs to
this row and not to the route.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any, cast

from shapely.errors import ShapelyError
from shapely.geometry import MultiPolygon, Polygon, box, mapping, shape
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

from .boundary import cell_shape, extension_zone, is_valid_cell
from .cells import cell_to_boundary, geo_to_cell
from .constants import WGS84_A, WGS84_E2
from .exceptions import (
    GeometryError,
    NonExistentCellError,
    UnsupportedGeometryTypeError,
)
from .geometry import Containment, polyfill
from .index import compose, decompose
from .resolutions import effective_cell_area, get_resolution, nominal_cell_area

if TYPE_CHECKING:
    import geopandas

__all__ = [
    "cells_to_geojson",
    "cell_to_wkt",
    "cells_to_wkt",
    "to_geodataframe",
    "from_geodataframe",
    "from_geojson",
    "recover_from_geojson",
]

#: Default target length of a densified segment along the base parallel of a
#: polar cap, in metres. Matches the package-wide densification default.
_DEFAULT_MAX_SEGMENT_M = 1000.0

#: Property names written by :func:`cells_to_geojson` when
#: ``include_metadata`` is true. These are public API: once emitted, renaming
#: one breaks every consumer that reads it.
INDEX_PROPERTY = "itacart_index"
RESOLUTION_PROPERTY = "itacart_resolution"
SHAPE_PROPERTY = "itacart_shape"
NOMINAL_AREA_PROPERTY = "itacart_nominal_area_m2"
EFFECTIVE_AREA_PROPERTY = "itacart_effective_area_m2"
EXTENSION_ZONE_PROPERTY = "itacart_extension_zone"


def _require_cell(cell: str) -> None:
    """Refuse an index that names no cell.

    Six public functions return geometry for a cell that does not exist. The
    exporters do not inherit that: a caller who asks for GeoJSON of an
    unoccupied address gets a refusal rather than a plausible polygon, and
    the refusal is a package exception so that one ``except ITACaRTError``
    still guards the pipeline.
    """
    if not is_valid_cell(cell):
        raise NonExistentCellError(f"index names no cell: {cell}")


def _holds_pole(ring: list[tuple[float, float]]) -> bool:
    """Whether a boundary ring reaches the pole.

    This is the operational test for a polar cap, and it is deliberately not
    a longitude-jump test. Minus 180 and plus 180 are the same physical
    meridian, and at the pole longitude is undefined, so a detector built on
    the size of a longitude step cannot tell a cap from a feature that
    genuinely occupies both sides of the seam.

    Nor is it a proximity test, and that distinction is not cosmetic. This
    was a comparison against 89.999 degrees, justified by the measurement
    that no cell which is not a cap reaches beyond 89.9824. That
    measurement was true and its scope was resolution one, which was never
    said. Refined, ordinary cells walk toward the pole without limit: three
    of the cap's grandchildren pass 89.999 at resolution four and three
    more at resolution five, none of them caps, and each was handed to the
    cap routine and exported as a band around the whole parallel. The
    inflation ran from ten to seventeen times the cell's own area, and
    three distinct cells came out as one identical footprint. No fixed
    threshold below ninety survives refinement, because the cells keep
    coming.

    The pole is not approximated anywhere in this package, so it does not
    have to be approximated here. ``itacart.boundary.to_geodetic`` answers
    ``copysign(90.0, y)`` as a literal whenever a plane vertex reaches the
    meridian quadrant, which is exactly the construction that produces a
    cap. Asking for that literal asks the same question the grid answers,
    rather than a question that resembles it.
    """
    return any(abs(lat) >= 90.0 for _, lat in ring)


def _parallel_circumference_m(latitude: float) -> float:
    """Length of the full parallel at ``latitude`` on the WGS84 ellipsoid."""
    phi = math.radians(latitude)
    prime_vertical = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(phi) ** 2)
    return 2.0 * math.pi * prime_vertical * math.cos(phi)


def _polar_cap_polygon(
    ring: list[tuple[float, float]],
    max_segment_m: float = _DEFAULT_MAX_SEGMENT_M,
) -> Polygon:
    """Build the interoperable polygon of a polar cap.

    The three-vertex ring that :func:`itacart.cell_to_boundary` returns for a
    cap is a coordinate-domain figure, not a triangle on the ellipsoid. Its
    two meridional sides at minus 180 and plus 180 are two copies of the same
    physical seam, its base is a complete parallel, and its apex is a point
    at which every longitude names the same place. A GeoJSON consumer has no
    way to know that, so the exporter converts the cap into a figure that
    means the same thing under RFC 7946 rules: the base parallel densified
    across the whole longitude range, the surface opened along the seam, and
    the apex carried as a degenerate edge at the pole.

    The result is one polygon with one part. A cap does not occupy two sides
    of the seam; it is produced by the seam, so cutting it would invent a
    division that the surface does not have.
    """
    latitudes = [lat for _, lat in ring]
    northern = max(latitudes, key=abs) > 0.0
    base = min(latitudes) if northern else max(latitudes)
    circumference = _parallel_circumference_m(base)
    steps = max(4, math.ceil(circumference / max_segment_m))
    lons = [-180.0 + 360.0 * step / steps for step in range(steps + 1)]
    if northern:
        parallel = [(lon, base) for lon in lons]
        polar_edge = [(180.0, 90.0), (-180.0, 90.0)]
    else:
        parallel = [(lon, -90.0) for lon in lons]
        polar_edge = [(180.0, base), (-180.0, base)]
    positions = parallel + polar_edge
    return Polygon(positions + [positions[0]])


def _unwrap_longitudes(
    ring: list[tuple[float, float]],
) -> list[tuple[float, float]]:
    """Make a ring continuous in longitude by lifting it off the seam."""
    unwrapped = [ring[0]]
    for lon, lat in ring[1:]:
        previous = unwrapped[-1][0]
        while lon - previous > 180.0:
            lon -= 360.0
        while previous - lon > 180.0:
            lon += 360.0
        unwrapped.append((lon, lat))
    return unwrapped


def _cut_at_antimeridian(geometry: BaseGeometry) -> BaseGeometry:
    """Cut an arbitrary geometry that occupies both sides of the seam.

    RFC 7946 section 3.1.9 says geometry crossing the antimeridian should be
    split so that no part crosses it. This is the general rule and it applies
    to arbitrary input geometry on the way out. It is not what handles a
    polar cap: see :func:`_polar_cap_polygon` for why those are a different
    problem with a different answer.

    No ITACaRT cell reaches this path. Measured across every column of the
    equatorial, mid-latitude and polar rows in all four quadrants, the widest
    ring spans exactly 180 degrees and the last column of every row stops on
    the seam without passing it. The routine exists for geometry that did not
    come from the grid.
    """
    if not isinstance(geometry, Polygon):
        raise GeometryError(
            "antimeridian cutting supports Polygon input only, "
            f"got {geometry.geom_type}"
        )
    exterior = list(geometry.exterior.coords)
    lifted = Polygon(_unwrap_longitudes(exterior))
    west = lifted.intersection(box(-540.0, -90.0, 180.0, 90.0))
    east = lifted.intersection(box(180.0, -90.0, 540.0, 90.0))
    parts: list[Polygon] = []
    for piece, offset in ((west, 0.0), (east, -360.0)):
        if piece.is_empty:
            continue
        for polygon in getattr(piece, "geoms", [piece]):
            shifted = [(lon + offset, lat) for lon, lat in polygon.exterior.coords]
            parts.append(Polygon(shifted))
    if len(parts) == 1:
        return parts[0]
    return MultiPolygon(parts)


def _cell_polygon(cell: str) -> Polygon:
    """The RFC 7946 polygon of one atomic cell, closed and right-handed.

    Rings come out explicitly closed, first position identical to last, as
    section 3.1.6 requires: a triangle in four positions, a parallelogram in
    five. Orientation is normalised rather than trusted. Measured,
    :func:`itacart.cell_to_boundary` already returns counterclockwise rings
    in all four quadrants, but section 3.1.6 makes right-handedness a
    requirement on production, so the exporter enforces it instead of
    relying on a measurement of today's behaviour. That is a deliberate
    divergence from the binary serialisation, which preserves whatever
    orientation it was given because it exists to reproduce identity, while
    this module exists to interoperate.
    """
    ring = cast("list[tuple[float, float]]", cell_to_boundary(cell, close=True))
    if _holds_pole(ring):
        return orient(_polar_cap_polygon(ring), sign=1.0)
    return orient(Polygon(ring), sign=1.0)


def _cell_metadata(cell: str) -> dict[str, Any]:
    """The declared metadata block of one cell."""
    return {
        INDEX_PROPERTY: cell,
        RESOLUTION_PROPERTY: get_resolution(cell),
        SHAPE_PROPERTY: str(cell_shape(cell)),
        NOMINAL_AREA_PROPERTY: nominal_cell_area(get_resolution(cell)),
        EFFECTIVE_AREA_PROPERTY: effective_cell_area(cell),
        EXTENSION_ZONE_PROPERTY: extension_zone(cell),
    }


def cells_to_geojson(
    index: str,
    properties: dict[str, Any] | None = None,
    include_metadata: bool = True,
) -> dict[str, Any]:
    """Export a compositional index as a GeoJSON FeatureCollection.

    One cell becomes one Feature, always, and a Feature carries a Polygon or
    a MultiPolygon rather than being split into two Features. The index goes
    into the Feature ``id`` member, which RFC 7946 section 3.2 reserves for
    an identifier already in common use, and is repeated in ``properties``
    when metadata is asked for, because some consumers drop ``id`` on
    import. Coordinates are WGS84 decimal degrees and nothing else: section 4
    fixes the coordinate reference system and records that alternatives were
    removed from the specification for interoperability reasons.

    Args:
        index: Compositional index; every cell it names becomes a Feature.
        properties: Extra members merged into every Feature's ``properties``.
        include_metadata: Whether to add the ITACaRT metadata block.

    Returns:
        A FeatureCollection dictionary.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    features: list[dict[str, Any]] = []
    for cell in decompose(index):
        _require_cell(cell)
        members: dict[str, Any] = {}
        if include_metadata:
            members.update(_cell_metadata(cell))
        if properties is not None:
            members.update(properties)
        features.append(
            {
                "type": "Feature",
                "id": cell,
                "geometry": mapping(_cell_polygon(cell)),
                "properties": members,
            }
        )
    return {"type": "FeatureCollection", "features": features}


def cell_to_wkt(cell: str) -> str | list[str]:
    """Export a cell as WKT, one string per cell.

    Follows the convention of the six functions that already answer a list
    when handed a compositional index naming several cells: the list is in
    the order :func:`itacart.decompose` returns.

    Args:
        cell: Atomic cell or compositional index.

    Returns:
        A WKT string, or a list of them in decompose order.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    cells = decompose(cell)
    for atom in cells:
        _require_cell(atom)
    if len(cells) == 1:
        return str(_cell_polygon(cells[0]).wkt)
    return [str(_cell_polygon(atom).wkt) for atom in cells]


def cells_to_wkt(index: str, dissolve: bool = False) -> str:
    """Export a whole index as a single WKT geometry.

    With ``dissolve`` false the result is a GEOMETRYCOLLECTION holding one
    polygon per cell. With ``dissolve`` true adjacent cells are merged, and
    the two seams behave differently: along the prime meridian the shared
    cell is a single cell claimed geometrically from both sides, so the union
    closes without a gap; across the antimeridian, cells that are neighbours
    on the ellipsoid are not neighbours in longitude and latitude, so the
    union leaves them apart and the result stays multi-part.

    Args:
        index: Compositional index.
        dissolve: Whether to merge adjacent cells into one geometry.

    Returns:
        A WKT string.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    cells = decompose(index)
    for atom in cells:
        _require_cell(atom)
    polygons = [_cell_polygon(atom) for atom in cells]
    if dissolve:
        return str(unary_union(polygons).wkt)
    return str(
        "GEOMETRYCOLLECTION (" + ", ".join(polygon.wkt for polygon in polygons) + ")"
    )


def geometry_to_geojson(geometry: BaseGeometry) -> dict[str, Any]:
    """Serialise arbitrary geometry, cutting it at the antimeridian.

    This is the export path for geometry that did not come from the grid.
    Anything that occupies both sides of the seam is split per RFC 7946
    section 3.1.9 so that no part crosses it. Polar caps are not handled
    here; they are a cell concern and have their own routine.

    Args:
        geometry: Any Shapely polygon in WGS84 degrees.

    Returns:
        A GeoJSON geometry dictionary.

    Raises:
        GeometryError: If the geometry type is not supported.
    """
    if not isinstance(geometry, Polygon):
        raise GeometryError(f"unsupported geometry type: {geometry.geom_type}")
    lons = [lon for lon, _ in geometry.exterior.coords]
    if max(lons) - min(lons) > 180.0:
        return dict(mapping(_cut_at_antimeridian(geometry)))
    return dict(mapping(geometry))


def recover_from_geojson(obj: dict[str, Any]) -> list[str]:
    """Read back the indices that :func:`cells_to_geojson` wrote.

    This is the recovery path and it is exact: the index is in the file, in
    the Feature ``id`` member, so nothing is inferred from geometry and
    nothing is lost. It is the inverse that makes the exporter checkable
    rather than merely plausible.

    **It refuses rather than returning an empty list.** A bare geometry,
    a single Feature or a collection without features is not something
    :func:`cells_to_geojson` writes -- the exporter refuses an empty index,
    so it never writes an empty collection -- and an empty list for any of
    them would read as a file that recovered to nothing.

    Args:
        obj: A FeatureCollection produced by :func:`cells_to_geojson`.

    Returns:
        One index per Feature, in file order, never empty.

    Raises:
        GeometryError: If ``obj`` is not a FeatureCollection, if it holds
            no Feature, or if a Feature carries no recoverable index.
    """
    kind = obj.get("type") if isinstance(obj, dict) else type(obj).__name__
    if kind != "FeatureCollection":
        raise GeometryError(
            f"recover_from_geojson reads the FeatureCollection that "
            f"cells_to_geojson writes, and got {kind!r}; a geometry carries no "
            "index to recover, and filling one is from_geojson"
        )
    if not obj.get("features"):
        raise GeometryError(
            "the FeatureCollection holds no Feature, and cells_to_geojson never "
            "writes an empty one, so there is no index to recover"
        )
    indices: list[str] = []
    for feature in obj["features"]:
        index = feature.get("id")
        if index is None:
            index = feature.get("properties", {}).get(INDEX_PROPERTY)
        if index is None:
            raise GeometryError("feature carries no ITACaRT index in id or properties")
        indices.append(str(index))
    return indices


def from_geojson(
    obj: dict[str, Any],
    resolution: int,
    containment: str = "center",
) -> list[str]:
    """Fill every geometry of a GeoJSON object into a compositional index.

    This is the filling path, and it is lossy by construction: the cells that
    cover a shape are not the shape. It is named apart from
    :func:`recover_from_geojson` on purpose, because the two answer different
    questions and sharing a name with a switching argument would hide which
    one a caller got.

    There is no ``crs`` argument. RFC 7946 section 4 fixes GeoJSON to WGS84
    and records that alternative coordinate reference systems were removed
    from the specification, so accepting one here would be accepting
    something the format does not have.

    **A geometry of dimension zero is quantized, not filled.** A Point or
    a MultiPoint has no area for a fill to cover, and each of its points
    has exactly one owning cell, which :func:`itacart.cells.geo_to_cell`
    answers. A Point becomes that cell; a MultiPoint becomes the region its
    owners compose. ``containment`` is validated but decides nothing here,
    since all three predicates relate a cell to an area. Lines and polygons
    go to :func:`itacart.geometry.polyfill`, which refuses a line under
    ``center`` or ``contains`` by saying that it has no area.

    Args:
        obj: A GeoJSON Feature, FeatureCollection or bare geometry.
        resolution: Target resolution level.
        containment: Predicate passed through to
            :func:`itacart.geometry.polyfill`; not consulted for a geometry
            of dimension zero.

    Returns:
        One compositional index per input geometry, in file order. A point
        inside the cell that holds a pole at resolution 10 is answered, at
        resolutions 11 and 12, with the spelling
        :func:`itacart.cells.geo_to_cell` gives it, which can name no
        cell: a known limitation, stated in :mod:`itacart.cells`.

    Raises:
        ValueError: If ``containment`` is not one of the three predicates.
        GeometryError: If a geometry of dimension zero is empty, or if a
            fill keeps no cell. And from resolution 11, if a line or a
            polygon reaches the cell that holds a pole at resolution 10,
            whose children cannot be enumerated: a known limitation,
            stated in :mod:`itacart.hierarchy`.
        UnsupportedGeometryTypeError: If a geometry is of a type the
            reader cannot draw. The geometry library raises its own error
            for that, which is outside this package's hierarchy, so one
            ``except ITACaRTError`` around a pipeline would not have
            caught it.
    """
    if containment not in ("center", "intersects", "contains"):
        raise ValueError(
            f"containment must be 'center', 'intersects' or 'contains', "
            f"got {containment!r}"
        )
    if obj.get("type") == "FeatureCollection":
        geometries = [feature["geometry"] for feature in obj.get("features", [])]
    elif obj.get("type") == "Feature":
        geometries = [obj["geometry"]]
    else:
        geometries = [obj]
    filled = []
    for geometry in geometries:
        try:
            drawn = shape(geometry)
        except ShapelyError as exc:
            # Shapely's own error for a type it does not know. Letting it
            # through means one ``except ITACaRTError`` no longer guards the
            # pipeline, and this route is why nobody noticed: the totality
            # sweep classifies consumers by their first parameter's name,
            # this one is called ``obj``, and the sweep never reached it.
            raise UnsupportedGeometryTypeError(
                f"geometry is not a GeoJSON type this package reads: {exc}"
            ) from exc
        if drawn.geom_type in ("Point", "MultiPoint"):
            filled.append(_quantize_points(drawn, resolution))
            continue
        filled.append(
            polyfill(
                drawn,
                resolution,
                containment=cast("Containment", containment),
            )
        )
    return filled


def _quantize_points(drawn: BaseGeometry, resolution: int) -> str:
    """The owning cells of a Point or MultiPoint, as one index."""
    if drawn.is_empty:
        raise GeometryError(f"an empty {drawn.geom_type} has no position to quantize")
    points = [drawn] if drawn.geom_type == "Point" else list(drawn.geoms)
    return compose([geo_to_cell(point.x, point.y, resolution) for point in points])


def to_geodataframe(index: str, crs: str = "EPSG:4326") -> "geopandas.GeoDataFrame":
    """Export a compositional index as a GeoDataFrame.

    Requires the ``geo`` extra, which brings ``pyproj`` in transitively. The
    package does not depend on it: its geodesy is computed directly, and the
    extra pulls in a reprojection library for the caller's benefit only.

    Args:
        index: Compositional index.
        crs: Target coordinate reference system; reprojected from EPSG:4326.

    Returns:
        One row per cell, carrying the declared metadata columns.

    The index is checked before the extra is imported, so an address
    naming no cell is refused with the package's own exception whether or
    not GeoPandas is installed. The other order made the answer depend on
    the environment: with the extra present the caller was told the cell
    does not exist, and without it the caller was told to install
    GeoPandas, which was true and was not the problem.

    Raises:
        NonExistentCellError: If a named cell does not exist.
        ImportError: If GeoPandas is not installed.
    """
    collection = cells_to_geojson(index, include_metadata=True)
    geopandas_module = _import_geopandas()
    frame = geopandas_module.GeoDataFrame.from_features(
        collection["features"], crs="EPSG:4326"
    )
    if crs != "EPSG:4326":
        frame = frame.to_crs(crs)
    return frame


def from_geodataframe(
    gdf: "geopandas.GeoDataFrame",
    resolution: int,
    containment: str = "center",
) -> list[str]:
    """Fill every geometry of a GeoDataFrame into a compositional index.

    Requires the ``geo`` extra. This is a thin adapter: it reprojects to
    EPSG:4326 and delegates to :func:`from_geojson`, so the filling logic
    exists in exactly one place, and that place is the path that needs no
    optional extra and is therefore testable in any environment.

    Args:
        gdf: Input frame in any CRS; reprojected to EPSG:4326 as needed.
        resolution: Target resolution level, 1 to 13.
        containment: Predicate passed through to
            :func:`itacart.geometry.polyfill`.

    Returns:
        One compositional index per row, in row order, as
        :func:`from_geojson` answers it, the known limitation of a point
        inside the cell that holds a pole at resolution 10 included.

    Raises:
        ImportError: If GeoPandas is not installed.
        GeometryError: As :func:`from_geojson` raises it, the refusal of a
            geometry that reaches the cell that holds a pole at resolution
            10 included.
    """
    _import_geopandas()
    frame = gdf if gdf.crs is None else gdf.to_crs("EPSG:4326")
    return from_geojson(
        dict(frame.__geo_interface__), resolution, containment=containment
    )


def _import_geopandas() -> Any:
    """Import GeoPandas, or say plainly which extra is missing."""
    try:
        import geopandas as geopandas_module
    except ImportError as exc:  # pragma: no cover - exercised by injection
        raise ImportError(
            "to_geodataframe and from_geodataframe need the 'geo' extra: "
            "pip install itacart[geo]"
        ) from exc
    return geopandas_module
