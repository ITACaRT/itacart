"""Vector geometry against the grid: filling, vertex mapping, canonical form.

Two ways to represent a vector feature, both described in section 4 of the
paper:

**Cell filling** - the feature is the set of cells it covers. Area follows
from a cell count with no projection distortion, but the index grows
verbose at fine resolutions.

**Vertex representation** - only the cells holding the polygon vertices are
kept, as in conventional vector data. Far more compact, and the form the
binary encodings in :mod:`itacart.serialization` are built on.

The descent
-----------

Filling never tests cells one by one. It projects the geometry onto the
parallels plane, shears the 45-degree tiling onto a square lattice exactly
as :mod:`itacart.cells` does, and walks the hierarchy: a node disjoint from
the geometry is dropped whole, a node wholly inside is accepted whole, and
only a node straddling the outline is subdivided. Cost therefore tracks the
outline, not the area.

The shear is the same one the quantizer uses, ``u = |x| + |y|`` and
``v = |y|``, applied per quadrant as a signed affine map so that the
mirrored quadrants keep their own handedness. Its Jacobian is 1, so a cell
count on the lattice is a cell count on the ellipsoid.

Scope
-----

Filling runs on the parallelogram interior and on every family that is
not the sheared square the ordinary descent tests. The prime-meridian
column holds triangles and has a descent of its own, see
:func:`_fill_meridian_node`; the last lattice column of a row, which
absorbs the border strip, the polar row and the two caps are walked by
:func:`_fill_border_root`. No family is refused.

Provenance: ``itacart_core/cell_filling.py``, ``densification.py``,
``geometry_blob.py`` (``canonicalize_rings``) and
``cadastral_processor/vertex_extractor.py``. Two of the seven functions
have no counterpart there: ``cells_to_geometry`` is absent from all four
files, and the origin's filling implements one containment mode of the
three.
"""

from __future__ import annotations

import math
import threading
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, Final, Literal, Sequence, TypeVar

from .boundary import _CLIP_EPSILON_M, crosses_antemeridian, last_lattice_column
from .cells import _ROW_LETTERS, cell_to_centroid, geo_to_cell
from .constants import (
    DESCENT_CLOSE,
    DESCENT_OPEN,
    MAX_RESOLUTION,
)
from .constants import MERIDIAN_QUADRANT as _MERIDIAN_QUADRANT
from .constants import QUADRANTS, RES1_DIGITS, RES1_SEPARATOR, SIBLING_SEPARATOR
from .exceptions import (
    AntemeridianError,
    DensificationError,
    DomainError,
    GeometryError,
    NonExistentCellError,
    ResolutionError,
    UnsupportedGeometryTypeError,
)
from .geodesy import direct_geodesic, geodetic_to_sinusoidal, inverse_geodesic
from .resolutions import cell_size, linear_refinement_ratio, refinement_ratio

if TYPE_CHECKING:
    from shapely.geometry import LineString, Point, Polygon
    from shapely.geometry.base import BaseGeometry
    from shapely.prepared import PreparedGeometry

__all__ = [
    "polyfill",
    "count_internal_cells",
    "vertex_to_cell",
    "cells_to_geometry",
    "densify_orthodromic",
    "densify_segment",
    "canonicalize_rings",
]

Containment = Literal["center", "intersects", "contains"]

_L1: float = cell_size(1)
"""Resolution-1 cell side, in metres."""

MAX_FILL_CELLS: int = 200_000_000
"""Ceiling on cells emitted by one :func:`polyfill` call.

A fill is bounded by the geometry's area divided by the cell area, which
at resolution 13 is one square centimetre: a single football pitch holds
seventy billion cells. Without a ceiling the call does not fail, it
simply never returns, and the caller cannot tell the two apart.

The number is measured, not chosen. The index a fill returns costs about
three bytes per cell, so two hundred million cells is a string of some
six hundred megabytes and a peak of roughly twice that while it is being
joined. That is large but survivable on an ordinary machine, and it is
the point at which the output itself, rather than any bookkeeping around
it, becomes the thing that does not fit.

Raising it further is a memory question, not a policy one: measure
``len(polyfill(...))`` against available memory before doing so.
"""

# --------------------------------------------------------------------------
# Plane and lattice helpers
# --------------------------------------------------------------------------


_QUADRANT_SHEAR: dict[str, tuple[float, float, float, float]] = {
    "NE": (1.0, 1.0, 0.0, 1.0),
    "NW": (-1.0, 1.0, 0.0, 1.0),
    "SE": (1.0, -1.0, 0.0, -1.0),
    "SW": (-1.0, -1.0, 0.0, -1.0),
}
"""Per-quadrant affine shear ``(a, b, d, e)`` taking plane ``(x, y)`` to ``(u, v)``.

``u = a*x + b*y`` and ``v = d*x + e*y``. Each is the composition of the
quadrant mirror with the unit shear of :mod:`itacart.cells`, so that
``u = |x| + |y|`` and ``v = |y|`` hold inside the quadrant. The determinant
is 1 in absolute value for all four, which is why a count on the lattice is
a count on the ellipsoid.

Written per quadrant rather than as ``abs()`` because ``abs()`` folds the
four quadrants onto one and loses which one a piece came from. The origin
folds, and recovers the quadrant afterwards; recovering it afterwards is
what makes a cross-quadrant geometry impossible to fill.
"""

_PLANE_LIMIT: float = 4.0e7
"""Half-width of the clipping box, in metres: past any plane coordinate."""

_NON_AREAL: Final[frozenset[str]] = frozenset(
    {"Point", "LineString", "MultiPoint", "MultiLineString"}
)
"""Geometry types a quadrant clip can produce that carry no area.

From an areal input such a piece means the geometry touched the axis
without crossing it. From a point or a line it means the feature itself,
which is why the two cases cannot share one rule.
"""


_EXTENSION_EDGE_SEGMENTS: Final[int] = 4096
"""Segments along a zone window's meridian edges.

A meridian is a curve on the sinusoidal plane, so the straight edge of a
box in degrees projects to a chord that misses it badly: measured, 26.6 km
over Fiji's latitude band and 20.2 km over Chukotka's, more than two
resolution-1 cells. The chord error falls with the square of the segment
count, so this many puts it under a millimetre.
"""


@lru_cache(maxsize=1)
def _extension_windows() -> tuple["BaseGeometry", ...]:
    """The region each extension zone lifts past the antemeridian.

    Built from the rows the zone is realized on, **not** from its declared
    degrees. ``ZONE_ROWS`` rounds the declared band outward to whole rows
    -- Fiji is declared -21.5 to -15.5 and realized on -21.5141 to
    -15.4610 -- and :func:`itacart.cells.geo_to_cell` shifts a position by
    the rows. A window drawn on the declared numbers would leave a strip
    up to a row wide where the point function shifts and the fill does not.

    Only the meridian edges are densified. A parallel is a straight line on
    the plane, because the ordinate depends on latitude alone.
    """
    from shapely.geometry import Polygon

    from .boundary import ZONE_ROWS
    from .constants import ANTEMERIDIAN_LON, EXTENSION_ZONES
    from .geodesy import sinusoidal_to_geodetic

    windows: list["BaseGeometry"] = []
    for name, (first, last) in ZONE_ROWS.items():
        spec = EXTENSION_ZONES[name]
        sign = 1.0 if spec.quadrant[0] == "N" else -1.0
        low, high = sorted(
            sinusoidal_to_geodetic(0.0, sign * ordinate)[1]
            for ordinate in (first * _L1, (last + 1) * _L1)
        )
        steps = _EXTENSION_EDGE_SEGMENTS
        span = high - low
        limit = [
            (spec.lon_limit, low + span * step / steps) for step in range(steps + 1)
        ]
        line = [
            (-ANTEMERIDIAN_LON, high - span * step / steps) for step in range(steps + 1)
        ]
        windows.append(Polygon(limit + line))
    return tuple(windows)


def _lift_extensions(geometry: "BaseGeometry") -> "BaseGeometry":
    """Move the extension-zone parts of a geometry past the antemeridian.

    Two zones let a row's domain run past 180 degrees so that land beyond
    the line is reached from the near side rather than by crossing it.
    :func:`itacart.cells.geo_to_cell` honours that by adding 360 degrees to
    a position inside a zone, which puts it in the eastern quadrant; the
    fill projected at the literal longitude instead, so the same ground
    landed in the **western** quadrant at a column past the end of its row,
    and every cell it named there failed :func:`itacart.is_valid_cell`.

    The lift is the areal form of the same rule. A geometry that only
    partly enters a zone is cut at the zone's own edge and lifted in part,
    because the limit is where the extension stops: a position at
    ``lon_limit`` is lifted and one just east of it is not.

    Runs after densification, so the lifted part keeps the vertices the
    densifier gave it, and the cut edge follows the meridian by taking the
    window's own vertices rather than a straight line in degrees.
    """
    from shapely.affinity import translate
    from shapely.ops import unary_union

    from .constants import ANTEMERIDIAN_LON

    lifted: list["BaseGeometry"] = []
    remainder = geometry
    areal = geometry.geom_type in {"Polygon", "MultiPolygon"}
    for window in _extension_windows():
        inside = remainder.intersection(window)
        if inside.is_empty:
            continue
        if areal and inside.geom_type in _NON_AREAL:
            continue  # the geometry touches the zone without entering it
        lifted.append(translate(inside, xoff=2.0 * ANTEMERIDIAN_LON))
        remainder = remainder.difference(window)
    if not lifted:
        return geometry
    if not remainder.is_empty:
        lifted.append(remainder)
    united = unary_union(lifted)
    if areal and united.geom_type == "GeometryCollection":
        # Cutting a densified boundary along the zone's own edge leaves
        # debris behind when the boundary lies *on* that edge: the union
        # comes back as the figure plus a train of zero-area linestrings,
        # one per densified vertex that fell on the cut. The rule is the
        # one already applied to the intersection above, in the other
        # direction: an areal geometry stays areal, and a part with no
        # area is not part of the figure. It is a rule about type rather
        # than about size, so no threshold has to be chosen and none can
        # go stale.
        areal_parts = [
            part
            for part in united.geoms
            if part.geom_type in {"Polygon", "MultiPolygon"}
        ]
        united = unary_union(areal_parts)
    return united


def _closes_at_a_pole(geometry: "BaseGeometry") -> bool:
    """Does this figure reach a pole, and so close around it?

    The distinction the antemeridian screen needs, and the smallest one
    that separates the two figures it was conflating.

    The polar cap runs from -180 to 180 because a cap has to: its two
    meridian edges are the *same* edge, and they meet at the pole, which
    is one point rather than a seam. An ordinary crossing has its two
    sides of the line at the same latitude, never meeting, and splitting
    it at 180 is exactly the right advice.

    Reaching the pole is what tells them apart, and it is a property of
    the figure rather than of its longitude span, so a caller cannot
    dress a crossing up as a cap by widening it.
    """
    from .constants import WGS84_A

    _, min_lat, _, max_lat = geometry.bounds
    slack = math.degrees(_CLIP_EPSILON_M / WGS84_A)
    return bool(max_lat >= 90.0 - slack or min_lat <= -90.0 + slack)


def _refuse_invalid_area(geometry: "BaseGeometry") -> None:
    """Refuse a polygon or multipolygon that is not a valid geometry.

    Validity is a property of the geometry as the caller wrote it, not of
    where it lies, so it is asked before any screen, clip or fill reads the
    geometry. Asked any later, an invalid outline is treated according to
    what the machinery happens to do with it: the quadrant split refuses a
    bowtie, while a polygon with a zero-width spike past the antemeridian
    passes the extension screen, which reads a part with no area as a
    polygon touching the far edge of a zone, and is filled from its valid
    remainder.

    Nothing is repaired. A zero-width buffer or ``make_valid`` would decide
    which area was meant -- one lobe of a bowtie or both, or a line where a
    polygon collapsed -- and that decision belongs to the caller.

    Geometries without area are not asked. A line that crosses itself is a
    valid line, and points and lines keep contracts of their own.
    """
    if geometry.geom_type in {"Polygon", "MultiPolygon"} and not geometry.is_valid:
        from shapely.validation import explain_validity

        raise GeometryError(
            f"the {geometry.geom_type} is not a valid geometry "
            f"({explain_validity(geometry)}); an invalid areal geometry is "
            "refused rather than repaired, because a repair would decide "
            "which area was meant"
        )


def _refuse_unzoned_extension(geometry: "BaseGeometry") -> None:
    """Refuse longitude past the line that no extension zone legitimises.

    The domain stops at 180 degrees except where a zone carries it
    further, and only in the rows that zone is realized on. A footprint
    written past the line anywhere else names ground the grid does not
    address, and it is refused here, by the antemeridian's own name.

    The rule has to be stated here rather than left to the fill. Such a
    footprint projects into the western quadrant at a column past the end
    of its row, which is the last lattice column, and the fill descends
    that column like any other: nothing further on would refuse it.

    The legitimate region is the lifted window itself, so the two cannot
    drift: :func:`_lift_extensions` moves a zone's ground east by a full
    turn, and what is allowed past the line is exactly where that lands.
    """
    from shapely.affinity import translate
    from shapely.geometry import box
    from shapely.ops import unary_union

    from .constants import ANTEMERIDIAN_LON, WGS84_A

    # A cell of the last column at the equator has its border *on* the
    # line, and the round trip through the projection puts its far
    # vertex 2.8e-14 degrees past it. That is arithmetic noise, not a
    # footprint reaching past the domain, so the line is given the
    # width the boundary module already clips with -- 1e-6 metres, or
    # 9e-12 degrees at the equator, three orders above the noise -- and
    # no new threshold is chosen here.
    slack = math.degrees(_CLIP_EPSILON_M / WGS84_A)
    edge = ANTEMERIDIAN_LON + slack
    min_lon, _, max_lon, _ = geometry.bounds
    if min_lon >= -edge and max_lon <= edge:
        return
    inside_the_line = box(-edge, -90.0, edge, 90.0)
    past = geometry.difference(inside_the_line)
    allowed = unary_union(
        [
            translate(window, xoff=2.0 * ANTEMERIDIAN_LON)
            for window in _extension_windows()
        ]
    )
    # The window's edges are the zone's own row boundaries, and a cell of
    # the first or last row of a zone lands on one of them, so the same
    # projection noise puts a hair of it outside. The allowance carries
    # the same slack the line does, for the same reason.
    stray = past.difference(allowed.buffer(slack))

    if stray.is_empty:
        return
    raise AntemeridianError(
        "geometry reaches past 180 degrees longitude where no extension "
        "zone carries the domain; only the Fiji and Chukotka rows are "
        "addressable past the line, and only as far as each one reaches"
    )


def _check_resolution(resolution: int) -> None:
    """Reject a resolution outside 1..13."""
    if not isinstance(resolution, int) or isinstance(resolution, bool):
        raise ResolutionError(
            f"resolution must be an int, got {type(resolution).__name__}"
        )
    if not 1 <= resolution <= MAX_RESOLUTION:
        raise ResolutionError(f"resolution {resolution} outside 1..{MAX_RESOLUTION}")


def _check_jobs(n_jobs: int) -> None:
    """Reject a worker count below one."""
    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool):
        raise ValueError(f"n_jobs must be an int, got {type(n_jobs).__name__}")
    if n_jobs < 1:
        raise ValueError(f"n_jobs must be >= 1, got {n_jobs}")


def _leaves_between(level: int, target: int) -> int:
    """How many ``target``-resolution cells sit under one ``level`` cell."""
    total = 1
    for step in range(level + 1, target + 1):
        total *= refinement_ratio(step)
    return total


def _child_code(row: int, column: int, level: int) -> str:
    """Refinement code of the child at ``(row, column)`` of a ``level`` grid.

    Rows run south to north and columns west to east *in the cell's own
    quadrant*, which is the frame the descent works in. The quaternary
    code is ``row * 2 + column + 1``; the quinary code is the row letter
    followed by the one-based column. Both are the inverse of the mapping
    :func:`itacart.cells._anchor_on_plane` applies when it ascends, and
    the letters come from that module so the north-south orientation
    keeps a single definition in the package.
    """
    if linear_refinement_ratio(level) == 2:
        return str(row * 2 + column + 1)
    return f"{_ROW_LETTERS[row]}{column + 1}"


def _project(geometry: "BaseGeometry") -> "BaseGeometry":
    """Map a geodetic geometry onto the parallels plane, vertex by vertex.

    ``shapely.ops.transform`` is not used: the projection is applied to
    the coordinates that are there, and no vertex is invented. A caller
    wanting the edges to follow geodesics densifies first, which is what
    :func:`polyfill` does on its behalf.
    """
    from shapely.geometry import LineString, MultiPolygon, Point, Polygon

    if isinstance(geometry, Point):
        return Point(*geodetic_to_sinusoidal(geometry.x, geometry.y))
    if isinstance(geometry, LineString):
        return LineString([geodetic_to_sinusoidal(x, y) for x, y in geometry.coords])
    if isinstance(geometry, Polygon):
        return Polygon(
            [geodetic_to_sinusoidal(x, y) for x, y in geometry.exterior.coords],
            [
                [geodetic_to_sinusoidal(x, y) for x, y in hole.coords]
                for hole in geometry.interiors
            ],
        )
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon([_project(part) for part in geometry.geoms])
    raise UnsupportedGeometryTypeError(
        f"cannot fill a {geometry.geom_type}; supported types are Point, "
        "LineString, Polygon and MultiPolygon"
    )


def _quadrant_pieces(plane: "BaseGeometry") -> list[tuple[str, "BaseGeometry"]]:
    """Split a plane geometry into its four quadrant parts, dropping empties.

    A geometry lying in one quadrant yields one piece, which is the
    ordinary case. A geometry straddling an axis yields two or four, each
    filled independently.

    **The windows are closed and the convention is not.** Shapely can only
    clip with a closed box, so all four windows carry their own boundary
    and an axis position lands in two of them, or in four at the origin.
    :func:`itacart.cells.geo_to_cell` awards that position to exactly one
    cell -- the eastern and northern side, measured -- so the clip has to
    be brought back to the same rule or the two functions disagree about
    who owns the axis.

    Two different corrections are needed because the leak takes two forms:

    * An areal geometry that only touches an axis leaves a degenerate
      piece on the far side, which is dropped.
    * A point or a line **on** the axis is not degenerate; it is the whole
      feature, and both windows keep all of it. Each quadrant therefore
      subtracts the axes it does not own. The subtraction is exact rather
      than approximate: a point on the line vanishes, a line along it
      vanishes, and a line merely crossing it keeps its full length.

    Areal pieces are left alone because subtracting a line from a polygon
    is a no-op on the closure, so the cost would buy nothing.
    """
    from shapely.geometry import LineString, box

    out: list[tuple[str, "BaseGeometry"]] = []
    limit = _PLANE_LIMIT
    windows = {
        "NE": box(0.0, 0.0, limit, limit),
        "NW": box(-limit, 0.0, 0.0, limit),
        "SE": box(0.0, -limit, limit, 0.0),
        "SW": box(-limit, -limit, 0.0, 0.0),
    }
    meridian = LineString([(0.0, -limit), (0.0, limit)])
    equator = LineString([(-limit, 0.0), (limit, 0.0)])
    disowned: dict[str, tuple["BaseGeometry", ...]] = {
        "NE": (),
        "NW": (meridian,),
        "SE": (equator,),
        "SW": (meridian, equator),
    }
    for quadrant, window in windows.items():
        try:
            piece = plane.intersection(window)
        except Exception as exc:
            # The geometry library refuses a self-intersecting plane
            # geometry, and its exception is not one this package owns.
            #
            # It used to arrive from inside: a polar-row cell whose own
            # densified boundary folded across the globe, because the
            # walk over the pole put the far half on the opposite
            # longitude branch from the vertex it was walking towards.
            # That is repaired in :func:`densify_segment`, and no cell of
            # the grid reaches here any more.
            #
            # What does reach here is a caller's outline that crosses
            # itself. Such a figure has no inside for a clip to keep, and
            # the package owes it a name rather than the engine's.
            raise GeometryError(
                "the boundary cannot be split by quadrant because it is "
                "self-intersecting; an outline that crosses itself has no "
                "interior for the split to keep"
            ) from exc
        if piece.is_empty:
            continue
        if piece.geom_type in _NON_AREAL:
            if plane.geom_type in {"Polygon", "MultiPolygon"}:
                continue  # a polygon touching the axis, not crossing it
            for axis in disowned[quadrant]:
                piece = piece.difference(axis)
            if piece.is_empty:
                continue
        out.append((quadrant, piece))
    return out


def _to_lattice(piece: "BaseGeometry", quadrant: str) -> "BaseGeometry":
    """Shear a quadrant piece onto the square lattice."""
    from shapely.affinity import affine_transform

    a, b, d, e = _QUADRANT_SHEAR[quadrant]
    return affine_transform(piece, [a, b, d, e, 0.0, 0.0])


#: The one row the pole falls inside. Rows are ``_L1`` tall from the
#: equator, so the pole at ``MERIDIAN_QUADRANT`` lands in this one and in
#: no other.
#:
#: It used to be derived, as "the row above addresses no cell", read off
#: ``last_lattice_column`` answering zero or less. That answer saturates:
#: it is zero for row 1000 and zero again for row 1001, and zero means
#: one cell in an eastern quadrant and none in a western one. So row 999
#: was called polar in all four quadrants, and its cells -- valid,
#: well-formed, six real children each, their highest point 1 966 metres
#: below the pole -- were refused for a property they do not have.
_POLAR_ROW: Final[int] = int(math.floor(_MERIDIAN_QUADRANT / _L1))


def _first_lattice_column(quadrant: str) -> int:
    """The lowest column index a quadrant addresses.

    Column zero is the meridian column and the meridian belongs to the
    east, so the western quadrants start at one.
    """
    return 0 if quadrant[1] == "E" else 1


def _check_addressable(quadrant: str, column: int, row: int) -> None:
    """Refuse a resolution-1 cell that is not an ordinary parallelogram.

    Two families are refused, and the reason is the same in both: the
    descent tests a sheared square, and in these families the cell is not
    that square.

    * the last lattice column of a row absorbs the strip between itself
      and the domain border, so it is wider than the square.
    * the polar row is clipped by the pole and carries a fraction of the
      nominal area.

    **Column 0 is not among them.** The prime-meridian column holds
    triangles, which are not squares either, but they are filled by
    :func:`_fill_meridian_node` rather than refused. Negative columns
    never arrive: :func:`_base_cells` drops them, because a square west
    of the meridian in an eastern frame is an artefact of a closed
    clipping window and names no cell.

    The screen is applied at resolution 1 and is therefore conservative:
    a geometry inside base column 500 cannot reach a border-absorbing
    cell at any finer resolution, so refusing the whole base column
    refuses more than strictly necessary and never less.

    **Order matters.** The polar row is tested before the last column,
    because the polar row *is* a row whose only column is its last one:
    at resolution 1 it holds exactly one, so a last-column test placed
    first would answer every polar position with the wrong diagnosis and
    the polar branch would never be reached at all. Measured: the branch
    was unreachable until the two were swapped.
    """
    if row > _POLAR_ROW:
        raise DomainError(
            f"row {row} of quadrant {quadrant} addresses no cell: it lies "
            "beyond the pole"
        )
    last = last_lattice_column(quadrant, row, _L1)
    if last < _first_lattice_column(quadrant):
        raise DomainError(
            f"row {row} of quadrant {quadrant} addresses no cell: its "
            "parallel circle is shorter than one cell side"
        )
    if row == _POLAR_ROW:
        raise DomainError(
            f"row {row} of quadrant {quadrant} is the polar row, which is "
            "clipped by the pole and does not carry the nominal cell area"
        )
    if column >= last:
        raise NonExistentCellError(
            f"row {row} column {column} of quadrant {quadrant} is the last "
            "lattice column of its row, which absorbs the border strip and "
            "is not a parallelogram"
        )


def _anomalous_band(
    quadrant: str, column: int, row: int, side: float
) -> "BaseGeometry":
    """The part of a resolution-1 cell that is still anomalous at ``side``.

    The two families are bounded by lines of constant ``w = u - v`` or
    of constant ``v``, so each one is a band and their union is what the
    screen has to protect. Returned in lattice coordinates.

    Derived, and validated by enumeration against
    :func:`_check_addressable` over 88,880 nodes in four quadrants. The
    lattice column of a node of side ``s`` at ``(u0, v0)`` is
    ``(u0 - v0) / s``, which reduces to ``u_index - row`` at resolution
    1. From it, ``column`` of a child is ``column * d + (c - r)`` with
    ``c`` and ``r`` in ``[0, d)``, so a node with a positive column never
    fathers one at column zero: neither family spreads toward the
    meridian.

    Each band is drawn one cell side wide, which is an outer bound on a
    family that is really a staircase of single cells. Refusing the
    bound refuses at most one extra cell of width and never one less,
    and one cell of width is the whole point of the exercise.
    """
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union

    u0 = (column + row) * _L1
    v0 = row * _L1
    cell = box(u0, v0, u0 + _L1, v0 + _L1)
    reach = 4.0 * _L1
    bands: list["BaseGeometry"] = []

    def by_column(lo: float, hi: float) -> "BaseGeometry":
        """The band ``lo <= u - v <= hi``, as a sheared quadrilateral."""
        return Polygon(
            [
                (v0 - reach + lo, v0 - reach),
                (v0 - reach + hi, v0 - reach),
                (v0 + reach + hi, v0 + reach),
                (v0 + reach + lo, v0 + reach),
            ]
        )

    top_row = int(math.floor((v0 + _L1 - side / 2.0) / side))
    last_at_top = last_lattice_column(quadrant, top_row, side)
    if last_at_top > 0:
        # A square of side ``s`` at column ``c`` spans ``u - v`` over
        # ``[(c - 1) s, (c + 1) s]``, so every square from the last column
        # outward lies at ``u - v >= (last - 1) s``. Derived rather than
        # chosen, and the two ends of it are what the screen needs:
        #
        # * at resolution 1 that threshold is the base cell's own lower
        #   edge, so the absorbing cell is refused whole rather than
        #   emitted whenever the geometry happens to sit on its inner
        #   side;
        # * at a finer target it is a strip beside the border, which is
        #   what lets the rest of the base cell still be filled.
        #
        # The outward end is drawn past the cell rather than one band
        # wide. A column beyond the last is not one cell out but as many
        # as the border strip is wide -- two, over Chukotka -- and a band
        # that stopped short of them left an empty intersection, so the
        # escape fired and the fill named cells that do not exist.
        lower = (last_at_top - 1) * side
        bands.append(by_column(lower, max(lower, (column + 2) * _L1)))
    else:
        bands.append(cell)

    polar_v = None
    probe = int(math.floor(v0 / side))
    if last_lattice_column(quadrant, probe + 1, side) <= 0:
        polar_v = v0
    elif last_lattice_column(quadrant, top_row + 1, side) <= 0:
        lo, hi = probe, top_row
        while lo < hi:
            mid = (lo + hi) // 2
            if last_lattice_column(quadrant, mid + 1, side) <= 0:
                hi = mid
            else:
                lo = mid + 1
        polar_v = lo * side
    if polar_v is not None:
        bands.append(box(u0 - reach, polar_v, u0 + reach, v0 + _L1))

    return unary_union([band.intersection(cell) for band in bands])


def _base_cells(
    prepared: "PreparedGeometry",
    lattice: "BaseGeometry",
    quadrant: str,
    resolution: int,
) -> list[tuple[int, int, bool]]:
    """Resolution-1 ``(column, row)`` pairs whose square meets ``lattice``.

    Candidates come from the bounding box, but only those a square
    actually meets are screened. Screening the whole box instead would
    refuse a geometry for a neighbouring column it never touches.

    The screen runs at ``resolution``, not at resolution 1. A base cell
    of one of the refused families holds mostly ordinary descendants --
    measured, 107 of 114 for the last column and 139 of 144 for the
    polar row, and 90 of 100 for the meridian before that family left
    the screen entirely -- and screening the family at
    resolution 1 refused all of them along with the few that deserved
    it. That is a band up to ten kilometres wide standing in for an
    anomaly one cell wide: a thousand times too much at resolution 7 and
    a million times at resolution 13.

    So an anomalous base cell is not refused for being one. Its
    anomalous band at the target side is computed and the geometry is
    refused only if it reaches it.
    """
    from shapely.geometry import box

    min_u, min_v, max_u, max_v = lattice.bounds
    row_lo = int(math.floor(min_v / _L1))
    row_hi = int(math.floor(max_v / _L1))
    u_lo = int(math.floor(min_u / _L1))
    u_hi = int(math.floor(max_u / _L1))
    out: list[tuple[int, int, bool]] = []
    for row in range(row_lo, row_hi + 1):
        for u_index in range(u_lo, u_hi + 1):
            u0 = u_index * _L1
            v0 = row * _L1
            if not prepared.intersects(box(u0, v0, u0 + _L1, v0 + _L1)):
                continue
            column = u_index - row
            if column <= 0:
                # A negative column is a square west of the meridian read
                # in an eastern frame: the clipping window is closed, so
                # such a square touches the piece at one corner and names
                # no cell. Column 0 names the meridian triangle, which
                # _fill_meridian_node walks in the plane rather than here.
                continue
            # A candidate position past the end of its row names no
            # cell. The enumeration produces such positions because it
            # walks a bounding box in lattice coordinates, which does not
            # know where each row stops, and they were reaching
            # ``_check_addressable`` and drawing that row's diagnosis --
            # so a valid cell of row 998 was cancelled by a position of
            # row 999 that is outside the grid. Dropping them is not a
            # weakening of the screen: there is nothing there to fill.
            reach = last_lattice_column(quadrant, row, _L1)
            if column > reach or reach < _first_lattice_column(quadrant):
                continue
            try:
                _check_addressable(quadrant, column, row)
            except NonExistentCellError:
                # Only the last lattice column is refused here. The other
                # refusals of _check_addressable concern a whole row -- past
                # the pole, a parallel shorter than one cell, the polar row
                # -- and none of those rows offers a positive column, so the
                # filter above has already dropped them. Anything else
                # propagates. The last column is walked rather than refused:
                # it absorbs the strip between its square and the domain
                # border, and _fill_border_node descends the tree the
                # hierarchy proved for it.
                band = _anomalous_band(quadrant, column, row, cell_size(resolution))
                if band.is_empty or not prepared.intersects(band):
                    out.append((column, row, False))
                    continue
                out.append((column, row, True))
                continue
            out.append((column, row, False))
    return out


# --------------------------------------------------------------------------
# Cell filling
# --------------------------------------------------------------------------


def _auto_segment(resolution: int) -> float:
    """Densification threshold implied by a target resolution, in metres.

    A straight line on the parallels plane is not a geodesic, and the two
    part company by a sagitta that grows with the square of the span. A
    chord of length ``d`` departs from its geodesic by roughly
    ``d^2 / (8 R)``; requiring that to stay under half a cell side gives
    ``d = sqrt(4 R l)``. At resolution 7 that is about sixteen kilometres,
    at resolution 13 about five hundred metres.

    The result is capped at one kilometre, the default segment of
    :func:`densify_orthodromic`, so this rule is never looser than that
    blanket threshold and is stricter wherever the cell asks for it.
    """
    from .constants import WGS84_A

    return min(math.sqrt(4.0 * WGS84_A * cell_size(resolution)), 1000.0)


def _accept_leaf(
    prepared: "PreparedGeometry",
    containment: Containment,
    u0: float,
    v0: float,
    side: float,
) -> bool:
    """Whether a target-resolution square is kept, given how it was reached.

    Called only for a square that the descent has found to intersect the
    geometry and *not* to be wholly inside it. Both facts are already
    established, so two of the three modes answer without a further
    predicate call, and the containment chain follows from the geometry
    of a square rather than from three independent tests:

    ``contains`` reached here means not wholly inside, so the cell is
    dropped. ``intersects`` reached here means touching, so it is kept.
    ``center`` asks the one question still open. Since the centre is a
    point of the square, wholly-inside implies centre-inside implies
    touching, and therefore ``contains`` is a subset of ``center`` is a
    subset of ``intersects`` by construction, for every square the
    descent visits.
    """
    from shapely.geometry import Point

    if containment == "contains":
        return False
    if containment == "intersects":
        return True
    return bool(prepared.contains(Point(u0 + side / 2.0, v0 + side / 2.0)))


@lru_cache(maxsize=None)
def _expansion(level: int, target: int) -> str:
    """Index suffix that expands a wholly-accepted node down to ``target``.

    Every node accepted whole at the same level expands to the same
    text, so it is built once and the callers share one string object.
    That is what keeps the interior of a fill cheap: a region wholly
    inside the geometry contributes one pointer per node rather than one
    string per cell.

    Empty when the node is already at the target resolution.
    """
    if level == target:
        return ""
    step = level + 1
    divisor = linear_refinement_ratio(step)
    inner = _expansion(step, target)
    parts = [
        _child_code(row, column, step) + inner
        for row in range(divisor)
        for column in range(divisor)
    ]
    return f"{DESCENT_OPEN}{SIBLING_SEPARATOR.join(parts)}{DESCENT_CLOSE}"


class _Budget:
    """Running cell count for one fill, checked before any expansion.

    The check has to come before :func:`_expansion` rather than after
    the fact, because expanding a resolution-1 node to resolution 13
    builds a string of a trillion cells and there is no counting that
    afterwards.
    """

    __slots__ = ("spent",)

    def __init__(self) -> None:
        self.spent = 0

    def charge(self, cells: int) -> None:
        """Add ``cells`` to the tally and refuse to pass the ceiling."""
        self.spent += cells
        if self.spent > MAX_FILL_CELLS:
            raise GeometryError(
                f"fill exceeded {MAX_FILL_CELLS} cells; pass compact=True where "
                "the geometry stays clear of the border-absorbing cells, which "
                "refuse it, or fill a smaller geometry. A coarser resolution is a last "
                "resort: in cadastral use the resolution is prescribed by the "
                "mapping scale and is not the caller's to lower"
            )


def _fill_node(
    prepared: "PreparedGeometry",
    u0: float,
    v0: float,
    side: float,
    level: int,
    target: int,
    containment: Containment,
    code: str,
    compact: bool,
    budget: _Budget,
) -> str | None:
    """Walk one node and return its index fragment, or ``None`` if empty.

    The fragment is the node's own code followed by its accepted
    descendants in brackets, which is the compositional index of the
    subtree. Building it during the descent rather than collecting
    atomic indices and composing them afterwards is not a
    micro-optimisation: the atomic form of a fill costs about two
    kilobytes of live memory per cell, against three bytes per cell for
    the composed form it is about to be turned into. Measured at 133,018
    cells: 303 MB peak for a 406 kB answer.

    Three outcomes. Disjoint: nothing below can be kept under any mode,
    because every descendant square is a subset of this one. Wholly
    inside: everything below is kept under every mode, so the node is
    emitted as itself when compacting and with a shared expansion
    otherwise. Straddling: subdivide, or apply the leaf predicate if
    this is already the target.
    """
    from shapely.geometry import box

    square = box(u0, v0, u0 + side, v0 + side)
    if not prepared.intersects(square):
        return None
    if prepared.contains(square):
        if compact or level == target:
            budget.charge(1)
            return code
        budget.charge(_leaves_between(level, target))
        return code + _expansion(level, target)
    if level == target:
        if _accept_leaf(prepared, containment, u0, v0, side):
            budget.charge(1)
            return code
        return None
    step = level + 1
    divisor = linear_refinement_ratio(step)
    child = side / divisor
    parts: list[str] = []
    for row in range(divisor):
        for column in range(divisor):
            fragment = _fill_node(
                prepared,
                u0 + column * child,
                v0 + row * child,
                child,
                step,
                target,
                containment,
                _child_code(row, column, step),
                compact,
                budget,
            )
            if fragment is not None:
                parts.append(fragment)
    if not parts:
        return None
    if compact and _folds_into_its_parent(parts, divisor):
        return code
    return f"{code}{DESCENT_OPEN}{SIBLING_SEPARATOR.join(parts)}{DESCENT_CLOSE}"


def _folds_into_its_parent(parts: list[str], divisor: int) -> bool:
    """Whether every child came back whole, so the parent stands for them.

    Compacting has to mean the same thing as compacting the uniform fill
    afterwards, so a node folds when **all of its children were kept**,
    not only when the node itself lies wholly inside the geometry.
    Wholly inside is a sufficient condition and the descent still takes
    it as a shortcut; it is not a necessary one, and treating it as one
    left the compaction incomplete -- measured, 15 per cent more cells
    than the same fill compacted afterwards under ``center`` and 80 per
    cent more under ``intersects``.

    A child came back whole when its fragment carries no descent of its
    own: either it is a leaf at the target resolution, or it is itself a
    fold. Sibling completeness is only decidable here, once every child
    has answered, which is why the test cannot be moved earlier.
    """
    return len(parts) == divisor * divisor and not any(
        DESCENT_OPEN in part for part in parts
    )


# --------------------------------------------------------------------------
# The border-absorbing family
# --------------------------------------------------------------------------


class _BorderWalk:
    """One fill's view of the absorbing family: the tree, cached once.

    The ordinary descent subdivides a sheared square into ``d x d``
    smaller ones and tests each. Neither half of that holds here: the
    cell is not the square, and the family runs two to seven children.
    So this walk borrows nothing from it -- no leaf count in closed
    form, no shared expansion string, no fold on ``d * d`` parts.

    ``hierarchy._border_children_of`` carries a bounded cache and one
    resolution-1 absorbing cell holds twenty thousand nodes by
    resolution 5, so a counting pass would evict what the fill is about
    to ask for and the tree would be paid for twice. This holds it for
    the length of one call and both passes read it.
    """

    __slots__ = ("_children", "_rings", "quadrant", "sheared")

    def __init__(self, quadrant: str, sheared: bool = True) -> None:
        self.quadrant = quadrant
        self.sheared = sheared
        self._children: dict[str, tuple[str, ...]] = {}
        self._rings: dict[str, "BaseGeometry | None"] = {}

    def children(self, cell: str) -> tuple[str, ...]:
        """The proved children of ``cell``, in canonical order."""
        from . import hierarchy

        got = self._children.get(cell)
        if got is None:
            got = tuple(hierarchy._children_of(cell))
            self._children[cell] = got
        return got

    def ring(self, cell: str) -> "BaseGeometry | None":
        """The cell's effective ring, sheared onto the lattice.

        :func:`~itacart.boundary.plane_ring` is the authority for
        effective -- it answers after the pole and after absorption --
        and the descent works sheared, so the ring is carried across
        with the matrix the query already went through. ``None`` when
        the cell names no surface.
        """
        if cell in self._rings:
            return self._rings[cell]
        from shapely.geometry import Polygon

        from .boundary import plane_ring

        _, ring = plane_ring(cell)
        if not self.sheared:
            # The meridian family straddles the line, so it lives in the
            # unsplit plane and the query was never sheared. Carrying the
            # ring across anyway would put the two in different frames.
            unsheared = Polygon(ring) if len(ring) >= 3 else None
            self._rings[cell] = unsheared
            return unsheared
        a, b, d, e = _QUADRANT_SHEAR[self.quadrant]
        shape: "BaseGeometry | None" = (
            Polygon([(a * x + b * y, d * x + e * y) for x, y in ring])
            if len(ring) >= 3
            else None
        )
        self._rings[cell] = shape
        return shape

    def square(self, cell: str, level: int) -> tuple[float, float, float]:
        """Lattice anchor and side of a cell that does not absorb.

        A cell outside the family is the sheared square the ordinary
        descent tests, so the moment a branch leaves the family it goes
        back to the closed-form walk. The anchor is read from the ring
        so the two frames cannot drift; the side comes from the
        resolution table so it cannot pick up rounding from bounds.
        """
        shape = self.ring(cell)
        assert shape is not None
        u0, v0 = shape.bounds[0], shape.bounds[1]
        return u0, v0, cell_size(level)


def _accept_border_leaf(
    prepared: "PreparedGeometry", containment: Containment, shape: "BaseGeometry"
) -> bool:
    """Whether an absorbing target-resolution cell is kept.

    Called once the cell is known to meet the query, by the cell's
    effective geometry rather than by a lattice square it is not.

    The three modes are asked independently, unlike :func:`_accept_leaf`,
    which infers two of them from the square's containment chain. That
    chain holds because a square's descendants are subsets of it. Here
    they are not: a child of an absorbing cell reaches past its parent
    by up to 2.411 per cent of its own area under the chord
    representation, so nothing may be inferred from the parent.
    """
    if containment == "intersects":
        return True
    if containment == "contains":
        return bool(prepared.contains(shape))
    return bool(prepared.contains(shape.centroid))


def _count_border_node(
    prepared: "PreparedGeometry",
    walk: _BorderWalk,
    cell: str,
    level: int,
    target: int,
    containment: Containment,
    remaining: int | None,
) -> int:
    """Count the target cells under one absorbing node, exactly.

    Exact rather than bounded, because no power of a bound describes
    this family. Measured over 572 resolution-1 parents, the real tree
    at resolution 5 runs from 749 nodes to 20 539 against a uniform
    10 000: a seven-to-the-depth ceiling refuses fills that fit, and a
    ``d`` -to-the-depth one accepts fills that cannot be materialised.

    ``remaining`` is what the fill's budget has left. The walk stops as
    soon as the running total passes it, so an oversized fill is refused
    without walking the rest of the tree. The point of the pre-count is to
    keep the output under the ceiling, not to keep the walk constant-time
    -- that was only ever possible because the uniform tree had a formula.

    ``None`` is for :func:`count_internal_cells`, which names no cell and
    so has no budget to protect. Stopping there would return the partial
    total at the moment of passing the ceiling, as if it were the count.
    """
    from .boundary import absorbs_border

    shape = walk.ring(cell)
    if shape is None or not prepared.intersects(shape):
        return 0
    if level == target:
        return 1 if _accept_border_leaf(prepared, containment, shape) else 0
    total = 0
    for child in walk.children(cell):
        if absorbs_border(child) or not walk.sheared:
            # Unsheared, the closed-form walk has no square to test: a
            # cell of the meridian family is a parallelogram or a triangle
            # in the plane, not an axis-aligned box. Every node is asked
            # about its own ring instead. The subtree is small -- 1 194
            # nodes under the cap at resolution 5 -- so that is affordable.
            total += _count_border_node(
                prepared,
                walk,
                child,
                level + 1,
                target,
                containment,
                None if remaining is None else remaining - total,
            )
        else:
            u0, v0, side = walk.square(child, level + 1)
            total += _count_node(prepared, u0, v0, side, level + 1, target, containment)
        if remaining is not None and total > remaining:
            return total
    return total


def _fill_border_node(
    prepared: "PreparedGeometry",
    walk: _BorderWalk,
    cell: str,
    level: int,
    target: int,
    containment: Containment,
) -> list[str]:
    """Every target-resolution cell the query keeps under one absorbing node.

    Atomic spellings rather than an index fragment. Of the 2 426
    resolution-2 children of the 572 lateral parents, 453 are spelled
    under a stem one or two columns east of their parent's, and that
    stem names no resolution-1 cell of its own. A fragment of the form
    ``code(children)`` cannot hold them, so the cells are returned and
    the composition is left to the caller, where the roots are grouped.

    There is no wholly-inside shortcut, for the reason
    :func:`_accept_border_leaf` gives: containing the parent says
    nothing about the children in this family.
    """
    from .boundary import absorbs_border

    shape = walk.ring(cell)
    if shape is None or not prepared.intersects(shape):
        return []
    if level == target:
        return [cell] if _accept_border_leaf(prepared, containment, shape) else []
    kept: list[str] = []
    for child in walk.children(cell):
        if absorbs_border(child) or not walk.sheared:
            kept.extend(
                _fill_border_node(prepared, walk, child, level + 1, target, containment)
            )
        else:
            u0, v0, side = walk.square(child, level + 1)
            kept.extend(
                _fill_ordinary_atoms(
                    prepared, u0, v0, side, level + 1, target, containment, child
                )
            )
    return kept


def _fill_border_root(
    prepared: "PreparedGeometry",
    quadrant: str,
    column: int,
    row: int,
    target: int,
    containment: Containment,
    compact: bool,
    budget: _Budget,
    sheared: bool = True,
) -> str | None:
    """Fill one resolution-1 absorbing cell and spell what it keeps.

    The pre-count runs first and is charged against what the budget has
    left, so an oversized fill is refused before any cell is named. Both
    passes read one ``_BorderWalk``, so the tree is derived once.

    The fragment is built from the atoms rather than during the descent,
    because 453 of the 2 426 resolution-2 children of the lateral family
    are spelled under a stem one or two columns east of their parent's.
    Those roots name no resolution-1 cell of their own -- they are past
    the last lattice column -- so they never collide with a base cell,
    and :func:`~itacart.index.compose` groups them.
    """
    from .index import compose

    if compact:
        raise GeometryError(
            "compact=True is not available over the border-absorbing "
            "family: the fold rule reads a uniform tree, and this family "
            "does not have one. Pass compact=False"
        )
    walk = _BorderWalk(quadrant, sheared=sheared)
    root = (
        f"{quadrant}{DESCENT_OPEN}{column:0{RES1_DIGITS}d}"
        f"{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}{DESCENT_CLOSE}"
    )
    remaining = MAX_FILL_CELLS - budget.spent
    budget.charge(
        _count_border_node(prepared, walk, root, 1, target, containment, remaining)
    )
    kept = _fill_border_node(prepared, walk, root, 1, target, containment)
    if not kept:
        return None
    composed = compose(kept)
    return composed[len(quadrant) + 1 : -1]


def _fill_ordinary_atoms(
    prepared: "PreparedGeometry",
    u0: float,
    v0: float,
    side: float,
    level: int,
    target: int,
    containment: Containment,
    cell: str,
) -> list[str]:
    """The closed-form descent, resumed once a branch leaves the family.

    Spelled as atoms rather than as a fragment because the caller is
    collecting cells across roots that do not share a stem. The walk
    itself is the ordinary one: a square, its ``d x d`` children, and
    the wholly-inside shortcut, which is sound again here because a
    square's descendants are subsets of it.
    """
    from shapely.geometry import box

    square = box(u0, v0, u0 + side, v0 + side)
    if not prepared.intersects(square):
        return []
    if level == target:
        if prepared.contains(square) or _accept_leaf(
            prepared, containment, u0, v0, side
        ):
            return [cell]
        return []
    step = level + 1
    divisor = linear_refinement_ratio(step)
    child = side / divisor
    out: list[str] = []
    for row in range(divisor):
        for column in range(divisor):
            out.extend(
                _fill_ordinary_atoms(
                    prepared,
                    u0 + column * child,
                    v0 + row * child,
                    child,
                    step,
                    target,
                    containment,
                    _descend_spelling(cell, _child_code(row, column, step)),
                )
            )
    return out


def _descend_spelling(cell: str, code: str) -> str:
    """``cell`` with one more refinement level appended, as an atom.

    A level goes *inside* the closing brackets, not beside them: a child
    of ``NE(1414/0500(1))`` is ``NE(1414/0500(1(A1)))``, while
    ``NE(1414/0500(1)(A1))`` is a sibling pair naming two cells.
    """
    depth = len(cell) - len(cell.rstrip(DESCENT_CLOSE))
    return f"{cell[:-depth]}{DESCENT_OPEN}{code}{DESCENT_CLOSE * (depth + 1)}"


def _meridian_triangle(base: float, y_sign: float, side: float) -> "BaseGeometry":
    """Plane ring of one prime-meridian triangle, from its base and side.

    The same three points :func:`itacart.boundary._nominal_ring` builds:
    a triangle of double width centred on the line, apex poleward. Base
    is unsigned and ``y_sign`` carries the hemisphere, so one expression
    serves north and south.
    """
    from shapely.geometry import Polygon

    y = y_sign * base
    return Polygon([(-side, y), (side, y), (0.0, y + y_sign * side)])


def _meridian_centre(base: float, y_sign: float, side: float) -> "BaseGeometry":
    """Centroid of that triangle, which always lies on the meridian.

    ``x`` is zero for every triangle at every resolution, so under
    ``center`` the point falls on the line the quadrant windows cut
    along. Tested against the geometry **before** the split, where the
    line is interior and the answer is one rather than two.
    """
    from shapely.geometry import Point

    return Point(0.0, y_sign * (base + side / 3.0))


def _accept_meridian_leaf(
    prepared: "PreparedGeometry",
    containment: Containment,
    base: float,
    y_sign: float,
    side: float,
) -> bool:
    """Whether a target-resolution triangle is kept, given how it was reached.

    The counterpart of :func:`_accept_leaf`, and the chain nests for the
    same reason: the centroid is a point of the triangle, so wholly
    inside implies centre inside implies touching.
    """
    if containment == "contains":
        return False
    if containment == "intersects":
        return True
    return bool(prepared.contains(_meridian_centre(base, y_sign, side)))


def _fill_meridian_node(
    plane: "PreparedGeometry",
    lattices: dict[str, "_LatticeView"],
    hemisphere: str,
    base: float,
    y_sign: float,
    side: float,
    level: int,
    target: int,
    containment: Containment,
    code: str,
    compact: bool,
    budget: _Budget,
) -> str | None:
    """Walk one prime-meridian triangle and return its index fragment.

    The three outcomes of :func:`_fill_node`, on a triangle instead of a
    square, and with one difference that is the whole reason this
    function exists: the predicates run against the geometry **before**
    the quadrant split. A triangle straddles the line, so its two halves
    live in different clipped pieces; testing the whole figure against
    the whole geometry gives ``intersects`` the union of the two halves
    and ``contains`` their intersection, in one predicate call each,
    without a second descent to reconcile.

    Refinement follows the fold of Figure 4(c), read through
    :func:`itacart.boundary.meridian_child`: the ordinary ``d x d`` grid
    is reflected onto its own diagonal, so ``(i, j)`` and ``(j, i)`` land
    in the same sub-row on opposite sides of the line. The diagonal is
    another triangle and recurses here; everything else is an ordinary
    parallelogram whose lattice anchor is ``(|offset| * sub + base,
    base)`` in the frame of the side it fell on, so its whole subtree is
    handed to :func:`_fill_node` unchanged.

    Codes come from :func:`itacart.boundary.child_code` on the grid
    position, so the fragment is emitted in index order and no merge is
    needed afterwards.
    """
    from .boundary import child_code, meridian_child

    triangle = _meridian_triangle(base, y_sign, side)
    if not plane.intersects(triangle):
        return None
    if plane.contains(triangle):
        if compact or level == target:
            budget.charge(1)
            return code
        budget.charge(_leaves_between(level, target))
        return code + _expansion(level, target)
    if level == target:
        if _accept_meridian_leaf(plane, containment, base, y_sign, side):
            budget.charge(1)
            return code
        return None

    step = level + 1
    divisor = linear_refinement_ratio(step)
    sub = side / divisor
    parts: list[str] = []
    for row in range(divisor):
        for column in range(divisor):
            sub_row, offset = meridian_child(row, column, divisor)
            child_base = base + sub_row * sub
            text = child_code(row, column, step)
            if offset == 0:
                fragment = _fill_meridian_node(
                    plane,
                    lattices,
                    hemisphere,
                    child_base,
                    y_sign,
                    sub,
                    step,
                    target,
                    containment,
                    text,
                    compact,
                    budget,
                )
            else:
                view = lattices.get(hemisphere + ("E" if offset > 0 else "W"))
                if view is None:
                    continue
                fragment = _fill_node(
                    view.prepared,
                    abs(offset) * sub + child_base,
                    child_base,
                    sub,
                    step,
                    target,
                    containment,
                    text,
                    compact,
                    budget,
                )
            if fragment is not None:
                parts.append(fragment)
    if not parts:
        return None
    if compact and _folds_into_its_parent(parts, divisor):
        return code
    return f"{code}{DESCENT_OPEN}{SIBLING_SEPARATOR.join(parts)}{DESCENT_CLOSE}"


def _count_meridian_node(
    plane: "PreparedGeometry",
    lattices: dict[str, "_LatticeView"],
    hemisphere: str,
    base: float,
    y_sign: float,
    side: float,
    level: int,
    target: int,
) -> int:
    """Accumulate the cell count under one meridian triangle.

    The counting twin of :func:`_fill_meridian_node`, with the
    wholly-inside case answered by arithmetic, exactly as
    :func:`_count_node` answers it for a square.
    """
    from .boundary import meridian_child

    triangle = _meridian_triangle(base, y_sign, side)
    if not plane.intersects(triangle):
        return 0
    if plane.contains(triangle):
        return _leaves_between(level, target)
    if level == target:
        return int(bool(plane.contains(_meridian_centre(base, y_sign, side))))

    step = level + 1
    divisor = linear_refinement_ratio(step)
    sub = side / divisor
    total = 0
    for row in range(divisor):
        for column in range(divisor):
            sub_row, offset = meridian_child(row, column, divisor)
            child_base = base + sub_row * sub
            if offset == 0:
                total += _count_meridian_node(
                    plane,
                    lattices,
                    hemisphere,
                    child_base,
                    y_sign,
                    sub,
                    step,
                    target,
                )
                continue
            view = lattices.get(hemisphere + ("E" if offset > 0 else "W"))
            if view is None:
                continue
            total += _count_node(
                view.prepared,
                abs(offset) * sub + child_base,
                child_base,
                sub,
                step,
                target,
            )
    return total


def _meridian_rows(plane: "BaseGeometry") -> list[tuple[str, int]]:
    """Hemisphere and resolution-1 row of every meridian triangle to try.

    Candidates, not answers: a row is offered when the strip
    ``|x| <= _L1`` reaches its band of latitudes, and the descent decides
    whether the triangle is actually met. The strip is clipped first so
    that a geometry nowhere near the line offers nothing at all.

    No row past the pole is offered, and none without a column, and the
    walk does not check for either: a latitude past 90 degrees is refused
    before anything is projected, the pole projects into the polar row,
    and every eastern row up to that one addresses the meridian column.

    The polar row is offered: ``_base_cells`` skips column zero on
    purpose, because that column names the meridian triangle this walk is
    here to reach, and the cap is column zero of its row -- so if this
    walk does not offer it, nothing does. No test on
    ``last_lattice_column`` could stand in for the range: that answer
    saturates at zero, for row 1000, which holds one cell in an eastern
    quadrant, and again for row 1001, which holds none.
    """
    from shapely.geometry import box

    strip = plane.intersection(box(-_L1, -_PLANE_LIMIT, _L1, _PLANE_LIMIT))
    if strip.is_empty:
        return []
    _west, low, _east, high = strip.bounds
    out: list[tuple[str, int]] = []
    for hemisphere, y_sign in (("N", 1.0), ("S", -1.0)):
        first = low if y_sign > 0 else -high
        last = high if y_sign > 0 else -low
        # The equator belongs to the north, which is the rule
        # geo_to_cell applies and _quadrant_pieces enforces on the
        # clipped pieces. This walk is driven by the unsplit plane, so it
        # has to apply the rule itself: a geometry that only reaches y=0
        # offers a northern row and no southern one. Without this a point
        # on the equator is filled twice, once from each hemisphere.
        if last < 0.0 or (y_sign < 0.0 and last == 0.0):
            continue
        for row in range(
            max(int(math.floor(first / _L1)), 0), int(math.floor(last / _L1)) + 1
        ):
            out.append((hemisphere, row))
    return out


def _count_node(
    prepared: "PreparedGeometry",
    u0: float,
    v0: float,
    side: float,
    level: int,
    target: int,
    containment: Containment = "center",
) -> int:
    """Accumulate the cell count under one node without naming any cell.

    The same three outcomes as :func:`_fill_node`, with the wholly-inside
    case answered by arithmetic rather than by enumeration: the number of
    target cells under a node is the product of the refinement ratios
    between the two levels, so a node inside the geometry contributes its
    whole subtree in one step and is never descended.

    Memory is the recursion stack, whose depth is the resolution
    difference and therefore at most twelve. Nothing accumulates per
    cell, which is what lets the count run at resolution 13 where naming
    the cells could not.
    """
    from shapely.geometry import box

    square = box(u0, v0, u0 + side, v0 + side)
    if not prepared.intersects(square):
        return 0
    if prepared.contains(square):
        return _leaves_between(level, target)
    if level == target:
        return 1 if _accept_leaf(prepared, containment, u0, v0, side) else 0
    step = level + 1
    divisor = linear_refinement_ratio(step)
    child = side / divisor
    total = 0
    for row in range(divisor):
        for column in range(divisor):
            total += _count_node(
                prepared,
                u0 + column * child,
                v0 + row * child,
                child,
                step,
                target,
                containment,
            )
    return total


def _prepare(
    geometry: "BaseGeometry", resolution: int, densify: bool
) -> tuple["BaseGeometry", list[tuple[str, "_LatticeView"]]]:
    """Screen, densify, project and shear a geometry, one part per quadrant.

    Returns the projected plane geometry **before** the split, and the
    ``(quadrant, lattice view)`` pairs ready for descent. The unsplit
    plane is not a convenience: a prime-meridian triangle straddles the
    cut, so the only figure that can answer a containment question about
    it is the one that has not been cut.
    """

    if geometry.is_empty:
        return geometry, []
    _refuse_invalid_area(geometry)
    if crosses_antemeridian(geometry) and not _closes_at_a_pole(geometry):
        raise AntemeridianError(
            "geometry crosses 180 degrees longitude outside an extension "
            "zone; split it at the antemeridian or express it with "
            "longitudes past 180 inside a defined zone"
        )
    # Judged on what the caller supplied, before densification. A
    # parallel is not a geodesic, so densifying the northern edge of a
    # cell in the last row of a zone bulges it four metres past the
    # zone's own latitude band -- measured, 0.000036 degrees at NE row
    # 799. Refusing on that would refuse a footprint the caller wrote
    # inside the zone, for something this package did to it afterwards.
    _refuse_unzoned_extension(geometry)
    if densify and geometry.geom_type in {"Polygon", "MultiPolygon"}:
        geometry = _densify_any(geometry, _auto_segment(resolution))
    geometry = _lift_extensions(geometry)
    plane = _project(geometry)
    out: list[tuple[str, "_LatticeView"]] = []
    for quadrant, piece in _quadrant_pieces(plane):
        out.append((quadrant, _LatticeView(_to_lattice(piece, quadrant))))
    return plane, out


_T = TypeVar("_T")
_R = TypeVar("_R")


class _LatticeView:
    """A per-thread prepared copy of one quadrant's lattice geometry.

    Shapely's prepared geometry is **not** safe to share across threads.
    A prepared geometry builds its spatial index lazily on the first
    predicate call and mutates the underlying GEOS object while doing so,
    so two threads entering it together corrupt that structure and the
    interpreter dies with a segmentation fault rather than an exception.
    Found by ``test_threads_and_serial_agree``, which is why that test
    compares the two paths instead of only exercising one.

    Each thread therefore rebuilds the geometry from its WKB and prepares
    its own. Nothing is shared but the bytes, which are immutable, and
    the copy costs one deserialisation per thread rather than per cell.
    """

    def __init__(self, lattice: "BaseGeometry") -> None:
        self._wkb = lattice.wkb
        self.geometry = lattice
        self._local = threading.local()

    @property
    def prepared(self) -> "PreparedGeometry":
        """This thread's prepared geometry, built on first use."""
        from shapely import wkb as _wkb
        from shapely.prepared import prep

        cached = getattr(self._local, "prepared", None)
        if cached is None:
            cached = prep(_wkb.loads(self._wkb))
            self._local.prepared = cached
        return cached


def _map_jobs(work: Sequence[_T], call: Callable[[_T], _R], n_jobs: int) -> list[_R]:
    """Run ``call`` over ``work``, on a thread pool when asked for one.

    Threads rather than processes: the geometry predicates are Shapely's,
    which release the interpreter lock, and the work items are small
    enough that pickling a polygon to a child would dominate. Thread
    safety comes from :class:`_LatticeView` giving each worker its own
    prepared copy, not from sharing one.
    """
    if n_jobs == 1:
        return [call(item) for item in work]
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        return list(pool.map(call, work))


def polyfill(
    geometry: "BaseGeometry",
    resolution: int,
    containment: Containment = "center",
    compact: bool = False,
    n_jobs: int = 1,
) -> str:
    """Rasterise a geometry into a compositional index.

    Descends the hierarchy over the sinusoidal plane rather than testing
    cells one by one, so the cost tracks the boundary rather than the
    area.

    ``center`` keeps cells whose centre falls inside, ``intersects``
    keeps every cell touching the geometry, and ``contains`` keeps only
    cells wholly inside. The three nest, ``contains`` inside ``center``
    inside ``intersects``, and they nest by construction rather than by
    coincidence: see :func:`_accept_leaf`.

    **The centre, not the anchor.** A cell's anchor is a vertex and lies
    on its own border, so an anchor test would award a cell by a point it
    shares with three neighbours, and ``contains`` would stop being a
    subset of ``center`` for every cell whose anchor sits on the outline.
    The centre is interior to the cell and has neither problem.

    **Densification is applied, not assumed.** A straight edge on the
    plane is not a geodesic, so an undensified long edge fills the wrong
    cells in between and the result is a plausible, wrong cell set with
    no error raised. Areal input is therefore densified first, at the
    threshold :func:`_auto_segment` derives from the target resolution.
    Pass an already densified geometry and the step is idempotent.

    **The prime-meridian column.** Its cells are triangles straddling
    the line, so the ordinary descent cannot test them and a separate
    walk does. They carry the eastern spelling, which is the only one
    the grammar admits, and their containment is decided against the
    geometry before the quadrant split -- the one figure that sees both
    halves of a cell the split cuts in two.

    **Every family is filled.** The last lattice column of a row, which
    absorbs the border strip, the polar row and the two caps are not the
    sheared square the descent tests, and each is walked by a descent of
    its own rather than refused.

    **A geometry without area fills nothing under an areal predicate.**
    ``center`` and ``contains`` keep a cell by its relation to an area, so
    a point or a line keeps no cell under them, and the refusal says that
    rather than guessing the geometry was empty. ``intersects`` keeps the
    cells such a geometry touches.

    **An areal geometry must be valid.** A polygon or multipolygon that is
    not a valid geometry is refused before anything reads it, with the
    reason the geometry engine gives, and is never repaired: validity
    belongs to the geometry as written, not to where it lies. Geometries
    without area are not asked.

    Args:
        geometry: A Shapely geometry in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        containment: Predicate deciding whether a cell is kept.
        compact: Fold every node whose children were all kept, giving the
            same index as compacting the uniform fill afterwards with
            :func:`itacart.hierarchy.compact_cells` under the same
            containment mode. Returns a mixed-resolution index instead of a
            uniform one. Refused where the fill walks the border-absorbing
            family -- where the geometry meets the strip the last lattice
            column of a row absorbs, or a polar cap -- because the fold
            reads a uniform tree and that family has none.
        n_jobs: Worker count; above 1 spreads base cells over threads.

    Returns:
        A compositional index string covering the geometry.

    Raises:
        AntemeridianError: If the geometry crosses 180 degrees outside an
            extension zone.
        UnsupportedGeometryTypeError: On unsupported geometry types.
        DomainError: If a position lies outside the latitude domain.
        GeometryError: If a polygon or multipolygon is not a valid
            geometry, if the fill exceeds :data:`MAX_FILL_CELLS`, if
            ``compact=True`` meets the border-absorbing family, or if it
            keeps no cell: the geometry is empty, it has no area under ``center``
            or ``contains``, or it is narrower than one cell under the
            chosen mode. The message names which.
    """
    _check_resolution(resolution)
    _check_jobs(n_jobs)
    if containment not in ("center", "intersects", "contains"):
        raise ValueError(
            f"containment must be 'center', 'intersects' or 'contains', "
            f"got {containment!r}"
        )

    budget = _Budget()
    plane, views = _prepare(geometry, resolution, densify=True)
    lattices = dict(views)
    ordered: dict[str, list[str]] = {}
    for quadrant, view in views:

        def _one(
            base: tuple[int, int, bool],
            _v: "_LatticeView" = view,
            _q: str = quadrant,
        ) -> str | None:
            column, row, border = base
            code = f"{column:0{RES1_DIGITS}d}{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}"
            if border:
                return _fill_border_root(
                    _v.prepared,
                    _q,
                    column,
                    row,
                    resolution,
                    containment,
                    compact,
                    budget,
                )
            return _fill_node(
                _v.prepared,
                (column + row) * _L1,
                row * _L1,
                _L1,
                1,
                resolution,
                containment,
                code,
                compact,
                budget,
            )

        # Sorted by column then row so that the siblings come out in the
        # order the index is read in. compose() preserves the order it is
        # given rather than imposing one, so the ordering has to be here.
        cells = sorted(_base_cells(view.prepared, view.geometry, quadrant, resolution))
        parts = [f for f in _map_jobs(cells, _one, n_jobs) if f is not None]
        if parts:
            ordered.setdefault(quadrant, []).extend(parts)

    # The prime-meridian column, walked once per hemisphere against the
    # unsplit plane and spelled from the east, which is the only side that
    # names it. Prepended rather than appended: column 0 sorts before every
    # ordinary column, and the eastern root may not exist yet when the
    # geometry lies wholly west of the line.
    meridian = _meridian_rows(plane)
    if meridian:
        from shapely.prepared import prep

        prepared_plane = prep(plane)
        seam: dict[str, list[str]] = {}
        for hemisphere, row in meridian:
            code = f"{0:0{RES1_DIGITS}d}{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}"
            if row == _POLAR_ROW:
                # The cap absorbs the domain border, so the uniform walk
                # below does not describe it: measured, that walk named
                # two spellings ``is_valid_cell`` rejects and missed two
                # it should have named. It is handed to the walker built
                # for that family instead.
                fragment = _fill_border_root(
                    prepared_plane,
                    hemisphere + "E",
                    0,
                    row,
                    resolution,
                    containment,
                    compact,
                    budget,
                    sheared=False,
                )
                if fragment is not None:
                    seam.setdefault(hemisphere + "E", []).append(fragment)
                continue
            fragment = _fill_meridian_node(
                prepared_plane,
                lattices,
                hemisphere,
                row * _L1,
                1.0 if hemisphere == "N" else -1.0,
                _L1,
                1,
                resolution,
                containment,
                code,
                compact,
                budget,
            )
            if fragment is not None:
                seam.setdefault(hemisphere + "E", []).append(fragment)
        for quadrant, fragments in seam.items():
            ordered[quadrant] = fragments + ordered.get(quadrant, [])

    roots = [
        f"{quadrant}{DESCENT_OPEN}"
        f"{SIBLING_SEPARATOR.join(ordered[quadrant])}{DESCENT_CLOSE}"
        for quadrant in QUADRANTS
        if ordered.get(quadrant)
    ]
    if not roots:
        raise GeometryError(_covers_no_cell(geometry, containment))
    return SIBLING_SEPARATOR.join(roots)


def _covers_no_cell(geometry: "BaseGeometry", containment: str) -> str:
    """Name why a fill kept no cell, from the geometry rather than a guess."""
    if geometry.is_empty:
        cause = "it is empty"
    elif geometry.area == 0.0 and containment != "intersects":
        cause = (
            f"a {geometry.geom_type} has no area, and {containment!r} keeps a "
            "cell only by its relation to an area; 'intersects' keeps the "
            "cells the geometry touches, and a point has exactly one owner, "
            "which geo_to_cell answers"
        )
    else:
        cause = f"it is narrower than one cell under {containment!r}"
    return f"geometry covers no cell at this resolution: {cause}"


def count_internal_cells(polygon: "Polygon", resolution: int, n_jobs: int = 1) -> int:
    """Count the cells a polygon covers, without naming any of them.

    Fast path over :func:`polyfill` with ``containment="center"``: the
    count is accumulated during descent, so nothing is materialised and
    the call runs at resolution 13 where the index could not be built.

    Since every cell carries the same area, the count times
    :func:`itacart.resolutions.nominal_cell_area` gives a
    distortion-free area, which is the property the paper builds
    tokenization on.

    **Why the centre and not strict containment.** Counting only the
    cells wholly inside undercounts by the ring of cells the outline
    crosses, about ``perimeter * side / 2`` in area terms. For a
    square kilometre at resolution 7 that is two per cent. Under the
    centre rule the cells the outline crosses are kept or dropped
    according to which side their centre falls, and the two errors
    cancel to first order. The origin counts the same way.

    **What the residual error is, and what it is not.** The count is
    exact for every cell the outline does not cross, so the area it
    implies can only be wrong about the crossed ones:

    ``abs(n * a - A) / A <= P * s / A`` for perimeter ``P``, cell side
    ``s`` and area ``A``. That bound is arithmetic and holds for every
    parcel.

    Where inside the bound a given parcel lands is not predictable.
    Measured over 24 phases at each of seven parcel sizes, three
    estimators of the scaling exponent disagree outright: 0.61 from the
    RMS, -0.47 from the median, 0.00 from the ninetieth percentile.
    There is no power law to fit, because the residual depends on how
    that outline happens to sit against that lattice.

    This matters for how the resolution is chosen. It is not chosen from
    the parcel: in cadastral work it follows from the mapping scale and
    the positional tolerance the applicable standard attaches to it, and
    it is fixed before any parcel is measured. The parcel's size and
    perimeter then determine the residual, which is an outcome to be
    reported against the bound above, not a target to refine toward.

    **Every family is counted.** The last lattice column of a row, the
    polar row and the two caps are not the sheared square the closed form
    counts, so their trees are walked -- ring by ring inside the family,
    in closed form wherever a branch leaves it -- and the prime-meridian
    column is counted by a walk of its own. None of those walks stops at
    :data:`MAX_FILL_CELLS`. That ceiling keeps a fill from naming more
    cells than can be materialised, and a count names none: stopping there
    would return the running total at the moment of passing it as though
    it were the count.

    Provenance: ``itacart_core/cell_filling.py``
    (``polygon_to_cells_count``).

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        n_jobs: Worker count; above 1 spreads base cells over threads.

    Returns:
        Number of cells whose centre lies inside the polygon.

    Raises:
        AntemeridianError: If the polygon crosses 180 degrees outside an
            extension zone.
        DomainError: If a position lies outside the latitude domain.
        GeometryError: If the polygon is not a valid geometry; it is refused,
            not repaired.
    """
    _check_resolution(resolution)
    _check_jobs(n_jobs)
    total = 0
    plane, views = _prepare(polygon, resolution, densify=True)
    lattices = dict(views)
    for quadrant, view in views:

        def _one(
            base: tuple[int, int, bool],
            _v: "_LatticeView" = view,
            _q: str = quadrant,
        ) -> int:
            column, row, border = base
            if border:
                walk = _BorderWalk(_q)
                root = f"{_q}{DESCENT_OPEN}{column:0{RES1_DIGITS}d}"
                root += f"{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}{DESCENT_CLOSE}"
                return _count_border_node(
                    _v.prepared,
                    walk,
                    root,
                    1,
                    resolution,
                    "center",
                    None,
                )
            return _count_node(
                _v.prepared, (column + row) * _L1, row * _L1, _L1, 1, resolution
            )

        cells = _base_cells(view.prepared, view.geometry, quadrant, resolution)
        total += sum(_map_jobs(cells, _one, n_jobs))

    meridian = _meridian_rows(plane)
    if meridian:
        from shapely.prepared import prep

        prepared_plane = prep(plane)
        for hemisphere, row in meridian:
            if row == _POLAR_ROW:
                walk = _BorderWalk(hemisphere + "E", sheared=False)
                root = (
                    f"{hemisphere}E{DESCENT_OPEN}{0:0{RES1_DIGITS}d}"
                    f"{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}{DESCENT_CLOSE}"
                )
                total += _count_border_node(
                    prepared_plane,
                    walk,
                    root,
                    1,
                    resolution,
                    "center",
                    None,
                )
                continue
            total += _count_meridian_node(
                prepared_plane,
                lattices,
                hemisphere,
                row * _L1,
                1.0 if hemisphere == "N" else -1.0,
                _L1,
                1,
                resolution,
            )
    return total


# --------------------------------------------------------------------------
# Vertex representation
# --------------------------------------------------------------------------


def _dedupe_consecutive(cells: Sequence[str], cyclic: bool) -> list[str]:
    """Collapse runs of the same index, optionally treating the list as a ring.

    Only runs. A repeat that is not adjacent is left alone: at the target
    resolution two consecutive vertices in the same cell are the same
    surveyed point written twice, while a vertex revisited later in the
    ring is a self-touching outline, which is a genuine pathology and the
    caller has to be able to see it.

    With ``cyclic`` the wrap-around pair is treated as adjacent too,
    since a ring's last vertex is followed by its first.
    """
    if not cells:
        return []
    out = [cells[0]]
    for cell in cells[1:]:
        if cell != out[-1]:
            out.append(cell)
    if cyclic and len(out) >= 2 and out[0] == out[-1]:
        out.pop()
    return out


def _ring_cells(
    coords: Sequence[tuple[float, ...]], resolution: int, dedupe: bool
) -> list[str]:
    """Map one closed ring to cells, dropping the repeated closing vertex."""
    points = list(coords)
    if len(points) >= 2 and points[0] == points[-1]:
        points = points[:-1]
    cells = [geo_to_cell(point[0], point[1], resolution) for point in points]
    return _dedupe_consecutive(cells, cyclic=True) if dedupe else cells


def vertex_to_cell(
    geometry: "Point | LineString | Polygon",
    resolution: int,
    dedupe_consecutive: bool = True,
) -> list[str]:
    """Map geometry vertices to cells, preserving sequence.

    Order is topological, not sorted: it is what allows the original
    geometry to be reconstructed. Rings keep their winding, and holes
    follow the exterior.

    Consecutive vertices landing in the same cell are collapsed by
    default, since at resolution 13 that reflects survey precision rather
    than distinct corners. Non-consecutive repeats are kept, being a
    genuine self-touching pathology the caller should see.

    A ring is treated as cyclic for the purposes of collapsing, so a ring
    whose last vertex falls in the same cell as its first loses the last
    one rather than closing on a duplicate. A LINESTRING is not cyclic
    and keeps both ends even when they coincide.

    **The cell is the precision.** A vertex written as a cell says how
    well it is known -- within its cell at the resolution chosen -- which
    a coordinate cannot say: section 3.1.10 of RFC 7946 rules out reading
    uncertainty from the number of digits. :func:`cells_to_geometry` is the
    inverse, and its docstring states what a rebuilt vertex promises.

    Provenance: ``cadastral_processor/vertex_extractor.py``.

    Args:
        geometry: A Shapely geometry in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        dedupe_consecutive: Collapse consecutive duplicates.

    Returns:
        Atomic index strings in traversal order.

    Raises:
        UnsupportedGeometryTypeError: On unsupported geometry types.
    """
    from shapely.geometry import LineString, Point, Polygon

    _check_resolution(resolution)
    if isinstance(geometry, Point):
        return [geo_to_cell(geometry.x, geometry.y, resolution)]
    if isinstance(geometry, Polygon):
        out = _ring_cells(
            list(geometry.exterior.coords), resolution, dedupe_consecutive
        )
        for hole in geometry.interiors:
            out.extend(_ring_cells(list(hole.coords), resolution, dedupe_consecutive))
        return out
    if isinstance(geometry, LineString):
        cells = [geo_to_cell(x, y, resolution) for x, y in geometry.coords]
        return _dedupe_consecutive(cells, cyclic=False) if dedupe_consecutive else cells
    raise UnsupportedGeometryTypeError(
        f"cannot take vertices of a {geometry.geom_type}; supported types "
        "are Point, LineString and Polygon"
    )


def cells_to_geometry(
    cells: list[str], geometry_type: str = "Polygon"
) -> "BaseGeometry":
    """Rebuild a geometry from an ordered vertex cell list.

    Inverse of :func:`vertex_to_cell`. **Every vertex is rebuilt at the
    representative position of its cell, the centroid**, which is the
    point :func:`itacart.cells.cell_to_centroid` returns and
    :func:`itacart.engine.describe` declares. There is no argument
    choosing another: the anchor is a vertex of the cell, on its boundary
    rather than within it, and it does not answer the question this
    function asks.

    **What a rebuilt vertex promises.** It lies within the cell that
    encoded the original vertex and quantizes back to that cell at the
    same resolution, in every family, extension zones and polar caps
    included. Its distance from the original is at most the cell's
    *reach*: the distance from the centroid to the cell's farthest vertex.
    For an ordinary parallelogram of side ``l`` the reach follows the
    shear of the projection::

        l * max(sqrt(1 + u**2) / 2, sqrt((1 - u / 2)**2 + 1 / 4))
        u = |lon * sin(lat)|, with lon in radians, at the centroid

    which is 1.118 l on the prime meridian, 0.707 l at its minimum where
    ``u = 1``, and 1.65 l as ``u`` reaches pi. That is the reach to first
    order in the side, and :func:`itacart.metrics.cell_max_uncertainty`
    returns the reach of any cell. A finite cell departs from the law in
    proportion to its side, and most near the pole: over the population
    that function documents, which runs through the last rows before the
    polar row and past ``u = pi`` in Chukotka, the reach of a regular
    parallelogram exceeds the law by at most ``4.2e-3 * l / 10 km`` below
    85 degrees of latitude and ``0.2 * l / 10 km`` from 85 degrees to row
    998, plus up to ``1e-6`` from the rounding of the coordinates.

    Three families are not bounded by that law, and are named rather than
    folded into it. A prime-meridian triangle reaches ``sqrt(10) / 3 * l``
    to first order, and a finite one falls short of it: by at most
    ``2.7e-3 * l / 10 km`` below 85 degrees and ``0.13 * l / 10 km`` from
    85 degrees to row 998, over the triangles of rows 930 to 998 and down
    the prime meridian to resolution 13. A cell that absorbs the domain
    border, a polar-row cell and a cap cell reach further: up to 2.06 l
    over the last existing column of every resolution-1 row in the four
    quadrants, together with the border, polar-row and cap cells of a
    lattice of the four quadrants at five resolutions from 10 km to 1 cm and
    both caps at every resolution. That is a measurement over the population
    named, not a bound beyond it.

    **Precision.** Coordinates are returned at full double precision and
    never rounded. How precisely a rebuilt vertex is known is declared by
    the resolution of its cell, not by its digits: section 3.1.10 of RFC
    7946 says the number of digits of a coordinate does not indicate its
    uncertainty, so a coordinate cannot carry it and an index can. Rounding
    would also destroy what it claims to summarise: resolution-13
    neighbours lie about 9e-8 degree of latitude apart, and at six decimal
    places their centroids become one point.

    **The inverse is exact for a hole-free geometry only.**
    :func:`vertex_to_cell` returns one flat sequence with the holes
    appended to the exterior, and a flat sequence does not say where the
    exterior ended. A caller needing holes keeps the rings apart and
    calls this once per ring.

    Args:
        cells: Atomic index strings in traversal order.
        geometry_type: OGC SFA type to build: ``Point``, ``LineString``,
            ``LinearRing``, ``Polygon`` or ``MultiPoint``.

    Returns:
        A Shapely geometry in EPSG:4326.

    Raises:
        UnsupportedGeometryTypeError: On unsupported types.
        GeometryError: If the cell count cannot make the type asked for.
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    from shapely.geometry import LinearRing, LineString, MultiPoint, Point, Polygon

    if not cells:
        raise GeometryError("cannot rebuild a geometry from an empty cell list")
    coords: list[tuple[float, float]] = []
    for cell in cells:
        position = cell_to_centroid(cell)
        if not isinstance(position, tuple):
            raise GeometryError(
                f"{cell!r} names more than one cell; cells_to_geometry takes "
                "atomic indices in traversal order"
            )
        coords.append((float(position[0]), float(position[1])))

    if geometry_type == "Point":
        if len(coords) != 1:
            raise GeometryError(f"a Point needs exactly one cell, got {len(coords)}")
        return Point(coords[0])
    if geometry_type == "MultiPoint":
        return MultiPoint(coords)
    if geometry_type == "LineString":
        if len(coords) < 2:
            raise GeometryError(
                f"a LineString needs at least two cells, got {len(coords)}"
            )
        return LineString(coords)
    if geometry_type in ("Polygon", "LinearRing"):
        if len(coords) < 3:
            raise GeometryError(
                f"a {geometry_type} needs at least three cells, got {len(coords)}"
            )
        ring = coords + [coords[0]]
        return LinearRing(ring) if geometry_type == "LinearRing" else Polygon(ring)
    raise UnsupportedGeometryTypeError(
        f"cannot rebuild a {geometry_type!r}; supported types are Point, "
        "MultiPoint, LineString, LinearRing and Polygon"
    )


# --------------------------------------------------------------------------
# Densification
# --------------------------------------------------------------------------


EDGE_MODELS: tuple[str, ...] = ("WGS84_GEODESIC",)
"""Edge interpretations :func:`densify_segment` understands."""


def _check_threshold(max_segment_m: float) -> None:
    """Reject a threshold that is not a positive finite length."""
    if not isinstance(max_segment_m, (int, float)) or isinstance(max_segment_m, bool):
        raise DensificationError(
            f"max_segment_m must be a number, got {type(max_segment_m).__name__}"
        )
    if not math.isfinite(max_segment_m) or max_segment_m <= 0.0:
        raise DensificationError(
            f"max_segment_m must be a positive finite length in metres, "
            f"got {max_segment_m!r}"
        )


def densify_segment(
    p1: tuple[float, float],
    p2: tuple[float, float],
    max_segment_m: float,
    edge_model: str = "WGS84_GEODESIC",
) -> list[tuple[float, float]]:
    """Densify a single segment.

    The building block of :func:`densify_orthodromic`, exposed for
    callers working segment by segment such as open LINESTRING handling.

    The span is measured with the Vincenty inverse and the intermediate
    points are placed with the Vincenty direct along the forward azimuth,
    equally spaced in geodesic arc length. Nothing is measured on the
    plane, because the projection preserves area and not direction.

    ``n_segments = floor(d / max_segment_m) + 1``, so every resulting leg
    is strictly shorter than the threshold whenever any subdivision
    happens at all. That is what makes the operation idempotent: a second
    pass measures each leg, finds it under the threshold, and leaves it.

    **Each interior point is put on the branch of the endpoint nearer to
    it in arc length**, the start for the first half and the destination
    for the second. On an ordinary segment the two references sit less
    than half a turn apart and agree, so the rule changes nothing. On a
    segment spanning exactly half a turn of longitude they do not agree,
    because such a segment is meridional: its geodesic passes through a
    pole and its longitude jumps a hundred and eighty degrees there.
    Referring the far half to the start would let it land on the opposite
    branch from the vertex it is walking towards, which folds the ring
    across the globe and leaves the projection with a self-intersection
    the caller cannot read. The midpoint itself, which is the pole when
    the segment is symmetric about it, is referred to the start; its
    longitude is immaterial because the parallels plane collapses every
    meridian at the pole onto one point.

    Provenance: ``itacart_core/geometry_blob.py`` (``densify_segment``).

    Args:
        p1: ``(lon, lat)`` of the start point.
        p2: ``(lon, lat)`` of the end point.
        max_segment_m: Longest segment to leave undensified, in metres.
        edge_model: Edge interpretation, currently ``"WGS84_GEODESIC"``.

    Returns:
        ``(lon, lat)`` pairs including both endpoints.

    Raises:
        DensificationError: If the threshold is not positive and finite,
            or the edge model is unknown.
    """
    _check_threshold(max_segment_m)
    if edge_model not in EDGE_MODELS:
        raise DensificationError(
            f"unknown edge model {edge_model!r}; supported models are "
            f"{', '.join(EDGE_MODELS)}"
        )
    lon1, lat1 = float(p1[0]), float(p1[1])
    lon2, lat2 = float(p2[0]), float(p2[1])
    distance, azimuth = inverse_geodesic(lon1, lat1, lon2, lat2)
    if distance <= 0.0:
        return [(lon1, lat1), (lon2, lat2)]
    pieces = int(math.floor(distance / max_segment_m)) + 1
    if pieces == 1:
        return [(lon1, lat1), (lon2, lat2)]
    step = distance / pieces
    out = [(lon1, lat1)]
    for index in range(1, pieces):
        longitude, latitude = direct_geodesic(lon1, lat1, azimuth, step * index)
        reference = lon1 if index * 2 <= pieces else lon2
        out.append((_on_the_branch_of(reference, longitude), latitude))
    out.append((lon2, lat2))
    return out


def _on_the_branch_of(reference: float, longitude: float) -> float:
    """Move ``longitude`` onto the 360-degree branch of ``reference``.

    :func:`itacart.geodesy.direct_geodesic` normalises what it returns to
    ``(-180, 180]``, which is the right answer to the question it is
    asked and the wrong one for densification. A segment written from
    179.9 to 180.3 -- the natural way to describe a footprint inside an
    extension zone, where the domain reaches past the antemeridian --
    has interior points at 180.1, and normalising those to -179.9 folds
    the ring back across the globe. The result self-intersects, and the
    quadrant clip downstream fails inside GEOS with a topology error
    rather than anywhere the package can explain.

    Densification therefore keeps the branch the caller wrote. Segments
    that genuinely wrap are refused upstream by
    :func:`itacart.boundary.crosses_antemeridian`, so no segment reaching
    here spans **more** than half the globe.

    At exactly half a globe the branch is not unambiguous, and this
    function cannot make it so from one reference. Two points half a turn
    of longitude apart lie on a single meridian circle: the geodesic
    joining them runs through the pole, and the longitude changes by a
    hundred and eighty degrees at the crossing. Both ``+180`` and
    ``-180`` name that change and the quotient is exactly one half, so
    the rounding decides, and it decides without knowing which of the two
    the segment is walking towards. :func:`densify_segment` therefore
    supplies the destination as the reference for the far half of the
    walk rather than asking this function to guess.
    """
    return longitude - 360.0 * round((longitude - reference) / 360.0)


def _densify_ring(
    coords: Sequence[tuple[float, ...]], max_segment_m: float
) -> list[tuple[float, float]]:
    """Densify one ring, returned open with no repeated closing vertex."""
    points = [(float(point[0]), float(point[1])) for point in coords]
    if len(points) >= 2 and points[0] == points[-1]:
        points = points[:-1]
    if len(points) < 2:
        return points
    out: list[tuple[float, float]] = []
    for index, start in enumerate(points):
        end = points[(index + 1) % len(points)]
        out.extend(densify_segment(start, end, max_segment_m)[:-1])
    return out


def densify_orthodromic(polygon: "Polygon", max_segment_m: float = 1000.0) -> "Polygon":
    """Insert intermediate vertices along geodesics.

    Applied to the exterior ring and every hole, with
    ``n_segments = floor(d_geo / max_segment_m) + 1`` and points equally
    spaced in geodesic distance.

    Needed because a straight line on the sinusoidal plane is not a
    geodesic on the ellipsoid; without densification a long edge would
    fill the wrong cells in between.

    Holes are densified too. The origin drops them, which is silent: the
    filled cell set comes back plausible and includes the hole.

    Provenance: ``itacart_core/densification.py``.

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        max_segment_m: Longest segment to leave undensified, in metres.

    Returns:
        A densified polygon in EPSG:4326.

    Raises:
        DensificationError: If a segment cannot be densified.
        TypeError: If the argument is not a Shapely polygon.
    """
    from shapely.geometry import Polygon

    if not isinstance(polygon, Polygon):
        raise TypeError(f"expected a shapely Polygon, got {type(polygon).__name__}")
    _check_threshold(max_segment_m)
    if polygon.is_empty:
        return polygon
    return Polygon(
        _densify_ring(list(polygon.exterior.coords), max_segment_m),
        [_densify_ring(list(hole.coords), max_segment_m) for hole in polygon.interiors],
    )


def _densify_any(geometry: "BaseGeometry", max_segment_m: float) -> "BaseGeometry":
    """Densify a polygon or a multipolygon, leaving anything else alone."""
    from shapely.geometry import MultiPolygon, Polygon

    if isinstance(geometry, Polygon):
        return densify_orthodromic(geometry, max_segment_m)
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon(
            [densify_orthodromic(part, max_segment_m) for part in geometry.geoms]
        )
    return geometry


# --------------------------------------------------------------------------
# Canonicalization
# --------------------------------------------------------------------------


def _min_rotation(ring: Sequence[str]) -> int:
    """Start index of the lexicographically smallest rotation of ``ring``.

    Booth's algorithm, linear in the ring length. A naive scan over the
    smallest vertex breaks on a ring that repeats its minimum, which a
    self-touching outline does; Booth compares whole rotations and orders
    those cases too.

    Provenance: ``itacart_core/geometry_blob.py``
    (``_booth_min_rotation``).
    """
    length = len(ring)
    if length == 0:
        return 0
    doubled = list(ring) * 2
    failure = [-1] * (2 * length)
    start = 0
    for position in range(1, 2 * length):
        candidate = doubled[position]
        offset = failure[position - start - 1]
        while offset != -1 and candidate != doubled[start + offset + 1]:
            if candidate < doubled[start + offset + 1]:
                start = position - offset - 1
            offset = failure[offset]
        if candidate != doubled[start + offset + 1]:
            if candidate < doubled[start + offset + 1]:
                start = position
            failure[position - start] = -1
        else:
            failure[position - start] = offset + 1
    return start


def canonicalize_rings(rings: list[list[str]]) -> list[list[str]]:
    """Normalise rings to a single spelling per geometry.

    Rings are rotated to their minimum lexicographic cyclic rotation, so
    the same ring starting at a different vertex canonicalises to the
    same sequence. Rotation applies to closed rings only; LINESTRING is
    directional and is left as given, which is why this function takes
    rings and not lines.

    Winding is **not** touched. A ring's direction distinguishes an
    exterior from a hole and is information the caller put there;
    rotating is a change of spelling, reversing would be a change of
    meaning. Ordering the rings themselves is likewise the caller's,
    since the exterior comes first by contract.

    The ordering key is the index string. Two cells of the same
    resolution have strings of the same length, so string order is a
    total order over the ring and the rotation is unique.

    This is what makes a geometry content-addressable, and therefore what
    makes hashing it meaningful.

    Provenance: ``itacart_core/geometry_blob.py``
    (``_canonicalize_polygon_rings``, ``_booth_min_rotation``).

    Args:
        rings: Rings as lists of atomic index strings, exterior first.

    Returns:
        The canonicalised rings, in input order.
    """
    out: list[list[str]] = []
    for ring in rings:
        if len(ring) < 2:
            out.append(list(ring))
            continue
        start = _min_rotation(ring)
        out.append(list(ring[start:]) + list(ring[:start]))
    return out
