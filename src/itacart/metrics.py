"""Cell shape metrics, for comparison against other DGGS implementations.

The two metrics of :func:`compactness` and :func:`normalized_cell_area`
are transcribed from Kmoch, Vasilyev, Virro and Uuemaa (2022), "Area and
shape distortions in open-source discrete global grid systems", *Big
Earth Data* 6(3), 256-275, doi:10.1080/20964471.2022.2094926. That paper
measures H3, S2, OpenEAGGR, rHEALPix and DGGRID on exactly these two
numbers, and the point of reproducing its conventions rather than better
ones is that the result has to be comparable with its figures. Where a
geodesic quantity would have been more accurate, the paper's planar
convention is kept and the difference is stated.

:func:`cell_base_angle` is not from that paper. It answers the question
the ITACaRT paper leaves qualitative: the cells are identical
parallelograms with a 45-degree acute angle on the sinusoidal plane, and
Figure 8 of Silva, Dietzsch and Shiguemori (2025) shows in pictures what
becomes of that angle on the ellipsoid. This function is that picture as
a number.

Three properties of the results are worth knowing before reading them.
Compactness is not monotone in latitude: it rises from 0.539 at the
equator to 0.785 at row 900, where the shear of the index cancels the
convergence of the meridians and the cell is a square to five decimal
places, and then falls again. The angle is a function of the lattice
position, not of latitude alone, and it is not monotone either. And
compactness is undefined on a cell whose ring wraps the whole parallel --
the cell that holds a pole, at every resolution, and from resolution 4
down many cells around it inside the polar caps -- for a reason documented
in :func:`compactness`; the angle is defined there, and is zero on the cell
that holds the pole.

Both metrics hold at every resolution, down to the one-centimetre cell of
resolution 13 and the thinnest trapezoid the border cuts. Nothing in them is
a tolerance in degrees: the cells span six orders of magnitude, and a
tolerance that suited the coarse ones once swallowed every cell from
resolution 11 down.

:func:`cell_max_uncertainty` is not a shape metric. It says how far a
position written as a cell can lie from the point the cell stands for, its
centroid, measured on the ground: the reach that
:func:`itacart.geometry.cells_to_geometry` states for a rebuilt vertex.
"""

from __future__ import annotations

import math
from typing import Literal, cast

from ._existence import require_existing_cells
from .boundary import _side_quadrant
from .cells import cell_to_anchor, cell_to_boundary, cell_to_centroid
from .constants import WGS84_A, WGS84_E2
from .exceptions import GeometryError
from .geodesy import geodetic_to_sinusoidal, inverse_geodesic
from .index import is_atomic, iter_cells, quadrant_of
from .resolutions import (
    cell_size,
    effective_cell_area,
    get_resolution,
    nominal_cell_area,
)

__all__ = [
    "compactness",
    "cell_base_angle",
    "normalized_cell_area",
    "cell_max_uncertainty",
]

AngleOutput = Literal["degrees", "radians", "sin"]
UncertaintyOutput = Literal["radius", "diameter", "ratio"]

_E = math.sqrt(WGS84_E2)

#: How far apart two longitudes may be and still count as a full turn, in
#: units in the last place of 360. This is the rounding of a longitude and
#: not a distance on the ground, so it does not shrink with the cell: the
#: tolerance this replaced was a micro-degree, about 0.11 m, and every cell
#: from resolution 11 down is smaller than that.
_FULL_TURN_ULPS = 4


# --------------------------------------------------------------------------
# Lambert Azimuthal Equal Area, oblique aspect, on the ellipsoid
# --------------------------------------------------------------------------


def _authalic_q(sin_phi: float) -> float:
    """Snyder's ``q``, the authalic area function.

    Equation (3-12) of Snyder, *Map Projections - A Working Manual*, USGS
    Professional Paper 1395 (1987), page 16. ``q`` is proportional to the
    area of the zone between the equator and the given latitude, which is
    what makes the projection built on it equal-area.
    """
    t = _E * sin_phi
    return (1.0 - WGS84_E2) * (
        sin_phi / (1.0 - t * t) - (1.0 / (2.0 * _E)) * math.log((1.0 - t) / (1.0 + t))
    )


_QP = _authalic_q(1.0)

#: Radius of the sphere of equal surface area, Snyder equation (3-13).
_RQ = WGS84_A * math.sqrt(_QP / 2.0)


def _sin_cos_deg(lat_deg: float) -> tuple[float, float]:
    """Sine and cosine of a latitude given in degrees, both to full precision.

    Near a pole the cosine is small, and ``cos(radians(lat))`` computes it
    from an argument already rounded to the precision of a number near
    ``pi / 2``: a relative error of 1e-10 at 89.9999 degrees. Taking it as
    the sine of the colatitude, which is exact in degrees, keeps every
    digit.
    """
    if abs(lat_deg) <= 45.0:
        phi = math.radians(lat_deg)
        return math.sin(phi), math.cos(phi)
    colatitude = math.radians(90.0 - abs(lat_deg))
    return math.copysign(math.cos(colatitude), lat_deg), math.sin(colatitude)


def _authalic_sin_cos(lat_deg: float) -> tuple[float, float]:
    """Sine and cosine of the authalic latitude, Snyder equation (3-11).

    ``sin(beta) = q / q_p``, and the authalic latitude was taken as the
    arcsine of that ratio. Near a pole the ratio is within a few parts in
    a trillion of one, the arcsine is ill-conditioned there, and the
    rounding of the ratio alone moved the projected point by about 0.4 mm
    at 89.9999 degrees -- which, for a one-centimetre cell, is most of the
    cell. The cosine is computed instead from ``q_p - q``, the area between
    the latitude and the pole, in a form with no cancellation:

        q_p - q = (1 - s)(1 + e2 s) / (1 - e2 s**2)
                  + (1 - e2) atanh(e (1 - s) / (1 - e2 s)) / e

    with ``s = |sin(phi)|`` and ``1 - s`` taken from the colatitude. Both
    terms carry the factor ``1 - s`` explicitly, so the difference is never
    formed by subtracting two nearly equal numbers. The sine is computed at
    ``s`` and given the latitude's sign, so south mirrors north to the bit.
    """
    sin_phi, cos_phi = _sin_cos_deg(lat_deg)
    s = abs(sin_phi)
    if abs(lat_deg) <= 45.0:
        one_minus_s = 1.0 - s
    else:
        half = math.radians(90.0 - abs(lat_deg)) / 2.0
        one_minus_s = 2.0 * math.sin(half) ** 2
    to_pole = (
        one_minus_s * (1.0 + WGS84_E2 * s) / (1.0 - WGS84_E2 * s * s)
        + (1.0 - WGS84_E2) * math.atanh(_E * one_minus_s / (1.0 - WGS84_E2 * s)) / _E
    )
    sin_beta = math.copysign(_authalic_q(s) / _QP, lat_deg)
    cos_beta = math.sqrt(to_pole * (2.0 * _QP - to_pole)) / _QP
    return sin_beta, cos_beta


def _laea(
    lon_deg: float, lat_deg: float, lon0_deg: float, lat0_deg: float
) -> tuple[float, float]:
    """Project one point into the oblique ellipsoidal LAEA of an origin.

    Snyder equations (24-13) through (24-20), pages 187-190. The origin
    is the cell's own centroid, which is what step 3 of section 2.2 of
    Kmoch et al. (2022) prescribes: each cell is reprojected individually
    with its centroid as the origin of the projection.

    The projection is equal-area, so the planar area of a cell's ring
    agrees with :func:`itacart.resolutions.effective_cell_area` to within
    the error of drawing the cell's sides as straight lines in the plane.
    That agreement is asserted in the test suite, and it is the control
    that says this transcription is the projection it claims to be.
    """
    lam = math.radians(lon_deg - lon0_deg)
    sin_beta, cos_beta = _authalic_sin_cos(lat_deg)
    sin_beta0, cos_beta0 = _authalic_sin_cos(lat0_deg)
    sin_phi0, cos_phi0 = _sin_cos_deg(lat0_deg)
    d = (WGS84_A * cos_phi0) / (
        math.sqrt(1.0 - WGS84_E2 * sin_phi0**2) * _RQ * cos_beta0
    )
    denominator = 1.0 + sin_beta0 * sin_beta + cos_beta0 * cos_beta * math.cos(lam)
    b = _RQ * math.sqrt(2.0 / denominator)
    x = b * d * cos_beta * math.sin(lam)
    y = (b / d) * (cos_beta0 * sin_beta - sin_beta0 * cos_beta * math.cos(lam))
    return x, y


# --------------------------------------------------------------------------
# Ring predicates and planar measures
# --------------------------------------------------------------------------


def _repeated_vertex(ring: list[tuple[float, float]]) -> tuple[int, int] | None:
    """Return the first pair of ring positions holding the same point.

    A ring that visits a point twice is not a simple closed curve, and
    the area enclosed by it is not the area of the cell. The cell that
    holds a pole is in exactly that state, at every resolution -- the polar
    cap and, below it, the one triangle per resolution that still spans the
    parallel: its ring runs pole, -180, +180 at one latitude, the two ends
    of the base are the same point, and the limiting parallel is missing
    from the boundary altogether. Some cells around it meet the same
    parallel at -180 and +180 as well, and many more, inside the caps from
    resolution 4 down, wrap the parallel several times over without
    repeating a vertex at all; :func:`_longitude_span` is what refuses
    those.

    The test is the construction, not a distance. Two vertices hold the
    same point here when they lie on the same parallel -- the same
    latitude, bit for bit, because both come from one ordinate -- with
    longitudes a full turn apart. The turn is measured to the rounding of
    a longitude, a few units in the last place of 360, because the
    antimeridian vertices come out at 180.00000000000003 at some
    resolutions. Two vertices with the *same* longitude are a short side,
    not a repeated point: a border trapezoid can have a base too short to
    survive in a longitude, and it is still a cell with an area.
    """
    turn = _FULL_TURN_ULPS * math.ulp(360.0)
    for i, (lon_i, lat_i) in enumerate(ring):
        for j in range(i + 1, len(ring)):
            lon_j, lat_j = ring[j]
            if lat_i == lat_j and abs(abs(lon_i - lon_j) - 360.0) <= turn:
                return i, j
    return None


def _longitude_span(ring: list[tuple[float, float]]) -> float:
    """How many degrees of longitude a ring reaches, end to end.

    A ring that reaches a full turn or more wraps the whole parallel, and
    inside the polar caps, from resolution 4 down, many rings wrap it more
    than once: their longitudes span from just over a turn to some eighty
    turns, and no two of their vertices need be a full turn apart for the
    ring to cross itself.
    """
    longitudes = [longitude for longitude, _ in ring]
    return max(longitudes) - min(longitudes)


def _planar_area(points: list[tuple[float, float]]) -> float:
    """Shoelace area of a closed planar ring, in square metres."""
    total = 0.0
    for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1]):
        total += x1 * y2 - x2 * y1
    return abs(total) / 2.0


def _planar_perimeter(points: list[tuple[float, float]]) -> float:
    """Perimeter of a closed planar ring, in metres."""
    return sum(
        math.hypot(x2 - x1, y2 - y1)
        for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1])
    )


def _ring_of(cell: str) -> list[tuple[float, float]]:
    return cast("list[tuple[float, float]]", cell_to_boundary(cell))


def _projected_ring(cell: str) -> list[tuple[float, float]]:
    """The cell's ring in its own LAEA, refusing a non-simple ring."""
    ring = _ring_of(cell)
    repeated = _repeated_vertex(ring)
    if repeated is not None:
        i, j = repeated
        raise GeometryError(
            f"{cell!r} has no compactness: its boundary visits the same point "
            f"at positions {i} and {j}, so the ring is not a simple closed "
            f"curve and the area it encloses is not the area of the cell"
        )
    span = _longitude_span(ring)
    if span >= 360.0 - _FULL_TURN_ULPS * math.ulp(360.0):
        raise GeometryError(
            f"{cell!r} has no compactness: its ring wraps the whole parallel, "
            f"its longitudes spanning {span:.1f} degrees, so the ring is not a "
            f"simple closed curve and the area it encloses is not the area of "
            f"the cell"
        )
    lon0, lat0 = cast("tuple[float, float]", cell_to_centroid(cell))
    return [_laea(lon, lat, lon0, lat0) for lon, lat in ring]


# --------------------------------------------------------------------------
# Public metrics
# --------------------------------------------------------------------------


def compactness(cell: str) -> float | list[float]:
    """Isoperimetric quotient of a cell, ``4 * pi * A / p ** 2``.

    A transcription of equation (1) of section 2.2 of Kmoch, Vasilyev,
    Virro and Uuemaa (2022), *Big Earth Data* 6(3), 256-275. The metric
    is unitless, a circle scores 1, and less compact shapes fall towards
    zero.

    ``A`` and ``p`` are **planar**, measured in the oblique ellipsoidal
    Lambert Azimuthal Equal Area projection centred on the cell's own
    centroid, which is step 3 of the paper's workflow. They are not
    geodesic and not measured on the sinusoidal plane. For a cell of 10
    km the difference against the geodesic quantities is of the order of
    ``(d / R) ** 2``, about ``2e-6``; the paper's convention is kept
    anyway, because the purpose of the number is to be comparable with
    the paper's figures.

    The projection runs on the authalic latitude, whose cosine is formed
    from the area between the latitude and the pole, so the quotient keeps
    its digits near the pole, where an arcsine of a ratio within parts in a
    trillion of one did not.

    Values to expect, all measured: a regular hexagon scores 0.9069, a
    square 0.7854 and an equilateral triangle 0.6046, which is the ladder
    the paper's Figure 11(b) puts H3, S2 and the triangular grids on.
    ITACaRT cells fall between the triangular and the square grids, from
    0.539 at the equator to 0.785 at row 900, where the cell is a square
    to five decimal places. That is the price of the 45-degree shear and
    it is a number to declare, not a defect to hide.

    Args:
        cell: Compositional index string.

    Returns:
        The quotient for a single terminal cell, or a positionally
        aligned list when the index holds several.

    Raises:
        GeometryError: If the cell's boundary is not a simple closed
            curve. The members of that class are the cells whose ring
            wraps the whole parallel, all of them inside the polar caps:
            the cell holding the pole, at every resolution -- the polar cap
            and, below it, the one triangle per resolution that still spans
            the parallel -- and, from resolution 4 down, many cells around
            that triangle whose rings wrap the parallel once or many times
            over. Such a ring crosses itself, so the area it encloses is not
            the cell's.
            Step 2 of the paper's own workflow removes cells whose
            geometry is invalid for exactly this kind of reason. A border
            trapezoid is not in it however short its base: a base too
            short to survive in a longitude is still a side, and the cell
            still has an area.
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    require_existing_cells(cell)
    values = []
    for atom in iter_cells(cell):
        points = _projected_ring(atom)
        perimeter = _planar_perimeter(points)
        values.append(4.0 * math.pi * _planar_area(points) / (perimeter * perimeter))
    return values[0] if is_atomic(cell) else values


def cell_base_angle(
    cell: str,
    output: AngleOutput = "degrees",
) -> float | list[float]:
    """Ellipsoidal image of the 45-degree angle of the sinusoidal plane.

    On the sinusoidal plane every ITACaRT cell is the same parallelogram:
    Figure 1 of Silva, Dietzsch and Shiguemori (2025) places its vertices
    at ``{x, y}``, ``{x + l, y}``, ``{x, y + l}`` and ``{x - l, y + l}``,
    so the base equals the height and the acute angle is exactly 45
    degrees. This function reports what that angle becomes once the cell
    is carried to the WGS84 ellipsoid.

    It is measured at the **anchor**, between two lines. One is the
    **parallel**, which is not a geodesic: its direction at the anchor is
    the due east and due west tangent, and the arm taken is the one the
    45 degrees opens against, opposite the cell's base. The other is the
    **geodesic** from the anchor to the vertex at the most negative X in
    the cell's own frame -- the corner the shear displaces -- in a
    parallelogram, a triangle or a trapezoid alike.

    Reading the base as a geodesic chord instead of as a parallel biases
    the result by half the convergence of the meridians across the base.
    That bias is zero at the equator, where a parallel is a geodesic,
    and it grows with latitude: 0.023 degree at row 300 and 0.284 at row
    900, where it is the whole difference between 89.426479 and
    89.710055.

    The angle is a function of ``(X, Y)``, the lattice position, not of
    latitude. Along row 0 it grows with longitude at 3.95e-4 degree per
    degree, from 45.000035 at column 1 to 45.070984 at column 2000, at
    constant latitude throughout. It is also not monotone: it reaches
    89.426479 at row 900, where the cell is a square, and keeps opening
    past a right angle, to 134.375799 at row 950 and 142.6159 in
    Kamchatka, which is the "severely obtuse angle exceeding 135
    degrees" that Figure 8 of the paper shows as a picture. The four
    quadrants agree to nine decimals at every one of those latitudes.

    The frame is the cell's own. Which way a cell leans is read from the
    side of the prime meridian its anchor lies on, on the plane, and not
    from its quadrant code: the western half of a meridian triangle and of
    a polar cap is spelled under an eastern code and leans the western way.
    Such a cell and its eastern mirror measure the same angle to a few
    units in the last place.

    The angle holds at every resolution. The geodesic is Vincenty's, which
    holds its azimuth to the rounding of its inputs at every length, so a
    centimetre cell is measured to about 4e-5 degree: the rounding of
    coordinates written in degrees, not the method.

    Triangular cells return 0, and the zero is measured rather than
    declared by a clause. The anchor of a meridian triangle is not a
    vertex of its ring: it sits at longitude 0 in the middle of the base,
    and the vertex at the most negative X is the far end of that base, at
    the same latitude. The side leaving the anchor towards it *is* the
    parallel, so the angle against the parallel is zero. The polar caps
    are triangles whose anchor is on Greenwich at the middle of their
    base, and they fall under the same rule without a clause of their
    own. The test is exact rather than a tolerance: a latitude depends on
    the ordinate alone, the leaning vertex of a triangle shares the
    anchor's ordinate, and every other shape puts it one side higher, at
    every resolution. Zero therefore means a triangle and nothing else.

    Args:
        cell: Compositional index string.
        output: ``"degrees"`` for the angle, ``"radians"`` for the same
            angle in radians, or ``"sin"`` for its sine. The sine is the
            multiplicand of ``base * side * sin(angle)``, so it is the
            sine of this angle and of no other quantity.

    Returns:
        The angle for a single terminal cell, or a positionally aligned
        list when the index holds several.

    Raises:
        GeometryError: If ``output`` is not one of the three named forms.
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    require_existing_cells(cell)
    if output not in ("degrees", "radians", "sin"):
        raise GeometryError(
            f"output must be 'degrees', 'radians' or 'sin', not {output!r}"
        )
    values = []
    for atom in iter_cells(cell):
        values.append(_angle_in(_atom_base_angle(atom), output))
    return values[0] if is_atomic(cell) else values


def _local_offsets(
    cell: str, ring: list[tuple[float, float]], anchor: tuple[float, float]
) -> list[float]:
    """Each ring vertex's X offset from the anchor, in the cell's own frame.

    Measured on the sinusoidal plane and not in longitude, because the
    45 degrees is a property of that plane: Figure 1 of the paper puts
    the vertices at ``{x, y}``, ``{x + l, y}``, ``{x, y + l}`` and
    ``{x - l, y + l}``, so the offsets are exactly ``0``, ``+l``, ``0``
    and ``-l`` at every latitude. In longitude they are not: the
    meridians converge, and above roughly row 900 the vertex that the
    shear puts to the west of the anchor comes out east of it, which
    would pick the wrong vertex for the very cells where the angle is
    most interesting.

    The sign follows the cell's chirality, so ``+X`` always runs the way
    the lattice counts columns on the cell's side of the prime meridian:
    east of it in the east, west of it in the west, where the shear falls
    out negative all the same. The side is the anchor's, read off its
    abscissa on the plane, and not the quadrant's letter: the western half
    of a meridian triangle, and of a polar cap, is tiled with the western
    chirality under an eastern quadrant code.
    """
    east = _chirality(cell, anchor)
    anchor_x, _ = geodetic_to_sinusoidal(*anchor)
    return [east * (geodetic_to_sinusoidal(*v)[0] - anchor_x) for v in ring]


def _chirality(cell: str, anchor: tuple[float, float]) -> float:
    """``+1`` where the cell leans like the east, ``-1`` like the west.

    The side of the prime meridian the anchor lies on, read on the plane
    from its unwrapped longitude. Not the quadrant letter: the western half
    of a meridian triangle and of a polar cap carries an eastern code and
    leans the western way, and reading the letter there picked the other
    end of the base and answered zero for exactly half of the cells under
    every triangle. Not a normalised longitude either: an extension cell
    east of the antimeridian has its anchor past 180 degrees, east of the
    meridian on the plane, and folding it to -175 would turn it west.
    """
    anchor_x, _ = geodetic_to_sinusoidal(*anchor)
    side = _side_quadrant(cast(str, quadrant_of(cell)), anchor_x)
    return 1.0 if side[1] == "E" else -1.0


def _atom_base_angle(cell: str) -> float:
    """Angle at the anchor, between the parallel and the leaning side."""
    ring = _ring_of(cell)
    anchor = cast("tuple[float, float]", cell_to_anchor(cell))
    offsets = _local_offsets(cell, ring, anchor)
    target = ring[min(range(len(ring)), key=lambda k: offsets[k])]
    # Exact, because the construction is: a latitude depends only on the
    # ordinate, and the leaning vertex of a triangle shares its ordinate
    # with the anchor while every other shape puts it one side higher.
    # A tolerance in degrees here once swallowed whole cells, and from
    # resolution 11 down every cell is shorter than the tolerance was.
    if target[1] == anchor[1]:
        return 0.0
    reference = 270.0 if _chirality(cell, anchor) > 0.0 else 90.0
    _, azimuth = inverse_geodesic(*anchor, *target)
    separation = abs(azimuth - reference) % 360.0
    return separation if separation <= 180.0 else 360.0 - separation


def _angle_in(degrees: float, output: AngleOutput) -> float:
    if output == "degrees":
        return degrees
    radians = math.radians(degrees)
    return radians if output == "radians" else math.sin(radians)


def normalized_cell_area(cell: str) -> float | list[float]:
    """Cell area divided by the nominal area of its resolution.

    The second metric of section 2.2 of Kmoch et al. (2022), with one
    documented substitution. The paper divides each cell's area by the
    **mean area of all cells at the same resolution**; that mean is not
    computable here, because resolution 1 alone holds about four million
    cells and every finer level multiplies that. The denominator used
    instead is :func:`itacart.resolutions.nominal_cell_area`, the area
    the resolution table assigns to the level.

    The substitution is safe to state and unsafe to leave unstated. It is
    safe because the sinusoidal plane is equal-area by construction, so
    every unclipped cell has exactly the nominal area and scores exactly
    1.0 -- against plus or minus 50 percent for H3 and S2, whose
    normalized area has a standard deviation of 0.13 and 0.15 in the
    paper. It has to be stated because it is a different denominator, and
    a claim that names its own substitution is worth more than one that
    hides it.

    The exceptions are the boundary families, and they are large.
    Enumerated over all four quadrants at resolution 1: 4000 trapezoids
    range from 0.032263 to 2.067292, and 4004 triangles from 0.121394 to
    1.0, the minimum of the latter being the four polar cap cells.

    Args:
        cell: Compositional index string.

    Returns:
        The ratio for a single terminal cell, or a positionally aligned
        list when the index holds several.

    Raises:
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    require_existing_cells(cell)
    values = []
    for atom in iter_cells(cell):
        area = cast(float, effective_cell_area(atom))
        values.append(area / nominal_cell_area(get_resolution(atom)))
    return values[0] if is_atomic(cell) else values


def cell_max_uncertainty(
    cell: str,
    output: UncertaintyOutput = "radius",
) -> float | list[float]:
    """How far a position written as a cell can lie from the cell's centre.

    A position recorded as a cell is known to within that cell, and is
    read back at the cell's representative position, its centroid -- the
    point :func:`itacart.cells.cell_to_centroid` returns and every
    single-point route of the package answers with. This function measures
    the worst case of that substitution on the ground: the geodesic
    distance, on the WGS84 ellipsoid, from the centroid to the farthest
    vertex of the cell's ring. It is the *reach* that
    :func:`itacart.geometry.cells_to_geometry` states for a rebuilt vertex.

    Every distance is a geodesic solved by
    :func:`itacart.geodesy.inverse_geodesic`. A ring has three or four
    vertices, so a radius takes four geodesics at most and a diameter six,
    whatever the cell, and the cost per cell is constant. Nothing is
    measured on the sinusoidal plane: there every regular parallelogram is
    the same figure, with diagonals ``l`` and ``sqrt(5) l``, and the
    measure would not see where the cell is.

    **On the ground the position shows.** For a regular parallelogram of
    side ``l`` the radius follows the shear of the projection, to first
    order in the side::

        l * max(sqrt(1 + u**2) / 2, sqrt((1 - u / 2)**2 + 1 / 4))
        u = |lon * sin(lat)|, lon in radians from Greenwich, at the centroid

    which is ``1.118 l`` on the prime meridian, ``l / sqrt(2)`` where
    ``u = 1``, which is where the cell is a square on the ground with a
    diameter of ``sqrt(2) l``, and ``1.65 l`` as ``u`` reaches pi. In
    Chukotka, east of the antimeridian, ``u`` passes pi, read on the
    unwrapped longitude the package returns past 180 degrees. A finite cell
    departs from the law in proportion to its side. Over the regular
    parallelograms of resolution-1 rows 943 and 990 to 998 in the four
    quadrants and of the cells past the antimeridian in three rows of
    Chukotka, and over a lattice in latitude and ``u`` at resolutions 2 to
    13 that reaches the top of row 998 and ``u`` past pi, the radius
    exceeds the law by at most ``4.2e-3 * l / 10 km`` below 85 degrees of
    latitude and ``0.2 * l / 10 km`` from 85 degrees to row 998, plus up to
    ``1e-6`` from the rounding of the coordinates, which is what sets the
    departure below 85 degrees from resolution 11 down. It never falls
    short of the law by more than that rounding. In a parallelogram the
    diameter is twice the radius within the same bound. In a triangle, a
    trapezoid or a cap it is not, and no relation between the two is
    claimed there.

    **The maxima are taken at the vertices**, not over the continuous
    boundary. The two agree, to the rounding of the coordinates, in every
    family measured -- interior, border, polar row, extension, the cells
    the quantizer answers with near each pole at every resolution, and
    every cell of the caps' physical refinement down to resolution 5 --
    with one exception at every resolution: the cell that holds a pole. Its
    base is a whole parallel and its two base vertices, at -180 and +180
    degrees, are one point, so its diameter over the vertices is half the
    diameter of its boundary, 0.196573 l against 0.393146 l at resolution
    1. Its radius agrees. Under a cap, a spelling outside the cap's
    physical refinement can also reach further along its boundary than at
    its vertices. Over the whole subtree at resolution 4, with every edge
    cut into eight, the boundary's radius reaches 1.156 times the
    vertices' and its diameter 1.825 times; in a sample at resolution 5,
    2.648 and 3.916 times. Finer cutting can only find more.

    The number is a property of the tessellation and not the accuracy of a
    survey: it says how far the cell's centre can be from a position the
    cell holds, not how well that position was measured.

    Args:
        cell: Compositional index string.
        output: ``"radius"``, the default, for the distance from the
            centroid to the farthest vertex, in metres; ``"diameter"`` for
            the largest distance between two vertices, in metres; or
            ``"ratio"`` for the radius over
            :func:`itacart.resolutions.cell_size` of the cell's
            resolution, unitless.

    Returns:
        The value for a single terminal cell, or a positionally aligned
        list when the index holds several.

    Raises:
        GeometryError: If ``output`` is not one of the three named forms.
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    require_existing_cells(cell)
    if output not in ("radius", "diameter", "ratio"):
        raise GeometryError(
            f"output must be 'radius', 'diameter' or 'ratio', not {output!r}"
        )
    values = []
    for atom in iter_cells(cell):
        values.append(_atom_uncertainty(atom, output))
    return values[0] if is_atomic(cell) else values


def _atom_uncertainty(cell: str, output: UncertaintyOutput) -> float:
    """Radius, diameter or ratio of one cell, over the vertices of its ring."""
    ring = _ring_of(cell)
    if output == "diameter":
        return max(
            inverse_geodesic(*ring[i], *ring[j])[0]
            for i in range(len(ring))
            for j in range(i + 1, len(ring))
        )
    centre = cast("tuple[float, float]", cell_to_centroid(cell))
    radius = max(inverse_geodesic(*centre, *vertex)[0] for vertex in ring)
    if output == "radius":
        return radius
    return radius / cell_size(get_resolution(cell))
