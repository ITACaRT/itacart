"""The centre rule where a centre lies on the outline.

Under ``center`` a cell is kept when its centre lies inside the geometry,
and ``contains`` is false on the outline, so a centre exactly on it used to
be kept by no region at all. The centre is now moved one micrometre east
and ten nanometres north on the parallels plane before it is asked about.
A centre on an edge then falls on one side of it: of two regions sharing
the edge exactly one keeps it, west and south edges keeping the centres on
them and east and north edges dropping them.

Most polygons here are drawn through cell centroids at resolution 13, where
the outline runs through the lattice of centres and every vertex is such a
tie. They are written in lattice steps -- ``(i, j)`` counts steps to the
east and to the north neighbour of an origin cell -- and the count is held
to an exact integer reference on the same steps: the lattice points a
polygon keeps when each is moved an infinitesimal distance in the direction
of the shift. On the lattice that direction is 101 steps east to one north
in the north-eastern and south-western quadrants and 99 to one in the other
two, which the first test reads off the grid rather than assuming.

The vertices of such a polygon reach the plane through the projection,
whose last bits differ between C libraries, so before anything is tested
the ``center`` rule moves every vertex lying within a micrometre of the
lattice onto it: onto the five-millimetre plane grid that holds every
corner and centre of an ordinary cell, or onto a prime-meridian
centroid. The tests that simulate another machine's bits check that the
move erases them, and one turns it off to show what it settles.

What the rule does not settle is pinned at the end, in two parts: what it
decides other than the exact count but alike on every machine, by number,
and what is still left to the machine's bits, by existence.
"""

from __future__ import annotations

import math
import random
from fractions import Fraction
from functools import lru_cache

import pytest
from shapely import get_coordinates
from shapely.affinity import affine_transform
from shapely.geometry import LineString, MultiPolygon, Point, Polygon, box

import itacart
from itacart import boundary
from itacart import geometry as G
from itacart.exceptions import GeometryError
from itacart.geodesy import geodetic_to_sinusoidal, sinusoidal_to_geodetic

RES = 13
QUADRANTS = ("NE", "NW", "SE", "SW")

Ring = list[tuple[int, int]]

#: Where each quadrant's lattice is taken: near (-47, -22), where the rule
#: was first measured, and its mirror images in the other three quadrants.
ORIGINS = {
    "NE": (47.0, 22.0),
    "NW": (-47.0, 22.0),
    "SE": (47.0, -22.0),
    "SW": (-47.0, -22.0),
}

#: The same, close to the poles, where the ordinate approaches ten thousand
#: kilometres and the north part of the shift is five units in its last place.
NEAR_THE_POLES = {
    "NE": (40.0, 88.0),
    "NW": (-40.0, 88.0),
    "SE": (40.0, -88.0),
    "SW": (-40.0, -85.0),
}

#: The shift's direction in lattice steps, east then north, by quadrant.
STEP_DIRECTION = {"NE": (101, 1), "NW": (99, 1), "SE": (99, 1), "SW": (101, 1)}


# --------------------------------------------------------------------------
# The lattice of centres, addressed in steps
# --------------------------------------------------------------------------


class _Lattice:
    """Resolution-13 cell centres around an origin, addressed in lattice steps.

    A step east or north is the plane offset from the origin cell's
    centroid to its neighbour's, so position ``(i, j)`` is a cell centre,
    and a polygon drawn through such positions has its vertices on the
    lattice, as a polygon drawn through cell centroids does.
    """

    def __init__(self, lon: float, lat: float) -> None:
        origin = itacart.geo_to_cell(lon, lat, RES)
        self.quadrant = origin[:2]
        self._o = self._plane(origin)
        east = self._plane(itacart.get_neighbor(origin, "E"))
        north = self._plane(itacart.get_neighbor(origin, "N"))
        self._e = (east[0] - self._o[0], east[1] - self._o[1])
        self._n = (north[0] - self._o[0], north[1] - self._o[1])

    @staticmethod
    def _plane(cell: str) -> tuple[float, float]:
        return geodetic_to_sinusoidal(*itacart.cell_to_centroid(cell))

    def cell(self, i: int, j: int) -> str:
        x = self._o[0] + i * self._e[0] + j * self._n[0]
        y = self._o[1] + i * self._e[1] + j * self._n[1]
        return itacart.geo_to_cell(*sinusoidal_to_geodetic(x, y), RES)

    def vertices(self, ring: Ring) -> list[tuple[float, float]]:
        return [itacart.cell_to_centroid(self.cell(i, j)) for i, j in ring]

    def polygon(self, ring: Ring) -> Polygon:
        return Polygon(self.vertices(ring))

    def shift_in_steps(self) -> tuple[float, float]:
        """The package's plane shift, written in east and north steps."""
        east, north = G._tie_shift(None)
        (ex, ey), (nx, ny) = self._e, self._n
        det = ex * ny - ey * nx
        return ((east * ny - north * nx) / det, (ex * north - ey * east) / det)


@lru_cache(maxsize=None)
def _lattice(lon: float, lat: float) -> _Lattice:
    return _Lattice(lon, lat)


def _count(lattice: _Lattice, ring: Ring) -> int:
    return itacart.count_internal_cells(lattice.polygon(ring), RES)


def _moved_by_one_ulp(
    points: list[tuple[float, float]], seed: int
) -> list[tuple[float, float]]:
    """Every coordinate moved by 0 or 1 unit in the last place, either way.

    The size of the difference between two C libraries' projection of the
    same vertex, so another platform is simulated without needing it.
    """
    rng = random.Random(seed)

    def move(value: float) -> float:
        if rng.random() < 2 / 3:
            return math.nextafter(value, rng.choice((-1.0, 1.0)) * math.inf)
        return value

    return [(move(x), move(y)) for x, y in points]


#: How far another machine's bits are simulated beyond the one unit above:
#: 64 units in the last place are 25 to 100 nanometres on the plane at these
#: latitudes, several times the 7.5 nanometres the projection's round trip
#: was measured at, and a tenth of the micrometre the vertices move by or less.
FAR_BITS = 64


def _moved_by_ulps(
    points: list[tuple[float, float]], seed: int, most: int
) -> list[tuple[float, float]]:
    """Every coordinate moved by up to ``most`` units in the last place, either way."""
    rng = random.Random(seed)

    def move(value: float) -> float:
        steps = rng.randint(-most, most)
        for _ in range(abs(steps)):
            value = math.nextafter(value, math.copysign(math.inf, steps))
        return value

    return [(move(x), move(y)) for x, y in points]


# --------------------------------------------------------------------------
# The exact reference, in integers
# --------------------------------------------------------------------------


def _floor_sum(n: int, m: int, a: int, b: int) -> int:
    """``sum(floor((a * x + b) / m) for x in range(n))``, in logarithmic time."""
    total = 0
    q, a = divmod(a, m)
    total += q * n * (n - 1) // 2
    q, b = divmod(b, m)
    total += q * n
    while True:
        if a >= m:
            q, a = divmod(a, m)
            total += q * n * (n - 1) // 2
        if b >= m:
            q, b = divmod(b, m)
            total += q * n
        y_max = a * n + b
        if y_max < m:
            return total
        n, b = divmod(y_max, m)
        m, a = a, m


def _east_then_north(ring: Ring) -> int:
    """Lattice points kept when each is moved infinitesimally east, then north.

    Summed edge by edge: the points of each row strictly left of the edge,
    signed by the edge's direction. A point on an edge that is not
    horizontal goes to the side east of it, and on a horizontal edge to the
    side north of it.
    """
    area2, half = 0, 0
    for k, (x0, y0) in enumerate(ring):
        x1, y1 = ring[(k + 1) % len(ring)]
        area2 += x0 * y1 - x1 * y0
        if y0 != y1:
            sign = 1 if y1 > y0 else -1
            (xl, yl), (xu, yu) = (
                ((x0, y0), (x1, y1)) if y1 > y0 else ((x1, y1), (x0, y0))
            )
            dy, dx = yu - yl, xu - xl
            half += sign * _floor_sum(dy, dy, dx, xl * dy + dy - 1)
    return -half if area2 < 0 else half


def _to_east(direction: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    """An integer matrix of determinant one taking ``direction`` to ``(1, 0)``."""
    p, q = direction
    old_r, r, old_x, x, old_y, y = p, q, 1, 0, 0, 1
    while r:
        k = old_r // r
        old_r, r = r, old_r - k * r
        old_x, x = x, old_x - k * x
        old_y, y = y, old_y - k * y
    if old_r < 0:
        old_r, old_x, old_y = -old_r, -old_x, -old_y
    assert old_r == 1, "the direction must be primitive"
    return ((old_x, old_y), (-q, p))


def _reference(ring: Ring, direction: tuple[int, int]) -> int:
    """Lattice points kept when each is moved infinitesimally along ``direction``.

    The unimodular map keeps the lattice and the orientation, so the count
    is the count of the image under :func:`_east_then_north`. An edge
    parallel to ``direction`` is settled by the second direction the map
    leaves, which is this reference's own convention: nothing below holds
    the package to it.
    """
    (a, b), (c, d) = _to_east(direction)
    return _east_then_north([(a * i + b * j, c * i + d * j) for i, j in ring])


def _inside(point: tuple[Fraction, Fraction], ring: Ring) -> bool:
    """Even-odd test in exact arithmetic, for a point off the outline."""
    x, y = point
    inside = False
    for k, (x0, y0) in enumerate(ring):
        x1, y1 = ring[(k + 1) % len(ring)]
        if (y0 > y) != (y1 > y):
            if x0 + (y - y0) * Fraction(x1 - x0) / (y1 - y0) > x:
                inside = not inside
    return inside


def _brute_force(ring: Ring, direction: tuple[int, int]) -> int:
    """The same count, point by point, each moved a nanostep along ``direction``."""
    eps = Fraction(1, 10**9)
    p, q = direction
    columns = [i for i, _ in ring]
    rows = [j for _, j in ring]
    return sum(
        _inside((i + eps * p, j + eps * q), ring)
        for i in range(min(columns), max(columns) + 1)
        for j in range(min(rows), max(rows) + 1)
    )


def _has_an_edge_along(ring: Ring, direction: tuple[int, int]) -> bool:
    p, q = direction
    return any(
        (x1 - x0) * q == (y1 - y0) * p
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1])
    )


def _random_ring(rng: random.Random) -> Ring | None:
    """A small star-shaped polygon with integer vertices, or None if degenerate."""
    size = rng.choice((6, 10, 16))
    angles = sorted(rng.uniform(0.0, 2 * math.pi) for _ in range(rng.randint(3, 8)))
    ring: Ring = []
    for angle in angles:
        radius = rng.uniform(0.3 * size, size)
        point = (round(radius * math.cos(angle)), round(radius * math.sin(angle)))
        if not ring or ring[-1] != point:
            ring.append(point)
    if len(ring) > 2 and ring[0] == ring[-1]:
        ring.pop()
    shape = Polygon(ring) if len(ring) >= 3 else Polygon()
    if shape.is_empty or not shape.is_valid or shape.area == 0:
        return None
    return ring


def _interior_and_boundary(ring: Ring) -> tuple[int, int]:
    """Lattice points strictly inside, and on the outline, by Pick's theorem."""
    area2 = abs(
        sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]))
    )
    on_outline = sum(
        math.gcd(abs(x1 - x0), abs(y1 - y0))
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1])
    )
    return (area2 - on_outline + 2) // 2, on_outline


# --------------------------------------------------------------------------
# Shapes, in lattice steps
# --------------------------------------------------------------------------


def _circle(radius: int, step: int) -> Ring:
    """Lattice points of a circle every ``step`` rows, in order around it."""
    points = (
        [(math.isqrt(radius**2 - j * j), j) for j in range(0, radius + 1, step)]
        + [
            (-math.isqrt(radius**2 - j * j), j)
            for j in range(radius, -radius - 1, -step)
        ]
        + [(math.isqrt(radius**2 - j * j), j) for j in range(-radius, 0, step)]
    )
    ring: Ring = []
    for point in points:
        if not ring or ring[-1] != point:
            ring.append(point)
    return ring[:-1] if ring[0] == ring[-1] else ring


def _arc(ring: Ring, first: int, last: int) -> Ring:
    return [
        ring[(first + k) % len(ring)] for k in range((last - first) % len(ring) + 1)
    ]


#: 706 758 480 cells; no edge runs near the shift's direction, and none is
#: long enough on the lattice for a point to sit closer to it than the shift.
CIRCLE = _circle(15_000, 60)
SQUARE: Ring = [(0, 0), (1_000, 0), (1_000, 1_000), (0, 1_000)]
STRIP: Ring = [(0, 0), (50_000, 0), (50_000, 1_000), (0, 1_000)]
#: Six edges in six directions, each with hundreds of lattice points on it.
HEXAGON: Ring = [
    (0, 0),
    (600, -300),
    (1_500, 0),
    (1_800, 900),
    (1_200, 1_500),
    (300, 1_200),
]
SMALL_HEXAGON: Ring = [(i // 10, j // 10) for i, j in HEXAGON]
SHAPES = {"circle": CIRCLE, "square": SQUARE, "strip": STRIP, "hexagon": HEXAGON}

#: An edge running 1:101, the shift's own direction on the lattice of the
#: north-eastern and south-western quadrants, with eleven lattice points on it.
ALONG_THE_SHIFT: Ring = [(0, 0), (1_010, 10), (1_010, 500), (0, 500)]
#: A primitive edge of 283 m on the plane whose nearest lattice point lies
#: 0.35 micrometres from it, against 0.70 of the shift across it.
NEAR_THE_EDGE: Ring = [(0, 0), (40_001, 20_000), (0, 20_000)]
#: A strip 400 m tall on the south-western lattice. On the plane its west
#: edge runs 1 m east for every 99.99 m south and its east edge 1 m east for
#: every 100.01 m north. Across either, the lattice lines parallel to it lie
#: exactly as far apart as the shift reaches across it, about a micrometre,
#: so the shift carries the points of the line beside the edge onto the edge
#: itself: four on each, one every 100 m. In exact arithmetic they land on
#: it; in floating point each is a near tie.
ONE_SHIFT_FROM_THE_EDGE: Ring = [
    (0, 0),
    (-39_596, -39_996),
    (-35_596, -37_996),
    (4_808, 2_008),
]


# --------------------------------------------------------------------------
# The instrument
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_the_reference_direction_is_the_packages_shift(quadrant: str) -> None:
    """The direction the reference moves points in is read off the grid.

    The package's plane shift, written in each quadrant's east and north
    steps, points 101 steps east to one north in two quadrants and 99 to
    one in the other two: the shear adds the north part to the east step
    in the first pair and takes it away in the second. If the shift or the
    lattice moved, this is where the reference would stop describing the
    rule.
    """
    lattice = _lattice(*ORIGINS[quadrant])
    assert lattice.quadrant == quadrant
    east_steps, north_steps = lattice.shift_in_steps()
    p, q = STEP_DIRECTION[quadrant]
    assert east_steps > 0 and north_steps > 0
    assert east_steps / north_steps == pytest.approx(p / q, rel=1e-5)


def test_the_reference_agrees_with_brute_force_on_small_polygons() -> None:
    """The reference is checked against a count it shares no code with.

    Random small polygons, every lattice point moved a billionth of a step
    along the direction and tested in exact rational arithmetic. Polygons
    with an edge parallel to the direction are left out: there the
    reference follows a convention of its own. The east-then-north rule of
    the reference is the special case ``(1, 0)``, and it reproduces the
    exact counts that were measured before the shift existed.
    """
    rng = random.Random(20261121)
    directions = [(101, 1), (99, 1), (1, 0), (-101, -1), (3, 7)]
    checked = 0
    while checked < 100:
        ring = _random_ring(rng)
        direction = directions[checked % len(directions)]
        if ring is None or _has_an_edge_along(ring, direction):
            continue
        assert _reference(ring, direction) == _brute_force(ring, direction), ring
        checked += 1
    assert _east_then_north(CIRCLE) == _reference(CIRCLE, (101, 1)) == 706_758_480
    assert _reference(SQUARE, (99, 1)) == 1_000_000


# --------------------------------------------------------------------------
# The count is the lattice count, and the bits do not move it
# --------------------------------------------------------------------------


@pytest.mark.parametrize("shape", sorted(SHAPES))
@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_the_count_is_the_integer_reference(quadrant: str, shape: str) -> None:
    """Drawn through cell centroids, the count is the exact lattice count.

    Every vertex and every lattice point along an edge is a centre on the
    outline, and each is decided by the direction of the shift rather than
    by the last bits of the projection. Before the shift every one of them
    was dropped: the square counted its interior, 998 001.
    """
    lattice = _lattice(*ORIGINS[quadrant])
    ring = SHAPES[shape]
    assert _count(lattice, ring) == _reference(ring, STEP_DIRECTION[quadrant])


@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_the_count_is_the_integer_reference_near_the_poles(quadrant: str) -> None:
    """The north part of the shift still decides, where it is five ulps.

    East-west edges are decided by the north part alone, ten nanometres,
    and near the pole that is five units in the last place of the
    ordinate. These squares sit at 85 and 88 degrees.
    """
    lattice = _lattice(*NEAR_THE_POLES[quadrant])
    for ring in (SQUARE, HEXAGON):
        assert _count(lattice, ring) == _reference(ring, STEP_DIRECTION[quadrant])


CASES_UNDER_THE_ULP = [
    ("SW", "circle", False),
    ("NW", "circle", False),
    *((quadrant, "hexagon", False) for quadrant in QUADRANTS),
    *((quadrant, "square", True) for quadrant in QUADRANTS),
]


@pytest.mark.parametrize(("quadrant", "shape", "polar"), CASES_UNDER_THE_ULP)
def test_one_ulp_on_the_vertices_does_not_move_the_count(
    quadrant: str, shape: str, polar: bool
) -> None:
    """Another C library's vertices, simulated: the count does not move.

    The rule used to depend on whether a centre landed exactly on the
    outline or one unit in the last place to either side, which is the
    size of the difference between two libraries' projection of the same
    vertex. Moving every coordinate by that much, eight times over, leaves
    the count where it was: with the vertices moved onto the lattice, which
    erases the difference, and without, which shows that on these shapes
    the shift holds on its own.
    """
    lattice = _lattice(*(NEAR_THE_POLES if polar else ORIGINS)[quadrant])
    points = lattice.vertices(SHAPES[shape])
    expected = itacart.count_internal_cells(Polygon(points), RES)
    for seed in range(8):
        moved = Polygon(_moved_by_one_ulp(points, seed))
        assert itacart.count_internal_cells(moved, RES) == expected, seed
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(G, "_SNAP_TOLERANCE", 0.0)
        for seed in range(-1, 8):
            moved = Polygon(points if seed < 0 else _moved_by_one_ulp(points, seed))
            assert itacart.count_internal_cells(moved, RES) == expected, seed


# --------------------------------------------------------------------------
# A shared edge gives each centre to one side
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ("SW", "NW"))
def test_cuts_through_the_circle_add_up_to_the_whole(quadrant: str) -> None:
    """Cut along a row, a column and a chord, the parts sum to the whole.

    Each cut is a lattice line with up to 29 999 centres on it. Before the
    shift every one was dropped by both halves: the row cut lost all of
    them.
    """
    lattice = _lattice(*ORIGINS[quadrant])
    whole = _count(lattice, CIRCLE)
    position = {point: k for k, point in enumerate(CIRCLE)}
    radius, rise = 15_000, 3_000
    run = math.isqrt(radius**2 - rise**2)
    cuts = {
        "row": ((radius, 0), (-radius, 0)),
        "column": ((0, radius), (0, -radius)),
        "chord": ((run, rise), (-run, -rise)),
    }
    for label, (start, stop) in cuts.items():
        first, last = position[start], position[stop]
        parts = _count(lattice, _arc(CIRCLE, first, last)) + _count(
            lattice, _arc(CIRCLE, last, first)
        )
        assert parts == whole, label


@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_a_shared_edge_gives_each_cell_to_one_side(quadrant: str) -> None:
    """Filled, the two halves of a cut hold every cell of the whole once."""
    lattice = _lattice(*ORIGINS[quadrant])
    ring = _circle(120, 12)
    position = {point: k for k, point in enumerate(ring)}
    first, last = position[(120, 0)], position[(-120, 0)]

    def filled(part: Ring) -> set[str]:
        return set(itacart.decompose(itacart.polyfill(lattice.polygon(part), RES)))

    north, south = filled(_arc(ring, first, last)), filled(_arc(ring, last, first))
    assert not north & south
    assert north | south == filled(ring)


@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_south_and_west_edges_keep_their_centres(quadrant: str) -> None:
    """The cells centred on the south and west edges are kept, the rest not."""
    lattice = _lattice(*ORIGINS[quadrant])
    side = 100
    ring = [(0, 0), (side, 0), (side, side), (0, side)]
    kept = set(itacart.decompose(itacart.polyfill(lattice.polygon(ring), RES)))
    assert len(kept) == side * side
    assert {lattice.cell(i, 0) for i in range(side)} <= kept
    assert {lattice.cell(0, j) for j in range(side)} <= kept
    assert not {lattice.cell(i, side) for i in range(side + 1)} & kept
    assert not {lattice.cell(side, j) for j in range(side + 1)} & kept


def test_the_count_lies_between_the_interior_and_the_closed_count() -> None:
    """No more than the closed polygon holds, and no less than its interior.

    Both bounds come from Pick's theorem, in integers. The two polygons the
    rule does not settle are included: whatever they are off by, they stay
    inside the bounds.
    """
    lattice = _lattice(*ORIGINS["SW"])
    for ring in (*SHAPES.values(), ALONG_THE_SHIFT, NEAR_THE_EDGE):
        interior, on_outline = _interior_and_boundary(ring)
        assert interior <= _count(lattice, ring) <= interior + on_outline


@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_polyfill_keeps_what_the_count_counts_where_every_vertex_is_a_tie(
    quadrant: str,
) -> None:
    """The fill and the count take every decision the same way."""
    lattice = _lattice(*ORIGINS[quadrant])
    for ring in (SMALL_HEXAGON, _circle(120, 12)):
        polygon = lattice.polygon(ring)
        counted = itacart.count_internal_cells(polygon, RES)
        assert itacart.count_cells(itacart.polyfill(polygon, RES)) == counted
        assert counted == _reference(ring, STEP_DIRECTION[quadrant])


# --------------------------------------------------------------------------
# The families that are not squares
# --------------------------------------------------------------------------


def _either_side_of_the_meridian(
    west: float, east: float, south: float, north: float
) -> tuple[Polygon, Polygon, Polygon]:
    """Two boxes sharing the meridian, and their union written with their edges.

    The union carries the two points where the meridian meets its outline.
    Written as one box instead, its northern and southern edges are single
    geodesics, longer than the parts' and bulging further, and the union
    holds cells that neither part holds.
    """
    return (
        Polygon([(0.0, south), (east, south), (east, north), (0.0, north)]),
        Polygon([(west, south), (0.0, south), (0.0, north), (west, north)]),
        Polygon(
            [
                (west, south),
                (0.0, south),
                (east, south),
                (east, north),
                (0.0, north),
                (west, north),
            ]
        ),
    )


def _filled(geometry: Polygon, resolution: int) -> set[str]:
    try:
        return set(itacart.decompose(itacart.polyfill(geometry, resolution)))
    except GeometryError:
        return set()


def _without_the_shift(patch: pytest.MonkeyPatch) -> None:
    """The rule the shift replaced: the centre asked about where it is.

    The vertices are read as projected, too: with a tolerance of zero the
    move leaves every vertex where it is, and the rule is the one released
    before, to the bit.
    """
    patch.setattr(G, "_TIE_EAST", 0.0)
    patch.setattr(G, "_TIE_NORTH", 0.0)
    patch.setattr(G, "_SNAP_TOLERANCE", 0.0)


@pytest.mark.parametrize(
    ("resolution", "west", "east", "south", "north"),
    [
        (3, -0.5, 0.5, 10.0, 11.0),
        (4, -0.3, 0.3, -20.3, -20.0),
        (5, -0.05, 0.05, 45.0, 45.05),
    ],
)
def test_the_meridian_triangles_go_to_the_east_side(
    resolution: int, west: float, east: float, south: float, north: float
) -> None:
    """A triangle's centroid lies on the meridian, and the east keeps it.

    The prime-meridian triangles are decided in the unsplit plane by their
    centroid, which has ``x`` zero at every resolution, so it lies on any
    outline that follows the meridian. The two boxes are mirror images, and
    the eastern one holds more cells by exactly the triangles; the parts
    add up to the union, in the count and in the fill. Without the shift
    both boxes drop the triangles and the union holds them.
    """
    parts = _either_side_of_the_meridian(west, east, south, north)
    counts = [itacart.count_internal_cells(part, resolution) for part in parts]
    assert counts[0] + counts[1] == counts[2]
    first, second, union = (_filled(part, resolution) for part in parts)
    assert not first & second
    assert first | second == union
    triangles = {cell for cell in union if itacart.cell_shape(cell) == "triangle"}
    assert triangles and triangles <= first
    assert counts[0] - counts[1] == len(triangles)

    with pytest.MonkeyPatch.context() as patch:
        _without_the_shift(patch)
        old = [itacart.count_internal_cells(part, resolution) for part in parts]
    assert old[0] == old[1] == counts[1]
    assert old[0] + old[1] < old[2] == counts[2]


@pytest.mark.parametrize("hemisphere", ("N", "S"))
def test_the_polar_caps_go_to_the_east_side(hemisphere: str) -> None:
    """The cap's centroid lies on the meridian, and the east keeps it.

    The polar row holds the two caps and nothing else, and each cap is
    decided by the centroid of its effective ring, in the unsplit plane.
    That centroid has ``x`` zero. A sliver east of the meridian reaching
    the cap keeps it, the mirror sliver to the west does not, and the
    union keeps it once. Without the shift neither sliver keeps it.
    """
    south, north = (89.98, 89.999) if hemisphere == "N" else (-89.999, -89.98)
    cap = f"{hemisphere}E(0000/1000)"
    parts = _either_side_of_the_meridian(-10.0, 10.0, south, north)
    for resolution in (1, 2):
        counts = [itacart.count_internal_cells(part, resolution) for part in parts]
        assert counts == [1, 0, 1], resolution
        first, second, union = (_filled(part, resolution) for part in parts)
        assert first == union and not second, resolution
        assert all(cell.startswith(cap[:-1]) for cell in first), resolution

    with pytest.MonkeyPatch.context() as patch:
        _without_the_shift(patch)
        assert [itacart.count_internal_cells(part, 1) for part in parts] == [0, 0, 1]


def _on_the_lattice(
    lattice_geometry: Polygon, quadrant: str, resolution: int, patch: pytest.MonkeyPatch
) -> tuple[int, set[str]]:
    """Count and fill a figure given on a quadrant's lattice, as the descent sees it.

    The public calls take a geodetic polygon, and no geodetic polygon puts a
    line exactly through the centroid of an absorbing cell's sheared ring.
    So the preparation is replaced by the figure itself, as the descent
    tests of the geometry module do, and everything after it is the public
    path: the screen, the absorbing walk and its leaves.
    """
    a, b, d, e = G._QUADRANT_SHEAR[quadrant]
    det = a * e - b * d
    plane = affine_transform(
        lattice_geometry, [e / det, -b / det, -d / det, a / det, 0.0, 0.0]
    )

    def prepared(
        _geometry: object, _resolution: int, densify: bool = True, *, snap: bool
    ) -> tuple[Polygon, list[tuple[str, G._LatticeView]]]:
        return plane, [(quadrant, G._LatticeView(lattice_geometry))]

    patch.setattr(G, "_prepare", prepared)
    anything = Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
    return itacart.count_internal_cells(anything, resolution), _filled(
        anything, resolution
    )


def test_an_absorbing_cell_goes_to_the_north_side() -> None:
    """The last column's cells are decided by their ring's centroid, shifted.

    The last cell of a row absorbs the strip between its square and the
    border, and its descendants that keep absorbing are decided by the
    centroid of their effective ring, on the quadrant's lattice. Two
    rectangles share a horizontal edge through that centroid: the north
    part keeps the cell, the south part does not, and their union keeps it
    once. Without the shift neither part keeps it.
    """
    row = 300
    last = boundary.last_lattice_column("NE", row, G._L1)
    root = f"NE({last:04d}/{row:04d})"
    walk = G._BorderWalk("NE")
    child = next(cell for cell in walk.children(root) if itacart.absorbs_border(cell))
    ring = walk.ring(child)
    assert ring is not None
    centre = ring.centroid
    reach, rise = 3_000.0, 400.0
    north = box(centre.x - reach, centre.y, centre.x + reach, centre.y + rise)
    south = box(centre.x - reach, centre.y - rise, centre.x + reach, centre.y)
    union = box(centre.x - reach, centre.y - rise, centre.x + reach, centre.y + rise)

    with pytest.MonkeyPatch.context() as patch:
        answers = [_on_the_lattice(part, "NE", 2, patch) for part in (north, south)]
        answers.append(_on_the_lattice(union, "NE", 2, patch))
    assert answers == [(1, {child}), (0, set()), (1, {child})]

    with pytest.MonkeyPatch.context() as patch:
        _without_the_shift(patch)
        old = [_on_the_lattice(part, "NE", 2, patch)[0] for part in (north, south)]
        old.append(_on_the_lattice(union, "NE", 2, patch)[0])
    assert old == [0, 0, 1]


# --------------------------------------------------------------------------
# One rule in every evaluation of it
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ("SW", "NW"))
def test_the_batched_and_the_recursive_count_agree_on_ties(quadrant: str) -> None:
    """The two descents shift the same centre the same way, to the last bit."""
    lattice = _lattice(*ORIGINS[quadrant])
    for ring in (_circle(3_000, 60), SQUARE, HEXAGON):
        polygon = lattice.polygon(ring)
        batched = itacart.count_internal_cells(polygon, RES)
        assert batched == G._count_internal_cells_recursive(polygon, RES)
        assert batched == _reference(ring, STEP_DIRECTION[quadrant])


def test_the_workers_count_ties_as_this_process_does(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two worker processes give the serial count, and take the shift with them.

    A worker is started by ``spawn`` and imports the package afresh, so it
    sees nothing this process changed after its own import. The shift
    travels with the work instead: with the rule it replaced set here, the
    workers follow it, and the pool still matches the serial count.
    """
    monkeypatch.setattr(G, "_PARALLEL_MIN_NODES", 0)
    monkeypatch.setattr(G, "_FRONTIER_PER_JOB", 4)
    ring = _circle(3_000, 60)
    polygon = _lattice(*ORIGINS["SW"]).polygon(ring)
    serial = itacart.count_internal_cells(polygon, RES)
    assert serial == _reference(ring, STEP_DIRECTION["SW"])
    assert itacart.count_internal_cells(polygon, RES, n_jobs=2) == serial

    _without_the_shift(monkeypatch)
    old = itacart.count_internal_cells(polygon, RES)
    assert old < serial
    assert itacart.count_internal_cells(polygon, RES, n_jobs=2) == old


def test_the_rule_it_replaced_fails_the_reference_and_the_cut(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control: without the shift, the same tests go red.

    Every centre on the outline is dropped, so the circle counts fewer
    cells than the lattice holds, the two halves of the row cut lose the
    centres on the cut between them, and a unit in the last place on the
    square's vertices moves its count. The figures depend on the platform,
    so only the direction is asserted.
    """
    lattice = _lattice(*ORIGINS["SW"])
    reference = _reference(CIRCLE, STEP_DIRECTION["SW"])
    _without_the_shift(monkeypatch)
    whole = _count(lattice, CIRCLE)
    assert whole < reference
    position = {point: k for k, point in enumerate(CIRCLE)}
    first, last = position[(15_000, 0)], position[(-15_000, 0)]
    parts = _count(lattice, _arc(CIRCLE, first, last)) + _count(
        lattice, _arc(CIRCLE, last, first)
    )
    assert parts < whole
    square = lattice.vertices(SQUARE)
    before = itacart.count_internal_cells(Polygon(square), RES)
    assert before < 1_000_000
    moved = [
        itacart.count_internal_cells(Polygon(_moved_by_one_ulp(square, seed)), RES)
        for seed in range(8)
    ]
    assert any(count != before for count in moved)


# --------------------------------------------------------------------------
# A vertex drawn on the lattice is moved onto it
# --------------------------------------------------------------------------


def _on_the_grid(x: float, y: float) -> bool:
    step = G._SNAP_STEP
    return round(x / step) * step == x and round(y / step) * step == y


def _landing(lon_lat: tuple[float, float]) -> set[tuple[float, float]]:
    """Where one vertex lands, projected and moved, under eight machines' bits."""
    return {
        G._snap_vertex(
            *geodetic_to_sinusoidal(*_moved_by_ulps([lon_lat], seed, FAR_BITS)[0])
        )
        for seed in range(8)
    }


@pytest.mark.parametrize("polar", (False, True))
@pytest.mark.parametrize("quadrant", QUADRANTS)
def test_a_corner_or_a_centre_lands_on_one_grid_point_whatever_the_bits(
    quadrant: str, polar: bool
) -> None:
    """Every resolution's corners and centroid, as the package returns them.

    Taken in degrees from ``cell_to_boundary`` and ``cell_to_centroid``,
    with each coordinate moved by up to ``FAR_BITS`` units in the last
    place, projected and moved onto the lattice: each lands on one point
    of the five-millimetre grid, the same under every simulated machine,
    and within ten nanometres of where the unmoved vertex projects.
    """
    lon, lat = (NEAR_THE_POLES if polar else ORIGINS)[quadrant]
    for resolution in range(1, RES + 1):
        cell = itacart.geo_to_cell(lon, lat, resolution)
        assert itacart.cell_shape(cell) == "parallelogram", cell
        for vertex in [itacart.cell_to_centroid(cell), *itacart.cell_to_boundary(cell)]:
            landed = _landing(vertex)
            assert len(landed) == 1, (cell, vertex)
            point = landed.pop()
            assert _on_the_grid(*point), (cell, point)
            assert math.dist(point, geodetic_to_sinusoidal(*vertex)) < 1e-8


@pytest.mark.parametrize("latitude", (45.0, -45.0))
def test_a_meridian_centroid_lands_on_its_formula_whatever_the_bits(
    latitude: float,
) -> None:
    """The triangle's centroid is not on the grid and is found by its formula.

    On the meridian, a third of a side above the triangle's base. At every
    resolution the centroid ``cell_to_centroid`` returns, with its bits
    moved, lands on one point: abscissa zero, and an ordinate that is
    ``k*side + side/3`` for a whole ``k``, computed as such.
    """
    for resolution in range(1, RES + 1):
        cell = itacart.geo_to_cell(0.0, latitude, resolution)
        assert itacart.is_triangular_cell(cell), cell
        centroid = itacart.cell_to_centroid(cell)
        landed = _landing(centroid)
        assert len(landed) == 1, cell
        x, y = landed.pop()
        assert x == 0.0 and math.copysign(1.0, y) == math.copysign(1.0, latitude)
        assert any(
            abs(y) == round((abs(y) - side / 3.0) / side) * side + side / 3.0
            for side in G._MERIDIAN_SIDES
        ), (cell, y)
        assert not _on_the_grid(x, y)
        assert math.dist((x, y), geodetic_to_sinusoidal(*centroid)) < 1e-8


def test_a_border_or_cap_centroid_is_left_where_it_is() -> None:
    """A vertex at the centroid of a border-absorbing cell or a cap stays put.

    A cell absorbing the domain border has a ring no formula describes, and
    so has a polar cap; the fill tests the centroid of that ring as GEOS
    computes it. A vertex drawn there, from the centroid the package
    returns, is not moved onto anything, so near it the bits still decide
    what they decided before.
    """
    cells = ["NE(0000/1000)", "SE(0000/1000)"]
    for quadrant, row in (("NE", 300), ("SW", 700)):
        last = boundary.last_lattice_column(quadrant, row, G._L1)
        root = f"{quadrant}({last:04d}/{row:04d})"
        children = G._BorderWalk(quadrant).children(root)
        cells.append(root)
        cells.extend(cell for cell in children if itacart.absorbs_border(cell))
    assert len(cells) > 4
    for cell in cells:
        x, y = geodetic_to_sinusoidal(*itacart.cell_to_centroid(cell))
        assert G._snap_vertex(x, y) == (x, y), cell


def test_the_move_reaches_a_micrometre_in_each_coordinate_and_no_further() -> None:
    """Within the tolerance on both coordinates a vertex moves; past it, it stays."""
    x, y = geodetic_to_sinusoidal(*ORIGINS["SW"])
    step, reach = G._SNAP_STEP, G._SNAP_TOLERANCE
    grid = (round(x / step) * step, round(y / step) * step)
    for dx, dy in ((0.9, 0.0), (0.0, -0.9), (-0.9, 0.9)):
        assert G._snap_vertex(grid[0] + dx * reach, grid[1] + dy * reach) == grid
    for dx, dy in ((1.1, 0.0), (0.0, 1.1), (-2.0, 2.0)):
        off = (grid[0] + dx * reach, grid[1] + dy * reach)
        assert G._snap_vertex(*off) == off

    side = itacart.cell_size(1)
    centroid = 498 * side + side / 3.0
    assert G._snap_vertex(0.9 * reach, -centroid + 0.5 * reach) == (0.0, -centroid)
    for off in ((1.1 * reach, centroid), (0.0, centroid + 1.1 * reach)):
        assert G._snap_vertex(*off) == off


def test_every_kind_of_geometry_has_its_lattice_vertices_moved() -> None:
    """Points, lines, polygons with holes and multipolygons alike.

    Each vertex half a micrometre off a grid point, in both coordinates,
    lands on it, and the geometry keeps its kind and its rings.
    """
    x, y = geodetic_to_sinusoidal(*ORIGINS["SW"])
    step, off = G._SNAP_STEP, 0.5e-6
    gx, gy = round(x / step) * step, round(y / step) * step

    def near(dx: int, dy: int) -> tuple[float, float]:
        return (gx + dx * step + off, gy + dy * step - off)

    shell = [near(0, 0), near(400, 0), near(400, 400), near(0, 400)]
    hole = [near(100, 100), near(100, 200), near(200, 200), near(200, 100)]
    other = [near(600, 0), near(700, 0), near(700, 100)]
    for geometry in (
        Point(near(3, 7)),
        LineString([near(0, 0), near(5, 9)]),
        Polygon(shell, [hole]),
        MultiPolygon([Polygon(shell, [hole]), Polygon(other)]),
    ):
        moved = G._snap_to_grid(geometry)
        assert moved.geom_type == geometry.geom_type
        before, after = get_coordinates(geometry), get_coordinates(moved)
        assert len(after) == len(before) and not (before == after).all()
        assert all(_on_the_grid(float(px), float(py)) for px, py in after)
    assert len(G._snap_to_grid(Polygon(shell, [hole])).interiors) == 1


def test_a_figure_the_move_makes_invalid_is_refused() -> None:
    """Two vertices 0.8 micrometres apart about one grid point become one.

    An hourglass whose waist is 0.8 micrometres wide is valid in degrees and
    on the plane. With both waist vertices moved onto the grid point between
    them its ring touches itself, and the ``center`` rule refuses it rather
    than repair it, in the count and in the fill. ``intersects`` and
    ``contains`` read the vertices as projected, and fill it.
    """
    x, y = geodetic_to_sinusoidal(*ORIGINS["SW"])
    step, half = G._SNAP_STEP, 0.4e-6
    gx, gy = round(x / step) * step, round(y / step) * step
    plane = [
        (gx + half, gy),
        (gx + 10.0, gy + 10.0),
        (gx - 10.0, gy + 10.0),
        (gx - half, gy),
        (gx - 10.0, gy - 10.0),
        (gx + 10.0, gy - 10.0),
    ]
    hourglass = Polygon([sinusoidal_to_geodetic(px, py) for px, py in plane])
    assert hourglass.is_valid
    with pytest.raises(GeometryError, match="once its vertices are moved onto"):
        itacart.count_internal_cells(hourglass, 9)
    with pytest.raises(GeometryError, match="once its vertices are moved onto"):
        itacart.polyfill(hourglass, 9)
    for containment in ("intersects", "contains"):
        assert itacart.polyfill(hourglass, 9, containment)


def test_vertices_on_the_lattice_decide_alike_on_every_machine() -> None:
    """Eight centres the shift puts on an edge, and the bits they used to follow.

    ``ONE_SHIFT_FROM_THE_EDGE`` has eight lattice points that the shift
    carries exactly onto an edge. With the vertices where the projection
    puts them, each one's side follows the vertices' last bits, and
    coordinates moved by up to ``FAR_BITS`` units in the last place --
    another machine, simulated -- move the count by up to four either way.
    Moved onto the lattice first, the vertices are the same numbers on
    every machine, and so is the count. It is not the exact count, which
    no finite shift gives here, but it is one count, within four of it.
    """
    lattice = _lattice(*ORIGINS["SW"])
    points = lattice.vertices(ONE_SHIFT_FROM_THE_EDGE)
    reference = _reference(ONE_SHIFT_FROM_THE_EDGE, STEP_DIRECTION["SW"])

    def counts() -> list[int]:
        return [
            itacart.count_internal_cells(
                Polygon(_moved_by_ulps(points, seed, FAR_BITS)), RES
            )
            for seed in range(12)
        ]

    moved = counts()
    assert len(set(moved)) == 1
    assert abs(moved[0] - reference) <= 4
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(G, "_SNAP_TOLERANCE", 0.0)
        as_projected = counts()
    assert len(set(as_projected)) > 1
    assert all(abs(count - reference) <= 4 for count in as_projected)


# --------------------------------------------------------------------------
# What the rule decides otherwise than the exact count, alike everywhere
# --------------------------------------------------------------------------


def test_an_edge_along_the_shift_is_left_to_the_arithmetic_on_its_points() -> None:
    """The shift cannot take a centre off an edge that runs its own way.

    The edge from ``(0, 0)`` to ``(1010, 10)`` runs 1:101, the shift's
    direction on this quadrant's lattice, and carries eleven lattice points
    counting its ends. The shift slides them along the edge rather than off
    it, so their side is settled by the rounding of each shifted point and
    not by the rule. With the vertices moved onto the lattice that rounding
    is the same on every machine, and vertices up to ``FAR_BITS`` units in
    the last place away give one count. No other point is near enough to
    any edge to move, so it stays within eleven of the exact count.
    """
    lattice = _lattice(*ORIGINS["SW"])
    reference = _reference(ALONG_THE_SHIFT, STEP_DIRECTION["SW"])
    on_the_edge = math.gcd(1_010, 10) + 1
    points = lattice.vertices(ALONG_THE_SHIFT)
    counted = {
        itacart.count_internal_cells(
            Polygon(_moved_by_ulps(points, seed, FAR_BITS)), RES
        )
        for seed in range(8)
    }
    assert len(counted) == 1
    assert abs(counted.pop() - reference) <= on_the_edge


def test_a_lattice_point_nearer_than_the_shift_changes_side() -> None:
    """The shift departs from the exact count by one here, on every platform.

    The long edge of this triangle is 283 m on the plane and primitive on
    the lattice, and its nearest lattice point lies 0.35 micrometres from
    it, inside the triangle. The shift's component across the edge is 0.70
    micrometres and points out of it, so that point is moved out: the
    count is one less than the lattice count. Nothing here is a tie, and
    vertices a unit in the last place away do not change it.
    """
    lattice = _lattice(*ORIGINS["SW"])
    points = lattice.vertices(NEAR_THE_EDGE)
    counted = itacart.count_internal_cells(Polygon(points), RES)
    assert counted - _reference(NEAR_THE_EDGE, STEP_DIRECTION["SW"]) == -1
    for seed in range(4):
        moved = Polygon(_moved_by_one_ulp(points, seed))
        assert itacart.count_internal_cells(moved, RES) == counted, seed


# --------------------------------------------------------------------------
# What is still left to the machine
# --------------------------------------------------------------------------


def test_the_densification_threshold_is_the_one_the_documentation_states() -> None:
    """505.1 m at resolution 13, and one kilometre at every other one."""
    assert round(G._auto_segment(13), 1) == 505.1
    assert all(G._auto_segment(resolution) == 1000.0 for resolution in range(1, 13))


def test_an_edge_the_count_densifies_leaves_the_lattice() -> None:
    """A square of 2 km, with its vertices every 400 m and with its corners only.

    The count reads a long edge as a geodesic and densifies it with points
    that are not on the lattice, so the same square on the lattice is two
    different figures to it. With a vertex every 400 m no edge reaches the
    threshold and the count is the lattice count. With the corners only it
    is not: it counts the figure bounded by geodesics, and near the points
    it inserted, which are not moved onto anything, its decisions are the
    platform's again. The difference is a property of the reading, so only
    its existence is asserted.
    """
    lattice = _lattice(*ORIGINS["SW"])
    side, step = 200_000, 40_000
    corners: Ring = [(0, 0), (side, 0), (side, side), (0, side)]
    dense: Ring = (
        [(i, 0) for i in range(0, side, step)]
        + [(side, j) for j in range(0, side, step)]
        + [(i, side) for i in range(side, 0, -step)]
        + [(0, j) for j in range(side, 0, -step)]
    )
    reference = _reference(corners, STEP_DIRECTION["SW"])
    assert _reference(dense, STEP_DIRECTION["SW"]) == reference == side * side
    assert _count(lattice, dense) == reference
    assert _count(lattice, corners) != reference


def test_a_vertex_off_the_lattice_is_still_left_to_the_bits() -> None:
    """Two micrometres off the lattice, a vertex is not moved, and the bits decide.

    ``ONE_SHIFT_FROM_THE_EDGE`` is slid two micrometres along its west edge:
    the same line, carrying the same four near ties, but every vertex now
    lies farther from the grid than the tolerance and stays where the
    projection puts it. Under simulated machines the count moves again, as
    it did with the move turned off.
    """
    lattice = _lattice(*ORIGINS["SW"])
    lon_lat = lattice.vertices(ONE_SHIFT_FROM_THE_EDGE)
    plane = [G._snap_vertex(*geodetic_to_sinusoidal(*vertex)) for vertex in lon_lat]
    along = 2e-6 / math.dist(plane[0], plane[1])
    slide_x = along * (plane[1][0] - plane[0][0])
    slide_y = along * (plane[1][1] - plane[0][1])
    slid = [sinusoidal_to_geodetic(x + slide_x, y + slide_y) for x, y in plane]
    for vertex in slid:
        projected = geodetic_to_sinusoidal(*vertex)
        assert G._snap_vertex(*projected) == projected
    counted = [
        itacart.count_internal_cells(Polygon(_moved_by_ulps(slid, seed, FAR_BITS)), RES)
        for seed in range(12)
    ]
    assert len(set(counted)) > 1
    assert max(counted) - min(counted) <= 8
