"""Cells as vertices: the inverse of ``vertex_to_cell``, declared and enumerated.

A vertex written as a cell is a position together with how well it is
known. Reading it back asks one question -- which point of the cell stands
for the vertex -- and this module holds the package to a single answer on
every route, the centroid, before it measures what that answer costs.

**The population is enumerated, not drawn.** Positions sit on a lattice of
latitudes and longitudes in the four quadrants, at five resolutions from
10 km to 1 cm, plus the two extension zones and a lattice over the polar
caps at every resolution. The families this grid breaks rules in -- the
prime-meridian column, the border-absorbing column, the southern equator
row, the polar row, the caps, the extensions -- are counted, and the counts
are asserted, so that a family cannot drop out of the population and leave
its tests passing over nothing.

**Distances are geodesic and the farthest point is a vertex.** A cell's
edges are straight on the sinusoidal plane. Densifying every edge six-fold
on the plane, inverting and measuring moved no reach in the enumerated
population, so the reach of a cell -- the distance from its centroid to
its farthest point -- is taken over its vertices.

**Longitude is the shear, not latitude.** The shape of an ordinary cell
on the ground is set by ``u = |lon * sin(lat)|``, lon in radians, so a
latitude sweep at one longitude covers a fraction of the range. The
lattice spans ``u`` from 0 to pi.

**The reach is public.** :func:`itacart.cell_max_uncertainty` returns it,
and the last part of this module measures that function: against the
definition here, against the law over a population that reaches the polar
approach and ``u`` past pi, and against the boundary between the vertices.
"""

from __future__ import annotations

import math
from functools import lru_cache

import pytest
from shapely.geometry import Point, Polygon

import itacart
from itacart import geodesy, metrics
from itacart.boundary import plane_ring, to_geodetic
from itacart.geodesy import (
    geodetic_to_sinusoidal,
    inverse_geodesic,
    sinusoidal_to_geodetic,
)

LATITUDES = (1e-7, 10.0, 20.0, 30.0, 45.0, 60.0, 70.0, 80.0, 85.0, 89.0, 89.9)
LONGITUDES = (1e-7, 0.5, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0, 179.0, 179.999)
RESOLUTIONS = (1, 2, 5, 9, 13)
SIGNS = ((1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0))

#: Fiji and Chukotka, each written both ways the package reads a position
#: past the antemeridian: west of it, and in the eastern spelling.
EXTENSION_POSITIONS = (
    (-179.5, -16.0),
    (-178.5, -18.5),
    (180.5, -21.0),
    (-179.0, 65.0),
    (-175.0, 68.0),
    (186.0, 71.0),
)

#: The caps begin at about 89.9823 degrees; the lattice runs to the pole.
CAP_LATITUDES = (89.9824, 89.9868, 89.9912, 89.9956, 90.0)
CAP_LONGITUDES = tuple(-180.0 + 45.0 * k for k in range(9))

#: Latitude below which the shear law holds within ORDINARY_LAW_TOLERANCE.
LAW_LATITUDE_LIMIT = 85.0


def _res1_row(cell: str) -> int:
    return int(itacart.split_components(cell)[1].split("/")[1])


def _family(cell: str) -> str:
    row = _res1_row(cell)
    column = int(itacart.split_components(cell)[1].split("/")[0])
    if row == 1000:
        return "cap"
    if row == 999:
        return "polar row"
    if itacart.is_extension_cell(cell):
        return "extension"
    if itacart.absorbs_border(cell):
        return "absorbing"
    if itacart.cell_shape(cell) == "triangle":
        return "meridian triangle"
    if itacart.is_quadrant_boundary_cell(cell):
        return "southern equator row"
    return "ordinary" if column > 0 else "meridian column"


@lru_cache(maxsize=None)
def _population() -> tuple[str, ...]:
    positions = [
        (sx * lon, sy * lat, resolution)
        for resolution in RESOLUTIONS
        for lat in LATITUDES
        for lon in LONGITUDES
        for sx, sy in SIGNS
    ]
    positions += [
        (lon, lat, resolution)
        for resolution in RESOLUTIONS
        for lon, lat in EXTENSION_POSITIONS
    ]
    positions += [
        (lon, sign * lat, resolution)
        for resolution in range(1, 14)
        for lat in CAP_LATITUDES
        for lon in CAP_LONGITUDES
        for sign in (1.0, -1.0)
    ]
    cells = {itacart.geo_to_cell(lon, lat, res) for lon, lat, res in positions}
    return tuple(sorted(cells))


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    if a == b:
        return 0.0
    return inverse_geodesic(a[0], a[1], b[0], b[1])[0]


def _reach(cell: str, point: tuple[float, float]) -> float:
    """Farthest vertex of the cell from ``point``, in metres."""
    ring = itacart.cell_to_boundary(cell)
    return max(_distance(point, (lon, lat)) for lon, lat in ring)


def _shear(cell: str) -> float:
    lon, lat = itacart.cell_to_centroid(cell)
    return abs(math.radians(lon) * math.sin(math.radians(lat)))


def shear_law(u: float) -> float:
    """Reach of an ordinary parallelogram over its side, to first order.

    On the ground the cell's plane vertices ``(0, 0)``, ``(l, 0)``,
    ``(0, l)`` and ``(-l, l)`` sit at ``(0, 0)``, ``(l, 0)``, ``(u l, l)``
    and ``((u - 1) l, l)``, because a plane displacement ``(dx, dy)`` is a
    ground displacement ``(dx + u dy, dy)``. The centroid is at
    ``(u l / 2, l / 2)``, and the farthest vertex gives the law.
    """
    return max(0.5 * math.sqrt(1.0 + u * u), math.sqrt((1.0 - u / 2.0) ** 2 + 0.25))


def anchor_law(u: float) -> float:
    """Reach of the same parallelogram measured from its anchor instead."""
    return max(math.sqrt(1.0 + (1.0 - u) ** 2), math.sqrt(1.0 + u * u))


# --------------------------------------------------------------------------
# The population is what it says it is
# --------------------------------------------------------------------------


def test_the_population_names_every_family_and_says_how_many() -> None:
    """Every family is present, in every quadrant it exists in, and counted."""
    counts: dict[str, int] = {}
    quadrants: dict[str, set[str]] = {}
    for cell in _population():
        family = _family(cell)
        counts[family] = counts.get(family, 0) + 1
        quadrants.setdefault(family, set()).add(itacart.quadrant_of(cell))
    assert counts == {
        "absorbing": 62,
        "cap": 520,
        "extension": 87,
        "meridian column": 56,
        "meridian triangle": 86,
        "ordinary": 1744,
        "polar row": 338,
        "southern equator row": 76,
    }
    for family in ("absorbing", "ordinary", "polar row"):
        assert quadrants[family] == {"NE", "NW", "SE", "SW"}, family
    assert quadrants["cap"] == {"NE", "SE"}
    assert all(itacart.is_valid_cell(cell) for cell in _population())


# --------------------------------------------------------------------------
# One declared point, whatever the route
# --------------------------------------------------------------------------


def test_every_single_point_route_answers_with_the_centroid() -> None:
    """The rebuilt vertex, the H3 alias, the engine and the declaration agree.

    The answer is compared for equality, not closeness: the routes are
    required to be one computation, and a second computation that agreed
    to a nanometre would still be a second answer.
    """
    engine = itacart.ITACaRT()
    assert itacart.describe()["cell_geometry"]["representative_position"] == (
        "centroid"
    )
    for cell in _population():
        centroid = itacart.cell_to_centroid(cell)
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        assert (rebuilt.x, rebuilt.y) == centroid, cell
        assert itacart.cell_to_latlng(cell) == (centroid[1], centroid[0]), cell
        assert engine.cell_to_centroid(cell) == centroid, cell


def test_a_rebuilt_ring_is_the_centroids_of_its_cells_in_order() -> None:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    cells = itacart.vertex_to_cell(Polygon(ring), 13)
    rebuilt = itacart.cells_to_geometry(cells, "Polygon")
    assert list(rebuilt.exterior.coords)[:-1] == [
        itacart.cell_to_centroid(cell) for cell in cells
    ]


def test_the_anchor_is_on_the_boundary_and_the_rebuilt_vertex_is_not() -> None:
    """Why the anchor is not offered: a vertex is not a point within the cell."""
    for cell in ("NE(0500/0300(3(C2)))", "SW(0476/0260(4(A2(2(E3)))))"):
        ring = plane_ring(cell)[1]
        polygon = Polygon([(x - ring[0][0], y - ring[0][1]) for x, y in ring])
        anchor = geodetic_to_sinusoidal(*itacart.cell_to_anchor(cell))
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        centre = geodetic_to_sinusoidal(rebuilt.x, rebuilt.y)
        shift = (-ring[0][0], -ring[0][1])
        on_edge = Point(anchor[0] + shift[0], anchor[1] + shift[1])
        inside = Point(centre[0] + shift[0], centre[1] + shift[1])
        assert polygon.exterior.distance(on_edge) < 1e-6
        assert polygon.contains(inside)


# --------------------------------------------------------------------------
# The round trip, enumerated
# --------------------------------------------------------------------------


def test_every_rebuilt_vertex_quantizes_back_to_the_cell_that_encoded_it() -> None:
    """Over the whole population, extension zones and caps included.

    The two families this could not claim before are the two the
    enumeration found: a cap position quantized to a spelling that names
    no cell, and a rebuilt extension vertex, written past 180 degrees as
    the package writes it, was refused by the quantizer.
    """
    for cell in _population():
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        resolution = itacart.get_resolution(cell)
        assert itacart.geo_to_cell(rebuilt.x, rebuilt.y, resolution) == cell, cell


def _probes(cell: str, fraction: float) -> list[tuple[float, float]]:
    """Positions just inside each vertex, moved toward the centroid on the plane."""
    centre = geodetic_to_sinusoidal(*itacart.cell_to_centroid(cell))
    probes = []
    for lon, lat in itacart.cell_to_boundary(cell):
        x, y = geodetic_to_sinusoidal(lon, lat)
        probes.append(
            sinusoidal_to_geodetic(
                x + fraction * (centre[0] - x), y + fraction * (centre[1] - y)
            )
        )
    return probes


def test_a_vertex_comes_back_within_its_cells_reach_and_near_it_at_a_corner() -> None:
    """The reach is a bound the round trip reaches, not a loose ceiling.

    Every cell of the ordinary and meridian families is probed just inside
    each of its vertices, one per cent of the way to the centroid. Each
    probe goes through ``vertex_to_cell`` and ``cells_to_geometry``; it must
    come back to the same cell, within the cell's reach, and the probe at
    the farthest corner must come back at least 98 per cent of the reach
    away, which is as close as a probe one per cent inside can get.
    """
    probed = 0
    for cell in _population():
        if _family(cell) not in ("ordinary", "meridian column", "meridian triangle"):
            continue
        resolution = itacart.get_resolution(cell)
        centroid = itacart.cell_to_centroid(cell)
        reach = _reach(cell, centroid)
        farthest = 0.0
        for lon, lat in _probes(cell, 0.01):
            cells = itacart.vertex_to_cell(Point(lon, lat), resolution)
            assert cells == [cell], (cell, lon, lat, cells)
            back = itacart.cells_to_geometry(cells, "Point")
            error = _distance((lon, lat), (back.x, back.y))
            assert error <= reach, cell
            farthest = max(farthest, error)
        assert farthest >= 0.98 * reach, cell
        probed += 1
    assert probed == 1886


def test_the_reach_of_an_ordinary_cell_follows_the_shear_law() -> None:
    """Measured reach over side against the closed form, below 85 degrees.

    The law is first order in the side. Closer to the pole a finite cell
    spans enough of the meridian's convergence to exceed it, and the
    population above the limit is asserted to do so, so that the limit is
    earned rather than conservative.
    """
    below = []
    above = []
    for cell in _population():
        if _family(cell) != "ordinary":
            continue
        ratio = _reach(cell, itacart.cell_to_centroid(cell)) / itacart.cell_size(
            itacart.get_resolution(cell)
        )
        law = shear_law(_shear(cell))
        lat = abs(itacart.cell_to_centroid(cell)[1])
        (below if lat < LAW_LATITUDE_LIMIT else above).append(ratio / law)
    assert (len(below), len(above)) == (1456, 288)
    assert min(below) > 1.0 - 1e-3
    assert max(below) < 1.0 + 3e-3
    assert 1.0 + 3e-3 < max(above) < 1.0 + 2.1e-2


def test_the_shear_law_spans_its_range_in_the_population() -> None:
    """Instrument check: the lattice reaches both ends of the law."""
    shears = [_shear(cell) for cell in _population() if _family(cell) == "ordinary"]
    assert min(shears) < 1e-3
    assert max(shears) > 3.0
    assert shear_law(0.0) == pytest.approx(math.sqrt(5.0) / 2.0)
    assert shear_law(1.0) == pytest.approx(math.sqrt(2.0) / 2.0)
    assert shear_law(math.pi) == pytest.approx(math.sqrt(1.0 + math.pi**2) / 2.0)


def test_a_meridian_triangle_reaches_root_ten_over_three_sides() -> None:
    """The triangle has base two sides and height one; its centroid is a third up.

    A finite triangle departs from the value for the reason the
    parallelogram departs from its law, its side not being small against
    the distance to the pole: within 0.3 per cent below 85 degrees, and up
    to 1.4 per cent above.
    """
    triangles = [c for c in _population() if _family(c) == "meridian triangle"]
    ratios = {"below": [], "above": []}
    for cell in triangles:
        side = itacart.cell_size(itacart.get_resolution(cell))
        reach = _reach(cell, itacart.cell_to_centroid(cell)) / side
        lat = abs(itacart.cell_to_centroid(cell)[1])
        band = "below" if lat < LAW_LATITUDE_LIMIT else "above"
        ratios[band].append(reach / (math.sqrt(10.0) / 3.0))
    assert (len(ratios["below"]), len(ratios["above"])) == (70, 16)
    assert all(abs(ratio - 1.0) < 3e-3 for ratio in ratios["below"])
    assert all(abs(ratio - 1.0) < 1.5e-2 for ratio in ratios["above"])
    assert max(abs(ratio - 1.0) for ratio in ratios["above"]) > 3e-3


@pytest.mark.slow
def test_the_border_families_are_not_bounded_by_the_law_and_are_measured() -> None:
    """The named exemption, earned and pinned.

    Population: the last existing column of every resolution-1 row in the
    four quadrants, and the border, polar-row, cap and extension cells of
    the enumerated population above. Some exceed the shear law, which is what makes the
    exemption necessary; the largest reach is pinned by value, so that a
    change in either direction is seen.
    """
    side = itacart.cell_size(1)
    border = set(c for c in _population() if _family(c) in _BORDER_FAMILIES)
    for quadrant in ("NE", "NW", "SE", "SW"):
        for row in range(1001):
            last = itacart.last_lattice_column(quadrant, row, side)
            cell = f"{quadrant}({max(last, 0):04d}/{row:04d})"
            if itacart.is_valid_cell(cell) and _family(cell) in _BORDER_FAMILIES:
                border.add(cell)
    exceeding = 0
    largest = 0.0
    for cell in border:
        ratio = _reach(cell, itacart.cell_to_centroid(cell)) / itacart.cell_size(
            itacart.get_resolution(cell)
        )
        largest = max(largest, ratio)
        if ratio > shear_law(_shear(cell)) * (1.0 + 5e-3):
            exceeding += 1
    assert (len(border), exceeding, round(largest, 4)) == (4967, 1969, 2.0532)
    assert largest < 2.06


_BORDER_FAMILIES = ("absorbing", "polar row", "cap", "extension")


# --------------------------------------------------------------------------
# The single point of the opening measurement is a position, not a law
# --------------------------------------------------------------------------


SAO_PAULO = (-46.633, -23.55)


def test_the_worst_case_favours_the_centroid_in_every_ordinary_cell() -> None:
    """Measured against the anchor, the centroid's reach is never larger.

    The ratio follows the two laws: 1.26 on the prime meridian, a minimum
    of ``sqrt(5) - 1`` between, and exactly 2 from ``u = 1``, where the
    anchor becomes an end of the long diagonal.
    """
    ratios = []
    for cell in _population():
        if _family(cell) != "ordinary":
            continue
        centroid = _reach(cell, itacart.cell_to_centroid(cell))
        anchor = _reach(cell, itacart.cell_to_anchor(cell))
        ratios.append(anchor / centroid)
    floor = min(anchor_law(k / 1000.0) / shear_law(k / 1000.0) for k in range(3200))
    assert floor == pytest.approx(math.sqrt(5.0) - 1.0, abs=1e-6)
    assert min(ratios) > floor * (1.0 - 1e-3)
    assert max(ratios) < 2.0 + 1e-3
    assert anchor_law(0.0) / shear_law(0.0) == pytest.approx(2.0 * math.sqrt(0.4))
    assert anchor_law(2.0) / shear_law(2.0) == pytest.approx(2.0)


def test_the_anchor_wins_a_single_point_only_inside_its_own_corner() -> None:
    """The opening table's reversal, reproduced and then taken apart.

    At resolutions 5 and 9 the São Paulo point is closer to the centroid;
    at 11 and 13 to the anchor. Moving the same point to the far corner of
    its resolution-13 cell, without changing the resolution, reverses it:
    what decided the table was where the point sits in its cell.
    """
    lon, lat = SAO_PAULO
    winners = {}
    for resolution in (5, 9, 11, 13):
        cell = itacart.geo_to_cell(lon, lat, resolution)
        to_anchor = _distance((lon, lat), itacart.cell_to_anchor(cell))
        to_centroid = _distance((lon, lat), itacart.cell_to_centroid(cell))
        winners[resolution] = "anchor" if to_anchor < to_centroid else "centroid"
    assert winners == {5: "centroid", 9: "centroid", 11: "anchor", 13: "anchor"}

    cell = itacart.geo_to_cell(lon, lat, 13)
    probes = _probes(cell, 0.01)
    anchor = itacart.cell_to_anchor(cell)
    far = max(probes, key=lambda probe: _distance(probe, anchor))
    assert itacart.geo_to_cell(far[0], far[1], 13) == cell
    assert _distance(far, itacart.cell_to_centroid(cell)) < _distance(far, anchor)


# --------------------------------------------------------------------------
# Precision is the index's to declare
# --------------------------------------------------------------------------


@pytest.mark.parametrize("lat", [-23.55, 0.5, 45.0, 70.0])
@pytest.mark.parametrize("sx, sy", SIGNS)
def test_resolution_13_neighbours_survive_only_at_full_precision(
    lat: float, sx: float, sy: float
) -> None:
    """Six decimal places merge neighbours; the rebuilt coordinates do not.

    Each disk is a resolution-13 cell and its ring of neighbours. Rebuilt,
    every vertex is distinct and quantizes back to its own cell. Rounded to
    six decimal places, two neighbours become one point, and requantizing
    the rounded coordinates no longer returns the cells.
    """
    centre = itacart.geo_to_cell(sx * 46.633, sy * abs(lat), 13)
    disk = sorted(itacart.grid_disk(centre, 1))
    assert len(disk) == 9
    rebuilt = itacart.cells_to_geometry(disk, "MultiPoint")
    points = [(point.x, point.y) for point in rebuilt.geoms]
    assert len(set(points)) == len(points)
    assert [itacart.geo_to_cell(x, y, 13) for x, y in points] == disk

    rounded = [(round(x, 6), round(y, 6)) for x, y in points]
    assert len(set(rounded)) < len(rounded)
    assert [itacart.geo_to_cell(x, y, 13) for x, y in rounded] != disk


# --------------------------------------------------------------------------
# The reach, made public
# --------------------------------------------------------------------------

#: How far a regular parallelogram's radius may depart from the shear law,
#: as a fraction of the law, and its diameter from twice the radius. The
#: departure is first order in the side and grows towards the pole, so the
#: bound is a slope per band of latitude, times the side in units of 10 km,
#: plus a floor for the rounding of the coordinates, which is what sets the
#: departure below 85 degrees from resolution 11 down. These are the figures
#: the documentation of ``cell_max_uncertainty`` states. Measured over the
#: population below: slopes of 4.07e-3 under 85 degrees and 0.190 from 85
#: degrees to row 998 for the radius, 4.10e-3 and 0.191 for the diameter,
#: and 6.4e-7 where the floor decides.
LAW_SLOPES = {"below 85": 4.2e-3, "from 85": 0.20}
LAW_FLOOR = 1e-6

#: The middle of each resolution-1 row from 990 to 998, the last rows of
#: regular cells before the polar row, where the law departs most.
POLAR_APPROACH = tuple(
    sinusoidal_to_geodetic(0.0, (row + 0.5) * 1e4)[1] for row in range(990, 999)
)

OUTPUTS = ("radius", "diameter", "ratio")


def _band(latitude: float) -> str:
    return "below 85" if abs(latitude) < LAW_LATITUDE_LIMIT else "from 85"


def _law_bound(resolution: int, latitude: float) -> float:
    side = itacart.cell_size(resolution) / itacart.cell_size(1)
    return LAW_SLOPES[_band(latitude)] * side + LAW_FLOOR


def _regular(cell: str) -> bool:
    """A parallelogram of the lattice: no triangle, no border, no polar row."""
    return (
        _res1_row(cell) < 999
        and itacart.cell_shape(cell) == "parallelogram"
        and not itacart.absorbs_border(cell)
    )


def _row_cells(quadrant: str, row: int) -> list[str]:
    """Every resolution-1 cell of a row, those past the antemeridian included."""
    last = itacart.last_lattice_column(quadrant, row, itacart.cell_size(1))
    cells = (f"{quadrant}({column:04d}/{row:04d})" for column in range(last + 1))
    return [cell for cell in cells if itacart.is_valid_cell(cell)]


def _law_positions(resolution: int) -> list[tuple[float, float]]:
    """A lattice in latitude and ``u``, in four quadrants, and Chukotka.

    The latitudes sit on both sides of 85 degrees and run through the
    polar approach to half a cell below the top of row 998. At each one
    ``u`` runs from 0 to just short of the antimeridian, where it is
    largest for that latitude. Chukotka, east of the antimeridian, carries
    ``u`` past pi.
    """
    side = itacart.cell_size(resolution)
    top = sinusoidal_to_geodetic(0.0, 9.99e6 - side / 2.0)[1]
    positions = []
    for lat in (0.5, 30.0, 60.0, 84.99, 85.01, *POLAR_APPROACH, top):
        sine = math.sin(math.radians(lat))
        far = math.pi * sine
        for u in (0.0005, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, far - 0.02, far - 0.005):
            lon = math.degrees(u / sine)
            if u < far and lon < 179.999:
                positions += [(sx * lon, sy * lat) for sx, sy in SIGNS]
    chukotka = [
        (lon, lat)
        for lat in (64.2, 66.5, 69.5, 71.9)
        for lon in (180.5, 184.0, 188.0, 189.5, 190.3)
    ]
    return positions + chukotka


#: Cells in the law population, all thirteen resolutions together.
LAW_POPULATION_SIZE = 7073


@lru_cache(maxsize=None)
def _law_population() -> tuple[str, ...]:
    """Where the law is tested: regular parallelograms, enumerated by construction.

    At resolution 1, every regular cell of row 943, the last one below 85
    degrees, and of rows 990 to 998, in the four quadrants, and the cells
    past the antimeridian in three rows of Chukotka. From resolution 2 to
    13, the regular cells that hold the positions of the lattice above.
    """
    cells: set[str] = set()
    for quadrant in ("NE", "NW", "SE", "SW"):
        for row in (943, *range(990, 999)):
            cells.update(cell for cell in _row_cells(quadrant, row) if _regular(cell))
    for row in (709, 754, 799):
        cells.update(
            cell
            for cell in _row_cells("NE", row)
            if _regular(cell) and itacart.cell_to_centroid(cell)[0] > 180.0
        )
    for resolution in range(2, 14):
        for lon, lat in _law_positions(resolution):
            cell = itacart.geo_to_cell(lon, lat, resolution)
            assert itacart.is_valid_cell(cell), (lon, lat, resolution)
            if _regular(cell):
                cells.add(cell)
    return tuple(sorted(cells))


@lru_cache(maxsize=None)
def _outputs(cell: str) -> tuple[float, float, float]:
    """Radius, diameter and ratio of a cell, from the public function, once."""
    radius, diameter, ratio = (
        itacart.cell_max_uncertainty(cell, output) for output in OUTPUTS
    )
    return radius, diameter, ratio


def _holds_a_pole(cell: str) -> bool:
    return any(abs(lat) == 90.0 for _, lat in itacart.cell_to_boundary(cell))


@lru_cache(maxsize=None)
def _law_measures() -> tuple[tuple[str, int, float, float, float, float], ...]:
    """Cell, resolution, latitude, ``u``, radius and diameter, measured once."""
    rows = []
    for cell in _law_population():
        lon, lat = itacart.cell_to_centroid(cell)
        rows.append(
            (
                cell,
                itacart.get_resolution(cell),
                lat,
                abs(math.radians(lon) * math.sin(math.radians(lat))),
                itacart.cell_max_uncertainty(cell),
                itacart.cell_max_uncertainty(cell, "diameter"),
            )
        )
    return tuple(rows)


@pytest.mark.slow
def test_the_uncertainty_is_the_reach_of_the_centroid_bit_for_bit() -> None:
    """The public name is the definition this module measures, to the bit.

    The radius is the farthest vertex from the centroid, over this
    module's population and the law population below; the diameter is the
    farthest pair of vertices and the ratio the radius over the side of
    the resolution, over this module's population.
    """
    for cell in _population():
        radius, diameter, ratio = _outputs(cell)
        assert radius == _reach(cell, itacart.cell_to_centroid(cell)), cell
        ring = itacart.cell_to_boundary(cell)
        assert diameter == max(
            _distance(first, second)
            for k, first in enumerate(ring)
            for second in ring[k + 1 :]
        ), cell
        assert ratio == radius / itacart.cell_size(itacart.get_resolution(cell)), cell
    for cell, _, _, _, radius, _ in _law_measures():
        assert radius == _reach(cell, itacart.cell_to_centroid(cell)), cell


@pytest.mark.slow
def test_the_law_population_reaches_what_it_is_named_for() -> None:
    """Both bands, the polar approach and ``u`` past pi, at every resolution."""
    rows: dict[int, list[tuple[str, float, float]]] = {}
    for cell, resolution, latitude, u, _, _ in _law_measures():
        rows.setdefault(resolution, []).append((cell, latitude, u))
    assert sorted(rows) == list(range(1, 14))
    for resolution, members in rows.items():
        assert {_band(latitude) for _, latitude, _ in members} == {
            "below 85",
            "from 85",
        }
        assert {_res1_row(cell) for cell, _, _ in members} >= set(range(990, 999))
        assert max(u for _, _, u in members) > math.pi, resolution
        assert {cell[:2] for cell, _, _ in members} == {"NE", "NW", "SE", "SW"}
    assert len(_law_population()) == LAW_POPULATION_SIZE


@pytest.mark.slow
def test_the_uncertainty_follows_the_shear_law_by_resolution_and_band() -> None:
    """Every regular parallelogram of the population, against the law at its own u.

    The radius never falls short of the law by more than the rounding
    floor. The bound is earned in both bands: where the side and not the
    rounding sets the departure, at resolutions 1 to 10, the measured slope
    comes within ten per cent of the one the bound allows.
    """
    slopes = {"below 85": 0.0, "from 85": 0.0}
    for cell, resolution, latitude, u, radius, _ in _law_measures():
        side = itacart.cell_size(resolution)
        departure = radius / (side * shear_law(u)) - 1.0
        assert abs(departure) <= _law_bound(resolution, latitude), cell
        assert departure >= -LAW_FLOOR, cell
        if resolution <= 10:
            slope = abs(departure) * itacart.cell_size(1) / side
            slopes[_band(latitude)] = max(slopes[_band(latitude)], slope)
    assert 4.0e-3 < slopes["below 85"] < LAW_SLOPES["below 85"]
    assert 0.18 < slopes["from 85"] < LAW_SLOPES["from 85"]


#: How far a meridian triangle's reach may fall short of ``sqrt(10) / 3 * l``,
#: as a fraction of it: a slope per band of latitude times the side in units
#: of 10 km, as for the parallelogram. Measured over the population below:
#: slopes of 2.68e-3 under 85 degrees and 0.125 from 85 degrees to row 998,
#: and never past the value by more than the rounding floor.
TRIANGLE_SLOPES = {"below 85": 2.7e-3, "from 85": 0.13}

#: Triangles in that population, all thirteen resolutions together.
TRIANGLE_POPULATION_SIZE = 494


@lru_cache(maxsize=None)
def _triangle_population() -> tuple[str, ...]:
    """Meridian triangles from row 930 to row 998, and down the meridian.

    At resolution 1, the triangle of every row from 930 to 998 in the two
    eastern quadrants. From resolution 2 to 13, the triangle that holds a
    position just east of the prime meridian at each latitude of the law
    lattice, in both hemispheres.
    """
    cells = {f"{q}(0000/{row:04d})" for q in ("NE", "SE") for row in range(930, 999)}
    for resolution in range(2, 14):
        side = itacart.cell_size(resolution)
        top = sinusoidal_to_geodetic(0.0, 9.99e6 - side / 2.0)[1]
        for lat in (0.5, 30.0, 60.0, 84.99, 85.01, *POLAR_APPROACH, top):
            for sign in (1.0, -1.0):
                cells.add(itacart.geo_to_cell(1e-9, sign * lat, resolution))
    assert all(itacart.cell_shape(cell) == "triangle" for cell in cells)
    return tuple(sorted(cells))


@pytest.mark.slow
def test_a_meridian_triangle_falls_short_of_its_value_near_the_pole() -> None:
    """The triangle's reach over a population that reaches row 998.

    ``sqrt(10) / 3 * l`` is first order in the side. A finite triangle
    falls short of it by a slope per band of latitude times the side, and
    never passes it by more than the rounding floor. The slopes are earned
    as the law's are, within ten per cent of the bound.
    """
    target = math.sqrt(10.0) / 3.0
    slopes = {"below 85": 0.0, "from 85": 0.0}
    population = _triangle_population()
    assert len(population) == TRIANGLE_POPULATION_SIZE
    for cell in population:
        resolution = itacart.get_resolution(cell)
        band = _band(itacart.cell_to_centroid(cell)[1])
        side = itacart.cell_size(resolution) / itacart.cell_size(1)
        departure = itacart.cell_max_uncertainty(cell, "ratio") / target - 1.0
        assert departure <= LAW_FLOOR, cell
        assert -departure <= TRIANGLE_SLOPES[band] * side + LAW_FLOOR, cell
        if resolution <= 10:
            slopes[band] = max(slopes[band], -departure / side)
    assert 2.5e-3 < slopes["below 85"] < TRIANGLE_SLOPES["below 85"]
    assert 0.12 < slopes["from 85"] < TRIANGLE_SLOPES["from 85"]


#: A resolution-13 cell at 30 degrees where ``u`` is one, and the line from
#: its centroid to its farthest vertex. The reach along it with this
#: release's solver, and with the solver of release 1.0.0.
SQUARE_13 = "NE(1105/0332(2(A2(2(B2(1(B3(2(D2(2(D2(4(E3)))))))))))))"
SQUARE_13_LINE = (
    114.59155901034882,
    29.999999973475177,
    114.59155906216967,
    30.00000001858018,
)
SQUARE_13_REACH = 0.007071068784917265
SQUARE_13_REACH_RELEASE_1_0_0 = 0.0070710463726038645


@pytest.mark.slow
def test_the_law_check_fails_with_the_solver_of_release_1_0_0(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control: restore the old stop rule and resolutions 9 to 13 fail.

    Release 1.0.0 stopped the Vincenty iteration once its longitude
    iterate moved by less than ``tol``, whatever the length of the line.
    This release stops below ``tol * min(1, sigma)``, or at four units in
    the last place of the iterate. Taking ``min(1, sigma)`` as one restores
    the old rule exactly: four units in the last place of any iterate here
    lie far below ``tol``, so the floor never stops the iteration first.
    The pinned line says the restoration is the old solver and not an
    approximation of it, to the sixteen units in the last place a Vincenty
    result is held to across platforms.
    """
    # Measured before the patch: the population's figures are cached for
    # the tests that read them, and must be this release's.
    measures = _law_measures()
    square = itacart.cell_max_uncertainty(SQUARE_13)
    assert square == pytest.approx(SQUARE_13_REACH, abs=16 * math.ulp(square))
    line = inverse_geodesic(*SQUARE_13_LINE)[0]
    assert line == pytest.approx(SQUARE_13_REACH, abs=16 * math.ulp(line))
    monkeypatch.setattr(geodesy, "min", lambda *_: 1.0, raising=False)
    old = inverse_geodesic(*SQUARE_13_LINE)[0]
    assert old == pytest.approx(SQUARE_13_REACH_RELEASE_1_0_0, abs=16 * math.ulp(old))
    for resolution in range(9, 14):
        failures = 0
        for cell, level, latitude, u, _, _ in measures:
            if level != resolution:
                continue
            side = itacart.cell_size(resolution)
            radius = itacart.cell_max_uncertainty(cell)
            departure = radius / (side * shear_law(u)) - 1.0
            failures += abs(departure) > _law_bound(resolution, latitude)
        assert failures > 0, resolution


@pytest.mark.slow
def test_in_a_parallelogram_the_diameter_is_twice_the_radius_within_the_bound() -> None:
    """To the same order as the law, and over the same population."""
    for cell, resolution, latitude, _, radius, diameter in _law_measures():
        departure = diameter / (2.0 * radius) - 1.0
        assert abs(departure) <= _law_bound(resolution, latitude), cell


def test_outside_the_parallelogram_the_diameter_is_not_tied_to_the_radius() -> None:
    """A meridian triangle, a border trapezoid and a cap, measured.

    The centroid of a parallelogram is the midpoint of both diagonals, so
    its farthest vertex is half the longest diagonal away. No other shape
    has that point, and no relation is claimed for them.
    """
    measured = {
        cell: round(
            itacart.cell_max_uncertainty(cell, "diameter")
            / (2.0 * itacart.cell_max_uncertainty(cell)),
            3,
        )
        for cell in ("NE(0000/0300)", "NE(2003/0000)", "NE(0000/1000)")
    }
    assert measured == {
        "NE(0000/0300)": 0.949,
        "NE(2003/0000)": 0.859,
        "NE(0000/1000)": 0.3,
    }


def test_where_u_is_one_the_cell_is_a_square_on_the_ground() -> None:
    """Radius ``l / sqrt(2)`` and diameter ``sqrt(2) l``, from real cells.

    The position where ``u`` is one is quantized at 30 and at 60 degrees,
    at four resolutions, and each cell is compared with the law at its
    own ``u``, which its size moves off one. The cell the specification
    calls a square, on row 900, is one of them too.
    """
    for lat in (30.0, 60.0):
        lon = math.degrees(1.0 / math.sin(math.radians(lat)))
        for resolution in (1, 5, 9, 13):
            cell = itacart.geo_to_cell(lon, lat, resolution)
            side = itacart.cell_size(resolution)
            law = shear_law(_shear(cell))
            bound = _law_bound(resolution, lat)
            ratio = itacart.cell_max_uncertainty(cell, "ratio")
            diameter = itacart.cell_max_uncertainty(cell, "diameter") / side
            assert ratio == pytest.approx(law, rel=bound), cell
            assert diameter == pytest.approx(2.0 * law, rel=bound), cell
            assert ratio == pytest.approx(math.sqrt(0.5), abs=2e-3), cell
            assert diameter == pytest.approx(math.sqrt(2.0), abs=4e-3), cell
    assert shear_law(1.0) == math.sqrt(2.0) / 2.0
    row_900 = "NE(0100/0900)"
    assert round(_shear(row_900), 4) == 0.9949
    assert round(itacart.cell_max_uncertainty(row_900, "ratio"), 6) == 0.709775


@pytest.mark.parametrize("resolution", [1, 7, 13])
def test_at_one_resolution_the_uncertainty_follows_the_position(
    resolution: int,
) -> None:
    """Same resolution, seven places along one parallel: the ratio follows ``u``.

    Enumerated in ``u`` and not in latitude, at 80 degrees, where the
    parallel reaches ``u`` of 3.09 before the antimeridian. The first
    place is a little east of the meridian triangle, so that the cell is a
    parallelogram at every resolution. The ratio runs from 1.11 there down
    to ``1 / sqrt(2)`` at ``u = 1`` and up past 1.58 near the antimeridian.
    """
    lat = 80.0
    sine = math.sin(math.radians(lat))
    ratios = {}
    for u in (0.02, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        cell = itacart.geo_to_cell(math.degrees(u / sine), lat, resolution)
        assert _regular(cell), cell
        ratio = itacart.cell_max_uncertainty(cell, "ratio")
        bound = _law_bound(resolution, lat)
        assert ratio == pytest.approx(shear_law(_shear(cell)), rel=bound), cell
        ratios[u] = ratio
    assert ratios[0.02] > 1.10
    assert min(ratios.values()) == ratios[1.0]
    assert ratios[1.0] == pytest.approx(math.sqrt(0.5), abs=2e-3)
    assert ratios[3.0] > 1.58


def _uncertainty_family(cell: str) -> str:
    """The family a cell's uncertainty is read in, from the cell itself.

    The extension is read off the centroid, east of 180 degrees on the
    unwrapped longitude, and a mirrored cell is one spelled under an
    eastern quadrant that leans the western way.
    """
    ring = itacart.cell_to_boundary(cell)
    longitudes = [lon for lon, _ in ring]
    if any(abs(lat) == 90.0 for _, lat in ring):
        return "pole"
    if max(longitudes) - min(longitudes) >= 360.0 - 4 * math.ulp(360.0):
        return "wrapping"
    row = _res1_row(cell)
    if row == 1000:
        return "cap"
    if row == 999:
        return "polar row"
    if itacart.absorbs_border(cell):
        return "absorbing"
    if itacart.cell_shape(cell) == "triangle":
        return "meridian triangle"
    if abs(itacart.cell_to_centroid(cell)[0]) > 180.0:
        return "extension"
    anchor = itacart.cell_to_anchor(cell)
    if itacart.quadrant_of(cell)[1] == "E" and metrics._chirality(cell, anchor) < 0:
        return "mirrored"
    return "ordinary"


def _cap_spellings(cap: str, resolution: int) -> list[str]:
    """Every spelling under a cap at one resolution that names a cell."""
    stem = itacart.split_components(cap)
    paths: list[list[str]] = [[]]
    for level in range(2, resolution + 1):
        paths = [
            path + [code]
            for path in paths
            for code in itacart.refinement_alphabet(level)
        ]
    spellings = (itacart.join_components(stem + path) for path in paths)
    return [spelling for spelling in spellings if itacart.is_valid_cell(spelling)]


def _under_the_pole(resolution: int) -> list[str]:
    """The cells spelled under the resolution-10 cell that holds each pole.

    Built rather than reached: the hierarchy does not descend from that
    cell, so its children are spelled code by code and kept where the
    predicate accepts them.
    """
    cells = []
    for lat in (90.0, -90.0):
        stem = itacart.split_components(itacart.geo_to_cell(0.0, lat, resolution - 1))
        for code in itacart.refinement_alphabet(resolution):
            spelling = itacart.join_components(stem + [code])
            if itacart.is_valid_cell(spelling):
                cells.append(spelling)
    return cells


@lru_cache(maxsize=None)
def _family_population() -> tuple[str, ...]:
    """This module's population, the caps' subtrees to 4, and under the pole."""
    cells = set(_population())
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        for resolution in (2, 3, 4):
            cells.update(_cap_spellings(cap, resolution))
    cells.update(_under_the_pole(11))
    return tuple(sorted(cells))


#: Per family: how many cells, and the smallest and largest ratio.
FAMILY_SUMMARY = {
    "absorbing": (94, 0.6422, 2.0029),
    "cap": (644, 0.784, 5.0881),
    "extension": (50, 0.7089, 1.6363),
    "meridian triangle": (86, 1.0396, 1.0541),
    "mirrored": (28, 1.1141, 1.118),
    "ordinary": (1853, 0.7124, 1.6482),
    "polar row": (338, 0.7234, 1.9294),
    "pole": (26, 0.2431, 1.6095),
    "wrapping": (36, 1.4882, 3.3415),
}


@pytest.mark.slow
def test_every_family_answers_with_a_finite_positive_number() -> None:
    """No family is refused, and each is counted with its extremes pinned.

    The population is this module's, the whole subtree of each cap down to
    resolution 4, with the cells around the pole whose rings wrap the
    parallel, and the cells spelled under the resolution-10 cell that holds
    each pole, which include the two trapezoids at resolution 11 that wrap
    it. The extremes are of the ratio.
    """
    ratios: dict[str, list[float]] = {}
    for cell in _family_population():
        values = _outputs(cell)
        assert all(math.isfinite(value) and value > 0.0 for value in values), cell
        ratios.setdefault(_uncertainty_family(cell), []).append(values[2])
    summary = {
        family: (len(values), round(min(values), 4), round(max(values), 4))
        for family, values in ratios.items()
    }
    assert summary == FAMILY_SUMMARY


def _contour(cell: str, pieces: int) -> tuple[float, float]:
    """Radius and diameter over the boundary, each edge cut into ``pieces``.

    The edges are straight on the sinusoidal plane, so they are cut there
    and every point is inverted, the vertices by the same function as the
    ring the package returns.
    """
    ring = plane_ring(cell)[1]
    points = [
        to_geodetic(x1 + k / pieces * (x2 - x1), y1 + k / pieces * (y2 - y1))
        for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1])
        for k in range(pieces)
    ]
    centroid = itacart.cell_to_centroid(cell)
    radius = max(_distance(centroid, point) for point in points)
    diameter = max(
        _distance(first, second)
        for k, first in enumerate(points)
        for second in points[k + 1 :]
    )
    return radius, diameter


#: Diameter over the vertices of the cell that holds a pole, over the side.
POLE_DIAMETERS = {1: 0.196573, 6: 0.314586, 8: 0.145863, 11: 0.293127, 13: 0.931272}


@pytest.mark.slow
def test_the_cell_holding_a_pole_measures_half_its_diameter_at_the_vertices() -> None:
    """The limitation the definition declares, at every resolution and pole.

    The base of that triangle is a whole parallel, and its two base
    vertices, at -180 and +180 degrees, are one point. The diameter over
    the vertices is therefore the pole to the base, while the boundary
    holds two points of the base half a turn apart: exactly twice as far,
    across the pole. The radius is not affected.
    """
    for lat in (90.0, -90.0):
        for resolution in range(1, 14):
            cell = itacart.geo_to_cell(0.0, lat, resolution)
            assert _uncertainty_family(cell) == "pole"
            radius = itacart.cell_max_uncertainty(cell)
            diameter = itacart.cell_max_uncertainty(cell, "diameter")
            on_boundary = _contour(cell, 32)
            assert on_boundary[0] == radius, cell
            assert diameter / on_boundary[1] == pytest.approx(0.5, abs=1e-7), cell
            if resolution in POLE_DIAMETERS:
                side = itacart.cell_size(resolution)
                assert round(diameter / side, 6) == POLE_DIAMETERS[resolution]


#: The cells of the next test: this module's population and the caps'
#: refinement to resolution 5, less the cell at each pole.
AWAY_FROM_THE_POLE = 5433


@pytest.mark.slow
def test_away_from_the_pole_the_vertices_hold_both_maxima() -> None:
    """Every family of this module, and the caps' physical refinement.

    The cells of the population other than the one at each pole, and every
    cell the hierarchy reaches under the two caps down to resolution 5.
    Cutting every edge into four on the plane moves neither the radius nor
    the diameter of any of them.
    """
    cells = list(_population())
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        for resolution in range(2, 6):
            cells += itacart.get_descendants(cap, resolution)
    cells = [cell for cell in cells if not _holds_a_pole(cell)]
    assert len(cells) == AWAY_FROM_THE_POLE
    for cell in cells:
        assert _contour(cell, 4) == _outputs(cell)[:2], cell


#: Per resolution: the cap spellings whose boundary reaches past their
#: vertices, and the largest ratio of boundary to vertices, in radius and
#: in diameter. One spelling in seven at resolution 5.
CAP_SPELLINGS_PAST_THEIR_VERTICES = {
    3: (8, 1.091, 1.178),
    4: (64, 1.156, 1.825),
    5: (316, 2.648, 3.916),
}


@pytest.mark.slow
def test_a_cap_spelling_outside_the_refinement_can_reach_past_its_vertices() -> None:
    """The other place the vertices fall short, measured and pinned.

    Every spelling under a cap that the predicate accepts at resolutions 3
    and 4, and one in seven at resolution 5, in lexical order. Those the
    hierarchy reaches agree with their boundary, and some of the others do
    not. Per resolution: how many reach past their vertices, and by how
    much at most, in radius and in diameter, with every edge cut into
    eight. Cutting finds a point of the boundary and can miss a farther
    one, so the counts and the ratios are what eight pieces see: lower
    bounds of what the continuous boundary reaches.
    """
    physical = {
        cell
        for cap in ("NE(0000/1000)", "SE(0000/1000)")
        for resolution in (3, 4, 5)
        for cell in itacart.get_descendants(cap, resolution)
    }
    found: dict[int, tuple[int, float, float]] = {}
    for resolution, step in ((3, 1), (4, 1), (5, 7)):
        count, radius_ratio, diameter_ratio = 0, 1.0, 1.0
        for cap in ("NE(0000/1000)", "SE(0000/1000)"):
            for cell in _cap_spellings(cap, resolution)[::step]:
                if _holds_a_pole(cell):
                    continue
                radius, diameter, _ = _outputs(cell)
                on_boundary = _contour(cell, 8)
                if on_boundary == (radius, diameter):
                    continue
                assert cell not in physical, cell
                count += 1
                radius_ratio = max(radius_ratio, on_boundary[0] / radius)
                diameter_ratio = max(diameter_ratio, on_boundary[1] / diameter)
        found[resolution] = (count, round(radius_ratio, 3), round(diameter_ratio, 3))
    assert found == CAP_SPELLINGS_PAST_THEIR_VERTICES
