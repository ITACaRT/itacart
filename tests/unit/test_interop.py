"""Tests for :mod:`itacart.interop`.

Every statement the module makes about its own behaviour is a named test
here, or is recorded, by name, where measurement bounds it.

Two statements hold only within a scope, and both are pinned with it. No
cell but a polar cap comes near the pole at resolution one; refined,
ordinary cells walk toward the pole without limit. A vertex at the pole
and a ring spanning the full circle are the same condition through
resolution four and not at five, where four cells hold the pole and one
spans the circle.

No tests are portable from itacart_core, which has no code for this
module.

Enumerations here walk a column range and test each address, rather than
counting up from zero until an address fails. Column zero does not exist
in the western quadrants, so the second form stops immediately there and
silently measures the eastern half of the grid while reporting a total
that looks whole. That mistake was made once while writing this file and
its cost was a table of 11 514 cells that should have read 23 014.
"""

from __future__ import annotations

import builtins
import math
import sys
from importlib.abc import MetaPathFinder
from typing import Any

import pytest
from shapely.geometry import LineString, MultiPolygon, Polygon, shape
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from shapely.wkt import loads as wkt_loads

import itacart
from itacart import interop
from itacart.boundary import plane_ring, to_geodetic
from itacart.exceptions import GeometryError, ITACaRTError, NonExistentCellError

# --------------------------------------------------------------------------
# Named cells, one per family, so a test says which family it is about
# --------------------------------------------------------------------------

INTERIOR = "NE(0500/0300)"
INTERIOR_FINE = "SW(0476/0259(4(E5(4(E2(2(B4(1(B4(2(D1(2(D3)))))))))))))"
MERIDIAN_TRIANGLE = "NE(0000/0300)"
TRAPEZOID = "NE(2003/0000)"
EQUATOR = "NE(0500/0000)"
POLAR_ROW = "NE(0001/0999)"
NORTHERN_CAP = "NE(0000/1000)"
SOUTHERN_CAP = "SE(0000/1000)"

#: The nine families: the four exception families, the four quadrants, and
#: the finest resolution the package addresses.
NINE_FAMILIES: tuple[tuple[str, str], ...] = (
    ("interior", INTERIOR),
    ("interior at resolution 13", INTERIOR_FINE),
    ("meridian column", MERIDIAN_TRIANGLE),
    ("trapezoid", TRAPEZOID),
    ("equator, row zero", EQUATOR),
    ("polar row", POLAR_ROW),
    ("polar cap", NORTHERN_CAP),
    ("south", "SE(0500/0300)"),
    ("west", "NW(0500/0300)"),
)

QUADRANTS = ("NE", "NW", "SE", "SW")

#: Wide enough to pass the last column of row zero, which is 2003.
_COLUMN_LIMIT = 2100


def _columns(quadrant: str, row: int) -> list[int]:
    """Every column of a row that names a cell, tested rather than counted."""
    return [
        column
        for column in range(_COLUMN_LIMIT)
        if itacart.is_valid_cell(f"{quadrant}({column:04d}/{row:04d})")
    ]


def _cells_of_row(row: int) -> list[str]:
    return [
        f"{quadrant}({column:04d}/{row:04d})"
        for quadrant in QUADRANTS
        for column in _columns(quadrant, row)
    ]


def _ring(cell: str) -> list[tuple[float, float]]:
    ring = itacart.cell_to_boundary(cell, close=True)
    assert isinstance(ring, list)
    return [(float(lon), float(lat)) for lon, lat in ring]


def _span(cell: str) -> float:
    longitudes = [lon for lon, _ in _ring(cell)]
    return max(longitudes) - min(longitudes)


def _apex(cell: str) -> float:
    return max(abs(lat) for _, lat in _ring(cell))


def _flatten(groups: Any) -> list[str]:
    out: list[str] = []
    for group in groups:
        out.extend(group if isinstance(group, list) else [group])
    return out


def test_the_row_enumerator_finds_the_columns_the_grid_has() -> None:
    """The enumerator is checked before anything is enumerated with it.

    Row zero has 2004 columns in the eastern quadrants and 2003 in the
    western ones, which start at column one rather than at column zero;
    the cap sits alone in row 1000 and only in the east. All four facts
    are established elsewhere, so an enumerator that disagrees with them
    is broken rather than interesting.
    """
    assert _columns("NE", 0)[0] == 0
    assert len(_columns("NE", 0)) == 2004
    assert _columns("NW", 0)[0] == 1
    assert len(_columns("NW", 0)) == 2003
    assert _columns("NE", 1000) == [0]
    assert _columns("NW", 1000) == []
    assert _columns("SW", 1000) == []


# --------------------------------------------------------------------------
# The polar cap
# --------------------------------------------------------------------------


def test_exporting_the_cap_ring_as_written_would_lose_exactly_half_of_it() -> None:
    """The cap is the one cell whose own ring is the wrong figure.

    ``cell_to_boundary`` answers three vertices: an apex at the pole and
    a base at the two ends of the seam. Read under RFC 7946, that is a
    triangle. Read on the ellipsoid, the base is a complete parallel and
    the apex is a point every longitude names, so the figure is a disc.
    A triangle and a rectangle over the same base and height stand in a
    fixed ratio, and here that ratio is the whole of the error: the
    exported polygon has exactly twice the coordinate-domain area of the
    ring it came from.

    What makes it worth a test rather than a comment is that the wrong
    answer is not malformed. The raw ring is a valid simple polygon and
    valid GeoJSON; a consumer receiving it has nothing to complain about,
    and half the cap is gone.
    """
    for cap in (NORTHERN_CAP, SOUTHERN_CAP):
        raw = Polygon(_ring(cap))
        exported = interop._cell_polygon(cap)
        assert raw.is_valid, "the wrong answer would not even look wrong"
        assert exported.is_valid
        assert exported.area == pytest.approx(2.0 * raw.area, rel=1e-12)
        assert raw.area / exported.area == pytest.approx(0.5, abs=1e-15)


def test_the_cap_exports_as_one_polygon_with_one_part() -> None:
    """A cap is produced by the seam, not divided by it.

    Cutting it at the antimeridian would invent a division the surface
    does not have, so the cap routine and the antimeridian cutter are
    deliberately different answers to deliberately different questions.
    """
    for cap in (NORTHERN_CAP, SOUTHERN_CAP):
        collection = interop.cells_to_geojson(cap)
        assert len(collection["features"]) == 1
        geometry = collection["features"][0]["geometry"]
        assert geometry["type"] == "Polygon"
        assert len(geometry["coordinates"]) == 1


def test_the_cap_is_the_only_triangle_that_is_not_equal_area() -> None:
    """``cell_shape`` does not tell a cap from a meridian triangle.

    Both answer ``triangle``. The separator is ``is_equal_area_cell``,
    and the cap is the only triangle for which it is false. Anything
    that needs to treat caps apart has to ask the second question too.
    """
    assert itacart.cell_shape(NORTHERN_CAP) == "triangle"
    assert itacart.cell_shape(MERIDIAN_TRIANGLE) == "triangle"
    assert itacart.is_equal_area_cell(NORTHERN_CAP) is False
    assert itacart.is_equal_area_cell(MERIDIAN_TRIANGLE) is True


def test_the_cap_exists_only_in_the_eastern_quadrants() -> None:
    """Row 1000 holds one cell, in NE and in SE, and nothing in the west."""
    assert itacart.is_valid_cell(NORTHERN_CAP)
    assert itacart.is_valid_cell(SOUTHERN_CAP)
    assert not itacart.is_valid_cell("NW(0000/1000)")
    assert not itacart.is_valid_cell("SW(0000/1000)")


def test_a_deep_cap_still_exports_as_a_cap() -> None:
    """The cap refines into a cap, and the routine follows it down.

    The base parallel of a resolution-five cap is a few hundred metres
    round, which puts its densification below the floor of four segments.
    The polygon still comes out closed, valid, and twice the area of the
    ring it came from.
    """
    cell = NORTHERN_CAP
    for _ in range(4):
        children = _flatten(itacart.get_children(cell))
        holders = [child for child in children if interop._holds_pole(_ring(child))]
        assert len(holders) == 1
        cell = holders[0]
    exported = interop._cell_polygon(cell)
    assert exported.is_valid
    assert exported.area == pytest.approx(2.0 * Polygon(_ring(cell)).area, rel=1e-12)
    assert len(exported.exterior.coords) >= 5


# --------------------------------------------------------------------------
# The pole detector, and the two claims about it that did not survive
# --------------------------------------------------------------------------


def test_no_cell_but_the_cap_reaches_the_pole_at_resolution_one() -> None:
    """The measurement that justified a threshold, with its scope stated.

    Every column of seven rows in all four quadrants, 23 014 cells: the
    largest latitude any of them reaches is 89.982400758..., and only the
    cap reaches ninety. The gap looks wide, and at this resolution it is.
    The claim is about resolution one and says so, which the original
    statement of it did not; what happens under refinement is the subject
    of the next test.
    """
    rows = {
        0: (8014, 0.0904369469508576),
        1: (8014, 0.18087388937666998),
        500: (5658, 45.22545419443361),
        900: (1246, 81.11818138973437),
        995: (58, 89.6242793429033),
        998: (18, 89.89287041700082),
        999: (6, 89.98240075856275),
    }
    total = 0
    highest = 0.0
    for row, (expected_n, expected_apex) in rows.items():
        cells = _cells_of_row(row)
        assert len(cells) == expected_n, f"row {row}"
        apex = max(_apex(cell) for cell in cells)
        assert apex == pytest.approx(expected_apex, rel=1e-12), f"row {row}"
        assert apex < 90.0
        highest = max(highest, apex)
        total += len(cells)
    assert total == 23014
    assert highest == pytest.approx(89.98240075856275, rel=1e-12)
    assert _apex(NORTHERN_CAP) == 90.0


def test_the_pole_detector_asks_for_the_pole_and_not_for_nearness() -> None:
    """Regression: a proximity threshold turns wedges into whole bands.

    The detector compared against 89.999 degrees. Three of the cap's
    descendants pass that at resolution four and three more at resolution
    five without being caps, and each was handed to the cap routine,
    which returns the band around the entire parallel. The three at
    resolution five came out as one identical footprint between nine and
    seventeen times their own area.

    The pole is exact everywhere else in the package:
    ``itacart.boundary.to_geodetic`` answers ``copysign(90.0, y)`` as a
    literal whenever a plane vertex reaches the meridian quadrant. Asking
    for that literal asks the grid's own question. Asking for nearness
    asks a question no fixed threshold answers, because refinement keeps
    producing cells that are nearer.
    """
    cap = NORTHERN_CAP
    witnesses = 0
    for _ in range(6):
        siblings = _flatten(itacart.get_children(cap))

        near_but_not_pole = [cell for cell in siblings if 89.999 <= _apex(cell) < 90.0]
        witnesses += len(near_but_not_pole)
        for cell in near_but_not_pole:
            assert not interop._holds_pole(_ring(cell))
            raw = Polygon(_ring(cell))
            assert interop._cell_polygon(cell).area == pytest.approx(
                raw.area, rel=1e-12
            )

        footprints = {interop._cell_polygon(cell).wkt for cell in near_but_not_pole}
        assert len(footprints) == len(near_but_not_pole), "near cells exported as one"

        holders = [cell for cell in siblings if _apex(cell) == 90.0]
        assert len(holders) == 1
        assert interop._holds_pole(_ring(holders[0]))
        cap = holders[0]

    assert witnesses > 0, "the defect's own witnesses are gone"


def test_a_pole_vertex_and_a_full_span_name_the_same_cells_where_measured() -> None:
    """An operational rule, with the population it was measured over.

    The rule is that a ring reaching the pole and a ring spanning three
    hundred and sixty degrees name the same cells. Walking the cap's own
    descent from resolution one to six, they do: exactly one child per
    level satisfies both, and the two sets are equal at every level.

    It is recorded as an empirical regularity over an enumerated
    population and not as a property of the grid, and the distinction
    earns its keep here. The *proxy* for the first condition -- a vertex
    within a thousandth of a degree of the pole -- does not agree. What
    is asserted about it is the property and not its count: the proxy
    never names the cap. How many cells it does name varies with the
    level and with the enumeration, and pinning that number pinned an
    accident: it read 0 or 3 while the descent lost children near the
    pole, and 0, 0, 0, 7, 3, 5, 3 once the enumeration was corrected. A
    rule that survives and a proxy that does not are different facts, and
    only the rule is stated.

    The descent used to stop at resolution six, where ``get_children``
    raised ``GeometryError`` on a self-intersecting refinement ring. It
    does not raise any more, so the population is claimed over seven
    refinements rather than five.
    """
    cell = NORTHERN_CAP
    levels = 0
    for _ in range(7):
        children = _flatten(itacart.get_children(cell))
        holders = [child for child in children if _apex(child) == 90.0]
        full = [child for child in children if abs(_span(child) - 360.0) < 1e-9]
        assert holders == full, f"they parted company one level below {cell}"
        assert len(holders) == 1
        near = [child for child in children if 89.999 <= _apex(child) < 90.0]
        assert not set(near) & set(holders), "the proxy must not name the cap"
        levels += 1
        cell = holders[0]
    assert levels == 7, "the rule is claimed over seven refinements and no more"


# --------------------------------------------------------------------------
# The antimeridian
# --------------------------------------------------------------------------


def test_the_widest_ring_outside_the_cap_spans_exactly_half_the_globe() -> None:
    """No cell reaches the antimeridian cutter, and the margin is exact.

    The last column of every row stops on the seam rather than passing
    it. The widest ring outside the cap spans exactly 180 degrees, in the
    polar trapezoid, so the cutter's own condition of more than 180 is
    never met by grid output.
    """
    rows = {0: 8014, 1: 8014, 500: 5658, 900: 1246, 995: 58, 998: 18, 999: 6}
    widest = 0.0
    widest_cell = None
    total = 0
    for row, expected_n in rows.items():
        cells = _cells_of_row(row)
        assert len(cells) == expected_n, f"row {row}"
        total += len(cells)
        for cell in cells:
            span = _span(cell)
            if span > widest:
                widest, widest_cell = span, cell
    assert total == 23014
    assert widest == 180.0
    assert widest_cell == "NE(0001/0999)"


def test_the_generic_cutter_is_reached_only_by_a_ring_no_cell_produces() -> None:
    """The cutter exists for geometry that did not come from the grid."""
    rectangle = Polygon([(170.0, 0.0), (-170.0, 0.0), (-170.0, 10.0), (170.0, 10.0)])
    cut = interop.geometry_to_geojson(rectangle)
    assert cut["type"] == "MultiPolygon"
    assert len(cut["coordinates"]) == 2
    ranges = sorted(
        (min(lon for lon, _ in part[0]), max(lon for lon, _ in part[0]))
        for part in cut["coordinates"]
    )
    assert ranges == [(-180.0, -170.0), (170.0, 180.0)]


def test_a_wide_ring_that_does_not_actually_cross_stays_one_polygon() -> None:
    """Width is the cutter's trigger; crossing is its subject.

    A ring whose longitudes span more than 180 degrees without occupying
    both sides of the seam lifts into one continuous piece, and comes
    back as a Polygon rather than a MultiPolygon of one part.
    """
    wide = Polygon([(-170.0, 0.0), (170.0, 0.0), (170.0, 10.0), (-170.0, 10.0)])
    longitudes = [lon for lon, _ in wide.exterior.coords]
    assert max(longitudes) - min(longitudes) > 180.0
    assert interop.geometry_to_geojson(wide)["type"] == "Polygon"


def test_geometry_that_does_not_span_the_seam_passes_through_untouched() -> None:
    narrow = Polygon([(10.0, 0.0), (11.0, 0.0), (11.0, 1.0), (10.0, 1.0)])
    result = interop.geometry_to_geojson(narrow)
    assert result["type"] == "Polygon"
    assert shape(result).equals(narrow)


def test_the_cutter_and_the_arbitrary_export_both_refuse_a_non_polygon() -> None:
    line = LineString([(0.0, 0.0), (1.0, 1.0)])
    with pytest.raises(GeometryError):
        interop.geometry_to_geojson(line)
    with pytest.raises(GeometryError):
        interop._cut_at_antimeridian(line)


# --------------------------------------------------------------------------
# Precision
# --------------------------------------------------------------------------


def test_neighbours_at_resolution_thirteen_are_lost_to_six_decimal_places() -> None:
    """Why the exporter does not round, stated as two distinct cells.

    Two adjacent cells at resolution thirteen stand 9.790384e-08 degrees
    apart in longitude. Rounded to the six places RFC 7946 suggests, both
    become the same number, and the file stays perfectly valid GeoJSON
    while two cadastral vertices become one. Section 3.1.10 says the
    digit count must not be read as an uncertainty; the index carries the
    uncertainty the coordinates cannot.
    """
    first = INTERIOR_FINE
    second = INTERIOR_FINE.replace("D3)", "D2)")
    assert first != second
    assert itacart.is_valid_cell(first) and itacart.is_valid_cell(second)
    assert itacart.get_resolution(first) == 13

    lon_first, _ = itacart.cell_to_centroid(first)
    lon_second, _ = itacart.cell_to_centroid(second)
    separation = abs(lon_first - lon_second)
    assert separation == pytest.approx(9.790384e-08, rel=1e-5)
    assert round(lon_first, 6) == round(lon_second, 6)
    assert round(lon_first, 9) != round(lon_second, 9)


def test_coordinates_leave_the_exporter_unrounded() -> None:
    """Full precision on the way out, for both exporters."""
    collection = interop.cells_to_geojson(INTERIOR_FINE)
    positions = collection["features"][0]["geometry"]["coordinates"][0]
    assert any(value != round(value, 6) for position in positions for value in position)
    text = interop.cell_to_wkt(INTERIOR_FINE)
    assert isinstance(text, str)
    tokens = text.replace(",", " ").replace("(", " ").replace(")", " ").split()
    assert any(len(token.partition(".")[2]) > 6 for token in tokens)


# --------------------------------------------------------------------------
# Refusal
# --------------------------------------------------------------------------

#: Addresses that parse but name no cell: column zero does not exist in
#: the western quadrants, and row 1000 exists only in the eastern ones.
NOT_CELLS = ("NW(0000/0300)", "SW(0000/0300)", "NW(0000/1000)", "SW(0000/1000)")

EXPORTERS = ("cells_to_geojson", "cell_to_wkt", "cells_to_wkt", "to_geodataframe")


@pytest.mark.parametrize("address", NOT_CELLS)
@pytest.mark.parametrize("exporter", EXPORTERS)
def test_every_exporter_refuses_an_index_that_names_no_cell(
    exporter: str, address: str
) -> None:
    """A refusal is what makes the defect visible.

    Six public functions answer geometry for an address that names no
    cell. An exporter that inherited that would emit a plausible polygon
    and raise nothing, which means it would pass the totality sweep while
    being wrong. Refusing is the behaviour that can be measured.
    """
    with pytest.raises(NonExistentCellError):
        getattr(interop, exporter)(address)


def test_the_refusal_is_a_package_exception_and_never_a_bare_value_error() -> None:
    """One ``except ITACaRTError`` still guards the pipeline."""
    assert issubclass(NonExistentCellError, ITACaRTError)
    assert not issubclass(NonExistentCellError, ValueError)
    with pytest.raises(ITACaRTError) as caught:
        interop.cells_to_geojson("NW(0000/0300)")
    assert isinstance(caught.value, NonExistentCellError)
    assert "names no cell" in str(caught.value)


def test_a_composed_index_is_refused_when_one_of_its_cells_is_absent() -> None:
    """The check is per cell, not per index.

    The mixed index is written out rather than built with ``compose``,
    which rewrites ``NW(0000/0300)`` into its eastern spelling and would
    hand the exporters a perfectly good pair of cells. Measured while
    writing this test; it belongs to the index module and is only noted
    here because it makes the obvious construction test nothing.
    """
    good = itacart.compose([INTERIOR, "NE(0501/0300)"])
    assert len(interop.cells_to_geojson(good)["features"]) == 2

    mixed = f"{INTERIOR},NW(0000/0300)"
    assert itacart.decompose(mixed) == [INTERIOR, "NW(0000/0300)"]
    with pytest.raises(NonExistentCellError):
        interop.cells_to_wkt(mixed)
    with pytest.raises(NonExistentCellError):
        interop.cell_to_wkt(mixed)
    with pytest.raises(NonExistentCellError):
        interop.cells_to_geojson(mixed)


# --------------------------------------------------------------------------
# Orientation and ring closure
# --------------------------------------------------------------------------


def test_the_boundary_already_turns_counterclockwise_in_all_four_quadrants() -> None:
    """Measured, and the measurement is the reason not to rely on it.

    Every quadrant's first row-zero cell has signed area +0.008124...,
    so the rings arrive right-handed already. Section 3.1.6 makes
    right-handedness a requirement on production rather than an
    observation about today's output, which is why the exporter
    normalises anyway.
    """
    areas = set()
    for quadrant in QUADRANTS:
        ring = _ring(f"{quadrant}(0001/0000)")
        signed = 0.0
        for index in range(len(ring) - 1):
            (x1, y1), (x2, y2) = ring[index], ring[index + 1]
            signed += x1 * y2 - x2 * y1
        signed /= 2.0
        assert signed > 0.0
        areas.add(round(signed, 12))
    assert areas == {round(0.008124094195761612, 12)}


def test_the_exporter_enforces_orientation_rather_than_inheriting_it() -> None:
    """A reversed ring comes back right-handed."""
    for _, cell in NINE_FAMILIES:
        assert interop._cell_polygon(cell).exterior.is_ccw, cell
    backwards = Polygon(list(reversed(_ring(INTERIOR))))
    assert not backwards.exterior.is_ccw
    assert orient(backwards, sign=1.0).exterior.is_ccw


def test_rings_close_with_four_positions_for_a_triangle_and_five_otherwise() -> None:
    """Closure is already the boundary's job, and it does it."""
    expected = {
        MERIDIAN_TRIANGLE: 4,
        NORTHERN_CAP: 4,
        INTERIOR: 5,
        TRAPEZOID: 5,
        POLAR_ROW: 5,
    }
    for cell, positions in expected.items():
        ring = _ring(cell)
        assert len(ring) == positions, cell
        assert ring[0] == ring[-1], cell


def test_the_five_families_each_export_and_are_told_apart() -> None:
    """Shape, equal-area flag and effective area, one line per family.

    The trapezoid is not equal-area, which is the whole reason the
    package separates ``nominal_cell_area`` from ``effective_cell_area``:
    row zero's trapezoid carries a quarter more area than nominal, the
    polar row's carries two thirds more, and the cap carries an eighth of
    it. Only the parallelogram and the meridian triangle come out at
    nominal exactly.
    """
    nominal = itacart.nominal_cell_area(1)
    families = {
        INTERIOR: ("parallelogram", True, 1.000000),
        MERIDIAN_TRIANGLE: ("triangle", True, 1.000000),
        TRAPEZOID: ("trapezoid", False, 1.249595),
        POLAR_ROW: ("trapezoid", False, 1.688347),
        NORTHERN_CAP: ("triangle", False, 0.121394),
    }
    for cell, (shape_name, equal_area, ratio) in families.items():
        assert itacart.cell_shape(cell) == shape_name, cell
        assert itacart.is_equal_area_cell(cell) is equal_area, cell
        assert itacart.effective_cell_area(cell) / nominal == pytest.approx(
            ratio, rel=1e-5
        ), cell
        polygon = interop._cell_polygon(cell)
        assert polygon.is_valid and polygon.area > 0.0, cell


# --------------------------------------------------------------------------
# The exported polygon is not the cell
# --------------------------------------------------------------------------

#: Points sampled along each plane edge. The departure is a sagitta and is
#: already resolved at four samples; sixty-four gives the same maximum to
#: every digit declared below, so eight is generosity rather than accuracy.
_EDGE_SAMPLES = 8


def _chord_distance(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    """Degrees from ``point`` to the straight lon/lat segment start-end."""
    (px, py), (ax, ay), (bx, by) = point, start, end
    dx, dy = bx - ax, by - ay
    if dx == 0.0 and dy == 0.0:  # pragma: no cover - no cell has a null edge
        return math.hypot(px - ax, py - ay)
    along = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    along = max(0.0, min(1.0, along))
    return math.hypot(px - (ax + along * dx), py - (ay + along * dy))


def _deviation(cell: str) -> float:
    """How far a cell's edges depart from the segments GeoJSON draws.

    A cell's edge is straight in the sinusoidal plane. RFC 7946 section
    3.1.1 makes the exported edge straight in longitude and latitude.
    Walking the first and measuring against the second is the difference
    between the two conventions, in degrees of the plane the RFC defines
    its segment in.
    """
    _, ring = plane_ring(cell)
    worst = 0.0
    count = len(ring)
    for index in range(count):
        x1, y1 = ring[index]
        x2, y2 = ring[(index + 1) % count]
        start, end = to_geodetic(x1, y1), to_geodetic(x2, y2)
        for step in range(1, _EDGE_SAMPLES):
            fraction = step / _EDGE_SAMPLES
            point = to_geodetic(x1 + fraction * (x2 - x1), y1 + fraction * (y2 - y1))
            worst = max(worst, _chord_distance(point, start, end))
    return worst


def test_the_deviation_instrument_reads_zero_on_an_edge_that_is_straight() -> None:
    """A predicate nobody has seen answer zero measures nothing.

    The lower edge of an ordinary cell runs along a parallel: latitude is
    constant, so every point of it lies on the straight lon/lat segment
    between its ends and the instrument must return exactly zero.
    """
    _, ring = plane_ring(INTERIOR)
    (x1, y1), (x2, y2) = ring[0], ring[1]
    assert y1 == y2, "the control edge is not a parallel"
    start, end = to_geodetic(x1, y1), to_geodetic(x2, y2)
    worst = max(
        _chord_distance(to_geodetic(x1 + step / 8 * (x2 - x1), y1), start, end)
        for step in range(1, 8)
    )
    assert worst == 0.0


#: The rows the departure is enumerated over, and the count of cells in
#: each. The counts are combinatorial and exact; the maxima are bracketed
#: rather than pinned because the matrix runs three operating systems.
_DEPARTURE_BY_ROW = {
    0: (8014, 5.5681626e-05),
    1: (8014, 5.5681623e-05),
    100: (7918, 5.5663028e-05),
    300: (7138, 5.7922932e-05),
    500: (5658, 5.8452497e-05),
    700: (3630, 8.0737605e-05),
    900: (1246, 2.2623814e-04),
    995: (58, 4.7658066e-03),
    999: (6, 3.7590492e-02),
    1000: (2, 1.5399336e-02),
}


@pytest.mark.slow
def test_the_exported_polygon_is_not_the_cell_by_latitude() -> None:
    """Enumerated by row, with the count of every row declared.

    41 684 cells, every column of ten rows in all four quadrants. The
    departure is never zero and never large in the middle latitudes: it
    sits near 5.6e-05 degrees from the equator to forty-five, grows by
    half again by sixty-three, and then runs away, reaching 3.8e-02 in
    the polar row.

    The unit is degrees in the longitude-latitude plane, which is the
    plane RFC 7946 defines its segment in. It is not a distance: near the
    pole a degree of longitude is short, so the figure overstates how far
    the ground moves.
    """
    total = 0
    for row, (count, worst) in _DEPARTURE_BY_ROW.items():
        cells = _cells_of_row(row)
        assert len(cells) == count, f"row {row}"
        measured = max(_deviation(cell) for cell in cells)
        assert measured > 0.0, f"row {row}: the polygon would be the cell"
        assert measured == pytest.approx(worst, rel=1e-3), f"row {row}"
        total += len(cells)
    assert total == 41684


@pytest.mark.slow
def test_the_exported_polygon_is_not_the_cell_by_shape() -> None:
    """The same enumeration bucketed by family, with each count declared.

    The ordinary cells, 41 628 of the 41 684, stay under 5e-03 degrees.
    The 56 that are not parallelograms hold the extremes, and the polar
    row's trapezoid and triangle tie for the largest departure in the
    grid.
    """
    expected = {
        "parallelogram": (41628, 4.7658066e-03),
        "polar cap": (2, 1.5399336e-02),
        "trapezoid": (36, 3.7590492e-02),
        "triangle": (18, 3.7590492e-02),
    }
    buckets: dict[str, list[float]] = {name: [] for name in expected}
    for row in _DEPARTURE_BY_ROW:
        for cell in _cells_of_row(row):
            name = str(itacart.cell_shape(cell))
            if name == "triangle" and not itacart.is_equal_area_cell(cell):
                name = "polar cap"
            buckets[name].append(_deviation(cell))
    for name, (count, worst) in expected.items():
        assert len(buckets[name]) == count, name
        assert min(buckets[name]) > 0.0, name
        assert max(buckets[name]) == pytest.approx(worst, rel=1e-3), name
    assert sum(len(values) for values in buckets.values()) == 41684


# --------------------------------------------------------------------------
# Two routes leave a cell, and only one of them is an inverse
# --------------------------------------------------------------------------

#: The row the declaration is measured over. It is the last row before the
#: cap, and it holds the largest departure in the grid.
_DECLARED_ROW = 999

#: The orthodromic step the declared coverage figures are read at. It is
#: the default of :func:`itacart.densify_orthodromic`, and the figures move
#: with it, so no coverage number in this file is written without it.
_DECLARED_STEP_M = 1000.0

#: Coverage of the reingested outline at ``_DECLARED_STEP_M``, by family,
#: with the count of cells in each. Enumerated, not sampled: six is the
#: whole of row 999 at resolution one.
_REINGESTED_COVERAGE = {
    "trapezoid": (4, 0.331111),
    "triangle": (2, 0.516779),
}

#: The same coverage read with the step refined until it stops moving.
#: The gap between this and the figures above is the point of pinning
#: both: neither number is a property of the grid.
_CONVERGED_COVERAGE = {
    "trapezoid": 0.330251,
    "triangle": 0.515698,
}


def _plane_area(cell: str) -> float:
    """Area of the cell on the plane its edges are straight in."""
    return float(Polygon(plane_ring(cell)[1]).area)


def _rebuilt_in_the_plane(cell: str) -> Polygon:
    """The exported ring read back as straight segments in the plane."""
    collection = interop.cells_to_geojson(cell)
    ring = collection["features"][0]["geometry"]["coordinates"][0]
    return Polygon([itacart.geodetic_to_sinusoidal(lon, lat) for lon, lat in ring])


def _reingested(cell: str, step_m: float = _DECLARED_STEP_M) -> Polygon:
    """The exported ring read back as a geodesic outline.

    This is not an inverse of the export. It is what a consumer does when
    it takes the exported polygon for a new geographic geometry: joins the
    positions along the ellipsoid, which RFC 7946 does not ask for and
    which draws a different figure.
    """
    outline = itacart.cell_to_polygon(cell)
    assert isinstance(outline, Polygon)
    dense = itacart.densify_orthodromic(outline, max_segment_m=step_m)
    return Polygon(
        [itacart.geodetic_to_sinusoidal(lon, lat) for lon, lat in dense.exterior.coords]
    )


def _reingested_coverage(cell: str, step_m: float = _DECLARED_STEP_M) -> float:
    """How much of the cell the reingested outline covers."""
    truth = Polygon(plane_ring(cell)[1])
    return float(_reingested(cell, step_m).intersection(truth).area) / _plane_area(cell)


def _family_of(cell: str) -> str:
    return str(itacart.cell_shape(cell))


def test_the_declared_row_holds_six_cells_in_two_families() -> None:
    """The population every figure below is measured over, enumerated.

    Four cells sit in the last column, one per quadrant, and each absorbs
    the border. Two are prime-meridian triangles and they exist in the
    eastern quadrants only, because column zero is the meridian column and
    the western quadrants do not have one. Six is the total, not a sample,
    and the arithmetic that reaches it is written out so a later reading
    cannot take four plus two for a coincidence.
    """
    cells = _cells_of_row(_DECLARED_ROW)
    assert len(cells) == 6

    by_family: dict[str, list[str]] = {}
    for cell in cells:
        by_family.setdefault(_family_of(cell), []).append(cell)
    assert {name: len(group) for name, group in by_family.items()} == {
        "trapezoid": 4,
        "triangle": 2,
    }

    assert all(itacart.absorbs_border(cell) for cell in by_family["trapezoid"])
    assert {cell[:2] for cell in by_family["trapezoid"]} == set(QUADRANTS)
    assert {cell[:2] for cell in by_family["triangle"]} == {"NE", "SE"}
    assert all(cell[3:7] == "0000" for cell in by_family["triangle"])
    assert _columns("NW", _DECLARED_ROW) == _columns("SW", _DECLARED_ROW) == [1]


def test_the_exported_ring_rebuilds_the_cell_in_the_plane() -> None:
    """The exact route, and what makes it exact.

    The exporter writes the cell's vertices and nothing between them: four
    positions for a triangle, five for a parallelogram. Joined again in the
    sinusoidal plane -- the plane the cell is defined in -- they rebuild
    the cell to the last bit, on every cell of the row.

    The ring length is asserted rather than assumed. If the exporter ever
    samples along the plane edge before inverting the projection, this test
    fails and says so, which is the only way that change becomes visible
    from here.
    """
    for cell in _cells_of_row(_DECLARED_ROW):
        collection = interop.cells_to_geojson(cell)
        ring = collection["features"][0]["geometry"]["coordinates"][0]
        assert len(ring) == (4 if _family_of(cell) == "triangle" else 5), cell
        rebuilt = _rebuilt_in_the_plane(cell)
        assert rebuilt.area == pytest.approx(_plane_area(cell), rel=1e-12), cell


def test_the_recovered_index_is_carried_rather_than_computed() -> None:
    """Identity survives the trip because it travels, not because it is read.

    ``recover_from_geojson`` takes the index from the Feature, so the round
    trip returns the same cell however the geometry was drawn. Saying the
    route preserves identity without saying where the identity rode would
    be a true claim with its scope left off.

    The control is the second half: give a Feature the geometry of a cell
    in another quadrant and the recovery is unchanged. A recovery that read
    the coordinates would answer the other cell, and this test would fail.
    """
    cells = _cells_of_row(_DECLARED_ROW)
    for cell in cells:
        assert interop.recover_from_geojson(interop.cells_to_geojson(cell)) == [cell]

    collection = interop.cells_to_geojson(cells[0])
    stranger = next(cell for cell in cells if cell[:2] != cells[0][:2])
    collection["features"][0]["geometry"] = interop.cells_to_geojson(stranger)[
        "features"
    ][0]["geometry"]
    assert interop.recover_from_geojson(collection) == [cells[0]]


def test_the_reingested_outline_covers_a_third_of_the_declared_row() -> None:
    """The figures the module docstring declares, pinned by family.

    Every number here is a reading of the instrument at
    ``_DECLARED_STEP_M``; the step is named in the constant and in the
    declaration. What does not depend on the step is the direction: the
    reingested outline never covers the cell, on any cell of this row.

    The four absorbers agree with each other exactly, and so do the two
    triangles, because each family is one figure mirrored into its
    quadrants. That agreement is asserted, so a defect that moved one
    quadrant alone could not hide inside a family average.
    """
    measured: dict[str, list[float]] = {}
    for cell in _cells_of_row(_DECLARED_ROW):
        measured.setdefault(_family_of(cell), []).append(_reingested_coverage(cell))

    for family, (count, declared) in _REINGESTED_COVERAGE.items():
        values = measured[family]
        assert len(values) == count, family
        assert max(values) - min(values) == pytest.approx(0.0, abs=1e-12), family
        assert values[0] == pytest.approx(declared, rel=1e-5), family
        assert values[0] < 1.0, family


def test_the_reingested_coverage_is_a_reading_of_its_step() -> None:
    """A number out of a numerical instrument belongs to the instrument.

    Refining the orthodromic step moves the coverage, so the declared
    figures are pinned together with the figures the refinement converges
    to, and the two are asserted to differ. Pinning only the first would
    let a reader take it for a property of the grid; pinning only the
    second would drop the number the declaration actually states.

    A coarse step reads high, because it skips the part of the geodesic
    that leaves the cell. Below one kilometre the reading settles.
    """
    for family, converged in _CONVERGED_COVERAGE.items():
        cell = next(
            cell for cell in _cells_of_row(_DECLARED_ROW) if _family_of(cell) == family
        )
        coarse = _reingested_coverage(cell, 10000.0)
        declared = _reingested_coverage(cell, _DECLARED_STEP_M)
        fine = _reingested_coverage(cell, 200.0)
        finer = _reingested_coverage(cell, 50.0)
        finest = _reingested_coverage(cell, 20.0)

        assert coarse > declared > fine > finer, family
        assert finer == pytest.approx(finest, abs=1e-4), family
        assert finer == pytest.approx(converged, rel=1e-3), family
        assert declared - converged > 5e-4, family


def test_center_does_not_return_the_declared_row_from_the_reingested_outline() -> None:
    """What ``center`` answers on this row, and why it is never the cell.

    Handed the reingested outline of one of these six, ``polyfill`` under
    ``center`` never answers the cell. For four of them it raises rather
    than answering an empty cover: the figure that comes back is narrower
    than a cell and holds no cell centre. For the two absorbers of the
    eastern quadrants it answers the cap of their own hemisphere, a cell
    the figure barely reaches. Under ``intersects`` the cell is in the
    answer, together with others the outline touches.

    The cap comes from the centre rule, and the test asserts the mechanism
    as well as the answer. The cap's centroid has ``x`` equal to zero, a
    negative zero as it happens, and every absorber of the row carries an
    edge along the meridian through it: the centroid lies exactly on the
    outline, where ``contains`` is false. The centre rule moves a centre a
    micrometre east before asking, which puts it inside the two eastern
    figures and outside the two western ones. Before that rule the six
    all raised; the outline did not change, the rule for a centre on it did.

    The control is an ordinary cell, where the same route loses less than a
    thousandth of the figure and ``center`` returns the cell itself. Without
    it the loss would read as a property of the route, and it is a property
    of this row.
    """
    from shapely.geometry import Point

    from itacart.geometry import _tie_shift

    east, north = _tie_shift(None)
    answered: dict[str, list[str]] = {}
    refused: list[str] = []
    for cell in _cells_of_row(_DECLARED_ROW):
        outline = itacart.cell_to_polygon(cell)
        assert isinstance(outline, Polygon)
        try:
            filled = itacart.polyfill(outline, 1, containment="center")
        except GeometryError as error:
            assert "narrower than one cell" in str(error), cell
            refused.append(cell)
        else:
            answered[cell] = itacart.decompose(filled)
        touched = itacart.decompose(
            itacart.polyfill(outline, 1, containment="intersects")
        )
        assert cell in touched, cell
        assert len(touched) > 1, cell

        if itacart.absorbs_border(cell):
            cap = f"{cell[0]}E(0000/1000)"
            centroid = Polygon(plane_ring(cap)[1]).centroid
            figure = _reingested(cell)
            assert centroid.x == 0.0, cell
            assert figure.boundary.distance(centroid) == 0.0, cell
            assert not figure.contains(centroid), cell
            shifted = Point(centroid.x + east, centroid.y + north)
            assert figure.contains(shifted) == (cell[1] == "E"), cell

    assert answered == {
        "NE(0001/0999)": [NORTHERN_CAP],
        "SE(0001/0999)": [SOUTHERN_CAP],
    }
    assert refused == [
        "NE(0000/0999)",
        "NW(0001/0999)",
        "SE(0000/0999)",
        "SW(0001/0999)",
    ]

    for label, cell in (("interior", INTERIOR), ("equator", EQUATOR)):
        outline = itacart.cell_to_polygon(cell)
        assert isinstance(outline, Polygon)
        resolution = itacart.get_resolution(cell)
        assert _reingested_coverage(cell) > 0.999, label
        filled = itacart.polyfill(outline, resolution, containment="center")
        assert itacart.decompose(filled) == [cell], label


# --------------------------------------------------------------------------
# WKT
# --------------------------------------------------------------------------


def test_one_cell_answers_a_string_and_several_answer_a_list() -> None:
    """The convention of the six functions that already do this."""
    single = interop.cell_to_wkt(INTERIOR)
    assert isinstance(single, str)
    assert single.startswith("POLYGON")

    index = itacart.compose([INTERIOR, "NE(0501/0300)", "NE(0502/0300)"])
    many = interop.cell_to_wkt(index)
    assert isinstance(many, list)
    assert len(many) == 3
    assert [wkt_loads(text).centroid.x for text in many] == [
        interop._cell_polygon(cell).centroid.x for cell in itacart.decompose(index)
    ]


def test_without_dissolve_the_answer_is_a_geometry_collection() -> None:
    index = itacart.compose([INTERIOR, "NE(0501/0300)"])
    text = interop.cells_to_wkt(index)
    assert text.startswith("GEOMETRYCOLLECTION (")
    collection = wkt_loads(text)
    assert collection.geom_type == "GeometryCollection"
    assert len(collection.geoms) == 2


def test_the_prime_meridian_union_closes_without_a_gap() -> None:
    """One cell claimed geometrically from both sides, so no gap opens."""
    index = itacart.compose(["NW(0001/0000)", "NE(0000/0000)", "NE(0001/0000)"])
    merged = wkt_loads(interop.cells_to_wkt(index, dissolve=True))
    assert merged.geom_type == "Polygon"
    assert merged.area == pytest.approx(0.024372277561029398, rel=1e-9)


def test_the_antimeridian_union_stays_in_two_parts() -> None:
    """Neighbours on the ellipsoid are not neighbours in longitude."""
    index = itacart.compose(["NE(2003/0000)", "NW(2003/0000)"])
    merged = wkt_loads(interop.cells_to_wkt(index, dissolve=True))
    assert merged.geom_type == "MultiPolygon"
    assert len(merged.geoms) == 2
    assert merged.area == pytest.approx(0.02030365330974189, rel=1e-9)


def test_the_interior_control_shows_the_union_can_produce_one_part() -> None:
    """A predicate nobody has seen succeed is a predicate measuring nothing.

    Without this, the two parts at the antimeridian could be an artefact
    of the instrument rather than a result. Two ordinary neighbours merge
    into one polygon, so the instrument can produce one part when one
    part is the answer.
    """
    index = itacart.compose([INTERIOR, "NE(0501/0300)"])
    merged = wkt_loads(interop.cells_to_wkt(index, dissolve=True))
    assert merged.geom_type == "Polygon"
    assert merged.area == pytest.approx(0.01821039179899931, rel=1e-9)
    polygons = [interop._cell_polygon(cell) for cell in itacart.decompose(index)]
    assert unary_union(polygons).geom_type == "Polygon"


# --------------------------------------------------------------------------
# GeoJSON structure, metadata and the id
# --------------------------------------------------------------------------


def test_the_property_names_are_the_declared_ones() -> None:
    """These are public API: renaming one breaks every consumer."""
    assert interop.INDEX_PROPERTY == "itacart_index"
    assert interop.RESOLUTION_PROPERTY == "itacart_resolution"
    assert interop.SHAPE_PROPERTY == "itacart_shape"
    assert interop.NOMINAL_AREA_PROPERTY == "itacart_nominal_area_m2"
    assert interop.EFFECTIVE_AREA_PROPERTY == "itacart_effective_area_m2"
    assert interop.EXTENSION_ZONE_PROPERTY == "itacart_extension_zone"

    properties = interop.cells_to_geojson(INTERIOR)["features"][0]["properties"]
    assert set(properties) == {
        "itacart_index",
        "itacart_resolution",
        "itacart_shape",
        "itacart_nominal_area_m2",
        "itacart_effective_area_m2",
        "itacart_extension_zone",
    }
    assert properties["itacart_index"] == INTERIOR
    assert properties["itacart_resolution"] == 1
    assert properties["itacart_shape"] == "parallelogram"
    assert properties["itacart_extension_zone"] is None


def test_the_index_goes_into_the_id_and_is_repeated_in_the_properties() -> None:
    """Section 3.2 reserves ``id``; some consumers drop it on import."""
    collection = interop.cells_to_geojson(itacart.compose([INTERIOR, TRAPEZOID]))
    assert collection["type"] == "FeatureCollection"
    assert [feature["id"] for feature in collection["features"]] == [
        INTERIOR,
        TRAPEZOID,
    ]
    for feature in collection["features"]:
        assert feature["type"] == "Feature"
        assert feature["properties"]["itacart_index"] == feature["id"]


def test_metadata_can_be_switched_off_and_extra_properties_merged() -> None:
    bare = interop.cells_to_geojson(INTERIOR, include_metadata=False)
    assert bare["features"][0]["properties"] == {}
    assert bare["features"][0]["id"] == INTERIOR

    tagged = interop.cells_to_geojson(
        INTERIOR, properties={"owner": "parcel-7"}, include_metadata=False
    )
    assert tagged["features"][0]["properties"] == {"owner": "parcel-7"}

    both = interop.cells_to_geojson(INTERIOR, properties={"owner": "parcel-7"})
    assert both["features"][0]["properties"]["owner"] == "parcel-7"
    assert both["features"][0]["properties"]["itacart_index"] == INTERIOR


def test_an_extension_cell_reports_its_zone() -> None:
    """The metadata block names which antimeridian zone a cell is in."""
    west, south, _, _ = itacart.extension_bounds("CHUKOTKA")
    footprint = Polygon(
        [
            (west + 0.01, south + 0.01),
            (west + 0.05, south + 0.01),
            (west + 0.05, south + 0.05),
            (west + 0.01, south + 0.05),
        ]
    )
    cell = itacart.decompose(itacart.polyfill(footprint, 4, compact=False))[0]
    properties = interop.cells_to_geojson(cell)["features"][0]["properties"]
    assert properties["itacart_extension_zone"] == "CHUKOTKA"


def test_one_cell_is_one_feature_even_when_the_geometry_is_a_cap() -> None:
    index = itacart.compose([NORTHERN_CAP, SOUTHERN_CAP, INTERIOR])
    collection = interop.cells_to_geojson(index)
    assert len(collection["features"]) == 3
    assert [f["id"] for f in collection["features"]] == itacart.decompose(index)


# --------------------------------------------------------------------------
# Round trip
# --------------------------------------------------------------------------


def _rounded(value: Any) -> Any:
    if isinstance(value, float):
        return round(value, 6)
    if isinstance(value, list):
        return [_rounded(item) for item in value]
    if isinstance(value, dict):
        return {key: _rounded(item) for key, item in value.items()}
    return value


def test_recovery_is_exact_and_survives_rounded_coordinates() -> None:
    """The recovery path reads the file rather than inferring from geometry.

    That is what makes it exact, and it is exact for the nine families at
    once. Rounding every coordinate to six decimal places, which
    collapses a resolution-thirteen cell into a degenerate polygon,
    leaves the recovery untouched: the index is not in the coordinates.
    """
    index = itacart.compose([cell for _, cell in NINE_FAMILIES])
    expected = itacart.decompose(index)
    collection = interop.cells_to_geojson(index)
    assert interop.recover_from_geojson(collection) == expected
    assert interop.recover_from_geojson(_rounded(collection)) == expected

    without_metadata = interop.cells_to_geojson(index, include_metadata=False)
    assert interop.recover_from_geojson(without_metadata) == expected


def test_recovery_falls_back_to_the_properties_when_the_id_is_dropped() -> None:
    collection = interop.cells_to_geojson(INTERIOR)
    del collection["features"][0]["id"]
    assert interop.recover_from_geojson(collection) == [INTERIOR]


def test_recovery_refuses_a_feature_that_carries_no_index() -> None:
    collection = interop.cells_to_geojson(INTERIOR, include_metadata=False)
    del collection["features"][0]["id"]
    with pytest.raises(GeometryError):
        interop.recover_from_geojson(collection)


@pytest.mark.parametrize(
    "obj, named",
    [
        ({"type": "Point", "coordinates": [-46.633, -23.55]}, "'Point'"),
        (
            {"type": "LineString", "coordinates": [[-46.633, -23.55], [-46.6, -23.5]]},
            "'LineString'",
        ),
        ({"type": "Feature", "id": INTERIOR, "geometry": None}, "'Feature'"),
        ({}, "None"),
    ],
)
def test_recovery_refuses_what_the_exporter_never_writes(
    obj: dict[str, object], named: str
) -> None:
    """An empty list was the answer to all four, and it read as success."""
    with pytest.raises(GeometryError, match=f"got {named}"):
        interop.recover_from_geojson(obj)


@pytest.mark.parametrize(
    "collection",
    [{"type": "FeatureCollection"}, {"type": "FeatureCollection", "features": []}],
)
def test_recovery_refuses_a_collection_without_features(
    collection: dict[str, object],
) -> None:
    """The exporter refuses an empty index, so it never writes this."""
    with pytest.raises(GeometryError, match="holds no Feature"):
        interop.recover_from_geojson(collection)
    with pytest.raises(itacart.InvalidIndexError):
        interop.cells_to_geojson("")


def test_the_geographic_round_trip_loses_the_row_below_the_pole() -> None:
    """**This measures the round trip through longitude, not the fill.**

    Read the name literally. Export a cell to GeoJSON and refill it, and
    what is being asked about is the trip: the cell leaves the sinusoidal
    plane, becomes longitude and latitude, and comes back in through the
    path meant for a caller's polygon, where every edge is rebuilt as a
    geodesic on the ellipsoid before the projection. An ITACaRT cell's
    edge is a straight line in the plane, so the figure that returns is
    not the figure that left.

    A failure here is therefore **not** evidence that ``polyfill`` is
    wrong, and reading it that way cost this project several sessions.
    The descent is exercised against the grid's own geometry in
    ``test_the_descent_recovers_the_real_tree_from_the_grids_own_geometry``,
    where all five families -- both caps included -- come back exactly.

    Measured, the trip returns 99.9269 per cent of a border-absorbing
    cell, 50 per cent of the southern cap and 33.1111 per cent of a cell
    in row 999. Eight families survive it anyway, which is what makes the
    loss look like an exception rather than the rule it is: the error is
    the same everywhere and only its size changes with latitude.

    The southern cap is one of the eight. Its cell is the northern cap's
    mirror and its round trip loses half the figure, but the half that
    survives still covers the one child there is, so the trip comes back
    exact. That is luck rather than correctness, and it is written down
    here so a later reading does not take it for the latter.

    The direction that would preserve the shape is the other one --
    sample the straight edge in the plane, then invert the projection --
    and it is an interoperability question rather than a fill one.

    The cell of row 999 comes back as a different cell, the northern cap,
    and that is the same loss showing through the centre rule. The figure
    that returns holds no centre of its own row, but it carries an edge
    along the meridian, and the cap's centroid lies exactly on it. The
    centre rule moves a centre on an outline a micrometre east before
    asking, so the eastern figure keeps the cap. Before that rule the trip
    covered nothing and the fill raised; the trip lost the cell either
    way.
    """
    exact: list[str] = []
    elsewhere: dict[str, list[str]] = {}
    refused: dict[str, str] = {}
    for label, cell in NINE_FAMILIES:
        collection = interop.cells_to_geojson(cell)
        try:
            filled = interop.from_geojson(collection, itacart.get_resolution(cell))
        except ITACaRTError as error:
            refused[label] = type(error).__name__
            continue
        cells = itacart.decompose(filled[0])
        if cells != [cell]:
            elsewhere[label] = cells
            continue
        exact.append(label)

    assert exact == [
        "interior",
        "interior at resolution 13",
        "meridian column",
        "trapezoid",
        "equator, row zero",
        "polar cap",
        "south",
        "west",
    ]
    # Not refused: reached the fill and answered a cell the trip landed on.
    # The distinction matters, because a screen turning it away would be a
    # limitation of the grid's addressing and this is a limitation of the
    # trip. Empty is the stronger pin: a refusal coming back would show.
    assert elsewhere == {"polar row": [NORTHERN_CAP]}
    assert refused == {}


def test_the_filling_loss_is_measured_rather_than_asserted() -> None:
    """Refilling at a finer resolution misses and adds, and both are small.

    At the cell's own resolution the cover is the cell and both errors
    are exactly zero, which is the control. Finer, the cover neither
    contains the exported polygon nor is contained by it: it misses about
    2.0e-07 square degrees and adds about the same, and the pair settles
    rather than shrinking. That the loss is two-sided matters, because a
    one-sided loss would let a caller argue the cover is conservative.
    """
    collection = interop.cells_to_geojson(INTERIOR)
    target = shape(collection["features"][0]["geometry"])

    measured: dict[int, tuple[int, float, float]] = {}
    for resolution in (1, 2, 3):
        filled = interop.from_geojson(collection, resolution)[0]
        cover = wkt_loads(interop.cells_to_wkt(filled, dissolve=True))
        measured[resolution] = (
            itacart.count_cells(filled),
            target.difference(cover).area,
            cover.difference(target).area,
        )

    assert measured[1] == (1, 0.0, 0.0)
    assert measured[2][0] == 4
    assert measured[3][0] == 100
    for resolution in (2, 3):
        _, missed, added = measured[resolution]
        assert 1e-07 < missed < 3e-07, resolution
        assert 1e-07 < added < 3e-07, resolution
        assert missed / target.area < 3e-05, resolution


def test_from_geojson_reads_a_collection_a_feature_and_a_bare_geometry() -> None:
    """Three shapes of input, one filling path."""
    collection = interop.cells_to_geojson(INTERIOR)
    feature = collection["features"][0]
    geometry = feature["geometry"]
    assert (
        interop.from_geojson(collection, 1)
        == interop.from_geojson(feature, 1)
        == interop.from_geojson(geometry, 1)
        == [INTERIOR]
    )


@pytest.mark.parametrize("containment", ["center", "intersects", "contains"])
def test_a_point_is_quantized_whatever_the_containment(containment: str) -> None:
    """Dimension zero has one destination: the owning cell, not a fill."""
    point = {"type": "Point", "coordinates": [-46.633, -23.55]}
    assert interop.from_geojson(point, 9, containment=containment) == [
        itacart.geo_to_cell(-46.633, -23.55, 9)
    ]


def test_a_point_on_a_cell_edge_gets_its_one_owner_not_its_neighbours() -> None:
    """``intersects`` would name every cell the point touches; the owner is one."""
    cell = itacart.geo_to_cell(-46.633, -23.55, 9)
    lon, lat = itacart.cell_to_anchor(cell)
    point = {"type": "Point", "coordinates": [lon, lat]}
    owner = itacart.geo_to_cell(lon, lat, 9)
    assert interop.from_geojson(point, 9, containment="intersects") == [owner]
    touching = itacart.polyfill(shape(point), 9, containment="intersects")
    assert itacart.count_cells(touching) > 1


def test_a_multipoint_becomes_the_region_its_owners_compose() -> None:
    coordinates = [[-46.633, -23.55], [-46.6, -23.5], [-46.633, -23.55]]
    owners = [itacart.geo_to_cell(lon, lat, 9) for lon, lat in coordinates]
    collection = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {},
                "geometry": {"type": "MultiPoint", "coordinates": coordinates},
            }
        ],
    }
    assert interop.from_geojson(collection, 9) == [itacart.compose(owners)]


@pytest.mark.parametrize("kind", ["Point", "MultiPoint"])
def test_an_empty_point_geometry_is_refused_by_name(kind: str) -> None:
    with pytest.raises(GeometryError, match=f"empty {kind} has no position"):
        interop.from_geojson({"type": kind, "coordinates": []}, 9)


@pytest.mark.parametrize(
    "obj",
    [
        {"type": "Point", "coordinates": [-46.633, -23.55]},
        {"type": "Nope", "coordinates": [0.0, 0.0]},
        {"type": "Polygon", "coordinates": "not coordinates"},
    ],
    ids=["readable point", "unsupported type", "unreadable coordinates"],
)
def test_a_bad_containment_is_refused_before_any_geometry_is_read(
    obj: dict[str, object],
) -> None:
    """The mode is judged first, so no reader error can stand in for it.

    A readable point is refused alike in either order. An unsupported type
    and coordinates that cannot be read are what tell the orders apart:
    read first, each would raise the reader's own error before the mode
    was ever looked at.
    """
    with pytest.raises(ValueError, match="containment must be"):
        interop.from_geojson(obj, 9, containment="touches")


@pytest.mark.parametrize("containment", ["center", "contains"])
def test_a_line_is_refused_for_having_no_area_not_for_being_empty(
    containment: str,
) -> None:
    line = {"type": "LineString", "coordinates": [[-46.633, -23.55], [-46.6, -23.5]]}
    with pytest.raises(GeometryError, match="LineString has no area") as caught:
        interop.from_geojson(line, 9, containment=containment)
    assert "may be empty" not in str(caught.value)
    assert interop.from_geojson(line, 9, containment="intersects")[0]


def test_from_geojson_passes_the_containment_predicate_through() -> None:
    collection = interop.cells_to_geojson(INTERIOR)
    centred = interop.from_geojson(collection, 2, containment="center")
    contained = interop.from_geojson(collection, 2, containment="contains")
    assert itacart.count_cells(centred[0]) == 4
    assert itacart.count_cells(contained[0]) <= itacart.count_cells(centred[0])


# --------------------------------------------------------------------------
# The GeoDataFrame adapters
# --------------------------------------------------------------------------
#
# These are measured against a stand-in rather than against GeoPandas, and
# the reason is a property the package holds on purpose: the suite
# has to count the same with and without the optional extras, and a skip
# keyed on an import undoes that. What the package owns here is the
# adapter -- which features it hands over, which CRS it declares, when it
# reprojects, and that a missing extra is named plainly. Those four are
# what the stand-in measures. What GeoPandas then does with the features
# is GeoPandas's contract, and asserting it here would test the dependency
# rather than the package.


class _FakeFrame:
    """The three members the adapters actually use on a frame."""

    def __init__(self, features: list[dict[str, Any]], crs: str | None) -> None:
        self.features = features
        self.crs = crs
        self.reprojected_to: list[str] = []

    def to_crs(self, crs: str) -> "_FakeFrame":
        other = _FakeFrame(self.features, crs)
        other.reprojected_to = self.reprojected_to + [crs]
        return other

    @property
    def __geo_interface__(self) -> dict[str, Any]:
        return {"type": "FeatureCollection", "features": self.features}


class _FakeGeoDataFrame:
    calls: list[tuple[int, str]] = []

    @staticmethod
    def from_features(features: list[dict[str, Any]], crs: str) -> _FakeFrame:
        _FakeGeoDataFrame.calls.append((len(features), crs))
        return _FakeFrame(list(features), crs)


class _FakeGeoPandas:
    GeoDataFrame = _FakeGeoDataFrame


@pytest.fixture()
def fake_geopandas(monkeypatch: pytest.MonkeyPatch) -> Any:
    _FakeGeoDataFrame.calls = []
    monkeypatch.setitem(sys.modules, "geopandas", _FakeGeoPandas)
    return _FakeGeoPandas


def test_to_geodataframe_hands_over_every_feature_in_wgs84(
    fake_geopandas: Any,
) -> None:
    """Metadata is always requested, and the declared CRS is EPSG:4326."""
    index = itacart.compose([INTERIOR, "NE(0501/0300)"])
    frame = interop.to_geodataframe(index)
    assert _FakeGeoDataFrame.calls == [(2, "EPSG:4326")]
    assert frame.crs == "EPSG:4326"
    assert frame.reprojected_to == []
    assert frame.features[0]["properties"]["itacart_index"] == INTERIOR


def test_to_geodataframe_reprojects_only_when_asked(fake_geopandas: Any) -> None:
    frame = interop.to_geodataframe(INTERIOR, crs="EPSG:3857")
    assert _FakeGeoDataFrame.calls == [(1, "EPSG:4326")]
    assert frame.reprojected_to == ["EPSG:3857"]
    assert frame.crs == "EPSG:3857"


def test_from_geodataframe_reprojects_before_delegating(fake_geopandas: Any) -> None:
    """The filling logic exists once, in the path that needs no extra."""
    features = interop.cells_to_geojson(INTERIOR)["features"]
    frame = _FakeFrame(features, "EPSG:3857")
    assert interop.from_geodataframe(frame, 1) == [INTERIOR]


def test_from_geodataframe_leaves_a_frame_without_a_crs_alone(
    fake_geopandas: Any,
) -> None:
    features = interop.cells_to_geojson(INTERIOR)["features"]
    frame = _FakeFrame(features, None)
    assert interop.from_geodataframe(frame, 1) == [INTERIOR]
    assert frame.reprojected_to == []


def test_a_missing_extra_is_named_plainly() -> None:
    """The error says which extra is missing and how to install it."""
    real_import = builtins.__import__

    def blocked(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "geopandas" or name.startswith("geopandas."):
            raise ImportError("blocked so the message has to be earned")
        return real_import(name, *args, **kwargs)

    saved = sys.modules.pop("geopandas", None)
    builtins.__import__ = blocked
    try:
        with pytest.raises(ImportError, match=r"itacart\[geo\]"):
            interop.to_geodataframe(INTERIOR)
        with pytest.raises(ImportError, match=r"itacart\[geo\]"):
            interop.from_geodataframe(_FakeFrame([], None), 1)
    finally:
        builtins.__import__ = real_import
        if saved is not None:
            sys.modules["geopandas"] = saved


# --------------------------------------------------------------------------
# The module's own surface
# --------------------------------------------------------------------------


def test_the_module_exports_seven_names_and_the_package_carries_them() -> None:
    assert set(interop.__all__) == {
        "cells_to_geojson",
        "cell_to_wkt",
        "cells_to_wkt",
        "to_geodataframe",
        "from_geodataframe",
        "from_geojson",
        "recover_from_geojson",
    }
    for name in interop.__all__:
        assert getattr(itacart, name) is getattr(interop, name)


def test_the_arbitrary_geometry_path_is_public_in_fact_and_not_in_name() -> None:
    """``geometry_to_geojson`` is reachable and undeclared, deliberately.

    It is not in the module's ``__all__``, so the surface invariant does
    not carry it to the package top and the census does not see it. That
    is a decision about the engine and the conformance surface, not a
    defect; what it must not be is untested, because the
    coverage threshold counts it either way.
    """
    assert "geometry_to_geojson" not in interop.__all__
    assert "geometry_to_geojson" not in itacart.__all__
    assert callable(interop.geometry_to_geojson)


def test_a_multipolygon_reaching_the_arbitrary_path_is_refused() -> None:
    """Only Polygon input is supported, and saying so is the contract."""
    parts = MultiPolygon(
        [
            Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            Polygon([(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]),
        ]
    )
    with pytest.raises(GeometryError, match="unsupported geometry type"):
        interop.geometry_to_geojson(parts)


def test_the_refusal_does_not_depend_on_whether_the_extra_is_installed() -> None:
    """Which problem gets reported must not turn on the environment.

    ``to_geodataframe`` used to import GeoPandas before it looked at the
    index, so an address naming no cell produced ``NonExistentCellError``
    where the extra was installed and ``ImportError`` where it was not.
    Both were true statements; only one was about what the caller did
    wrong. The index is checked first now, and this test makes the
    missing extra visible without needing an environment that lacks it.
    """

    class _Missing(MetaPathFinder):
        def find_spec(self, fullname: str, path: Any = None, target: Any = None) -> Any:
            if fullname == "geopandas" or fullname.startswith("geopandas."):
                raise ModuleNotFoundError(f"No module named {fullname!r}")
            return None

    saved = sys.modules.pop("geopandas", None)
    finder = _Missing()
    sys.meta_path.insert(0, finder)
    try:
        with pytest.raises(NonExistentCellError):
            interop.to_geodataframe("NW(0000/0300)")
        with pytest.raises(ImportError):
            interop.to_geodataframe(INTERIOR)
    finally:
        sys.meta_path.remove(finder)
        if saved is not None:
            sys.modules["geopandas"] = saved
