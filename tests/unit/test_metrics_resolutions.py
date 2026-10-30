"""Shape metrics at every resolution, and the controls that see them fail.

``test_metrics.py`` measures the metrics at resolution 1, where every one
of its cells lives. Below it, three things broke without a test noticing.
A tolerance of a micro-degree, about 0.11 m, was wider than every cell
from resolution 11 down, so the angle answered zero for all of them and
compactness refused nearly all of them, and from resolution 3 on it also
refused border trapezoids whose base was shorter than the tolerance. The
western half of every meridian triangle and polar cap leans the western
way under an eastern code, and reading the lean off the quadrant letter
answered zero for exactly half of the cells under each triangle. And the
geodesic solver stopped on an absolute threshold that left 0.03 degree
of azimuth in a one-centimetre line. The totality sweep called the
metrics over produced indices and accepted any answer, which is the
right instrument for its own question and the wrong one for this.

This module asks the question at every resolution, over populations
built to reach the cases that broke: chains through physical children
from resolution 1 to 13 in four quadrants, the meridian triangles, the
caps and the last column of each row; an enumeration of border cells,
because a chain leaves the border by resolution 5 or 6 and never finds a
sliver; and the descendants of the triangles and caps, enumerated to the
depth that pairs every western cell with its eastern twin. Each claim
ships with the control that watches it fail.
"""

from __future__ import annotations

import decimal
import math
from typing import Any, Iterator

import pytest

import itacart
from itacart import geodesy, metrics
from itacart.boundary import _L1, _x_border, plane_ring, ring_centroid
from itacart.constants import (
    CELL_SIZE_M,
    MERIDIAN_QUADRANT,
    WGS84_E2,
    refinement_alphabet,
)
from itacart.exceptions import ConvergenceError, GeometryError
from itacart.index import split_components
from itacart.resolutions import cell_size
from tests.unit.test_geodesy import (
    CASE_AZIMUTHS,
    CASE_LATITUDES,
    CASE_LONGITUDES,
    integrate_geodesic,
)

QUADRANTS = ("NE", "NW", "SE", "SW")
SIDE_1 = CELL_SIZE_M[1]
assert SIDE_1 is not None

#: Radius used to turn a cell side into an angle in the closed-form bound.
MEAN_RADIUS_M = 6371000.0

#: Rows the chains start from: the equator, mid latitudes, the square
#: latitude near row 900, the obtuse rows above it, and the last row below
#: the cap.
START_ROWS = (0, 1, 100, 300, 500, 700, 850, 900, 950, 990, 999)


def _starts() -> list[str]:
    """Resolution-1 cells the chains descend from, enumerated.

    Four quadrants; the meridian triangle where it exists; columns 1, 500
    and 1000; and the last two columns of every start row, which are the
    border trapezoids. The two caps close the list.
    """
    starts: list[str] = []
    for quadrant in QUADRANTS:
        for row in START_ROWS:
            top = itacart.last_lattice_column(quadrant, row, SIDE_1)
            columns = {1, 500, 1000, max(1, top - 1), top}
            if quadrant[1] == "E":
                columns.add(0)
            for column in sorted(columns):
                cell = f"{quadrant}({column:04d}/{row:04d})"
                if itacart.is_valid_cell(cell):
                    starts.append(cell)
    return starts + ["NE(0000/1000)", "SE(0000/1000)"]


def _chains() -> list[str]:
    """One chain per start, through physical children, to resolution 13.

    The child taken at each level is a fixed function of the start's
    position and the depth, so the population is the same on every run
    without a random generator, and spreads across the children instead
    of always taking the first.
    """
    cells: list[str] = []
    for position, start in enumerate(_starts()):
        cell = start
        cells.append(cell)
        for depth in range(2, 14):
            children = list(itacart.get_children(cell, flatten=True))
            cell = children[(5 * position + 3 * depth) % len(children)]
            cells.append(cell)
    return cells


def _triangle_chains() -> list[str]:
    """From each meridian triangle, down its triangle children to 13.

    A chain through arbitrary children leaves the triangles by resolution
    4, so the claim that a triangle still answers zero at resolution 13
    needs a chain of its own. The middle triangle child is taken at each
    level: two of four at an even level, five of twenty-five at an odd one.
    """
    cells: list[str] = []
    for start in _starts():
        if start[3:7] != "0000" or start.endswith("1000)"):
            continue
        cell = start
        for _ in range(2, 14):
            children = list(itacart.get_children(cell, flatten=True))
            triangles = [c for c in children if itacart.cell_shape(c) == "triangle"]
            cell = triangles[len(triangles) // 2]
            cells.append(cell)
    return cells


def _western(cell: str) -> bool:
    _, ring = plane_ring(cell)
    return ring_centroid(ring)[0] < 0.0


def _western_chains() -> list[str]:
    """From each meridian triangle and cap, down its western half to 13.

    The western half leans the western way under an eastern code, which is
    the case the quadrant letter got wrong, and below resolution 4 no
    arbitrary chain stays in it. While the chain is on a triangle it steps
    to a western child that is not one, or, under a cap at resolution 2,
    where the only child is the triangle at the pole, down that triangle.
    Every descendant of a western cell is western, so from there the chain
    continues by the same rule as the others.
    """
    cells: list[str] = []
    for position, start in enumerate(_starts()):
        if start[3:7] != "0000":
            continue
        cell = start
        for depth in range(2, 14):
            children = list(itacart.get_children(cell, flatten=True))
            if itacart.cell_shape(cell) == "triangle":
                western = [
                    c
                    for c in children
                    if _western(c) and itacart.cell_shape(c) != "triangle"
                ]
                children = western or children
            cell = children[(5 * position + 3 * depth) % len(children)]
            cells.append(cell)
    return cells


CHAINS = _chains()
TRIANGLE_CHAINS = _triangle_chains()
WESTERN_CHAINS = _western_chains()
POPULATION = CHAINS + TRIANGLE_CHAINS + WESTERN_CHAINS


def _holds_a_pole(cell: str) -> bool:
    ring = itacart.cell_to_boundary(cell)
    assert isinstance(ring, list)
    return any(abs(latitude) == 90.0 for _, latitude in ring)


def _spans_a_full_turn(cell: str) -> bool:
    """Whether the ring's longitudes reach both ends of a full turn.

    Read off the extremes of the ring, not off pairs of vertices as the
    function under test reads it: a ring that wraps the whole parallel has
    a vertex on each side of the antimeridian, a turn apart.
    """
    ring = itacart.cell_to_boundary(cell)
    assert isinstance(ring, list)
    longitudes = [longitude for longitude, _ in ring]
    return max(longitudes) - min(longitudes) >= 360.0 - 4 * math.ulp(360.0)


def _from_the_pole(cell: str) -> float:
    """Metres from the nearer pole to the closest vertex of the ring."""
    ring = itacart.cell_to_boundary(cell)
    assert isinstance(ring, list)
    closest = min(90.0 - abs(latitude) for _, latitude in ring)
    return math.radians(closest) * MEAN_RADIUS_M


def _lexical_siblings(cell: str) -> list[str]:
    """The other cells spelled under the same parent prefix."""
    components = split_components(cell)
    level = len(components) - 1
    siblings = []
    for code in refinement_alphabet(level):
        if code == components[-1]:
            continue
        codes = components[2:-1] + [code]
        spelling = (
            f"{components[0]}({components[1]}"
            + "".join(f"({c}" for c in codes)
            + ")" * (len(codes) + 1)
        )
        if itacart.is_valid_cell(spelling):
            siblings.append(spelling)
    return siblings


def _around_the_poles(resolution: int) -> dict[str, list[str]]:
    """A sample of the cells within four sides of each pole.

    Circles around the pole out to four sides, four samples per side of
    circumference, and then every cell spelled under the same parent as one
    already found, while it stays within four sides of the pole. It is a
    sample and not an enumeration: inside the caps many spellings are
    neither the answer at a sampled point nor a sibling of one, which is
    why the caps are also enumerated, to resolution 5. Near either pole at
    resolutions 11 and 12, :func:`itacart.geo_to_cell` can even return a
    spelling that names no cell, which the sample leaves out.
    """
    side = cell_size(resolution)
    fractions = [0.0005, 0.001, 0.005, 0.01, 0.02, 0.05] + [
        0.1 * k for k in range(1, 41)
    ]
    cells: dict[str, list[str]] = {}
    for name, sign in (("north", 1.0), ("south", -1.0)):
        found: set[str] = set()
        for fraction in fractions:
            radius = fraction * side
            latitude = sign * (90.0 - math.degrees(radius / MEAN_RADIUS_M))
            samples = max(32, math.ceil(4 * 2 * math.pi * radius / side))
            for k in range(samples):
                longitude = -180.0 + 360.0 * (k + 0.5) / samples
                cell = str(itacart.geo_to_cell(longitude, latitude, resolution))
                if itacart.is_valid_cell(cell):
                    found.add(cell)
        frontier = sorted(found)
        while frontier:
            fresh = [
                sibling
                for cell in frontier
                if resolution > 1
                for sibling in _lexical_siblings(cell)
                if sibling not in found and _from_the_pole(sibling) <= 4 * side
            ]
            found.update(fresh)
            frontier = sorted(set(fresh))
        cells[name] = sorted(c for c in found if _from_the_pole(c) <= 4 * side)
    return cells


def _refused(cell: str) -> bool:
    try:
        itacart.compactness(cell)
    except GeometryError:
        return True
    return False


def _closed_form(anchor: tuple[float, float]) -> float:
    """The angle of the infinitesimal cell at ``anchor``, in degrees.

    On the ellipsoidal sinusoidal plane ``x = lambda N cos(phi)``, and
    since ``d(N cos phi) / d(phi) = -M sin(phi)`` a plane displacement
    ``(dx, dy)`` moves ``dx + lambda sin(phi) dy`` east and ``dy`` north.
    The leaning side ``(-l, l)`` therefore meets the parallel at
    ``atan2(1, 1 - |lambda| sin|phi|)``, with longitude in radians from
    Greenwich: 45 degrees on the equator and on the meridian, a right
    angle where ``|lambda| sin|phi| = 1``, and obtuse beyond it.
    """
    u = abs(math.radians(anchor[0])) * math.sin(abs(math.radians(anchor[1])))
    return math.degrees(math.atan2(1.0, 1.0 - u))


def _closed_form_tolerance(anchor: tuple[float, float], side: float) -> float:
    """How far a cell of ``side`` metres may sit from the closed form, degrees.

    Two terms, both measured before they were written. The first is the
    discretisation residual: half the cell's side, as an angle, times the
    rate at which ``|lambda| sin|phi|`` changes across the plane, which is
    ``tan|phi|`` in one direction and ``|lambda| / cos(phi)`` in the other.
    Over 1 500 uniform points per resolution, the chains below and samples
    near 89 degrees and in the two extension zones, the largest measured
    residual was 0.84 of this term. The second is the rounding floor of
    coordinates written in degrees, which grows as the cell shrinks and is
    what remains at resolution 13: 4e-5 degree at most, against 6e-5 here.
    """
    phi = math.radians(abs(anchor[1]))
    lam = abs(math.radians(anchor[0]))
    residual = 0.5 * side / MEAN_RADIUS_M * (math.tan(phi) + lam / math.cos(phi))
    return math.degrees(residual + 1.1e-8 / side)


def _letter_chirality(cell: str, anchor: tuple[float, float]) -> float:
    """The rule this module exists to rule out, for the controls."""
    return 1.0 if itacart.quadrant_of(cell)[1] == "E" else -1.0


# --------------------------------------------------------------------------
# The populations, declared
# --------------------------------------------------------------------------


def _census(cells: list[str]) -> dict[int, dict[str, int]]:
    counts: dict[int, dict[str, int]] = {}
    for cell in cells:
        shapes = counts.setdefault(itacart.get_resolution(cell), {})
        shape = str(itacart.cell_shape(cell))
        shapes[shape] = shapes.get(shape, 0) + 1
    return counts


def test_the_populations_reach_the_families_they_name() -> None:
    """A green sweep over an empty bucket is still green, so count them.

    The arbitrary chains leave the triangles by resolution 4 and the
    border by resolution 6, which is why the other two chains and the
    border enumeration exist: from resolution 6 down, every one of the 192
    arbitrary chains is on a parallelogram.
    """
    assert len(_starts()) == 192
    assert {cell[:2] for cell in CHAINS} == set(QUADRANTS)
    chains = _census(CHAINS)
    assert chains[1] == {"triangle": 24, "parallelogram": 124, "trapezoid": 44}
    assert chains[5] == {"parallelogram": 191, "trapezoid": 1}
    assert all(chains[r] == {"parallelogram": 192} for r in range(6, 14))
    assert all(_census(TRIANGLE_CHAINS)[r] == {"triangle": 22} for r in range(2, 14))
    western = _census(WESTERN_CHAINS)
    assert western[2] == {"parallelogram": 22, "triangle": 2}
    assert all(western[r] == {"parallelogram": 24} for r in range(5, 14))
    assert all(
        _western(cell)
        for cell in WESTERN_CHAINS[1:]
        if itacart.get_resolution(cell) > 2
    )
    assert len(POPULATION) == 2496 + 264 + 288


# --------------------------------------------------------------------------
# Zero only on a triangle, refusal only where the ring wraps the parallel
# --------------------------------------------------------------------------


def test_zero_means_triangle_and_refusal_means_a_full_turn() -> None:
    """The two claims the metrics make about themselves, from 1 to 13.

    The angle is zero on a triangle, whose leaning vertex is the far end
    of its own base, and nowhere else. Compactness is refused for a cell
    whose ring wraps the whole parallel, and so visits the same point
    twice, and for no other. Over these chains the cells that wrap are the
    ones holding a pole. Both claims used to hold at resolution 1 only.
    """
    zero_off_triangle = [
        cell
        for cell in POPULATION
        if (itacart.cell_base_angle(cell) == 0.0)
        != (itacart.cell_shape(cell) == "triangle")
    ]
    assert zero_off_triangle == []
    refused_off_turn = [
        cell for cell in POPULATION if _refused(cell) != _spans_a_full_turn(cell)
    ]
    assert refused_off_turn == []
    assert all(_holds_a_pole(cell) for cell in POPULATION if _spans_a_full_turn(cell))


def _cap_subtree(cap: str, resolution: int) -> list[str]:
    """Every cell spelled under a cap at one resolution, its whole subtree."""
    paths: list[list[str]] = [[]]
    for level in range(2, resolution + 1):
        paths = [path + [code] for path in paths for code in refinement_alphabet(level)]
    spellings = [
        cap[:-1] + "".join(f"({code}" for code in path) + ")" * (len(path) + 1)
        for path in paths
    ]
    return [spelling for spelling in spellings if itacart.is_valid_cell(spelling)]


#: Per cap and resolution: the cells of its lexical subtree, and how many of
#: them have a ring that wraps the parallel, the triangle at the pole
#: included. Enumerated, not sampled.
CAP_CENSUS = {2: (1, 1), 3: (16, 1), 4: (98, 17), 5: (3360, 731)}


def test_compactness_refuses_exactly_the_rings_that_wrap_the_parallel() -> None:
    """Inside the caps, every cell to resolution 5, and a sample below it.

    The cell that holds the pole wraps at every resolution: the parallel
    through it collapses onto one triangle, and its ring runs pole, -180,
    +180. From resolution 4 down many trapezoids around it wrap as well,
    most of them several times over, with rings that reach some nineteen
    turns of longitude at resolution 5. Compactness refuses exactly the
    cells that wrap and measures every other one, over the whole subtree
    of each cap to resolution 5 and over the sample of cells around the
    poles at every resolution.
    """
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        for resolution, counts in CAP_CENSUS.items():
            cells = _cap_subtree(cap, resolution)
            wrapping = [cell for cell in cells if _spans_a_full_turn(cell)]
            assert (len(cells), len(wrapping)) == counts, (cap, resolution)
            assert [cell for cell in cells if _refused(cell)] == wrapping
            assert sum(_holds_a_pole(cell) for cell in wrapping) == 1
    for resolution in range(1, 14):
        for cells in _around_the_poles(resolution).values():
            wrapping = [cell for cell in cells if _spans_a_full_turn(cell)]
            assert [cell for cell in cells if _refused(cell)] == wrapping
            assert sum(_holds_a_pole(cell) for cell in wrapping) == 1


def test_a_ring_that_wraps_more_than_once_is_refused_by_its_span(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control: without the span, such a ring gets a number.

    ``NE(0000/1000(1(B3(4(D5)))))`` comes back with longitudes spanning
    691.7 degrees and no two vertices a full turn apart, so the repeated
    vertex never fires on it. Release 1.0.0 measured it at 0.0096, which is
    the quotient of a ring crossing itself and not of the cell.
    """
    cell = "NE(0000/1000(1(B3(4(D5)))))"
    ring = itacart.cell_to_boundary(cell)
    assert isinstance(ring, list)
    assert metrics._repeated_vertex(ring) is None
    assert round(metrics._longitude_span(ring), 1) == 691.7
    with pytest.raises(GeometryError, match="wraps the whole parallel"):
        itacart.compactness(cell)
    monkeypatch.setattr(metrics, "_longitude_span", lambda ring: 0.0)
    assert round(itacart.compactness(cell), 4) == 0.0096


def test_the_cell_holding_a_pole_is_refused_at_every_resolution() -> None:
    """One such cell per pole per resolution.

    Below the cap the parallel through the pole still collapses onto one
    triangle per resolution, and its ring keeps the cap's defect. At some
    resolutions its two antimeridian vertices come out at
    180.00000000000003 rather than 180, so the full turn is read to the
    rounding of a longitude, and this is what pins that.
    """
    for latitude in (90.0, -90.0):
        for resolution in range(1, 14):
            cell = itacart.geo_to_cell(0.0, latitude, resolution)
            assert _holds_a_pole(cell)
            assert _refused(cell)
            assert itacart.cell_base_angle(cell) == 0.0
    for cap in ("NE(0000/1000)", "SE(0000/1000)"):
        under = list(itacart.get_descendants(cap, 4))
        assert len(under) == 44
        assert [cell for cell in under if _refused(cell)] == [
            cell for cell in under if _holds_a_pole(cell)
        ]
        assert len([cell for cell in under if _refused(cell)]) == 1


def test_the_zero_instrument_catches_the_rule_it_replaced(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control: put the quadrant letter back and the sweep reports.

    Reading the lean off the letter is what answered zero under every
    meridian triangle. Restoring it here turns the western half of the
    triangle's children into zeros, which is the shape of the defect.
    """
    west, east = "NE(0000/0900(3))", "NE(0000/0900(2))"
    assert itacart.cell_base_angle(west) == pytest.approx(45.212720312, abs=1e-9)
    monkeypatch.setattr(metrics, "_chirality", _letter_chirality)
    assert itacart.cell_base_angle(west) == 0.0
    assert itacart.cell_base_angle(east) != 0.0
    assert itacart.cell_shape(west) == "parallelogram"


# --------------------------------------------------------------------------
# The mirrored twin
# --------------------------------------------------------------------------


def _triangle_family() -> Iterator[str]:
    """Every descendant of six meridian triangles and the two caps.

    Both to resolution 4. A western cell's twin is found on the plane
    rather than in this list, so the depth is a cost and not a condition.
    """
    roots = (
        "NE(0000/0000)",
        "NE(0000/0300)",
        "NE(0000/0900)",
        "NE(0000/0999)",
        "SE(0000/0300)",
        "SE(0000/0650)",
        "NE(0000/1000)",
        "SE(0000/1000)",
    )
    for root in roots:
        for resolution in (2, 3, 4):
            yield from itacart.get_descendants(root, resolution)


def test_the_mirrored_twin_leans_by_the_same_angle() -> None:
    """A western cell and its eastern mirror measure the same angle.

    The twin is found on the plane, at the mirror image of the cell's own
    centroid, and its ring is asserted to be the mirror image bit for bit
    before the angles are compared. The angles then agree to a few units
    in the last place: they are the same lean measured from opposite
    sides, through a geodesic solver whose last bits depend on the order
    of operations. Exact equality would be the wrong claim, and it fails
    on more than a third of the pairs.
    """
    family = list(_triangle_family())
    assert len(family) == 3134
    pairs = 0
    exact = 0
    worst_ulps = 0.0
    for cell in family + WESTERN_CHAINS:
        if itacart.cell_shape(cell) == "triangle":
            continue
        angle = itacart.cell_base_angle(cell)
        assert isinstance(angle, float)
        assert angle != 0.0
        _, ring = plane_ring(cell)
        x, y = ring_centroid(ring)
        if x >= 0.0:
            continue
        twin = itacart.sinusoidal_to_cell(-x, y, itacart.get_resolution(cell))
        _, twin_ring = plane_ring(twin)
        assert sorted((-u, v) for u, v in ring) == sorted(twin_ring)
        twin_angle = itacart.cell_base_angle(twin)
        assert isinstance(twin_angle, float)
        pairs += 1
        exact += angle == twin_angle
        worst_ulps = max(worst_ulps, abs(angle - twin_angle) / math.ulp(angle))
    # Half of the enumerated family lies west of the meridian, and every
    # cell of the western chains does: 1 464 and 286.
    assert pairs == 1464 + 286
    assert 0 < exact < pairs
    assert worst_ulps <= 16


# --------------------------------------------------------------------------
# Convergence to the closed form
# --------------------------------------------------------------------------


def test_the_angle_converges_to_its_closed_form_below_89_degrees() -> None:
    """Every chain cell below 89 degrees sits within the measured bound.

    The closed form is the angle of the infinitesimal cell at the anchor,
    so it is an oracle that needs no library. The bound shrinks with the
    side until the rounding of the coordinates takes over at resolution
    13. Above 89 degrees the finite size of the cell grows towards the
    pole faster than this bound allows, and the oracle would have to be
    another one there.
    """
    measured = 0
    worst = 0.0
    for cell in POPULATION:
        if itacart.cell_shape(cell) == "triangle":
            continue
        anchor = itacart.cell_to_anchor(cell)
        assert isinstance(anchor, tuple)
        if abs(anchor[1]) >= 89.0:
            continue
        side = cell_size(itacart.get_resolution(cell))
        angle = itacart.cell_base_angle(cell)
        assert isinstance(angle, float)
        departure = abs(angle - _closed_form(anchor))
        worst = max(worst, departure / _closed_form_tolerance(anchor, side))
        measured += 1
    assert measured == 2396
    assert worst <= 1.0


def test_the_closed_form_bound_is_not_vacuous(monkeypatch: pytest.MonkeyPatch) -> None:
    """The bound refuses the zero the letter rule produced, at resolution 9."""
    cell = WESTERN_CHAINS[7]
    assert itacart.get_resolution(cell) == 9
    assert cell.startswith("NE(0000/0000(")
    anchor = itacart.cell_to_anchor(cell)
    assert isinstance(anchor, tuple)
    tolerance = _closed_form_tolerance(anchor, cell_size(9))
    angle = itacart.cell_base_angle(cell)
    assert isinstance(angle, float)
    assert abs(angle - _closed_form(anchor)) <= tolerance
    monkeypatch.setattr(metrics, "_chirality", _letter_chirality)
    assert abs(itacart.cell_base_angle(cell) - _closed_form(anchor)) > tolerance


# --------------------------------------------------------------------------
# The border slivers
# --------------------------------------------------------------------------


def _border_cells(resolution: int, per_quadrant: int) -> list[str]:
    """The last cell of a row at the antimeridian border, on an ordinate grid.

    Enumerated, because a chain leaves the border within a few levels and
    never meets the slivers the border cuts at depth. The ordinates are
    evenly spaced from the equator to the row below the cap, in all four
    quadrants, and the cell taken is the one holding the point a hair
    inside the border at that ordinate.
    """
    side = cell_size(resolution)
    cells: list[str] = []
    for quadrant in QUADRANTS:
        x_sign = 1.0 if quadrant[1] == "E" else -1.0
        y_sign = 1.0 if quadrant[0] == "N" else -1.0
        for step in range(per_quadrant):
            y = (step + 0.5) / per_quadrant * (MERIDIAN_QUADRANT - 2.0 * _L1)
            border = _x_border(quadrant, int(y // _L1), y)
            cell = itacart.sinusoidal_to_cell(
                x_sign * (border - 1e-6 * side), y_sign * y, resolution
            )
            if cell not in cells:
                cells.append(cell)
    return cells


#: Slivers named because they are the thinnest found. The first two are a
#: base of 6.1 cm at resolution 3, north and south; the third a base of
#: 2.9 micrometres at resolution 13, found by enumerating 16 000 border
#: cells there.
NAMED_SLIVERS = (
    "NE(1984/0090(3(E3)))",
    "SE(1984/0090(3(E3)))",
    "SE(0854/0719(2(D2(1(C4(3(A1(2(E4(1(A3(4(A3)))))))))))))",
)


def test_the_border_slivers_are_measured_rather_than_refused() -> None:
    """A trapezoid has an area however short its base, so it has a shape.

    A sliver has no minimum thickness: the border cuts the lattice
    wherever it falls, and the thinner the cell the closer its base comes
    to a single point. A tolerance on distance only moves the problem to
    a finer resolution, so the predicate is the construction instead.
    """
    cells = list(NAMED_SLIVERS)
    for resolution in (3, 5, 7, 9, 11, 13):
        cells += _border_cells(resolution, 40)
    assert len(cells) == 3 + 6 * 160
    for cell in cells:
        assert itacart.cell_shape(cell) == "trapezoid", cell
        assert not _refused(cell), cell
        assert itacart.cell_base_angle(cell) != 0.0, cell
    _, ring = plane_ring(NAMED_SLIVERS[2])
    bases = [
        abs(ring[k][0] - ring[k - 1][0])
        for k in range(len(ring))
        if ring[k][1] == ring[k - 1][1]
    ]
    assert round(min(bases) * 1e6, 1) == 2.9


def test_the_old_predicate_refused_a_sliver_this_one_measures() -> None:
    """The control: the tolerance in degrees refuses the 6.1 cm base."""
    ring = itacart.cell_to_boundary(NAMED_SLIVERS[0])
    assert isinstance(ring, list)

    def by_tolerance(points: list[tuple[float, float]]) -> tuple[int, int] | None:
        for i, (lon_i, lat_i) in enumerate(points):
            for j in range(i + 1, len(points)):
                lon_j, lat_j = points[j]
                turned = ((lon_i - lon_j + 180.0) % 360.0) - 180.0
                if abs(lat_i - lat_j) <= 1e-6 and abs(turned) <= 1e-6:
                    return i, j
        return None

    assert by_tolerance(ring) is not None
    assert metrics._repeated_vertex(ring) is None
    assert itacart.compactness(NAMED_SLIVERS[0]) == pytest.approx(0.4774, abs=1e-4)


# --------------------------------------------------------------------------
# The azimuth of a short line
# --------------------------------------------------------------------------


def _azimuth_departures(**options: Any) -> dict[float, float]:
    """Worst azimuth error per length, against the integrated geodesic.

    The integrator solves the equations that define a geodesic, from a
    start point along a known azimuth, and the inverse solution is asked
    for the azimuth from the start to where the integrator arrived. Four
    steps are plenty at these lengths and keep the integrator's own
    rounding below what is being measured.
    """
    worst: dict[float, float] = {}
    for length in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
        worst[length] = 0.0
        for latitude in CASE_LATITUDES:
            for longitude in CASE_LONGITUDES:
                for azimuth in CASE_AZIMUTHS:
                    end = integrate_geodesic(longitude, latitude, azimuth, length, 4)
                    assert end is not None
                    _, solved = itacart.inverse_geodesic(
                        longitude, latitude, end[0], end[1], **options
                    )
                    error = abs((solved - azimuth + 540.0) % 360.0 - 180.0)
                    worst[length] = max(worst[length], error)
    return worst


def _azimuth_bound(length: float) -> float:
    """A few times the rounding of a coordinate, over the line's length."""
    return 4e-7 / length


def test_a_short_line_keeps_its_azimuth() -> None:
    """From one centimetre to a kilometre, at the rounding of the inputs.

    The departure an absolute stopping rule left behind was 3e-3 degree
    at one centimetre on this grid and 1.1e-4 at one metre; the bound
    here is 4e-5 and 4e-7.
    """
    for length, departure in _azimuth_departures().items():
        assert departure <= _azimuth_bound(length), (length, departure)


def test_the_azimuth_bound_sees_an_unconverged_solution() -> None:
    """The control: stop the solver early and the bound reports it."""
    loose = _azimuth_departures(tol=1e-3)
    assert loose[0.01] > _azimuth_bound(0.01)


def _near_antipodal_solutions() -> list[tuple[float, float, int] | None]:
    """The lines near the antipodes, each solved or refused.

    6 724 lines, from four starting points to within 0.6 degree of each
    one's antipode in steps of 0.03 degree, solved by the instrumented
    solver at its default tolerance and cap; ``None`` where it refuses.
    """
    tol, cap = geodesy.VINCENTY_TOL_RAD, geodesy.VINCENTY_MAX_ITER
    solutions: list[tuple[float, float, int] | None] = []
    for latitude in (0.0, 0.3, 10.0, 45.0):
        for i in range(-20, 21):
            for j in range(-20, 21):
                try:
                    solution = geodesy._vincenty_inverse(
                        0.0,
                        latitude,
                        180.0 + i * 0.03,
                        -latitude + j * 0.03,
                        tol,
                        cap,
                    )
                except ConvergenceError:
                    solutions.append(None)
                else:
                    solutions.append(solution)
    return solutions


#: Two lines the sine-scaled rule refused, and their lengths in metres.
REFUSED_BY_THE_SINE_RULE = (
    ((0.0, 0.0, 180.0, 0.02), 20001719.97310625),
    ((0.0, 0.0, 179.64, 0.6), 19931638.117068194),
)

#: A line of a centimetre and a half, where the two stopping rules part.
SHORT_LINE = (0.0, 30.0, 1e-7, 30.0000001)


def test_a_long_line_converges_as_it_always_did(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Past one radian the stopping rule is the absolute one, unchanged.

    A threshold scaled by the sine of the arc shrinks again towards the
    antipode, where the iteration is slowest, and asks it for digits it
    cannot reach within its cap: such a rule refused lines that the
    absolute one had answered to a tenth of a millimetre. Capped at one
    radian, the relative rule leaves every line longer than that exactly
    as release 1.0.0 solved it.

    The solver differs from release 1.0.0 only in its stopping rule, and
    taking ``min(1, sigma)`` as one restores the old rule. So the lines
    near the antipodes are solved twice on the same machine, with this
    release's rule and with the old one, and must agree to the last bit,
    iterations included, with every refusal shared. A digest pinned on one
    machine would not travel: the last bits of a near-antipodal solution
    follow the platform's trigonometric functions, and some of these
    lines converge at the last iteration the cap allows. The short line
    shows that the restoration takes hold. The two named lines, which the
    sine-scaled rule refused, are answered alike by both rules, at lengths
    pinned to a tenth of a millimetre.
    """
    tol, cap = geodesy.VINCENTY_TOL_RAD, geodesy.VINCENTY_MAX_ITER
    this_release = _near_antipodal_solutions()
    named = [
        geodesy._vincenty_inverse(*line, tol, cap)
        for line, _ in REFUSED_BY_THE_SINE_RULE
    ]
    short = geodesy._vincenty_inverse(*SHORT_LINE, tol, cap)
    monkeypatch.setattr(geodesy, "min", lambda *_: 1.0, raising=False)
    assert geodesy._vincenty_inverse(*SHORT_LINE, tol, cap) != short
    assert _near_antipodal_solutions() == this_release
    assert [
        geodesy._vincenty_inverse(*line, tol, cap)
        for line, _ in REFUSED_BY_THE_SINE_RULE
    ] == named
    for (_, length), (distance, _, _) in zip(REFUSED_BY_THE_SINE_RULE, named):
        assert distance == pytest.approx(length, abs=1e-4)
    answered = [solution for solution in this_release if solution is not None]
    assert len(answered) < len(this_release)
    assert max(iterations for _, _, iterations in answered) > 150


# --------------------------------------------------------------------------
# The authalic latitude near the pole
# --------------------------------------------------------------------------

_PI = decimal.Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459"
)


def _decimal_cos_beta(latitude: float) -> decimal.Decimal:
    """``cos`` of the authalic latitude at 60 digits, for latitudes near 90.

    Every operation is decimal, so the arcsine's ill-conditioning costs
    digits that 60 of them can spare. The sine of the latitude is the
    cosine of the colatitude, by its series; the colatitude is formed from
    the very float the function under test receives.
    """
    with decimal.localcontext() as context:
        context.prec = 60
        colatitude = (decimal.Decimal(90) - decimal.Decimal(latitude)) * _PI / 180
        sine = decimal.Decimal(0)
        term = decimal.Decimal(1)
        for k in range(1, 30):
            sine += term
            term *= -colatitude * colatitude / ((2 * k - 1) * (2 * k))
        e2 = decimal.Decimal(WGS84_E2)
        e = e2.sqrt()

        def q(s: decimal.Decimal) -> decimal.Decimal:
            t = e * s
            atanh = ((1 + t) / (1 - t)).ln() / 2
            return (1 - e2) * (s / (1 - t * t) + atanh / e)

        ratio = q(sine) / q(decimal.Decimal(1))
        return (1 - ratio * ratio).sqrt()


def test_the_authalic_latitude_keeps_its_digits_near_the_pole() -> None:
    """The cosine of the authalic latitude, to fourteen digits, to the pole.

    Compactness is measured in a Lambert projection centred on the cell,
    and that projection runs on the authalic latitude. Taken as the
    arcsine of ``q / q_p``, its cosine lost a digit for every decade the
    latitude moved towards the pole: at a micro-degree of colatitude it
    kept seven. Formed from the area between the latitude and the pole, it
    keeps them all. The control is the arcsine form itself, measured
    against the same reference.
    """
    for colatitude in (30.0, 1.0, 1e-2, 1e-4, 1e-6, 1e-7):
        latitude = 90.0 - colatitude
        reference = _decimal_cos_beta(latitude)
        _, cos_beta = metrics._authalic_sin_cos(latitude)
        assert abs(decimal.Decimal(cos_beta) / reference - 1) < decimal.Decimal("1e-14")

        sin_phi, _ = metrics._sin_cos_deg(latitude)
        by_arcsine = math.cos(math.asin(metrics._authalic_q(sin_phi) / metrics._QP))
        error = abs(decimal.Decimal(by_arcsine) / reference - 1)
        if colatitude <= 1e-4:
            assert error > decimal.Decimal("1e-10")


def test_the_authalic_latitude_is_odd_and_bounded() -> None:
    """South mirrors north to the bit, and the pair stays on the unit circle.

    The sine and the cosine are computed by different routes, so their
    squares sum to one only to a few units in the last place: 1.3e-15 at
    worst over these latitudes, against eight units here.
    """
    for latitude in (0.0, 12.5, 45.0, 45.000001, 60.0, 89.99, 90.0):
        north = metrics._authalic_sin_cos(latitude)
        south = metrics._authalic_sin_cos(-latitude)
        assert north[0] == -south[0] and north[1] == south[1]
        assert abs(north[0] ** 2 + north[1] ** 2 - 1.0) <= 8 * math.ulp(1.0)
    assert metrics._authalic_sin_cos(90.0) == (1.0, 0.0)
