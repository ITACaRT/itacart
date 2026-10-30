"""Shape metrics: the instrument first, then the grid.

Every number in this file was measured before it was written down. The
order of the file is the order of the argument: the compactness formula
is checked against shapes whose quotient is known in closed form, then
the projection it runs in is checked against the package's own areas,
and only then is a cell measured. A metric that disagrees with the grid
after all three have passed is evidence about the grid; one that has not
passed them is evidence about nothing.

Two claims here are negative -- that no cell scores above 1, and that
only the polar caps are refused -- and each ships with the control that
sees it fail, because a predicate nobody has watched fail is a predicate
nobody knows measures anything.
"""

from __future__ import annotations

import math
import statistics

import pytest

import itacart
from itacart.constants import CELL_SIZE_M
from itacart.exceptions import (
    GeometryError,
    ITACaRTError,
    NonExistentCellError,
)
from itacart.metrics import (
    _laea,
    _planar_area,
    _planar_perimeter,
    _projected_ring,
    _repeated_vertex,
)

QUADRANTS = ("NE", "NW", "SE", "SW")
#: The quadrants that carry a meridian column. Column zero is the
#: prime-meridian column and it belongs to the east: west of the
#: meridian the count starts at column one, so a western zero column
#: is a spelling that names no cell.
EASTERN = ("NE", "SE")
SIDE_1 = CELL_SIZE_M[1]
assert SIDE_1 is not None
ROWS = 1001

#: The window the enumerated claims below run over: the meridian column
#: and the last four columns of every row, in all four quadrants. Named
#: because a bound is only as wide as the set it was measured on.
BORDER_WINDOW = "column 0 and the last four columns of each row"


def _quotient(area: float, perimeter: float) -> float:
    return 4.0 * math.pi * area / (perimeter * perimeter)


def _border_window(quadrant: str, row: int) -> list[int]:
    top = itacart.last_lattice_column(quadrant, row, SIDE_1)
    return [0] + [c for c in range(max(1, top - 3), top + 1)]


def _enumerate_window() -> list[str]:
    cells = []
    for quadrant in QUADRANTS:
        for row in range(ROWS):
            for column in _border_window(quadrant, row):
                cell = f"{quadrant}({column:04d}/{row:04d})"
                try:
                    itacart.cell_shape(cell)
                except ITACaRTError:
                    continue
                cells.append(cell)
    return cells


# --------------------------------------------------------------------------
# The instrument, before any cell
# --------------------------------------------------------------------------


def test_the_quotient_scores_the_shapes_the_paper_puts_on_its_ladder() -> None:
    """Figure 11(b) of Kmoch et al. (2022) read as three closed forms.

    The paper places the hexagonal grids near 0.9, the square-based ones
    near 0.77 to 0.79 and the triangular ones near 0.6. Those are the
    quotients of the regular hexagon, the square and the equilateral
    triangle, and if this formula does not return them the disagreement
    that follows is the instrument's and not the grid's.
    """
    side = 1.0
    hexagon = _quotient(3 * math.sqrt(3) / 2 * side * side, 6 * side)
    square = _quotient(side * side, 4 * side)
    triangle = _quotient(math.sqrt(3) / 4 * side * side, 3 * side)
    assert round(hexagon, 4) == 0.9069
    assert round(square, 4) == 0.7854
    assert round(triangle, 4) == 0.6046
    assert hexagon > square > triangle


def test_the_quotient_falls_when_the_shape_does() -> None:
    """The control that sees the ladder fail.

    A quotient that returned 0.9 for everything would pass the test
    above on the hexagon alone. Squeezing a square into a sliver has to
    move the number, and towards zero.
    """
    fat = _quotient(1.0, 4.0)
    thin = _quotient(0.01, 2 * (1.0 + 0.01))
    assert thin < fat
    assert thin < 0.05


def test_the_planar_helpers_measure_a_square_of_known_size() -> None:
    ring = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    assert _planar_area(ring) == pytest.approx(100.0)
    assert _planar_perimeter(ring) == pytest.approx(40.0)


def test_the_projection_puts_its_own_origin_at_the_origin() -> None:
    assert _laea(12.0, 42.0, 12.0, 42.0) == pytest.approx((0.0, 0.0), abs=1e-9)


def test_the_projection_is_equal_area_where_the_cell_is_narrow() -> None:
    """The control that says the LAEA transcription is that projection.

    An equal-area projection has to return the package's own area for
    the same ring. It does, to better than one part in a hundred, over
    the 12 000 cells of the window below rows 900 -- and it stops doing
    so above them, which is a property of drawing a wide cell's sides as
    straight lines, not of the projection. The exception is named with
    its number rather than left for someone to find.
    """
    checked = 0
    worst = 0.0
    for cell in _enumerate_window():
        row = int(cell[8:12])
        if row >= 900:
            continue
        try:
            projected = _projected_ring(cell)
        except GeometryError:
            continue
        measured = _planar_area(projected)
        stated = itacart.effective_cell_area(cell)
        assert isinstance(stated, float)
        worst = max(worst, abs(measured - stated) / stated)
        checked += 1
    assert checked == 16200
    assert round(worst, 6) == 0.007467

    widest = "NE(0001/0999)"
    ring = itacart.cell_to_boundary(widest)
    assert isinstance(ring, list)
    longitudes = [point[0] for point in ring]
    assert max(longitudes) - min(longitudes) == pytest.approx(180.0)
    stated = itacart.effective_cell_area(widest)
    assert isinstance(stated, float)
    departure = abs(_planar_area(_projected_ring(widest)) - stated) / stated
    assert round(departure, 4) == 0.6338


# --------------------------------------------------------------------------
# Compactness of cells
# --------------------------------------------------------------------------


def test_compactness_is_pinned_on_the_four_quadrants_at_the_equator() -> None:
    for quadrant in QUADRANTS:
        cell = f"{quadrant}(0500/0000)"
        assert round(itacart.compactness(cell), 6) == 0.539207


def test_the_square_latitude_is_pinned_by_number() -> None:
    """Row 900 is where the shear cancels the meridian convergence.

    The cell is a square to five decimal places, and the two numbers are
    written out rather than compared with a tolerance so that correcting
    either of them also has to fail here.
    """
    perfect = math.pi / 4.0
    assert round(perfect, 6) == 0.785398
    for quadrant in QUADRANTS:
        cell = f"{quadrant}(0100/0900)"
        assert round(itacart.compactness(cell), 6) == 0.785388
    assert abs(0.785388 - perfect) < 2e-5


def test_compactness_is_not_monotone_in_latitude() -> None:
    """It rises to the square and comes back down.

    Sweeping too few rows shows half the curve and reports a monotone
    rise that is not there.
    """
    ladder = {
        "NE(0500/0000)": 0.539207,
        "NE(0100/0300)": 0.564744,
        "NE(0100/0600)": 0.609529,
        "NE(0100/0900)": 0.785388,
        "NE(0100/0950)": 0.536390,
    }
    for cell, expected in ladder.items():
        assert round(itacart.compactness(cell), 6) == expected
    values = list(ladder.values())
    assert values[3] == max(values)
    assert values[4] < values[0]


def test_compactness_of_the_border_shapes_is_pinned() -> None:
    assert round(itacart.compactness("NE(0000/0300)"), 6) == 0.538923
    assert round(itacart.compactness("NE(2003/0000)"), 6) == 0.650915


def test_compactness_is_bounded_over_the_enumerated_window() -> None:
    """``0 < C <= 1`` for every cell the metric accepts.

    Enumerated over the window, not sampled from it: the extremes of the
    quotient sit on the border families, which a random draw is most
    likely to miss.
    """
    values = []
    refused = []
    for cell in _enumerate_window():
        try:
            value = itacart.compactness(cell)
        except GeometryError:
            refused.append(cell)
            continue
        assert isinstance(value, float)
        assert 0.0 < value <= 1.0, f"{cell} scored {value}"
        values.append(value)
    assert len(values) == 17988
    assert round(min(values), 6) == 0.094971
    assert round(max(values), 6) == 0.785398
    assert refused == [f"{quadrant}(0000/1000)" for quadrant in EASTERN]


def test_the_bound_is_not_vacuous() -> None:
    """The control that sees ``C <= 1`` fail.

    Feeding the quotient the polar cap's stated area against the length
    of its own ring returns 9.8696, which is impossible for a simple
    closed curve and is what the refusal is protecting the bound from.
    The number is here so that the bound above is known to be a claim
    about the grid rather than about the cells nobody fed it.
    """
    cap = "NE(0000/1000)"
    ring = itacart.cell_to_boundary(cap, close=True)
    assert isinstance(ring, list)
    perimeter = 0.0
    for start, end in zip(ring, ring[1:]):
        distance, _ = itacart.inverse_geodesic(start[0], start[1], end[0], end[1])
        perimeter += distance
    area = itacart.effective_cell_area(cap)
    assert isinstance(area, float)
    assert round(perimeter, 1) == 3931.5
    assert round(_quotient(area, perimeter), 4) == 9.8696


def test_only_the_polar_caps_lack_a_simple_ring() -> None:
    """Two cells, two caps, and the reason named at the vertex pair.

    This said four cells and two caps covered twice over, and that was
    the leniency being read as a property of the grid. Only the eastern
    quadrants carry a cap: the cap sits in the meridian column, column
    zero is the meridian column, and west of the meridian the count
    starts at one. ``NW(0000/1000)`` and ``SW(0000/1000)`` name no cell,
    and what they used to report was the eastern cap's own geometry,
    reached by folding. The counts are pinned against the two that exist,
    and the two that do not are pinned as refusals so the fold cannot
    come back unnoticed.
    """
    northern = {"NE": (1, 2), "SE": (0, 2)}
    for quadrant, positions in northern.items():
        cap = f"{quadrant}(0000/1000)"
        ring = itacart.cell_to_boundary(cap)
        assert isinstance(ring, list)
        assert _repeated_vertex(ring) == positions
        first, second = ring[positions[0]], ring[positions[1]]
        assert first[1] == second[1]
        assert abs(first[0]) == abs(second[0]) == 180.0
        with pytest.raises(GeometryError, match="not a simple closed curve"):
            itacart.compactness(cap)

    north = itacart.effective_cell_area("NE(0000/1000)")
    assert isinstance(north, float)
    assert round(north, 3) == 12139402.004
    assert north == itacart.effective_cell_area("SE(0000/1000)")

    for quadrant in ("NW", "SW"):
        with pytest.raises(NonExistentCellError):
            itacart.effective_cell_area(f"{quadrant}(0000/1000)")


def test_the_repeated_vertex_predicate_passes_an_ordinary_ring() -> None:
    """The control that sees the refusal not fire."""
    ring = itacart.cell_to_boundary("NE(0500/0000)")
    assert isinstance(ring, list)
    assert _repeated_vertex(ring) is None
    assert _repeated_vertex([(180.0, 10.0), (-180.0, 10.0), (0.0, 0.0)]) == (0, 1)


# --------------------------------------------------------------------------
# The angle
# --------------------------------------------------------------------------


def test_the_angle_is_pinned_at_two_columns_of_row_zero() -> None:
    """Fixed by number at both ends of the sweep.

    The angle is a function of ``(X, Y)``: both cells below sit on the
    same row, at the same latitude, and disagree.
    """
    assert round(itacart.cell_base_angle("NE(0001/0000)"), 6) == 45.000035
    assert round(itacart.cell_base_angle("NE(2000/0000)"), 6) == 45.070984


def test_the_angle_grows_with_longitude_at_a_constant_rate() -> None:
    rates = []
    for column in (100, 500, 1000, 2000):
        cell = f"NE({column:04d}/0000)"
        anchor = itacart.cell_to_anchor(cell)
        assert isinstance(anchor, tuple)
        angle = itacart.cell_base_angle(cell)
        assert isinstance(angle, float)
        rates.append((angle - 45.0) / anchor[0])
    assert all(round(rate, 6) == 3.95e-4 for rate in rates)


def test_the_reference_is_the_parallel_and_not_a_geodesic_chord() -> None:
    """The control that separates the parallel from its chord.

    A cell's base runs along a parallel, and a parallel is not a
    geodesic: the geodesic between the two ends of the base leaves the
    anchor north of due east by half the convergence of the meridians
    over the base. Reading the base as a geodesic therefore biases the
    angle by that half-convergence, which is zero at the equator, where
    the parallel *is* a geodesic, and grows to almost three tenths of a
    degree at row 900.

    The predicate is the residual itself: the difference between the two
    readings has to equal ``(delta lambda / 2) * sin(phi)``, computed
    from the cell's own base. A test that only asserted the two differ
    would pass on any bias at all.
    """
    for cell, expected in (("NE(0500/0000)", 0.0), ("NE(0100/0900)", 0.283576)):
        ring = itacart.cell_to_boundary(cell)
        assert isinstance(ring, list)
        anchor, base = ring[0], ring[1]
        _, chord = itacart.inverse_geodesic(*anchor, *base)
        half_convergence = abs(90.0 - chord)
        span = base[0] - anchor[0]
        predicted = abs(span / 2.0 * math.sin(math.radians(anchor[1])))
        assert round(half_convergence, 6) == round(predicted, 6) == expected

    biased = 89.426478519 + 0.283576381
    assert round(biased, 6) == 89.710055
    assert round(itacart.cell_base_angle("NE(0100/0900)"), 6) == 89.426479


def test_the_angle_agrees_across_the_four_quadrants() -> None:
    """Enumerated, and now they do agree.

    They only agree because the reference is the parallel. Reading it
    off the ring's own edges splits the quadrants into two pairs,
    because mirroring reverses the ring and lands on the other of the
    two corners that carry 45 degrees on the plane.
    """
    # Pinned against an independent geodesic solution over the same two
    # points, not against this function's own output: the first two pins
    # once carried three units in the ninth decimal of the solver's
    # convergence residue, and agreed with themselves all the same.
    for column, row, expected in (
        (500, 0, 45.017729479),
        (100, 900, 89.426478519),
        (100, 950, 134.375798680),
    ):
        cell = f"({column:04d}/{row:04d})"
        values = {q: itacart.cell_base_angle(q + cell) for q in QUADRANTS}
        assert set(round(v, 9) for v in values.values()) == {expected}


def test_the_leaning_vertex_is_chosen_on_the_sinusoidal_plane() -> None:
    """And it has to be, because longitude reverses the sign near the pole.

    On the plane the four offsets are exactly ``0``, ``+l``, ``0`` and
    ``-l`` at every latitude. In longitude the meridians converge until
    the vertex the shear puts west of the anchor comes out east of it,
    which would pick the anchor's own position instead of the leaning
    side and report zero for a cell that leans 134 degrees.
    """
    from itacart.metrics import _local_offsets, _ring_of

    for row in (0, 300, 900, 950):
        cell = f"NE(0100/{row:04d})"
        ring = _ring_of(cell)
        anchor = itacart.cell_to_anchor(cell)
        assert isinstance(anchor, tuple)
        offsets = _local_offsets(cell, ring, anchor)
        assert [round(value) for value in offsets] == [0, 10000, 0, -10000]

    ring = _ring_of("NE(0100/0950)")
    anchor = itacart.cell_to_anchor("NE(0100/0950)")
    assert isinstance(anchor, tuple)
    in_longitude = [point[0] - anchor[0] for point in ring]
    assert min(in_longitude) == 0.0
    assert round(itacart.cell_base_angle("NE(0100/0950)"), 6) == 134.375799


def test_the_angle_reaches_the_obtuse_the_paper_shows_in_a_picture() -> None:
    """Figure 8 of the ITACaRT paper, as three numbers.

    Near the equidistant lines the cell keeps its 45 degrees; at mid
    latitude and mid longitude it is orthogonal; far from the origin it
    is severely obtuse. The paper says "exceeding 135 degrees" and the
    grid says 142.62.
    """
    italy = itacart.cell_base_angle(itacart.geo_to_cell(12.5, 42.0, 1))
    china = itacart.cell_base_angle(itacart.geo_to_cell(125.0, 30.0, 1))
    kamchatka = itacart.cell_base_angle(itacart.geo_to_cell(160.0, 56.0, 1))
    assert round(italy, 4) == 49.4921
    assert round(china, 4) == 95.2679
    assert round(kamchatka, 4) == 142.6159
    assert kamchatka > 135.0


def test_the_angle_passes_through_the_right_angle_and_keeps_opening() -> None:
    """Not monotone in latitude, and not bounded by 90.

    The cell is square near row 900 and obtuse above it. A reading
    folded into the acute half would report the same number on both
    sides of the square and lose the reversal entirely.
    """
    ladder = [
        ("NE(0100/0000)", 45.003545),
        ("NE(0100/0300)", 47.406155),
        ("NE(0100/0600)", 51.937956),
        ("NE(0100/0900)", 89.426479),
        ("NE(0100/0950)", 134.375799),
    ]
    values = []
    for cell, expected in ladder:
        angle = itacart.cell_base_angle(cell)
        assert isinstance(angle, float)
        assert round(angle, 6) == expected
        values.append(angle)
    assert values == sorted(values)
    assert values[-1] > 90.0


def test_the_trapezoid_is_measured_on_its_own_figure() -> None:
    """And the agreement with the neighbour is the result, not the rule.

    The clip cuts the two corners on the border side. Neither of them is
    the leaning vertex, so the angle is read off the corner the clip
    leaves standing. That it then matches the unclipped neighbour to
    four decimals is a measured consequence and not an analogy.
    """
    from itacart.metrics import _local_offsets, _ring_of

    clipped = "NE(2003/0000)"
    assert itacart.cell_shape(clipped) == "trapezoid"
    ring = _ring_of(clipped)
    anchor = itacart.cell_to_anchor(clipped)
    assert isinstance(anchor, tuple)
    offsets = _local_offsets(clipped, ring, anchor)
    assert [round(value) for value in offsets] == [0, 7508, 7484, -10000]

    assert round(itacart.cell_base_angle(clipped), 4) == 45.0711
    assert round(itacart.cell_base_angle("NE(2002/0000)"), 4) == 45.0711
    assert itacart.cell_base_angle(clipped) != itacart.cell_base_angle("NE(2002/0000)")


def test_the_triangle_zero_is_measured_not_declared() -> None:
    """The leaning vertex of a triangle is on the anchor's own parallel.

    The anchor of a meridian triangle is not a vertex: it sits at
    longitude 0 in the middle of the base, and the vertex at the most
    negative X is the other end of that base, at the same latitude. The
    side leaving the anchor in that direction *is* the parallel, so the
    angle against it is zero by measurement rather than by a clause.
    """
    from itacart.metrics import _local_offsets, _ring_of

    cell = "NE(0000/0300)"
    assert itacart.cell_shape(cell) == "triangle"
    ring = _ring_of(cell)
    anchor = itacart.cell_to_anchor(cell)
    assert isinstance(anchor, tuple)
    assert anchor[0] == 0.0
    assert anchor not in ring

    offsets = _local_offsets(cell, ring, anchor)
    assert [round(value) for value in offsets] == [-10000, 10000, 0]
    leaning = ring[min(range(len(ring)), key=lambda k: offsets[k])]
    assert leaning[1] == anchor[1]
    assert itacart.cell_base_angle(cell) == 0.0


def test_the_polar_cap_angle_needs_no_clause_of_its_own() -> None:
    """The cap is a triangle whose anchor is Greenwich, mid-base.

    And the anchor is not on the ring at all, which is the same defect
    that :func:`compactness` refuses: the boundary omits the limiting
    parallel and the anchor with it. The angle survives it, because the
    vertex at the most negative X is still on the anchor's parallel.
    """
    from itacart.metrics import _local_offsets, _ring_of

    for quadrant in EASTERN:
        cap = f"{quadrant}(0000/1000)"
        assert itacart.cell_shape(cap) == "triangle"
        anchor = itacart.cell_to_anchor(cap)
        assert isinstance(anchor, tuple)
        assert anchor[0] == 0.0
        assert round(abs(anchor[1]), 8) == 89.98240076
        ring = _ring_of(cap)
        assert anchor not in ring

        offsets = _local_offsets(cap, ring, anchor)
        leaning = ring[min(range(len(ring)), key=lambda k: offsets[k])]
        assert leaning[1] == anchor[1]
        assert itacart.cell_base_angle(cap) == 0.0


def test_the_output_forms_are_the_same_angle() -> None:
    cell = "NE(0500/0000)"
    degrees = itacart.cell_base_angle(cell)
    radians = itacart.cell_base_angle(cell, "radians")
    sine = itacart.cell_base_angle(cell, "sin")
    assert isinstance(degrees, float) and isinstance(radians, float)
    assert isinstance(sine, float)
    assert radians == pytest.approx(math.radians(degrees))
    assert sine == pytest.approx(math.sin(radians))
    assert round(sine, 9) == 0.707325553


def test_an_unknown_output_form_is_refused() -> None:
    with pytest.raises(GeometryError, match="must be 'degrees'"):
        itacart.cell_base_angle("NE(0500/0000)", "gradians")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# Normalized area
# --------------------------------------------------------------------------


def test_normalized_area_is_exactly_one_on_the_equal_area_interior() -> None:
    """Exactly, not nearly. The plane is equal-area by construction."""
    for quadrant in QUADRANTS:
        for column, row in ((500, 0), (100, 300), (100, 600), (100, 900)):
            cell = f"{quadrant}({column:04d}/{row:04d})"
            assert itacart.normalized_cell_area(cell) == 1.0


def test_normalized_area_uses_the_nominal_denominator() -> None:
    """The substitution is visible in the arithmetic, not only in prose."""
    cell = "NE(2003/0000)"
    effective = itacart.effective_cell_area(cell)
    assert isinstance(effective, float)
    nominal = itacart.nominal_cell_area(1)
    assert itacart.normalized_cell_area(cell) == effective / nominal
    assert round(effective / nominal, 6) == 1.249595


def test_the_border_families_carry_the_spread_and_are_enumerated() -> None:
    """The exceptions the equal-area claim has to name, with their counts.

    Over the window in all four quadrants at resolution 1. Every row of
    this table declares its ``n``, including the one that scores exactly
    1 and therefore contributes nothing to the spread.
    """
    families: dict[str, list[float]] = {}
    members: dict[str, list[str]] = {}
    interior = 0
    for cell in _enumerate_window():
        shape = itacart.cell_shape(cell)
        assert isinstance(shape, str)
        ratio = itacart.normalized_cell_area(cell)
        assert isinstance(ratio, float)
        if shape == "parallelogram":
            interior += 1
            assert ratio == 1.0
            continue
        families.setdefault(shape, []).append(ratio)
        members.setdefault(shape, []).append(cell)

    assert interior == 11988
    assert len(families["trapezoid"]) == 4000
    assert len(families["triangle"]) == 2002

    # The triangle count was 4004, and half of it was spelling rather than
    # grid: the western zero column folded onto its eastern twin and was
    # counted a second time. The number is decomposed here rather than
    # left bare, because a count whose scope is not said is a trap set for
    # whoever reads it next. Two thousand meridian triangles over rows 0
    # to 999 in the two eastern quadrants, plus the two polar caps at row
    # 1000, all of them in column zero.
    triangles_seen = members["triangle"]
    assert {cell[:2] for cell in triangles_seen} == set(EASTERN)
    assert {cell[3:7] for cell in triangles_seen} == {"0000"}
    caps = [cell for cell in triangles_seen if int(cell[8:12]) == 1000]
    assert len(caps) == 2
    assert len(triangles_seen) - len(caps) == 2000

    trapezoids = families["trapezoid"]
    assert round(min(trapezoids), 6) == 0.032263
    assert round(max(trapezoids), 6) == 2.067292
    assert round(statistics.fmean(trapezoids), 6) == 1.108982
    assert round(statistics.stdev(trapezoids), 6) == 0.455217

    triangles = families["triangle"]
    assert round(min(triangles), 6) == 0.121394
    assert max(triangles) == 1.0
    # 0.027760 was this statistic over a sample that counted every
    # meridian triangle twice. The spread itself did not move: every
    # phantom carried a ratio identical to its twin's, so the minimum and
    # the maximum are unchanged and the whole difference is the n-1
    # correction, sqrt(4003/4002) applied to the old figure. Measured
    # rather than argued, and it agrees to nine decimals.
    assert round(statistics.stdev(triangles), 6) == 0.027763


# --------------------------------------------------------------------------
# Compositional semantics
# --------------------------------------------------------------------------


def test_the_three_accept_a_compositional_index_in_decompose_order() -> None:
    index = "NE(0500/0000,0100/0300,0100/0900)"
    parts = itacart.decompose(index)
    assert len(parts) == 3
    for function in (
        itacart.compactness,
        itacart.cell_base_angle,
        itacart.normalized_cell_area,
    ):
        many = function(index)
        assert isinstance(many, list)
        assert many == [function(part) for part in parts]


def test_the_three_return_a_scalar_for_an_atomic_index() -> None:
    for function in (
        itacart.compactness,
        itacart.cell_base_angle,
        itacart.normalized_cell_area,
    ):
        assert isinstance(function("NE(0500/0000)"), float)


# --------------------------------------------------------------------------
# The uncertainty of a position written as a cell
# --------------------------------------------------------------------------

#: One cell of each family the uncertainty is read in: interior, the square
#: of row 900, a meridian triangle, a border trapezoid, the polar row, the
#: cap, an extension cell in Chukotka and a resolution-13 cell.
UNCERTAINTY_CELLS = (
    "NE(0500/0000)",
    "NE(0100/0900)",
    "NE(0000/0300)",
    "NE(2003/0000)",
    "SW(0001/0999)",
    "SE(0000/1000)",
    "NE(0931/0709)",
    "SW(0476/0260(4(A2(2(E3(4(B4(3(E2(1(C3(2(D4)))))))))))))",
)


def test_the_uncertainty_answers_in_metres_and_in_sides() -> None:
    """Radius and diameter in metres, the ratio as the radius over the side."""
    cell = "NE(0100/0900)"
    radius = itacart.cell_max_uncertainty(cell)
    assert radius == itacart.cell_max_uncertainty(cell, "radius")
    diameter = itacart.cell_max_uncertainty(cell, "diameter")
    ratio = itacart.cell_max_uncertainty(cell, "ratio")
    assert isinstance(radius, float) and isinstance(diameter, float)
    assert isinstance(ratio, float)
    assert ratio == radius / SIDE_1
    assert (round(radius, 3), round(diameter, 3)) == (7097.751, 14178.2)


def test_the_uncertainty_follows_decompose_over_several_roots() -> None:
    index = "NE(0500/0000,0100/0900),SE(0000/1000),NW(0003/0004(2))"
    parts = itacart.decompose(index)
    assert len(parts) == 4
    for output in ("radius", "diameter", "ratio"):
        many = itacart.cell_max_uncertainty(index, output)  # type: ignore[arg-type]
        assert isinstance(many, list)
        assert many == [
            itacart.cell_max_uncertainty(part, output)  # type: ignore[arg-type]
            for part in parts
        ]


def test_the_uncertainty_refuses_an_unknown_output_and_a_missing_cell() -> None:
    refusal = "must be 'radius', 'diameter' or 'ratio'"
    uncertainty = itacart.cell_max_uncertainty
    with pytest.raises(GeometryError, match=refusal):
        uncertainty("NE(0500/0000)", "diagonal")  # type: ignore[arg-type]
    for missing in ("NW(0000/1000)", "NE(0625/1000)", "NE(0500/0000),NW(0000/0110)"):
        with pytest.raises(NonExistentCellError):
            itacart.cell_max_uncertainty(missing)


def test_the_uncertainty_is_the_same_twice_and_in_a_fresh_process() -> None:
    """No state decides it: a second call, and a second interpreter.

    The second interpreter imports this same package and answers the cells
    in reverse order, so that a cache warmed in one order cannot hide a
    dependence on it.
    """
    import json
    import os
    import subprocess
    import sys

    index = itacart.compose(list(UNCERTAINTY_CELLS))
    assert sorted(itacart.decompose(index)) == sorted(UNCERTAINTY_CELLS)
    here = {
        output: [
            itacart.cell_max_uncertainty(cell, output) for cell in UNCERTAINTY_CELLS
        ]
        for output in ("radius", "diameter", "ratio")
    }
    again = {
        output: [
            itacart.cell_max_uncertainty(cell, output) for cell in UNCERTAINTY_CELLS
        ]
        for output in ("radius", "diameter", "ratio")
    }
    assert again == here
    source = os.path.dirname(os.path.dirname(itacart.__file__))
    script = (
        "import json, sys\n"
        "import itacart\n"
        "cells = json.load(sys.stdin)\n"
        "answers = {output: [itacart.cell_max_uncertainty(cell, output)\n"
        "                    for cell in reversed(cells)][::-1]\n"
        "           for output in ('radius', 'diameter', 'ratio')}\n"
        "json.dump({'file': itacart.__file__, 'answers': answers}, sys.stdout)\n"
    )
    path = os.pathsep.join(filter(None, (source, os.environ.get("PYTHONPATH"))))
    run = subprocess.run(
        [sys.executable, "-c", script],
        input=json.dumps(list(UNCERTAINTY_CELLS)),
        capture_output=True,
        text=True,
        check=True,
        env=dict(os.environ, PYTHONPATH=path),
    )
    there = json.loads(run.stdout)
    assert there["file"] == itacart.__file__
    assert there["answers"] == here
