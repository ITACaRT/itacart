"""Unit tests for :mod:`itacart.cells`.

Two measurement notes, both learned the hard way while writing these.

**Areas are measured locally.** The shoelace of a resolution-13 ring in
absolute plane coordinates is pure cancellation: the coordinates are
around 1e7 m and the side is 0.01 m, so the products are around 1e14 and
the answer around 1e-4. Every ring is therefore shifted to its own anchor
before its area is taken. Measuring wrong here would have looked exactly
like a broken cell.

**Orientation is measured on the plane, not on lon/lat.** The plane is
where the cell is a true parallelogram; the sign of the geodetic shoelace
happens to agree for small cells but that is a coincidence of scale, not
a property.
"""

from __future__ import annotations

import math

import pytest

from itacart import boundary, cells
from itacart.cells import FLOOR_EPSILON_M, _anchor_on_plane
from itacart.constants import RES1_MAX_INDEX
from itacart.exceptions import (
    DomainError,
    InvalidIndexError,
    NonExistentCellError,
    ResolutionError,
)
from itacart.geodesy import geodetic_to_sinusoidal
from itacart.index import decompose, join_components, split_components
from itacart.resolutions import cell_size, nominal_cell_area

RESOLUTIONS = tuple(range(1, 14))

# One interior position per quadrant, far from every axis.
QUADRANT_POINTS = {
    "NE": (13.0, 42.0),  # Italian peninsula
    "NW": (-73.9665, 40.7812),  # Central Park, Figure 7
    "SE": (151.2150784, -33.8567529),  # Sydney Opera House
    "SW": (-46.6328862, -23.5508962),  # Praca da Se
}


def local_shoelace(ring: list[tuple[float, float]]) -> float:
    """Signed area of a ring, computed relative to its first vertex."""
    origin_x, origin_y = ring[0]
    shifted = [(x - origin_x, y - origin_y) for x, y in ring]
    count = len(shifted)
    return 0.5 * math.fsum(
        shifted[i][0] * shifted[(i + 1) % count][1]
        - shifted[(i + 1) % count][0] * shifted[i][1]
        for i in range(count)
    )


def plane_ring(cell: str) -> list[tuple[float, float]]:
    """Geodetic boundary of a cell, projected back onto the plane."""
    ring = cells.cell_to_boundary(cell)
    return [geodetic_to_sinusoidal(lon, lat) for lon, lat in ring]


# ==========================================================================
# Round trip through the anchor
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_the_anchor_requantizes_to_its_own_cell(resolution: int, quadrant: str) -> None:
    """The property the 1 um floor epsilon exists to buy.

    Without the epsilon the anchor, which sits exactly on two grid lines,
    floors back into the previous cell and the identity fails at every
    resolution.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    anchor_lon, anchor_lat = cells.cell_to_anchor(cell)
    assert cells.geo_to_cell(anchor_lon, anchor_lat, resolution) == cell


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_projection_noise_stays_far_below_the_epsilon(
    resolution: int,
) -> None:
    """The margin, measured rather than assumed.

    The epsilon is only safe while the projection round trip is much
    quieter than it. This pins the ratio so that a future change to the
    geodesy that degrades the round trip fails here, loudly, instead of
    silently eating the margin.
    """
    worst = 0.0
    for lon, lat in QUADRANT_POINTS.values():
        cell = cells.geo_to_cell(lon, lat, resolution)
        exact_x, exact_y = _anchor_on_plane(cell)[:2]
        anchor_lon, anchor_lat = cells.cell_to_anchor(cell)
        back_x, back_y = geodetic_to_sinusoidal(anchor_lon, anchor_lat)
        worst = max(worst, math.hypot(back_x - exact_x, back_y - exact_y))
    assert worst < FLOOR_EPSILON_M / 10.0


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_the_reference_cells_are_not_on_a_quadrant_axis(
    resolution: int, quadrant: str
) -> None:
    """The scope of the anchor round trip, asserted rather than assumed.

    The round trip is stated for cells off the quadrant axes. If a
    reference point ever drifted onto an axis, the round-trip test above
    would start passing vacuously; this makes that impossible.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    assert cells.is_quadrant_boundary_cell(cell) is False


# --------------------------------------------------------------------------
# The quadrant axes, enumerated
# --------------------------------------------------------------------------

QUADRANT_CODES = ("NE", "NW", "SE", "SW")
ENUMERATED_COLUMNS = range(0, 400)


def enumerate_row(template: str) -> dict[str, int]:
    """Requantize the anchor of every cell in a base-index family.

    Enumeration, not sampling: the boundary families have effectively
    zero measure, so drawing random positions finds them by luck and
    reports a count that means nothing.
    """
    tally = {"round_trip": 0, "mirror_collision": 0, "non_existent": 0, "domain": 0}
    for quadrant in QUADRANT_CODES:
        for column in ENUMERATED_COLUMNS:
            cell = template.format(quadrant=quadrant, k=column)
            if boundary.is_valid_cell(cell) is False:
                tally["non_existent"] += 1
                continue
            try:
                anchor = cells.cell_to_anchor(cell)
                back = cells.geo_to_cell(anchor[0], anchor[1], 1)
            except DomainError:
                tally["domain"] += 1
            else:
                key = "round_trip" if back == cell else "mirror_collision"
                tally[key] += 1
    return tally


def test_the_equator_row_collides_only_in_the_south() -> None:
    """The anchor round trip holds except on the southern equator row.

    A vertex belongs to every cell that meets there, and the half-open
    convention awards it to one of them. The equator goes to the north,
    so a southern row ``0000`` cell holds an anchor it does not own and
    requantizes to its northern twin. Every southern cell in the row
    collides and every northern one round-trips; the count is exact
    because the family is enumerated rather than sampled.

    The western ``X = 0000`` cells are absent from both tallies: they do
    not exist, which is the paper's own remedy for the same collision on
    the prime meridian.
    """
    tally = enumerate_row("{quadrant}({k:04d}/0000)")
    assert tally == {
        "round_trip": 799,
        "mirror_collision": 799,
        "non_existent": 2,
        "domain": 0,
    }


def test_the_mirror_collision_is_the_meridian_case_without_its_remedy() -> None:
    """Why the equator cannot be repaired the way the meridian was.

    On the meridian the eastern triangle covers both sides, so deleting
    the western column loses no ground and the collision disappears. On
    the equator the two rows cover disjoint territory: deleting either
    would erase a ten-kilometre band of a hemisphere. The property is
    unsatisfiable for one of the two rows, and the predicate names which.
    """
    northern, southern = "NE(0518/0000)", "SE(0518/0000)"
    assert cells.cell_to_anchor(northern) == cells.cell_to_anchor(southern)
    assert cells.is_quadrant_boundary_cell(southern) is True
    assert cells.is_quadrant_boundary_cell(northern) is False

    triangle, absent = "NE(0000/0400)", "NW(0000/0400)"
    assert boundary.is_valid_cell(triangle) is True
    assert boundary.is_valid_cell(absent) is False
    anchor = cells.cell_to_anchor(triangle)
    assert cells.geo_to_cell(anchor[0], anchor[1], 1) == triangle


def test_the_prime_meridian_column_round_trips_at_every_row() -> None:
    """``X = 0000`` is a triangle in the east and nothing in the west.

    Its anchor is the midpoint of its base, which lies on the meridian,
    and the meridian is awarded to the east -- so the anchor requantizes
    to the cell that published it. The single collision is the origin
    cell of the southern quadrant, which sits on the equator as well and
    loses that tie for the reason the equator test gives.
    """
    tally = enumerate_row("{quadrant}(0000/{k:04d})")
    assert tally["non_existent"] == 800
    assert tally["round_trip"] == 799
    assert tally["mirror_collision"] == 1


def test_the_polar_row_holds_only_the_meridian_triangle() -> None:
    """Row ``Y = 1000`` is one cell wide, not empty and not a rectangle.

    At 89.98 degrees the parallel circle is about 6 km long, shorter than
    one 10 km cell, so no column but the meridian one has any ground
    under it. The triangle survives, clipped; every ``X >= 1`` is refused
    as non-existent rather than answered with a longitude off the planet.
    """
    tally = enumerate_row("{quadrant}({k:04d}/1000)")
    assert tally["non_existent"] == 1598
    assert tally["round_trip"] == 2
    assert tally["mirror_collision"] == 0
    assert boundary.is_valid_cell("NE(0000/1000)") is True
    assert boundary.is_valid_cell("NE(0001/1000)") is False
    assert boundary.is_equal_area_cell("NE(0000/1000)") is False


def test_the_resolution_one_index_space_is_not_a_rectangle() -> None:
    """``RES1_MAX_INDEX`` names a bounding-box corner, not a cell.

    The addressable ``X`` range shrinks with ``cos(phi)``: 2003 columns
    at the equator, 313 at row 900, none at row 1000. Table 1 gives the
    range as ``0000/0000`` to ``2003/1000``, which is the bounding box of
    a sinusoidal region rather than the region itself.

    The anchor is still *computable* -- the inverse shear and the
    inverse projection are total functions -- and it is the longitude
    that leaves the planet. What used to make this quiet was that the
    public entry point answered anyway; it no longer does, so both halves
    are asserted here. The interior computes, because nothing polices a
    non-canonical coordinate on its way to an answer, and the surface
    refuses, because the cell does not exist.
    """
    corner = f"NE({RES1_MAX_INDEX})"
    corner_lon, _ = boundary.to_geodetic(*_anchor_on_plane(corner)[:2])
    assert abs(corner_lon) > 180.0
    with pytest.raises(DomainError, match="outside"):
        cells.geo_to_cell(corner_lon, 89.98, 1)

    with pytest.raises(NonExistentCellError):
        cells.cell_to_anchor(corner)


def test_the_equator_restriction_is_a_declared_rule_not_a_pinned_bug() -> None:
    """The former ``xfail`` retired, replaced by the rule that explains it.

    It pinned a southern equator anchor failing to requantize to itself,
    as a defect awaiting repair. It is not a defect: a vertex is shared,
    the half-open convention awards it, and one of the two rows must
    lose. What was open was the choice of rule, and the rule is now
    stated -- equator to the north, meridian to the east -- with a public
    predicate naming exactly where it bites.
    """
    cell = "SW(0518/0000)"
    assert cells.is_quadrant_boundary_cell(cell) is True
    anchor = cells.cell_to_anchor(cell)
    assert anchor[1] == 0.0
    assert cells.geo_to_cell(anchor[0], anchor[1], 1) == "NW(0518/0000)"


def test_is_quadrant_boundary_cell_flags_the_southern_equator_row_alone() -> None:
    """The predicate narrowed once the meridian stopped being a problem.

    It used to flag both axes, because both carried a shared anchor. The
    meridian no longer does: the western column is gone and the eastern
    triangle owns its own base midpoint. What is left is the equator, and
    only its southern side.
    """
    assert cells.is_quadrant_boundary_cell("SW(0518/0000)") is True
    assert cells.is_quadrant_boundary_cell("SE(0518/0000)") is True
    assert cells.is_quadrant_boundary_cell("NE(0518/0000)") is False
    assert cells.is_quadrant_boundary_cell("NE(0000/0400)") is False
    assert cells.is_quadrant_boundary_cell("NE(0518/0400)") is False


def test_is_quadrant_boundary_cell_looks_at_the_anchor_not_the_base_pair() -> None:
    """A refined cell inherits the axis only if its own anchor sits on it.

    Child ``1`` of a ``Y = 0000`` cell keeps the parent's anchor and is
    still on the equator; child ``4`` is offset upward and is not.
    """
    assert cells.is_quadrant_boundary_cell("SW(0518/0000(1))") is True
    assert cells.is_quadrant_boundary_cell("SW(0518/0000(4))") is False


def test_is_quadrant_boundary_cell_is_vectorised() -> None:
    flags = cells.is_quadrant_boundary_cell("SW(0518/0000(1,4))")
    assert flags == [True, False]


def test_resolution_zero_has_no_anchor() -> None:
    """``"NE"`` is a valid index; the operation is what fails.

    Claiming ``InvalidIndexError`` would assert malformation, which
    the index grammar says is false. ``ResolutionError`` names exactly this case
    -- invalid for the operation -- and matches how
    :mod:`itacart.resolutions` already answers resolution 0.
    """
    with pytest.raises(ResolutionError, match="whole quadrant"):
        cells.cell_to_anchor("NE")
    assert not isinstance(
        pytest.raises(ResolutionError, cells.cell_to_anchor, "NE").value,
        InvalidIndexError,
    )


# ==========================================================================
# Four vertices, counter-clockwise
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_boundary_has_four_vertices(resolution: int, quadrant: str) -> None:
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    assert len(cells.cell_to_boundary(cell)) == 4


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_ring_is_counter_clockwise_in_every_quadrant(
    resolution: int, quadrant: str
) -> None:
    """The half of the counter-clockwise rule that mirroring can break.

    Reflection reverses orientation, so NW and SE -- one reflection each
    -- are where a naive re-signing yields a clockwise ring. SW composes
    two reflections and comes back correct on its own, which is why
    testing two quadrants would have missed it.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    assert local_shoelace(plane_ring(cell)) > 0.0


@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_the_ring_starts_at_the_anchor(quadrant: str) -> None:
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, 9)
    assert cells.cell_to_boundary(cell)[0] == cells.cell_to_anchor(cell)


def test_close_repeats_the_first_vertex() -> None:
    cell = cells.geo_to_cell(13.0, 42.0, 9)
    ring = cells.cell_to_boundary(cell, close=True)
    assert len(ring) == 5
    assert ring[0] == ring[-1]


@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_vertices_match_the_closed_form_of_the_paper(
    quadrant: str,
) -> None:
    """Section 3, Figure 1: ``(x,y)``, ``(x+l,y)``, ``(x,y+l)``, ``(x-l,y+l)``.

    Transcription, not convention: the paper gives the four vertices in
    closed form and its order is already counter-clockwise.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, 9)
    anchor_x, anchor_y, side, _ = _anchor_on_plane(cell)
    ring = [(x - anchor_x, y - anchor_y) for x, y in plane_ring(cell)]
    expected = [(0.0, 0.0), (side, 0.0), (0.0, side), (-side, side)]
    sign_x = -1.0 if quadrant[1] == "W" else 1.0
    sign_y = -1.0 if quadrant[0] == "S" else 1.0
    mirrored = [(sign_x * vx, sign_y * vy) for vx, vy in expected]
    if (quadrant[1] == "W") ^ (quadrant[0] == "S"):
        mirrored = [mirrored[0], *reversed(mirrored[1:])]
    for got, want in zip(ring, mirrored):
        assert got[0] == pytest.approx(want[0], abs=1e-6)
        assert got[1] == pytest.approx(want[1], abs=1e-6)


# ==========================================================================
# Mirroring across all four quadrants
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_mirrored_positions_give_mirrored_indices(
    resolution: int,
) -> None:
    """The rule stated literally in the paper, tested on all four.

    Reflection in the x axis for the southern quadrants, in the y axis
    for the western ones. A position and its three mirrors must give the
    same base pair and the same refinement codes, differing only in the
    quadrant code.
    """
    lon, lat = 20.0, 10.0
    reference = split_components(cells.geo_to_cell(lon, lat, resolution))
    for signed_lon, signed_lat, quadrant in (
        (lon, lat, "NE"),
        (-lon, lat, "NW"),
        (lon, -lat, "SE"),
        (-lon, -lat, "SW"),
    ):
        components = split_components(
            cells.geo_to_cell(signed_lon, signed_lat, resolution)
        )
        assert components[0] == quadrant
        assert components[1:] == reference[1:]


def test_the_southwest_composes_both_reflections() -> None:
    """SW is the case the origin's ``abs()`` gets right by accident.

    Two reflections restore orientation, so a bug that only checks
    orientation would pass here and fail in NW and SE. Stated as its own
    test so the composition is on the record.
    """
    cell = cells.geo_to_cell(-20.0, -10.0, 9)
    assert cell.startswith("SW")
    assert local_shoelace(plane_ring(cell)) > 0.0
    x, y, _, _ = _anchor_on_plane(cell)
    assert x < 0.0 and y < 0.0


# ==========================================================================
# Single traversal
# ==========================================================================


def test_descent_visits_each_level_exactly_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Structural, as the amendment asks. Not a clock reading.

    The origin's recursive version asked each ancestor for its
    coordinates from inside a loop over the ancestors, which is where the
    quadratic behaviour came from. Counting per-level work proves the
    traversal is single without measuring a machine.
    """
    calls: list[int] = []
    original = cells.linear_refinement_ratio

    def counting(resolution: int) -> int:
        calls.append(resolution)
        return original(resolution)

    monkeypatch.setattr(cells, "linear_refinement_ratio", counting)
    cells.geo_to_cell(13.0, 42.0, 13)
    assert calls == list(range(2, 14))


def test_ascent_visits_each_level_exactly_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[int] = []
    original = cells.linear_refinement_ratio

    def counting(resolution: int) -> int:
        calls.append(resolution)
        return original(resolution)

    cell = cells.geo_to_cell(13.0, 42.0, 13)
    monkeypatch.setattr(cells, "linear_refinement_ratio", counting)
    _anchor_on_plane(cell)
    assert calls == list(range(2, 14))

    # What the entry point costs on top, declared rather than folded into
    # the count above. Resolving whether a deep cell exists walks the same
    # ladder the ascent walks, so the public function pays for two passes
    # where the ascent itself pays for one. Linear in the depth either
    # way, which is the claim these counts make.
    calls.clear()
    cells.cell_to_anchor(cell)
    assert calls == list(range(2, 14)) * 2


def test_the_count_catches_a_descent_redone_from_the_root(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control for the descent count: the quadratic case fails it.

    Redoing the descent from the root at every level returns the same
    cell and asks for every earlier ratio again on the way down, which is
    how the origin's recursion grew with the square of the depth. Counting
    sees it in the calls. A clock cannot be relied on to: the fixed cost
    of each call dilutes the growth, and a ratio bound loose enough to
    survive the CI matrix admits the quadratic case.
    """
    calls: list[int] = []
    original_ratio = cells.linear_refinement_ratio
    original_descent = cells._descend_parallelogram

    def counting(resolution: int) -> int:
        calls.append(resolution)
        return original_ratio(resolution)

    def redone_from_the_root(
        residual_x: float,
        residual_y: float,
        side: float,
        first: int,
        resolution: int,
    ) -> list[str]:
        codes: list[str] = []
        for level in range(first, resolution + 1):
            codes = original_descent(residual_x, residual_y, side, first, level)
        return codes

    expected = cells.geo_to_cell(13.0, 42.0, 13)
    monkeypatch.setattr(cells, "linear_refinement_ratio", counting)
    monkeypatch.setattr(cells, "_descend_parallelogram", redone_from_the_root)
    assert cells.geo_to_cell(13.0, 42.0, 13) == expected
    assert len(calls) == sum(level - 1 for level in range(2, 14))


# ==========================================================================
# The centroid exists and is geodetic
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_centroid_is_a_geodetic_position(resolution: int, quadrant: str) -> None:
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    centroid_lon, centroid_lat = cells.cell_to_centroid(cell)
    assert -180.0 <= centroid_lon <= 180.0
    assert -90.0 <= centroid_lat <= 90.0


@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_centroid_sits_half_a_side_above_the_anchor(
    quadrant: str,
) -> None:
    """For a parallelogram the four-vertex mean collapses to this.

    The abscissa is unchanged because the shear leans the top edge back
    by exactly one side; the ordinate rises by half a side.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, 5)
    anchor_x, anchor_y, side, _ = _anchor_on_plane(cell)
    sign_y = -1.0 if quadrant[0] == "S" else 1.0
    centroid_x, centroid_y = geodetic_to_sinusoidal(*cells.cell_to_centroid(cell))
    assert centroid_x == pytest.approx(anchor_x, abs=1e-6)
    assert centroid_y == pytest.approx(anchor_y + sign_y * side / 2.0, abs=1e-6)


def test_centroid_differs_from_the_anchor() -> None:
    """The anchor and the centroid are different points, stated once."""
    cell = cells.geo_to_cell(13.0, 42.0, 5)
    assert cells.cell_to_centroid(cell) != cells.cell_to_anchor(cell)


# ==========================================================================
# Equal area, measured from the reprojected boundary
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
@pytest.mark.parametrize("quadrant", sorted(QUADRANT_POINTS))
def test_measured_area_matches_the_nominal_area(resolution: int, quadrant: str) -> None:
    """Under 0.01%, from the geodetic ring projected back onto the plane.

    This is the equal-area property end to end: the boundary goes out
    through the inverse projection and comes back through the forward
    one, and the area survives the trip.
    """
    lon, lat = QUADRANT_POINTS[quadrant]
    cell = cells.geo_to_cell(lon, lat, resolution)
    measured = local_shoelace(plane_ring(cell))
    nominal = nominal_cell_area(resolution)
    assert abs(measured / nominal - 1.0) < 1e-4


def test_area_is_the_same_at_the_equator_and_near_the_pole() -> None:
    """Equal area means equal area, not equal-area-where-convenient."""
    equator = cells.geo_to_cell(30.0, 0.5, 5)
    polar = cells.geo_to_cell(15.6469, 78.2232, 5)  # Longyearbyen
    area_equator = local_shoelace(plane_ring(equator))
    area_polar = local_shoelace(plane_ring(polar))
    assert area_equator == pytest.approx(area_polar, rel=1e-4)


# ==========================================================================
# The three orientations, pinned by name
# ==========================================================================


def test_orientation_1_the_y_index_grows_northward() -> None:
    """Measured, against the paper's prose, which reads the other way.

    Two positions on the same meridian: the northern one must carry the
    larger Y. The geodesy suite measures the same thing on the projection;
    this pins
    it on the index, which is what the boundary and topology modules consume.
    """
    south = split_components(cells.geo_to_cell(20.0, 10.0, 1))[1]
    north = split_components(cells.geo_to_cell(20.0, 20.0, 1))[1]
    assert int(north.split("/")[1]) > int(south.split("/")[1])


def test_orientation_2_codes_1_and_2_sit_on_the_southern_row() -> None:
    """The 2x2 quaternary layout of Figure 3(c).

    Two positions inside one resolution-1 cell, differing only in
    latitude: the southern one takes a code from ``{1, 2}`` and the
    northern one from ``{3, 4}``.
    """
    parent = cells.geo_to_cell(20.0, 10.0, 1)
    base_x, base_y = _anchor_on_plane(parent)[:2]
    side = cell_size(2)
    south_code = split_components(
        cells.sinusoidal_to_cell(base_x + side / 2, base_y + side / 4, 2)
    )[-1]
    north_code = split_components(
        cells.sinusoidal_to_cell(base_x + side / 2, base_y + side * 1.25, 2)
    )[-1]
    assert south_code in {"1", "2"}
    assert north_code in {"3", "4"}


def test_orientation_3_the_letter_grows_from_a_south_to_e_north() -> None:
    """The 5x5 quinary layout of Figure 3(d).

    Row ``A`` is the southern row. The vertical wrap-around ``E <-> A``
    that :mod:`itacart.topology` relies on inherits its sign from this.
    """
    parent = cells.geo_to_cell(20.0, 10.0, 2)
    base_x, base_y = _anchor_on_plane(parent)[:2]
    side = cell_size(3)
    codes = [
        split_components(
            cells.sinusoidal_to_cell(
                base_x + side * 2.5, base_y + side * (row + 0.5), 3
            )
        )[-1]
        for row in range(5)
    ]
    assert [code[0] for code in codes] == ["A", "B", "C", "D", "E"]


# ==========================================================================
# Vectorised semantics
# ==========================================================================


def test_atomic_index_returns_a_scalar_not_a_list_of_one() -> None:
    cell = cells.geo_to_cell(13.0, 42.0, 9)
    anchor = cells.cell_to_anchor(cell)
    assert isinstance(anchor, tuple)
    assert len(anchor) == 2


def test_compositional_index_returns_a_list_aligned_with_decompose() -> None:
    region = "NE(0500/0400(1(A1,A2,B3)))"
    anchors = cells.cell_to_anchor(region)
    assert isinstance(anchors, list)
    assert len(anchors) == len(decompose(region))
    for atom, anchor in zip(decompose(region), anchors):
        assert cells.cell_to_anchor(atom) == anchor


def test_boundary_of_a_compositional_index_is_a_list_of_rings() -> None:
    region = "NE(0500/0400(1(A1,A2)))"
    rings = cells.cell_to_boundary(region)
    assert len(rings) == 2
    assert all(len(ring) == 4 for ring in rings)


def test_cell_to_polygon_round_trips_through_shapely() -> None:
    cell = cells.geo_to_cell(13.0, 42.0, 9)
    polygon = cells.cell_to_polygon(cell)
    assert polygon.is_valid
    assert len(polygon.exterior.coords) == 5


def test_cell_to_polygon_of_a_compositional_index_is_a_list() -> None:
    polygons = cells.cell_to_polygon("NE(0500/0400(1(A1,A2)))")
    assert len(polygons) == 2
    assert all(polygon.is_valid for polygon in polygons)


# ==========================================================================
# The plane entry point
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_sinusoidal_to_cell_agrees_with_geo_to_cell(resolution: int) -> None:
    lon, lat = 13.0, 42.0
    x, y = geodetic_to_sinusoidal(lon, lat)
    assert cells.sinusoidal_to_cell(x, y, resolution) == cells.geo_to_cell(
        lon, lat, resolution
    )


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_cell_to_sinusoidal_inverts_the_plane_entry_point(resolution: int) -> None:
    cell = cells.geo_to_cell(13.0, 42.0, resolution)
    x, y = cells.cell_to_sinusoidal(cell)
    assert cells.sinusoidal_to_cell(x, y, resolution) == cell


# ==========================================================================
# Boundary scope and input validation
# ==========================================================================


def test_the_prime_meridian_itself_addresses_a_triangle() -> None:
    """The meridian is covered, not refused, and it is covered from the east."""
    cell = cells.geo_to_cell(0.0, 42.0, 9)
    assert boundary.cell_shape(cell) == "triangle"
    assert cell.startswith("NE(0000/")


@pytest.mark.parametrize("lon", [-0.01, 0.0, 0.01])
def test_both_sides_of_the_meridian_address_the_same_eastern_column(
    lon: float,
) -> None:
    """No western cell is created along the meridian, at any resolution.

    The triangle spans both sides, so a position just west of the line is
    part of the hierarchical subdivision of the eastern cell rather than
    a cell of its own.
    """
    cell = cells.geo_to_cell(lon, 42.0, 9)
    assert cell.startswith("NE(0000/")
    assert boundary.cell_shape(cell) in ("triangle", "parallelogram")


@pytest.mark.parametrize("lon", [179.9, -179.9, 180.0, -180.0])
def test_the_antemeridian_is_addressable_rather_than_refused(lon: float) -> None:
    """At 42 degrees north there is no extension, so the line is a hard edge.

    A half-degree band around it is not an exclusion zone: half a degree is
    the precision the extension limits are quoted to, and the paper asks
    for clipped cells here rather than none.
    """
    cell = cells.geo_to_cell(lon, 42.0, 9)
    assert boundary.is_valid_cell(cell) is True


def test_both_signs_of_one_hundred_and_eighty_name_the_same_cell() -> None:
    """One meridian, two spellings, awarded to the east like the prime one."""
    assert cells.geo_to_cell(-180.0, 42.0, 9) == cells.geo_to_cell(180.0, 42.0, 9)


@pytest.mark.parametrize("resolution", [0, 14, -1])
def test_quantization_rejects_resolutions_outside_one_to_thirteen(
    resolution: int,
) -> None:
    with pytest.raises(ResolutionError, match="outside 1"):
        cells.geo_to_cell(13.0, 42.0, resolution)


@pytest.mark.parametrize("bogus", [True, 9.0, "9", None])
def test_quantization_rejects_non_integer_resolutions(bogus: object) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        cells.geo_to_cell(13.0, 42.0, bogus)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "lon,lat",
    [
        (13.0, 91.0),
        (13.0, -91.0),
        (181.0, 42.0),
        (float("nan"), 42.0),
        (13.0, float("inf")),
    ],
)
def test_quantization_rejects_positions_outside_the_domain(
    lon: float, lat: float
) -> None:
    with pytest.raises(DomainError):
        cells.geo_to_cell(lon, lat, 9)


#: A position inside each extension zone, written west of the antemeridian.
ZONE_POSITIONS = {"FIJI": (-179.5, -18.5), "CHUKOTKA": (-172.0, 68.0)}


@pytest.mark.parametrize("resolution", [1, 5, 13])
@pytest.mark.parametrize("zone", sorted(ZONE_POSITIONS))
def test_a_centroid_written_past_180_quantizes_back_to_its_own_cell(
    zone: str, resolution: int
) -> None:
    """The package writes extension positions past 180; it reads them back.

    The first assertion keeps the test from passing over nothing: the
    centroid really is written past the antemeridian.
    """
    lon, lat = ZONE_POSITIONS[zone]
    cell = cells.geo_to_cell(lon, lat, resolution)
    centroid = cells.cell_to_centroid(cell)
    assert centroid[0] > 180.0
    assert cells.geo_to_cell(centroid[0], centroid[1], resolution) == cell


@pytest.mark.parametrize("zone", sorted(ZONE_POSITIONS))
def test_both_spellings_of_an_extension_position_name_one_cell(zone: str) -> None:
    lon, lat = ZONE_POSITIONS[zone]
    for resolution in (1, 9, 13):
        western = cells.geo_to_cell(lon, lat, resolution)
        assert cells.geo_to_cell(lon + 360.0, lat, resolution) == western


@pytest.mark.parametrize("zone, limit", [("FIJI", 182.0), ("CHUKOTKA", 190.5)])
def test_the_eastern_spelling_stops_at_the_meridian_the_zone_reaches(
    zone: str, limit: float
) -> None:
    lat = ZONE_POSITIONS[zone][1]
    assert boundary.is_valid_cell(cells.geo_to_cell(limit, lat, 5)) is True
    with pytest.raises(DomainError, match="outside"):
        cells.geo_to_cell(limit + 1e-4, lat, 5)


@pytest.mark.parametrize("lon, lat", [(181.0, 42.0), (181.0, -30.0), (181.0, 18.5)])
def test_a_longitude_past_180_outside_every_zone_band_is_refused(
    lon: float, lat: float
) -> None:
    """Not a general widening: Fiji's band on the wrong hemisphere is refused too."""
    with pytest.raises(DomainError, match="outside"):
        cells.geo_to_cell(lon, lat, 5)


@pytest.mark.parametrize("lon, lat", [(-181.0, -18.5), (-190.0, 68.0)])
def test_a_longitude_west_of_minus_180_is_refused_inside_a_zone_band(
    lon: float, lat: float
) -> None:
    """Only the eastern spelling is the package's; there is no western one."""
    with pytest.raises(DomainError, match="outside"):
        cells.geo_to_cell(lon, lat, 5)


#: A lattice over both caps, which begin at about 89.9823 degrees.
CAP_LATITUDES = (89.9824, 89.9868, 89.9912, 89.9956, 90.0)
CAP_LONGITUDES = tuple(-180.0 + 22.5 * k for k in range(17))


def _cap_positions() -> list[tuple[float, float, int]]:
    return [
        (lon, sign * lat, resolution)
        for resolution in range(2, 14)
        for lat in CAP_LATITUDES
        for lon in CAP_LONGITUDES
        for sign in (1.0, -1.0)
    ]


def test_quantization_in_the_caps_names_only_cells_that_exist() -> None:
    """Where the border carry meets the triangle descent, every answer exists.

    The carry moved the sheared column but the triangle descent reads the
    position, so a cap position in the strip beyond the last column was
    descended as if the strip had a cell of its own, and in a fine row that
    holds the pole as if a column other than the meridian survived there.
    """
    answers = [cells.geo_to_cell(lon, lat, res) for lon, lat, res in _cap_positions()]
    assert len(answers) == 2040
    assert all(boundary.is_valid_cell(answer) for answer in answers)


def test_the_cap_lattice_reaches_the_carry_at_every_resolution() -> None:
    """The control: the lattice does contain positions the carry acts on."""
    from itacart.cells import _resolve_extension
    from itacart.geodesy import geodetic_to_sinusoidal

    carried = {resolution: 0 for resolution in range(2, 14)}
    for lon, lat, resolution in _cap_positions():
        lon, quadrant = _resolve_extension(lon, lat)
        x, y = geodetic_to_sinusoidal(lon, lat)
        side = cells.cell_size(resolution)
        fine_row = int(math.floor((abs(y) + FLOOR_EPSILON_M) / side))
        fine_column = int(math.floor((abs(x) + abs(y) + FLOOR_EPSILON_M) / side))
        last = boundary.last_lattice_column(quadrant, fine_row, side)
        carried[resolution] += fine_column - fine_row > last
    assert all(count > 0 for count in carried.values()), carried


def test_a_cap_answer_holds_its_position_within_a_chord() -> None:
    """The answer contains the position, or lies within 0.1 mm of it.

    The tolerance is the chord: the cell's border side is a straight
    segment on the plane standing in for a curved meridian. Which cell is
    canonical where two existing cells overlap is not asserted here.
    """
    from shapely.geometry import Point, Polygon

    from itacart.cells import _resolve_extension
    from itacart.geodesy import geodetic_to_sinusoidal

    for lon, lat, resolution in _cap_positions():
        answer = cells.geo_to_cell(lon, lat, resolution)
        lon, _ = _resolve_extension(lon, lat)
        x, y = geodetic_to_sinusoidal(lon, lat)
        ring = boundary.plane_ring(answer)[1]
        origin = ring[0]
        polygon = Polygon([(a - origin[0], b - origin[1]) for a, b in ring])
        point = Point(x - origin[0], y - origin[1])
        assert polygon.contains(point) or polygon.exterior.distance(point) <= 1e-4, (
            answer,
            lon,
            lat,
        )


#: The population the cap repair was measured on: 25 latitudes from 89.9824
#: degrees to the pole, 49 longitudes 7.5 degrees apart, both hemispheres and
#: every resolution from 2 to 13 -- 29 400 positions. Before the repair, 2 064
#: of its answers named no cell, 414 distinct spellings among them, at every
#: resolution; every answer the repair changed was one of those.
CAP_GRID_LATITUDES = tuple(89.9824 + k * (90.0 - 89.9824) / 24 for k in range(25))
CAP_GRID_LONGITUDES = tuple(-180.0 + 7.5 * k for k in range(49))


def _cap_grid() -> list[tuple[float, float, int]]:
    return [
        (lon, sign * lat, resolution)
        for resolution in range(2, 14)
        for sign in (1.0, -1.0)
        for lat in CAP_GRID_LATITUDES
        for lon in CAP_GRID_LONGITUDES
    ]


@pytest.mark.slow
def test_every_cap_answer_on_the_measured_population_names_a_cell() -> None:
    """The repair holds over the whole population it was measured on."""
    answers = [cells.geo_to_cell(lon, lat, res) for lon, lat, res in _cap_grid()]
    assert len(answers) == 29_400
    assert all(boundary.is_valid_cell(answer) for answer in answers)


@pytest.mark.slow
def test_the_cap_answers_are_the_same_in_a_fresh_process() -> None:
    """Where existing cells overlap, the answer does not depend on the process.

    The canonical choice could in principle follow from state -- a cache
    warmed in another order, a set iterated in another order. A second
    interpreter, importing this same package, answers the population in
    reverse order and must give back every answer this one gives.
    """
    import json
    import os
    import subprocess
    import sys

    import itacart

    positions = _cap_grid()
    here = [cells.geo_to_cell(lon, lat, res) for lon, lat, res in positions]
    source = os.path.dirname(os.path.dirname(itacart.__file__))
    script = (
        "import json, sys\n"
        "import itacart\n"
        "positions = json.load(sys.stdin)\n"
        "answers = [itacart.geo_to_cell(lon, lat, res)\n"
        "           for lon, lat, res in reversed(positions)]\n"
        "json.dump({'file': itacart.__file__, 'answers': answers[::-1]}, sys.stdout)\n"
    )
    path = os.pathsep.join(filter(None, (source, os.environ.get("PYTHONPATH"))))
    run = subprocess.run(
        [sys.executable, "-c", script],
        input=json.dumps(positions),
        capture_output=True,
        text=True,
        check=True,
        env=dict(os.environ, PYTHONPATH=path),
    )
    there = json.loads(run.stdout)
    assert there["file"] == itacart.__file__
    assert there["answers"] == here


#: Spellings the cap quantizer answered on that population before the repair,
#: none of which names a cell: one per resolution, the middle one of that
#: resolution's list, the hemispheres alternating.
CORRECTED_CAP_SPELLINGS = (
    (-180.0, 89.9846, 2, "NE(0000/1000(2))"),
    (127.5, -89.99486666666667, 3, "SE(0000/1000(1(B3)))"),
    (-157.5, 89.99706666666667, 4, "NE(0000/1000(1(C2(4))))"),
    (-172.5, -89.99413333333334, 5, "SE(0000/1000(1(D2(1(D3)))))"),
    (-180.0, 89.99413333333334, 6, "NE(0000/1000(1(B4(1(D4(2))))))"),
    (-180.0, -89.99193333333334, 7, "SE(0000/1000(1(B4(2(A4(4(B5)))))))"),
    (180.0, 89.99119999999999, 8, "NE(0000/1000(1(A5(3(E1(4(D3(1))))))))"),
    (180.0, -89.99046666666666, 9, "SE(0000/1000(1(A5(3(E3(1(A5(2(A2)))))))))"),
    (-180.0, 89.99193333333334, 10, "NE(0000/1000(1(B4(2(A4(4(B5(2(E1(3))))))))))"),
    (180.0, -89.99193333333334, 11, "SE(0000/1000(1(B4(2(A4(4(B5(2(E1(3(C4)))))))))))"),
    (
        180.0,
        89.99046666666666,
        12,
        "NE(0000/1000(1(A5(3(E3(1(A5(2(A2(3(E2(1))))))))))))",
    ),
    (
        -180.0,
        -89.99193333333334,
        13,
        "SE(0000/1000(1(B4(2(A4(4(B5(2(E1(3(C4(1(D1)))))))))))))",
    ),
)


@pytest.mark.parametrize(
    "lon, lat, resolution, spelling",
    CORRECTED_CAP_SPELLINGS,
    ids=[f"resolution-{case[2]}" for case in CORRECTED_CAP_SPELLINGS],
)
def test_a_spelling_the_cap_quantizer_used_to_answer_is_not_answered_again(
    lon: float, lat: float, resolution: int, spelling: str
) -> None:
    """At the position that produced it, the answer now exists and is not it."""
    assert not boundary.is_valid_cell(spelling)
    answer = cells.geo_to_cell(lon, lat, resolution)
    assert boundary.is_valid_cell(answer)
    assert answer != spelling
    assert cells.geo_to_cell(lon, lat, resolution) == answer


@pytest.mark.parametrize("bogus", [float("nan"), float("inf")])
def test_the_plane_entry_point_rejects_non_finite_coordinates(bogus: float) -> None:
    with pytest.raises(DomainError, match="finite"):
        cells.sinusoidal_to_cell(bogus, 1000.0, 9)


# ==========================================================================
# The index assembler
# ==========================================================================


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_the_assembler_inverts_split_components(resolution: int) -> None:
    """The invariant that keeps the private helper honest."""
    cell = cells.geo_to_cell(13.0, 42.0, resolution)
    assert join_components(split_components(cell)) == cell


def test_the_assembler_renders_a_bare_quadrant() -> None:
    assert join_components(["NE"]) == "NE"


@pytest.mark.parametrize("resolution", RESOLUTIONS)
def test_every_emitted_index_is_valid_and_atomic(resolution: int) -> None:
    for lon, lat in QUADRANT_POINTS.values():
        cell = cells.geo_to_cell(lon, lat, resolution)
        assert cells._is_valid_atomic(cell)


def test_the_paper_example_parses_to_the_documented_shape(
    paper_example_index: str,
) -> None:
    """Section 3.1: ``SE(1400/0374(3(C2(3))))``, a resolution-4 cell."""
    assert split_components(paper_example_index) == [
        "SE",
        "1400/0374",
        "3",
        "C2",
        "3",
    ]
    assert cells._is_valid_atomic(paper_example_index)
