"""Tests for :mod:`itacart.geometry`.

Nothing is portable from ``itacart_core``: it has no coverage for this
module. The 406 figure that used to stand here was the size of its whole
suite, not a count of anything reusable.

Two families of test carry most of the weight. The containment chain is
asserted as a set inclusion rather than as a count comparison, because
three counts can nest while the sets do not. And the border families are
enumerated rather than sampled: the last column of a row is one cell in
two thousand, and a random draw finds it only by luck.
"""

from __future__ import annotations

import itertools
import math
from typing import Callable, Iterable

import pytest
from shapely.geometry import (
    GeometryCollection,
    LinearRing,
    LineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
    box,
)

import itacart
from itacart import boundary, geometry, hierarchy
from itacart.constants import MERIDIAN_QUADRANT, QUADRANTS
from itacart.exceptions import (
    AntemeridianError,
    DensificationError,
    DomainError,
    GeometryError,
    NonExistentCellError,
    ResolutionError,
    UnsupportedGeometryTypeError,
)

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

FIGURE_7A = (
    "NW(0625/0451(1(E1(3(B2(4(A2,B2,B3,B4,C2,C3,C4,C5,D1,D2,D3,D4,"
    "D5,E1,E2,E3,E4,E5)),B3(3(C1,D1,D2,E1,E2,E3,E4)),C2(1(A5,B5,C5,"
    "D4,D5,E4,E5),2,3(A4,A5,B3,B4,B5,C3,C4,C5,D2,D3,D4,D5,E2,E3,E4,"
    "E5),4),C3(1,2(A1,B1,B2,B3,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,"
    "E5),3,4),D2(1(A1,A2,A3,A4,A5,B2,B3,B4,B5,C3,C4,C5,D5),2,4(A3,"
    "A4,A5,B5)),D3(1,2,3(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,C2,C3,C4,C5,"
    "D4,D5),4),D4(1(D1,E1),3(A1,B1,C1,D1,E1)))))))"
)
"""Alternative (a) of Figure 7: a compositional fill at resolutions 6 and 7."""

FIGURE_7B = (
    "NW(0625/0451(1(E1(3(C3(2(C5(1(A2(3(C3(3(C4)))))))),"
    "B2(4(A2(1(D5(4(A3(1(D4)))))))),"
    "D2(1(B1(1(B2(4(C2(4(E5))))))),4(B3(2(A4(4(B4(2(D5)))))))),"
    "D3(3(E5(2(D5(1(D5(2(E3)))))))),"
    "D4(3(E1(4(D2(3(C3(4(B1))))))),1(E1(2(D2(4(C5(1(B1)))))))))))))"
)
"""Alternative (b) of Figure 7: vertex representation at resolution 13.

Seven vertices, each a chain from resolution 1 to 13 without a gap, so
each walks all twelve quaternary/quinary alternations. It is the only
artefact in the paper that does, which is why it is the cheapest
end-to-end check of the alphabet the package has.

Transcribed across a page break: the chain opens on the last line of
page 12 and closes on page 13, and a reader stopping at the page
boundary gets a truncated string that still balances its parentheses.
"""


@pytest.fixture
def parcel() -> Polygon:
    """A parcel inside the NW quadrant, clear of every border family."""
    return Polygon(
        [
            (-73.9812, 40.7681),
            (-73.9581, 40.8005),
            (-73.9497, 40.7968),
            (-73.9730, 40.7644),
        ]
    )


@pytest.fixture
def small_square() -> Polygon:
    """A square of a few hundred metres, cheap to fill at resolution 5."""
    return Polygon([(10.00, 45.00), (10.01, 45.00), (10.01, 45.01), (10.00, 45.01)])


# --------------------------------------------------------------------------
# The paper's Figure 7
# --------------------------------------------------------------------------


def test_figure_7a_is_a_valid_mixed_resolution_index() -> None:
    """Alternative (a) parses and holds cells at resolutions 6 and 7 only."""
    assert itacart.is_valid_index(FIGURE_7A)
    cells = itacart.decompose(FIGURE_7A)
    assert len(cells) == 114
    assert sorted({itacart.get_resolution(cell) for cell in cells}) == [6, 7]


def test_figure_7a_refills_from_its_own_footprint() -> None:
    """Filling the footprint of Figure 7(a) returns Figure 7(a)'s cells.

    The paper publishes an index, not the polygon behind it, so this is
    a fixed point -- ``fill(union(published)) == published`` -- and not
    a check that the author's parcel produces the author's index. The
    stronger claim is not available from the paper, and saying so is
    part of the evidence.
    """
    from shapely.ops import unary_union

    cells = itacart.decompose(FIGURE_7A)
    footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
    refilled = set(itacart.decompose(itacart.polyfill(footprint, 7)))
    expected = set(itacart.uncompact_cells(FIGURE_7A, 7))
    assert refilled == expected


def test_figure_7a_is_reproduced_as_a_string() -> None:
    """The compacted fill reproduces the published index character for character.

    Stronger than the set comparison above, and the strongest check the
    paper allows: the fill, run through the package's own compaction,
    lands on the exact 415-character string the article prints, so the
    package agrees with the authors on where to compact and not merely
    on which ground is covered.
    """
    from shapely.ops import unary_union

    cells = itacart.decompose(FIGURE_7A)
    footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
    assert itacart.compact_cells(itacart.polyfill(footprint, 7)) == FIGURE_7A


def test_the_two_compactions_agree_on_a_real_parcel() -> None:
    """``compact=True`` and ``compact_cells`` reach the same index.

    This test used to assert the opposite -- that the fill collapsed
    *less* than ``compact_cells`` -- and called the difference an open
    contract question. It was neither: the fill folded a node only when
    its square lay wholly inside the geometry, which is sufficient and
    not necessary, so the compaction was incomplete. Recording that
    shortfall as a property is how a limitation of the algorithm becomes
    a rule of the grid, and the old wording is kept here so that the
    substitution is visible rather than silent.

    What the old version measured is worth keeping: the excess it pinned
    was 495 square metres over 24 nodes, none of them near the
    floating-point noise floor, so the divergence was real rather than a
    sliver. It is now zero, and the round trips say the ground covered
    is unchanged.
    """
    parcel = Polygon(
        [
            (-73.9812, 40.7681),
            (-73.9581, 40.8005),
            (-73.9497, 40.7968),
            (-73.9730, 40.7644),
        ]
    )
    plain = itacart.polyfill(parcel, 7)
    geometric = itacart.polyfill(parcel, 7, compact=True)
    lexical = itacart.compact_cells(plain)

    assert geometric == lexical
    assert itacart.count_cells(geometric) < itacart.count_cells(plain)
    assert set(itacart.uncompact_cells(geometric, 7)) == set(itacart.decompose(plain))
    assert set(itacart.uncompact_cells(lexical, 7)) == set(itacart.decompose(plain))


@pytest.mark.parametrize("containment", ("intersects", "center", "contains"))
def test_compacting_equals_compacting_the_uniform_fill(containment: str) -> None:
    """``compact=True`` is post-hoc compaction, not a second semantics.

    A node folds when all of its children were kept, not only when it
    lies wholly inside the geometry. Wholly inside is sufficient and the
    descent still takes it as a shortcut; treating it as necessary left
    the compaction incomplete, at 15 per cent more cells under ``center``
    and 80 per cent more under ``intersects``, measured at resolution 3.

    Equality of the index, not of the cell count. Two indices can hold
    the same number of cells and describe different sets, so counting
    would measure a neighbour of the property.
    """
    regions = {
        "away from the axes": box(29.7, 9.7, 30.3, 10.3),
        "meridian, north": box(-0.22, 9.5, 0.22, 9.9),
        "meridian, south": box(-0.22, -9.9, 0.22, -9.5),
        "equator": box(29.7, -0.2, 30.3, 0.2),
        "the origin": box(-0.2, -0.2, 0.2, 0.2),
        "inside an extension zone": box(-179.8, -20.5, -179.2, -19.5),
    }
    for label, region in regions.items():
        for resolution in (2, 3):
            direct = itacart.polyfill(
                region, resolution, containment=containment, compact=True
            )
            afterwards = itacart.compact_cells(
                itacart.polyfill(region, resolution, containment=containment)
            )
            assert direct == afterwards, f"{label} at resolution {resolution}"


def test_compacting_actually_folds_across_a_discontinuity() -> None:
    """The control. Equality would also hold if nothing folded at all.

    Both sides of the equality above run through the same descent, so a
    fold that never fired would satisfy it and leave the fill uniform.
    """
    region = box(-0.22, 9.5, 0.22, 9.9)
    uniform = itacart.polyfill(region, 3, containment="intersects")
    compacted = itacart.polyfill(region, 3, containment="intersects", compact=True)
    assert itacart.count_cells(compacted) < itacart.count_cells(uniform)
    resolutions = {
        itacart.get_resolution(cell) for cell in itacart.decompose(compacted)
    }
    assert len(resolutions) > 1


def test_figure_7a_lies_in_the_nw_quadrant() -> None:
    """The fixture is mirrored, and the fill has to honour that."""
    assert all(cell.startswith("NW(") for cell in itacart.decompose(FIGURE_7A))


def test_figure_7b_walks_all_twelve_alternations() -> None:
    """Alternative (b) is a single chain from resolution 1 to 13."""
    cells = itacart.decompose(FIGURE_7B)
    assert cells
    assert all(itacart.get_resolution(cell) == 13 for cell in cells)


def test_figure_7b_alternates_quaternary_and_quinary() -> None:
    """Each level of the chain uses the alphabet its resolution declares."""
    for cell in itacart.decompose(FIGURE_7B):
        codes = itacart.split_components(cell)[2:]
        for offset, code in enumerate(codes):
            level = 2 + offset
            if itacart.linear_refinement_ratio(level) == 2:
                assert code in {"1", "2", "3", "4"}, (cell, level, code)
            else:
                assert code[0] in "ABCDE" and code[1] in "12345", (cell, level, code)


def test_figure_7b_round_trips_through_vertex_to_cell() -> None:
    """Rebuilding and requantizing (b)'s vertices returns the same cells."""
    cells = itacart.decompose(FIGURE_7B)
    ring = itacart.cells_to_geometry(cells, "Polygon")
    assert itacart.vertex_to_cell(ring, 13) == cells


# --------------------------------------------------------------------------
# The containment chain
# --------------------------------------------------------------------------


def test_containment_chain_is_a_set_inclusion(parcel: Polygon) -> None:
    """contains is inside center is inside intersects, as sets."""
    sets = {
        mode: set(itacart.decompose(itacart.polyfill(parcel, 6, containment=mode)))
        for mode in ("contains", "center", "intersects")
    }
    assert sets["contains"] < sets["center"] < sets["intersects"]


@pytest.mark.parametrize("resolution", [4, 5, 6])
def test_containment_chain_holds_at_several_resolutions(
    parcel: Polygon, resolution: int
) -> None:
    """The chain is a property of the descent, not of one resolution."""
    inner = set(
        itacart.decompose(itacart.polyfill(parcel, resolution, containment="contains"))
    )
    middle = set(
        itacart.decompose(itacart.polyfill(parcel, resolution, containment="center"))
    )
    outer = set(
        itacart.decompose(
            itacart.polyfill(parcel, resolution, containment="intersects")
        )
    )
    assert inner <= middle <= outer


def test_unknown_containment_mode_is_rejected(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="containment must be"):
        itacart.polyfill(parcel, 5, containment="touches")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# Count times nominal area against the equal-area plane
# --------------------------------------------------------------------------


def test_area_agreement_converges_under_refinement(parcel: Polygon) -> None:
    """The area error falls as the cell shrinks, for a fixed parcel.

    The reference is the polygon's area on the parallels plane, which is
    the equal-area image of its ellipsoidal area: that identity is the
    paper's central claim and the reason a cell count is an area at all.

    This is a regression property of the package, not a statement about
    how a resolution should be chosen. It varies the resolution against
    a fixed parcel, which is the opposite of the cadastral situation,
    where the resolution follows from the mapping scale and its
    positional tolerance and is fixed before the parcel is measured. The
    operational property is the bound in
    ``test_area_error_is_bounded_by_the_outline_band``.
    """
    reference = geometry._project(parcel).area
    errors = []
    for resolution in (5, 6, 7, 8):
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        errors.append(abs(area - reference) / reference)
    assert errors[-1] < errors[0] / 100.0


def test_area_error_is_bounded_by_the_outline_band(parcel: Polygon) -> None:
    """The error never exceeds the area of the cells the outline crosses.

    A cell the outline crosses is counted or dropped by which side its
    centre falls; every other cell is decided exactly. The total error is
    therefore bounded by the band of crossed cells, whose area is at most
    the perimeter times the cell side.

    This is the property a cadastral user needs, because it is the one
    that holds at a resolution they did not choose.
    """
    plane = geometry._project(parcel)
    reference = plane.area
    for resolution in (5, 6, 7):
        side = itacart.cell_size(resolution)
        band = plane.exterior.length * side
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        assert abs(area - reference) <= band


def test_the_bound_holds_at_a_prescribed_resolution_across_parcel_sizes() -> None:
    """At a fixed resolution the bound holds for every parcel size.

    The cadastral direction: the resolution is prescribed by the mapping
    scale, and the parcel's size and perimeter then determine the
    residual. Enumerated over five sizes spanning a factor of sixteen,
    at one resolution, because that is the axis the practitioner varies
    and the earlier tests do not.

    Only the bound is asserted. Where inside it a parcel lands is not
    predictable: measured over 24 phases at each of seven sizes, the
    RMS, median and ninetieth-percentile estimators of the scaling
    exponent give 0.61, -0.47 and 0.00. There is no power law to pin.
    """
    resolution = 7
    side = itacart.cell_size(resolution)
    for length in (200, 400, 800, 1600, 3200):
        delta_lat = length / 111_320.0
        delta_lon = delta_lat / math.cos(math.radians(45.0))
        parcel = Polygon(
            [
                (10.0, 45.0),
                (10.0 + delta_lon, 45.0),
                (10.0 + delta_lon, 45.0 + delta_lat),
                (10.0, 45.0 + delta_lat),
            ]
        )
        plane = geometry._project(parcel)
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        assert abs(area - plane.area) <= plane.exterior.length * side


def _plane_perimeter(plane: object) -> float:
    """Total ring length of a projected polygon, holes included."""
    return plane.exterior.length + sum(  # type: ignore[attr-defined]
        hole.length for hole in plane.interiors  # type: ignore[attr-defined]
    )


_BOUND_RADIUS_IN_SIDES = math.sqrt(5.0)
"""Radius of the outline neighbourhood, in cell sides.

The radius is the diameter of a cell, and the cell is not a square. On
the plane of the parallels it has base ``(s, 0)`` and side ``(-s, s)``,
so its longer diagonal is ``s*sqrt(5)`` and not the ``s*sqrt(2)`` a
square would give. The value is asserted against a measured cell by
``test_the_bound_radius_is_the_measured_cell_diameter`` so that it
cannot drift back to a hand-written number.
"""


def _residual_within_derived_bound(parcel: Polygon, resolution: int) -> bool:
    """Whether the count's area error respects the proved bound.

    The proof, in two lines. A cell whose closure does not meet the
    outline is decided exactly, so only cells meeting the outline can be
    misclassified. Every such cell lies inside the neighbourhood of the
    outline whose radius ``r`` is the cell diameter, and that
    neighbourhood has area at most ``2*r*P + pi*r*r``. Substituting
    ``r = s*sqrt(5)`` gives

        |n*a - A| <= 2*sqrt(5)*P*s + 5*pi*s*s

    which holds for any rectifiable outline, convex or not, holed or not.
    """
    plane = geometry._project(parcel)
    side = itacart.cell_size(resolution)
    perimeter = _plane_perimeter(plane)
    radius = _BOUND_RADIUS_IN_SIDES * side
    bound = 2.0 * radius * perimeter + math.pi * radius * radius
    counted = itacart.count_internal_cells(parcel, resolution)
    area = counted * itacart.nominal_cell_area(resolution)
    return abs(area - plane.area) <= bound


def _metres(value: float) -> float:
    """Degrees of latitude for a length in metres."""
    return value / 111_320.0


def _degrees_lon(value: float, latitude: float = 45.0) -> float:
    """Degrees of longitude for a length in metres at ``latitude``."""
    return _metres(value) / math.cos(math.radians(latitude))


def _extended_parcel_family() -> list[tuple[str, Polygon]]:
    """Parcels the rectangle family does not cover.

    The 168 parcels behind the residual measurement were all convex
    axis-aligned rectangles far larger than a cell, which is a narrow
    family to draw a bound from. These add a hole, a reentrant outline, a
    sliver, a parcel smaller than one cell, a triangle and a star.
    """
    shell = [
        (10.0, 45.0),
        (10.0 + _degrees_lon(800), 45.0),
        (10.0 + _degrees_lon(800), 45.0 + _metres(800)),
        (10.0, 45.0 + _metres(800)),
    ]
    hole = [
        (10.0 + _degrees_lon(300), 45.0 + _metres(300)),
        (10.0 + _degrees_lon(500), 45.0 + _metres(300)),
        (10.0 + _degrees_lon(500), 45.0 + _metres(500)),
        (10.0 + _degrees_lon(300), 45.0 + _metres(500)),
    ]
    star = []
    for step in range(12):
        angle = 2.0 * math.pi * step / 12.0
        radius = 400.0 if step % 2 == 0 else 150.0
        star.append(
            (
                10.0 + _degrees_lon(radius * math.cos(angle)),
                45.0 + _metres(radius * math.sin(angle)),
            )
        )
    return [
        ("holed", Polygon(shell, [hole])),
        (
            "reentrant",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(800), 45.0),
                    (10.0 + _degrees_lon(800), 45.0 + _metres(300)),
                    (10.0 + _degrees_lon(300), 45.0 + _metres(300)),
                    (10.0 + _degrees_lon(300), 45.0 + _metres(800)),
                    (10.0, 45.0 + _metres(800)),
                ]
            ),
        ),
        (
            "sliver",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(2000), 45.0),
                    (10.0 + _degrees_lon(2000), 45.0 + _metres(20)),
                    (10.0, 45.0 + _metres(20)),
                ]
            ),
        ),
        (
            "smaller than a cell",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(3), 45.0),
                    (10.0 + _degrees_lon(3), 45.0 + _metres(3)),
                    (10.0, 45.0 + _metres(3)),
                ]
            ),
        ),
        (
            "triangle",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(600), 45.0),
                    (10.0, 45.0 + _metres(600)),
                ]
            ),
        ),
        ("star", Polygon(star)),
    ]


@pytest.mark.parametrize(
    "name, parcel",
    _extended_parcel_family(),
    ids=[name for name, _ in _extended_parcel_family()],
)
def test_the_derived_bound_holds_beyond_the_rectangle_family(
    name: str, parcel: Polygon
) -> None:
    """The proved bound holds for holes, reentrant outlines and slivers."""
    assert _residual_within_derived_bound(parcel, 7), name


@pytest.mark.parametrize(
    "name, parcel",
    _extended_parcel_family(),
    ids=[name for name, _ in _extended_parcel_family()],
)
def test_the_tighter_observed_bound_also_holds(name: str, parcel: Polygon) -> None:
    """``P*s`` is tighter than the proof gives, and is not a theorem.

    The derived bound carries a factor of ``2*sqrt(5)`` and an ``s*s``
    term. ``P*s`` is about 4.5 times tighter and every parcel measured
    respects it -- rectangles, holes, slivers, and combs of four to
    sixty-four teeth one cell wide, worst ratio 0.34. No counterexample
    was found, which is not the same as a proof, and the distinction is
    recorded rather than blurred.
    """
    plane = geometry._project(parcel)
    side = itacart.cell_size(7)
    counted = itacart.count_internal_cells(parcel, 7)
    area = counted * itacart.nominal_cell_area(7)
    assert abs(area - plane.area) <= _plane_perimeter(plane) * side, name


def test_the_bound_radius_is_the_measured_cell_diameter() -> None:
    """The neighbourhood radius is measured off a cell, not off a square.

    ``s*sqrt(2)`` is the diagonal of a square and the cell is not one, so
    a bound built on it is understated by a factor of 1.58 and is not the
    bound the proof gives. Measuring the diameter here keeps the constant
    honest: writing it by hand again fails this test.
    """
    plane = geometry._project(Polygon(itacart.cell_to_boundary("NE(0500/0300)")))
    corners = list(plane.exterior.coords)[:-1]
    diameter = max(math.dist(a, b) for a, b in itertools.combinations(corners, 2))
    assert diameter / itacart.cell_size(1) == pytest.approx(
        _BOUND_RADIUS_IN_SIDES, rel=1e-9
    )
    assert _BOUND_RADIUS_IN_SIDES > math.sqrt(2.0)


def test_the_bound_says_nothing_about_a_parcel_narrower_than_a_cell() -> None:
    """Below about one cell across, the bound is true and useless.

    A three-metre parcel at ten-metre cells has a bound of thirteen times
    its own area, so it constrains nothing. The bound informs only while
    ``P*s < A``, which is the same as saying the parcel is more than
    roughly one cell wide. Worth pinning, because a bound quoted outside
    its useful range is how a specification becomes decoration.
    """
    tiny = dict(_extended_parcel_family())["smaller than a cell"]
    plane = geometry._project(tiny)
    side = itacart.cell_size(7)
    assert _plane_perimeter(plane) * side > plane.area
    assert _residual_within_derived_bound(tiny, 7)


def test_a_positional_tolerance_selects_one_resolution() -> None:
    """The ladder a tolerance picks from is decimal and unambiguous.

    The selection rule is the *coarsest* cell that does not exceed the
    tolerance, not the finest: every resolution below the tolerance
    satisfies it, and taking the finest would always land on resolution
    13 and over-resolve the survey. The tolerance itself comes from the
    applicable standard and is deliberately not encoded here.
    """
    for tolerance, expected in (
        (0.10, 11),
        (0.25, 11),
        (0.50, 10),
        (1.00, 9),
        (2.50, 9),
        (5.00, 8),
    ):
        chosen = min(
            r
            for r in range(1, itacart.MAX_RESOLUTION + 1)
            if itacart.cell_size(r) <= tolerance
        )
        assert chosen == expected, (tolerance, chosen)
        assert itacart.cell_size(chosen) <= tolerance
        assert itacart.cell_size(chosen - 1) > tolerance


def test_count_matches_the_centre_mode_fill(parcel: Polygon) -> None:
    """The fast path counts exactly what the slow path names.

    This is the test that keeps ``count_internal_cells`` honest: it is a
    different code path with a different accumulator, and the only thing
    tying it to ``polyfill`` is that they agree.
    """
    for resolution in (4, 5, 6):
        named = itacart.count_cells(
            itacart.polyfill(parcel, resolution, containment="center")
        )
        assert itacart.count_internal_cells(parcel, resolution) == named


# --------------------------------------------------------------------------
# The count does not materialise indices
# --------------------------------------------------------------------------


_PREVIOUSLY_LISTED_DOORS: tuple[tuple[str, str], ...] = (
    ("itacart.cells", "_quantize"),
    ("itacart.cells", "geo_to_cell"),
    ("itacart.cells", "sinusoidal_to_cell"),
    ("itacart.hierarchy", "_descend"),
    ("itacart.hierarchy", "_shift_column"),
    ("itacart.hierarchy", "_parent_cell"),
    ("itacart.hierarchy", "common_ancestor"),
    ("itacart.hierarchy", "compact_cells"),
    ("itacart.index", "_canonical_base"),
    ("itacart.index", "join_components"),
    ("itacart.index", "_render_node"),
    ("itacart.index", "_render_tree"),
    ("itacart.index", "compose"),
    ("itacart.index", "normalize"),
    ("itacart.topology", "_contact"),
    ("itacart.topology", "_eastern"),
    ("itacart.topology", "_spell"),
    ("itacart.topology", "_encode"),
)
"""The twenty doors the previous rule found, kept as a regression floor.

The rule that produced them matched a bare ``str`` return annotation in
four named modules. Widening it must not lose any of these, which is
what ``test_the_widening_dropped_no_previous_door`` checks.
"""


_DOORS_ON_THE_COUNT_PATH: frozenset[tuple[str, str]] = frozenset(
    {
        ("itacart.boundary", "_zone_of_row"),
        ("itacart.geometry", "_meridian_rows"),
        ("itacart.geometry", "_prepare"),
        ("itacart.geometry", "_quadrant_pieces"),
    }
)
"""The only routes to a ``str`` the count legitimately takes.

Measured rather than reasoned: everything the rule finds was armed at
once, and whatever fired was moved here until the count survived. All
four return a label -- an extension zone name, a quadrant name or a
hemisphere name -- and none of them names a cell, so the property under
test survives their being called.

Pinned in both directions, because a list of exemptions is exactly the
kind of thing that grows quietly. ``test_count_never_builds_an_index_fragment``
arms everything else and survives, so nothing else is on the path;
``test_every_exempted_door_is_really_on_the_count_path`` arms each of
these alone and requires it to fire, so nothing here is idle.
"""


def _index_rendering_doors() -> list[tuple[str, str]]:
    """Every route in the package from components to an index string.

    Discovered, not listed, and widened on the three axes where the
    previous rule was narrow.

    First, it matched ``hints.get("return") in ("str", str)``, so only a
    bare ``str`` counted. ``list[str]``, ``Iterator[str]``,
    ``str | None`` and ``tuple[str, ...]`` all escaped it, and among the
    escapees were ``index.decompose``, which materialises the whole list
    of index strings, and ``index.iter_cells``, which is by name the
    door a descent would use.

    Second, it skipped anything that was a class, and with the class
    went every method it holds.

    Third, the universe it swept was the set of module names appearing
    in its own list, so a module missing from the list was invisible to
    the check whose only job was to find omissions. The universe is now
    the package, walked with ``pkgutil``.
    """
    import importlib
    import inspect
    import pkgutil

    def mentions_str(annotation: object) -> bool:
        if annotation is str:
            return True
        return isinstance(annotation, str) and "str" in annotation

    modules = [itacart]
    for info in pkgutil.walk_packages(itacart.__path__, "itacart."):
        modules.append(importlib.import_module(info.name))

    found: set[tuple[str, str]] = set()
    for module in modules:
        module_name = module.__name__
        for attribute, value in list(vars(module).items()):
            members: list[tuple[str, object]] = []
            if isinstance(value, type):
                if getattr(value, "__module__", None) != module_name:
                    continue
                members.extend(
                    (f"{attribute}.{name}", member)
                    for name, member in vars(value).items()
                )
            elif callable(value) and getattr(value, "__module__", None) == module_name:
                members.append((attribute, value))
            for label, member in members:
                target: object = member
                if isinstance(member, property):
                    target = member.fget
                elif isinstance(member, (staticmethod, classmethod)):
                    target = member.__func__
                if not callable(target):
                    continue
                try:
                    hints = inspect.get_annotations(target, eval_str=False)
                except (TypeError, ValueError):  # pragma: no cover - defensive
                    continue
                if mentions_str(hints.get("return")):
                    found.add((module_name, label))
    return sorted(found)


def _arm_doors(doors: Iterable[tuple[str, str]], patch: pytest.MonkeyPatch) -> None:
    """Replace every named door with something that fails if it is called."""
    import importlib

    def detonator(module_name: str, label: str) -> Callable[..., str]:
        def fired(*args: object, **kwargs: object) -> str:
            raise AssertionError(f"the count reached {module_name}.{label}")

        return fired

    for module_name, label in doors:
        module = importlib.import_module(module_name)
        if "." in label:
            class_name, attribute = label.split(".", 1)
            owner = getattr(module, class_name)
            if attribute not in vars(owner):  # pragma: no cover - defensive
                continue
            patch.setattr(owner, attribute, detonator(module_name, label))
        else:
            patch.setattr(module, label, detonator(module_name, label))


def test_the_widening_dropped_no_previous_door() -> None:
    """Widening the rule may add doors; it may not lose one."""
    doors = set(_index_rendering_doors())
    missing = set(_PREVIOUSLY_LISTED_DOORS) - doors
    assert not missing, sorted(missing)


def test_the_door_rule_sees_what_a_bare_str_rule_cannot() -> None:
    """The rule fails if it narrows back to a bare ``str`` return.

    Each of these escaped the previous rule, and between them they pin
    all three widenings: the annotation shapes, the module universe, and
    the sheer size of the result. Twenty was never the number.
    """
    doors = set(_index_rendering_doors())
    for door in (
        ("itacart.index", "decompose"),
        ("itacart.index", "iter_cells"),
        ("itacart.index", "split_components"),
        ("itacart.hierarchy", "get_children"),
        ("itacart.hierarchy", "uncompact_cells"),
        ("itacart.hierarchy", "_descendants_of"),
        ("itacart.topology", "grid_disk"),
        ("itacart.topology", "_atoms"),
        ("itacart.boundary", "child_code"),
        ("itacart.boundary", "_side_quadrant"),
    ):
        assert door in doors, door
    assert len(doors) >= 100, len(doors)


def test_count_never_builds_an_index_fragment(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Counting names no cell, proved against every route to a string.

    Patching the descent's own helpers would prove nothing: the property
    is that no *package* route to an index string is reached, and a
    descent that rendered through another module would satisfy a guard
    aimed at itself.

    So every door the rule finds is armed at once, less the three that
    the count is measured to take and that return a label rather than a
    cell name. That is 106 detonators against the 20 this test used to
    arm, and the count has to survive with all of them live.
    """
    doors = [
        door
        for door in _index_rendering_doors()
        if door not in _DOORS_ON_THE_COUNT_PATH
    ]
    assert len(doors) >= 100, len(doors)
    _arm_doors(doors, monkeypatch)
    assert itacart.count_internal_cells(parcel, 7) == 33_258


def test_every_exempted_door_is_really_on_the_count_path(parcel: Polygon) -> None:
    """No exemption is idle: armed alone, each one fires.

    Without this, the exemption list is a place to hide a door that the
    count should not be taking.
    """
    for door in sorted(_DOORS_ON_THE_COUNT_PATH):
        with pytest.MonkeyPatch.context() as patch:
            _arm_doors([door], patch)
            with pytest.raises(AssertionError, match="the count reached"):
                itacart.count_internal_cells(parcel, 7)


def test_count_stack_depth_is_the_resolution_difference(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Recursion depth is bounded by the resolution, not by the cell count.

    The batched descent takes one frame per level and never one per node
    or per batch: a batch is descended to the target before the next is
    read. Batches of seven nodes make every level many batches wide, and
    the depth still stops at the resolution. The count does not move
    with the batch either, since it is a sum of integers.
    """
    levels: list[int] = []
    depth = [0, 0]  # current, deepest
    original = geometry._count_nodes_batched

    def _traced(
        lattice: object,
        u0: object,
        v0: object,
        side: float,
        level: int,
        target: int,
        strict: bool,
        shift: tuple[float, float],
    ) -> int:
        levels.append(level)
        depth[0] += 1
        depth[1] = max(depth[1], depth[0])
        try:
            return original(  # type: ignore[arg-type]
                lattice, u0, v0, side, level, target, strict, shift
            )
        finally:
            depth[0] -= 1

    whole = itacart.count_internal_cells(parcel, 6)
    monkeypatch.setattr(geometry, "_COUNT_BATCH", 7)
    monkeypatch.setattr(geometry, "_count_nodes_batched", _traced)
    assert itacart.count_internal_cells(parcel, 6) == whole
    assert min(levels) == 1
    assert max(levels) == 6
    assert depth[1] == 6
    assert len(levels) > 10 * depth[1]


def test_ordinary_cells_are_not_counted_one_node_at_a_time(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The public count takes the batched descent, and not the recursive one.

    Both give the same count, so a return to the node-by-node path would
    pass every count in the suite and show only on the clock. The path is
    asserted by calls instead of seconds. The parcel is clear of every
    border family, whose walks still descend node by node where a branch
    leaves the family; the reference descent, which takes the recursive
    path on purpose, is the control that the spy sees the calls it is
    there to see.
    """
    calls = [0]
    original = geometry._count_node

    def _spy(*args: object, **kwargs: object) -> int:
        calls[0] += 1
        return original(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(geometry, "_count_node", _spy)
    public = itacart.count_internal_cells(parcel, 6)
    assert calls[0] == 0
    reference = geometry._count_internal_cells_recursive(parcel, 6)
    assert calls[0] > 1_000
    assert public == reference


# --------------------------------------------------------------------------
# Densification
# --------------------------------------------------------------------------


def test_densify_segment_respects_the_threshold() -> None:
    points = itacart.densify_segment((0.0, 0.0), (1.0, 0.0), 10_000.0)
    for start, end in zip(points, points[1:]):
        assert itacart.inverse_geodesic(*start, *end)[0] <= 10_000.0


def test_densify_segment_spaces_points_equally_in_geodesic_distance() -> None:
    points = itacart.densify_segment((0.0, 0.0), (2.0, 0.0), 50_000.0)
    legs = [
        itacart.inverse_geodesic(*start, *end)[0]
        for start, end in zip(points, points[1:])
    ]
    assert len(legs) > 1
    assert max(legs) - min(legs) < 1e-6


def test_densify_segment_leaves_a_short_span_alone() -> None:
    span = ((0.0, 0.0), (0.001, 0.0))
    assert itacart.densify_segment(*span, 10_000.0) == [span[0], span[1]]


def test_densify_segment_handles_coincident_endpoints() -> None:
    assert itacart.densify_segment((3.0, 4.0), (3.0, 4.0), 10.0) == [
        (3.0, 4.0),
        (3.0, 4.0),
    ]


def test_densify_segment_rejects_an_unknown_edge_model() -> None:
    with pytest.raises(DensificationError, match="unknown edge model"):
        itacart.densify_segment((0.0, 0.0), (1.0, 0.0), 100.0, edge_model="LOXODROME")


@pytest.mark.parametrize("threshold", [0.0, -1.0, float("inf"), float("nan")])
def test_densify_segment_rejects_a_bad_threshold(threshold: float) -> None:
    with pytest.raises(DensificationError, match="positive finite"):
        itacart.densify_segment((0.0, 0.0), (1.0, 0.0), threshold)


def test_densify_segment_rejects_a_non_numeric_threshold() -> None:
    with pytest.raises(DensificationError, match="must be a number"):
        itacart.densify_segment(
            (0.0, 0.0), (1.0, 0.0), "1000"  # type: ignore[arg-type]
        )


def test_densify_orthodromic_is_idempotent() -> None:
    polygon = Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    once = itacart.densify_orthodromic(polygon, 20_000.0)
    twice = itacart.densify_orthodromic(once, 20_000.0)
    assert list(once.exterior.coords) == list(twice.exterior.coords)


def test_densify_orthodromic_densifies_holes_too() -> None:
    """The origin drops interior rings, and drops them silently."""
    shell = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    hole = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)]
    polygon = Polygon(shell, [hole])
    dense = itacart.densify_orthodromic(polygon, 20_000.0)
    assert len(dense.interiors) == 1
    assert len(dense.interiors[0].coords) > len(polygon.interiors[0].coords)


def test_densify_orthodromic_keeps_every_side_under_the_threshold() -> None:
    polygon = Polygon([(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (0.0, 2.0)])
    dense = itacart.densify_orthodromic(polygon, 30_000.0)
    coords = list(dense.exterior.coords)
    for start, end in zip(coords, coords[1:]):
        assert itacart.inverse_geodesic(*start, *end)[0] <= 30_000.0


def test_densify_orthodromic_rejects_a_non_polygon() -> None:
    with pytest.raises(TypeError, match="expected a shapely Polygon"):
        itacart.densify_orthodromic(
            LineString([(0, 0), (1, 1)])  # type: ignore[arg-type]
        )


def test_densify_orthodromic_passes_an_empty_polygon_through() -> None:
    empty = Polygon()
    assert itacart.densify_orthodromic(empty).is_empty


def test_densify_any_leaves_a_line_alone() -> None:
    line = LineString([(0.0, 0.0), (10.0, 0.0)])
    assert geometry._densify_any(line, 1000.0) is line


def test_densify_any_handles_a_multipolygon() -> None:
    parts = MultiPolygon(
        [
            Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            Polygon([(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]),
        ]
    )
    dense = geometry._densify_any(parts, 20_000.0)
    assert isinstance(dense, MultiPolygon)
    assert len(dense.geoms) == 2


def test_densify_ring_passes_a_degenerate_ring_through() -> None:
    assert geometry._densify_ring([(1.0, 2.0)], 100.0) == [(1.0, 2.0)]


def test_auto_segment_is_monotone_and_capped() -> None:
    thresholds = [geometry._auto_segment(r) for r in range(1, 14)]
    assert thresholds == sorted(thresholds, reverse=True)
    assert max(thresholds) <= 1000.0
    assert geometry._auto_segment(13) < 1000.0


# --------------------------------------------------------------------------
# Vertex order and winding
# --------------------------------------------------------------------------


def test_vertex_to_cell_preserves_ring_winding() -> None:
    """A ring and its reverse produce reversed cell sequences."""
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    forward = itacart.vertex_to_cell(Polygon(ring), 9)
    backward = itacart.vertex_to_cell(Polygon(ring[::-1]), 9)
    assert forward != backward
    assert set(forward) == set(backward)


def test_vertex_to_cell_puts_holes_after_the_exterior() -> None:
    shell = [(10.0, 45.0), (10.01, 45.0), (10.01, 45.01), (10.0, 45.01)]
    hole = [(10.004, 45.004), (10.006, 45.004), (10.006, 45.006), (10.004, 45.006)]
    cells = itacart.vertex_to_cell(Polygon(shell, [hole]), 9)
    assert len(cells) == 8
    assert cells[:4] == itacart.vertex_to_cell(Polygon(shell), 9)


def test_vertex_to_cell_handles_a_point() -> None:
    assert itacart.vertex_to_cell(Point(10.0, 45.0), 9) == [
        itacart.geo_to_cell(10.0, 45.0, 9)
    ]


def test_vertex_to_cell_handles_a_linestring() -> None:
    line = LineString([(10.0, 45.0), (10.01, 45.0), (10.02, 45.0)])
    assert len(itacart.vertex_to_cell(line, 9)) == 3


def test_vertex_to_cell_rejects_an_unsupported_type() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot take vertices"):
        itacart.vertex_to_cell(
            MultiPoint([(0, 0), (1, 1)]), 9  # type: ignore[arg-type]
        )


def test_vertex_to_cell_rejects_a_bad_resolution() -> None:
    with pytest.raises(ResolutionError):
        itacart.vertex_to_cell(Point(10.0, 45.0), 0)


def test_cells_to_geometry_inverts_vertex_to_cell_for_a_ring() -> None:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    cells = itacart.vertex_to_cell(Polygon(ring), 13)
    rebuilt = itacart.cells_to_geometry(cells, "Polygon")
    assert itacart.vertex_to_cell(rebuilt, 13) == cells


@pytest.mark.parametrize(
    "kind, count, expected",
    [
        ("Point", 1, Point),
        ("MultiPoint", 3, MultiPoint),
        ("LineString", 3, LineString),
        ("LinearRing", 3, LinearRing),
        ("Polygon", 3, Polygon),
    ],
)
def test_cells_to_geometry_builds_each_supported_type(
    kind: str, count: int, expected: type
) -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(count)]
    assert isinstance(itacart.cells_to_geometry(cells, kind), expected)


def test_cells_to_geometry_rejects_an_empty_list() -> None:
    with pytest.raises(GeometryError, match="empty cell list"):
        itacart.cells_to_geometry([])


def test_cells_to_geometry_rejects_a_composed_index() -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(2)]
    with pytest.raises(GeometryError, match="more than one cell"):
        itacart.cells_to_geometry([itacart.compose(cells)], "Point")


@pytest.mark.parametrize(
    "kind, count",
    [("Point", 2), ("LineString", 1), ("Polygon", 2), ("LinearRing", 2)],
)
def test_cells_to_geometry_rejects_a_count_the_type_cannot_use(
    kind: str, count: int
) -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(count)]
    with pytest.raises(GeometryError):
        itacart.cells_to_geometry(cells, kind)


def test_cells_to_geometry_rejects_an_unsupported_type() -> None:
    cells = [itacart.geo_to_cell(10.0, 45.0, 9)]
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot rebuild"):
        itacart.cells_to_geometry(cells, "Tetrahedron")


# --------------------------------------------------------------------------
# Consecutive collapse, non-consecutive kept
# --------------------------------------------------------------------------


def test_consecutive_duplicates_collapse() -> None:
    line = LineString([(10.0, 45.0), (10.0000001, 45.0), (10.01, 45.0)])
    assert len(itacart.vertex_to_cell(line, 9)) == 2


def test_non_consecutive_repeat_is_preserved() -> None:
    """A revisited vertex is a self-touching outline, not survey noise."""
    line = LineString([(10.0, 45.0), (10.01, 45.0), (10.0, 45.0)])
    cells = itacart.vertex_to_cell(line, 9)
    assert len(cells) == 3
    assert cells[0] == cells[2]


def test_dedupe_can_be_switched_off() -> None:
    line = LineString([(10.0, 45.0), (10.0000001, 45.0), (10.01, 45.0)])
    kept = itacart.vertex_to_cell(line, 9, dedupe_consecutive=False)
    assert len(kept) == 3


def test_ring_closing_duplicate_is_dropped_cyclically() -> None:
    """A ring whose last vertex shares the first cell loses the last one."""
    ring = [
        (10.0, 45.0),
        (10.01, 45.0),
        (10.01, 45.01),
        (10.0, 45.01),
        (10.0000001, 45.0),
    ]
    assert len(itacart.vertex_to_cell(Polygon(ring), 9)) == 4


def test_dedupe_consecutive_handles_an_empty_sequence() -> None:
    assert geometry._dedupe_consecutive([], cyclic=True) == []


def test_ring_cells_handles_a_degenerate_ring() -> None:
    assert geometry._ring_cells([], 9, True) == []


# --------------------------------------------------------------------------
# Canonical rings
# --------------------------------------------------------------------------


def _ring_cells_for_canonical() -> list[str]:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    return itacart.vertex_to_cell(Polygon(ring), 9)


def test_canonicalize_rings_is_idempotent() -> None:
    rings = [_ring_cells_for_canonical()]
    once = itacart.canonicalize_rings(rings)
    assert itacart.canonicalize_rings(once) == once


def test_canonicalize_rings_agrees_across_starting_vertices() -> None:
    ring = _ring_cells_for_canonical()
    rotations = [ring[k:] + ring[:k] for k in range(len(ring))]
    canonical = {tuple(itacart.canonicalize_rings([r])[0]) for r in rotations}
    assert len(canonical) == 1


def test_canonicalize_rings_starts_at_the_smallest_index() -> None:
    ring = _ring_cells_for_canonical()
    assert itacart.canonicalize_rings([ring])[0][0] == min(ring)


def test_canonicalize_rings_does_not_reverse_winding() -> None:
    """Direction distinguishes an exterior from a hole; rotation does not."""
    ring = _ring_cells_for_canonical()
    forward = itacart.canonicalize_rings([ring])[0]
    backward = itacart.canonicalize_rings([ring[::-1]])[0]
    assert forward != backward


def test_canonicalize_rings_orders_a_ring_that_repeats_its_minimum() -> None:
    """Booth compares whole rotations, so a repeated minimum still orders."""
    ring = ["A", "B", "A", "C"]
    assert itacart.canonicalize_rings([ring])[0] == ["A", "B", "A", "C"]


def test_canonicalize_rings_leaves_a_short_ring_alone() -> None:
    assert itacart.canonicalize_rings([["A"], []]) == [["A"], []]


def test_min_rotation_of_an_empty_ring_is_zero() -> None:
    assert geometry._min_rotation([]) == 0


def test_canonicalize_rings_preserves_ring_order() -> None:
    a = _ring_cells_for_canonical()
    b = list(reversed(a))
    out = itacart.canonicalize_rings([a, b])
    assert len(out) == 2
    assert set(out[0]) == set(a)


# --------------------------------------------------------------------------
# The antemeridian
# --------------------------------------------------------------------------


def test_geometry_crossing_the_antemeridian_is_refused() -> None:
    crossing = Polygon([(179.9, 10.0), (-179.9, 10.0), (-179.9, 10.1), (179.9, 10.1)])
    with pytest.raises(AntemeridianError, match="180 degrees"):
        itacart.polyfill(crossing, 5)


def test_a_geometry_near_but_not_across_the_line_is_not_refused() -> None:
    """Strictly astride, not merely close."""
    near = Polygon([(179.0, 10.0), (179.2, 10.0), (179.2, 10.2), (179.0, 10.2)])
    assert not itacart.crosses_antemeridian(near)


def test_the_antemeridian_test_is_strict_and_not_a_tolerance() -> None:
    """A geometry a thousandth of a degree from the line still does not cross.

    This is the case that separates the strict reading from a tolerance
    reading: a predicate that refused anything within some epsilon of 180
    would refuse this polygon, whose far edge sits 0.001 degrees short of
    the line. The strict predicate asks whether the outline is astride
    the meridian at the ordinate in question, and it is not.

    The earlier test does not discriminate: at 179.0 both readings agree.
    """
    hugging = Polygon(
        [
            (179.990, 10.0),
            (179.999, 10.0),
            (179.999, 10.01),
            (179.990, 10.01),
        ]
    )
    assert not itacart.crosses_antemeridian(hugging)

    astride = Polygon(
        [(179.999, 10.0), (-179.999, 10.0), (-179.999, 10.01), (179.999, 10.01)]
    )
    assert itacart.crosses_antemeridian(astride)


# --------------------------------------------------------------------------
# The refused border families, enumerated
# --------------------------------------------------------------------------


def test_prime_meridian_column_is_filled_from_the_east() -> None:
    """Column zero holds triangles, and the triangles are walked.

    They are not the square the ordinary descent tests, which is why
    they have a descent of their own. What this pins is the outcome:
    every cell the fill emits there exists, carries the eastern
    spelling that is the only one ``is_valid_cell`` admits, and is a
    triangle rather than something the square path happened to name.
    """
    astride = Polygon([(-0.01, 45.0), (0.01, 45.0), (0.01, 45.01), (-0.01, 45.01)])
    cells = itacart.decompose(itacart.polyfill(astride, 5))
    column_zero = [c for c in cells if c[3:7] == "0000"]
    assert column_zero
    assert all(itacart.is_valid_cell(c) for c in column_zero)
    assert all(c.startswith(("NE", "SE")) for c in column_zero)


def test_every_row_refuses_its_own_last_column() -> None:
    """Enumerated, not sampled: one cell in two thousand per row.

    The polar row is skipped because it is a row whose only column is
    its last one, so it has no interior column to accept and it raises
    for a different and more specific reason.
    """
    side = itacart.cell_size(1)
    checked = 0
    for row in range(0, 1000):
        last = itacart.last_lattice_column("NE", row, side)
        if last <= 1 or itacart.last_lattice_column("NE", row + 1, side) <= 0:
            continue
        with pytest.raises(NonExistentCellError, match="last lattice column"):
            geometry._check_addressable("NE", last, row)
        with pytest.raises(NonExistentCellError, match="last lattice column"):
            geometry._check_addressable("NE", last + 5, row)
        geometry._check_addressable("NE", last - 1, row)
        checked += 1
    assert checked == 999


def test_the_polar_row_is_the_row_the_pole_falls_in() -> None:
    """Which row is polar, derived from the pole rather than from a count.

    It used to be read off ``last_lattice_column`` as the last row whose
    answer is above zero, which made it row 999. That derivation
    saturates: the answer is zero for row 1000 and zero again for 1001,
    and zero means one cell in an eastern quadrant and none in a western
    one. Row 999 was called polar in all four quadrants on the strength
    of it, and its cells -- valid, six real children each, their highest
    point 1 966 metres below the pole -- were refused for a property they
    do not have.
    """
    side = itacart.cell_size(1)
    polar = int(math.floor(MERIDIAN_QUADRANT / side))
    assert polar == 1000

    with pytest.raises(DomainError, match="polar row"):
        geometry._check_addressable("NE", 0, polar)

    # The row below it is not polar, and the refusal it draws is the one
    # its last column earns, which the fill knows how to descend.
    with pytest.raises(NonExistentCellError, match="last lattice column"):
        geometry._check_addressable("NE", 1, polar - 1)
    geometry._check_addressable("NE", 0, polar - 1)

    # A western quadrant has no column zero, so its polar row addresses
    # nothing at all and says so before the polar branch is reached.
    with pytest.raises(DomainError, match="addresses no cell"):
        geometry._check_addressable("NW", 1, polar)


def test_a_row_above_the_pole_addresses_nothing() -> None:
    side = itacart.cell_size(1)
    beyond = int(math.floor(MERIDIAN_QUADRANT / side)) + 1
    with pytest.raises(DomainError, match="addresses no cell"):
        geometry._check_addressable("NE", 1, beyond)


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_screen_is_the_same_in_all_four_quadrants(quadrant: str) -> None:
    """The mirror flips the shear; it does not move the border families.

    Read on the last lattice column, which is a family the screen still
    owns. The prime-meridian column used to serve here and no longer
    can: it is filled rather than refused, and by a descent that never
    consults this screen.
    """
    side = itacart.cell_size(1)
    last = itacart.last_lattice_column(quadrant, 100, side)
    with pytest.raises(NonExistentCellError, match="last lattice column"):
        geometry._check_addressable(quadrant, last, 100)
    geometry._check_addressable(quadrant, last - 1, 100)


# --------------------------------------------------------------------------
# Quadrants and the shear
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_shear_has_unit_jacobian(quadrant: str) -> None:
    """Area on the lattice is area on the ellipsoid, in all four quadrants."""
    a, b, d, e = geometry._QUADRANT_SHEAR[quadrant]
    assert abs(a * e - b * d) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "lon, lat, quadrant",
    [
        (10.0, 45.0, "NE"),
        (-10.0, 45.0, "NW"),
        (10.0, -45.0, "SE"),
        (-10.0, -45.0, "SW"),
    ],
)
def test_a_fill_lands_in_its_own_quadrant(
    lon: float, lat: float, quadrant: str
) -> None:
    parcel = Polygon(
        [(lon, lat), (lon + 0.01, lat), (lon + 0.01, lat + 0.01), (lon, lat + 0.01)]
    )
    cells = itacart.decompose(itacart.polyfill(parcel, 5))
    assert all(cell.startswith(f"{quadrant}(") for cell in cells)


def test_the_lattice_shear_agrees_with_the_quantizer() -> None:
    """The shear used for filling is the one ``cells`` uses for addressing."""
    lon, lat = -10.0, 45.0
    x, y = itacart.geodetic_to_sinusoidal(lon, lat)
    a, b, d, e = geometry._QUADRANT_SHEAR["NW"]
    assert a * x + b * y == pytest.approx(abs(x) + abs(y))
    assert d * x + e * y == pytest.approx(abs(y))


def test_quadrant_pieces_drops_a_polygon_that_only_touches_an_axis() -> None:
    plane = geometry._project(
        Polygon([(0.0, 45.0), (0.01, 45.0), (0.01, 45.01), (0.0, 45.01)])
    )
    quadrants = {quadrant for quadrant, _ in geometry._quadrant_pieces(plane)}
    assert quadrants == {"NE"}


# --------------------------------------------------------------------------
# Descent internals
# --------------------------------------------------------------------------


def test_child_code_follows_the_quaternary_and_quinary_alphabets() -> None:
    assert geometry._child_code(0, 0, 2) == "1"
    assert geometry._child_code(1, 1, 2) == "4"
    assert geometry._child_code(0, 0, 3) == "A1"
    assert geometry._child_code(4, 4, 3) == "E5"


def test_child_code_inverts_the_quantizer_descent() -> None:
    """The code the fill writes is the code the quantizer would read."""
    cell = itacart.geo_to_cell(10.0, 45.0, 7)
    codes = itacart.split_components(cell)[2:]
    for offset, code in enumerate(codes):
        level = 2 + offset
        divisor = itacart.linear_refinement_ratio(level)
        found = [
            (row, column)
            for row in range(divisor)
            for column in range(divisor)
            if geometry._child_code(row, column, level) == code
        ]
        assert len(found) == 1


def test_leaves_between_is_the_product_of_the_ratios() -> None:
    assert geometry._leaves_between(6, 6) == 1
    assert geometry._leaves_between(6, 7) == 25
    assert geometry._leaves_between(5, 7) == 100


def test_expansion_names_every_descendant_once() -> None:
    """The shared expansion holds each descendant of a whole node once."""
    suffix = geometry._expansion(5, 7)
    codes = suffix.strip("()").replace("(", ",").replace(")", ",").split(",")
    leaves = [c for c in codes if c]
    assert len(leaves) == geometry._leaves_between(5, 7) + 4


def test_expansion_is_empty_at_the_target() -> None:
    assert geometry._expansion(7, 7) == ""


def test_expansion_is_shared_between_callers() -> None:
    """Two whole nodes at the same level reuse one string object.

    This is what keeps the interior of a fill cheap: a region wholly
    inside the geometry costs one pointer per node rather than one
    string per cell. Identity, not equality, is the property.
    """
    assert geometry._expansion(9, 11) is geometry._expansion(9, 11)


def test_budget_refuses_past_the_ceiling(monkeypatch: pytest.MonkeyPatch) -> None:
    """The tally is charged before expansion, never after the fact.

    Charging afterwards cannot work: expanding a resolution-1 node to
    resolution 13 builds a string of a trillion cells, and there is no
    counting that once it exists.
    """
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", 5)
    budget = geometry._Budget()
    budget.charge(5)
    assert budget.spent == 5
    with pytest.raises(GeometryError, match="exceeded 5 cells"):
        budget.charge(1)


def test_project_handles_each_supported_type() -> None:
    assert geometry._project(Point(10.0, 45.0)).geom_type == "Point"
    assert geometry._project(LineString([(10.0, 45.0), (11.0, 45.0)])).geom_type == (
        "LineString"
    )
    shell = [(10.0, 45.0), (11.0, 45.0), (11.0, 46.0)]
    hole = [(10.4, 45.2), (10.6, 45.2), (10.6, 45.4)]
    assert len(geometry._project(Polygon(shell, [hole])).interiors) == 1
    assert (
        geometry._project(
            MultiPolygon(
                [Polygon(shell), Polygon([(20.0, 45.0), (21.0, 45.0), (21.0, 46.0)])]
            )
        ).geom_type
        == "MultiPolygon"
    )


def test_project_rejects_an_unsupported_type() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot fill"):
        geometry._project(GeometryCollection([Point(0, 0)]))


# --------------------------------------------------------------------------
# Guards, arguments and edge cases
# --------------------------------------------------------------------------


@pytest.mark.parametrize("resolution", [0, 14, -1])
def test_polyfill_rejects_a_resolution_out_of_range(
    parcel: Polygon, resolution: int
) -> None:
    with pytest.raises(ResolutionError, match="outside"):
        itacart.polyfill(parcel, resolution)


def test_polyfill_rejects_a_non_integer_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        itacart.polyfill(parcel, 5.0)  # type: ignore[arg-type]


def test_polyfill_rejects_a_boolean_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        itacart.polyfill(parcel, True)  # type: ignore[arg-type]


@pytest.mark.parametrize("jobs", [0, -3])
def test_polyfill_rejects_a_bad_worker_count(parcel: Polygon, jobs: int) -> None:
    with pytest.raises(ValueError, match="n_jobs must be >= 1"):
        itacart.polyfill(parcel, 5, n_jobs=jobs)


def test_polyfill_rejects_a_non_integer_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs must be an int"):
        itacart.polyfill(parcel, 5, n_jobs=2.0)  # type: ignore[arg-type]


def test_polyfill_rejects_a_boolean_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs must be an int"):
        itacart.polyfill(parcel, 5, n_jobs=True)  # type: ignore[arg-type]


def test_polyfill_refuses_an_empty_geometry() -> None:
    """The refusal names emptiness, not a guess among possible causes."""
    with pytest.raises(
        GeometryError, match="covers no cell at this resolution: it is empty"
    ):
        itacart.polyfill(Polygon(), 5)


def test_polyfill_refuses_a_geometry_narrower_than_a_cell() -> None:
    """A sliver under ``contains`` names nothing, and says so."""
    sliver = Polygon(
        [(10.0, 45.0), (10.000001, 45.0), (10.000001, 45.000001), (10.0, 45.000001)]
    )
    with pytest.raises(GeometryError, match="narrower than one cell under 'contains'"):
        itacart.polyfill(sliver, 3, containment="contains")


def test_polyfill_honours_the_cell_ceiling(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", 10)
    with pytest.raises(GeometryError, match="exceeded 10 cells"):
        itacart.polyfill(parcel, 7)


def test_compact_returns_a_mixed_resolution_index(parcel: Polygon) -> None:
    compacted = itacart.polyfill(parcel, 7, compact=True)
    uniform = itacart.polyfill(parcel, 7)
    resolutions = {
        itacart.get_resolution(cell) for cell in itacart.decompose(compacted)
    }
    assert len(resolutions) > 1
    assert itacart.count_cells(compacted) < itacart.count_cells(uniform)


def test_compact_covers_the_same_ground(parcel: Polygon) -> None:
    compacted = itacart.polyfill(parcel, 7, compact=True)
    uniform = set(itacart.decompose(itacart.polyfill(parcel, 7)))
    expanded = set(itacart.uncompact_cells(compacted, 7))
    assert expanded == uniform


def test_the_worker_count_never_changes_the_answer(parcel: Polygon) -> None:
    """The worker count changes the schedule, never the answer."""
    serial = itacart.polyfill(parcel, 6, n_jobs=1)
    four = itacart.polyfill(parcel, 6, n_jobs=4)
    assert serial == four
    assert itacart.count_internal_cells(parcel, 6, n_jobs=4) == (
        itacart.count_internal_cells(parcel, 6, n_jobs=1)
    )


def test_polyfill_starts_no_thread_whatever_the_worker_count(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``n_jobs`` above 1 fills in the calling thread, with the serial answer.

    Two threads made the fill slower than one, so it runs in series for
    every worker count until it has a batched descent of its own. Every
    thread start is recorded, so a pool coming back fails here; the
    recorder is first shown to see the threads of a pool.
    """
    import threading
    from concurrent.futures import ThreadPoolExecutor

    started: list[str] = []
    real_start = threading.Thread.start

    def recording(self: threading.Thread) -> None:
        started.append(self.name)
        real_start(self)

    monkeypatch.setattr(threading.Thread, "start", recording)
    with ThreadPoolExecutor(max_workers=2) as pool:
        assert list(pool.map(abs, [-1, -2, -3])) == [1, 2, 3]
    assert started
    started.clear()
    serial = itacart.polyfill(parcel, 6)
    assert itacart.polyfill(parcel, 6, n_jobs=4) == serial
    assert started == []


def test_count_rejects_a_bad_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError):
        itacart.count_internal_cells(parcel, 0)


def test_count_rejects_a_bad_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs"):
        itacart.count_internal_cells(parcel, 5, n_jobs=0)


def test_count_of_an_empty_polygon_is_zero() -> None:
    assert itacart.count_internal_cells(Polygon(), 5) == 0


def test_polyfill_accepts_a_point(small_square: Polygon) -> None:
    """A point is a degenerate fill and still names its cell."""
    filled = itacart.polyfill(Point(10.0, 45.0), 5, containment="intersects")
    assert itacart.count_cells(filled) >= 1


def test_polyfill_accepts_a_linestring() -> None:
    line = LineString([(10.0, 45.0), (10.02, 45.02)])
    filled = itacart.polyfill(line, 5, containment="intersects")
    assert itacart.count_cells(filled) > 1


def test_polyfill_accepts_a_multipolygon() -> None:
    parts = MultiPolygon(
        [
            box(10.00, 45.00, 10.01, 45.01),
            box(10.05, 45.05, 10.06, 45.06),
        ]
    )
    assert itacart.count_cells(itacart.polyfill(parts, 5)) > 1


def test_polyfill_rejects_an_unsupported_geometry() -> None:
    with pytest.raises(UnsupportedGeometryTypeError):
        itacart.polyfill(GeometryCollection([Point(10.0, 45.0)]), 5)


def test_polyfill_output_parses_and_is_atomic_per_cell(parcel: Polygon) -> None:
    filled = itacart.polyfill(parcel, 6)
    assert itacart.is_valid_index(filled)
    for cell in itacart.decompose(filled):
        assert itacart.is_atomic(cell)
        assert itacart.get_resolution(cell) == 6


def test_polyfill_cells_are_real_cells(parcel: Polygon) -> None:
    """Every cell the fill names survives the boundary validator."""
    for cell in itacart.decompose(itacart.polyfill(parcel, 5)):
        assert itacart.is_valid_cell(cell)
        assert not itacart.absorbs_border(cell)


def test_filled_cells_requantize_to_themselves(parcel: Polygon) -> None:
    """A cell's own centroid addresses that cell, at the fill resolution."""
    for cell in itacart.decompose(itacart.polyfill(parcel, 6)):
        lon, lat = itacart.cell_to_centroid(cell)  # type: ignore[misc]
        assert itacart.geo_to_cell(lon, lat, 6) == cell


def test_intersects_covers_the_whole_geometry(small_square: Polygon) -> None:
    """Nothing of the polygon falls outside the intersects fill."""
    from shapely.ops import unary_union

    filled = itacart.polyfill(small_square, 5, containment="intersects")
    cover = unary_union(
        [itacart.cell_to_polygon(cell) for cell in itacart.decompose(filled)]
    )
    assert cover.buffer(1e-9).contains(small_square)


def test_contains_stays_inside_the_geometry(small_square: Polygon) -> None:
    filled = itacart.polyfill(small_square, 5, containment="contains")
    for cell in itacart.decompose(filled):
        assert small_square.buffer(1e-9).contains(itacart.cell_to_polygon(cell))


# --------------------------------------------------------------------------
# Public surface
# --------------------------------------------------------------------------


def test_every_public_name_is_exported() -> None:
    for name in geometry.__all__:
        assert hasattr(itacart, name), name
        assert name in itacart.__all__, name


def test_module_exports_exactly_what_the_phase_promised() -> None:
    assert sorted(geometry.__all__) == sorted(
        [
            "polyfill",
            "count_internal_cells",
            "vertex_to_cell",
            "cells_to_geometry",
            "densify_orthodromic",
            "densify_segment",
            "canonicalize_rings",
        ]
    )


def test_a_row_with_one_column_is_the_polar_row() -> None:
    """The polar row and the single-column row are the same row.

    This is why the polar test comes first in the screen: at resolution
    1 the polar row holds exactly one column, so a last-column test
    placed ahead of it would answer every polar position with the wrong
    diagnosis and the polar branch would be unreachable.
    """
    side = itacart.cell_size(1)
    single = [
        row for row in range(1100) if itacart.last_lattice_column("NE", row, side) == 1
    ]
    polar = [
        row
        for row in range(1100)
        if itacart.last_lattice_column("NE", row, side) > 0
        and itacart.last_lattice_column("NE", row + 1, side) <= 0
    ]
    assert single == polar


def test_min_rotation_orders_a_ring_with_a_repeated_prefix() -> None:
    """Booth's failure function is exercised by a ring with a period."""
    assert geometry._min_rotation(["B", "A", "B", "A"]) == 1
    assert geometry._min_rotation(["A", "A", "B"]) == 0
    assert geometry._min_rotation(["C", "A", "A", "B"]) == 1
    assert geometry._min_rotation(["B", "A", "A", "A"]) == 1


def test_min_rotation_restarts_mid_comparison() -> None:
    """A smaller vertex found part way through a match restarts the scan.

    ``ABA`` and ``BBA`` are the shortest rings that reach that branch:
    the scan is inside a partial match when it meets a vertex smaller
    than the one it was matching against, and has to rewind to it.
    """
    assert geometry._min_rotation(["A", "B", "A"]) == 2
    assert geometry._min_rotation(["B", "B", "A"]) == 2


def test_fill_is_canonical_without_a_sort(parcel: Polygon) -> None:
    """The descent emits siblings in the order the index is read in.

    ``compose`` preserves the order it is handed rather than imposing
    one, so the canonical ordering has to come from the descent itself.
    Asserted against the sorted-and-composed form the fill used to build
    by materialising every atomic index first.
    """
    for resolution in (4, 5, 6):
        for mode in ("contains", "center", "intersects"):
            filled = itacart.polyfill(parcel, resolution, containment=mode)
            assert filled == itacart.compose(sorted(set(itacart.decompose(filled))))


def test_fill_is_canonical_across_base_cells_and_quadrants() -> None:
    """Ordering holds where it is easiest to get wrong: between roots."""
    wide = Polygon([(9.99, 44.99), (10.30, 44.99), (10.30, 45.30), (9.99, 45.30)])
    filled = itacart.polyfill(wide, 3)
    assert filled == itacart.compose(sorted(set(itacart.decompose(filled))))
    assert filled.count("/") > 1  # more than one resolution-1 base cell


def test_a_wholly_contained_region_costs_far_less_than_its_cells(
    parcel: Polygon,
) -> None:
    """The index of a fill is a fraction of the size of its cell list.

    The measured ratio is what sets :data:`MAX_FILL_CELLS`: the composed
    index runs about three bytes per cell, while naming every cell as an
    atomic string runs about sixty. A fill that materialised the atomic
    form before composing it would exhaust memory two orders of
    magnitude earlier.
    """
    filled = itacart.polyfill(parcel, 7)
    cells = itacart.decompose(filled)
    atomic_bytes = sum(len(cell) for cell in cells)
    assert len(filled) * 10 < atomic_bytes


# --------------------------------------------------------------------------
# Sweep: measurements the fill's contract rests on
# --------------------------------------------------------------------------


def test_polyfill_never_consults_the_topology_caches(parcel: Polygon) -> None:
    """Filling asks topology nothing, at any resolution.

    Whether the 4,096-entry contact cache is large enough under a
    resolution-13 fill is an empty question: the descent subdivides in the
    sheared lattice and never asks for a neighbour, so both caches stay at
    zero hits, zero misses and zero entries.

    It is also why no asymmetry in ``get_neighbor`` can reach the fill: a
    fill that composes no step never meets one.
    """
    from itacart import topology

    caches = [
        value
        for value in vars(topology).values()
        if hasattr(value, "cache_info") and hasattr(value, "cache_clear")
    ]
    assert caches, "topology has no caches to check; the test has gone stale"
    for cache in caches:
        cache.cache_clear()
    itacart.polyfill(parcel, 7)
    itacart.count_internal_cells(parcel, 9)
    for cache in caches:
        info = cache.cache_info()
        assert (info.hits, info.misses, info.currsize) == (0, 0, 0), info


def test_children_union_matches_the_parent_only_away_from_the_border() -> None:
    """The union of children agrees with the parent in the interior only.

    A single figure for the whole grid -- an excess of 1.17e-3 -- says
    nothing about where it arises. Enumerated by family rather than
    sampled, the four separate cleanly: an interior
    parent sits at about -1.5e-7, a deficit rather than an excess,
    because each child approximates the curved parallel by the chord of
    its own height and the chords fall inside the parent's. The
    meridian column and the border-absorbing column run positive and
    grow toward the pole.

    This is why the containment chain is asserted over squares that carry
    this agreement; the border families are filled by walks of their own.
    """
    from shapely.ops import unary_union

    from itacart import hierarchy

    side = itacart.cell_size(1)
    for row in (0, 10, 100, 300, 500):
        last = itacart.last_lattice_column("NE", row, side)
        interior = f"NE({last // 2:04d}/{row:04d})"
        assert not itacart.absorbs_border(interior)
        parent = itacart.cell_to_polygon(interior)
        children = list(itacart.get_children(interior))[0]
        union = unary_union([itacart.cell_to_polygon(k) for k in children])
        relative = (union.area - parent.area) / parent.area
        assert -1e-5 < relative <= 0.0, (interior, relative)
        assert len(children) == itacart.refinement_ratio(2)
    assert hierarchy is not None


def test_border_children_survive_every_refinement() -> None:
    """``get_children`` refines a border cell down to resolution 13.

    The overlap filter behind the validity gate is not redundant, and its
    threshold has to shrink with the cell. A fixed square metre is the whole
    nominal area of a resolution-9 cell and larger than every cell below
    it; the comparison is strict, so a child wholly inside its parent would
    be rejected for being the size it is supposed to be, and the eastern
    border of the grid would have no refinement below resolution 8. The
    threshold is a fraction of the nominal area of the level.

    Two assertions, and the second is the guard against overcorrection.
    Every level has to yield children, and the counts have to keep
    deviating from the canonical refinement ratio where the geometry
    says they should: a threshold pushed too far down would admit
    candidates the parent does not hold and iron that deviation flat.
    The measured series for row 100, resolutions 1 to 12, is
    ``2, 21, 4, 26, 5, 11, 5, 9, 2, 26, 4, 21``. It was re-measured when
    the child enumeration stopped assuming that a child is spelled under
    the parent's own stem or the column immediately east: the first three
    entries are unchanged and the rest moved, because the old walk lost
    children the deeper it went. The series before that correction read
    ``2, 21, 4, 24, 4, 23, 4, 25, 4, 25, 4, 24``, and its own first seven
    entries were what a still earlier fix had produced.
    """
    side = itacart.cell_size(1)
    row = 100
    cell = f"NE({itacart.last_lattice_column('NE', row, side):04d}/{row:04d})"
    assert itacart.absorbs_border(cell)
    counts: dict[int, int] = {}
    while itacart.get_resolution(cell) < 13:
        resolution = itacart.get_resolution(cell)
        children = list(itacart.get_children(cell))[0]
        counts[resolution] = len(children)
        absorbing = [k for k in children if itacart.absorbs_border(k)]
        assert absorbing, (resolution, children)
        cell = absorbing[-1]
    assert sorted(counts) == list(range(1, 13)), counts
    assert [counts[r] for r in sorted(counts)] == [
        2,
        21,
        4,
        26,
        5,
        11,
        5,
        9,
        2,
        26,
        4,
        21,
    ], counts
    assert all(count > 0 for count in counts.values()), counts
    deviating = [r for r in counts if counts[r] != itacart.refinement_ratio(r + 1)]
    assert deviating, counts


# --------------------------------------------------------------------------
# Filling inside an extension zone
# --------------------------------------------------------------------------


_ZONE_FOOTPRINTS: tuple[tuple[str, float, float], ...] = (
    ("Fiji, southern hemisphere", -18.0, -17.8),
    ("Chukotka, northern hemisphere", 65.0, 65.2),
)
"""One latitude band inside each extension zone, as ``(name, low, high)``.

Fiji is realised on resolution-1 rows 171 to 237 and reaches 182 degrees;
Chukotka on rows 709 to 799, reaching 190.5. Both bands sit well inside
their zone, so a footprint written past 180 there is inside the domain
and has to fill rather than be refused.
"""


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_a_footprint_astride_the_antemeridian_fills_inside_a_zone(
    name: str, low: float, high: float
) -> None:
    """The positive case of the antemeridian: a zone footprint astride it.

    A ring written from 179.9 to 180.3 lies on one side of the line
    already, so it does not cross: inside an extension zone the domain
    reaches past 180 and the whole footprint is addressable. It used to
    raise a bare ``GEOSException`` from the geometry engine instead,
    because densification normalised its interior vertices to
    ``(-180, 180]`` and folded the ring back across the globe.
    """
    footprint = Polygon([(179.9, low), (180.3, low), (180.3, high), (179.9, high)])
    assert not itacart.crosses_antemeridian(footprint)
    cells = itacart.decompose(itacart.polyfill(footprint, 3))
    assert cells, name
    assert all(itacart.get_resolution(cell) == 3 for cell in cells)


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_a_footprint_wholly_past_180_fills_inside_a_zone(
    name: str, low: float, high: float
) -> None:
    """Written entirely past the line, and still inside the domain.

    The strict reading of the antemeridian test takes the longitude
    literally, so this footprint does not cross either. It broke for the
    same reason and is fixed by the same change.
    """
    footprint = Polygon([(180.2, low), (180.4, low), (180.4, high), (180.2, high)])
    assert not itacart.crosses_antemeridian(footprint)
    assert itacart.decompose(itacart.polyfill(footprint, 3)), name


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_the_wrapped_spelling_of_the_same_footprint_is_still_refused(
    name: str, low: float, high: float
) -> None:
    """The other direction of travel, which is a crossing and stays one.

    Spelled with a jump from 179.9 to -179.7 the ring does cross, and the
    refusal is the package's own exception. Filling one spelling must not
    quietly start filling the other: they describe different rings.
    """
    wrapped = Polygon([(179.9, low), (-179.7, low), (-179.7, high), (179.9, high)])
    assert itacart.crosses_antemeridian(wrapped)
    with pytest.raises(AntemeridianError):
        itacart.polyfill(wrapped, 3)


def test_a_footprint_past_180_outside_every_zone_is_refused_by_name() -> None:
    """Outside a zone the domain stops at the line, and so does the fill.

    What matters here is not that it is refused but that the refusal is
    the package's, not the geometry engine's. A raw ``GEOSException``
    tells the caller nothing they can act on.

    The name used to be ``NonExistentCellError``, and it was an accident:
    such a footprint projects into the western quadrant at a column past
    the end of its row, and the fill refused that column outright, so the
    guarantee about the line rested on a limitation about absorbing
    cells. The fill descends that column now. The refusal is the
    antemeridian's own, and the zoned rows are asserted alongside it so
    that a guard which refused everything past the line would not pass.
    """
    footprint = Polygon([(179.9, 5.0), (180.3, 5.0), (180.3, 5.2), (179.9, 5.2)])
    assert not itacart.crosses_antemeridian(footprint)
    with pytest.raises(AntemeridianError):
        itacart.polyfill(footprint, 3)

    beyond_the_reach = Polygon(
        [(190.0, 68.0), (191.0, 68.0), (191.0, 68.2), (190.0, 68.2)]
    )
    with pytest.raises(AntemeridianError):
        itacart.polyfill(beyond_the_reach, 3)

    for zoned in (
        Polygon([(179.9, 68.0), (182.0, 68.0), (182.0, 68.2), (179.9, 68.2)]),
        Polygon([(179.9, -18.0), (181.0, -18.0), (181.0, -17.8), (179.9, -17.8)]),
    ):
        cells = itacart.decompose(itacart.polyfill(zoned, 3))
        assert cells
        assert all(itacart.is_valid_cell(cell) for cell in cells)


def test_densification_keeps_the_longitude_branch_it_was_given() -> None:
    """The mechanism behind the four tests above, isolated.

    ``direct_geodesic`` normalises to ``(-180, 180]``, which is correct
    for the question it answers and wrong for densification: the
    interior vertices of a segment crossing the line came back on the
    other branch and the ring self-intersected. Measured on the ring
    itself rather than through the fill, so that a future change to the
    fill cannot hide a regression here.
    """
    ring = Polygon([(179.9, -18.0), (180.3, -18.0), (180.3, -17.8), (179.9, -17.8)])
    dense = itacart.densify_orthodromic(ring, 1000.0)
    longitudes = [point[0] for point in dense.exterior.coords]
    assert min(longitudes) > 179.0, min(longitudes)
    assert max(longitudes) <= 180.3 + 1e-9, max(longitudes)
    assert dense.is_valid


def test_an_outline_that_crosses_itself_is_refused_with_our_own_name() -> None:
    """A fold is named by this package, at the entry and at the split.

    A caller's bowtie is refused at the entrance, as an invalid areal
    geometry, by both public entry points: they share the machinery, and a
    future split of it should not quietly leave one of them unguarded.

    The quadrant split keeps its own name for a fold as well. The engine
    refuses a plane figure with no interior with an exception this package
    does not own, and that exception once reached callers through a public
    name, carried by a polar-row cell whose densified boundary folded over
    the pole -- a defect of this package, repaired at the walk. A valid
    outline that folds once projected would reach the split the same way,
    so the split is handed a folded figure directly.
    """
    bowtie = Polygon([(10.0, 10.0), (10.2, 10.2), (10.0, 10.2), (10.2, 10.0)])
    assert not bowtie.is_valid

    for call in (itacart.polyfill, itacart.count_internal_cells):
        with pytest.raises(GeometryError, match="not a valid geometry"):
            call(bowtie, 7)

    folded = geometry._project(bowtie)
    assert not folded.is_valid
    with pytest.raises(GeometryError, match="crosses itself"):
        geometry._quadrant_pieces(folded)


# --------------------------------------------------------------------------
# The radius of the boundary screen
# --------------------------------------------------------------------------


def _parcel_east_of_greenwich(degrees: float, metres: float = 60.0) -> Polygon:
    """A small square parcel that many degrees east of the prime meridian."""
    side = metres / 111_319.0
    lat = 51.4776
    return Polygon(
        [
            (degrees, lat),
            (degrees + side, lat),
            (degrees + side, lat + side),
            (degrees, lat + side),
        ]
    )


def test_a_parcel_beside_the_meridian_is_filled_at_every_resolution() -> None:
    """What used to be a narrowing is now an absence.

    The screen ran at resolution 1 and refused a band up to ten
    kilometres wide to protect against a family one cell wide, so this
    parcel was refused coarse and filled fine. The meridian is no longer
    a refused family at all, so the parcel is filled at both, and the
    band that used to stand in front of it is gone rather than smaller.

    The property that band replaced -- a screen measured at the target
    side rather than at resolution 1 -- is still held, on the family
    that still has one, by
    ``test_the_band_of_a_polar_base_cell_is_the_whole_cell``.

    Read under ``intersects``. The parcel is sixty metres across and a
    resolution-3 cell is a kilometre, so under the default mode it holds
    no centre and is refused for its size rather than for its place --
    a different refusal, which would hide the one this test is about.
    """
    parcel = _parcel_east_of_greenwich(0.01)
    for resolution in (3, 7):
        cells = itacart.decompose(
            itacart.polyfill(parcel, resolution, containment="intersects")
        )
        assert cells
        assert all(itacart.get_resolution(c) == resolution for c in cells)


def test_the_absorbing_cell_is_filled_at_every_resolution() -> None:
    """The last lattice column is walked, and answers what names it.

    This is the equality the refusing version of this test said it would
    become. The fill used to refuse the whole column because its cell is
    not the sheared square the descent tests; it now descends the tree
    the hierarchy proves for that cell, whose family runs two to seven
    children rather than ``d`` squared.

    The assertion is against ``geo_to_cell`` rather than against a
    literal: a parcel small enough to sit inside one target cell must be
    answered by the cell that names its own interior. The neighbouring
    ordinary column is filled alongside, so a fill that had stopped
    working everywhere would not pass.
    """
    side = itacart.cell_size(1)
    metres = 60.0 / 111_319.0
    for row in (100, 300, 573):
        last = itacart.last_lattice_column("NE", row, side)
        assert itacart.is_valid_cell(f"NE({last:04d}/{row:04d})")
        ring = itacart.cell_to_boundary(f"NE({last:04d}/{row:04d})")
        longitudes = [x for x, _ in ring]
        latitudes = [y for _, y in ring]
        latitude = (min(latitudes) + max(latitudes)) / 2.0
        middle = (min(longitudes) + max(longitudes)) / 2.0
        absorbing = Polygon(
            [
                (middle, latitude),
                (middle + metres, latitude),
                (middle + metres, latitude + metres),
                (middle, latitude + metres),
            ]
        )
        inner = itacart.cell_to_centroid(f"NE({last - 1:04d}/{row:04d})")
        ordinary = Polygon(
            [
                (inner[0], inner[1]),
                (inner[0] + metres, inner[1]),
                (inner[0] + metres, inner[1] + metres),
                (inner[0], inner[1] + metres),
            ]
        )
        for resolution in (1, 5):
            for parcel in (absorbing, ordinary):
                cells = itacart.decompose(
                    itacart.polyfill(parcel, resolution, containment="intersects")
                )
                assert cells
                assert all(itacart.is_valid_cell(cell) for cell in cells)

            named = itacart.geo_to_cell(
                middle + metres / 2.0, latitude + metres / 2.0, resolution
            )
            covering = itacart.decompose(
                itacart.polyfill(absorbing, resolution, containment="intersects")
            )
            assert named in covering, (row, resolution, named)


def test_the_row_below_the_pole_no_longer_needs_the_escape() -> None:
    """The narrowing was carrying a row that should never have needed it.

    The screen is applied at resolution 1 and is conservative, so a
    refused base cell used to keep an escape: measured at the target side
    instead, most of its descendants are ordinary and can be let through.
    Row 999 lived on that escape, refused outright at resolution 1 and
    filled only from some resolution down.

    It was refused because it was miscounted as the polar row. It is a
    row of ordinary border-absorbing trapezoids, and it fills at every
    resolution now, resolution 1 included, without the escape being
    involved at all.

    The cells are checked against ``geo_to_cell`` rather than only for
    existence: a fill that named cells nothing else would name would be a
    hole, not a repair.
    """
    side = itacart.cell_size(1)
    polar = int(math.floor(MERIDIAN_QUADRANT / side))
    ring = itacart.cell_to_boundary(f"NE(0001/{polar - 1:04d})")
    latitudes = [y for _, y in ring]
    longitudes = [x for x, _ in ring]
    latitude = min(latitudes) + 0.95 * (max(latitudes) - min(latitudes))
    longitude = (min(longitudes) + max(longitudes)) / 2.0
    metres = 60.0 / 111_319.0
    parcel = Polygon(
        [
            (longitude, latitude),
            (longitude + metres, latitude),
            (longitude + metres, latitude + metres),
            (longitude, latitude + metres),
        ]
    )

    for resolution in (1, 5, 7):
        cells = itacart.decompose(
            itacart.polyfill(parcel, resolution, containment="intersects")
        )
        assert cells, resolution
        assert all(itacart.is_valid_cell(cell) for cell in cells)
        assert (
            itacart.geo_to_cell(
                longitude + metres / 2.0, latitude + metres / 2.0, resolution
            )
            in cells
        ), resolution


def test_the_meridian_walk_reaches_the_polar_row_and_stops_there() -> None:
    """The seam walk runs to the pole's own row, and no further.

    Column zero is skipped by ``_base_cells`` on purpose, because it
    names the meridian triangle this walk descends, so the cap is
    reachable through this walk or not at all. It used to stop two rows
    short: the walk asked ``last_lattice_column`` about the row and the
    row above and read an answer of zero as *no cell*, which is one cell
    in an eastern quadrant and none in a western one, and which is zero
    for row 1000 and zero again for 1001.

    It asks the pole now. A region reaching ninety degrees fills every
    row of the seam up to and including the polar one, and the cap is
    among the cells it names.
    """
    region = Polygon([(-0.2, 89.0), (0.2, 89.0), (0.2, 90.0), (-0.2, 90.0)])
    plane, _views = geometry._prepare(region, 1, densify=True, snap=True)
    rows = geometry._meridian_rows(plane)
    assert rows
    side = itacart.cell_size(1)
    polar = int(math.floor(MERIDIAN_QUADRANT / side))
    last = max(row for _hemisphere, row in rows)
    assert last == polar
    assert all(row <= polar for _hemisphere, row in rows), rows

    cells = itacart.decompose(itacart.polyfill(region, 1))
    assert cells
    assert all(itacart.is_valid_cell(cell) for cell in cells)
    assert f"NE({0:0{4}d}/{polar:04d})" in cells, cells


def test_a_geometry_touching_an_extension_zone_is_not_lifted() -> None:
    """Touching the zone is not entering it.

    A polygon whose edge lies exactly on a zone window leaves a
    zero-width piece inside it, and lifting that piece would move a line
    the geometry only borders past the antemeridian while the body of it
    stayed behind. Read on both zones and on both edges the geometry can
    meet -- the longitude limit and the poleward end of the row band --
    because the two are cut by different kinds of line.
    """
    for window in geometry._extension_windows():
        west, south, east, north = window.bounds
        beyond_the_limit = Polygon(
            [
                (east, south + 0.5),
                (east + 0.5, south + 0.5),
                (east + 0.5, north - 0.5),
                (east, north - 0.5),
            ]
        )
        beyond_the_band = Polygon(
            [
                (west + 0.5, north),
                (east - 0.5, north),
                (east - 0.5, north + 0.5),
                (west + 0.5, north + 0.5),
            ]
        )
        for outside in (beyond_the_limit, beyond_the_band):
            assert geometry._lift_extensions(outside).equals(outside)


def test_no_column_past_the_last_ever_escapes_the_screen() -> None:
    """The escape must not reach a column the row does not have.

    Near the end of a row the absorbing cell is wider than the square, so
    a position inside the border strip is given a raw column past the
    end. The band used to stop short of those columns, leaving an empty
    intersection, and an empty band is how the escape reads "the geometry
    does not reach the family" -- so the fill named cells that do not
    exist. Read on the rows where the strip is widest, and two columns
    out, because over Chukotka it is two and not one.
    """
    side = itacart.cell_size(1)
    for quadrant in ("NE", "SE"):
        for row in (100, 300, 573, 732):
            last = itacart.last_lattice_column(quadrant, row, side)
            for column in (last, last + 1, last + 2):
                band = geometry._anomalous_band(quadrant, column, row, side)
                u0 = (column + row) * geometry._L1
                v0 = row * geometry._L1
                square = box(u0, v0, u0 + geometry._L1, v0 + geometry._L1)
                assert band.covers(square), (
                    f"{quadrant} column {column} of row {row} is not covered, "
                    "so a geometry inside it would escape the screen"
                )


def test_the_geometric_screen_never_admits_an_anomalous_cell() -> None:
    """The screen is validated against the index side, by enumeration.

    Complete to level 3, over the three families and the ordinary
    interior in four quadrants. Every node is compared both ways, and
    the geometric answer is never allowed to be the more permissive of
    the two.

    Divergences do exist and are all in the safe direction. Inside the
    meridian family they sit at a negative column -- west of the line,
    outside the quadrant, unreachable because the quadrant clip leaves
    ``u - v`` positive -- and that much is asserted. Inside the
    border-absorbing family they sit at positive columns, because
    absorption gives the parent more children than the square
    subdivision has and the two are not the same tree there. Only
    conservativeness is claimed for that family; equivalence is not.

    Three levels is not thirteen, and the depth is the whole point of
    the screen, so the shallow sweep is the floor and not the argument.
    ``test_the_screen_agrees_all_the_way_down_the_frontier`` carries it
    to resolution 13.
    """
    side = itacart.cell_size(1)
    checked = 0
    skipped = 0
    for quadrant in QUADRANTS:
        for row in (0, 1, 100, 451, 700):
            last = itacart.last_lattice_column(quadrant, row, side)
            for column in sorted({0, 1, 2, last // 2, last - 1, last, last + 1}):
                stack = [
                    (
                        f"{quadrant}({column:04d}/{row:04d})",
                        (column + row) * side,
                        row * side,
                        side,
                        1,
                    )
                ]
                while stack:
                    cell, u0, v0, cell_side, level = stack.pop()
                    # The two sides are compared only where there is
                    # something to compare. The column set below includes
                    # column zero in every quadrant, and west of the prime
                    # meridian that column names no cell, so the index side
                    # has no truth to state there and the plane side would
                    # be answering about the eastern twin. Skipped and
                    # counted rather than silently walked.
                    if not itacart.is_valid_cell(cell):
                        skipped += 1
                        continue
                    refused = _screen_refuses(quadrant, u0, v0, cell_side)
                    anomalous = _is_anomalous(cell)
                    checked += 1
                    if anomalous and not refused:
                        raise AssertionError(f"screen admits {cell}")
                    if refused and not anomalous and column == 0:
                        assert round((u0 - v0) / cell_side) < 0, cell
                    if level < 3:
                        stack.extend(
                            _children_in_lockstep(cell, u0, v0, cell_side, level)
                        )
    # Pinned by equality rather than by a floor. A containment assertion
    # lets the number rot in both directions: a walk that quietly stops
    # reaching cells still passes it. The skipped population is every
    # spelling the walk builds that names no cell -- at resolution 1, ten
    # western zero columns and twenty columns past their row; below it,
    # 164 refined descendants that fall outside the border -- and it is
    # pinned too, so the leniency cannot return by widening what the walk
    # is willing to admit.
    assert (checked, skipped) == (11086, 194)


def test_the_screen_agrees_all_the_way_down_the_frontier() -> None:
    """The same agreement, carried from resolution 1 to resolution 13.

    A complete sweep to resolution 13 is not enumerable -- a single base
    cell holds a million cells per axis down there -- so this walks the
    frontier instead, which is where the derivation says the only
    anomalies can be. At every step *all* children of the current node
    are compared both ways, and the walk then continues into an
    anomalous one. Ordinary nodes are covered exhaustively by the
    shallow sweep; what this adds is depth on the one path where the
    families actually live.

    The shape of the shortcut is deliberate. The refinement is a fixed
    point only where the enumeration is complete, and the enumeration is
    complete at every level of the walk -- it is the breadth across
    levels that is traded away, not the breadth within one.

    This test exists because a bound measured at two levels and asserted
    for thirteen is how a fixed overlap threshold survives unnoticed.
    """
    base = itacart.cell_size(1)
    deepest = 0
    checked = 0
    skipped = 0
    for quadrant in QUADRANTS:
        last = itacart.last_lattice_column(quadrant, 451, base)
        polar = max(
            row
            for row in range(700, 1000)
            if itacart.last_lattice_column(quadrant, row, base) > 0
        )
        for column, row in ((0, 451), (last, 451), (1, polar)):
            cell = f"{quadrant}({column:04d}/{row:04d})"
            # Column zero west of the meridian names no cell, so there is
            # no descent to make and no agreement to check.
            if not itacart.is_valid_cell(cell):
                skipped += 1
                continue
            u0, v0, side, level = (column + row) * base, row * base, base, 1
            while level < 13:
                children = _children_in_lockstep(cell, u0, v0, side, level)
                anomalous_children = []
                for name, cu, cv, cs, step in children:
                    refused = _screen_refuses(quadrant, cu, cv, cs)
                    anomalous = _is_anomalous(name)
                    checked += 1
                    if anomalous and not refused:
                        raise AssertionError(f"screen admits {name}")
                    if anomalous:
                        anomalous_children.append((name, cu, cv, cs, step))
                if not anomalous_children:
                    break
                cell, u0, v0, side, level = anomalous_children[-1]
                deepest = max(deepest, level)
            assert level == 13, (quadrant, column, row, level)
    assert deepest == 13, deepest
    # The two skipped are column zero at row 451 in NW and SW, where the
    # meridian column does not exist. Equality for the same reason as the
    # screen test above.
    assert (checked, skipped) == (1740, 2)


def test_the_two_remaining_families_are_refused_when_actually_touched() -> None:
    """Narrowing the radius must not narrow the refusal itself.

    Two families, not three. The prime-meridian column left the screen
    when it gained a descent of its own, and a parcel sitting on the
    line is now filled; the last lattice column and the polar row are
    refused exactly as before, in all four quadrants.

    Read at resolution 7, where the old version read at 13. The old one
    could afford 13 because the screen refused before counting anything;
    now the parcel is actually filled, and at a one-centimetre cell a
    sixty-metre parcel is thirty-six million of them.
    """
    on_the_line = itacart.decompose(itacart.polyfill(_parcel_east_of_greenwich(0.0), 7))
    assert on_the_line
    assert all(itacart.is_valid_cell(cell) for cell in on_the_line)
    side = itacart.cell_size(1)
    for quadrant in QUADRANTS:
        polar = max(
            row
            for row in range(700, 1000)
            if itacart.last_lattice_column(quadrant, row, side) > 0
        )
        for column, row in (
            (itacart.last_lattice_column(quadrant, 100, side), 100),
            (1, polar),
        ):
            with pytest.raises((NonExistentCellError, DomainError)):
                geometry._check_addressable(quadrant, column, row)


def _screen_refuses(quadrant: str, u0: float, v0: float, side: float) -> bool:
    """The geometric screen, expressed only on the node's own coordinates.

    The lattice column of a node of side ``s`` at ``(u0, v0)`` is
    ``(u0 - v0) / s``, which reduces to ``u_index - row`` at resolution
    1. Nothing here names a cell, which is the point: the fill may not
    reach an index-rendering door, so it cannot ask ``absorbs_border``.
    """
    column = round((u0 - v0) / side)
    row = round(v0 / side)
    last_here = itacart.last_lattice_column(quadrant, row, side)
    beyond = itacart.last_lattice_column(quadrant, row + 1, side)
    return column <= 0 or last_here <= 0 or beyond <= 0 or column >= last_here


def _is_anomalous(cell: str) -> bool:
    """The index-side truth the screen has to agree with."""
    try:
        if itacart.absorbs_border(cell):
            return True
        return itacart.cell_shape(cell) != "parallelogram"
    except Exception:  # pragma: no cover - defensive
        return True


def _children_in_lockstep(
    cell: str, u0: float, v0: float, side: float, level: int
) -> list[tuple[str, float, float, float, int]]:
    """Every child of a node, as index and as square, paired by position."""
    step = level + 1
    divisor = itacart.linear_refinement_ratio(step)
    child = side / divisor
    return [
        (
            hierarchy._descend(cell, geometry._child_code(row, column, step)),
            u0 + column * child,
            v0 + row * child,
            child,
            step,
        )
        for row in range(divisor)
        for column in range(divisor)
    ]


def test_what_is_filled_beside_the_meridian_is_ordinary_and_nominal() -> None:
    """The cells the narrower screen lets through are real cells.

    Inside an anomalous base cell the index descent and the square
    descent are not the same tree, so admitting geometry there is only
    safe if the indices the fill emits name the squares it accepted.
    Measured on the round trip rather than assumed.
    """
    parcel = _parcel_east_of_greenwich(0.001)
    cells = itacart.decompose(itacart.polyfill(parcel, 7))
    assert cells
    nominal = itacart.nominal_cell_area(7)
    for cell in cells:
        assert itacart.is_valid_index(cell), cell
        assert not itacart.absorbs_border(cell), cell
        assert itacart.cell_shape(cell) == "parallelogram", cell
        area = Polygon(boundary.plane_ring(cell)[1]).area
        assert area == pytest.approx(nominal, rel=1e-9), cell
        assert itacart.cell_to_polygon(cell).intersects(parcel), cell


# --------------------------------------------------------------------------
# The grid's own geometry as the query
# --------------------------------------------------------------------------


def _prepared_from_the_grid(cell: str):
    """Plane geometry and lattice views built from the cell itself.

    The fill takes a geodetic polygon, so every measurement of it has to
    hand it one. Building that polygon with ``cell_to_boundary`` sends
    the cell out through longitude and back, and the trip does not
    return the cell: an edge of an ITACaRT cell is a straight line in
    the sinusoidal plane, and densifying it as a geodesic on the
    ellipsoid rebuilds it as a different curve. Measured, the round trip
    returns 99.9269 per cent of a border-absorbing cell, 50 per cent of
    the southern cap and 33.1111 per cent of a cell in row 999.

    So a test that wants to ask about the *descent* has to hand it the
    figure the grid produces, not one that has been through longitude.
    This builds exactly what ``_prepare`` builds -- projected plane,
    then one sheared lattice view per quadrant -- from ``plane_ring``.
    """
    plane = Polygon(boundary.plane_ring(cell)[1])
    views = [
        (quadrant, geometry._LatticeView(geometry._to_lattice(piece, quadrant)))
        for quadrant, piece in geometry._quadrant_pieces(plane)
    ]
    return plane, views


def _real_family(cell: str, resolution: int) -> set[str]:
    """Every descendant of ``cell`` at ``resolution``, from the hierarchy."""
    family = {cell}
    for _level in range(resolution - itacart.get_resolution(cell)):
        family = {child for node in family for child in hierarchy._children_of(node)}
    return family


def test_the_meridian_walk_offers_the_polar_row_in_both_hemispheres() -> None:
    """The cap is column zero, and only this walk can reach it.

    ``_base_cells`` skips column zero deliberately, because that column
    names the meridian triangle this walk exists to descend. So if the
    walk does not offer the polar row, nothing offers the cap and the
    fill returns empty without refusing anything.

    The row is decided by the pole. It used to be decided by a pair of
    ``last_lattice_column(...) <= 0`` tests, one on the row and one on
    the row above, and that answer saturates: zero means one cell in an
    eastern quadrant and none in a western one, and it is zero for both
    row 1000 and row 1001. The pair skipped row 999 on one test and row
    1000 on the other, so the walk stopped two rows short.
    """
    polar = int(math.floor(MERIDIAN_QUADRANT / itacart.cell_size(1)))

    for hemisphere, cap in (("N", "NE(0000/1000)"), ("S", "SE(0000/1000)")):
        plane, _views = _prepared_from_the_grid(cap)
        rows = geometry._meridian_rows(plane)
        assert (hemisphere, polar) in rows, (cap, rows)
        assert all(row <= polar for _hemisphere, row in rows), rows


@pytest.mark.parametrize(
    ("cell", "label"),
    [
        ("NE(0000/1000)", "northern cap"),
        ("SE(0000/1000)", "southern cap"),
        ("NE(0001/0999)", "row below the pole"),
        ("NE(1414/0500)", "border-absorbing, lateral"),
        ("NE(0500/0300)", "interior"),
    ],
)
@pytest.mark.parametrize("resolution", [2, 3])
def test_the_descent_recovers_the_real_tree_from_the_grids_own_geometry(
    monkeypatch: pytest.MonkeyPatch, cell: str, label: str, resolution: int
) -> None:
    """Asked about the cell itself, the fill answers with the cell's children.

    Set against set, not cardinality against cardinality: a walk that
    named the right number of the wrong cells would pass a count and
    fail here, and one did -- before the polar row was routed to the
    proved tree, the uniform meridian walk named two spellings
    ``is_valid_cell`` rejects and missed two it should have named.

    ``count_internal_cells`` is asserted alongside, because the two
    reach the same answer by different code and a disagreement between
    them is a defect in whichever is wrong.
    """
    monkeypatch.setattr(
        geometry,
        "_prepare",
        lambda _geometry, _resolution, densify=True, *, snap: _prepared_from_the_grid(
            cell
        ),
    )
    ignored = Polygon([(0.0, 0.0), (0.1, 0.0), (0.1, 0.1)])

    filled = set(
        itacart.decompose(itacart.polyfill(ignored, resolution, compact=False))
    )
    counted = itacart.count_internal_cells(ignored, resolution)

    assert filled == _real_family(cell, resolution), label
    assert counted == len(filled), label
    assert all(itacart.is_valid_cell(spelling) for spelling in filled), label


def test_the_band_of_a_polar_base_cell_is_the_whole_cell() -> None:
    """At the pole the anomaly is not a strip, and the band says so.

    ``_anomalous_band`` draws the three families as bands one cell side
    wide, which is the point of narrowing the screen. The polar family
    is the one case where that narrowing buys nothing: every target row
    in the last base row of a quadrant is clipped by the pole, so the
    band covers the base cell entirely and the screen refuses exactly
    what it refused before.

    Two rows and three resolutions, because the branch that detects the
    polar reach differs between them. Row 999 is the last row that
    addresses a cell at resolution 1; row 1000 is straddled by the pole,
    so it addresses none, and its topmost target rows lie beyond the
    meridian quadrant.
    """
    base = itacart.cell_size(1)
    for row in (999, 1000):
        for resolution in (1, 3, 5):
            band = geometry._anomalous_band("NE", 1, row, itacart.cell_size(resolution))
            assert not band.is_empty, (row, resolution)
            if row == 1000:
                assert band.area == pytest.approx(base * base, rel=1e-9), (
                    row,
                    resolution,
                )
            else:
                assert band.area > 0.7 * base * base, (row, resolution)


def test_the_polar_reach_is_found_by_search_when_it_starts_mid_cell() -> None:
    """The pole can fall inside a base cell, not only below it.

    Row 1000 at resolution 5 is the case: its lowest target rows still
    address cells, its highest do not, and the first row that does not
    has to be located rather than assumed. A linear walk would visit up
    to a million rows at resolution 13, so it is a search, and the
    search is what this covers.
    """
    base = itacart.cell_size(1)
    side = itacart.cell_size(5)
    lowest = int(math.floor(1000 * base / side))
    highest = int(math.floor((1000 * base + base - side / 2.0) / side))
    assert itacart.last_lattice_column("NE", lowest + 1, side) > 0
    assert itacart.last_lattice_column("NE", highest + 1, side) <= 0
    crossings = [
        row
        for row in range(lowest, highest + 1)
        if (itacart.last_lattice_column("NE", row + 1, side) <= 0)
        != (itacart.last_lattice_column("NE", row, side) <= 0)
    ]
    assert len(crossings) == 1, crossings
    band = geometry._anomalous_band("NE", 1, 1000, side)
    assert band.area == pytest.approx(base * base, rel=1e-9)


def test_the_polar_row_narrows_with_the_target_like_the_others() -> None:
    """The polar family narrows too, and the first version of this lied.

    It asserted that a parcel in the last base row is refused at every
    resolution, on the reasoning that the pole clips the whole row. That
    is true of base row 1000, which the pole straddles, and false of row
    999: most of its target rows address cells perfectly well, and only
    the topmost do not. The screen refuses the top and admits the rest,
    which is the whole intent of measuring the band at the target side.

    So the assertion is the one that survives measurement. Base row 999
    is refused entirely at resolution 1, because at that side the row is
    the polar row and holds one cell. Below that side the band shrinks,
    and it never grows back.
    """
    base = itacart.cell_size(1)
    whole = base * base
    assert geometry._anomalous_band("NE", 1, 999, base).area == pytest.approx(
        whole, rel=1e-9
    )
    areas = [
        geometry._anomalous_band("NE", 1, 999, itacart.cell_size(r)).area
        for r in (1, 3, 5, 7)
    ]
    assert areas[0] == pytest.approx(whole, rel=1e-9)
    assert all(area < whole for area in areas[1:]), areas
    assert max(areas[1:]) < 0.9 * whole, areas
    for row in (999, 1000):
        for resolution in (1, 3, 5):
            band = geometry._anomalous_band("NE", 1, row, itacart.cell_size(resolution))
            assert not band.is_empty, (row, resolution)


# --------------------------------------------------------------------------
# The border-absorbing roots: counting, the ceiling, and what is refused
# --------------------------------------------------------------------------


def _last_column_footprint(row: int) -> Polygon:
    """The drawn outline of the last lattice column of an eastern row."""
    last = itacart.last_lattice_column("NE", row, itacart.cell_size(1))
    return Polygon(itacart.cell_to_boundary(f"NE({last:04d}/{row:04d})"))


#: A band inside the northern cap, clear of its row boundary: the polar row
#: starts near 89.9824 degrees, and both edges sit inside it.
CAP_BAND = box(1.0, 89.985, 89.0, 89.999)


@pytest.mark.parametrize("family, resolution", (("last column", 3), ("polar cap", 5)))
def test_a_count_over_the_absorbing_family_ignores_the_fill_ceiling(
    monkeypatch: pytest.MonkeyPatch, family: str, resolution: int
) -> None:
    """A count names no cell, so the fill's ceiling is not its budget.

    Both absorbing walkers can stop once a running total passes a budget.
    The fill needs that, because it must refuse before it names anything.
    A count materialises nothing and is meant to run where the index could
    not be built, yet it was handed the fill's ceiling as its budget and
    returned the partial total at the moment of passing it. With the
    ceiling lowered to 5 it answered 25 for the 63 cells under the last
    column of row 300, and 6 for the 140 under this band of the cap. The
    real ceiling hides the same truncation until a single absorbing root
    holds more target cells than it.

    The answer is read from the fill, which counts by naming, before the
    ceiling is lowered; the fill itself would be refused under it.
    """
    parcel = _last_column_footprint(300) if family == "last column" else CAP_BAND
    expected = len(itacart.decompose(itacart.polyfill(parcel, resolution)))
    assert itacart.count_internal_cells(parcel, resolution) == expected
    ceiling = 5
    assert expected > ceiling
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", ceiling)
    assert itacart.count_internal_cells(parcel, resolution) == expected


def test_contains_keeps_a_cap_cell_exactly_when_its_ring_lies_inside() -> None:
    """Under ``contains`` an absorbing leaf is judged by its own ring.

    The cap is walked ring by ring, because none of its cells is the
    sheared square the ordinary descent tests. The rule is checked against
    the query the fill itself builds, cell by cell, over every cell the
    band touches: kept exactly when the effective ring lies inside. On this
    band at resolution 5 that keeps 106 of the 192, so the predicate is
    seen answering both ways, and the three modes still nest.
    """
    resolution = 5
    plane, _ = geometry._prepare(CAP_BAND, resolution, densify=True, snap=False)
    kept = set(itacart.decompose(itacart.polyfill(CAP_BAND, resolution, "contains")))
    centre = set(itacart.decompose(itacart.polyfill(CAP_BAND, resolution)))
    touched = set(
        itacart.decompose(itacart.polyfill(CAP_BAND, resolution, "intersects"))
    )
    assert kept <= centre <= touched
    assert kept and touched - kept
    for cell in touched:
        inside = plane.contains(Polygon(boundary.plane_ring(cell)[1]))
        assert (cell in kept) == inside, cell


def test_contains_refuses_an_ordinary_cell_an_absorbing_root_only_crosses() -> None:
    """Below an absorbing root the ordinary cells are squares again.

    Once a branch of the absorbing tree leaves the family, the closed-form
    descent resumes and its cells are tested as sheared squares. The query
    here covers the last column of row 300 from a parallel cut through the
    middle of a resolution-3 row upwards, and that parallel is a straight
    line in the plane, so whether an ordinary cell lies wholly inside is
    read from its southern edge alone. Of the 33 ordinary cells under the
    root that the query touches, 27 are kept and the 6 it crosses are not.
    Absorbing leaves are left out of the reckoning on purpose: whether
    their rings can lie inside a query that stops at the antemeridian is a
    question about their chords, not about this descent.
    """
    root = _last_column_footprint(300)
    last = itacart.last_lattice_column("NE", 300, itacart.cell_size(1))
    lons = [x for x, _ in root.exterior.coords]
    lats = [y for _, y in root.exterior.coords]
    cut = min(lats) + 0.35 * (max(lats) - min(lats))
    query = box(min(lons) - 0.05, cut, 180.0, max(lats) + 0.05)
    under = set(itacart.get_children(f"NE({last:04d}/0300)", 3, flatten=True))
    ordinary = {cell for cell in under if not itacart.absorbs_border(cell)}
    kept = set(itacart.decompose(itacart.polyfill(query, 3, "contains")))
    touched = set(itacart.decompose(itacart.polyfill(query, 3, "intersects")))
    crossed = touched & ordinary
    assert crossed & kept and crossed - kept
    for cell in crossed:
        south = min(lat for _, lat in itacart.cell_to_boundary(cell))
        assert (cell in kept) == (south > cut), cell


def test_an_absorbing_pre_count_stops_once_the_fill_ceiling_is_passed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An oversized fill over an absorbing root is refused mid-walk.

    The absorbing tree has no formula, so the fill pre-counts it by walking
    it, against what the budget has left. The walk stops at the first child
    that takes the running total past that. With the ceiling lowered to 5,
    the root of the last column of row 300 stops at 25 of its 63 cells and
    the fill is refused, and the refusal does not recommend the compaction
    this family refuses without saying so.
    """
    parcel = _last_column_footprint(300)
    whole = itacart.count_internal_cells(parcel, 3)
    real = geometry._count_border_node
    roots: list[tuple[int | None, int]] = []

    def watched(
        prepared: object,
        walk: object,
        cell: str,
        level: int,
        target: int,
        containment: str,
        remaining: int | None,
    ) -> int:
        total = real(prepared, walk, cell, level, target, containment, remaining)
        if level == 1:
            roots.append((remaining, total))
        return total

    monkeypatch.setattr(geometry, "_count_border_node", watched)
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", 5)
    with pytest.raises(GeometryError, match="exceeded 5 cells") as refusal:
        itacart.polyfill(parcel, 3)
    assert "border-absorbing cells, which refuse it" in str(refusal.value)
    passed = [
        (left, total) for left, total in roots if left is not None and total > left
    ]
    assert len(passed) == 1, roots
    assert passed[0][1] < whole


@pytest.mark.parametrize("family", ("last column", "polar cap"))
def test_compact_is_refused_by_name_over_an_absorbing_root(family: str) -> None:
    """``compact=True`` reaches a tree with no uniform fold, and says so.

    The fold rule reads a uniform tree, and neither absorbing walker has
    one: the last lattice column holds two to seven children per node and
    the cap its own staircase. The same fill without compaction succeeds,
    so the refusal is about the option and not about the geometry.
    """
    parcel, resolution = (
        (_last_column_footprint(300), 3) if family == "last column" else (CAP_BAND, 5)
    )
    assert itacart.decompose(itacart.polyfill(parcel, resolution))
    with pytest.raises(GeometryError, match="not available over the border-absorbing"):
        itacart.polyfill(parcel, resolution, compact=True)


def test_the_base_screen_meets_only_the_last_column_refusal() -> None:
    """Every whole-row refusal is filtered out before the screen asks.

    ``_base_cells`` asks ``_check_addressable`` only about a positive
    column no further out than its row's last one, and catches just the
    refusal of that last column. The refusals that concern a whole row --
    past the pole, a parallel shorter than one cell, the polar row -- are
    never asked: enumerated over the four quadrants and every row to one
    past the pole, each row that offers such a column is ordinary below
    its last column and refuses the last column as absorbing, and each
    row that offers none is refused whole. That second half is the control:
    it is what the screen would meet if the filter ever let a row through.
    """
    side = itacart.cell_size(1)
    whole_rows = 0
    for quadrant in QUADRANTS:
        first = max(geometry._first_lattice_column(quadrant), 1)
        for row in range(geometry._POLAR_ROW + 2):
            last = itacart.last_lattice_column(quadrant, row, side)
            if last < first:
                with pytest.raises(DomainError) as refusal:
                    geometry._check_addressable(quadrant, first, row)
                assert not isinstance(refusal.value, NonExistentCellError)
                whole_rows += 1
                continue
            if last - 1 >= first:
                geometry._check_addressable(quadrant, last - 1, row)
            with pytest.raises(NonExistentCellError):
                geometry._check_addressable(quadrant, last, row)
    assert whole_rows == 2 * len(QUADRANTS)


def test_the_meridian_walk_is_offered_no_row_past_the_pole() -> None:
    """The rows offered to the meridian walk stop at the polar row.

    Three facts make that so without a test inside the walk. A latitude
    past 90 degrees is refused before anything is projected. The pole
    projects into the polar row, in both hemispheres. And every eastern row
    up to it addresses the meridian column. The walk still reaches the
    polar row when the geometry does.
    """
    from itacart.geodesy import geodetic_to_sinusoidal

    past_the_pole = box(-0.05, 89.99, 0.05, 90.2)
    for call in (itacart.polyfill, itacart.count_internal_cells):
        with pytest.raises(DomainError, match=r"\[-90, 90\]"):
            call(past_the_pole, 1)
    for latitude in (90.0, -90.0):
        _, y = geodetic_to_sinusoidal(0.0, latitude)
        assert math.floor(abs(y) / geometry._L1) == geometry._POLAR_ROW
    side = itacart.cell_size(1)
    for quadrant in ("NE", "SE"):
        first = geometry._first_lattice_column(quadrant)
        rows = range(geometry._POLAR_ROW + 1)
        assert all(
            itacart.last_lattice_column(quadrant, r, side) >= first for r in rows
        )
    plane, _ = geometry._prepare(CAP_BAND, 5, densify=True, snap=True)
    assert geometry._meridian_rows(plane) == [("N", geometry._POLAR_ROW)]


# --------------------------------------------------------------------------
# An areal input must be valid, and is never repaired
# --------------------------------------------------------------------------

#: Invalid areal inputs. The first two have a part past the antemeridian
#: with no area, which is what used to carry them past the extension screen;
#: the bowtie and the overlapping parts are invalid anywhere.
SPIKE_PAST_THE_LINE = Polygon(
    [
        (179.0, 10.0),
        (179.9, 10.0),
        (179.9, 10.5),
        (180.5, 10.5),
        (179.9, 10.5),
        (179.9, 11.0),
        (179.0, 11.0),
    ]
)
COLLAPSED_PAST_THE_LINE = Polygon([(179.0, 10.0), (181.0, 10.0), (180.5, 10.0)])
BOWTIE = Polygon([(10.0, 10.0), (10.2, 10.2), (10.0, 10.2), (10.2, 10.0)])
OVERLAPPING_PARTS = MultiPolygon(
    [box(10.0, 10.0, 10.2, 10.2), box(10.1, 10.1, 10.3, 10.3)]
)
INVALID_AREAS = {
    "spike past the line": SPIKE_PAST_THE_LINE,
    "collapsed past the line": COLLAPSED_PAST_THE_LINE,
    "bowtie": BOWTIE,
    "overlapping parts": OVERLAPPING_PARTS,
}


@pytest.mark.parametrize(
    "valid",
    (
        box(10.0, 10.0, 10.2, 10.2),
        MultiPolygon([box(10.0, 10.0, 10.1, 10.1), box(10.2, 10.2, 10.3, 10.3)]),
    ),
    ids=("polygon", "multipolygon"),
)
def test_a_valid_areal_geometry_is_still_filled(valid: Polygon | MultiPolygon) -> None:
    """The guard on validity asks nothing of a geometry that is valid."""
    assert valid.is_valid
    for containment in ("center", "intersects", "contains"):
        cells = itacart.decompose(itacart.polyfill(valid, 5, containment))
        assert cells and all(itacart.is_valid_cell(cell) for cell in cells)
    centre = itacart.decompose(itacart.polyfill(valid, 5))
    assert itacart.count_internal_cells(valid, 5) == len(centre)


@pytest.mark.parametrize("name", sorted(INVALID_AREAS))
def test_an_invalid_areal_geometry_is_refused_by_its_reason_and_not_repaired(
    name: str,
) -> None:
    """Validity belongs to the geometry the caller wrote, not to where it lies.

    A polygon or multipolygon that is not valid is refused before anything
    reads it -- before the antemeridian screens, the clip or the fill --
    and the refusal carries the engine's own account of what is wrong,
    with where. Nothing is repaired: a zero-width buffer or ``make_valid``
    would decide which area was meant, and the two halves of a bowtie are
    not the caller's to lose. Both entry points refuse alike, under every
    containment mode, because the fill and the count share the entrance.
    """
    from shapely.validation import explain_validity

    invalid = INVALID_AREAS[name]
    assert not invalid.is_valid
    reason = explain_validity(invalid)
    for containment in ("center", "intersects", "contains"):
        with pytest.raises(GeometryError, match="not a valid geometry") as refusal:
            itacart.polyfill(invalid, 3, containment)
        assert type(refusal.value) is GeometryError
        assert reason in str(refusal.value)
    with pytest.raises(GeometryError, match="not a valid geometry"):
        itacart.count_internal_cells(invalid, 3)


@pytest.mark.parametrize("name", ("spike past the line", "collapsed past the line"))
def test_an_invalid_outline_past_the_antemeridian_dies_at_the_validity_guard(
    name: str,
) -> None:
    """The outline that slipped past the extension screen now stops earlier.

    Its part past the line has no area, and the extension screen let such
    a part through as if it were an areal geometry touching the far edge of
    a zone. Filled under ``intersects`` it named the cells of its valid
    remainder, while the same line as a ``LineString`` was refused. Now it
    is refused for what it is, and not for where it reaches: the refusal is
    the validity one, not the antemeridian's.
    """
    invalid = INVALID_AREAS[name]
    past = invalid.difference(box(-180.0, -90.0, 180.0, 90.0))
    assert not past.is_empty and past.area == 0.0
    with pytest.raises(GeometryError, match="not a valid geometry") as refusal:
        itacart.polyfill(invalid, 1, containment="intersects")
    assert not isinstance(refusal.value, AntemeridianError)


def test_the_valid_twin_ending_at_the_antemeridian_is_filled_and_its_rule_kept() -> (
    None
):
    """Valid geometry meets the same antemeridian rules as before.

    A valid box that stops at 180 degrees is filled under every mode, and
    counted as the centre fill is. The same box carried half a degree past
    the line, outside any extension zone, is still refused by the
    antemeridian's own name rather than by the validity guard.
    """
    twin = box(179.0, 10.0, 180.0, 11.0)
    past = box(179.0, 10.0, 180.5, 11.0)
    assert twin.is_valid and past.is_valid
    for containment in ("center", "intersects", "contains"):
        cells = itacart.decompose(itacart.polyfill(twin, 3, containment))
        assert cells and all(itacart.is_valid_cell(cell) for cell in cells)
    centre = itacart.decompose(itacart.polyfill(twin, 3))
    assert itacart.count_internal_cells(twin, 3) == len(centre)
    assert itacart.geo_to_cell(179.5, 10.5, 3) in centre
    for call in (itacart.polyfill, itacart.count_internal_cells):
        with pytest.raises(AntemeridianError):
            call(past, 3)


def test_a_line_keeps_its_own_contract_whatever_the_areal_guard_says() -> None:
    """Having no area is not being invalid.

    A line that crosses itself is a valid geometry, and under ``intersects``
    it is filled as before. Under an areal predicate it keeps no cell and
    says so. The line matching the invalid outline past 180 degrees is
    refused by the antemeridian's name, exactly as it was.
    """
    crossing = LineString([(10.0, 10.0), (10.2, 10.2), (10.0, 10.2), (10.2, 10.0)])
    assert crossing.is_valid and not crossing.is_simple
    assert itacart.decompose(itacart.polyfill(crossing, 3, "intersects"))
    with pytest.raises(GeometryError, match="has no area"):
        itacart.polyfill(crossing, 3)
    with pytest.raises(AntemeridianError):
        itacart.polyfill(LineString([(179.0, 10.0), (181.0, 10.0)]), 1, "intersects")


#: Valid representations the package supports on either side of the screens
#: the validity guard now precedes: the seam, the antemeridian cut, the pole
#: and both extension zones. Every one is valid as written.
VALID_BEFORE_THE_SCREENS = {
    "polygon touching 180": box(179.0, 10.0, 180.0, 11.0),
    "polygon touching -180": box(-180.0, 10.0, -179.0, 11.0),
    "multipolygon cut at the antimeridian": MultiPolygon(
        [box(179.5, 10.0, 180.0, 11.0), box(-180.0, 10.0, -179.5, 11.0)]
    ),
    "rectangle up to the pole": box(-180.0, 89.9, 180.0, 90.0),
    "band inside the northern cap": CAP_BAND,
    "fiji past 180": Polygon(
        [(179.9, -18.0), (181.0, -18.0), (181.0, -17.8), (179.9, -17.8)]
    ),
    "chukotka past 180": Polygon(
        [(179.9, 68.0), (182.0, 68.0), (182.0, 68.2), (179.9, 68.2)]
    ),
    "parts touching at a vertex on 180": MultiPolygon(
        [box(179.5, -18.5, 180.0, -18.0), box(180.0, -18.0, 180.5, -17.5)]
    ),
    "triangle with a vertex on 180": Polygon(
        [(179.0, 10.0), (180.0, 10.5), (179.0, 11.0)]
    ),
}


@pytest.mark.parametrize("name", sorted(VALID_BEFORE_THE_SCREENS))
def test_a_valid_geographic_representation_passes_the_validity_guard(
    name: str,
) -> None:
    """The guard's place, before every screen, refuses no supported shape.

    Validity is asked of the geometry as written, before the antemeridian
    and extension screens and before anything is projected, so a shape
    that only becomes legitimate once those screens have read it would be
    refused there. None of the ways the package lets a caller write the
    seam, the cut at the antimeridian, the pole or an extension zone is
    such a shape: each is valid as written, passes the guard, and is filled.
    """
    written = VALID_BEFORE_THE_SCREENS[name]
    assert written.is_valid
    geometry._refuse_invalid_area(written)
    cells = itacart.decompose(itacart.polyfill(written, 1, "intersects"))
    assert cells and all(itacart.is_valid_cell(cell) for cell in cells)


@pytest.mark.parametrize(
    "cell", ("NE(0000/1000)", "SE(0000/1000)", "SE(1911/0199(2(A2)))")
)
def test_an_outline_the_package_exports_passes_the_validity_guard(cell: str) -> None:
    """What the exporter writes, the fill accepts as written.

    Both caps and a cell of the Fiji zone, whose exported outline lies past
    180 degrees, come back through ``from_geojson`` without the guard
    objecting to any feature.
    """
    from shapely.geometry import shape

    assert itacart.is_valid_cell(cell)
    collection = itacart.cells_to_geojson(cell)
    for feature in collection["features"]:
        outline = shape(feature["geometry"])
        assert outline.is_valid, cell
        geometry._refuse_invalid_area(outline)
    assert itacart.from_geojson(collection, 3, containment="intersects")
