"""The parent-to-child contract on cells that absorb the domain border.

The contract is one sentence: the effective rings of the children cover
the effective ring of the parent, and their interiors do not overlap.
``boundary.plane_ring`` is the authority for "effective" -- it answers
after the pole and after absorption -- and every claim here is measured
against it rather than against the nominal lattice.

Why the nominal lattice is not enough is itself measured, in
:func:`test_the_nominal_footprint_does_not_cover_an_absorbing_parent`.
On the polar cap it covers 78.77% of the parent, and the shortfall does
not shrink with depth, because the absorbed surface is annexed at one
level and dropped by the next descent. Anything that enumerates children
by descending the footprint loses that surface permanently.

Three parents carry the regression. ``A5`` and ``E1`` are the two
trapezoids of the cap's first sub-row, mirror images of one another.
``B2`` is the polar triangle, whose effective ring is 2.93 times its
nominal one and which is therefore the hardest case in the grid.
"""

from __future__ import annotations

import statistics
from collections.abc import Iterator
from functools import lru_cache

import pytest
from shapely.geometry import Polygon

import itacart
from itacart import boundary
from itacart import hierarchy as hy
from itacart.constants import CELL_SIZE_M, QUADRANTS, refinement_alphabet
from itacart.exceptions import DomainError, GeometryError, NonExistentCellError
from itacart.index import join_components, split_components
from itacart.resolutions import refinement_ratio

POLAR_CAP = "NE(0000/1000)"
SOUTH_CAP = "SE(0000/1000)"
CAP_SUBDIVISION = "NE(0000/1000(1))"
EAST_TRAPEZOID = "NE(0000/1000(1(A5)))"
WEST_TRAPEZOID = "NE(0000/1000(1(E1)))"
POLAR_TRIANGLE = "NE(0000/1000(1(B2)))"

ABSORBING_PARENTS = [
    POLAR_CAP,
    CAP_SUBDIVISION,
    EAST_TRAPEZOID,
    WEST_TRAPEZOID,
    POLAR_TRIANGLE,
    SOUTH_CAP,
    "NE(1414/0500)",
    "NE(0311/0900)",
    "SW(1845/0237)",
    "NE(0817/0747)",
]

# Both match the tolerances the enumeration itself proves against, so a
# test that passes here is not passing on a slacker rule than the code
# applies.
#
# SHORTFALL_TOLERANCE is one-sided, and the name says so because the
# reading it replaced did not. It bounds the share of a parent that its
# children fail to reach, and it bounds nothing else: it does not say the
# cover is exact, and it does not say the summed child area is any
# particular multiple of the parent's. Both of those are false in this
# family, because a child's chord reaches past its parent's -- the
# subject of TestContainmentIsNotYetAProperty. Measured, the shortfall it
# bounds runs to 1.1e-12 over the 572 lateral parents and 2.1e-16 under
# the caps, six to fourteen orders of magnitude inside the bound; the
# excess it ignores runs to 2.4 per cent. Widening this into a two-sided area
# tolerance would let the excess pass as slack and is not what it is for.
SHORTFALL_TOLERANCE = 1e-6
OVERLAP_TOLERANCE = 1e-6


def _body(cell: str) -> Polygon:
    """The effective ring of ``cell`` as a polygon."""
    _, ring = boundary.plane_ring(cell)
    assert ring, f"{cell} has no effective ring"
    body = Polygon(ring)
    assert body.is_valid, f"{cell} has a self-intersecting effective ring"
    return body


def _descendants(cell: str, resolution: int) -> set[str]:
    """Every descendant of ``cell`` at ``resolution``, by refinement."""
    family = {cell}
    for _level in range(resolution - itacart.get_resolution(cell)):
        family = {child for node in family for child in hy._children_of(node)}
    return family


def _fill_from_the_grid(cell: str, resolution: int, containment: str) -> set[str]:
    """Fill ``cell`` from its own plane ring, through the real screen.

    The fill takes a geodetic polygon, and sending the cell out through
    longitude and back does not return the cell. So the prepared geometry
    is built from ``plane_ring`` directly, which is what the descent tests
    in the geometry module do for the same reason.
    """
    from itacart import geometry as gm

    def prepared(
        _geometry: object, _resolution: int, densify: bool = True, *, snap: bool
    ) -> object:
        plane = Polygon(boundary.plane_ring(cell)[1])
        views = [
            (quadrant, gm._LatticeView(gm._to_lattice(piece, quadrant)))
            for quadrant, piece in gm._quadrant_pieces(plane)
        ]
        return plane, views

    original = gm._prepare
    gm._prepare = prepared  # type: ignore[assignment]
    try:
        ignored = Polygon([(0.0, 0.0), (0.1, 0.0), (0.1, 0.1)])
        filled = itacart.polyfill(
            ignored, resolution, containment=containment, compact=False
        )
        return set(itacart.decompose(filled)) if filled else set()
    finally:
        gm._prepare = original  # type: ignore[assignment]


def _outside_share(child: Polygon, parent: Polygon) -> float:
    """How much of ``child`` lies outside ``parent``, as a share of itself.

    By areas rather than by ``difference``. Sibling and parent rings share
    edges exactly, and the difference of such rings is where the geometry
    engine reports a side location conflict; the areas ask the same
    question without building the shape.
    """
    return (child.area - child.intersection(parent).area) / child.area


def _every_absorbing_cell_of(quadrant: str, step: int) -> list[str]:
    """The last-column cell of every ``step``-th row of ``quadrant``."""
    side = CELL_SIZE_M[1]
    assert side is not None
    family = []
    for row in range(0, 1000, step):
        column = itacart.last_lattice_column(quadrant, row, side)
        if column < 0:
            continue
        cell = f"{quadrant}({column:04d}/{row:04d})"
        if boundary.is_valid_cell(cell):
            family.append(cell)
    return family


class TestTheContract:
    """The five properties that hold, asserted over every named parent."""

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_every_child_is_a_cell(self, parent: str) -> None:
        children = hy._children_of(parent)
        assert children
        for child in children:
            assert boundary.is_valid_cell(child), child

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_sibling_interiors_do_not_overlap(self, parent: str) -> None:
        rings = [_body(child) for child in hy._children_of(parent)]
        shared = sum(
            first.intersection(second).area
            for index, first in enumerate(rings)
            for second in rings[index + 1 :]
        )
        assert shared <= OVERLAP_TOLERANCE * _body(parent).area

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_the_children_cover_the_whole_effective_ring(self, parent: str) -> None:
        """The half of the contract the nominal descent cannot meet."""
        body = _body(parent)
        covered = sum(
            _body(child).intersection(body).area for child in hy._children_of(parent)
        )
        assert body.area - covered <= SHORTFALL_TOLERANCE * body.area

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_every_child_resolves_back_to_this_parent(self, parent: str) -> None:
        """The relation is a function both ways, spelling notwithstanding.

        A child of an absorbing cell can be spelled under a neighbour of
        its parent, by a column of the resolution-1 lattice or by a
        sibling code at an intermediate level. Neither spelling changes
        who the parent is, and this is what says so.
        """
        for child in hy._children_of(parent):
            assert hy._parent_cell(child) == parent, child

    def test_an_ordinary_cell_still_refines_into_exactly_the_alphabet(self) -> None:
        """Property seven: nothing outside the border family moved."""
        for cell in (
            "NE(0100/0100)",
            "SW(0500/0250)",
            "NE(0000/0500)",
            "SE(0007/0007)",
        ):
            for level in (2, 3):
                parent = cell if level == 2 else hy._children_of(cell)[0]
                if boundary.absorbs_border(parent):
                    continue
                children = hy._children_of(parent)
                assert len(children) == refinement_ratio(level), parent
                assert children == [
                    hy._descend(parent, code) for code in refinement_alphabet(level)
                ], parent


class TestThePolarCap:
    """The cap, where the nominal footprint fails hardest."""

    def test_the_nominal_footprint_does_not_cover_an_absorbing_parent(self) -> None:
        """Why children cannot be the descendants of the nominal lattice.

        Measured on the cap's first subdivision: the nominal rings of its
        own children cover 78.77% of it. The remaining fifth is surface
        the parent holds by absorption, and a descent that starts from
        the footprint never proposes a cell for it. The shortfall is a
        property of the footprint, not of depth, so it does not shrink
        one level down.
        """
        body = _body(CAP_SUBDIVISION)
        nominal = [
            Polygon(boundary._nominal_ring(child)[1])
            for child in hy._children_of(CAP_SUBDIVISION)
        ]
        covered = sum(ring.intersection(body).area for ring in nominal if ring.is_valid)
        assert covered / body.area < 0.79
        assert covered / body.area > 0.78

    def test_the_polar_triangle_is_the_hardest_parent_and_still_closes(self) -> None:
        """``B2``: effective ring 2.93 times the nominal one, tiled exactly.

        Its children are spelled under three different intermediate
        codes, two of which name no cell of their own. Four further
        spellings under that parent carry folded rings and are not
        children: the tiling closes without them, which is what says
        that a folded ring is an artefact of absorption rather than a
        cell.
        """
        body = _body(POLAR_TRIANGLE)
        nominal = Polygon(boundary._nominal_ring(POLAR_TRIANGLE)[1])
        assert body.area / nominal.area > 2.9

        children = hy._children_of(POLAR_TRIANGLE)
        assert len(children) == 8
        prefixes = {child[: child.rindex("(")] for child in children}
        assert len(prefixes) == 3

        covered = sum(_body(child).intersection(body).area for child in children)
        assert body.area - covered <= SHORTFALL_TOLERANCE * body.area

    def test_the_cap_enumerates_symmetrically_east_and_west(self) -> None:
        """Property six, asserted on extents rather than on names.

        The two sides cannot be checked by swapping a quadrant letter or
        transposing a code, because a child's spelling can change stem
        during refinement. What is symmetric is the geometry, so that is
        what is compared: the outer reach of the western children is the
        outer reach of the eastern ones, mirrored.
        """
        for parent in (CAP_SUBDIVISION, POLAR_TRIANGLE):
            west, east = [], []
            for child in hy._children_of(parent):
                xs = [x for x, _ in boundary.plane_ring(child)[1]]
                if max(xs) <= 1e-9:
                    west.append(round(-min(xs), 6))
                elif min(xs) >= -1e-9:
                    east.append(round(max(xs), 6))
            assert west, parent
            assert sorted(west) == sorted(east), parent

    def test_the_two_trapezoids_of_the_first_sub_row_are_mirrors(self) -> None:
        """``A5`` and ``E1``, the regression pair, closing identically."""
        counts, areas = [], []
        for parent in (EAST_TRAPEZOID, WEST_TRAPEZOID):
            body = _body(parent)
            children = hy._children_of(parent)
            covered = sum(_body(child).intersection(body).area for child in children)
            assert body.area - covered <= SHORTFALL_TOLERANCE * body.area
            counts.append(len(children))
            areas.append(round(body.area, 3))
        assert counts[0] == counts[1]
        assert areas[0] == areas[1]


class TestTheFamilyAsAWhole:
    """Totality instruments: the family, not an example from it."""

    @pytest.mark.slow
    def test_the_contract_holds_over_every_seventh_row_of_every_quadrant(self) -> None:
        """Coverage and disjointness across the whole antemeridian family.

        Every seventh row rather than every row, because the property is
        being asserted over a family and the cost is polygonal. The
        stride is coprime with nothing in the grid's periods, so it does
        not sample one phase of the border.
        """
        checked = 0
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 7):
                body = _body(parent)
                rings = [_body(child) for child in hy._children_of(parent)]
                covered = sum(ring.intersection(body).area for ring in rings)
                shared = sum(
                    first.intersection(second).area
                    for index, first in enumerate(rings)
                    for second in rings[index + 1 :]
                )
                assert body.area - covered <= SHORTFALL_TOLERANCE * body.area, parent
                assert shared <= OVERLAP_TOLERANCE * body.area, parent
                checked += 1
        assert checked > 500

    @pytest.mark.slow
    def test_the_child_count_of_the_family_runs_from_two_to_seven(self) -> None:
        """The bound, re-measured after the enumeration was corrected.

        Seven is not an instance. Over the four quadrants it occurs in 45
        cells of 4 000, and the old bound of six was measured against an
        enumeration that reached one column east and no further -- which
        missed children spelled two columns out.
        """
        counts: dict[int, int] = {}
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 1):
                count = len(hy._children_of(parent))
                counts[count] = counts.get(count, 0) + 1
        assert min(counts) == 2
        assert max(counts) == 7
        assert counts[7] == 45
        assert sum(counts.values()) == 4000


class TestContainmentIsNotYetAProperty:
    """Property two, which does not hold, measured rather than assumed.

    A child of a border-absorbing cell reaches past its parent. Not by an
    epsilon, and not by accident: re-measured over every seventh row of
    every quadrant, all 572 parents are affected and exactly two children
    of each leave, which is the two trapezoids. The share of a child's own
    area that falls outside runs to a median of 0.0381%, a mean of
    0.0569%, a ninetieth percentile of 0.1041% and a worst case of 2.4110%
    at ``SE(1930/0196)``.

    The cause is a chord, not a staircase. ``_x_border`` does not depend
    on the lattice row at all -- measuring that is what killed the first
    explanation. The border is a smooth convex curve on the plane, and an
    absorbing cell replaces its outer side by the chord of that curve
    across its own height. A parent spans one cell height and takes one
    chord; its children span half that each and take two, and a chord of a
    convex curve over a shorter span lies further out.

    The model is exact, which is what
    :func:`test_the_excess_is_the_chord_model_to_four_places` asserts: the
    excess area equals the sagitta times half the parent's height, to four
    decimal places, over parents whose sagitta runs from 0.95 to 6.12
    metres.

    None of this is the enumeration's doing. The child rings are what
    ``plane_ring`` returns and the child set is the one the contract
    proves. This is pinned, not asserted away: when the absorbing side
    stops being a per-cell chord, these numbers go to zero and this test
    is the one that should fail.
    """

    def test_a_child_reaches_past_its_parent_on_the_lateral_border(self) -> None:
        for parent in ("NE(1414/0500)", "SW(1845/0237)"):
            body = _body(parent)
            shares = [
                _outside_share(_body(child), body) for child in hy._children_of(parent)
            ]
            assert max(shares) > 1e-6, parent

    @pytest.mark.slow
    def test_it_happens_across_the_family_and_stays_under_three_percent(self) -> None:
        parents = 0
        reaching = 0
        shortfall = 0.0
        leaving: list[float] = []
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 7):
                body = _body(parent)
                shares = [
                    _outside_share(_body(child), body)
                    for child in hy._children_of(parent)
                ]
                parents += 1
                if max(shares) > 1e-6:
                    reaching += 1
                leaving.extend(share for share in shares if share > 1e-6)
                covered = sum(
                    _body(child).intersection(body).area
                    for child in hy._children_of(parent)
                )
                shortfall = max(shortfall, (body.area - covered) / body.area)
        # The population, declared beside the numbers that rest on it. The
        # docstring quoted all four and only the bound was asserted, so a
        # family that started losing a third child, or a stride that stopped
        # reaching a quadrant, would have left the prose standing.
        assert parents == 572
        assert reaching == parents
        assert len(leaving) == 2 * parents == 1_144

        leaving.sort()
        assert statistics.median(leaving) == pytest.approx(0.00038142, rel=1e-3)
        assert statistics.fmean(leaving) == pytest.approx(0.00056931, rel=1e-3)
        assert leaving[int(0.9 * len(leaving))] == pytest.approx(0.00104096, rel=1e-3)
        assert leaving[-1] == pytest.approx(0.02411015, rel=1e-3)
        assert leaving[-1] < 0.03

        # The other side of the same family, measured here so the excess
        # above and the shortfall the criterion bounds are read off the
        # same population. The reading belongs to that population: over a
        # coarser stride of 160 parents it is 2.4e-13, and over these 572
        # it is four times that. The bound is SHORTFALL_TOLERANCE, and
        # both readings sit six orders of magnitude inside it.
        assert shortfall < 1.2e-12
        assert shortfall <= SHORTFALL_TOLERANCE

    def test_the_excess_is_the_chord_model_to_four_places(self) -> None:
        """The closed form, which is what a fix has to remove.

        One chord over a full cell height against two over half heights
        encloses two triangles, each of base half the height and height
        the sagitta at the midpoint. Their total is the sagitta times half
        the parent's height, and that is the measured excess exactly.
        """
        for quadrant, row in (("SE", 196), ("NE", 500), ("NE", 100), ("SW", 237)):
            side = CELL_SIZE_M[1]
            assert side is not None
            column = itacart.last_lattice_column(quadrant, row, side)
            parent = f"{quadrant}({column:04d}/{row:04d})"
            if not boundary.is_valid_cell(parent):
                continue
            whole = _body(parent)
            measured = sum(
                _body(child).area - _body(child).intersection(whole).area
                for child in hy._children_of(parent)
            )
            ys = [y for _, y in boundary.plane_ring(parent)[1]]
            low, high = min(ys), max(ys)
            edges = boundary._x_border
            sagitta = abs(
                edges(quadrant, row, (low + high) / 2.0)
                - (edges(quadrant, row, low) + edges(quadrant, row, high)) / 2.0
            )
            predicted = sagitta * (high - low) / 2.0
            assert sagitta > 0.9, parent
            assert measured == pytest.approx(predicted, rel=1e-4), parent

    def test_the_excess_is_quantized_to_one_cell_and_to_the_same_one_twice(
        self,
    ) -> None:
        """Where a child leaves its parent, one cell still answers, always.

        The population is the excess itself: the three parents this class
        already names, the two children of each that leave, and one
        position inside each sliver. It is deliberately not a sweep --
        the question is whether a position that no parent ring holds has
        a single settled owner, and six positions in three quadrants
        answer it.

        Two answers are asked for. At the child's own resolution the
        owner is the child whose sliver it is, including the one spelled
        under a stem a column east of its parent's. At the parent's
        resolution the owner is the parent, although the position lies
        outside the parent's effective ring: the absorbing cell claims
        the strip its chord cuts off, and that claim is what makes the
        choice total rather than merely consistent.

        The second half is process state. A canonical choice could follow
        from a cache warmed in one order or a set iterated in another, so
        a fresh interpreter walks the same positions in reverse and has
        to return the same owners. That is the form the cap quantizer's
        determinism is already pinned in, reused here for the lateral
        family rather than re-derived.
        """
        import json
        import os
        import subprocess
        import sys

        positions: list[tuple[float, float, int, str, str]] = []
        for parent in ("SE(1930/0196)", "NE(1414/0500)", "SW(1845/0237)"):
            whole = _body(parent)
            for child in hy._children_of(parent):
                excess = _body(child).difference(whole)
                if excess.area <= 1e-6:
                    continue
                spot = excess.representative_point()
                lon, lat = boundary.to_geodetic(spot.x, spot.y)
                positions.append((lon, lat, 2, child, parent))

        assert len(positions) == 6, "two children of each parent leave"

        here: list[tuple[str, str]] = []
        for lon, lat, resolution, child, parent in positions:
            fine = itacart.geo_to_cell(lon, lat, resolution)
            coarse = itacart.geo_to_cell(lon, lat, 1)
            assert itacart.is_valid_cell(fine), (lon, lat)
            assert itacart.is_valid_cell(coarse), (lon, lat)
            assert fine == child
            assert coarse == parent
            assert not _body(parent).contains(Polygon(_body(child).exterior))
            here.append((fine, coarse))

        source = os.path.dirname(os.path.dirname(itacart.__file__))
        script = (
            "import json, sys\n"
            "import itacart\n"
            "positions = json.load(sys.stdin)\n"
            "answers = [(itacart.geo_to_cell(lon, lat, fine),\n"
            "            itacart.geo_to_cell(lon, lat, 1))\n"
            "           for lon, lat, fine in reversed(positions)]\n"
            "json.dump({'file': itacart.__file__, 'answers': answers[::-1]},\n"
            "          sys.stdout)\n"
        )
        path = os.pathsep.join(filter(None, (source, os.environ.get("PYTHONPATH"))))
        run = subprocess.run(
            [sys.executable, "-c", script],
            input=json.dumps([position[:3] for position in positions]),
            capture_output=True,
            text=True,
            check=True,
            env=dict(os.environ, PYTHONPATH=path),
        )
        there = json.loads(run.stdout)
        assert there["file"] == itacart.__file__
        assert [tuple(pair) for pair in there["answers"]] == here

    def test_the_measured_shares_stay_where_they_were_re_measured(self) -> None:
        """The four numbers the docstring quotes, over a smaller stride."""
        shares = []
        for quadrant in ("NE", "SE"):
            for parent in _every_absorbing_cell_of(quadrant, 29):
                whole = _body(parent)
                shares.extend(
                    _outside_share(_body(child), whole)
                    for child in hy._children_of(parent)
                    if _outside_share(_body(child), whole) > 1e-9
                )
        assert shares
        shares.sort()
        median = shares[len(shares) // 2]
        assert 0.0003 < median < 0.0005
        assert max(shares) < 0.03

    def test_a_query_bounded_by_the_border_drops_leaves_under_contains(
        self,
    ) -> None:
        """The same approximation, seen from the query side.

        ``contains`` tests the effective chordal ring, so a leaf whose ring
        protrudes past the piece the query was clipped to is not kept, even
        though the leaf is a descendant of the cell the query is. This is
        not a second defect and it is not a defect of ``contains``: it is
        the excess above, asked about from outside instead of from the
        parent.

        The query here is the parent's own effective ring, which is the
        figure a caller gets when the region of interest stops at the
        domain border. Row 300 in four quadrants, at resolutions 2 and 3,
        with the count of leaves declared beside every figure so a green
        line cannot stand on an empty one.

        Under ``intersects`` the whole family comes back, which is the
        control: it separates "the walk cannot reach these" from "this
        acceptance rule declines them".
        """
        for quadrant in QUADRANTS:
            side = CELL_SIZE_M[1]
            assert side is not None
            column = itacart.last_lattice_column(quadrant, 300, side)
            parent = f"{quadrant}({column:04d}/{300:04d})"
            assert boundary.is_valid_cell(parent) and boundary.absorbs_border(parent)

            query = _body(parent)
            for resolution, total, kept_count in ((2, 3, 1), (3, 63, 53)):
                family = _descendants(parent, resolution)
                assert len(family) == total, (parent, resolution)

                kept = _fill_from_the_grid(parent, resolution, "contains")
                assert len(kept) == kept_count, (parent, resolution)
                assert kept <= family, (parent, resolution)

                assert _fill_from_the_grid(parent, resolution, "intersects") == family

                dropped = family - kept
                protrusions = sorted(
                    _body(leaf).area - _body(leaf).intersection(query).area
                    for leaf in dropped
                )
                assert len(protrusions) == total - kept_count
                assert min(protrusions) > 0.0, (parent, resolution)

    def test_the_polar_cap_does_not_have_this_defect(self) -> None:
        """The cap is clipped by the pole, not by the staircase border."""
        for parent in (CAP_SUBDIVISION, POLAR_TRIANGLE, EAST_TRAPEZOID, WEST_TRAPEZOID):
            body = _body(parent)
            for child in hy._children_of(parent):
                assert _outside_share(_body(child), body) < 1e-6, child
            covered = sum(
                _body(child).intersection(body).area
                for child in hy._children_of(parent)
            )
            assert (body.area - covered) / body.area < 2.2e-16, parent


class TestFoldedSpellingsAreNotCells:
    """The regression for what a folded ring is, and is not.

    Absorbing a border steeper than the lattice shear can fold the
    quadrilateral the absorption builds. The shoelace area of such a ring
    is the difference of its two lobes rather than nothing, so an area
    test accepts it, and the existence predicate used to.

    Nothing repairs one. The pair below reports 104 724.4 square metres
    from rings that cover 322 621.8 once rebuilt, and the surface they
    seem to claim already belongs to the trapezoid that absorbed it.
    """

    FOLDED = ("NE(0000/1000(3(A1)))", "NE(0000/1000(2(A1)))")
    ABSORBERS = (WEST_TRAPEZOID, EAST_TRAPEZOID)

    @pytest.mark.parametrize("cell", FOLDED)
    def test_a_folded_ring_does_not_name_a_cell(self, cell: str) -> None:
        assert not boundary.is_valid_cell(cell)

    @pytest.mark.parametrize("cell", FOLDED)
    def test_the_area_alone_would_have_accepted_it(self, cell: str) -> None:
        """Which is why the predicate cannot rest on area."""
        ring = boundary.plane_ring(cell)[1]
        assert boundary.ring_area(ring) > 100_000.0
        assert not Polygon(ring).is_valid

    @pytest.mark.parametrize("absorber", ABSORBERS)
    def test_the_absorbing_twin_is_a_cell_and_holds_that_surface(
        self, absorber: str
    ) -> None:
        assert boundary.is_valid_cell(absorber)
        effective = _body(absorber)
        nominal = Polygon(boundary._nominal_ring(absorber)[1])
        assert effective.is_valid
        assert round(effective.area - nominal.area, 1) == 104_724.4

    def test_no_child_of_an_absorbing_parent_carries_a_folded_ring(self) -> None:
        for parent in ABSORBING_PARENTS:
            for child in hy._children_of(parent):
                assert _body(child).is_valid, child

    def test_a_folded_spelling_is_never_proposed_as_a_child(self) -> None:
        """``B2`` closes exactly, and none of the four folds is in the set."""
        folded = {
            "NE(0000/1000(1(D2(1))))",
            "NE(0000/1000(1(D2(2))))",
            "NE(0000/1000(1(B4(1))))",
            "NE(0000/1000(1(B4(2))))",
        }
        children = set(hy._children_of(POLAR_TRIANGLE))
        assert not children & folded
        for cell in folded:
            assert not boundary.is_valid_cell(cell)


class TestTheSimplicityClauseIsDifferential:
    """What the new clause in the existence predicate changed, enumerated.

    Adding "and the ring closes without crossing itself" to
    :func:`itacart.boundary.is_valid_cell` is a widening of a refusal, and
    a widening has to be measured against the population it acts on
    rather than against the examples that motivated it. Two things are
    asserted over that population: the change is one-way, and every
    spelling it refuses is refused for the stated reason.

    Measured over 32 656 spellings -- the full alphabet of both polar caps
    to resolution 5, and the lateral family every thirty-seventh row to
    resolution 3 -- 196 spellings changed verdict, none in the other
    direction, and none of them carries a simple ring.
    """

    @staticmethod
    def _population() -> list[str]:
        side = CELL_SIZE_M[1]
        assert side is not None
        population: list[str] = []
        for quadrant in ("NE", "SE"):
            level_cells = [f"{quadrant}(0000/1000)"]
            for level in (2, 3, 4, 5):
                level_cells = [
                    hy._descend(cell, code)
                    for cell in level_cells
                    for code in refinement_alphabet(level)
                ]
                population.extend(level_cells)
        for quadrant in ("NE", "NW", "SE", "SW"):
            for row in range(0, 1000, 37):
                column = itacart.last_lattice_column(quadrant, row, side)
                if column < 0:
                    continue
                level_cells = [f"{quadrant}({column:04d}/{row:04d})"]
                for level in (2, 3):
                    level_cells = [
                        hy._descend(cell, code)
                        for cell in level_cells
                        for code in refinement_alphabet(level)
                    ]
                    population.extend(level_cells)
        return population

    @pytest.mark.slow
    def test_the_clause_only_ever_refuses_and_only_ever_for_a_fold(self) -> None:
        population = self._population()
        assert len(population) > 30_000

        refused = []
        for cell in population:
            ring = boundary._safe_ring(cell)[1]
            had_area = boundary.ring_area(ring) > boundary._AREA_EPSILON_M2
            exists = bool(boundary.is_valid_cell(cell))
            # One-way: nothing the old rule refused is accepted now.
            assert not (exists and not had_area), cell
            if had_area and not exists:
                refused.append(cell)

        assert refused
        for cell in refused:
            assert not boundary._is_simple_ring(
                boundary._safe_ring(cell)[1]
            ), f"{cell} was refused although its ring is simple"

    def test_an_ordinary_ring_and_a_degenerate_one_are_told_apart(self) -> None:
        """The clause itself, on rings written here rather than found."""
        square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        bowtie = [(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)]
        sliver = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]
        assert boundary._is_simple_ring(square)
        assert not boundary._is_simple_ring(bowtie)
        # Collinear and enclosing nothing is degenerate, not folded: the
        # area test is what refuses it, and this clause must not.
        assert boundary._is_simple_ring(sliver)
        assert not boundary._is_simple_ring(square[:2])


def test_a_probe_answer_that_names_no_cell_is_dropped_rather_than_kept(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The discovery keeps its own existence check, and the check still works.

    The probes quantize positions inside the parent. Inside the polar cap
    the quantizer used to answer some of them with spellings that name no
    cell, and those reached this check; it no longer does, so the check is
    exercised here by handing it one such spelling directly. The children
    must come out unchanged.
    """
    denied = "NE(0000/1000(2))"
    assert boundary.is_valid_cell(denied) is False
    hy._border_children_of.cache_clear()
    expected = list(itacart.get_children(POLAR_CAP))
    import itacart.cells as quantizer

    real = quantizer.sinusoidal_to_cell
    calls: list[int] = []

    def with_one_denied_answer(x: float, y: float, resolution: int) -> str:
        calls.append(resolution)
        return denied if len(calls) == 1 else real(x, y, resolution)

    monkeypatch.setattr(quantizer, "sinusoidal_to_cell", with_one_denied_answer)
    hy._border_children_of.cache_clear()
    try:
        assert list(itacart.get_children(POLAR_CAP)) == expected
    finally:
        hy._border_children_of.cache_clear()
    assert calls and calls[0] == 2


# --------------------------------------------------------------------------
# The discovery's own guards, exercised by injection
# --------------------------------------------------------------------------
#
# Every parent the suite visits is enumerated on the first pass, with every
# probe answered and nothing left uncovered. The guards below exist for the
# parents nobody has visited, so each is exercised by handing the discovery
# the one condition it guards against and asserting what it does with it.


@pytest.fixture
def fresh_discovery() -> Iterator[None]:
    """Discovery is cached per parent; each injection must reach it anew."""
    hy._border_children_of.cache_clear()
    yield
    hy._border_children_of.cache_clear()


def _counted_passes(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """The probe spacing of every discovery pass, in order."""
    spacings: list[float] = []
    real = hy._probe_points

    def counted(body: Polygon, spacing: float) -> list[tuple[float, float]]:
        spacings.append(spacing)
        return real(body, spacing)

    monkeypatch.setattr(hy, "_probe_points", counted)
    return spacings


def test_a_pass_that_keeps_no_candidate_halves_the_spacing_and_retries(
    monkeypatch: pytest.MonkeyPatch, fresh_discovery: None
) -> None:
    """Finding nothing on one pass is a reason to look closer, not an error."""
    expected = hy._border_children_of(POLAR_CAP, 2)
    hy._border_children_of.cache_clear()
    spacings = _counted_passes(monkeypatch)
    real = hy._shared_area
    monkeypatch.setattr(
        hy, "_shared_area", lambda a, b: 0.0 if len(spacings) == 1 else real(a, b)
    )
    assert hy._border_children_of(POLAR_CAP, 2) == expected
    assert spacings == [spacings[0], spacings[0] / 2.0]


def test_a_pass_that_leaves_the_parent_uncovered_halves_the_spacing_and_retries(
    monkeypatch: pytest.MonkeyPatch, fresh_discovery: None
) -> None:
    """Candidates that do not cover the ring are not accepted as the children.

    The first pass is made to see half of every area, so its candidates
    survive the overlap filter but cover half the parent. The coverage
    check sends it round again at half the spacing, and the second pass,
    seeing true areas, returns what an undisturbed discovery returns.
    """
    expected = hy._border_children_of(POLAR_CAP, 2)
    hy._border_children_of.cache_clear()
    spacings = _counted_passes(monkeypatch)
    real = hy._shared_area
    monkeypatch.setattr(
        hy,
        "_shared_area",
        lambda a, b: real(a, b) / 2.0 if len(spacings) == 1 else real(a, b),
    )
    assert hy._border_children_of(POLAR_CAP, 2) == expected
    assert spacings == [spacings[0], spacings[0] / 2.0]


def test_a_discovery_that_never_covers_the_parent_is_refused_after_every_pass(
    monkeypatch: pytest.MonkeyPatch, fresh_discovery: None
) -> None:
    """Only an exhausted budget is an error, and it is spent before refusing."""
    spacings = _counted_passes(monkeypatch)
    monkeypatch.setattr(hy, "_shared_area", lambda a, b: 0.0)
    with pytest.raises(GeometryError, match="could not be enumerated"):
        hy._border_children_of(POLAR_CAP, 2)
    assert len(spacings) == hy._DISCOVERY_PASSES
    assert all(later == earlier / 2.0 for earlier, later in zip(spacings, spacings[1:]))


def test_a_candidate_that_overlaps_the_others_is_refused_not_kept(
    monkeypatch: pytest.MonkeyPatch, fresh_discovery: None
) -> None:
    """The proof refuses a set of cells whose interiors overlap.

    One probe is answered with the parent itself, a real cell at the wrong
    resolution that covers everything the true children cover. Coverage
    then balances exactly -- the surplus area is the overlap -- so the only
    check that can object is the overlap check, and it must.
    """
    import itacart.cells as quantizer

    real = quantizer.sinusoidal_to_cell
    calls: list[int] = []

    def with_the_parent_once(x: float, y: float, resolution: int) -> str:
        calls.append(resolution)
        return POLAR_CAP if len(calls) == 1 else real(x, y, resolution)

    monkeypatch.setattr(quantizer, "sinusoidal_to_cell", with_the_parent_once)
    with pytest.raises(GeometryError, match="overlap by"):
        hy._border_children_of(POLAR_CAP, 2)
    assert calls and calls[0] == 2


def test_a_probe_the_quantizer_refuses_is_skipped_not_fatal(
    monkeypatch: pytest.MonkeyPatch, fresh_discovery: None
) -> None:
    """A refused probe proposes nothing, and the proof decides as before."""
    expected = hy._border_children_of(POLAR_CAP, 2)
    hy._border_children_of.cache_clear()
    import itacart.cells as quantizer

    real = quantizer.sinusoidal_to_cell
    calls: list[int] = []

    def with_one_refusal(x: float, y: float, resolution: int) -> str:
        calls.append(resolution)
        if len(calls) == 1:
            raise DomainError("a probe the quantizer will not answer")
        return real(x, y, resolution)

    monkeypatch.setattr(quantizer, "sinusoidal_to_cell", with_one_refusal)
    assert hy._border_children_of(POLAR_CAP, 2) == expected
    assert len(calls) > 1 and calls[0] == 2


def test_a_shared_edge_the_engine_cannot_resolve_is_retried_on_rebuilt_rings() -> None:
    """The area falls back to rings rebuilt by a zero-width buffer.

    The engine has reported a side location conflict on sibling rings that
    share an edge. None of the rings the suite builds provokes one with the
    engine installed here, so the conflict is raised by a stand-in for the
    first ring, and the answer must be the area the rebuilt rings share --
    here a whole child, which is inside its parent.
    """
    from shapely.errors import GEOSException

    parent = Polygon(boundary.plane_ring(POLAR_CAP)[1])
    child = Polygon(boundary.plane_ring(CAP_SUBDIVISION)[1])
    refused: list[str] = []

    class Conflicted:
        def intersection(self, other: Polygon) -> Polygon:
            refused.append("intersection")
            raise GEOSException("TopologyException: side location conflict")

        def buffer(self, distance: float) -> Polygon:
            return parent.buffer(distance)

    shared = hy._shared_area(Conflicted(), child)  # type: ignore[arg-type]
    assert refused == ["intersection"]
    assert shared == pytest.approx(child.area, rel=1e-9)
    assert shared > 0.0


def test_a_prefix_that_does_not_father_the_cell_is_not_taken_for_its_parent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The lexical prefix is trusted only as far as its child relation.

    Where the prefix names an absorbing cell it is the parent only if its
    own children include the cell. No child found so far fails that, so the
    check is exercised by withholding one child from its prefix's answer.
    The resolution must then look past the prefix -- west along the row,
    then among the prefix's siblings, without asking the prefix again --
    and, finding no cell that fathers the child, refuse it by name.
    """
    last = itacart.last_lattice_column("NE", 300, itacart.cell_size(1))
    root = f"NE({last:04d}/0300)"
    prefix = next(
        cell
        for cell in hy._children_of(root)
        if boundary.absorbs_border(cell)
        and split_components(cell)[:-1] == split_components(root)
    )
    child = next(
        cell
        for cell in hy._children_of(prefix)
        if split_components(cell)[:-1] == split_components(prefix)
    )
    assert hy._parent_cell(child) == prefix
    real = hy._children_of
    asked: list[str] = []

    def withholding(cell: str) -> list[str]:
        asked.append(cell)
        return [c for c in real(cell) if not (cell == prefix and c == child)]

    monkeypatch.setattr(hy, "_children_of", withholding)
    with pytest.raises(NonExistentCellError, match="has no parent cell"):
        hy._parent_cell(child)
    siblings = [cell for cell in real(root) if cell != prefix]
    assert siblings
    assert asked.count(prefix) == 1
    assert set(siblings) <= set(asked)


def test_a_resolution_1_cell_is_answered_by_its_quadrant_without_a_search(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The neighbour search is reached only below resolution 1.

    It reads a column from the second component, which a resolution-1 cell
    has, but it never gets that far with one: the prefix of a resolution-1
    cell is its quadrant, and every quadrant names a cell that absorbs
    nothing. Checked on one ordinary and one absorbing resolution-1 cell
    per quadrant, with the search replaced by a refusal.
    """
    side = itacart.cell_size(1)
    for quadrant in QUADRANTS:
        assert boundary.is_valid_cell(quadrant)
        assert not boundary.absorbs_border(quadrant)

    def no_search(*_: object) -> Iterator[str]:
        raise AssertionError("the neighbour search was consulted")

    monkeypatch.setattr(hy, "_neighbouring_parents", no_search)
    for quadrant in QUADRANTS:
        last = itacart.last_lattice_column(quadrant, 300, side)
        for column in (last - 1, last):
            cell = f"{quadrant}({column:04d}/0300)"
            assert boundary.is_valid_cell(cell)
            assert hy._parent_cell(cell) == quadrant
        assert boundary.absorbs_border(f"{quadrant}({last:04d}/0300)")


# --------------------------------------------------------------------------
# Two relations: lexical ancestry and physical refinement
# --------------------------------------------------------------------------
#
# Away from the border a cell's children are the codes of the next alphabet
# spelled under it, and the two relations are one. An absorbing cell also
# fathers children spelled under the next column's stem, so there the string
# prefix and the refinement part: ``contains``, ``is_ancestor`` and
# ``get_parent`` answer the prefix, and ``get_children``, compaction,
# expansion and the tree blob answer the refinement. ``normalize`` is a
# syntactic rewrite that must not change the region an index denotes.


@lru_cache(maxsize=1)
def _lateral_families() -> tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...]:
    """Absorbing resolution-1 cells with their children, split by spelling.

    The last column of every row in the four quadrants, and near the pole
    also the column before it and column zero -- the population the
    divergence was measured on. Each entry is the parent, the children
    spelled under it, and the children spelled under another stem.
    """
    side = itacart.cell_size(1)
    families = []
    for quadrant in QUADRANTS:
        first = 0 if quadrant.endswith("E") else 1
        for row in range(1001):
            last = itacart.last_lattice_column(quadrant, row, side)
            columns = {last, max(first, last - 1), 0} if row >= 999 else {last}
            for column in sorted(columns):
                parent = f"{quadrant}({column:04d}/{row:04d})"
                if column < first or not boundary.is_valid_cell(parent):
                    continue
                if not boundary.absorbs_border(parent):
                    continue
                children = tuple(hy._children_of(parent))
                stem = split_components(parent)
                own = tuple(c for c in children if split_components(c)[:-1] == stem)
                other = tuple(c for c in children if c not in own)
                families.append((parent, own, other))
    return tuple(families)


@pytest.mark.slow
def test_the_two_relations_part_on_an_enumerated_border_population() -> None:
    """Where the prefix and the refinement disagree, and how far.

    Of 4 002 absorbing resolution-1 cells, 2 392 father 3 327 children
    spelled under a stem that names no cell. For every one of them the
    lexical predicates answer the prefix and the physical ones the
    refinement, and no child spelled under a denied stem is left without
    the parent that fathers it.
    """
    families = _lateral_families()
    assert len(families) == 4_002
    parting = [family for family in families if family[2]]
    assert len(parting) == 2_392
    others = [(parent, child) for parent, _, other in parting for child in other]
    assert len(others) == 3_327
    assert len({itacart.get_parent(child) for _, child in others}) == 2_437
    for parent, child in others:
        assert itacart.contains(parent, child) is False
        assert itacart.is_ancestor(parent, child) is False
        assert not itacart.is_valid_cell(itacart.get_parent(child))
        assert hy._parent_cell(child) == parent


@pytest.mark.slow
def test_normalize_does_not_grow_a_region_under_an_absorbing_parent() -> None:
    """Lexical completeness is not physical completeness at the border.

    In 1 624 of the parting families the children spelled under the parent
    already spell the whole alphabet. Collapsing them into the parent would
    name the children spelled under the next column as well, which the
    index did not hold. ``normalize`` leaves such a family as written, so
    the region it denotes is the region it was given.
    """
    alphabet = set(refinement_alphabet(2))
    complete = [
        (parent, own)
        for parent, own, other in _lateral_families()
        if other and {split_components(c)[-1] for c in own} == alphabet
    ]
    assert len(complete) == 1_624
    for parent, own in complete:
        normal = itacart.normalize(itacart.compose(list(own)))
        assert set(itacart.uncompact_cells(normal, 2)) == set(own), parent


@pytest.mark.slow
def test_compaction_by_the_children_a_parent_has_stays_with_compact_cells() -> None:
    """The physical fold is compaction's, and it still happens."""
    for parent, own, other in _lateral_families():
        if not other:
            continue
        assert itacart.compact_cells(itacart.compose(list(own))) != parent
        assert itacart.compact_cells(itacart.compose(list(own + other))) == parent


def test_normalize_still_collapses_a_complete_node_outside_the_absorbing_family() -> (
    None
):
    """The guard concerns absorbing nodes only.

    The column west of every absorbing one, in the four quadrants and every
    row below the polar pair, holds ordinary cells, and a complete alphabet
    under each still collapses to it. So does a complete alphabet under an
    ordinary child of an absorbing cell, while the same alphabet under an
    absorbing child does not.
    """
    side = itacart.cell_size(1)
    collapsed = 0
    for quadrant in QUADRANTS:
        for row in range(999):
            column = itacart.last_lattice_column(quadrant, row, side) - 1
            parent = f"{quadrant}({column:04d}/{row:04d})"
            if not boundary.is_valid_cell(parent):
                continue
            assert not boundary.absorbs_border(parent)
            children = [f"{parent[:-1]}({code}))" for code in refinement_alphabet(2)]
            assert itacart.normalize(itacart.compose(children)) == parent
            collapsed += 1
    assert collapsed > 3_900
    ordinary, absorbing = "NE(2003/0000(1))", "NE(1784/0300(2))"
    assert not boundary.absorbs_border(ordinary) and boundary.absorbs_border(absorbing)
    for cell, folds in ((ordinary, True), (absorbing, False)):
        stem = split_components(cell)
        children = [join_components([*stem, code]) for code in refinement_alphabet(3)]
        normal = itacart.normalize(itacart.compose(children))
        assert (normal == cell) is folds, cell


@pytest.mark.slow
def test_under_the_caps_the_relations_part_by_the_measured_counts() -> None:
    """The caps part the two relations deeper than the lateral border does.

    Walked by refinement from both caps to resolution 5, every cell's
    physical parent is the cell that returned it. By resolution from 2 to
    5 there are 2, 20, 88 and 2 388 cells; 0, 0, 12 and 228 of them are
    spelled under a prefix that names no cell, and for those ``contains``
    of the physical parent answers false; 0, 0, 0 and 284 more carry such a
    prefix further up while being spelled under their physical parent.
    """
    cells = [0, 0, 0, 0]
    own_prefix_denied = [0, 0, 0, 0]
    denied_further_up = [0, 0, 0, 0]
    for cap in (POLAR_CAP, "SE(0000/1000)"):
        frontier = [cap]
        for level in range(4):
            following = []
            for parent in frontier:
                for child in hy._children_of(parent):
                    following.append(child)
                    cells[level] += 1
                    assert hy._parent_cell(child) == parent
                    stem = split_components(child)
                    prefixes = [join_components(stem[:k]) for k in range(2, len(stem))]
                    if boundary.is_valid_cell(prefixes[-1]):
                        denied = not all(boundary.is_valid_cell(p) for p in prefixes)
                        denied_further_up[level] += denied
                        assert prefixes[-1] == parent or not denied
                    else:
                        own_prefix_denied[level] += 1
                        assert itacart.contains(parent, child) is False
            frontier = following
    assert cells == [2, 20, 88, 2_388]
    assert own_prefix_denied == [0, 0, 12, 228]
    assert denied_further_up == [0, 0, 0, 284]


@pytest.mark.slow
def test_the_quantizer_nests_by_refinement_in_the_caps() -> None:
    """Each answer in the caps refines the coarser answer at the same place.

    At 2 450 positions -- 25 latitudes from 89.9824 degrees to the pole, 49
    longitudes, both hemispheres -- the physical parent of the answer at
    every resolution from 2 to 6 is the answer one resolution coarser. Its
    lexical prefix is not, at 468, 468 and 92 positions from resolution 4
    to 6, which is the lexical relation parting from the physical one.
    """
    latitudes = [89.9824 + k * (90.0 - 89.9824) / 24 for k in range(25)]
    longitudes = [-180.0 + 7.5 * k for k in range(49)]
    lexical_misses = {resolution: 0 for resolution in range(2, 7)}
    for sign in (1.0, -1.0):
        for lat in latitudes:
            for lon in longitudes:
                answers = {
                    r: itacart.geo_to_cell(lon, sign * lat, r) for r in range(1, 7)
                }
                for resolution in range(2, 7):
                    finer, coarser = answers[resolution], answers[resolution - 1]
                    assert hy._parent_cell(finer) == coarser, (lon, sign * lat, finer)
                    if join_components(split_components(finer)[:-1]) != coarser:
                        lexical_misses[resolution] += 1
    assert lexical_misses == {2: 0, 3: 0, 4: 468, 5: 468, 6: 92}
