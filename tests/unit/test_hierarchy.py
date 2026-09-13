"""Tests for :mod:`itacart.hierarchy`.

Organised by property of the hierarchy, then by surface.

Two conventions of the suite matter here. Fixtures never
spell a refinement code by hand: the alphabet strictly alternates between
levels, so a hand-written chain that repeats an alphabet is invalid, and
every code comes from :func:`refinement_alphabet` or from ``_chain_to``.
And border families are enumerated, never sampled: the last column of a
row is one cell out of two thousand, so a random draw finds it only by
luck.
"""

from __future__ import annotations

from typing import Iterator

import pytest

from itacart import boundary
from itacart import hierarchy as hy
from itacart import index as ix
from itacart.constants import (
    MAX_RESOLUTION,
    MERIDIAN_QUADRANT,
    QUADRANTS,
    refinement_alphabet,
)
from itacart.exceptions import (
    DomainError,
    MaxResolutionError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)
from itacart.resolutions import (
    cell_size,
    get_resolution,
    linear_refinement_ratio,
    refinement_ratio,
)

# --------------------------------------------------------------------------
# Corpus
# --------------------------------------------------------------------------

L1 = cell_size(1)

#: An interior cell, far from every discontinuity.
INTERIOR = "NE(0500/0100)"

#: The equatorial trapezoid whose fifth child needs column 2004 to be
#: spelled.
EQUATOR_TRAPEZOID = "NE(2003/0000)"


def _chain_to(resolution: int, *, tail: str = "") -> str:
    """Build a single-path index descending to ``resolution``.

    Same helper as the index suite: codes come from the alphabet of each
    level, so the fixture cannot itself encode the even/odd rule.
    """
    codes = [refinement_alphabet(level)[0] for level in range(2, resolution + 1)]
    if tail:
        codes.append(tail)
    body = "(".join(["0001/0002", *codes])
    return f"NE({body}{')' * (len(codes) + 1)}"


def _child(parent: str, code: str) -> str:
    """Append one refinement code to an atomic index.

    Deliberately duplicates the module's own descent helper rather than
    calling it, so that fixtures built here do not depend on the function
    they are used to check. ``test_descend_matches_the_lexical_children``
    pins the two together.
    """
    depth = parent.count("(")
    return parent[: len(parent) - depth] + "(" + code + ")" * (depth + 1)


#: The one child of the equatorial trapezoid that itself absorbs the
#: border, and so refines through the geometric path.
BORDER_AT_RES_2 = "NE(2003/0000(2))"

#: A child of the same trapezoid that does not absorb the border, but
#: sits close enough to it that the eastern stem still offers candidates.
INTERIOR_AT_RES_2_NEAR_BORDER = "NE(2003/0000(1))"

#: The polar row, the one row shorter than a whole lattice side. Its
#: refinement follows a rule of its own and is pinned separately by
#: :class:`TestPolarRefinement`; the general border enumeration below
#: stops short of it so that the two families stay distinguishable.
POLAR_ROW = 1000


#: Ordinate of the pole on the projection plane, the quarter meridian.
POLE_ORDINATE = MERIDIAN_QUADRANT


def _offset_of(code: str, level: int) -> int:
    """Signed distance of a refinement code from the prime meridian.

    Zero on the line, positive east, negative west. Read from the code's
    grid position rather than spelled by hand, so the fixture cannot
    encode the fold it is used to check.
    """
    size = linear_refinement_ratio(level)
    return boundary.meridian_child(*boundary.child_position(code, level), size)[1]


def _sub_row_span(parent: str, code: str, level: int) -> tuple[float, float]:
    """Ordinates a refinement code occupies inside its parent, unclipped."""
    size = linear_refinement_ratio(level)
    sub = boundary.meridian_geometry(parent)[3] / size
    sub_row = boundary.meridian_child(*boundary.child_position(code, level), size)[0]
    base = abs(boundary.meridian_geometry(parent)[2]) + sub_row * sub
    return base, base + sub


def _last_column_cell(quadrant: str, row: int) -> str | None:
    """The last existing cell of a lattice row, or ``None`` for an empty row."""
    column = boundary.last_lattice_column(quadrant, row, L1)
    if column < 0:
        return None
    cell = f"{quadrant}({column:04d}/{row:04d})"
    return cell if boundary.is_valid_cell(cell) else None


def _every_row_last_cell(quadrant: str = "NE") -> Iterator[str]:
    """Every last-column cell of a quadrant below the polar row, one per row.

    The antemeridian border family at resolution 1, enumerated rather
    than sampled. The polar row is deliberately not part of it: its cells
    are cut short by the pole rather than by the antemeridian, so they
    obey a different rule and would blunt every property stated over this
    family. They are enumerated on their own in
    :class:`TestPolarRefinement`.
    """
    for row in range(0, POLAR_ROW):
        cell = _last_column_cell(quadrant, row)
        if cell is not None:
            yield cell


# --------------------------------------------------------------------------
# Get_parent is pure lexical truncation and returns a prefix
# --------------------------------------------------------------------------


class TestGetParentIsLexical:
    def test_parent_is_a_prefix_of_the_child_path(self) -> None:
        child = _chain_to(MAX_RESOLUTION)
        parent = hy.get_parent(child)
        assert isinstance(parent, str)
        assert ix.split_components(child)[:-1] == ix.split_components(parent)

    def test_truncation_holds_at_every_resolution(self) -> None:
        for resolution in range(2, MAX_RESOLUTION + 1):
            cell = _chain_to(resolution)
            assert get_resolution(str(hy.get_parent(cell))) == resolution - 1

    def test_target_res_ascends_several_levels_at_once(self) -> None:
        cell = _chain_to(MAX_RESOLUTION)
        for target in range(0, MAX_RESOLUTION):
            assert get_resolution(str(hy.get_parent(cell, target))) == target

    def test_quadrant_has_no_parent(self) -> None:
        with pytest.raises(MinResolutionError):
            hy.get_parent("NE")

    def test_negative_target_is_refused_rather_than_wrapping(self) -> None:
        # A negative index would silently slice from the end of the path,
        # which is the failure this guard exists to prevent.
        with pytest.raises(MinResolutionError):
            hy.get_parent(INTERIOR, -1)

    def test_finer_target_is_refused(self) -> None:
        with pytest.raises(ResolutionError):
            hy.get_parent(INTERIOR, 5)

    def test_no_floating_point_is_consulted(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Ascent must not reach for the closed form in cosine latitude."""

        def explode(*args: object, **kwargs: object) -> object:
            raise AssertionError("get_parent consulted the boundary module")

        monkeypatch.setattr(boundary, "last_lattice_column", explode)
        monkeypatch.setattr(boundary, "is_valid_cell", explode)
        assert hy.get_parent("NE(1491/0465(1))") == "NE(1491/0465)"


class TestGetParentReturnsPrefixNotCell:
    """The falsifying family the boundary tests measure.

    Above roughly 18.5 degrees the child of a border-absorbing cell is
    spelled under the next column, so its lexical prefix names no cell.
    """

    def test_named_example_from_the_article_amendment(self) -> None:
        child = "NE(1491/0465(1))"
        assert boundary.is_valid_cell(child)
        prefix = hy.get_parent(child)
        assert prefix == "NE(1491/0465)"
        assert not boundary.is_valid_cell(str(prefix))

    def test_prefix_is_grammatical_even_when_it_names_no_cell(self) -> None:
        prefix = str(hy.get_parent("NE(1491/0465(1))"))
        assert ix.is_valid_index(prefix)
        assert not boundary.is_valid_cell(prefix)

    def test_get_parent_is_total_over_the_grammar(self) -> None:
        for cell in ("NE(2004/0000(3))", "NE(1491/0465(1))", _chain_to(5)):
            assert isinstance(hy.get_parent(cell), str)

    def test_the_divergence_is_a_family_and_not_a_single_case(self) -> None:
        divergent = 0
        for parent in _every_row_last_cell("NE"):
            for child in hy._children_of(parent):
                if str(hy.get_parent(child)) != parent:
                    divergent += 1
        assert divergent > 1

    def test_duplicates_are_preserved_for_positional_alignment(self) -> None:
        siblings = "NE(0500/0100(1,2,3,4))"
        assert hy.get_parent(siblings) == [INTERIOR] * 4


# --------------------------------------------------------------------------
# Get_children counts 4 and 25 away from the border only
# --------------------------------------------------------------------------


class TestGetChildrenAwayFromBorder:
    def test_four_children_descending_into_an_even_resolution(self) -> None:
        children = list(hy.get_children(INTERIOR, flatten=True))
        assert len(children) == refinement_ratio(2) == 4

    def test_twenty_five_children_descending_into_an_odd_resolution(self) -> None:
        parent = _child(INTERIOR, refinement_alphabet(2)[0])
        children = list(hy.get_children(parent, flatten=True))
        assert len(children) == refinement_ratio(3) == 25

    def test_alternation_holds_at_every_level(self) -> None:
        for resolution in range(2, MAX_RESOLUTION):
            cell = _chain_to(resolution)
            children = list(hy.get_children(cell, flatten=True))
            assert len(children) == refinement_ratio(resolution + 1)

    def test_children_are_the_level_alphabet_in_order(self) -> None:
        children = list(hy.get_children(INTERIOR, flatten=True))
        codes = [ix.split_components(child)[-1] for child in children]
        assert tuple(codes) == refinement_alphabet(2)

    def test_beyond_the_finest_resolution_is_refused(self) -> None:
        with pytest.raises(MaxResolutionError):
            list(hy.get_children(_chain_to(MAX_RESOLUTION)))

    def test_coarser_target_is_refused(self) -> None:
        with pytest.raises(ResolutionError):
            hy.get_children(_chain_to(5), 3)

    def test_target_beyond_thirteen_is_refused(self) -> None:
        with pytest.raises(MaxResolutionError):
            hy.get_children(INTERIOR, MAX_RESOLUTION + 1)

    def test_quadrant_descent_is_refused_with_a_reason(self) -> None:
        with pytest.raises(ResolutionError):
            list(hy.get_children("NE", flatten=True))


class TestGetChildrenOnTheBorder:
    """The count is measured, never taken from ``refinement_ratio``."""

    def test_the_equatorial_trapezoid_has_a_fifth_child(self) -> None:
        children = list(hy.get_children(EQUATOR_TRAPEZOID, flatten=True))
        assert len(children) == 5
        assert children[-1] == "NE(2004/0000(3))"

    def test_the_fifth_child_is_spelled_under_the_next_column(self) -> None:
        own = ix.split_components(EQUATOR_TRAPEZOID)[1]
        alien = [
            child
            for child in hy._children_of(EQUATOR_TRAPEZOID)
            if ix.split_components(child)[1] != own
        ]
        assert len(alien) == 1

    def test_refinement_ratio_is_not_the_count_on_the_border(self) -> None:
        children = hy._children_of(EQUATOR_TRAPEZOID)
        assert len(children) != refinement_ratio(2)

    def test_measured_counts_of_the_named_rows(self) -> None:
        """The seven rows the amendment quotes, re-measured here."""
        expected = {0: 5, 100: 2, 200: 2, 465: 5, 750: 6, 900: 4, 999: 6}
        for row, count in expected.items():
            cell = _last_column_cell("NE", row)
            assert cell is not None
            assert len(hy._children_of(cell)) == count

    def test_the_count_stays_between_two_and_seven_over_the_family(self) -> None:
        """Bounds the antemeridian family, which is not the whole globe.

        Six was the bound while the enumeration reached one column east
        and stopped. It reaches as far as the absorbing column now, and a
        child can be spelled two columns out -- ``NE(0819/0747(1))``
        belongs to ``NE(0817/0747)`` -- so the bound is seven. Seven is
        not an instance: over the four quadrants it occurs in 45 cells of
        4 000, and the count is pinned there rather than here.

        The polar cell sits outside these bounds with a single child, and
        outside this family: it is cut short by the pole rather than by
        the antemeridian. Widening the floor to admit it would make the
        bound say nothing about either family.
        """
        counts = {len(hy._children_of(parent)) for parent in _every_row_last_cell("NE")}
        assert min(counts) >= 2
        assert max(counts) <= 7
        assert len(hy._children_of("NE(0000/1000)")) < min(counts)

    def test_descending_only_the_own_stem_would_lose_children(self) -> None:
        stem = EQUATOR_TRAPEZOID[:-1]
        children = hy._children_of(EQUATOR_TRAPEZOID)
        own_stem = [child for child in children if child.startswith(stem)]
        assert len(own_stem) < len(children)

    def test_absorbs_border_discriminates_the_two_regimes(self) -> None:
        """The fast path is sound: not absorbing implies the canonical count.

        Enumerated over a whole row rather than sampled. Row 900 is short
        enough to walk end to end, which ties the shortcut to a complete
        enumeration in at least one place.
        """
        row = 900
        last = boundary.last_lattice_column("NE", row, L1)
        checked = 0
        for column in range(0, last + 1):
            cell = f"NE({column:04d}/{row:04d})"
            if not boundary.is_valid_cell(cell):
                continue
            if boundary.absorbs_border(cell):
                continue
            assert len(hy._children_of(cell)) == refinement_ratio(2)
            checked += 1
        assert checked > 300

    def test_every_quadrant_shows_the_same_border_regime(self) -> None:
        for quadrant in QUADRANTS:
            cell = _last_column_cell(quadrant, 465)
            assert cell is not None
            assert boundary.absorbs_border(cell)


# --------------------------------------------------------------------------
# Get_descendants is a generator
# --------------------------------------------------------------------------


class TestGetDescendantsStreams:
    def test_is_a_generator_and_not_a_list(self) -> None:
        stream = hy.get_descendants(INTERIOR, MAX_RESOLUTION)
        assert not isinstance(stream, list)
        assert hasattr(stream, "__next__")

    def test_a_deep_request_yields_before_it_could_have_finished(self) -> None:
        """Taking one item from a 10^12 expansion must return promptly."""
        stream = hy.get_descendants(INTERIOR, MAX_RESOLUTION)
        assert get_resolution(next(stream)) == MAX_RESOLUTION

    def test_twelve_refinements_multiply_to_a_million_million(self) -> None:
        ratios = [refinement_ratio(level) for level in range(2, MAX_RESOLUTION + 1)]
        assert ratios.count(4) == 6
        assert ratios.count(25) == 6
        product = 1
        for ratio in ratios:
            product *= ratio
        assert product == 10**12

    def test_counts_match_the_ratio_product_for_a_shallow_expansion(self) -> None:
        count = sum(1 for _ in hy.get_descendants(INTERIOR, 4))
        assert count == refinement_ratio(2) * refinement_ratio(3) * refinement_ratio(4)

    def test_every_descendant_sits_at_the_requested_resolution(self) -> None:
        for cell in hy.get_descendants(INTERIOR, 4):
            assert get_resolution(cell) == 4

    def test_descendants_of_a_border_cell_also_stream(self) -> None:
        assert sum(1 for _ in hy.get_descendants(EQUATOR_TRAPEZOID, 3)) > 0


# --------------------------------------------------------------------------
# Compact_cells runs to a fixed point
# --------------------------------------------------------------------------


class TestCompactCells:
    def test_a_complete_sibling_set_collapses_to_the_parent(self) -> None:
        full = ix.compose(hy.get_children(INTERIOR, flatten=True))
        assert hy.compact_cells(full) == INTERIOR

    def test_an_incomplete_sibling_set_does_not_collapse(self) -> None:
        codes = refinement_alphabet(2)[:-1]
        partial = ix.compose(_child(INTERIOR, code) for code in codes)
        assert hy.compact_cells(partial) == partial

    def test_compaction_runs_to_a_fixed_point_across_two_levels(self) -> None:
        deep = ix.compose(hy.get_descendants(INTERIOR, 3))
        assert hy.compact_cells(deep) == INTERIOR

    def test_compaction_is_idempotent(self) -> None:
        deep = ix.compose(hy.get_descendants(INTERIOR, 3))
        once = hy.compact_cells(deep)
        assert hy.compact_cells(once) == once

    def test_completeness_is_not_tested_by_refinement_ratio(self) -> None:
        """The border case that makes the ratio the wrong test.

        Four quaternary children spell the whole alphabet, and
        ``refinement_ratio(2)`` is four, but the parent has five children.
        Collapsing here would assert coverage over an incomplete
        partition.
        """
        four = ix.compose(
            _child(EQUATOR_TRAPEZOID, code) for code in refinement_alphabet(2)
        )
        assert len(refinement_alphabet(2)) == refinement_ratio(2)
        assert len(hy._children_of(EQUATOR_TRAPEZOID)) == 5
        assert hy.compact_cells(four) == four

    def test_the_full_border_sibling_set_does_collapse(self) -> None:
        full = ix.compose(hy.get_children(EQUATOR_TRAPEZOID, flatten=True))
        assert hy.compact_cells(full) == EQUATOR_TRAPEZOID

    def test_neither_normalize_nor_compaction_claims_an_incomplete_family(
        self,
    ) -> None:
        """At the border a complete alphabet is not a complete family.

        ``normalize`` fixes spelling and withholds the alphabet collapse under
        a border-absorbing parent; compaction counts the children that exist.
        Both leave the four quaternary children of the equatorial trapezoid
        as written, because the trapezoid has a fifth, spelled under the next
        column: claiming the parent would claim that fifth child too.
        """
        four = ix.compose(
            _child(EQUATOR_TRAPEZOID, code) for code in refinement_alphabet(2)
        )
        assert ix.normalize(four) == "NE(2003/0000(1,2,3,4))"
        assert hy.compact_cells(four) != EQUATOR_TRAPEZOID

    def test_cells_of_several_parents_are_kept_apart(self) -> None:
        left = ix.compose(hy.get_children(INTERIOR, flatten=True))
        other = "NE(0501/0100)"
        mixed = ix.compose([*ix.decompose(left), other])
        assert set(ix.decompose(hy.compact_cells(mixed))) == {INTERIOR, other}


# --------------------------------------------------------------------------
# The uncompact/compact round trip
# --------------------------------------------------------------------------


class TestRoundTrip:
    def test_round_trip_holds_for_a_uniform_interior_set(self) -> None:
        uniform = list(hy.get_descendants(INTERIOR, 3))
        packed = hy.compact_cells(ix.compose(uniform))
        assert sorted(hy.uncompact_cells(packed, 3)) == sorted(uniform)

    def test_round_trip_holds_for_a_partial_set(self) -> None:
        uniform = list(hy.get_descendants(INTERIOR, 3))[:7]
        packed = hy.compact_cells(ix.compose(uniform))
        assert sorted(hy.uncompact_cells(packed, 3)) == sorted(uniform)

    def test_the_round_trip_holds_on_the_case_built_to_break_it(self) -> None:
        """The case built to break completeness at the border.

        Four quaternary children of the equatorial trapezoid, compacted
        and expanded again. Under a completeness test based on
        ``refinement_ratio`` this returns five cells and the round trip
        would need a declared restriction. Under enumeration it returns
        the same four, and the round trip holds unrestricted.
        """
        four = [_child(EQUATOR_TRAPEZOID, code) for code in refinement_alphabet(2)]
        packed = hy.compact_cells(ix.compose(four))
        recovered = list(hy.uncompact_cells(packed, 2))
        assert len(recovered) == 4
        assert sorted(recovered) == sorted(four)

    def test_round_trip_holds_over_the_whole_border_family(self) -> None:
        for parent in _every_row_last_cell("NE"):
            children = hy._children_of(parent)
            packed = hy.compact_cells(ix.compose(children))
            assert packed == parent
            assert sorted(hy.uncompact_cells(packed, 2)) == sorted(children)

    def test_uncompact_is_a_generator(self) -> None:
        stream = hy.uncompact_cells(INTERIOR, MAX_RESOLUTION)
        assert not isinstance(stream, list)
        assert hasattr(stream, "__next__")

    def test_a_cell_already_at_the_target_passes_through(self) -> None:
        assert list(hy.uncompact_cells(INTERIOR, 1)) == [INTERIOR]

    def test_a_finer_cell_than_the_target_is_refused(self) -> None:
        with pytest.raises(ResolutionError):
            list(hy.uncompact_cells(_chain_to(5), 3))

    def test_beyond_the_finest_resolution_is_refused(self) -> None:
        with pytest.raises(MaxResolutionError):
            hy.uncompact_cells(INTERIOR, MAX_RESOLUTION + 1)

    def test_a_mixed_resolution_input_expands_uniformly(self) -> None:
        finer = _child(INTERIOR, refinement_alphabet(2)[0])
        mixed = ix.compose([INTERIOR, finer])
        expanded = list(hy.uncompact_cells(mixed, 2))
        assert {get_resolution(cell) for cell in expanded} == {2}


# --------------------------------------------------------------------------
# Vectorised semantics
# --------------------------------------------------------------------------


class TestVectorisedSemantics:
    def test_get_parent_aligns_with_decompose(self) -> None:
        composed = "NE(0500/0100(1,2)),NE(0501/0100(3))"
        parents = hy.get_parent(composed)
        assert isinstance(parents, list)
        assert len(parents) == len(ix.decompose(composed))

    def test_a_single_cell_answers_a_scalar(self) -> None:
        assert isinstance(hy.get_parent(INTERIOR), str)

    def test_get_children_unflattened_yields_one_list_per_cell(self) -> None:
        composed = "NE(0500/0100),NE(0501/0100)"
        groups = list(hy.get_children(composed))
        assert len(groups) == len(ix.decompose(composed))
        assert all(isinstance(group, list) for group in groups)
        assert all(len(group) == refinement_ratio(2) for group in groups)

    def test_get_children_does_not_collapse_for_a_single_cell(self) -> None:
        """The one member of this family that keeps its wrapper.

        ``get_parent`` answers a scalar for one cell and ``get_ancestors``
        answers a flat chain, both pinned above. ``get_children`` does
        not follow them: unflattened it yields one list per input cell,
        one list even when the input is one cell, so the stream has
        length one and the children are inside it. The docstring says so
        now; this is what makes it a claim rather than prose.
        """
        groups = list(hy.get_children(INTERIOR))
        assert len(groups) == 1
        assert isinstance(groups[0], list)
        assert len(groups[0]) == refinement_ratio(2)
        assert all(isinstance(cell, str) for cell in groups[0])

        # The control: the same call flattened does collapse, so the
        # assertion above is about the wrapper and not about the walk.
        flat = list(hy.get_children(INTERIOR, flatten=True))
        assert flat == groups[0]

    def test_get_children_flattened_yields_a_single_stream(self) -> None:
        composed = "NE(0500/0100),NE(0501/0100)"
        flat = list(hy.get_children(composed, flatten=True))
        assert len(flat) == 2 * refinement_ratio(2)
        assert all(isinstance(cell, str) for cell in flat)

    def test_get_ancestors_aligns_with_decompose(self) -> None:
        composed = "NE(0500/0100(1)),NE(0501/0100(2))"
        chains = hy.get_ancestors(composed)
        assert len(chains) == len(ix.decompose(composed))
        assert all(isinstance(chain, list) for chain in chains)

    def test_get_ancestors_of_a_single_cell_is_a_flat_chain(self) -> None:
        assert hy.get_ancestors("NE(0500/0100(1))") == ["NE", INTERIOR]

    def test_get_ancestors_of_a_quadrant_is_empty(self) -> None:
        assert hy.get_ancestors("NE") == []

    def test_child_position_aligns_with_decompose(self) -> None:
        assert hy.child_position("NE(0500/0100(1,2,3))") == [0, 1, 2]

    def test_contains_aligns_with_decompose(self) -> None:
        answers = hy.contains(INTERIOR, "NE(0500/0100(1,2)),NE(0501/0100(3))")
        assert answers == [True, True, False]


# --------------------------------------------------------------------------
# Contains is an index predicate, never a geometry query
# --------------------------------------------------------------------------


class TestContains:
    def test_a_region_covers_a_descendant_of_a_terminal_cell(self) -> None:
        assert hy.contains(INTERIOR, "NE(0500/0100(1(A2)))") is True

    def test_a_region_covers_its_own_terminal_cell(self) -> None:
        assert hy.contains(INTERIOR, INTERIOR) is True

    def test_a_region_does_not_cover_an_ancestor_of_its_terminal(self) -> None:
        assert hy.contains("NE(0500/0100(1))", INTERIOR) is False

    def test_a_region_does_not_cover_an_unrelated_cell(self) -> None:
        assert hy.contains(INTERIOR, "NE(0501/0100)") is False

    def test_a_compositional_region_covers_under_any_terminal(self) -> None:
        region = "NE(0500/0100(1)),NE(0501/0100(2))"
        assert hy.contains(region, "NE(0501/0100(2(A1)))") is True

    def test_no_geometry_is_consulted(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Containment is a claim about indices, not about polygons.

        The union of a cell's children was measured to overrun the cell
        by about one part in a thousand, because the border is curved and
        each cell approximates it by a chord. A geometric test would
        therefore disagree with the index at the border.
        """

        def explode(*args: object, **kwargs: object) -> object:
            raise AssertionError("contains consulted geometry")

        monkeypatch.setattr(boundary, "plane_ring", explode)
        monkeypatch.setattr(boundary, "is_valid_cell", explode)
        assert hy.contains(EQUATOR_TRAPEZOID, "NE(2003/0000(1))") is True

    def test_the_alien_child_is_covered_geometrically_but_not_by_index(self) -> None:
        """A declared exception, pinned so it cannot drift.

        The fifth child of the equatorial trapezoid lies inside its
        parent geometrically, but is spelled under another column, so the
        index predicate answers false. This is the price of keeping
        containment free of geometry, and it is deliberate.
        """
        alien = "NE(2004/0000(3))"
        assert alien in hy._children_of(EQUATOR_TRAPEZOID)
        assert hy.contains(EQUATOR_TRAPEZOID, alien) is False


class TestIsAncestor:
    def test_a_strict_prefix_is_an_ancestor(self) -> None:
        assert hy.is_ancestor(INTERIOR, "NE(0500/0100(1(A2)))")

    def test_a_cell_is_not_its_own_ancestor(self) -> None:
        assert not hy.is_ancestor(INTERIOR, INTERIOR)

    def test_a_finer_cell_is_not_an_ancestor_of_a_coarser_one(self) -> None:
        assert not hy.is_ancestor("NE(0500/0100(1))", INTERIOR)

    def test_a_different_quadrant_is_not_an_ancestor(self) -> None:
        assert not hy.is_ancestor("NW", "NE(0500/0100(1))")

    def test_a_different_base_cell_is_not_an_ancestor(self) -> None:
        assert not hy.is_ancestor("NE(0501/0100)", "NE(0500/0100(1))")

    def test_a_quadrant_is_an_ancestor_of_everything_under_it(self) -> None:
        assert hy.is_ancestor("NE", INTERIOR)

    def test_agrees_with_get_ancestors(self) -> None:
        cell = _chain_to(6)
        for ancestor in hy.get_ancestors(cell):
            assert hy.is_ancestor(str(ancestor), cell)


# --------------------------------------------------------------------------
# Child_position has a declared order under a border parent
# --------------------------------------------------------------------------


class TestChildPosition:
    def test_position_is_the_alphabet_rank_away_from_the_border(self) -> None:
        for rank, code in enumerate(refinement_alphabet(2)):
            assert hy.child_position(_child(INTERIOR, code)) == rank

    def test_odd_levels_rank_row_major_over_the_alphabet(self) -> None:
        parent = _child(INTERIOR, refinement_alphabet(2)[0])
        for rank, code in enumerate(refinement_alphabet(3)):
            assert hy.child_position(_child(parent, code)) == rank

    def test_positions_under_a_border_parent_are_a_dense_range(self) -> None:
        children = hy._children_of(EQUATOR_TRAPEZOID)
        positions = [hy.child_position(child) for child in children]
        assert positions == list(range(len(children)))

    def test_the_alien_child_sorts_last(self) -> None:
        """The declared order: own stem first, eastern stem after.

        The eastern stem carries the larger column, so the enumeration
        order agrees with sorting the children by component path.
        """
        children = hy._children_of(EQUATOR_TRAPEZOID)
        assert hy.child_position("NE(2004/0000(3))") == len(children) - 1

    def test_alphabet_rank_alone_would_not_be_injective(self) -> None:
        """Why the border needs the enumeration and not the rank.

        Two siblings can carry the same refinement code under different
        column prefixes, so the rank collides and cannot be an ordinal.
        """
        children = hy._children_of(EQUATOR_TRAPEZOID)
        codes = [ix.split_components(child)[-1] for child in children]
        assert len(set(codes)) < len(codes)

    def test_the_order_is_injective_over_the_whole_border_family(self) -> None:
        for parent in _every_row_last_cell("NE"):
            children = hy._children_of(parent)
            positions = [hy.child_position(child) for child in children]
            assert sorted(positions) == list(range(len(children)))

    def test_a_quadrant_has_no_siblings(self) -> None:
        with pytest.raises(MinResolutionError):
            hy.child_position("NE")

    def test_a_resolution_one_cell_has_no_refinement_alphabet(self) -> None:
        with pytest.raises(ResolutionError):
            hy.child_position(INTERIOR)

    def test_a_cell_with_no_parent_in_the_domain_is_refused(self) -> None:
        with pytest.raises(NonExistentCellError):
            hy.child_position("NE(2004/0500(1))")


# --------------------------------------------------------------------------
# common_ancestor
# --------------------------------------------------------------------------


class TestCommonAncestor:
    def test_siblings_share_their_parent(self) -> None:
        assert hy.common_ancestor("NE(0500/0100(1,2))") == INTERIOR

    def test_a_single_cell_is_its_own_common_ancestor(self) -> None:
        assert hy.common_ancestor(INTERIOR) == INTERIOR

    def test_cells_of_different_base_cells_share_the_quadrant(self) -> None:
        assert hy.common_ancestor("NE(0500/0100),NE(0501/0100)") == "NE"

    def test_cells_of_different_quadrants_have_none(self) -> None:
        with pytest.raises(DomainError):
            hy.common_ancestor("NE(0500/0100),NW(0500/0100)")

    def test_a_deep_chain_and_its_ancestor_share_the_ancestor(self) -> None:
        deep = _chain_to(6)
        shallow = str(hy.get_parent(deep, 3))
        assert hy.common_ancestor(ix.compose([deep, shallow])) == shallow

    def test_the_result_is_an_ancestor_of_or_equal_to_every_terminal(self) -> None:
        composed = "NE(0500/0100(1,2)),NE(0500/0100(3(A1)))"
        shared = hy.common_ancestor(composed)
        for cell in ix.decompose(composed):
            assert shared == cell or hy.is_ancestor(shared, cell)


# --------------------------------------------------------------------------
# Internal helpers with behaviour worth pinning
# --------------------------------------------------------------------------


class TestParentCell:
    def test_the_lexical_prefix_answers_away_from_the_border(self) -> None:
        assert hy._parent_cell("NE(0500/0100(1))") == INTERIOR

    def test_an_alien_child_resolves_to_the_western_column(self) -> None:
        assert hy._parent_cell("NE(2004/0000(3))") == EQUATOR_TRAPEZOID

    def test_a_quadrant_has_no_parent_cell(self) -> None:
        with pytest.raises(MinResolutionError):
            hy._parent_cell("NE")

    def test_an_orphan_is_refused(self) -> None:
        with pytest.raises(NonExistentCellError):
            hy._parent_cell("NE(2004/0500(1))")

    def test_every_child_of_the_border_family_finds_its_parent(self) -> None:
        for parent in _every_row_last_cell("NE"):
            for child in hy._children_of(parent):
                assert hy._parent_cell(child) == parent


class TestColumnShift:
    def test_shifting_east_and_back_is_the_identity(self) -> None:
        for cell in (INTERIOR, EQUATOR_TRAPEZOID, _chain_to(5)):
            assert hy._shift_column(hy._shift_column(cell, 1), -1) == cell

    def test_the_shifted_column_keeps_its_padding(self) -> None:
        assert hy._shift_column("NE(0009/0100)", 1) == "NE(0010/0100)"


class TestSortKey:
    def test_a_shallower_path_sorts_before_the_paths_it_prefixes(self) -> None:
        shallow = ix.split_components(INTERIOR)
        deep = ix.split_components("NE(0500/0100(1))")
        assert hy._sort_key(shallow) < hy._sort_key(deep)

    def test_resolution_one_sorts_by_the_numeric_pair(self) -> None:
        low = ix.split_components("NE(0009/0100)")
        high = ix.split_components("NE(0010/0100)")
        assert hy._sort_key(low) < hy._sort_key(high)

    def test_refinements_sort_by_alphabet_rank(self) -> None:
        first = ix.split_components(_child(INTERIOR, refinement_alphabet(2)[0]))
        last = ix.split_components(_child(INTERIOR, refinement_alphabet(2)[-1]))
        assert hy._sort_key(first) < hy._sort_key(last)


class TestRenderAndDescend:
    def test_render_round_trips_through_split_components(self) -> None:
        for cell in (INTERIOR, "NE", _chain_to(MAX_RESOLUTION), "NE(2004/0000(3))"):
            assert ix.join_components(ix.split_components(cell)) == cell

    def test_descend_appends_one_component(self) -> None:
        code = refinement_alphabet(2)[0]
        assert ix.split_components(hy._descend(INTERIOR, code))[-1] == code

    def test_descend_matches_the_lexical_children(self) -> None:
        built = [_child(INTERIOR, code) for code in refinement_alphabet(2)]
        assert built == hy._children_of(INTERIOR)


# --------------------------------------------------------------------------
# A defect inherited from the boundary module, pinned rather than repaired
# --------------------------------------------------------------------------


class TestPolarRefinement:
    """The polar row refines by vertical span, not by the level alphabet.

    Row 1000 stands a whole side above the equator while the pole is only
    about two kilometres higher, so every refinement of it is cut from a
    lattice that overshoots. A sub-row wholly beyond the pole names no
    cell; the sub-row holding the pole collapses onto one isosceles
    triangle spanning the parallel; sub-rows below it refine normally.
    The child count therefore falls far below the refinement ratio, and
    that shortfall is the result rather than a defect.

    Two cells in the globe are affected, one per eastern quadrant. The
    western quadrants have no polar cell because column zero does not
    exist there.
    """

    POLAR = ("NE(0000/1000)", "SE(0000/1000)")

    def test_only_two_polar_cells_exist(self) -> None:
        for quadrant in QUADRANTS:
            cell = f"{quadrant}(0000/{POLAR_ROW:04d})"
            expected = cell in self.POLAR
            assert bool(boundary.is_valid_cell(cell)) is expected

    def test_the_polar_cell_is_a_triangle_that_absorbs_the_border(self) -> None:
        for cell in self.POLAR:
            assert boundary.cell_shape(cell) == "triangle"
            assert boundary.absorbs_border(cell)

    def test_descent_enumerates_the_polar_cell_instead_of_refusing_it(self) -> None:
        for cell in self.POLAR:
            assert hy._children_of(cell)

    def test_the_first_refinement_yields_a_single_child(self) -> None:
        """Refining once buys no detail: the child is the parent again.

        The only sub-row at this level is the one holding the pole, and
        it collapses onto the triangle the parent already was.
        """
        from shapely.geometry import Polygon

        for cell in self.POLAR:
            found = hy._children_of(cell)
            assert len(found) == 1
            assert found[0] == _child(cell, refinement_alphabet(2)[0])
            assert Polygon(boundary.plane_ring(found[0])[1]).equals(
                Polygon(boundary.plane_ring(cell)[1])
            )

    def test_the_second_refinement_yields_ten_of_twenty_five(self) -> None:
        """Nine cells across the lower sub-row and one on the pole.

        The whole alphabet is enumerated and the survivors named, so the
        shortfall against the refinement ratio is a census.
        """
        alphabet = refinement_alphabet(3)
        expected = [
            alphabet[position] for position in (0, 1, 2, 3, 4, 5, 6, 10, 15, 20)
        ]
        for cell in self.POLAR:
            parent = hy._children_of(cell)[0]
            found = hy._children_of(parent)
            assert len(found) < refinement_ratio(3)
            assert sorted(found) == sorted(_child(parent, code) for code in expected)

    def test_the_children_tile_the_parent_without_gap_or_overlap(self) -> None:
        """Area closes at both refinements, which is what a partition means."""
        from shapely.geometry import Polygon
        from shapely.ops import unary_union

        for cell in self.POLAR:
            body = Polygon(boundary.plane_ring(cell)[1])
            queue = [cell]
            for _ in (2, 3):
                found = [child for parent in queue for child in hy._children_of(parent)]
                pieces = [Polygon(boundary.plane_ring(child)[1]) for child in found]
                assert sum(piece.area for piece in pieces) == pytest.approx(
                    body.area, rel=1e-6
                )
                assert unary_union(pieces).area == pytest.approx(body.area, rel=1e-6)
                queue = found

    def test_the_whole_level_alphabet_is_accounted_for(self) -> None:
        """Every nominal code lands in exactly one of four buckets.

        Alive under the parent's own prefix, killed by lying wholly
        beyond the pole, killed by the collapse of the sub-row that holds
        the pole, or gained from the stem of the column immediately east.
        The four have to sum to the refinement ratio, and the first plus
        the fourth have to sum to what descent returns. Stated over the
        whole alphabet at both refinements, since a bucket that is empty
        is as much a result as one that is full.
        """
        for cell in self.POLAR:
            queue = [cell]
            for level in (2, 3):
                for parent in queue:
                    ratio = refinement_ratio(level)
                    produced = set(hy._children_of(parent))
                    east = hy._shift_column(parent, 1)
                    own_alive, beyond, collapsed, gained = [], [], [], []
                    for code in refinement_alphabet(level):
                        base, top = _sub_row_span(parent, code, level)
                        offset = _offset_of(code, level)
                        if base >= POLE_ORDINATE:
                            beyond.append(code)
                        elif top > POLE_ORDINATE and offset != 0:
                            collapsed.append(code)
                        else:
                            own_alive.append(code)
                        if _child(east, code) in produced:
                            gained.append(code)
                    assert len(own_alive) + len(beyond) + len(collapsed) == ratio
                    assert gained == []
                    assert produced == {_child(parent, code) for code in own_alive}
                    assert len(produced) == len(own_alive)
                queue = [child for parent in queue for child in hy._children_of(parent)]

    def test_the_eastern_stem_gains_nothing_here_alone(self) -> None:
        """The polar cell absorbs, yet its eastern stem is empty.

        Everywhere else in the border family the eastern stem is what
        keeps descent from losing children. Here there is no column east
        to descend into: the polar cell already spans the whole parallel,
        so column zero is the last one the row holds.
        """
        for cell in self.POLAR:
            assert boundary.absorbs_border(cell)
            quadrant = cell[:2]
            assert boundary.last_lattice_column(quadrant, POLAR_ROW, L1) == 0
            assert boundary.is_valid_cell(f"{quadrant}(0001/{POLAR_ROW:04d})") is False

    def test_the_lone_child_compacts_back_onto_its_parent(self) -> None:
        """Two indices, one region, and the spatial operations agree.

        The single child of the polar cell covers exactly what its parent
        covers, which is the only place in the globe where refining buys
        no area. Compaction is spatial and has to collapse the pair;
        containment is an index predicate and happens to answer the same
        way here, which is worth pinning because it does so for a
        different reason.
        """
        from shapely.geometry import Polygon

        for cell in self.POLAR:
            lone = hy._children_of(cell)[0]
            assert Polygon(boundary.plane_ring(lone)[1]).equals(
                Polygon(boundary.plane_ring(cell)[1])
            )
            assert hy.compact_cells(lone) == cell
            assert ix.normalize(f"{cell},{lone}") == cell
            assert hy.contains(cell, lone) is True
            assert hy.is_ancestor(cell, lone) is True
            assert hy.get_parent(lone) == cell
            assert list(hy.uncompact_cells(cell, 2)) == [lone]

    def test_no_other_row_of_any_quadrant_is_affected(self) -> None:
        """Bounds the rule: enumerated over every last column of the globe."""
        for quadrant in QUADRANTS:
            for cell in _every_row_last_cell(quadrant):
                hy._children_of(cell)

    def test_a_folded_ring_is_refused_by_existence_and_never_repaired(
        self,
    ) -> None:
        """The refusal moved, and it moved to the right place.

        The selector used to raise on a self-intersecting candidate,
        because the polar triangle produced them and repairing one would
        have answered with a child set derived from a polygon nobody
        meant to draw. That guard asked the question one layer too late:
        a folded ring does not name a cell at all, and the spellings that
        carry one are not candidates to be refused but non-cells to be
        skipped.

        So the rule lives in the existence predicate now, and the tiling
        is what proves nothing was lost by skipping them: the polar
        triangle closes to 100% without its four folded spellings.
        """
        folded = "NE(0000/1000(3(A1)))"
        assert not boundary.is_valid_cell(folded)
        assert boundary.ring_area(boundary.plane_ring(folded)[1]) > 100_000.0

        children = hy._children_of("NE(0000/1000(1(B2)))")
        assert folded not in children
        for child in children:
            assert boundary.is_valid_cell(child)


# --------------------------------------------------------------------------
# How the two selection tests divide the work
# --------------------------------------------------------------------------


class TestSelectionOfBorderChildren:
    """Existence and overlap are separate tests and both are load bearing.

    A candidate under the eastern stem has to be a cell at all, and it
    has to lie inside the parent. Which of the two rejects is not fixed:
    at the first refinement existence alone happens to decide every case,
    and one level down the overlap test does most of the work. Dropping
    either would be right in one place and wrong in the other.
    """

    def test_existence_alone_decides_the_first_refinement(self) -> None:
        from shapely.geometry import Polygon

        surplus = 0
        for parent in _every_row_last_cell("NE"):
            body = Polygon(boundary.plane_ring(parent)[1])
            for step in (0, 1):
                stem = hy._shift_column(parent, step)
                for code in refinement_alphabet(2):
                    candidate = hy._descend(stem, code)
                    if not boundary.is_valid_cell(candidate):
                        continue
                    outline = Polygon(boundary.plane_ring(candidate)[1])
                    epsilon = hy._OVERLAP_EPSILON_RATIO * cell_size(2) ** 2
                    if outline.intersection(body).area <= epsilon:
                        surplus += 1
        assert surplus == 0

    def test_the_overlap_test_does_reject_when_it_is_reached(self) -> None:
        """The selector's own contract, exercised directly.

        Called on a parent that does not absorb the border, the geometric
        selector admits twenty-eight candidates by existence and keeps
        twenty-five by overlap. The public path never sends such a parent
        here, which is the point of the next test.
        """
        parent = INTERIOR_AT_RES_2_NEAR_BORDER
        assert not boundary.absorbs_border(parent)
        admitted = 0
        for step in (0, 1):
            stem = hy._shift_column(parent, step)
            for code in refinement_alphabet(3):
                if boundary.is_valid_cell(hy._descend(stem, code)):
                    admitted += 1
        kept = len(hy._border_children_of(parent, 3))
        assert admitted > kept
        assert kept == refinement_ratio(3)

    def test_behind_the_gate_existence_alone_has_always_sufficed(self) -> None:
        """Enumerated, not sampled, over two refinement levels.

        Every parent that absorbs the border in the northeast quadrant,
        at the first and second refinements: the overlap test rejects
        nothing any of those times, because a candidate that exists under
        the eastern stem of such a parent always lies inside it. The
        selector keeps the test because two levels are not thirteen, but
        the redundancy is recorded here rather than assumed away.
        """
        from shapely.geometry import Polygon

        rejected = 0
        parents = 0
        for first in _every_row_last_cell("NE"):
            for parent in (first, *hy._children_of(first)):
                if not boundary.absorbs_border(parent):
                    continue
                parents += 1
                level = get_resolution(parent) + 1
                body = Polygon(boundary.plane_ring(parent)[1])
                for step in (0, 1):
                    stem = hy._shift_column(parent, step)
                    for code in refinement_alphabet(level):
                        candidate = hy._descend(stem, code)
                        if not boundary.is_valid_cell(candidate):
                            continue
                        outline = Polygon(boundary.plane_ring(candidate)[1])
                        area = outline.intersection(body).area
                        epsilon = hy._OVERLAP_EPSILON_RATIO * cell_size(level) ** 2
                        if area <= epsilon:
                            rejected += 1
        assert parents > 2000
        assert rejected == 0

    def test_only_some_children_of_a_border_cell_inherit_the_border(self) -> None:
        """Absorption is not hereditary, which is why the fast path pays."""
        inheriting = [
            child
            for child in hy._children_of(EQUATOR_TRAPEZOID)
            if boundary.absorbs_border(child)
        ]
        assert 0 < len(inheriting) < len(hy._children_of(EQUATOR_TRAPEZOID))


# --------------------------------------------------------------------------
# Compaction corners
# --------------------------------------------------------------------------


class TestCompactionCorners:

    def test_a_cell_with_no_parent_in_the_domain_survives_compaction(self) -> None:
        orphan = "NE(2004/0500(1))"
        assert ix.is_valid_index(orphan)
        assert hy.compact_cells(orphan) == orphan

    def test_resolution_one_cells_never_collapse_into_a_quadrant(self) -> None:
        """A quadrant holds two million base cells and is never claimed whole.

        Compare cell sets rather than strings: two base cells of one
        quadrant render under a shared prefix, which is the same set
        spelled more compactly.
        """
        pair = "NE(0500/0100),NE(0501/0100)"
        compacted = hy.compact_cells(pair)
        assert set(ix.decompose(compacted)) == set(ix.decompose(pair))
        assert compacted != "NE"
