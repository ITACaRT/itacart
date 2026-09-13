"""Hierarchical navigation: ascent, descent, compaction.

Part of OGC DGGS Core requirement 17 (topological query functions), with
the neighbour half living in :mod:`itacart.topology`.

Every function here accepts a compositional index. Where the input holds
several cells, the return is positionally aligned with
:func:`itacart.index.decompose`.

Two relations live here, and they coincide everywhere except at the
domain border.

*Lexical ancestry* relates spellings: one index is an ancestor of another
when it is a prefix of it, which costs a slice and never consults the
domain. :func:`get_parent`, :func:`is_ancestor` and :func:`contains`
answer it, and :func:`itacart.boundary.is_valid_cell` is the only arbiter
of whether a prefix they return names a cell.

*Physical refinement* relates cells: the children of a cell are the cells
that partition it. :func:`get_children`, :func:`get_descendants`,
:func:`child_position`, :func:`compact_cells` and :func:`uncompact_cells`
answer it, and so does the tree blob, whose compaction folds by the same
rule. :func:`itacart.index.normalize` answers neither: it rewrites
spelling, and withholds the one rewrite that would change the region.

Away from the border the children of a cell are the codes of the next
alphabet spelled under it, and the two relations are one. A cell whose
outer side has been carried onto the domain border reaches east past its
own nominal column, so some of its children are spelled under the next
column's resolution-1 prefix: they exist and are addressed uniquely, but
their prefix names no cell. Enumerated between resolutions 1 and 2, 2 392
of the 4 002 absorbing cells of the lateral border father 3 327 such
children under 2 437 prefixes that name no cell. Under the two caps the
same happens deeper: 12 of the 88 cells at resolution 4 and 228 of the
2 388 at resolution 5 are spelled under a prefix that names no cell. For
each of those the lexical parent is a well-formed string that names no
cell, :func:`contains` of the physical parent answers false, and the
physical parent is found by refinement. Their descendants carry the prefix
as a lexical ancestor; under the caps, the other 284 cells of resolution 5
with such an ancestor are spelled under their physical parent.

Functions that descend consult :func:`itacart.boundary.absorbs_border`
first and take a purely lexical path when it answers false, which it does
for every cell of the lattice except the border-absorbing family.

**The absorbing side is a chord, and a child's chord is not its parent's.**
This is the geometric half of the same parting, and it is a declared
limitation of this version rather than a defect awaiting repair. The
domain border is a smooth convex curve on the plane, and an absorbing cell
replaces its outer side by the chord of that curve across its own height.
A parent spans one cell height and takes one chord; each child spans half
that and takes two, and a chord of a convex curve over a shorter span lies
further out. So a child of an absorbing cell reaches past its parent's
effective ring, by construction rather than by accident, and no
enumeration order changes it.

Measured over 572 lateral parents -- the last column of every seventh row
in four quadrants -- every parent is affected and exactly two children of
each leave, which is 1 144 children. The share of a child's own area
falling outside its parent has a median of 0.038142 per cent, a mean of
0.056931, a ninetieth percentile of 0.104096 and a worst case of 2.411015
at ``SE(1930/0196)``. The excess has a closed form: the sagitta of the
border across the parent's height, times half that height, to four decimal
places. Every figure in this paragraph is fixed by a named test, and those
tests fail if the absorbing side stops being a per-cell chord.

What the approximation does not cost is stated with the same care, because
a limitation whose edges are not drawn is read as larger than it is. The
children still cover the parent: the shortfall stays under one part in a
million of the parent's area, measured at 1.1e-12 over those 572 parents
and 2.1e-16 under the caps. The criterion that bounds it is one-sided by
design -- it bounds the shortfall and says nothing about the excess, which
is this section's subject, so it must not be read as saying the cover is
exact or the summed child area a fixed multiple of the parent's. Neither
is true here. Quantisation stays total and deterministic: a position
inside the excess is answered by the child whose excess it is at the
child's resolution, and by the parent at the parent's, in a fresh
interpreter as in this one. No public name has to emit a spelling
:func:`itacart.boundary.is_valid_cell` denies because of any of it.
Physical refinement preserves the region, so :func:`itacart.index.normalize`
does not grow one and compaction still folds by the children a parent has.

What it does cost is :func:`contains`, which tests the effective chordal
ring. A query bounded by the domain border can leave out a leaf whose
effective geometry protrudes past the query's own piece, by a few square
metres or by a few thousand, and under ``contains`` that leaf is not kept.
Measured on row 300 in four quadrants, with the parent's effective ring as
the query, one leaf of three is kept at resolution 2 and 53 of 63 at
resolution 3, while ``intersects`` returns the whole family. That is the
excess asked about from outside instead of from the parent, not a second
defect and not a defect of ``contains``.

Representing the border exactly would mean carrying a side, a resolution-1
row and an interval instead of a ring. That was investigated and is not in
this version.
"""

from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING, Iterable, Iterator, Sequence, cast

from ._existence import require_existing_cells
from .constants import (
    CELL_SIZE_M,
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    QUADRANTS,
    RES1_DIGITS,
    RES1_SEPARATOR,
    refinement_alphabet,
)
from .exceptions import (
    DomainError,
    GeometryError,
    ITACaRTError,
    MaxResolutionError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)
from .index import (
    BASE_CELL_RESOLUTION,
    QUADRANT_RESOLUTION,
    compose,
    decompose,
    join_components,
    split_components,
)

if TYPE_CHECKING:  # pragma: no cover - import cycle guard for annotations
    from shapely.geometry import Polygon

__all__ = [
    "get_parent",
    "get_children",
    "get_ancestors",
    "get_descendants",
    "is_ancestor",
    "common_ancestor",
    "contains",
    "compact_cells",
    "uncompact_cells",
    "child_position",
]


_DESCENT_OPEN = "("
_DESCENT_CLOSE = ")"

_OVERLAP_EPSILON_RATIO = 1e-6

# The share of a parent's own area that may stay uncovered by the children
# found for it before the enumeration is declared incomplete. Measured, an
# exact tiling leaves parts in the 1e-9 range of the parent, so this sits
# three orders above the noise and far below one cell.
_COVERAGE_EPSILON_RATIO = 1e-6

# How many times the discovery lattice is halved before giving up. Each
# pass quadruples the number of probes, so six passes reach a spacing of
# one sixty-fourth of a cell side.
_DISCOVERY_PASSES = 6
"""Fraction of a child's nominal area below which overlap is only contact.

Border children are selected by intersecting plane rings, and two cells
that merely share an edge intersect in a sliver of rounding noise. The
threshold has to sit above that sliver and below the smallest genuine
overlap. Both of those scale with the cell, so the threshold scales with
it too, and the constant is a ratio rather than an area.

Measured, not assumed. The sliver between lexical neighbours is exactly
zero up to level 11 and at most 2.33e-9 of the nominal area at levels 12
and 13. The smallest genuine overlap, enumerated over all 4000
border-absorbing cells of resolution 1 and over 60 chains carried down
to resolution 13, is 2.05e-2 of the nominal area. This ratio sits 430
times above the first and 20500 times below the second, close to the
geometric mean of the two.

A fixed area cannot do this job. One square metre, which is what stood
here, is the entire nominal area of a resolution-9 cell and is larger
than every cell below it, so a child wholly inside its parent was
rejected for being the size it is supposed to be and the eastern border
of the grid had no refinement at all below resolution 8.
"""


# --------------------------------------------------------------------------
# Component paths
# --------------------------------------------------------------------------


def _resolution_of(components: Sequence[str]) -> int:
    """Resolution addressed by a component path.

    One component is a bare quadrant at resolution zero, two reach the
    resolution-1 lattice, and each further component descends one level.
    """
    return len(components) - 1


def _sort_key(components: Sequence[str]) -> tuple[object, ...]:
    """Deterministic ordering key for a component path.

    Quadrants order as they are declared, resolution-1 cells by the
    numeric pair, and refinements by their rank in the level alphabet.
    """
    key: list[object] = [QUADRANTS.index(components[0])]
    for position, component in enumerate(components[1:], start=BASE_CELL_RESOLUTION):
        if position == BASE_CELL_RESOLUTION:
            column, _, row = component.partition(RES1_SEPARATOR)
            key.append((int(column), int(row)))
        else:
            key.append((refinement_alphabet(position).index(component), 0))
    return tuple(key)


def _paths(index: str) -> list[list[str]]:
    """Component path of every cell addressed by an index, in write order."""
    return [split_components(cell) for cell in decompose(index)]


def _is_single(index: str) -> bool:
    """Whether an index addresses exactly one cell."""
    return len(decompose(index)) == 1


def _scalar_or_list(index: str, values: list[str]) -> str | list[str]:
    """Answer a scalar for a single-cell index and a list otherwise.

    Follows the convention of :func:`itacart.index.quadrant_of`: whoever
    asks about one cell wants one answer, not a list holding one.
    """
    if len(values) == 1 and _is_single(index):
        return values[0]
    return values


# --------------------------------------------------------------------------
# Ascent
# --------------------------------------------------------------------------


def get_parent(index: str, target_res: int | None = None) -> str | list[str]:
    """Ascend the hierarchy by lexical truncation.

    No floating-point work and no negative indices: the parent is a
    prefix of the child string, which is the property that makes the
    compositional index cheap to navigate.

    The return is a **prefix**, not necessarily a cell. For most of the
    grid the two coincide, but a cell that absorbs the domain border
    reaches past its own column and fathers children spelled under the
    next column's prefix; truncating such a child yields a well-formed
    string that names no cell. Callers that need the stronger claim pass
    the result to :func:`itacart.boundary.is_valid_cell`.

    Args:
        index: Compositional index string.
        target_res: Resolution to ascend to. Defaults to one level up.

    Returns:
        The ancestor prefix for a single cell, or a positionally aligned
        list. Duplicates are preserved so alignment holds; call
        :func:`itacart.index.normalize` on the composed result to dedupe.

    Raises:
        MinResolutionError: If ``target_res`` is above resolution 0.
        ResolutionError: If ``target_res`` is finer than the input.
    """
    if target_res is not None and target_res < MIN_RESOLUTION:
        raise MinResolutionError(
            f"resolution {target_res} is above the quadrant level {MIN_RESOLUTION}"
        )
    values: list[str] = []
    for components in _paths(index):
        current = _resolution_of(components)
        wanted = current - 1 if target_res is None else target_res
        if wanted < MIN_RESOLUTION:
            raise MinResolutionError(
                f"{join_components(components)!r} is at resolution {current} "
                "and has no coarser prefix"
            )
        if wanted > current:
            raise ResolutionError(
                f"resolution {wanted} is finer than {join_components(components)!r} "
                f"at resolution {current}"
            )
        values.append(join_components(components[: wanted + 1]))
    return _scalar_or_list(index, values)


def get_ancestors(index: str) -> list[str] | list[list[str]]:
    """Every ancestor of a cell, coarse to fine.

    Strict: a cell is not among its own ancestors. The reflexive relation
    is :func:`contains`.

    Args:
        index: Compositional index string.

    Returns:
        An ancestor chain for a single cell, or a positionally aligned
        list of chains. Entries are prefixes under the same caveat as
        :func:`get_parent`.
    """
    chains: list[list[str]] = []
    for components in _paths(index):
        chains.append(
            [
                join_components(components[: depth + 1])
                for depth in range(len(components) - 1)
            ]
        )
    if len(chains) == 1 and _is_single(index):
        return chains[0]
    return cast("list[str] | list[list[str]]", chains)


def common_ancestor(index: str) -> str:
    """Deepest cell containing every terminal cell of the index.

    Underlies :func:`compact_cells` and :func:`itacart.index.compose`, and
    is independently useful for finding the grouping node of a vertex set.

    Args:
        index: Compositional index string, or a composed list of cells.

    Returns:
        The deepest common ancestor index. Returns the quadrant when the
        cells share only that.

    Raises:
        DomainError: If the cells span more than one quadrant.
    """
    paths = _paths(index)
    quadrants = {components[0] for components in paths}
    if len(quadrants) > 1:
        raise DomainError(
            f"{index!r} spans quadrants {sorted(quadrants)} "
            "and has no common ancestor"
        )
    shared = list(paths[0])
    for components in paths[1:]:
        depth = 0
        limit = min(len(shared), len(components))
        while depth < limit and shared[depth] == components[depth]:
            depth += 1
        shared = shared[:depth]
    return join_components(shared)


# --------------------------------------------------------------------------
# Descent
# --------------------------------------------------------------------------


def _descend(cell: str, code: str) -> str:
    """Append one refinement code to an atomic index string."""
    depth = cell.count(_DESCENT_OPEN)
    head = cell[: len(cell) - depth]
    return head + _DESCENT_OPEN + code + _DESCENT_CLOSE * (depth + 1)


def _shift_column(cell: str, step: int) -> str:
    """The same index with its resolution-1 column moved east by ``step``."""
    quadrant, _, rest = cell.partition(_DESCENT_OPEN)
    column, _, tail = rest.partition(RES1_SEPARATOR)
    shifted = f"{int(column) + step:0{RES1_DIGITS}d}"
    return f"{quadrant}{_DESCENT_OPEN}{shifted}{RES1_SEPARATOR}{tail}"


def _children_of(cell: str) -> list[str]:
    """Every existing cell one resolution below ``cell``, in canonical order.

    Two regimes. A cell that does not absorb the domain border refines
    into exactly the level alphabet, and the children are read off the
    string with no geometry at all. A cell that does absorb it holds
    surface its nominal footprint does not, and that surface is spelled
    elsewhere, so the children are found by geometry against the ring the
    cell actually has.

    The order is canonical in both regimes: the component sort, which is
    the alphabet order under one prefix and a total order across
    prefixes. It does not depend on where a child was spelled, because
    under an absorbing parent there is no single stem to order by.
    """
    from . import boundary

    components = split_components(cell)
    current = _resolution_of(components)
    if current < BASE_CELL_RESOLUTION:
        raise ResolutionError(
            f"{cell!r} is a quadrant; its children are the resolution-1 "
            "lattice, addressed by a coordinate pair rather than by a "
            "refinement alphabet, and are not enumerated here"
        )
    if current >= MAX_RESOLUTION:
        raise MaxResolutionError(
            f"{cell!r} is at resolution {MAX_RESOLUTION} and has no children"
        )
    level = current + 1
    if not boundary.absorbs_border(cell):
        return [_descend(cell, code) for code in refinement_alphabet(level)]
    return list(_border_children_of(cell, level))


def _probe_points(body: "Polygon", spacing: float) -> list[tuple[float, float]]:
    """Points inside ``body`` on a lattice of ``spacing``, plus one anchor.

    The anchor guarantees at least one probe for a part thinner than the
    spacing; the lattice does the rest. Completeness is not claimed here
    and does not have to be: the caller proves it by area.
    """
    from shapely.geometry import Point

    anchor = body.representative_point()
    points = [(anchor.x, anchor.y)]
    min_x, min_y, max_x, max_y = body.bounds
    steps_x = int((max_x - min_x) / spacing) + 2
    steps_y = int((max_y - min_y) / spacing) + 2
    for i in range(steps_x):
        x = min_x + (i + 0.5) * spacing
        for j in range(steps_y):
            y = min_y + (j + 0.5) * spacing
            if body.contains(Point(x, y)):
                points.append((x, y))
    return points


def _shared_area(first: "Polygon", second: "Polygon") -> float:
    """Area common to two rings, robust to a shared edge.

    Sibling rings meet along an edge both of them carry, and the geometry
    engine reports a side location conflict on some of those pairs rather
    than an empty overlap. A zero-width buffer rebuilds each ring from
    its own edges and the second attempt succeeds; it is only reached
    when the first fails, so the ordinary pair pays nothing for it.

    This is not the repair of a folded ring. Both arguments have already
    passed :func:`itacart.boundary.is_valid_cell`, which refuses a ring
    that crosses itself, so what is being resolved here is a disagreement
    about a shared edge between two rings that each close cleanly.
    """
    from shapely.errors import GEOSException

    try:
        return float(first.intersection(second).area)
    except GEOSException:
        return float(first.buffer(0).intersection(second.buffer(0)).area)


@lru_cache(maxsize=4096)
def _border_children_of(cell: str, level: int) -> tuple[str, ...]:
    """Children of a border-absorbing cell, found geometrically and proved.

    The ring the package returns for the parent is the authority. It is
    the effective ring, after the pole and after absorption, and the
    contract this function meets is stated against it: the effective
    rings of the children cover the parent exactly and their interiors do
    not overlap.

    The children are not the descendants of the parent's nominal
    footprint. Measured on the polar cap, that footprint covers 78.77% of
    the parent and the shortfall does not shrink with depth, because the
    absorbed surface is annexed at one level and then dropped by the next
    descent. The surface the footprint misses is spelled under sibling
    codes at an intermediate level -- spellings that name no cell of their
    own -- so no shift of the resolution-1 column reaches it. The eastern
    stem this function used to try is one such spelling among many, and
    only for part of the lateral border.

    So the candidates come from the point resolver instead, which answers
    with the canonical spelling of whatever cell covers a place. A probe
    lattice over the parent's ring proposes; the area test disposes. A
    candidate whose ring folds is dropped rather than refused, because a
    folded ring is an artefact of absorbing a border steeper than the
    lattice and is not a child at all: measured on the polar triangle,
    the tiling closes exactly without the four folded spellings.

    Discovery is deliberately naive. If the probes miss surface, the
    lattice is halved and the search runs again, and only an exhausted
    budget is an error. Correctness is in the proof, not in the search,
    so the search may be replaced by a cheaper one without touching what
    this function means.
    """
    from shapely.geometry import Polygon

    from . import boundary, cells
    from .resolutions import cell_size

    body = Polygon(boundary.plane_ring(cell)[1])
    side = cell_size(level)
    overlap_epsilon = _OVERLAP_EPSILON_RATIO * side**2
    coverage_epsilon = _COVERAGE_EPSILON_RATIO * body.area
    spacing = side / 2.0

    for _ in range(_DISCOVERY_PASSES):
        # Two sources, added rather than chosen between. The probes reach
        # spellings no edit of the string produces, which is the polar
        # case. The alphabet under the own and eastern stems reaches
        # children whose overlap with the parent is a sliver too thin for
        # a probe to land in, which is the lateral case and is where the
        # old rule was right. Neither is trusted on its own; the contract
        # below is what decides.
        proposed: set[str] = set()
        for x, y in _probe_points(body, spacing):
            try:
                proposed.add(cells.sinusoidal_to_cell(x, y, level))
            except ITACaRTError:
                continue
        for step in (0, 1):
            stem = _shift_column(cell, step)
            for code in refinement_alphabet(level):
                candidate = _descend(stem, code)
                if boundary.is_valid_cell(candidate):
                    proposed.add(candidate)

        kept: list[str] = []
        rings: list["Polygon"] = []
        for candidate in proposed:
            # One authority for whether a spelling names a cell, and it
            # already refuses a folded ring. Nothing here repairs one: a
            # fold is an artefact of absorbing a border steeper than the
            # lattice, and the tiling closes without the spellings that
            # carry one.
            if not boundary.is_valid_cell(candidate):
                continue
            outline = Polygon(boundary.plane_ring(candidate)[1])
            if _shared_area(outline, body) <= overlap_epsilon:
                continue
            kept.append(candidate)
            rings.append(outline)

        if not rings:
            spacing /= 2.0
            continue

        # Areas rather than a union. Sibling rings share edges exactly,
        # and the union of such rings is where the geometry engine
        # reports a side location conflict; pairwise areas ask the same
        # question without building the shape. The two agree only while
        # the interiors are disjoint, which is the other half of the
        # contract and is measured right here.
        overlap = sum(
            _shared_area(first, second)
            for index, first in enumerate(rings)
            for second in rings[index + 1 :]
        )
        covered = sum(_shared_area(ring, body) for ring in rings)
        if body.area - covered + overlap > coverage_epsilon:
            spacing /= 2.0
            continue

        # The tolerance for a shared edge cannot be a share of the cell
        # area alone. Plane coordinates run to ten million metres, so one
        # unit in the last place is about two ten-billionths of a metre,
        # and two rings meeting along an edge disagree by a sliver that
        # wide however small the cells are. Deep in the hierarchy that
        # sliver outgrows a tolerance tied to the cell side: measured at
        # resolution 13, the sliver is 2.6e-10 square metres and the
        # side-based tolerance is 1.0e-10. The floor is therefore the
        # total edge length times that unit, with room to spare, and it
        # stays orders below any overlap a real pair of cells would show.
        unit = 2.0**-52 * max(abs(value) for value in body.bounds)
        noise = 8.0 * sum(ring.length for ring in rings) * unit
        if overlap > max(overlap_epsilon, noise):
            raise GeometryError(
                f"the children found for {cell!r} overlap by {overlap:.6f} "
                "square metres, so they do not partition it"
            )
        return tuple(sorted(kept, key=lambda c: _sort_key(split_components(c))))

    raise GeometryError(
        f"the children of {cell!r} could not be enumerated: the cells found "
        "leave part of its ring uncovered, so they do not partition it"
    )


def get_children(
    index: str, target_res: int | None = None, flatten: bool = False
) -> Iterator[str] | Iterator[list[str]]:
    """Project descendants of an index down to a target resolution.

    Yields 4 children when descending into an even resolution and 25 when
    descending into an odd one — away from the domain border. A cell
    whose outer side has been carried onto the border was measured to
    hold between two and six children instead, so the count comes from
    enumeration and never from
    :func:`itacart.resolutions.refinement_ratio`.

    Args:
        index: Compositional index string.
        target_res: Resolution to descend to. Defaults to one level down.
        flatten: ``False``, the default, yields one list per cell of the
            input, preserving positional alignment. ``True`` yields a
            single flat stream. The shape does not depend on how many
            cells the input holds.

    Yields:
        Lists of child index strings, one list per cell of the input,
        when ``flatten`` is ``False``; child index strings when it is
        ``True``.

        The grouping does not collapse for a single cell, and that is
        worth stating because the neighbours do collapse:
        :func:`get_parent` answers a scalar and :func:`get_ancestors`
        answers a flat chain when handed one cell. Here
        ``list(get_children(cell))`` is a list holding one list, whose
        length is one and not the child count. Pass ``flatten=True`` to
        read the children directly.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
        ResolutionError: If ``target_res`` is not finer than the input.
        NonExistentCellError: If any cell of the index names no cell.
            The predicate is the arbiter and the contract ends there:
            a spelling it denies is refused rather than answered for
            the cell it would otherwise fold onto.
    """
    if target_res is not None and target_res > MAX_RESOLUTION:
        raise MaxResolutionError(
            f"resolution {target_res} is finer than the maximum {MAX_RESOLUTION}"
        )
    plans: list[tuple[str, int]] = []
    for cell in decompose(index):
        current = _resolution_of(split_components(cell))
        wanted = current + 1 if target_res is None else target_res
        if wanted <= current:
            raise ResolutionError(
                f"resolution {wanted} is not finer than {cell!r} "
                f"at resolution {current}"
            )
        plans.append((cell, wanted))
    if flatten:
        return _flat_children(plans)
    return _grouped_children(plans)


def _flat_children(plans: Iterable[tuple[str, int]]) -> Iterator[str]:
    for cell, wanted in plans:
        yield from _descendants_of(cell, wanted)


def _grouped_children(plans: Iterable[tuple[str, int]]) -> Iterator[list[str]]:
    for cell, wanted in plans:
        yield list(_descendants_of(cell, wanted))


def _descendants_of(cell: str, target_res: int) -> Iterator[str]:
    """Stream descendants of one atomic cell, depth first."""
    for child in _children_of(cell):
        if _resolution_of(split_components(child)) == target_res:
            yield child
        else:
            yield from _descendants_of(child, target_res)


def get_descendants(index: str, target_res: int) -> Iterator[str]:
    """Stream every descendant at a target resolution.

    Cardinality grows fast: from resolution 1 to 13 a base cell away from
    the border undergoes twelve refinements, six of ratio 4 and six of
    ratio 25, hence ``(4 * 25) ** 6`` — one million million cells. Always
    a generator, never a list.

    Args:
        index: Compositional index string.
        target_res: Resolution to expand to.

    Yields:
        Atomic index strings at ``target_res``.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
        NonExistentCellError: If any addressed cell names no cell.
        ResolutionError: If ``target_res`` is not finer than any input cell.
    """
    # Descent answers about a cell, so it is on the entering side of the
    # existence contract. The check is here rather than left to the
    # descent itself because the refusal has to happen at the boundary:
    # a spelling whose resolution already equals the target never reaches
    # the machinery that would have consulted the border.
    require_existing_cells(index)
    return cast("Iterator[str]", get_children(index, target_res, flatten=True))


def _parent_cell(cell: str) -> str:
    """The cell that actually fathers ``cell``, border cases included.

    The lexical prefix is the answer wherever it names a cell, and away
    from an absorbing parent it is the answer with no geometry at all:
    the alphabet defines the relation.

    Where the prefix names nothing, the child was spelled under a
    neighbour of its parent, and there are two ways that happens. On the
    lateral border the neighbour is the resolution-1 column immediately
    east, which never holds a cell of its own in that row; stepping one
    column west finds the parent. On the polar cap the neighbour is a
    *sibling code at an intermediate level* -- ``...(1(C2(1)))`` is a
    child of ``...(1(B2))`` -- and no column shift reaches that spelling
    from that child.

    The old rule knew the first and not the second. Both are tried here,
    and neither is trusted: a candidate is the parent only if the
    corrected child relation returns this cell. The relation is
    consumed, not restated.
    """
    from . import boundary

    components = split_components(cell)
    if _resolution_of(components) <= QUADRANT_RESOLUTION:
        raise MinResolutionError(f"{cell!r} is a quadrant and has no parent cell")

    # The cheap path is only as good as the predicate under it. It used
    # to accept a prefix whose effective ring folded -- ``...(3(A1))``
    # reported an area and named nothing anyone drew -- and answered with
    # it, which is how a child of the polar cap came to have two parents.
    # The predicate refuses such a spelling now, so the path can be
    # trusted; it is written this way to say that the trust is borrowed.
    prefix = join_components(components[:-1])
    if boundary.is_valid_cell(prefix):
        if not boundary.absorbs_border(prefix) or cell in _children_of(prefix):
            return prefix

    for candidate in _neighbouring_parents(components, prefix):
        if cell in _children_of(candidate):
            return candidate
    raise NonExistentCellError(f"{cell!r} has no parent cell in the domain")


def _neighbouring_parents(components: Sequence[str], prefix: str) -> Iterator[str]:
    """Cells at the prefix's own resolution that could father the cell.

    Two families, in the order they cost. The column one west of the
    prefix, which is the lateral border case and is a string edit. Then
    the siblings of the prefix, reached through the nearest ancestor that
    names a cell, which is the polar case and needs the child relation of
    that ancestor.

    Only a cell below resolution 1 gets here. The prefix of a
    resolution-1 cell is its quadrant, and every quadrant names a cell
    that absorbs nothing, so :func:`_parent_cell` answers with it before
    any neighbour is looked for.
    """
    from . import boundary

    column_text, _, row_text = components[1].partition(RES1_SEPARATOR)
    column, row = int(column_text), int(row_text)
    side = CELL_SIZE_M[BASE_CELL_RESOLUTION]
    assert side is not None
    last = boundary.last_lattice_column(components[0], row, side)
    # West as far as the absorbing column, not one step. A child can be
    # spelled two columns past its parent -- ``NE(0819/0747(1))`` belongs
    # to ``NE(0817/0747)`` -- and a single step lands on a column that names
    # no cell either.
    for target in range(column - 1, last - 1, -1):
        candidate = _shift_column(prefix, target - column)
        if boundary.is_valid_cell(candidate):
            yield candidate

    wanted = len(components) - 1
    for depth in range(wanted - 1, 1, -1):
        ancestor = join_components(components[:depth])
        if not boundary.is_valid_cell(ancestor):
            continue
        # The nearest ancestor that names a cell can sit several levels
        # above the prefix, because an intermediate spelling names nothing
        # either. Descend back down through the corrected relation until
        # the prefix's own resolution is reached; those cells are the
        # siblings the prefix should have had.
        frontier = [ancestor]
        for _ in range(depth, wanted):
            frontier = [child for cell in frontier for child in _children_of(cell)]
        for sibling in frontier:
            if sibling != prefix:
                yield sibling
        return


def child_position(cell: str) -> int | list[int]:
    """Ordinal of a cell among its siblings.

    Zero-based against the refinement alphabet of the cell's resolution:
    ``1``-``4`` for even levels, ``A1``-``E5`` in row-major order for odd.

    Under a border-absorbing parent the alphabet rank is not an ordinal,
    because two siblings can carry the same code under different column
    prefixes. There the ordinal is the position in the enumeration of
    :func:`get_children`: own stem before eastern stem, alphabet order
    within each. Away from the border the two definitions agree, and the
    cheap one is used.

    Args:
        cell: Compositional index string.

    Returns:
        An ordinal for a single cell, or a positionally aligned list.

    Raises:
        MinResolutionError: If any addressed cell is a whole quadrant.
        NonExistentCellError: If any addressed cell names no cell, or has
            no parent cell.
        ResolutionError: If any addressed cell is a resolution-1 cell.
    """
    from . import boundary

    # A sibling ordinal is an answer about a cell. Without this the
    # resolution complaint arrived first for a denied resolution-1
    # spelling, which reads as a refusal but answers a different
    # question: the quadrant case below still raises its own error,
    # because a quadrant is a cell the predicate accepts.
    require_existing_cells(cell)

    positions: list[int] = []
    for atom in decompose(cell):
        components = split_components(atom)
        current = _resolution_of(components)
        if current <= QUADRANT_RESOLUTION:
            raise MinResolutionError(f"{atom!r} is a quadrant and has no siblings")
        if current == BASE_CELL_RESOLUTION:
            raise ResolutionError(
                f"{atom!r} is a resolution-1 cell; its siblings are the "
                "lattice of its quadrant, which carries no refinement alphabet"
            )
        prefix = join_components(components[:-1])
        if boundary.is_valid_cell(prefix) and not boundary.absorbs_border(prefix):
            positions.append(refinement_alphabet(current).index(components[-1]))
            continue
        positions.append(_children_of(_parent_cell(atom)).index(atom))
    if len(positions) == 1 and _is_single(cell):
        return positions[0]
    return positions


# --------------------------------------------------------------------------
# Containment
# --------------------------------------------------------------------------


def is_ancestor(parent_index: str, child_index: str) -> bool:
    """Whether one cell contains another in the hierarchy.

    Strict, so a cell is not its own ancestor; the reflexive relation is
    :func:`contains`. A predicate over index strings that never consults
    geometry: the union of a cell's children was measured to overrun the
    cell itself by about one part in a thousand, because the domain
    border is curved and each cell approximates it by a chord, so a
    geometric test would disagree at the border with the index that
    addresses the cell.

    Args:
        parent_index: Candidate ancestor, a single atomic index.
        child_index: Candidate descendant, a single atomic index.

    Returns:
        ``True`` if ``parent_index`` is a strict prefix of
        ``child_index``.
    """
    ancestor = split_components(parent_index)
    descendant = split_components(child_index)
    if len(ancestor) >= len(descendant):
        return False
    return descendant[: len(ancestor)] == ancestor


def contains(region: str, cell: str) -> bool | list[bool]:
    """Whether a compositional region covers a cell.

    Complements :func:`is_ancestor`, which only tests the vertical
    relation between two single cells. A region covers a cell when any of
    its terminal cells is that cell or an ancestor of it.

    Like :func:`is_ancestor`, a predicate over index strings only.

    Args:
        region: Compositional index string denoting the region.
        cell: Compositional index string to test.

    Returns:
        A boolean for a single test cell, or a positionally aligned list.
    """
    terminals = [tuple(components) for components in _paths(region)]
    answers: list[bool] = []
    for candidate in _paths(cell):
        answers.append(
            any(
                len(terminal) <= len(candidate)
                and tuple(candidate[: len(terminal)]) == terminal
                for terminal in terminals
            )
        )
    if len(answers) == 1 and _is_single(cell):
        return answers[0]
    return answers


# --------------------------------------------------------------------------
# Compaction
# --------------------------------------------------------------------------


def _is_complete(parent: str, present: set[tuple[str, ...]]) -> bool:
    """Whether ``present`` holds every child of ``parent``.

    No guard against a childless or over-deep parent: candidates come
    from :func:`_candidate_parents`, which only proposes the prefix of a
    cell, so the parent sits between resolution 1 and 12 and always has
    children. A defective ring still propagates, deliberately.
    """
    children = _children_of(parent)
    return all(tuple(split_components(child)) in present for child in children)


def _candidate_parents(present: set[tuple[str, ...]]) -> set[str]:
    """Cells that could absorb one of the paths in ``present``.

    Asks :func:`_parent_cell` rather than editing the string, for the
    reason given there: under an absorbing parent a child can be spelled
    under a sibling code at an intermediate level, which no column shift
    reaches. A path whose parent cannot be resolved proposes nothing,
    since a set that holds it was never going to collapse.
    """
    candidates: set[str] = set()
    for components in present:
        if _resolution_of(components) <= BASE_CELL_RESOLUTION:
            continue
        try:
            candidates.add(_parent_cell(join_components(components)))
        except (NonExistentCellError, MinResolutionError):
            continue
    return candidates


def compact_cells(index: str) -> str:
    """Replace exhaustive sibling sets by their parent, recursively.

    Runs to a fixed point, so a fully covered base cell collapses all the
    way to resolution 1.

    Distinct from :func:`itacart.index.normalize`. Normalization rewrites
    spelling and never changes the region: it collapses a complete
    alphabet away from the border and leaves a border-absorbing node as
    written, because four quaternary children spell the whole alphabet
    while such a parent may have five. Compaction counts the children the
    parent actually has, so it folds an absorbing cell exactly when all of
    them are present, whichever stem they are spelled under, and leaves a
    partition that is incomplete at the border uncompacted rather than
    claimed whole.

    Args:
        index: Compositional index string.

    Returns:
        The compacted index, mixed-resolution in general, with cells in
        deterministic order.
    """
    present = {tuple(components) for components in _paths(index)}
    changed = True
    while changed:
        changed = False
        for parent in sorted(_candidate_parents(present), key=len, reverse=True):
            path = tuple(split_components(parent))
            if not _is_complete(parent, present):
                continue
            for child in _children_of(parent):
                present.discard(tuple(split_components(child)))
            present.add(path)
            changed = True
    ordered = sorted(present, key=_sort_key)
    return compose(join_components(components) for components in ordered)


def uncompact_cells(index: str, target_res: int) -> Iterator[str]:
    """Expand a compacted index back to uniform resolution.

    Args:
        index: Compositional index string, possibly mixed-resolution.
        target_res: Resolution to expand every cell to.

    Yields:
        Atomic index strings at ``target_res``.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
        NonExistentCellError: If any cell of the index names no cell.
        ResolutionError: If any terminal cell is finer than ``target_res``.
    """
    if target_res > MAX_RESOLUTION:
        raise MaxResolutionError(
            f"resolution {target_res} is finer than the maximum {MAX_RESOLUTION}"
        )
    # A cell already at the target resolution is yielded straight back by
    # the expansion below, which is how a denied spelling used to leave
    # through this name without anything having asked whether it named a
    # cell. Refusing where the index enters covers both branches.
    require_existing_cells(index)
    plans: list[tuple[str, bool]] = []
    for cell in decompose(index):
        current = _resolution_of(split_components(cell))
        if current > target_res:
            raise ResolutionError(
                f"{cell!r} is at resolution {current}, finer than the "
                f"requested {target_res}"
            )
        plans.append((cell, current == target_res))
    return _expand(plans, target_res)


def _expand(plans: Iterable[tuple[str, bool]], target_res: int) -> Iterator[str]:
    for cell, already_there in plans:
        if already_there:
            yield cell
        else:
            yield from _descendants_of(cell, target_res)
