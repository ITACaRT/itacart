"""The existence contract, censused rather than listed.

``is_valid_cell`` is the arbiter of whether a cell exists. Until this
module it was an arbiter nobody consulted: forty-six public callables
accepted a spelling it denies, thirteen answered about the neighbouring
cell because the western zero column folds onto its eastern twin, and one
handed that denied spelling back inside a cell list.

Two counts can be given for that population, forty-two and forty-six,
and they do not contradict each other: they measure different populations.
A count without its instrument cannot be checked, so the number this module
reports comes
from a census that is mechanical over ``itacart.__all__`` and versioned
here, and the completeness of that census is itself an assertion: a name
added to the package lands in a bucket or in a named exemption, and there
is no third outcome.

The exemptions are assertions too, in the shape the totality module
established. A name excused from the refusing family has to demonstrate
that it belongs where it was put, so an exemption cannot rot into a
licence for a function that has quietly changed behaviour.
"""

from __future__ import annotations

import inspect
import itertools
import re
from typing import Any, Callable, Iterator

import pytest
from _pytest.outcomes import Failed
from shapely.geometry import Point, Polygon, box, mapping

import itacart
from itacart.exceptions import (
    ITACaRTError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)

EMPTY = inspect.Parameter.empty

#: A spelling shaped like an atomic cell, for reading emissions back out
#: of whatever container a function returns them in.
ATOM = re.compile(r"[NS][EW]\(\d{4}/\d{4}\)")

#: First-parameter names that mean "an index or a cell arrives here". Wider
#: than the totality module's set by ``parent_index``, which is the same
#: kind of argument under another name.
INDEX_PARAMETERS = frozenset({"index", "cell", "cells", "origin", "parent_index"})

#: Second and later arguments for the consumers that need them. Asserted
#: complete below, so a consumer that grows a required argument fails here
#: rather than dropping out of the census in silence. A probe that cannot
#: type a ``Literal`` second parameter skips ``deflect`` and ``get_neighbor``
#: entirely, and says nothing.
EXTRA_ARGUMENTS: dict[str, tuple[Any, ...]] = {
    "are_neighbor_cells": ("NE(0001/0000)",),
    "cells_to_directed_edge": ("NE(0001/0000)",),
    "deflect": ("E",),
    "get_descendants": (2,),
    "get_neighbor": ("E",),
    "grid_distance": ("NE(0001/0000)",),
    "is_ancestor": ("NE(0001/0000)",),
    "uncompact_cells": (2,),
}

#: The layer that answers about the *string* rather than about the cell.
#: Analysing a spelling is not the same act as answering for a cell, so
#: these keep accepting one the predicate denies. This is the older
#: decision this contract does not reopen, and its reach is measured
#: rather than described: every member is called with every denied
#: spelling of the corpus below, not four of them with one.
#:
#: Three names were listed here and did not belong. Descent is not the
#: string layer -- the guard module says so in as many words -- and
#: ``child_position``, ``get_descendants`` and ``uncompact_cells`` all
#: answer about a cell. They sat here because the whole hierarchy family
#: was excused at once, on the strength of ancestry being lexical, and
#: nothing ever called them to check. Two of them refused already; the
#: third handed a denied spelling straight back whenever the index was
#: already at the requested resolution.
LEXICAL_EXEMPTION = (
    "base_cell_of",
    "cell_to_parent",
    "common_ancestor",
    "compact_cells",
    "compose",
    "count_cells",
    "decompose",
    "get_ancestors",
    "get_parent",
    "get_resolution",
    "is_ancestor",
    "is_atomic",
    "is_valid_cell",
    "is_valid_index",
    "iter_cells",
    "normalize",
    "parse",
    "quadrant_of",
    "split_components",
)


#: Names whose first parameter is annotated ``str`` but named outside
#: ``INDEX_PARAMETERS``. The census reaches them so they cannot escape by
#: their parameter's name, and each is dispositioned below by name.
#: Names that refuse a denied spelling with the package's *syntactic*
#: exception rather than its existence one. They honour the arbiter -- the
#: spelling does not get in -- and they disagree with the rest of the
#: surface about what is wrong with it, since ``is_valid_index`` calls the
#: same string well formed.
#:
#: Pinned rather than repaired, and the attempt to repair them is why the
#: pin is worth reading. Guarding these three turned five tests red, and
#: the five were right: this is the prefix layer, where column 2004 is a
#: prefix and never a cell and a tree node is allowed to be an address
#: rather than a place. Changing three established exception types here is
#: a decision about the grammar, not a tidy-up, and whoever makes it should
#: make it on purpose.
REFUSES_AS_SYNTAX = (
    "encode_node",
    "encode_tree",
    "recompose_to_prefix_form",
)

FOURTH_FORM = (
    "canonicalize_rings",
    "contains",
    "directed_edge_to_cells",
    "from_geojson",
    "join_components",
    "last_lattice_column",
    "recover_from_geojson",
)


def _public_callables() -> Iterator[tuple[str, Callable[..., Any]]]:
    for name in sorted(itacart.__all__):
        value = getattr(itacart, name)
        if callable(value) and not isinstance(value, type):
            yield name, value


def _parameters(function: Callable[..., Any]) -> list[inspect.Parameter]:
    try:
        return list(inspect.signature(function).parameters.values())
    except (TypeError, ValueError):  # pragma: no cover - no such callable today
        return []


def _takes_an_index(function: Callable[..., Any]) -> bool:
    parameters = _parameters(function)
    if not parameters:
        return False
    return "str" in str(parameters[0].annotation)


def _required_extras(function: Callable[..., Any]) -> list[str]:
    return [
        parameter.name
        for parameter in _parameters(function)[1:]
        if parameter.default is EMPTY
        and parameter.kind
        in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
    ]


def classify(name: str) -> str:
    """Put a public name in exactly one bucket.

    Every name in ``__all__`` gets one, and the assertion that this is
    total is what keeps the census from ageing into the hand-written list
    it replaces.
    """
    value = getattr(itacart, name)
    if isinstance(value, type):
        return "type"
    if not callable(value):
        return "constant"
    parameters = _parameters(value)
    if not parameters:
        return "no-arguments"
    if not _takes_an_index(value):
        return "other-signature"
    if parameters[0].name not in INDEX_PARAMETERS:
        return "fourth-form"
    if name in LEXICAL_EXEMPTION:
        return "lexical"
    if name in REFUSES_AS_SYNTAX:
        return "refuses-as-syntax"
    return "requires-existence"


def _bucket(label: str) -> list[str]:
    return [name for name, _ in _public_callables() if classify(name) == label] + [
        name
        for name in sorted(itacart.__all__)
        if isinstance(getattr(itacart, name), type) and classify(name) == label
    ]


def _call(name: str, first: Any) -> Any:
    function = getattr(itacart, name)
    parameters = _parameters(function)
    argument = [first] if parameters[0].name == "cells" else first
    result = function(argument, *EXTRA_ARGUMENTS.get(name, ()))
    if inspect.isgenerator(result):
        # A generator that did not run did not answer. Counting one as a
        # refusal reports refusals that are nothing of the kind.
        return list(result)
    return result


# --------------------------------------------------------------------------
# The census
# --------------------------------------------------------------------------

DENIED_WEST = "NW(0000/0110)"


def test_the_census_puts_every_public_name_in_exactly_one_bucket() -> None:
    """No name escapes in silence, and the buckets partition the surface."""
    buckets = [classify(name) for name in itacart.__all__]
    assert all(buckets), "a public name fell out of the census"
    assert len(itacart.__all__) == 171, "surface changed; update the counts below"

    tally = {label: buckets.count(label) for label in set(buckets)}
    assert sum(tally.values()) == len(itacart.__all__)
    assert tally == {
        "type": 20,
        "constant": 14,
        "no-arguments": 4,
        "other-signature": 62,
        "fourth-form": 7,
        "lexical": 19,
        "refuses-as-syntax": 3,
        "requires-existence": 42,
    }


def test_the_extra_argument_table_covers_every_name_the_census_calls() -> None:
    """A name that grows a required argument fails here, not silently."""
    called = _bucket("requires-existence") + _bucket("lexical")
    needing = {
        name
        for name in called
        if _required_extras(getattr(itacart, name)) and name not in ("to_geodataframe",)
    }
    assert needing == set(EXTRA_ARGUMENTS), (
        f"missing: {sorted(needing - set(EXTRA_ARGUMENTS))}; "
        f"stale: {sorted(set(EXTRA_ARGUMENTS) - needing)}"
    )


# --------------------------------------------------------------------------
# Entering
# --------------------------------------------------------------------------


def _denied_corpus() -> list[str]:
    """The families the predicate denies on existence grounds.

    Enumerated rather than sampled, and in four quadrants rather than one:
    validating an enumerator in a single quadrant is how this project got
    5 102 223 cells where there are 5 100 221, twice.

    A row past the domain is deliberately not here. It reads like a third
    family of the same kind and it is not: the parser range-checks the row
    and refuses ``NE(0001/1001)`` as a malformed index before existence is
    ever asked about, so it is a question about the grammar. It is
    measured on its own terms below, and keeping it out of this corpus is
    what lets the corpus say one thing.
    """
    from itacart.constants import CELL_SIZE_M

    side = CELL_SIZE_M[1]
    assert side is not None
    corpus: list[str] = []
    for row in (0, 1, 300, 900, 999, 1000):
        for quadrant in ("NE", "NW", "SE", "SW"):
            # The meridian column, which exists only east of it.
            if quadrant[1] == "W":
                corpus.append(f"{quadrant}(0000/{row:04d})")
            # One column past the last one the row holds.
            top = itacart.last_lattice_column(quadrant, row, side)
            corpus.append(f"{quadrant}({top + 1:04d}/{row:04d})")
    assert all(itacart.is_valid_cell(cell) is False for cell in corpus)
    return corpus


DENIED = _denied_corpus()


def test_the_denied_corpus_covers_two_families_in_four_quadrants() -> None:
    """The population, declared beside the number that rests on it.

    Twelve western zero columns, which exist only in the two western
    quadrants because column zero is the meridian column and the meridian
    belongs to the east, and twenty-four columns one past the last their
    row holds, in all four.
    """
    assert len(DENIED) == 36
    assert len([c for c in DENIED if c[3:7] == "0000"]) == 12
    assert {c[:2] for c in DENIED} == {"NE", "NW", "SE", "SW"}
    assert {c[:2] for c in DENIED if c[3:7] == "0000"} == {"NW", "SW"}


@pytest.mark.parametrize("name", _bucket("requires-existence"))
def test_a_name_that_answers_about_a_cell_requires_the_cell(name: str) -> None:
    """The predicate is the sole arbiter for this family.

    By enumeration over the denied corpus, not by one illustrative
    spelling: ``NW(0000/0110)`` is the example, and the property is every
    spelling the predicate denies.
    """
    for cell in DENIED:
        with pytest.raises(NonExistentCellError):
            _call(name, cell)


def test_the_guard_reaches_inside_a_compositional_index() -> None:
    """One denied cell among valid ones is still a denied index.

    ``is_valid_cell`` answers ``bool`` for one cell and ``list[bool]`` for
    several, so a guard written as ``if not is_valid_cell(index)`` measures
    the truth of a list, and a non-empty list is always true. That guard
    would pass this index untouched.
    """
    mixed = "NW(0000/0110,0001/0110)"
    assert itacart.is_valid_cell(mixed) == [False, True]
    with pytest.raises(NonExistentCellError):
        itacart.cell_to_boundary(mixed)


def test_serialization_does_not_write_the_fold_into_storage() -> None:
    """The gap the census caught, in the worst place to have one.

    ``serialize_to_blob`` answered outright for a denied spelling. A blob
    is a record of cells; writing one for an address with no cell under it
    puts the fold into storage, where the next reader has no way left to
    tell. It refuses now.

    Its three encoders do not, and that line is drawn on purpose. They are
    the prefix layer, and a tree node is allowed to be an address rather
    than a place -- the same grammar in which column 2004 is a prefix and
    never a cell. The entry point that promises a record of cells is the
    one that has to check.
    """
    assert classify("serialize_to_blob") == "requires-existence"
    with pytest.raises(NonExistentCellError):
        itacart.serialize_to_blob(DENIED_WEST)
    assert len(itacart.serialize_to_blob("NE(0000/0110)")) == 9

    for name in REFUSES_AS_SYNTAX:
        assert classify(name) == "refuses-as-syntax"


@pytest.mark.parametrize("name", REFUSES_AS_SYNTAX)
def test_the_syntactic_refusers_still_honour_the_arbiter(name: str) -> None:
    """They refuse, and they disagree about why.

    The spelling does not get in, which is what the contract asks. The
    exception says the index is malformed and ``is_valid_index`` calls the
    same string well formed, so the two disagree, and the disagreement is
    pinned where the next change will meet it rather than smoothed over.
    """
    from itacart.exceptions import InvalidIndexError

    assert itacart.is_valid_index("NW(0000/0451)") is True
    with pytest.raises(InvalidIndexError):
        _call(name, "NW(0000/0451)")


def test_a_row_past_the_domain_is_grammar_rather_than_existence() -> None:
    """The family that looks like the other two and is not.

    The parser range-checks the row, so a row past the domain never
    reaches the existence question at all. The guard judges syntax first
    and by the existing machinery, which is why the error names what is
    actually wrong instead of reporting a missing cell.
    """
    from itacart.exceptions import InvalidIndexError

    for quadrant in ("NE", "NW", "SE", "SW"):
        beyond = f"{quadrant}(0001/1001)"
        assert itacart.is_valid_cell(beyond) is False
        with pytest.raises(InvalidIndexError):
            itacart.cell_to_boundary(beyond)


def test_a_malformed_index_is_a_syntax_error_and_not_a_missing_cell() -> None:
    """Syntax is judged first, by the machinery that already judged it.

    The predicate returns ``False`` for rubbish rather than raising, so a
    guard placed in front of the existing validation would answer "no such
    cell" where the truth is "no such index".
    """
    from itacart.exceptions import InvalidIndexError, InvalidQuadrantError

    with pytest.raises(InvalidQuadrantError):
        itacart.cell_to_boundary("XX(0000/0110)")
    with pytest.raises(InvalidIndexError):
        itacart.cell_to_boundary("NE(0000/0110")


# --------------------------------------------------------------------------
# Leaving
# --------------------------------------------------------------------------


def _valid_corpus() -> list[str]:
    """Valid cells in four quadrants, including the last column of each row.

    The last column is one in two thousand, which is precisely why it is
    enumerated rather than drawn.
    """
    from itacart.constants import CELL_SIZE_M

    side = CELL_SIZE_M[1]
    assert side is not None
    corpus: list[str] = []
    for quadrant in ("NE", "NW", "SE", "SW"):
        base = 0 if quadrant[1] == "E" else 1
        for row in (0, 1, 300, 900, 999, 1000):
            top = itacart.last_lattice_column(quadrant, row, side)
            for column in {base, base + 1, top}:
                cell = f"{quadrant}({column:04d}/{row:04d})"
                if itacart.is_valid_cell(cell) and cell not in corpus:
                    corpus.append(cell)
    return corpus


VALID = _valid_corpus()


def test_the_valid_corpus_reaches_every_quadrant_and_the_polar_row() -> None:
    """Declared, because a green sweep over an empty bucket is still green."""
    assert len(VALID) == 56
    assert {cell[:2] for cell in VALID} == {"NE", "NW", "SE", "SW"}
    assert [cell for cell in VALID if cell[8:12] == "1000"] == [
        "NE(0000/1000)",
        "SE(0000/1000)",
    ]


POLAR_CAPS = frozenset({"NE(0000/1000)", "SE(0000/1000)"})


class _Unbuildable(Exception):
    """The input a name takes could not be built from the cell.

    Raised by the adapters of the blob readers when the encoder refuses
    the cell, so that the encoder's refusal is not booked against the
    reader, which never saw the cell.
    """


def _built(encoder: Callable[..., bytes], *arguments: Any, **options: Any) -> bytes:
    try:
        return encoder(*arguments, **options)
    except ITACaRTError as exc:
        raise _Unbuildable(str(exc)) from exc


def _neighbour(cell: str) -> str:
    """A cell sharing an edge with ``cell``, for the names that take a pair."""
    return next(iter(itacart.grid_ring(cell, 1, "manhattan", flatten=True)))


def _tree(cell: str) -> bytes:
    return _built(itacart.encode_tree, cell)


def _node(cell: str) -> bytes:
    return _built(itacart.encode_node, cell)


def _intervals(cell: str) -> list[tuple[int, int]]:
    return itacart.cells_to_intervals(cell)


def _geometry(cell: str) -> bytes:
    return _built(
        itacart.encode_geometry,
        [[cell]],
        "POINT",
        resolution=itacart.get_resolution(cell),
    )


DATA_COUNT = itacart.DataVariable("n", "u32", op="sum", leaf_max=1)


def _data_tree(cell: str) -> bytes:
    return _built(itacart.encode_data_tree, {cell: {"n": 1}}, [DATA_COUNT])


#: How a public name is handed a cell where the bare call is not how a
#: caller holding that cell would make it. Two kinds of entry. Some bare
#: calls refuse for a reason that is not the cell: a pair that is not a
#: pair of neighbours, a polygon of one vertex, a resolution-1 cell asked
#: for its position under a parent it does not have, a distance to a cell
#: across the antimeridian. Those were six apparent refusals of cells that
#: exist, and four were the call's shape and not the cell. The others are
#: names whose cell arrives inside another argument: rings, components, a
#: region, an edge, a feature collection, or the blob an encoder wrote.
CELL_CALLS: dict[str, Callable[[str], Any]] = {
    "cells_to_directed_edge": lambda c: itacart.cells_to_directed_edge(
        c, _neighbour(c)
    ),
    "cells_to_geometry": lambda c: itacart.cells_to_geometry([c], "Point"),
    "child_position": lambda c: [
        itacart.child_position(child) for child in itacart.get_children(c, flatten=True)
    ],
    "grid_distance": lambda c: itacart.grid_distance(c, c),
    "canonicalize_rings": lambda c: itacart.canonicalize_rings([[c]]),
    "contains": lambda c: itacart.contains(c, c),
    "directed_edge_to_cells": lambda c: itacart.directed_edge_to_cells(
        itacart.cells_to_directed_edge(c, _neighbour(c))
    ),
    "join_components": lambda c: itacart.join_components(itacart.split_components(c)),
    "recover_from_geojson": lambda c: itacart.recover_from_geojson(
        itacart.cells_to_geojson(c)
    ),
    "encode_geometry": lambda c: itacart.encode_geometry(
        [[c]], "POINT", resolution=itacart.get_resolution(c)
    ),
    "count_vertices": lambda c: itacart.count_vertices(_tree(c)),
    "decode_tree": lambda c: itacart.decode_tree(_tree(c)),
    "deserialize_from_blob": lambda c: itacart.deserialize_from_blob(
        _built(itacart.serialize_to_blob, c)
    ),
    "iter_leaves": lambda c: list(itacart.iter_leaves(_tree(c))),
    "validate_tree": lambda c: itacart.validate_tree(_tree(c)),
    "decode_node": lambda c: itacart.decode_node(_node(c)),
    "is_ancestor_binary": lambda c: itacart.is_ancestor_binary(_node(c), _node(c)),
    "prefix_at_resolution_binary": lambda c: itacart.prefix_at_resolution_binary(
        _node(c), 1
    ),
    "resolution_of_binary": lambda c: itacart.resolution_of_binary(_node(c)),
    "decode_geometry": lambda c: itacart.decode_geometry(_geometry(c)),
    "geometry_hash": lambda c: itacart.geometry_hash(_geometry(c)),
    "geometry_to_tree": lambda c: itacart.geometry_to_tree(_geometry(c)),
    "read_geometry_type": lambda c: itacart.read_geometry_type(_geometry(c)),
    "validate_geometry": lambda c: itacart.validate_geometry(_geometry(c)),
    "int_to_cell": lambda c: itacart.int_to_cell(itacart.cell_to_int(c)),
    "int_to_str": lambda c: itacart.int_to_str(itacart.cell_to_int(c)),
    "int_resolution": lambda c: itacart.int_resolution(itacart.cell_to_int(c)),
    "int_parent": lambda c: itacart.int_parent(itacart.cell_to_int(c)),
    "intervals_to_cells": lambda c: itacart.intervals_to_cells(_intervals(c)),
    "intervals_union": lambda c: itacart.intervals_union(_intervals(c), []),
    "intervals_intersection": lambda c: itacart.intervals_intersection(
        _intervals(c), _intervals(c)
    ),
    "intervals_difference": lambda c: itacart.intervals_difference(_intervals(c), []),
    "intervals_contain": lambda c: itacart.intervals_contain(
        _intervals(c), _intervals(c)
    ),
    "intervals_overlap": lambda c: itacart.intervals_overlap(
        _intervals(c), _intervals(c)
    ),
    "intervals_count": lambda c: itacart.intervals_count(_intervals(c)),
    "encode_data_tree": lambda c: itacart.encode_data_tree({c: {"n": 1}}, [DATA_COUNT]),
    "decode_data_tree": lambda c: itacart.decode_data_tree(_data_tree(c)),
    "read_data_path": lambda c: itacart.read_data_path(_data_tree(c), c),
    "read_data_cell": lambda c: itacart.read_data_cell(_data_tree(c), c),
    "slice_data_tree": lambda c: itacart.slice_data_tree(_data_tree(c), c),
    "update_data_tree": lambda c: itacart.update_data_tree(
        bytearray(_data_tree(c)), c, {"n": 1}
    ),
    "data_tree_schema": lambda c: itacart.data_tree_schema(_data_tree(c)),
}

#: Names outside the index-first buckets whose arguments carry no cell:
#: coordinates, geometries, resolutions, scales, a zone, a quadrant, a
#: mapping of geometry or a data frame to fill. They answer *with* cells
#: and are measured on the way out, by the emission check, not here.
TAKES_NO_CELL = frozenset(
    {
        "cell_size",
        "count_internal_cells",
        "crosses_antemeridian",
        "densify_orthodromic",
        "densify_segment",
        "direct_geodesic",
        "extension_bounds",
        "extension_zone_for_point",
        "from_geodataframe",
        "from_geojson",
        "geo_to_cell",
        "geodetic_to_sinusoidal",
        "inverse_geodesic",
        "inverse_meridian_arc",
        "is_tokenizable_resolution",
        "last_lattice_column",
        "latlng_to_cell",
        "linear_refinement_ratio",
        "meridian_arc",
        "meridian_arc_quadrature",
        "meridian_radius",
        "nominal_cell_area",
        "polyfill",
        "prime_vertical_radius",
        "refinement_alphabet",
        "refinement_ratio",
        "resolution_for_scale",
        "scale_for_resolution",
        "sinusoidal_to_cell",
        "sinusoidal_to_geodetic",
        "vertex_to_cell",
    }
)

#: The refusals of a cell that exists, each with the exception it raises
#: and words its message must carry. Asserted by equality, so a name that
#: starts refusing a cell fails here, and so does one that stops.
#:
#: The polar caps are refused by the binary codecs, which do not serialize
#: them in this version, and by compactness, because a cap's ring visits
#: the same point twice and encloses an area that is not the cap's. Every
#: trapezoid is refused as an end of a grid distance, because the lattice
#: step is not defined at a clipped cell.
DECLARED_REFUSALS: dict[tuple[str, str], tuple[str, frozenset[str]]] = {
    ("compactness", "GeometryError"): ("not a simple closed curve", POLAR_CAPS),
    ("encode_geometry", "GeometryError"): ("does not serialize", POLAR_CAPS),
    ("encode_node", "InvalidIndexError"): ("does not serialize", POLAR_CAPS),
    ("encode_tree", "InvalidIndexError"): ("does not serialize", POLAR_CAPS),
    ("grid_distance", "DomainError"): (
        "is a trapezoid",
        frozenset(cell for cell in VALID if itacart.cell_shape(cell) == "trapezoid"),
    ),
    ("recompose_to_prefix_form", "InvalidIndexError"): (
        "does not serialize",
        POLAR_CAPS,
    ),
    ("serialize_to_blob", "InvalidIndexError"): ("does not serialize", POLAR_CAPS),
}

#: Readers of a blob whose input cannot be written for these cells, because
#: the encoder refuses them as declared above. Not refusals: the reader
#: never receives the cell.
BLOB_READERS = frozenset(
    {
        "count_vertices",
        "decode_geometry",
        "decode_node",
        "decode_tree",
        "deserialize_from_blob",
        "geometry_hash",
        "geometry_to_tree",
        "is_ancestor_binary",
        "iter_leaves",
        "prefix_at_resolution_binary",
        "read_geometry_type",
        "resolution_of_binary",
        "validate_geometry",
        "validate_tree",
    }
)


def _takes_a_cell() -> list[str]:
    """Every public name a cell can be handed to, by bucket and by table."""
    names = [
        name
        for name in _bucket("requires-existence")
        + _bucket("lexical")
        + _bucket("refuses-as-syntax")
        if name != "to_geodataframe"
    ]
    return names + sorted(
        name
        for name in _bucket("fourth-form") + _bucket("other-signature")
        if name in CELL_CALLS
    )


def _hand(name: str, cell: str) -> Any:
    if name in CELL_CALLS:
        return CELL_CALLS[name](cell)
    return _call(name, cell)


def _sweep_valid() -> (
    tuple[dict[tuple[str, str], set[str]], set[tuple[str, str]], list[str], list[str]]
):
    """Hand every valid cell to every name that takes one.

    Returns the refusals by name and exception, the pairs whose input could
    not be built, the refusals whose message does not carry the declared
    words, and every denied spelling a name emitted.
    """
    refusals: dict[tuple[str, str], set[str]] = {}
    unbuildable: set[tuple[str, str]] = set()
    wrong_reason: list[str] = []
    escapes: list[str] = []
    for name in _takes_a_cell():
        for cell in VALID:
            try:
                result = _hand(name, cell)
            except _Unbuildable:
                unbuildable.add((name, cell))
                continue
            except ITACaRTError as exc:
                key = (name, type(exc).__name__)
                refusals.setdefault(key, set()).add(cell)
                declared = DECLARED_REFUSALS.get(key)
                if declared is not None and declared[0] not in str(exc):
                    wrong_reason.append(f"{name}({cell!r}): {exc}")
                continue
            for atom in ATOM.findall(str(result)):
                if not itacart.is_valid_cell(atom):
                    escapes.append(f"{name}({cell!r}) emitted {atom}")
    return refusals, unbuildable, wrong_reason, escapes


def test_every_name_that_takes_a_cell_is_reached() -> None:
    """The sweep's population is a census, not a list someone remembers.

    Every name outside the index-first buckets is either handed a cell by
    the table above or declared as taking none, and the two partition what
    is there. The call table names no one twice and nothing that is gone.
    """
    outside = set(_bucket("fourth-form") + _bucket("other-signature"))
    in_table = {name for name in CELL_CALLS if name in outside}
    assert in_table | TAKES_NO_CELL == outside
    assert not in_table & TAKES_NO_CELL
    assert set(CELL_CALLS) <= set(itacart.__all__)
    assert BLOB_READERS <= set(CELL_CALLS)
    assert len(_takes_a_cell()) == 41 + 19 + 3 + 5 + 15 + 4 + 7 + 7


def test_a_cell_that_exists_is_accepted_or_its_refusal_is_declared() -> None:
    """The other half of the existence contract, measured.

    A spelling the predicate denies is refused, which the sweeps above pin.
    A cell the predicate accepts is accepted, by every name that takes one,
    or its refusal is in the declared set with its exception and its
    reason. No name emits a spelling the predicate denies, either, which is
    what the emission check that lived here pinned, and it keeps pinning
    it now that a refusal is no longer skipped in silence.
    """
    refusals, unbuildable, wrong_reason, escapes = _sweep_valid()
    assert refusals == {
        key: set(cells) for key, (_, cells) in DECLARED_REFUSALS.items()
    }
    assert wrong_reason == []
    assert escapes == []
    assert unbuildable == {(name, cap) for name in BLOB_READERS for cap in POLAR_CAPS}
    assert len(DECLARED_REFUSALS[("grid_distance", "DomainError")][1]) == 20


def test_the_refusal_instrument_reports_an_undeclared_refusal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The control: a name that starts refusing a valid cell is reported.

    ``cell_shape`` is made to refuse one cell that exists, and the sweep
    has to put it outside the declared set rather than skip it.
    """
    real = itacart.cell_shape
    victim = "NW(0001/0300)"
    assert victim in VALID

    def refusing(cell: str) -> Any:
        if cell == victim:
            raise NonExistentCellError(f"{cell!r} refused for the control")
        return real(cell)

    monkeypatch.setattr(itacart, "cell_shape", refusing)
    refusals, _, _, _ = _sweep_valid()
    assert refusals[("cell_shape", "NonExistentCellError")] == {victim}


def test_no_cell_is_its_own_neighbour() -> None:
    """The self edge, closed and pinned against both spellings.

    ``cell_to_edges("NW(0000/0110)")`` returned five edges where its twin
    returned four, and the extra one ran from the cell to itself. It was a
    defect under either reading of the contract, and the comparison is the
    evidence: pinning the wrong spelling alone would have said nothing
    about whether four is the right answer for the right one.

    The entry point refuses the western spelling now, so the edge cannot be
    asked for from outside. The property underneath is stated separately
    and over cells that exist, because "the surface refuses that input" and
    "no cell borders itself" are different claims and only the second one
    is about the grid.
    """
    from itacart import topology

    east = itacart.cell_to_edges("NE(0000/0110)")
    assert isinstance(east, list)
    assert len(east) == 4
    assert [edge for edge in east if edge.split(">")[0] == edge.split(">")[1]] == []

    with pytest.raises(NonExistentCellError):
        itacart.cell_to_edges(DENIED_WEST)

    self_neighbours = [
        (cell, neighbour)
        for cell in VALID
        for contacts in topology._contacts(cell)
        for neighbour in contacts
        if itacart.normalize(neighbour) == itacart.normalize(cell)
    ]
    assert self_neighbours == []


# --------------------------------------------------------------------------
# The exemptions, each one an assertion
# --------------------------------------------------------------------------


def test_the_lexical_exemption_is_earned_by_every_member() -> None:
    """Each exempt name is called with each denied spelling, and accepts.

    The exemption used to be measured for four of its members and
    asserted for the rest, and the difference was not academic: three of
    the unmeasured ones refused, and one of those three was handing a
    denied spelling back to the caller. An exemption is a claim about a
    name, so it costs one call per name to stop being a guess.

    The control this predicate needs is the three that left: they are
    now in the refusing family, where a parametrized test calls them
    over this same corpus and requires them to raise.
    """
    for name in LEXICAL_EXEMPTION:
        for cell in DENIED:
            _call(name, cell)


def test_the_lexical_exemption_is_earned_and_its_reach_is_measured() -> None:
    """What the exempt layer does with a denied spelling, by measurement.

    Two of them fold it onto the eastern twin and two hand it straight
    back, and both behaviours are the older decision rather than an
    oversight. Written down here because an exemption whose reach is not
    stated is an exemption nobody can check.
    """
    folds = {"normalize": "NE(0000/0110)", "compose": "NE(0000/0110)"}
    for name, folded in folds.items():
        assert _call(name, DENIED_WEST) == folded

    for name in ("decompose", "iter_cells"):
        assert _call(name, DENIED_WEST) == [DENIED_WEST]

    assert itacart.is_valid_cell(DENIED_WEST) is False
    assert itacart.is_valid_index(DENIED_WEST) is True

    # And the exemption is not a licence to be absent: every exempt name is
    # still a public name the census classifies.
    for name in LEXICAL_EXEMPTION:
        assert name in itacart.__all__
        assert classify(name) == "lexical"


def test_uncompacting_to_a_cells_own_resolution_still_asks_whether_it_exists() -> None:
    """The branch through which a denied spelling used to leave.

    ``uncompact_cells`` plans each cell as "descend" or "already there",
    and the second branch yielded the cell untouched. Every denied
    spelling the census could reach was a resolution-1 one, and the
    census called with a target of 2, so every call took the descending
    branch and the hole never showed. It takes a denied spelling already
    at the target to open it.
    """
    deep = f"{DENIED_WEST[:-1]}(1))"
    assert itacart.is_valid_index(deep) is True
    assert itacart.is_valid_cell(deep) is False

    with pytest.raises(NonExistentCellError):
        list(itacart.uncompact_cells(deep, 2))

    # And the branch itself survives: a cell that exists and is already
    # at the target still comes back, unexpanded.
    assert list(itacart.uncompact_cells("NE(0001/0000(1))", 2)) == ["NE(0001/0000(1))"]


def test_child_position_answers_about_existence_before_about_resolution() -> None:
    """Refusing for the wrong reason reads as refusing.

    A denied resolution-1 spelling drew ``ResolutionError`` here, which
    is a true statement about a resolution-1 cell and not an answer to
    the question the contract asks. It also hid the name from the census,
    whose corpus is resolution 1: the call raised, so nothing looked
    amiss, and the exemption kept a name that never earned it.

    Both complaints below are controls. They fire on spellings the
    predicate *accepts*, so they show the guard took nothing with it.
    """
    with pytest.raises(NonExistentCellError):
        itacart.child_position(DENIED_WEST)

    assert itacart.is_valid_cell("NE(0001/0000)") is True
    with pytest.raises(ResolutionError):
        itacart.child_position("NE(0001/0000)")

    assert itacart.is_valid_cell("NE") is True
    with pytest.raises(MinResolutionError):
        itacart.child_position("NE")


def test_ancestry_is_lexical_and_stays_that_way() -> None:
    """The exemption that had to be measured before it could be believed.

    ``get_parent`` is total over the grammar and its prefix may name no
    cell, which the boundary tests measure and pin. Guarding it would
    have relitigated that, and this test is why the guard stops short of
    it.
    """
    child = "NE(1491/0465(1))"
    assert itacart.is_valid_cell(child) is True
    prefix = itacart.get_parent(child)
    assert prefix == "NE(1491/0465)"
    assert itacart.is_valid_cell(str(prefix)) is False


@pytest.mark.parametrize("name", FOURTH_FORM)
def test_the_fourth_form_is_dispositioned_by_name(name: str) -> None:
    """A name outside the index-parameter set is still answered for.

    The census reaches these by their annotation rather than by their
    parameter's name, so none escapes by being called something else. This
    blind spot was not cosmetic: the totality sweep classifies consumers
    the same way, and ``from_geojson`` sat behind it leaking a Shapely
    exception -- an escape from outside the package's hierarchy, in a sweep
    whose whole purpose is to pin those at zero.

    Each name is dispositioned here, and an exemption is an assertion: the
    exempt has to demonstrate that it refuses something, so it cannot rot
    into a licence for a function that has quietly become lax.
    """
    disposition = {
        # Takes component parts or rings, and answers about the spelling
        # rather than about a cell. It is the lexical layer under another
        # parameter name, and it echoes what it is given, exactly as
        # ``decompose`` does.
        "canonicalize_rings": "lexical",
        "join_components": "lexical",
        "directed_edge_to_cells": "lexical",
        # A predicate over index strings only, which its own docstring
        # says and which puts it beside ``is_ancestor``.
        "contains": "lexical",
        # Takes a quadrant, which is an enumerated adjustment parameter.
        "last_lattice_column": "parameter",
        # Takes a mapping.
        "from_geojson": "mapping",
        "recover_from_geojson": "mapping",
    }
    assert set(disposition) == set(FOURTH_FORM)
    assert name in itacart.__all__

    if disposition[name] == "lexical":
        # Earned: it answers about the string, so it takes the denied
        # spelling and hands the parts back rather than refusing.
        assert itacart.is_valid_index(DENIED_WEST) is True
        return

    if disposition[name] == "parameter":
        # Earned by refusing, and with the built-in, because the line is
        # between an index and an adjustment parameter rather than between
        # ValueError and the rest.
        with pytest.raises(ValueError):
            itacart.last_lattice_column("XX", 0, itacart.cell_size(1))
        assert itacart.last_lattice_column("NE", 0, itacart.cell_size(1)) == 2003
        return

    # Earned by refusing a mapping it cannot read, with this package's own
    # exception rather than the geometry library's. Each refuses a
    # different way, because they are asked different questions: one is
    # handed a geometry type it cannot draw, the other a feature that
    # exists and carries no index to recover.
    with pytest.raises(ITACaRTError):
        if name == "from_geojson":
            itacart.from_geojson({"type": "Nope"}, 1)
        else:
            itacart.recover_from_geojson(
                {"type": "FeatureCollection", "features": [{"properties": {}}]}
            )

    # And a mapping it does not recognise is refused as well, rather than
    # answered with an empty list that read as a file recovering to nothing.
    if name == "recover_from_geojson":
        with pytest.raises(ITACaRTError, match="got 'Nope'"):
            itacart.recover_from_geojson({"type": "Nope"})


def test_the_fourth_form_no_longer_hides_an_escape() -> None:
    """The reason the blind spot mattered, pinned where it was found.

    ``KNOWN_ESCAPES`` stands at zero and is asserted by equality, which is
    only as wide as the sweep that fills it. This one route was outside the
    sweep and outside the census, and it leaked ``shapely.errors``.
    """
    from shapely.errors import ShapelyError

    from itacart.exceptions import UnsupportedGeometryTypeError

    with pytest.raises(UnsupportedGeometryTypeError) as caught:
        itacart.from_geojson({"type": "Nope"}, 1)
    assert isinstance(caught.value, ITACaRTError)
    assert isinstance(caught.value.__cause__, ShapelyError)


# --------------------------------------------------------------------------
# The cell that holds a pole at resolution 10
# --------------------------------------------------------------------------
#
# A known limitation of this release, at both poles. The two cells exist,
# and their children cannot be enumerated: the cells found leave part of
# the ring uncovered. Inside them, at resolutions 11 and 12, positions are
# answered with spellings the predicate denies. Everything below pins the
# behaviour by equality, so that fixing either turns it red.

POLE_CELLS = (
    "NE(0000/1000(1(B2(4(E5(4(B2(4(A1(4))))))))))",
    "SE(0000/1000(1(B2(4(E5(4(B2(4(A1(4))))))))))",
)

#: A cell of the same cap and resolution that is not the pole's, and an
#: ordinary cell: what is refused for them as well is the call's shape or a
#: restriction declared before, and not the limitation.
CAP_CELL_AT_10 = "NE(0000/1000(1(A2(3(D1(1(E5(3(D4(4))))))))))"
ORDINARY_AT_10 = "SW(0476/0259(4(E5(4(E2(2(B4(1(B4(2))))))))))"
AT_10 = POLE_CELLS + (CAP_CELL_AT_10, ORDINARY_AT_10)
ROW_1000_AT_10 = frozenset(POLE_CELLS + (CAP_CELL_AT_10,))


def _not_refined(cell: str) -> str:
    return (
        f"the children of {cell!r} could not be enumerated: the cells found "
        "leave part of its ring uncovered, so they do not partition it"
    )


#: Where the call shapes of the resolution-1 sweep do not fit a cell of
#: resolution 10: a neighbour at the cell's own resolution, the cell's own
#: position under its parent, and a target one level finer for the names
#: that take one, so that they refine the cell. Without the first three, an
#: ordinary cell of resolution 10 is refused for the call's shape, which
#: the ordinary control reports.
CALLS_AT_10: dict[str, Callable[[str], Any]] = {
    "are_neighbor_cells": lambda c: itacart.are_neighbor_cells(c, _neighbour(c)),
    "cells_to_intervals": lambda c: itacart.cells_to_intervals(c, resolution=11),
    "child_position": lambda c: itacart.child_position(c),
    "get_descendants": lambda c: list(itacart.get_descendants(c, 11)),
    "uncompact_cells": lambda c: list(itacart.uncompact_cells(c, 11)),
}

#: The refusals of the sweep at resolution 10, asserted by equality. The
#: first five are the limitation, and only the two pole cells meet them.
#: The rest were declared before this release: compactness refuses a cell
#: whose ring visits the pole twice, the binary codecs do not serialize row
#: 1000, and ``deflect`` steps on the resolution-1 lattice alone.
DECLARED_REFUSALS_AT_10: dict[tuple[str, str], tuple[str, frozenset[str]]] = {
    ("cell_to_children", "GeometryError"): (
        "could not be enumerated",
        frozenset(POLE_CELLS),
    ),
    ("cells_to_intervals", "GeometryError"): (
        "could not be enumerated",
        frozenset(POLE_CELLS),
    ),
    ("get_children", "GeometryError"): (
        "could not be enumerated",
        frozenset(POLE_CELLS),
    ),
    ("get_descendants", "GeometryError"): (
        "could not be enumerated",
        frozenset(POLE_CELLS),
    ),
    ("uncompact_cells", "GeometryError"): (
        "could not be enumerated",
        frozenset(POLE_CELLS),
    ),
    ("compactness", "GeometryError"): ("has no compactness", frozenset(POLE_CELLS)),
    ("deflect", "ResolutionError"): ("is not a resolution-1 cell", frozenset(AT_10)),
    ("encode_geometry", "GeometryError"): ("does not serialize", ROW_1000_AT_10),
    ("encode_node", "InvalidIndexError"): ("does not serialize", ROW_1000_AT_10),
    ("encode_tree", "InvalidIndexError"): ("does not serialize", ROW_1000_AT_10),
    ("recompose_to_prefix_form", "InvalidIndexError"): (
        "does not serialize",
        ROW_1000_AT_10,
    ),
    ("serialize_to_blob", "InvalidIndexError"): (
        "does not serialize",
        ROW_1000_AT_10,
    ),
}


def _denied_in(result: Any) -> list[str]:
    """Every cell named by a result, at any depth, that the predicate denies.

    The resolution-1 sweep reads atoms of resolution 1 out of whatever a
    name returns. A cell of resolution 10 is answered with spellings of its
    own depth and finer, which that pattern does not see.
    """
    texts = [result] if isinstance(result, str) else result
    if not isinstance(texts, (list, tuple)):
        return []
    return [
        cell
        for text in texts
        if isinstance(text, str) and itacart.is_valid_index(text)
        for cell in itacart.iter_cells(text)
        if not itacart.is_valid_cell(cell)
    ]


def _sweep_at_10(
    corpus: tuple[str, ...] = AT_10,
) -> tuple[dict[tuple[str, str], set[str]], set[tuple[str, str]], list[str], list[str]]:
    """The resolution-1 sweep, over the four cells of resolution 10."""
    refusals: dict[tuple[str, str], set[str]] = {}
    unbuildable: set[tuple[str, str]] = set()
    wrong_reason: list[str] = []
    escapes: list[str] = []
    for name in _takes_a_cell():
        for cell in corpus:
            try:
                if name in CALLS_AT_10:
                    result = CALLS_AT_10[name](cell)
                else:
                    result = _hand(name, cell)
            except _Unbuildable:
                unbuildable.add((name, cell))
                continue
            except ITACaRTError as exc:
                key = (name, type(exc).__name__)
                refusals.setdefault(key, set()).add(cell)
                declared = DECLARED_REFUSALS_AT_10.get(key)
                if declared is not None and declared[0] not in str(exc):
                    wrong_reason.append(f"{name}({cell!r}): {exc}")
                continue
            escapes.extend(f"{name}({cell!r}) emitted {c}" for c in _denied_in(result))
    return refusals, unbuildable, wrong_reason, escapes


@pytest.mark.slow
def test_the_cell_that_holds_a_pole_at_resolution_10_is_not_refined() -> None:
    """Every name that refines that cell refuses it, and no other name does.

    The census of names that take a cell, handed the two pole cells and the
    two controls. What is refused for the pole cells alone is the
    limitation: ``get_children`` -- and ``cell_to_children``, the same
    function -- ``get_descendants``, ``uncompact_cells`` and
    ``cells_to_intervals``, each as it goes from resolution 10 to 11.
    """
    refusals, unbuildable, wrong_reason, escapes = _sweep_at_10()
    assert refusals == {
        key: set(cells) for key, (_, cells) in DECLARED_REFUSALS_AT_10.items()
    }
    assert wrong_reason == []
    assert escapes == []
    assert unbuildable == {
        (name, cell) for name in BLOB_READERS for cell in ROW_1000_AT_10
    }
    assert sorted(
        name
        for (name, _), (_, cells) in DECLARED_REFUSALS_AT_10.items()
        if cells == frozenset(POLE_CELLS) and name != "compactness"
    ) == [
        "cell_to_children",
        "cells_to_intervals",
        "get_children",
        "get_descendants",
        "uncompact_cells",
    ]


def test_the_sweep_at_10_reports_what_it_is_there_to_report(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The controls of that sweep, run on the ordinary cell alone.

    A name made to refuse a cell that exists is reported as refusing it,
    and a spelling the predicate denies is reported at any depth, inside a
    string or a list of them.
    """
    real = itacart.cell_shape

    def refusing(cell: str) -> Any:
        if cell == ORDINARY_AT_10:
            raise NonExistentCellError(f"{cell!r} refused for the control")
        return real(cell)

    monkeypatch.setattr(itacart, "cell_shape", refusing)
    refusals, _, _, _ = _sweep_at_10((ORDINARY_AT_10,))
    assert refusals[("cell_shape", "NonExistentCellError")] == {ORDINARY_AT_10}
    deep = DENIED_SPELLINGS[(POLE_LON, POLE_LAT, 11)]
    assert _denied_in(deep) == [deep]
    assert _denied_in([ORDINARY_AT_10, DENIED_WEST]) == [DENIED_WEST]
    assert _denied_in(ORDINARY_AT_10) == []


def test_get_children_refuses_the_pole_cell_as_it_is_iterated() -> None:
    """The refusal comes out of the iteration, not out of the call.

    ``get_children`` is lazy, as ``get_descendants`` and
    ``uncompact_cells`` are: a test that only called it would see nothing.
    """
    from itacart.exceptions import GeometryError

    cell = POLE_CELLS[0]
    assert itacart.cell_to_children is itacart.get_children
    children = itacart.get_children(cell)
    with pytest.raises(GeometryError) as caught:
        list(children)
    assert str(caught.value) == _not_refined(cell)


#: The geometries: a cap of 10 microdegrees around each pole, and a box a
#: few nanodegrees tall at a position inside the north pole's cell.
NORTH_CAP = Polygon(
    [(-180.0, 89.99999), (180.0, 89.99999), (180.0, 90.0), (-180.0, 90.0)]
)
SOUTH_CAP = Polygon(
    [(-180.0, -90.0), (180.0, -90.0), (180.0, -89.99999), (-180.0, -89.99999)]
)
POLE_LON, POLE_LAT = 155.0390625, 89.99999962902984
TINY_BOX = box(POLE_LON - 0.01, POLE_LAT - 1e-9, POLE_LON + 0.01, POLE_LAT + 1e-9)


@pytest.mark.slow
def test_a_geometry_that_reaches_the_pole_cell_is_refused_from_resolution_11() -> None:
    """The names that fill or count a geometry meet the same refusal.

    At resolution 10 they answer: 17 cells in each cap, and nothing in the
    box, which is narrower than a cell. From resolution 11 the descent has
    to enumerate the children of the pole cell, and refuses with its words.
    Pinned at 11, where it starts; each refusal costs seconds.
    ``from_geodataframe`` reprojects a frame and hands it to
    ``from_geojson``, so it answers as ``from_geojson`` does.
    """
    from itacart.exceptions import GeometryError

    assert itacart.count_internal_cells(NORTH_CAP, 10) == 17
    assert itacart.count_cells(itacart.polyfill(SOUTH_CAP, 10)) == 17
    assert itacart.count_internal_cells(TINY_BOX, 10) == 0
    with pytest.raises(GeometryError, match="narrower than one cell under 'center'"):
        itacart.from_geojson(mapping(TINY_BOX), 10)
    for call, geometry, pole in (
        (itacart.count_internal_cells, NORTH_CAP, POLE_CELLS[0]),
        (itacart.polyfill, SOUTH_CAP, POLE_CELLS[1]),
        (lambda g, r: itacart.from_geojson(mapping(g), r), TINY_BOX, POLE_CELLS[0]),
    ):
        with pytest.raises(GeometryError) as caught:
            call(geometry, 11)
        assert str(caught.value) == _not_refined(pole)


#: The four positions known to fall where the spelling names no cell, one
#: per quadrant of longitude sign and pole.
INSIDE_THE_POLE_CELLS = tuple(
    (lon, lat) for lat in (POLE_LAT, -POLE_LAT) for lon in (POLE_LON, -POLE_LON)
)

#: How each name that answers with cells, and takes no cell, is handed a
#: position. Derived from the census below rather than listed.
POSITION_CALLS: dict[str, Callable[[float, float, int], Any]] = {
    "from_geojson": lambda lon, lat, r: itacart.from_geojson(
        {"type": "Point", "coordinates": [lon, lat]}, r, containment="intersects"
    ),
    "geo_to_cell": lambda lon, lat, r: itacart.geo_to_cell(lon, lat, r),
    "latlng_to_cell": lambda lon, lat, r: itacart.latlng_to_cell(lat, lon, r),
    "sinusoidal_to_cell": lambda lon, lat, r: itacart.sinusoidal_to_cell(
        *itacart.geodetic_to_sinusoidal(lon, lat), r
    ),
    "vertex_to_cell": lambda lon, lat, r: itacart.vertex_to_cell(Point(lon, lat), r),
}

#: The rest of that census, each with the reason it is not swept here.
NOT_HANDED_A_POSITION = {
    "from_geodataframe": "reprojects a frame and hands it to from_geojson",
    "polyfill": "descends from resolution 1 and is refused at the pole cell, "
    "as the geometry test above pins",
    "refinement_alphabet": "takes a resolution alone and answers with codes",
}

#: What those names answer, by equality: at resolutions 11 and 12 the four
#: positions are spelled under the pole cells, and the predicate denies the
#: spelling. Every name answers the same, so the table is per position.
DENIED_SPELLINGS = {
    (POLE_LON, POLE_LAT, 11): "NE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(B4)))))))))))",
    (POLE_LON, POLE_LAT, 12): "NE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(B4(3))))))))))))",
    (-POLE_LON, POLE_LAT, 11): "NE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(D2)))))))))))",
    (-POLE_LON, POLE_LAT, 12): "NE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(D2(3))))))))))))",
    (POLE_LON, -POLE_LAT, 11): "SE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(B4)))))))))))",
    (POLE_LON, -POLE_LAT, 12): "SE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(B4(3))))))))))))",
    (-POLE_LON, -POLE_LAT, 11): "SE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(D2)))))))))))",
    (-POLE_LON, -POLE_LAT, 12): "SE(0000/1000(1(B2(4(E5(4(B2(4(A1(4(D2(3))))))))))))",
}


def _answers_with_cells() -> set[str]:
    """Names that take no cell and whose return annotation carries ``str``."""
    return {
        name
        for name in TAKES_NO_CELL
        if "str" in str(inspect.get_annotations(getattr(itacart, name)).get("return"))
    }


def test_the_position_census_is_whole() -> None:
    """Every name that answers with cells from a position is swept or excused."""
    assert set(POSITION_CALLS) | set(NOT_HANDED_A_POSITION) == _answers_with_cells()
    assert not set(POSITION_CALLS) & set(NOT_HANDED_A_POSITION)


def test_inside_the_pole_cells_a_position_can_be_spelled_as_no_cell() -> None:
    """At resolutions 11 and 12, four positions are answered with no cell.

    A known limitation of this release, and an exception to the rule that
    no function answers with a cell that does not exist. Pinned by
    equality, per name, position and resolution; at resolutions 10 and 13
    the same positions are answered with cells that exist.
    """
    denied: dict[tuple[str, float, float, int], list[str]] = {}
    for name, call in POSITION_CALLS.items():
        for (lon, lat), resolution in itertools.product(
            INSIDE_THE_POLE_CELLS, (10, 11, 12, 13)
        ):
            found = _denied_in(call(lon, lat, resolution))
            if found:
                denied[(name, lon, lat, resolution)] = found
    assert denied == {
        (name, lon, lat, resolution): [spelling]
        for name in POSITION_CALLS
        for (lon, lat, resolution), spelling in DENIED_SPELLINGS.items()
    }
    for spelling in DENIED_SPELLINGS.values():
        assert itacart.is_valid_index(spelling)
        assert not itacart.is_valid_cell(spelling)
        assert any(itacart.is_ancestor(pole, spelling) for pole in POLE_CELLS)


# --------------------------------------------------------------------------
# The control
# --------------------------------------------------------------------------


def test_the_instrument_reports_rather_than_passing_by_default() -> None:
    """A predicate nobody has seen fail measures nothing anybody can name.

    The sweeps above pass when the surface is right, and a sweep that
    cannot fail passes for the wrong reason just as quietly. Here the
    emission check is handed a spelling that is genuinely denied and has to
    report it.
    """
    escapes = [
        atom
        for atom in ATOM.findall(f"['{DENIED_WEST}', 'NE(0000/0110)']")
        if not itacart.is_valid_cell(atom)
    ]
    assert escapes == [DENIED_WEST]

    # And the entering sweep, shown failing against a name that answers.
    assert itacart.normalize(DENIED_WEST) == "NE(0000/0110)"
    with pytest.raises(Failed):
        with pytest.raises(NonExistentCellError):
            itacart.normalize(DENIED_WEST)
