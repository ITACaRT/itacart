"""The package's public contract, stated as an invariant rather than a list.

A name a module declares public can fail to reach the top of the package,
and a list of names cannot catch it: a list goes stale the moment a module
gains one; an invariant does not.

The rule is: whatever a module declares in its own ``__all__`` is
importable from ``itacart`` and appears in ``itacart.__all__``. Nothing
here enumerates names, so a module that grows a public function is
covered the day it grows it.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import pytest

import itacart
from itacart.exceptions import InvalidIndexError

#: Submodules that declare a public contract. Named rather than
#: discovered, so that a new module cannot join the package and escape
#: the invariant by being overlooked; ``test_every_submodule_is_listed``
#: is what keeps this list honest.
CONTRACT_MODULES = (
    "boundary",
    "cells",
    "constants",
    "engine",
    "exceptions",
    "geodesy",
    "geometry",
    "hierarchy",
    "index",
    "integer_index",
    "interop",
    "intervals",
    "metrics",
    "resolutions",
    "serialization",
    "serialization.data_tree_blob",
    "serialization.geometry_blob",
    "serialization.tree_blob",
    "topology",
)

#: Empty since ``constants`` gained its ``__all__``. Kept rather than
#: deleted, because ``test_every_submodule_is_listed`` needs somewhere to
#: account for a module that deliberately declares no contract, and an
#: empty tuple says "there is none" where a missing name would say nothing.
MODULES_WITHOUT_CONTRACT: tuple[str, ...] = ()


def _module(name: str) -> ModuleType:
    return importlib.import_module(f"itacart.{name}")


@pytest.mark.parametrize("name", CONTRACT_MODULES)
def test_every_declared_name_reaches_the_package_top(name: str) -> None:
    module = _module(name)
    declared = getattr(module, "__all__", None)
    assert declared, f"itacart.{name} declares no __all__"
    unreachable = [entry for entry in declared if not hasattr(itacart, entry)]
    unlisted = [entry for entry in declared if entry not in itacart.__all__]
    assert unreachable == []
    assert unlisted == []


@pytest.mark.parametrize("name", CONTRACT_MODULES)
def test_the_package_re_exports_the_same_object(name: str) -> None:
    """Re-export, not a second definition.

    Importing a name from the top has to give back the object the module
    defines, or the package grows a private copy that drifts.
    """
    module = _module(name)
    for entry in getattr(module, "__all__", ()):
        assert getattr(itacart, entry) is getattr(module, entry), entry


def test_the_package_contract_holds_no_duplicates() -> None:
    assert len(itacart.__all__) == len(set(itacart.__all__))


def test_every_name_in_the_package_contract_resolves() -> None:
    assert [entry for entry in itacart.__all__ if not hasattr(itacart, entry)] == []


def test_every_submodule_is_listed() -> None:
    """No submodule may sit outside both lists.

    Without this, adding a module and forgetting to name it above would
    make the invariant silently narrower rather than fail.
    """
    import pkgutil

    found = {
        info.name
        for info in pkgutil.iter_modules(itacart.__path__)
        if not info.name.startswith("_")
    }
    accounted = set(CONTRACT_MODULES) | set(MODULES_WITHOUT_CONTRACT)
    assert found - accounted == set()


@pytest.mark.parametrize("name", MODULES_WITHOUT_CONTRACT)
def test_the_uncovered_module_is_still_uncovered(name: str) -> None:
    """Any module left outside the invariant really declares no contract."""
    assert getattr(_module(name), "__all__", None) is None


def test_the_declared_contract_of_constants_is_what_the_package_re_exports() -> None:
    """The gap closed without deciding anything new about the surface.

    ``constants`` was the one module whose public names were whatever the
    package happened to re-export. Its ``__all__`` is exactly that set,
    measured rather than chosen, so the contract became written without
    becoming different. The assertion is equality: promoting a constant
    later has to come here and be seen, and quietly widening the module's
    contract without widening the package's turns this red.
    """
    from itacart import constants

    imported = {
        name for name in constants.__all__ if getattr(itacart, name, None) is not None
    }
    assert set(constants.__all__) == imported
    assert set(constants.__all__) <= set(itacart.__all__)
    assert len(itacart.__all__) == 171


# --------------------------------------------------------------------------
# H3-compatible aliases
# --------------------------------------------------------------------------


class TestH3Aliases:
    """The aliases swap argument order, and nothing else.

    They exist so that code written against H3 v4 reads unchanged, which
    means the only thing worth checking is that the swap happens and that
    the underlying answer is untouched. They were the last three lines of
    the package without a test.
    """

    def test_latlng_to_cell_swaps_the_argument_order(
        self, praca_da_se: tuple[float, float]
    ) -> None:
        lon, lat = praca_da_se
        assert itacart.latlng_to_cell(lat, lon, 13) == itacart.geo_to_cell(lon, lat, 13)

    def test_cell_to_latlng_swaps_the_returned_order(
        self, praca_da_se: tuple[float, float]
    ) -> None:
        lon, lat = praca_da_se
        cell = itacart.geo_to_cell(lon, lat, 13)
        centroid = itacart.cell_to_centroid(cell)
        assert itacart.cell_to_latlng(cell) == (centroid[1], centroid[0])

    def test_the_two_aliases_invert_each_other(
        self, praca_da_se: tuple[float, float]
    ) -> None:
        lon, lat = praca_da_se
        cell = itacart.latlng_to_cell(lat, lon, 13)
        back_lat, back_lng = itacart.cell_to_latlng(cell)
        assert itacart.latlng_to_cell(back_lat, back_lng, 13) == cell

    def test_cell_to_parent_is_get_parent(self) -> None:
        assert itacart.cell_to_parent is itacart.get_parent


def test_cell_to_latlng_refuses_a_composed_index() -> None:
    """One cell in, one point out, or a loud refusal.

    The alias exists to offer the H3 contract, and a composed index has
    no single point to give. Refusing is the only answer that does not
    require breaking the signature.

    The two-cell case is the one that mattered: the old cast to a pair
    was a lie the interpreter could not catch, so the call returned a
    pair of coordinate *pairs*, swapped, with no error at all. Three or
    more raised an unpacking error naming nothing useful.

    It lives here because the function lives in ``__init__``: a test of the
    public surface belongs with the public surface, and putting it elsewhere
    makes a partially applied change fail in a file that has nothing to do
    with the cause.
    """
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(3)]
    assert len(set(cells)) == 3, "the sample cells have to be distinct"

    single = itacart.cell_to_latlng(cells[0])
    centroid = itacart.cell_to_centroid(cells[0])
    assert single == (centroid[1], centroid[0])  # type: ignore[index]

    for count in (2, 3):
        composed = itacart.compose(cells[:count])
        assert not itacart.is_atomic(composed)
        with pytest.raises(InvalidIndexError, match=f"names {count} cells"):
            itacart.cell_to_latlng(composed)


def test_cell_to_centroid_still_answers_for_a_composed_index() -> None:
    """The refusal is the alias's, not the underlying function's."""
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(2)]
    centroids = itacart.cell_to_centroid(itacart.compose(cells))
    assert isinstance(centroids, list)
    assert len(centroids) == 2
