"""The version is written in two places, and the two must agree.

``pyproject.toml`` gives the distribution its version and
``itacart.__version__`` is what the package reports about itself, through
:func:`itacart.describe` among others. A release once moved the first and
left the second behind, and seventeen green checks let it through because
nothing compared them. These two tests compare them.
"""

from __future__ import annotations

import importlib.metadata
import pathlib
import re

import pytest

import itacart


def test_the_installed_metadata_carries_the_version_the_package_reports() -> None:
    """What ``pip`` recorded is what the package says it is.

    Skipped when nothing is installed, as when the suite runs from a
    checkout through ``PYTHONPATH`` alone: there is no metadata to read
    then. The continuous integration installs the package before testing,
    so the comparison runs there on every push. An installed copy other
    than the working tree fails here too, which is a true finding: the
    metadata would then describe a different package from the one tested.
    """
    try:
        installed = importlib.metadata.version("itacart")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("itacart is not installed, so it has no metadata to read")
    assert installed == itacart.__version__


def test_the_project_file_carries_the_version_the_package_reports() -> None:
    """The same comparison from the source tree, with nothing installed.

    The project table is read out of the file text rather than parsed: the
    only TOML parser in the standard library arrived in Python 3.11, and
    this package supports 3.10. Skipped when the package is tested from an
    installed wheel, which carries no project file.
    """
    project = pathlib.Path(itacart.__file__).parent.parent.parent / "pyproject.toml"
    if not project.is_file():
        pytest.skip("no pyproject.toml beside the package; not a source tree")
    text = project.read_text(encoding="utf-8")
    table = text[text.index("[project]") :]
    table = table[: table.index("\n[", 1)]
    declared = re.findall(r'^version = "([^"]+)"$', table, flags=re.MULTILINE)
    assert declared == [itacart.__version__]
