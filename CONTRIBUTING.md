# Contributing

## Environment

```bash
git clone https://github.com/ICartCWB/itacart
cd itacart
pip install -e ".[dev,docs,geo]"
```

## Before opening a PR

```bash
black src tests
isort src tests
flake8 src tests --max-line-length=88 --extend-ignore=E203,W503
mypy src
pytest
```

## Documentation

The build imports the package from `src/`, not from `site-packages`. **No
installation is required**, and an installed copy that shadows the working tree
makes `conf.py` fail loudly rather than document the wrong code.

The `dev,docs,geo` install above already covers it. On its own, the `docs`
extra provides Sphinx, the Read the Docs theme, MyST-NB and the Sphinx
plugins:

```bash
pip install -e ".[docs]"
```

The worked examples under Getting Started execute during the build, so the
environment has to be able to import the package it documents. An example
whose code raises stops the build.

### Build it

```bash
cd docs
make html          # make.bat html on Windows
```

The result is `docs/_build/html/index.html`.

### Build it strictly

The CI job treats warnings as errors, and so should a local check before
pushing:

```bash
cd docs
make strict        # make.bat strict on Windows
```

A broken cross-reference, a missing image or a malformed docstring becomes a
build failure rather than a silently ugly page. Run this, not `make html`,
before opening a PR that touches `docs/` or a docstring.

### Live preview

`sphinx-autobuild` rebuilds on save and reloads the browser, which is the
useful mode while writing docstrings. It is not in the `docs` extra:

```bash
pip install sphinx-autobuild
sphinx-autobuild docs/source docs/_build/html --watch src
```

The `--watch src` is the point: a docstring edited in `src/itacart/` rebuilds
the page that documents it. Note that each rebuild re-executes the code of
any executable page that changed.

### Checking for undocumented code

```bash
cd docs
make coverage
cat _build/coverage/python.txt
```

`sphinx.ext.coverage` lists every public name reached by autodoc that carries
no docstring. On a package whose docstrings are the specification, that list
should stay empty.

### What is generated, and what is written

Almost nothing under `docs/source/` is maintained by hand.

| Path | Origin |
|---|---|
| `api/generated/` | `sphinx.ext.autosummary`, walking the package recursively |
| `_generated/figures/` | the `phase_figures` extension, scanning `docs/_static/fN/` |
| `_build/` | Sphinx output |
| everything else | written by hand: `conf.py`, `index.rst`, and every page under `getting_started/` and `concepts/` |

Both generated directories are rewritten on every build and are ignored by
Git. A stale one survives a rebuild, so delete them when a build behaves
oddly.

Adding a module requires **no edit here**: `autosummary` finds it through
`itacart`'s own namespace. Adding a set of figures requires no edit either —
drop the PNGs in `docs/_static/fN/`, named `fN_MM_slug.png`, and the gallery
page appears with captions derived from the filenames.

## Code discipline

- Docstrings, identifiers and comments in English.
- Type hints required on every public signature.
- Coverage target >= 85% per module.
- A fixed bug earns a permanent regression test.
- An architectural decision is recorded with its rationale **and the
  alternatives that were rejected**. Recording the rejection is what keeps a
  later change from relitigating a question that was already settled.

## Suggested order of implementation

The dependency graph between modules sets the order:

1. `constants`, `exceptions` — no dependencies
2. `geodesy` — projection and Vincenty
3. `index` — parser, compose, normalize
4. `resolutions` — table and scales
5. `cells` — quantization and inverse geometry
6. `boundary` — border families (depends on cells)
7. `hierarchy` — navigation
8. `topology` — neighbourhood (depends on boundary for deflection)
9. `geometry` — polyfill, densification, canonicalization
10. `serialization` — TreeBlob and GeometryBlob
11. `interop`, `engine`

## How the work is sequenced

Development proceeds in numbered phases, each one closed by a handoff that
records what was delivered, which decisions were taken, and what was left
open. The phase documents are coordination artefacts and are not versioned
here; the CHANGELOG is the record this repository keeps.
