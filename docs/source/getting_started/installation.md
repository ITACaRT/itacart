# Installation

```console
$ pip install itacart
```

Python 3.10 or newer. The runtime dependencies are `shapely` and `numpy`, and
Shapely requires `numpy` as well, so nothing is installed beyond what Shapely
brings. The ellipsoidal sinusoidal projection is computed from the paper's
equations, so PROJ is never called at run time and there is no compiled
geospatial stack to install.

## Optional extras

One extra adds capability:

```console
$ pip install "itacart[geo]"       # GeoDataFrame import and export
```

`geo` pulls in `geopandas` and enables the GeoDataFrame side of
[GeoJSON and interoperability](../concepts/data_and_interoperability/index.md).

Workers need no extra. `count_internal_cells` spreads a large count over
processes from the standard library, and
[polyfill](../concepts/spatial_operations/cell_filling_polyfill.md) runs
serially whatever `n_jobs` says. Release 1.2.0 removed the `parallel` extra,
which installed `joblib` and nothing in the package used.

## From source

```console
$ git clone https://github.com/itacart/itacart
$ cd itacart
$ pip install .
```

Setting up a development environment, running the test suite and building this
documentation are covered in
[CONTRIBUTING.md](https://github.com/itacart/itacart/blob/main/CONTRIBUTING.md).

## Checking the installation

```console
$ python -c "import itacart; print(itacart.__version__)"
```

If that prints a version and
[Quick Start](quick_start.md) runs to the end, the installation is complete.
