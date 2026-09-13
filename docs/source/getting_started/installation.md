# Installation

```console
$ pip install itacart
```

Python 3.10 or newer. The only hard runtime dependency is `shapely`; the
ellipsoidal sinusoidal projection is computed from the paper's equations, so
PROJ is never called at run time and there is no compiled geospatial stack to
install.

## Optional extras

Two extras add capability rather than replacing anything:

```console
$ pip install "itacart[geo]"       # GeoDataFrame import and export
$ pip install "itacart[parallel]"  # cell filling across several processes
```

`geo` pulls in `geopandas` and enables the GeoDataFrame side of
[GeoJSON and interoperability](../concepts/data_and_interoperability/index.md).
`parallel` pulls in `joblib` and is what lets
[polyfill](../concepts/spatial_operations/cell_filling_polyfill.md) use more
than one worker. Both can be combined:

```console
$ pip install "itacart[geo,parallel]"
```

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
