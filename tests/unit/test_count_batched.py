"""The batched count (1.2.0) against the recursive descent it replaces (1.1.0).

The recursive descent is kept in the package as ``_count_node`` and
``_count_internal_cells_recursive``; every test here holds the batched
evaluation to it, or to counts frozen from 1.1.0 itself.

``ITACART_EQUIV_N`` scales the adversarial suite (polygons per kind and
seed, default 40). The release check runs it at 1000.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import shapely
from shapely.geometry import MultiPolygon, Polygon, box
from shapely.prepared import prep

import itacart
from itacart import geometry as G
from itacart.resolutions import cell_size, linear_refinement_ratio

L1 = G._L1
N_PER_KIND = int(os.environ.get("ITACART_EQUIV_N", "40"))

# --------------------------------------------------------------------------
# Lattice-plane cases: the descent itself, without projection
# --------------------------------------------------------------------------


def _bases(geom):
    """Base squares around ``geom`` in the lattice plane, one ring of margin."""
    minx, miny, maxx, maxy = geom.bounds
    return [
        (u, r)
        for u in range(int(minx // L1) - 1, int(maxx // L1) + 2)
        for r in range(int(miny // L1) - 1, int(maxy // L1) + 2)
    ]


def _recursive(geom, target):
    p = prep(geom)
    return sum(G._count_node(p, u * L1, r * L1, L1, 1, target) for u, r in _bases(geom))


def _batched(geom, target, strict=False):
    g = shapely.from_wkb(geom.wkb)
    shapely.prepare(g)
    b = _bases(geom)
    u0 = np.array([u * L1 for u, _ in b], dtype=float)
    v0 = np.array([r * L1 for _, r in b], dtype=float)
    return G._count_nodes_batched(g, u0, v0, L1, 1, target, strict)


def _random_node(rng, level):
    """A descent node at ``level``, reached with the descent's own arithmetic.

    Lattice magnitudes of Brazil (u ~ 7.4e6 m, v ~ 2.6e6 m), so that the
    rounding of node coordinates at resolutions 11 to 13 is the real one.
    """
    u, r = int(rng.integers(735, 742)), int(rng.integers(258, 263))
    u0, v0, side = u * L1, r * L1, L1
    for lv in range(2, level + 1):
        d = linear_refinement_ratio(lv)
        child = side / d
        col, row = int(rng.integers(d)), int(rng.integers(d))
        u0, v0, side = u0 + col * child, v0 + row * child, child
    return u0, v0, side


def _star(rng, cx, cy, radius, k):
    ang = np.sort(rng.uniform(0, 2 * np.pi, k))
    rad = rng.uniform(0.35 * radius, radius, k)
    return Polygon(np.c_[cx + rad * np.cos(ang), cy + rad * np.sin(ang)])


def _case(rng, kind):
    """``(target, geometry)``; the geometry may be empty or invalid (skipped)."""
    target = int(rng.choice([3, 5, 6, 7, 8, 9, 10, 11, 12, 13]))
    s = cell_size(target)
    if kind == "star":
        u, v, _ = _random_node(rng, min(target, 6))
        return target, _star(
            rng,
            u + rng.uniform(0, 1) * s * 40,
            v + rng.uniform(0, 1) * s * 40,
            rng.uniform(20, 150) * s,
            int(rng.integers(6, 41)),
        )
    if kind == "holes":
        u, v, _ = _random_node(rng, min(target, 6))
        radius = rng.uniform(60, 150) * s
        shell = _star(rng, u, v, radius, int(rng.integers(12, 41)))
        holes = [
            _star(
                rng,
                u + rng.uniform(-0.2, 0.2) * radius,
                v + rng.uniform(-0.2, 0.2) * radius,
                rng.uniform(0.05, 0.2) * radius,
                int(rng.integers(5, 12)),
            )
            for _ in range(int(rng.integers(1, 4)))
        ]
        holes = [h for h in holes if h.is_valid]
        if not shell.is_valid or not holes:
            return target, Polygon()
        try:
            return target, shell.difference(shapely.union_all(holes))
        except shapely.errors.GEOSException:
            return target, Polygon()
    if kind == "grid":
        # Squares built exactly as the descent builds them, so the outline
        # runs along the very squares the descent asks about.
        level = int(rng.integers(max(2, target - 3), target + 1))
        u, v, side = _random_node(rng, level - 1)
        d = linear_refinement_ratio(level)
        child = side / d
        squares = [
            box(
                u + c * child,
                v + r * child,
                u + c * child + child,
                v + r * child + child,
            )
            for r in range(d)
            for c in range(d)
            if rng.random() < 0.6
        ] or [box(u, v, u + child, v + child)]
        return target, shapely.union_all(squares)
    if kind == "centres":
        # Left and bottom edges through leaf centres, as the descent computes
        # them: the case the centre rule decides.
        u, v, side = _random_node(rng, target)
        cx, cy = u + side / 2.0, v + side / 2.0
        w, h = rng.uniform(1, 60) * side, rng.uniform(1, 60) * side
        return target, Polygon([(cx, cy), (cx + w, cy), (cx + w, cy + h), (cx, cy + h)])
    # sliver: a fraction of a cell wide, crossing many cells
    u, v, _ = _random_node(rng, min(target, 6))
    ang = rng.uniform(0, np.pi)
    length, width = rng.uniform(50, 300) * s, rng.uniform(0.05, 0.9) * s
    dx, dy = np.cos(ang), np.sin(ang)
    nx, ny = -dy * width, dx * width
    return target, Polygon(
        [
            (u, v),
            (u + dx * length, v + dy * length),
            (u + dx * length + nx, v + dy * length + ny),
            (u + nx, v + ny),
        ]
    )


KINDS = ("star", "holes", "grid", "centres", "sliver")


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("kind", KINDS)
def test_batched_descent_matches_recursive(kind, seed):
    rng = np.random.default_rng(20261115 + 101 * seed + KINDS.index(kind))
    checked = 0
    while checked < N_PER_KIND:
        target, geom = _case(rng, kind)
        if geom.is_empty or not geom.is_valid or geom.area == 0:
            continue
        assert _batched(geom, target) == _recursive(geom, target), (kind, target)
        checked += 1


@pytest.mark.parametrize("kind", KINDS)
def test_strict_leaf_matches_shortcut_on_valid_geometry(kind):
    rng = np.random.default_rng(7 + KINDS.index(kind))
    checked = 0
    while checked < 10:
        target, geom = _case(rng, kind)
        if geom.is_empty or not geom.is_valid or geom.area == 0:
            continue
        assert _batched(geom, target, strict=True) == _batched(geom, target)
        checked += 1


def test_the_worker_counts_a_chunk_as_this_process_does():
    """The worker's body, run in this process, where coverage can see it.

    In a pool it runs in a child process, which the suite's coverage does
    not measure. Handed one chunk of a real frontier, it has to give what
    the batched descent gives on the same chunk, under either leaf rule,
    and the chunks have to add up to the count of the whole descent.
    """
    rng = np.random.default_rng(5)
    geom = _star(rng, 737 * L1 + 4_000.0, 260 * L1 + 6_000.0, 3_000.0, 24)
    target = 7
    g = shapely.from_wkb(geom.wkb)
    shapely.prepare(g)
    b = _bases(geom)
    u0 = np.array([u * L1 for u, _ in b], dtype=float)
    v0 = np.array([r * L1 for _, r in b], dtype=float)
    head, fu, fv, side, level = G._expand_frontier(g, u0, v0, L1, 1, target, 64)
    assert level < target and len(fu) >= 64
    chunks = 4
    for strict in (False, True):
        parts = []
        for k in range(chunks):
            cu, cv = fu[k::chunks], fv[k::chunks]
            part = G._count_chunk_in_process(
                (geom.wkb, cu, cv, side, level, target, strict)
            )
            assert part == G._count_nodes_batched(
                g, cu, cv, side, level, target, strict
            )
            parts.append(part)
        assert head + sum(parts) == _recursive(geom, target)


def test_contains_xy_matches_contains_on_points_near_the_outline():
    rng = np.random.default_rng(3)
    shell = _star(rng, 7_385_000.0, 2_605_000.0, 4_000.0, 400)
    g = shapely.from_wkb(shell.wkb)
    shapely.prepare(g)
    ring = np.asarray(shell.exterior.coords)
    k = rng.integers(0, len(ring) - 1, 200_000)
    t = rng.random(200_000)
    x = ring[k, 0] + t * (ring[k + 1, 0] - ring[k, 0]) + rng.normal(0, 0.05, 200_000)
    y = ring[k, 1] + t * (ring[k + 1, 1] - ring[k, 1]) + rng.normal(0, 0.05, 200_000)
    assert np.array_equal(
        shapely.contains_xy(g, x, y), shapely.contains(g, shapely.points(x, y))
    )


# --------------------------------------------------------------------------
# Invalid lattice pieces: the guard
# --------------------------------------------------------------------------


def _invalid_pieces():
    """Invalid figures whose outline runs through a level-3 leaf's centre."""
    u0, v0 = 737 * L1, 260 * L1
    cx, cy = u0 + 2500.0, v0 + 2500.0
    x0, y0, x1, y1 = u0 - 100, v0 - 100, u0 + L1 + 100, v0 + L1 + 100
    outer = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    return {
        "spike_to_centre": Polygon(
            [(x0, y0), (cx, y0), (cx, cy), (cx, y0), (x1, y0), (x1, y1), (x0, y1)]
        ),
        "degenerate_hole": Polygon(
            outer, [[(cx - 10, cy), (cx + 10, cy), (cx - 10, cy)]]
        ),
        "bow_tie": Polygon(
            [
                (cx - 300, cy - 300),
                (cx + 300, cy + 300),
                (cx + 300, cy - 300),
                (cx - 300, cy + 300),
            ]
        ),
        "overlapping_parts": MultiPolygon(
            [box(x0, y0, cx, y1), box(cx - 200, y0, x1, y1)]
        ),
    }


@pytest.mark.parametrize("name", list(_invalid_pieces()))
def test_invalid_piece_counts_as_the_recursive_descent(name):
    geom = _invalid_pieces()[name]
    assert not geom.is_valid
    for target in (3, 5):
        assert _batched(geom, target, strict=True) == _recursive(geom, target)


def test_validity_of_the_lattice_piece_selects_the_leaf_rule(monkeypatch):
    seen = []
    real = G._count_leaves

    def spy(geometry, u0, v0, side, strict):
        seen.append(strict)
        return real(geometry, u0, v0, side, strict)

    monkeypatch.setattr(G, "_count_leaves", spy)
    bases = [(737 - 260, 260)]  # (column, row) with column + row = 737
    G._count_ordinary(_invalid_pieces()["spike_to_centre"], bases, 3, 1)
    assert seen and all(seen)
    seen.clear()
    G._count_ordinary(
        box(737 * L1 + 10, 260 * L1 + 10, 737 * L1 + 9000, 260 * L1 + 7000), bases, 3, 1
    )
    assert seen and not any(seen)


# --------------------------------------------------------------------------
# Public API, in EPSG:4326
# --------------------------------------------------------------------------

SQ_LON0, SQ_LAT0 = -46.6, -23.5
SQ_DLON, SQ_DLAT = 0.009790404196363954, 0.009029832224583227  # ~1 km at lat -23.5

GEO_CASES = {
    "square_1km_sp": box(SQ_LON0, SQ_LAT0, SQ_LON0 + SQ_DLON, SQ_LAT0 + SQ_DLAT),
    "equator_two_quadrants": box(-50.02, -0.01, -50.00, 0.01),
    "prime_meridian": box(-0.01, 51.47, 0.01, 51.49),
    "antimeridian_last_column": box(179.97, 10.00, 179.999, 10.03),
    "polygon_with_hole": Polygon(
        [(-47.40, -22.03), (-47.36, -22.04), (-47.35, -22.00), (-47.39, -21.99)],
        [[(-47.38, -22.02), (-47.37, -22.02), (-47.37, -22.01), (-47.38, -22.01)]],
    ),
    "multipolygon": MultiPolygon(
        [box(-48.00, -15.80, -47.99, -15.79), box(-47.98, -15.80, -47.975, -15.785)]
    ),
}

# Counted with itacart 1.1.0 (n_jobs=1), Python 3.11.18, Shapely 2.2.0, GEOS 3.14.1.
FROZEN_1_1_0 = {
    ("square_1km_sp", 5): 100,
    ("square_1km_sp", 7): 10_000,
    ("square_1km_sp", 9): 1_000_036,
    ("square_1km_sp", 11): 100_003_606,
    ("square_1km_sp", 13): 10_000_994_628,
    ("equator_two_quadrants", 5): 506,
    ("equator_two_quadrants", 7): 49_506,
    ("equator_two_quadrants", 9): 4_923_912,
    ("prime_meridian", 5): 286,
    ("prime_meridian", 7): 30_996,
    ("prime_meridian", 9): 3_090_525,
    ("antimeridian_last_column", 5): 1_050,
    ("antimeridian_last_column", 7): 105_559,
    ("antimeridian_last_column", 9): 10_549_234,
    ("polygon_with_hole", 5): 1_833,
    ("polygon_with_hole", 7): 182_905,
    ("polygon_with_hole", 9): 18_293_209,
    ("multipolygon", 5): 204,
    ("multipolygon", 7): 20_681,
    ("multipolygon", 9): 2_075_369,
}


@pytest.mark.parametrize(("name", "resolution"), sorted(FROZEN_1_1_0))
def test_count_matches_frozen_1_1_0(name, resolution):
    assert (
        itacart.count_internal_cells(GEO_CASES[name], resolution)
        == FROZEN_1_1_0[(name, resolution)]
    )


@pytest.mark.parametrize("name", sorted(GEO_CASES))
def test_count_matches_recursive_oracle(name):
    geom = GEO_CASES[name]
    for resolution in (5, 7):
        assert itacart.count_internal_cells(
            geom, resolution
        ) == G._count_internal_cells_recursive(geom, resolution)


@pytest.mark.parametrize(
    "name", ["square_1km_sp", "equator_two_quadrants", "polygon_with_hole"]
)
def test_processes_match_one_process(name, monkeypatch):
    # Force the pool even for small counts, and split the frontier early.
    monkeypatch.setattr(G, "_PARALLEL_MIN_NODES", 0)
    monkeypatch.setattr(G, "_FRONTIER_PER_JOB", 4)
    geom = GEO_CASES[name]
    assert itacart.count_internal_cells(
        geom, 9, n_jobs=2
    ) == itacart.count_internal_cells(geom, 9)


def test_the_pool_starts_its_workers_by_spawn(monkeypatch):
    """The start method is asserted on the pool, not inferred from the count.

    The count comes out the same whatever started the workers, so a test
    that only compared counts would pass under ``fork`` as well. The
    workers are spawned on every system: forked from a process with
    another thread alive -- a notebook kernel, a server -- a child can
    deadlock, and Python warns that it may.
    """
    seen = []
    real = G.ProcessPoolExecutor

    def recording(*args, **kwargs):
        seen.append(kwargs.get("mp_context"))
        return real(*args, **kwargs)

    monkeypatch.setattr(G, "ProcessPoolExecutor", recording)
    monkeypatch.setattr(G, "_PARALLEL_MIN_NODES", 0)
    monkeypatch.setattr(G, "_FRONTIER_PER_JOB", 4)
    assert (
        itacart.count_internal_cells(GEO_CASES["square_1km_sp"], 9, n_jobs=2)
        == FROZEN_1_1_0[("square_1km_sp", 9)]
    )
    assert [context.get_start_method() for context in seen] == ["spawn"]


def test_small_counts_stay_in_one_process(monkeypatch):
    def refuse(*args, **kwargs):
        raise AssertionError("a pool was started for a small count")

    monkeypatch.setattr(G, "ProcessPoolExecutor", refuse)
    assert (
        itacart.count_internal_cells(GEO_CASES["square_1km_sp"], 7, n_jobs=4) == 10_000
    )


@pytest.mark.parametrize("n_jobs", [0, -1, 1.5, True])
def test_n_jobs_is_validated(n_jobs):
    with pytest.raises(ValueError):
        itacart.count_internal_cells(GEO_CASES["square_1km_sp"], 5, n_jobs=n_jobs)
