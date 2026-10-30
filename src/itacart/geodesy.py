"""Ellipsoidal geodesy: the parallels-plane projection and WGS84 geodesics.

Cartographic foundation of ITACaRT, from Silva, Dietzsch & Shiguemori (2025),
section 2.3, equations (1) and (2)::

    x = f1(lambda, phi) = lambda * a cos(phi) / sqrt(1 - e^2 sin^2 phi)   (1)
    y = f2(phi)         = int_0^phi a(1 - e^2) / (1 - e^2 sin^2 P)^{3/2} dP   (2)

Provenance: ``itacart_core/geodesy.py``, from the itacart-app work. Ported
rather than rewritten.

Read structurally, the two equations say something simple. The factor in
Eq. (1) is ``N(phi) cos phi``, the radius of the parallel circle, so ``x``
is true distance along the parallel. The integrand of Eq. (2) is
``M(phi)``, the meridional radius of curvature, so ``y`` is true distance
along the meridian. The projection is therefore equidistant on both axes
from the origin and equivalent by construction -- which is what carries
the equal-area property of ITACaRT.

Two consequences shape this module:

* :func:`meridian_arc` **is** the ``y`` component of
  :func:`geodetic_to_sinusoidal`. The integral is implemented once.
* The derivative of :func:`meridian_arc` is its own integrand,
  :func:`meridian_radius`, in closed form. That is what makes the Newton
  inversion in :func:`inverse_meridian_arc` cheap and unconditionally
  convergent -- ``M(phi) > 0`` everywhere, so the arc is strictly
  monotonic and has no saddle.

Units. **Every public function takes and returns degrees for angles and
metres for lengths.** There is no degree/radian boundary inside the public
surface: the radian domain lives in the ``_rad``-suffixed private helpers
and nowhere else. A unit error at this layer raises nothing --
it produces a plausible wrong answer -- so the rule is uniform rather than
convenient.

Argument order. The core is ``(lon, lat)``, cartesian ``(x, y)`` order,
throughout, geodesics included. The H3-compatible aliases take
``(lat, lng)`` and live in the package facade, not here.

**Porting warning.** ``itacart_core`` uses the opposite order:
``geodesy.forward(lat, lon)`` and ``geodesy.inverse(x, y) -> (lat, lon)``,
and ``cells.point_to_cell(lat, lon, resolution)`` on top of it. Whoever
ports ``cells.py`` has to swap both the call and the unpacking. A
swap raises nothing -- it lands on a valid coordinate somewhere else on
the globe -- so it is pinned by ``test_argument_order_is_lon_lat_not_lat_lon``
rather than left to a comment.

No PROJ anywhere: both equations are implemented directly, and the
geodesic solutions are held to a Runge-Kutta integration of the geodesic
equations and to Clairaut's relation rather than to a second library.
"""

from __future__ import annotations

import math
from typing import Final

from .constants import (
    WGS84_A,
    WGS84_B,
    WGS84_E2,
    WGS84_F,
)
from .exceptions import ConvergenceError, DomainError

__all__ = [
    "prime_vertical_radius",
    "meridian_radius",
    "meridian_arc",
    "meridian_arc_quadrature",
    "inverse_meridian_arc",
    "geodetic_to_sinusoidal",
    "sinusoidal_to_geodetic",
    "inverse_geodesic",
    "direct_geodesic",
]

# --------------------------------------------------------------------------
# Derived quantities of the meridian-arc series
# --------------------------------------------------------------------------
# The series in the third flattening n = (a - b) / (a + b) converges in
# five terms to below a nanometre for WGS84 (n ~ 1.68e-3), which is seven
# orders of magnitude finer than the 1 cm cell of resolution 13. It is the
# implementation the itacart-app validated and the one ported here; a
# composite-Simpson quadrature ships alongside it as an independent
# reference.

_N3: Final[float] = (WGS84_A - WGS84_B) / (WGS84_A + WGS84_B)
"""Third flattening of the WGS84 ellipsoid."""

_C0: Final[float] = 1.0 + _N3**2 / 4.0 + _N3**4 / 64.0
_C2: Final[float] = -1.5 * (_N3 - _N3**3 / 8.0)
_C4: Final[float] = (15.0 / 16.0) * (_N3**2 - _N3**4 / 4.0)
_C6: Final[float] = -(35.0 / 48.0) * _N3**3
_C8: Final[float] = (315.0 / 512.0) * _N3**4
_MERIDIAN_SCALE: Final[float] = WGS84_A / (1.0 + _N3)

_HALF_PI: Final[float] = math.pi / 2.0

# --------------------------------------------------------------------------
# Numerical defaults
# --------------------------------------------------------------------------

NEWTON_TOL_RAD: Final[float] = 1e-14
"""Convergence tolerance of the Newton step in :func:`inverse_meridian_arc`.

Radians. ``1e-14 rad`` is 5.7e-13 degrees, two orders below the 1e-10
degree the inversion is held to, and still above the double-precision
floor near the pole.
"""

NEWTON_MAX_ITER: Final[int] = 16
"""Iteration cap of the Newton inversion. Typical convergence is 3 steps."""

QUADRATURE_TOL_M: Final[float] = 1e-9
"""Absolute tolerance of :func:`meridian_arc_quadrature`, in metres.

A tolerance of 1e-10 m is unreachable in float64: one ULP of the
meridian quadrant (1.0e7 m) is 1.86e-9 m, so a tolerance ten times
finer than that lies below the representable
resolution of the result and the refinement loop would never exit on the
comparison. The tolerance is therefore floored at four ULP of the running
estimate.
"""

QUADRATURE_MAX_INTERVALS: Final[int] = 4096
"""Interval cap of the adaptive refinement. WGS84 converges by 8."""

VINCENTY_TOL_RAD: Final[float] = 1e-12
"""Convergence tolerance of the Vincenty iterations.

The direct solution holds its arc iterate to it in radians. The inverse
solution holds the change in its longitude iterate to this times the
line's angular length, up to one radian: a short line's azimuth inherits
that residue divided by its length, and an absolute 1e-12 radian left
0.03 degree of azimuth in a line of one centimetre. Past one radian the
threshold is this value itself, as it always was.
"""

VINCENTY_MAX_ITER: Final[int] = 200
"""Iteration cap of the Vincenty iterations.

Typical convergence is 5 to 10 steps. The cap is generous so that
exceeding it means the pair really is (near-)antipodal, which the inverse
formula does not resolve, rather than merely slow.
"""

# --------------------------------------------------------------------------
# Input validation
# --------------------------------------------------------------------------


def _check_latitude(lat_deg: float, name: str = "lat_deg") -> float:
    """Return ``lat_deg`` as a float, rejecting non-finite or out-of-range."""
    value = float(lat_deg)
    if not math.isfinite(value):
        raise DomainError(f"{name} must be finite, got {lat_deg!r}")
    if not -90.0 <= value <= 90.0:
        raise DomainError(f"{name} must lie in [-90, 90] degrees, got {value}")
    return value


def _check_finite(value: float, name: str) -> float:
    """Return ``value`` as a float, rejecting NaN and infinity."""
    result = float(value)
    if not math.isfinite(result):
        raise DomainError(f"{name} must be finite, got {value!r}")
    return result


# --------------------------------------------------------------------------
# Radii of curvature
# --------------------------------------------------------------------------


def _prime_vertical_radius_rad(phi: float) -> float:
    """Prime-vertical radius ``N(phi)``, latitude in radians."""
    return WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(phi) ** 2)


def _meridian_radius_rad(phi: float) -> float:
    """Meridional radius ``M(phi)``, latitude in radians."""
    # math.pow rather than ** so the return type is float and not Any
    # under mypy --strict; the C call underneath is identical.
    return (
        WGS84_A * (1.0 - WGS84_E2) / math.pow(1.0 - WGS84_E2 * math.sin(phi) ** 2, 1.5)
    )


def prime_vertical_radius(lat_deg: float) -> float:
    """Return the prime-vertical radius of curvature ``N`` at a latitude.

    ``N(phi) = a / sqrt(1 - e^2 sin^2 phi)``. Together with ``cos phi`` it
    forms the radius of the parallel circle, which is the scale factor of
    Eq. (1).

    Args:
        lat_deg: Geodetic latitude, in degrees.

    Returns:
        Radius in metres, between ``a`` at the equator and ``a / sqrt(1 -
        e^2)`` at the pole.

    Raises:
        DomainError: If ``lat_deg`` is not finite or lies outside
            ``[-90, 90]``.
    """
    return _prime_vertical_radius_rad(math.radians(_check_latitude(lat_deg)))


def meridian_radius(lat_deg: float) -> float:
    """Return the meridional radius of curvature ``M`` at a latitude.

    ``M(phi) = a(1 - e^2) / (1 - e^2 sin^2 phi)^{3/2}``. This is both the
    integrand of Eq. (2) and, by the fundamental theorem of calculus, the
    derivative of :func:`meridian_arc` -- which is the derivative the
    Newton iteration of :func:`inverse_meridian_arc` uses in closed form.

    Args:
        lat_deg: Geodetic latitude, in degrees.

    Returns:
        Radius in metres. Strictly positive everywhere, hence the meridian
        arc is strictly monotonic.

    Raises:
        DomainError: If ``lat_deg`` is not finite or lies outside
            ``[-90, 90]``.
    """
    return _meridian_radius_rad(math.radians(_check_latitude(lat_deg)))


# --------------------------------------------------------------------------
# Meridian arc -- Eq. (2)
# --------------------------------------------------------------------------


def _meridian_arc_rad(phi: float) -> float:
    """Meridian arc from the equator, latitude in radians."""
    return _MERIDIAN_SCALE * (
        _C0 * phi
        + _C2 * math.sin(2.0 * phi)
        + _C4 * math.sin(4.0 * phi)
        + _C6 * math.sin(6.0 * phi)
        + _C8 * math.sin(8.0 * phi)
    )


def meridian_arc(lat_deg: float) -> float:
    """Return the meridian arc length from the equator to a latitude.

    Evaluates Eq. (2) of the paper through its series expansion in the
    third flattening: five terms, closed form, no iteration. The function
    is odd, so southern latitudes return negative lengths, and it is the
    ``y`` component of :func:`geodetic_to_sinusoidal`.

    :func:`meridian_arc_quadrature` evaluates the same integral by
    composite Simpson quadrature and agrees to a few nanometres. It exists
    as an independent check on this one; this is the implementation on the
    hot path.

    Args:
        lat_deg: Geodetic latitude, in degrees.

    Returns:
        Signed arc length in metres, ``0`` at the equator and
        ``+-10 001 965.7293 m`` at the poles.

    Raises:
        DomainError: If ``lat_deg`` is not finite or lies outside
            ``[-90, 90]``.

    Example:
        >>> round(meridian_arc(90.0), 4)
        10001965.7293
        >>> meridian_arc(0.0)
        0.0
    """
    return _meridian_arc_rad(math.radians(_check_latitude(lat_deg)))


def meridian_arc_quadrature(
    lat_deg: float,
    *,
    tol: float = QUADRATURE_TOL_M,
    max_intervals: int = QUADRATURE_MAX_INTERVALS,
) -> float:
    """Return the meridian arc by composite Simpson quadrature of Eq. (2).

    Integrates :func:`meridian_radius` directly, doubling the interval
    count until two successive estimates agree to ``tol``. Reference
    implementation: it does not share a line of arithmetic with the series
    of :func:`meridian_arc`, so agreement between the two is evidence
    about the integral itself rather than about one coefficient table.

    The effective tolerance is ``max(tol, 4 * ulp(estimate))``. At the
    quadrant one ULP is 1.86e-9 m, so a tolerance finer than that cannot
    be met in float64 and asking for it would only spend intervals
    at no gain in accuracy.

    Args:
        lat_deg: Geodetic latitude, in degrees.
        tol: Absolute tolerance in metres, floored as described above.
        max_intervals: Cap on the number of Simpson intervals.

    Returns:
        Signed arc length in metres.

    Raises:
        DomainError: If ``lat_deg`` is not finite or lies outside
            ``[-90, 90]``, or if ``tol`` is not positive.
        ConvergenceError: If the refinement reaches ``max_intervals``
            without meeting the tolerance.
    """
    phi = math.radians(_check_latitude(lat_deg))
    if _check_finite(tol, "tol") <= 0.0:
        raise DomainError(f"tol must be positive, got {tol}")
    if phi == 0.0:
        return 0.0

    previous = math.inf
    intervals = 2
    while intervals <= max_intervals:
        step = phi / intervals
        terms = [_meridian_radius_rad(0.0), _meridian_radius_rad(phi)]
        terms.extend(
            (4.0 if i % 2 else 2.0) * _meridian_radius_rad(i * step)
            for i in range(1, intervals)
        )
        # Compensated summation: a naive loop accumulates round-off that
        # grows with the interval count, and the refinement then chases a
        # difference dominated by that noise instead of by the quadrature.
        estimate = math.fsum(terms) * step / 3.0
        floor_tol = 4.0 * math.ulp(abs(estimate))
        if abs(estimate - previous) < max(tol, floor_tol):
            return estimate
        previous = estimate
        intervals *= 2

    raise ConvergenceError(
        f"meridian arc quadrature did not reach tol={tol} m within "
        f"{max_intervals} intervals at lat={lat_deg} deg"
    )


def _inverse_meridian_arc_rad(
    arc_m: float, tol: float, max_iter: int
) -> tuple[float, int]:
    """Invert the meridian arc, returning ``(phi_rad, iterations)``."""
    phi = arc_m / (_MERIDIAN_SCALE * _C0)
    for iteration in range(1, max_iter + 1):
        delta = (arc_m - _meridian_arc_rad(phi)) / _meridian_radius_rad(phi)
        phi += delta
        if abs(delta) < tol:
            return phi, iteration
    raise ConvergenceError(
        f"inverse meridian arc did not converge in {max_iter} Newton "
        f"iterations for arc={arc_m} m"
    )


def inverse_meridian_arc(
    arc_m: float,
    *,
    tol: float = NEWTON_TOL_RAD,
    max_iter: int = NEWTON_MAX_ITER,
) -> float:
    """Return the latitude whose meridian arc from the equator is ``arc_m``.

    Newton iteration on ``meridian_arc(phi) - arc_m``, with
    :func:`meridian_radius` as the derivative in closed form. Since
    ``M(phi) > 0`` everywhere the arc is strictly monotonic, so the
    iteration has no saddle and the rectifying-latitude first guess
    converges in about three steps.

    Args:
        arc_m: Signed meridian arc length from the equator, in metres.
        tol: Convergence tolerance on the latitude increment, in radians.
        max_iter: Iteration cap.

    Returns:
        Geodetic latitude in degrees, exactly inverting
        :func:`meridian_arc`.

    Raises:
        DomainError: If ``arc_m`` is not finite or exceeds the meridian
            quadrant in absolute value -- there is no latitude beyond the
            pole.
        ConvergenceError: If Newton fails to converge within ``max_iter``.

    Example:
        >>> round(inverse_meridian_arc(meridian_arc(-33.75)), 12)
        -33.75
    """
    arc = _check_finite(arc_m, "arc_m")
    quadrant = _meridian_arc_rad(_HALF_PI)
    if abs(arc) > quadrant:
        raise DomainError(
            f"arc_m={arc} m lies beyond the meridian quadrant "
            f"({quadrant:.3f} m); no latitude corresponds to it"
        )
    phi, _ = _inverse_meridian_arc_rad(arc, tol, max_iter)
    return math.degrees(phi)


# --------------------------------------------------------------------------
# Parallels-plane projection -- Eq. (1) and (2)
# --------------------------------------------------------------------------


def geodetic_to_sinusoidal(lon_deg: float, lat_deg: float) -> tuple[float, float]:
    """Project a geodetic position onto the parallels plane.

    Applies Eq. (1) to the abscissa and Eq. (2) to the ordinate. The
    result is signed: eastern longitudes give positive ``x``, northern
    latitudes positive ``y``. Quadrant mirroring, which is what keeps the
    ITACaRT index free of negative components, happens in
    :mod:`itacart.cells`, not here.

    Args:
        lon_deg: Longitude in degrees, relative to the prime meridian.
        lat_deg: Geodetic latitude in degrees.

    Returns:
        ``(x, y)`` in metres.

    Raises:
        DomainError: If either argument is not finite, or ``lat_deg`` is
            outside ``[-90, 90]``.

    Note:
        Longitude is not range-checked. Eq. (1) is linear in lambda and
        the addressable ITACaRT domain -- including the antemeridian
        extension zones of section 3.2 -- is decided in ``boundary``.
        Enforcing ``[-180, 180]`` here would pre-empt that decision and
        reject the extension zones.

    Example:
        >>> x, y = geodetic_to_sinusoidal(0.0, 0.0)
        >>> (x, y)
        (0.0, 0.0)
    """
    lat = _check_latitude(lat_deg)
    lam = math.radians(_check_finite(lon_deg, "lon_deg"))
    phi = math.radians(lat)
    x = lam * _prime_vertical_radius_rad(phi) * math.cos(phi)
    y = _meridian_arc_rad(phi)
    return x, y


def sinusoidal_to_geodetic(
    x_m: float,
    y_m: float,
    *,
    tol: float = NEWTON_TOL_RAD,
    max_iter: int = NEWTON_MAX_ITER,
) -> tuple[float, float]:
    """Invert the parallels-plane projection back to geodetic coordinates.

    Recovers the latitude from ``y`` through :func:`inverse_meridian_arc`,
    then divides ``x`` by the radius of the parallel circle at that
    latitude to recover the longitude.

    Args:
        x_m: Abscissa on the projection plane, in metres.
        y_m: Ordinate on the projection plane, in metres.
        tol: Newton tolerance, passed to the arc inversion, in radians.
        max_iter: Newton iteration cap.

    Returns:
        ``(lon_deg, lat_deg)`` in degrees.

    Raises:
        DomainError: If either argument is not finite, or ``y_m`` lies
            beyond the meridian quadrant.
        ConvergenceError: Propagated from the Newton iteration.

    Note:
        At a pole the parallel circle collapses to a point and every
        longitude projects to ``x = 0``, so longitude is not recoverable.
        The function returns ``0.0`` there by convention rather than
        raising: the pole is a legitimate position, only its longitude is
        undefined.
    """
    _check_finite(x_m, "x_m")
    lat_deg = inverse_meridian_arc(y_m, tol=tol, max_iter=max_iter)
    phi = math.radians(lat_deg)
    cos_phi = math.cos(phi)
    if abs(cos_phi) < 1e-15:
        return 0.0, lat_deg
    lam = x_m / (_prime_vertical_radius_rad(phi) * cos_phi)
    return math.degrees(lam), lat_deg


# --------------------------------------------------------------------------
# Geodesics -- Vincenty (1975)
# --------------------------------------------------------------------------
#
# Vincenty's nested equations solve the inverse and direct geodesic
# problems on the ellipsoid. Ported from itacart_core/geodesy.py. The
# series is checked against a numerical solution of the differential
# equations a geodesic satisfies, which agrees to under a nanodegree over
# an enumerated grid of 568 cases, and against Clairaut's relation, which
# closes at machine precision. The
# package needs them for orthodromic densification: a straight line
# on the sinusoidal plane is not a geodesic on the ellipsoid, so a long
# edge has to be resampled along the true geodesic before it is filled.
#
# Reference: T. Vincenty, "Direct and inverse solutions of geodesics on
# the ellipsoid with application of nested equations", Survey Review
# XXIII (176), April 1975, pp. 88-93.


def _vincenty_inverse(
    lon1: float,
    lat1: float,
    lon2: float,
    lat2: float,
    tol: float,
    max_iter: int,
) -> tuple[float, float, int]:
    """Vincenty inverse, instrumented.

    Returns ``(distance_m, azimuth_deg, iterations)``. The iteration count
    locates the near-antipodal frontier, where convergence slows; the public
    wrapper drops it.

    The iteration stops when the change in the longitude iterate falls
    below ``tol`` times ``min(1, sigma)``, the line's angular length up to
    one radian, and not below ``tol`` itself. The azimuth is read off the
    converged iterate, and a residue ``d`` in it moves the azimuth by
    about ``d / sin(sigma)``: an absolute threshold of 1e-12 radian left
    3.4e-2 degree in a line of one centimetre and 2.8e-4 in a line of one
    metre, against an independent solution, while the distance stayed
    within micrometres. Scaled by the length, a short line's azimuth lands
    on the rounding floor of its own inputs.

    Past one radian, about 6 400 km, the threshold is ``tol`` itself, as it
    always was. The sine shrinks again towards the antipode, and a
    threshold scaled by it there asked the slow near-antipodal iteration
    for digits it cannot reach within its cap: it refused lines it had
    answered, to a tenth of a millimetre. With the default tolerance a line
    longer than a radian is solved exactly as before, near-antipodal
    frontier included, and so is the pair the port is pinned against,
    which keeps its five iterations and its last bit.

    The threshold never goes below four units in the last place of the
    iterate, which is the most a float can resolve there, so a line the
    iterate cannot tell from its neighbour -- a centimetre that crosses
    many meridians near a pole -- still stops rather than circling.
    """

    if lat1 == lat2 and lon1 == lon2:
        return 0.0, 0.0, 0

    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_lon = math.radians(lon2 - lon1)

    one_minus_f = 1.0 - WGS84_F
    tan_u1 = one_minus_f * math.tan(phi1)
    tan_u2 = one_minus_f * math.tan(phi2)
    cos_u1 = 1.0 / math.sqrt(1.0 + tan_u1 * tan_u1)
    sin_u1 = tan_u1 * cos_u1
    cos_u2 = 1.0 / math.sqrt(1.0 + tan_u2 * tan_u2)
    sin_u2 = tan_u2 * cos_u2

    lam = delta_lon
    iterations = 0
    cos_sq_alpha = 0.0
    sin_sigma = 0.0
    cos_sigma = 0.0
    sigma = 0.0
    cos_2sigma_m = 0.0
    sin_lam = 0.0
    cos_lam = 0.0

    for iterations in range(1, max_iter + 1):
        sin_lam = math.sin(lam)
        cos_lam = math.cos(lam)
        sin_sigma = math.sqrt(
            (cos_u2 * sin_lam) ** 2 + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam) ** 2
        )
        if sin_sigma == 0.0:  # pragma: no cover - defensive, see note below
            # Coincident after reduction to the auxiliary sphere. Kept from
            # the origin implementation; the equality guard above already
            # absorbs every input that can reach it in float64, so the
            # branch is defensive rather than live.
            return 0.0, 0.0, iterations
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
        cos_sq_alpha = 1.0 - sin_alpha * sin_alpha
        if cos_sq_alpha == 0.0:
            cos_2sigma_m = 0.0
        else:
            cos_2sigma_m = cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_sq_alpha
        c = WGS84_F / 16.0 * cos_sq_alpha * (4.0 + WGS84_F * (4.0 - 3.0 * cos_sq_alpha))
        lam_prev = lam
        lam = delta_lon + (1.0 - c) * WGS84_F * sin_alpha * (
            sigma
            + c
            * sin_sigma
            * (
                cos_2sigma_m
                + c * cos_sigma * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m)
            )
        )
        change = abs(lam - lam_prev)
        if change < tol * min(1.0, sigma) or change <= 4.0 * math.ulp(lam):
            break
    else:
        raise ConvergenceError(
            f"Vincenty inverse did not converge in {max_iter} iterations "
            f"(near-antipodal pair): ({lon1}, {lat1}) -> ({lon2}, {lat2})"
        )

    u_sq = cos_sq_alpha * (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_B * WGS84_B)
    a_v = 1.0 + u_sq / 16384.0 * (
        4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq))
    )
    b_v = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)))
    delta_sigma = (
        b_v
        * sin_sigma
        * (
            cos_2sigma_m
            + b_v
            / 4.0
            * (
                cos_sigma * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m)
                - b_v
                / 6.0
                * cos_2sigma_m
                * (-3.0 + 4.0 * sin_sigma * sin_sigma)
                * (-3.0 + 4.0 * cos_2sigma_m * cos_2sigma_m)
            )
        )
    )
    distance = WGS84_B * a_v * (sigma - delta_sigma)

    azimuth = math.atan2(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam)
    return distance, math.degrees(azimuth) % 360.0, iterations


def inverse_geodesic(
    lon1_deg: float,
    lat1_deg: float,
    lon2_deg: float,
    lat2_deg: float,
    *,
    tol: float = VINCENTY_TOL_RAD,
    max_iter: int = VINCENTY_MAX_ITER,
) -> tuple[float, float]:
    """Return the geodesic distance and initial azimuth between two points.

    Vincenty inverse solution on the WGS84 ellipsoid.

    Args:
        lon1_deg: Longitude of the first point, in degrees.
        lat1_deg: Latitude of the first point, in degrees.
        lon2_deg: Longitude of the second point, in degrees.
        lat2_deg: Latitude of the second point, in degrees.
        tol: Convergence tolerance on the longitude iterate, relative to
            the line's angular length up to one radian: the iteration
            stops when the iterate moves by less than
            ``tol * min(1, sigma)`` radians, or by no more than four units
            in its last place. Relative because a short line's azimuth
            inherits the residue divided by its length, so an absolute
            threshold that suits a long line leaves a short one with a
            wrong azimuth; capped at one radian, so that a long line,
            near-antipodal ones included, converges as it always did.
        max_iter: Iteration cap.

    Returns:
        ``(distance_m, azimuth_deg)``, the geodesic length in metres and
        the forward azimuth at the first point in degrees, clockwise from
        north, normalised to ``[0, 360)``.

    Raises:
        DomainError: If any argument is not finite, or a latitude lies
            outside ``[-90, 90]``.
        ConvergenceError: If the iteration does not converge, which for
            this formula means a (near-)antipodal pair. Returning an
            unconverged iterate would look exactly like a distance and
            silently poison whatever consumes it.

    Note:
        Coincident points return ``(0.0, 0.0)``. The azimuth is a
        convention there, not a bearing -- the geodesic is degenerate and
        has no direction. Raising instead would force every caller to
        guard the duplicated-vertex case that vertex deduplication in
        :mod:`itacart.geometry` exists to absorb.

    Example:
        >>> d, az = inverse_geodesic(0.0, 0.0, 0.0, 1.0)
        >>> round(d, 3), round(az, 6)
        (110574.389, 0.0)
    """
    distance, azimuth, _ = _vincenty_inverse(
        _check_finite(lon1_deg, "lon1_deg"),
        _check_latitude(lat1_deg, "lat1_deg"),
        _check_finite(lon2_deg, "lon2_deg"),
        _check_latitude(lat2_deg, "lat2_deg"),
        tol,
        max_iter,
    )
    return distance, azimuth


def direct_geodesic(
    lon1_deg: float,
    lat1_deg: float,
    azimuth_deg: float,
    distance_m: float,
    *,
    tol: float = VINCENTY_TOL_RAD,
    max_iter: int = VINCENTY_MAX_ITER,
) -> tuple[float, float]:
    """Return the end point of a geodesic given start, azimuth and length.

    Vincenty direct solution on the WGS84 ellipsoid. The end point lies on
    the unique geodesic leaving the start point along ``azimuth_deg``, at
    arc length ``distance_m``.

    Args:
        lon1_deg: Longitude of the start point, in degrees.
        lat1_deg: Latitude of the start point, in degrees.
        azimuth_deg: Forward azimuth at the start point, in degrees
            clockwise from north.
        distance_m: Geodesic distance, in metres. Zero returns the start
            point unchanged.

    Returns:
        ``(lon2_deg, lat2_deg)`` in degrees, longitude normalised to
        ``(-180, 180]``.

    Raises:
        DomainError: If any argument is not finite, ``lat1_deg`` lies
            outside ``[-90, 90]``, or ``distance_m`` is negative. A
            negative distance is a sign error at the call site, not a
            reverse traverse; accepting it would return a plausible point
            on the opposite side.
        ConvergenceError: If the iteration does not converge.
    """
    lat1 = _check_latitude(lat1_deg, "lat1_deg")
    lon1 = _check_finite(lon1_deg, "lon1_deg")
    azimuth = _check_finite(azimuth_deg, "azimuth_deg")
    distance = _check_finite(distance_m, "distance_m")
    if distance < 0.0:
        raise DomainError(
            f"distance_m must be non-negative, got {distance}; to travel "
            "the other way, add 180 to the azimuth"
        )
    if distance == 0.0:
        return lon1, lat1

    phi1 = math.radians(lat1)
    alpha1 = math.radians(azimuth)

    one_minus_f = 1.0 - WGS84_F
    tan_u1 = one_minus_f * math.tan(phi1)
    cos_u1 = 1.0 / math.sqrt(1.0 + tan_u1 * tan_u1)
    sin_u1 = tan_u1 * cos_u1

    sigma_1 = math.atan2(tan_u1, math.cos(alpha1))
    sin_alpha = cos_u1 * math.sin(alpha1)
    cos_sq_alpha = 1.0 - sin_alpha * sin_alpha
    u_sq = cos_sq_alpha * (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_B * WGS84_B)
    a_v = 1.0 + u_sq / 16384.0 * (
        4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq))
    )
    b_v = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)))

    sigma = distance / (WGS84_B * a_v)
    sin_sigma = 0.0
    cos_sigma = 0.0
    cos_2sigma_m = 0.0
    for _ in range(max_iter):
        cos_2sigma_m = math.cos(2.0 * sigma_1 + sigma)
        sin_sigma = math.sin(sigma)
        cos_sigma = math.cos(sigma)
        delta_sigma = (
            b_v
            * sin_sigma
            * (
                cos_2sigma_m
                + b_v
                / 4.0
                * (
                    cos_sigma * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m)
                    - b_v
                    / 6.0
                    * cos_2sigma_m
                    * (-3.0 + 4.0 * sin_sigma * sin_sigma)
                    * (-3.0 + 4.0 * cos_2sigma_m * cos_2sigma_m)
                )
            )
        )
        sigma_prev = sigma
        sigma = distance / (WGS84_B * a_v) + delta_sigma
        if abs(sigma - sigma_prev) < tol:
            break
    else:
        raise ConvergenceError(
            f"Vincenty direct did not converge in {max_iter} iterations "
            f"from ({lon1}, {lat1}) on azimuth {azimuth} over {distance} m"
        )

    cos_alpha1 = math.cos(alpha1)
    phi2 = math.atan2(
        sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alpha1,
        one_minus_f
        * math.sqrt(
            sin_alpha * sin_alpha
            + (sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alpha1) ** 2
        ),
    )
    lam = math.atan2(
        sin_sigma * math.sin(alpha1),
        cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alpha1,
    )
    c = WGS84_F / 16.0 * cos_sq_alpha * (4.0 + WGS84_F * (4.0 - 3.0 * cos_sq_alpha))
    lon_diff = lam - (1.0 - c) * WGS84_F * sin_alpha * (
        sigma
        + c
        * sin_sigma
        * (cos_2sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m))
    )
    lon2 = ((lon1 + math.degrees(lon_diff) + 540.0) % 360.0) - 180.0
    return lon2, math.degrees(phi2)
