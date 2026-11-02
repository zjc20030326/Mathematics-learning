"""Independent checks of the Jost kernel and the two spectral infinities.

Run from the thesis directory with ``python calculations/verify_jost_asymptotics.py``.
The ODE experiment uses a smooth localized test potential, not the explicit
one-orbit solution, and integrates the spatial Lax equation directly.
"""

from __future__ import annotations

import json

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm


def symbolic_checks() -> dict[str, bool]:
    z, q0 = sp.symbols("z q0", nonzero=True)
    q, qc, qx, qcx, j = sp.symbols("q qc qx qcx j")
    sigma = sp.diag(1, -1)
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    eye = sp.eye(2)
    k = (z - q0**2 / z) / 2
    lam = (z + q0**2 / z) / 2
    potential = sp.Matrix([[0, q], [-qc, 0]])
    background = sp.Matrix([[0, q0], [-q0, 0]])
    b = sigma * potential
    bx = sp.Matrix([[0, qx], [qcx, 0]])
    u = -sp.I * k * sigma + potential
    ub = -sp.I * k * sigma + background
    transform = eye - sp.I * b / z
    residual = u * transform + sp.I * lam * transform * sigma + sp.I * bx / z
    expected = sp.I * ((q0**2 - q * qc) * sigma + bx) / z
    a = -sp.I * b + sp.I * j * sigma
    small_constant = sp.Matrix([[q, -j], [j, qc]]) / q0

    def zero(matrix: sp.Matrix) -> bool:
        return all(sp.simplify(entry) == 0 for entry in matrix)

    checks = {
        "background_characteristic_identity": zero(ub**2 + lam**2 * eye),
        "transformed_spatial_equation": zero(residual - expected),
        "large_z_off_diagonal_coefficient": zero(-sp.I * (sigma * a - a * sigma) / 2 + potential),
        "small_z_constant_from_sheet_change": zero(sp.I * a * sigma1 / q0 - small_constant),
    }
    assert all(checks.values()), checks
    return checks


SIGMA = np.diag([1.0, -1.0]).astype(complex)
SIGMA1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
EYE = np.eye(2, dtype=complex)
Q0 = 1.0
ENDPOINT = 8.0


def potential(x: float) -> complex:
    return Q0 + 0.22 * np.exp(-(x - 0.27) ** 2) * (1.0 + 1j * (0.11 + 0.35 * x))


def spectral(z: complex) -> tuple[complex, complex]:
    return (z - Q0**2 / z) / 2, (z + Q0**2 / z) / 2


def qmatrix(q: complex) -> np.ndarray:
    return np.array([[0.0, q], [-q.conjugate(), 0.0]], dtype=complex)


def jost_column(z: complex, endpoint_sign: int, column: int) -> np.ndarray:
    """Integrate the original modified spatial equation, with no asymptotic ansatz."""
    k, lam = spectral(z)
    s = 1 if column == 0 else -1
    y = EYE - 1j * Q0 * SIGMA1 / z

    def equation(x: float, value: np.ndarray) -> np.ndarray:
        return (-1j * k * SIGMA + qmatrix(potential(x)) + 1j * s * lam * EYE) @ value

    result = solve_ivp(
        equation,
        (endpoint_sign * ENDPOINT, 0.0),
        y[:, column],
        method="DOP853",
        rtol=2e-12,
        atol=2e-13,
    )
    assert result.success, result.message
    return result.y[:, -1]


def phase_integral(endpoint_sign: int) -> float:
    value, _ = quad(lambda x: Q0**2 - abs(potential(x)) ** 2, -ENDPOINT, 0.0, epsabs=1e-13)
    if endpoint_sign < 0:
        return value
    value, _ = quad(lambda x: Q0**2 - abs(potential(x)) ** 2, 0.0, ENDPOINT, epsabs=1e-13)
    return -value


def asymptotic_checks() -> dict[str, object]:
    rows = []
    q = potential(0.0)
    b = SIGMA @ qmatrix(q)
    # D_+ has mu_{-,1} and mu_{+,2}; D_- has the opposite columns.
    for endpoint_sign, column, direction in [(-1, 0, 0.8 + 0.6j), (1, 1, 0.8 + 0.6j), (1, 0, 0.8 - 0.6j), (-1, 1, 0.8 - 0.6j)]:
        j = phase_integral(endpoint_sign)
        a = -1j * b + 1j * j * SIGMA
        small_constant = np.array([[q, -j], [j, q.conjugate()]], dtype=complex) / Q0
        large_errors = []
        small_errors = []
        for radius in [12.0, 24.0, 48.0]:
            z = radius * direction
            large = jost_column(z, endpoint_sign, column)
            large_errors.append(float(np.linalg.norm(large - EYE[:, column] - a[:, column] / z)))
            small_z = -Q0**2 / z
            # Sheet change swaps the column and its analytic domain.
            small = jost_column(small_z, endpoint_sign, 1 - column)
            prediction = -1j * Q0 * SIGMA1[:, 1 - column] / small_z + small_constant[:, 1 - column]
            small_errors.append(float(np.linalg.norm(small - prediction)))
        large_orders = np.log2(np.array(large_errors[:-1]) / np.array(large_errors[1:]))
        small_orders = np.log2(np.array(small_errors[:-1]) / np.array(small_errors[1:]))
        assert min(large_orders) > 1.7, (endpoint_sign, column, large_errors)
        assert min(small_orders) > 0.7, (endpoint_sign, column, small_errors)
        rows.append({"endpoint": endpoint_sign, "column": column + 1, "large_errors": large_errors, "large_orders": large_orders.tolist(), "small_errors": small_errors, "small_orders": small_orders.tolist()})
    return {"potential": "1 + 0.22 exp(-(x-0.27)^2) (1 + i (0.11+0.35x))", "radii": [12, 24, 48], "columns": rows}


def branch_checks() -> dict[str, float]:
    max_identity_error = 0.0
    endpoint_errors = []
    for sign in [-1, 1]:
        for eps in [0.1, 0.01, 0.001]:
            z = sign * 1j * (Q0 + eps)
            k, lam = spectral(z)
            ub = -1j * k * SIGMA + qmatrix(Q0)
            r = float(sign)  # Im(lambda) * r > 0 for column 1.
            e = np.exp(2j * lam * r)
            kernel = (1 + e) * EYE / 2 + np.expm1(2j * lam * r) * ub / (2j * lam)
            exact = np.exp(1j * lam * r) * expm(r * ub)
            max_identity_error = max(max_identity_error, float(np.linalg.norm(kernel - exact)))
            assert np.linalg.norm(kernel, 1) <= 1 + abs(r) * np.linalg.norm(ub, 1) + 1e-12
        zb = sign * 1j * Q0
        branch = jost_column(zb, -1, 0)
        # Approach from the admissible closed domain for this negative-end column.
        for eps in [0.01, 0.001, 0.0001]:
            z = sign * 1j * Q0 + 1j * eps
            endpoint_errors.append(float(np.linalg.norm(jost_column(z, -1, 0) - branch)))
        assert endpoint_errors[-1] < endpoint_errors[-2] < endpoint_errors[-3]
    assert max_identity_error < 1e-12
    return {"max_kernel_identity_error": max_identity_error, "max_last_branch_column_error": max(endpoint_errors[2::3])}


def scattering_checks() -> dict[str, object]:
    jtotal = phase_integral(-1) - phase_integral(1)
    large_errors = []
    small_errors = []
    for radius in [12.0, 24.0, 48.0]:
        for z, errors in [(radius, large_errors), (-Q0**2 / radius, small_errors)]:
            minus = np.column_stack([jost_column(z, -1, j) for j in range(2)])
            plus = np.column_stack([jost_column(z, 1, j) for j in range(2)])
            scattering = np.linalg.solve(minus, plus)  # x=t=0 makes theta=0.
            first_order = -1j * jtotal * SIGMA / z if abs(z) > 1 else -1j * jtotal * z * SIGMA / Q0**2
            errors.append(float(np.linalg.norm(scattering - EYE - first_order)))
            assert abs(np.linalg.det(scattering) - 1) < 1e-9
    large_orders = np.log2(np.array(large_errors[:-1]) / np.array(large_errors[1:]))
    small_orders = np.log2(np.array(small_errors[:-1]) / np.array(small_errors[1:]))
    assert min(large_orders) > 1.7
    assert min(small_orders) > 1.7
    return {"J_total": jtotal, "large_errors": large_errors, "large_orders": large_orders.tolist(), "small_errors": small_errors, "small_orders": small_orders.tolist()}


if __name__ == "__main__":
    print(json.dumps({"symbolic": symbolic_checks(), "branch": branch_checks(), "asymptotic": asymptotic_checks(), "scattering": scattering_checks()}, ensure_ascii=False, indent=2))
