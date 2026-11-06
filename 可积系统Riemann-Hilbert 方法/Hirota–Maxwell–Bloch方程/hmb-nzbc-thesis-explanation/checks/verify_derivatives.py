"""Independently check the added projection/field derivative calculations.

Run from any directory:
    python checks/verify_derivatives.py

Requires SymPy. The default JSON report is written only inside this explanation
project. No original-thesis files, TeX files, or long numerical checks are run.

T and Tc represent T and its complex conjugate. They are algebraically
independent during simplification, as are F and Fc. The derivative operators
encode the actual exponential evolution; checking identities with independent
symbols avoids relying on a selected parameter sample.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "chapters" / "09-derivative-details.tex"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "calculations" / "derivative-checks.json",
    )
    args = parser.parse_args()
    checks: list[dict[str, object]] = []

    def residual(value):
        return sp.simplify(sp.cancel(sp.together(value)))

    def check(name: str, value) -> None:
        if isinstance(value, sp.MatrixBase):
            reduced = value.applyfunc(residual)
            passed = all(entry == 0 for entry in reduced)
            detail = [[str(entry) for entry in reduced.row(j)]
                      for j in range(reduced.rows)]
        else:
            reduced = residual(value)
            passed = reduced == 0
            detail = str(reduced)
        checks.append({"name": name, "passed": bool(passed),
                       "residual": detail})

    h, T, Tc, kappa, F, Fc = sp.symbols(
        "h_0 T Tc kappa F Fc", nonzero=True,
    )
    a, b, q0, eta0, Omega = sp.symbols(
        "a b q_0 eta_0 Omega", real=True, nonzero=True,
    )
    phase = sp.symbols("phase", nonzero=True)
    delta = (1 + h**2) * (1 + T * Tc) - 2 * h * (T + Tc)
    numerator = (T - h) * (1 - h * Tc)
    P = sp.Matrix([
        [(T - h) * (Tc - h), numerator],
        [(1 - h * T) * (Tc - h), (1 - h * T) * (1 - h * Tc)],
    ]) / delta

    def dx(expr):
        return 2 * kappa * (T * sp.diff(expr, T)
                            + Tc * sp.diff(expr, Tc))

    def dt(expr):
        return (-2 * kappa * (F * T * sp.diff(expr, T)
                              + Fc * Tc * sp.diff(expr, Tc))
                - sp.I * Omega * phase * sp.diff(expr, phase))

    delta_x = 4 * kappa * (1 + h**2) * T * Tc - 4 * kappa * h * (T + Tc)
    delta_xx = 16 * kappa**2 * (1 + h**2) * T * Tc - 8 * kappa**2 * h * (T + Tc)
    delta_xxx = 64 * kappa**3 * (1 + h**2) * T * Tc - 16 * kappa**3 * h * (T + Tc)
    delta_t = (-2 * kappa * (1 + h**2) * (F + Fc) * T * Tc
               + 4 * kappa * h * (F * T + Fc * Tc))
    check("Delta_x", dx(delta) - delta_x)
    check("Delta_xx", dx(dx(delta)) - delta_xx)
    check("Delta_xxx", dx(dx(dx(delta))) - delta_xxx)
    check("Delta_t", dt(delta) - delta_t)
    check("absT2_x", dx(T * Tc) - 4 * kappa * T * Tc)
    check("absT2_xx", dx(dx(T * Tc)) - 16 * kappa**2 * T * Tc)
    check("absT2_xxx", dx(dx(dx(T * Tc))) - 64 * kappa**3 * T * Tc)
    check("absT2_t", dt(T * Tc) + 2 * kappa * (F + Fc) * T * Tc)

    numerator_x = 2 * kappa * (T + h**2 * Tc - 2 * h * T * Tc)
    numerator_xx = 4 * kappa**2 * (T + h**2 * Tc - 4 * h * T * Tc)
    numerator_xxx = 8 * kappa**3 * (T + h**2 * Tc - 8 * h * T * Tc)
    numerator_t = -2 * kappa * (F * T + h**2 * Fc * Tc
                                - h * (F + Fc) * T * Tc)
    check("P12_numerator_x", dx(numerator) - numerator_x)
    check("P12_numerator_xx", dx(dx(numerator)) - numerator_xx)
    check("P12_numerator_xxx", dx(dx(dx(numerator))) - numerator_xxx)
    check("P12_numerator_t", dt(numerator) - numerator_t)

    p12_x = numerator_x / delta - numerator * delta_x / delta**2
    p12_xx = (numerator_xx / delta
              - 2 * numerator_x * delta_x / delta**2
              - numerator * delta_xx / delta**2
              + 2 * numerator * delta_x**2 / delta**3)
    p12_xxx = (numerator_xxx / delta
               - 3 * numerator_xx * delta_x / delta**2
               - 3 * numerator_x * delta_xx / delta**2
               + 6 * numerator_x * delta_x**2 / delta**3
               - numerator * delta_xxx / delta**2
               + 6 * numerator * delta_x * delta_xx / delta**3
               - 6 * numerator * delta_x**3 / delta**4)
    p12_t = numerator_t / delta - numerator * delta_t / delta**2
    check("P12_x_expanded", dx(P[0, 1]) - p12_x)
    check("P12_xx_expanded", dx(dx(P[0, 1])) - p12_xx)
    check("P12_xxx_expanded", dx(dx(dx(P[0, 1]))) - p12_xxx)
    check("P12_t_expanded", dt(P[0, 1]) - p12_t)
    check("P12_x_product_recurrence",
          delta * p12_x + delta_x * P[0, 1] - numerator_x)
    check("P12_xx_product_recurrence",
          delta * p12_xx + 2 * delta_x * p12_x
          + delta_xx * P[0, 1] - numerator_xx)
    check("P12_xxx_product_recurrence",
          delta * p12_xxx + 3 * delta_x * p12_xx
          + 3 * delta_xx * p12_x + delta_xxx * P[0, 1] - numerator_xxx)

    other_numerators = {
        (0, 0): (
            (T - h) * (Tc - h),
            4 * kappa * T * Tc - 2 * kappa * h * (T + Tc),
            -2 * kappa * (F + Fc) * T * Tc
            + 2 * kappa * h * (F * T + Fc * Tc),
        ),
        (1, 0): (
            (1 - h * T) * (Tc - h),
            2 * kappa * (Tc + h**2 * T - 2 * h * T * Tc),
            -2 * kappa * (Fc * Tc + h**2 * F * T
                          - h * (F + Fc) * T * Tc),
        ),
        (1, 1): (
            (1 - h * T) * (1 - h * Tc),
            -2 * kappa * h * (T + Tc) + 4 * kappa * h**2 * T * Tc,
            2 * kappa * h * (F * T + Fc * Tc)
            - 2 * kappa * h**2 * (F + Fc) * T * Tc,
        ),
    }
    for (j, k), (num, num_x, num_t) in other_numerators.items():
        name = f"P{j + 1}{k + 1}"
        check(name + "_numerator_x", dx(num) - num_x)
        check(name + "_numerator_t", dt(num) - num_t)
        check(name + "_x", dx(P[j, k])
              - (num_x / delta - num * delta_x / delta**2))
        check(name + "_t", dt(P[j, k])
              - (num_t / delta - num * delta_t / delta**2))
    check("trace_P_x", dx(P[0, 0]) + dx(P[1, 1]))
    check("trace_P_t", dt(P[0, 0]) + dt(P[1, 1]))
    conjugate_formal = lambda expr: expr.xreplace({T: Tc, Tc: T, F: Fc, Fc: F})
    check("P21_x_conjugation", dx(P[1, 0]) - conjugate_formal(p12_x))
    check("P21_t_conjugation", dt(P[1, 0]) - conjugate_formal(p12_t))

    qtilde = -phase * (q0 + 4 * b * P[0, 1])
    check("qtilde_x", dx(qtilde) + 4 * b * phase * p12_x)
    check("qtilde_xx", dx(dx(qtilde)) + 4 * b * phase * p12_xx)
    check("qtilde_xxx", dx(dx(dx(qtilde))) + 4 * b * phase * p12_xxx)
    check("qtilde_t", dt(qtilde)
          - phase * (sp.I * Omega * (q0 + 4 * b * P[0, 1]) - 4 * b * p12_t))

    N = (a + sp.I * b) * sp.eye(2) - 2 * sp.I * b * P
    Nd = (a - sp.I * b) * sp.eye(2) + 2 * sp.I * b * P
    Nx = N.applyfunc(dx)
    Ndx = Nd.applyfunc(dx)
    check("N_x_all_entries", Nx + 2 * sp.I * b * P.applyfunc(dx))
    check("N_dagger_x_all_entries", Ndx - 2 * sp.I * b * P.applyfunc(dx))
    H = sp.diag(sp.I, -sp.I)
    Rbg = eta0 * sp.Matrix([[1, sp.I * q0 / a], [-sp.I * q0 / a, -1]])
    medium = N * Rbg * Nd
    medium_x = Nx * Rbg * Nd + N * Rbg * Ndx
    check("medium_product_derivative_all_entries", medium.applyfunc(dx) - medium_x)
    check("medium_H_diagonal_sign", (H * medium_x * H.inv())[0, 0] - medium_x[0, 0])
    check("medium_H_offdiagonal_sign", (H * medium_x * H.inv())[0, 1] + medium_x[0, 1])

    # Differentiate the printed component formulas with independent matrix
    # entries, then substitute n_x = -2 i b P_x and its conjugate.
    n11, n12, n21, n22, c11, c12, c21, c22 = sp.symbols(
        "n11 n12 n21 n22 nc11 nc12 nc21 nc22",
    )
    px11, px12, px21, px22, pc11, pc12, pc21, pc22 = sp.symbols(
        "Px11 Px12 Px21 Px22 Pxc11 Pxc12 Pxc21 Pxc22",
    )
    entries = [n11, n12, n21, n22, c11, c12, c21, c22]
    rates = [-2 * sp.I * b * px11, -2 * sp.I * b * px12,
             -2 * sp.I * b * px21, -2 * sp.I * b * px22,
             2 * sp.I * b * pc11, 2 * sp.I * b * pc12,
             2 * sp.I * b * pc21, 2 * sp.I * b * pc22]
    differentiate_entries = lambda expr: sum(
        sp.diff(expr, var) * rate for var, rate in zip(entries, rates)
    )
    eta_component = n11 * c11 - n12 * c12 + sp.I * q0 / a * (n11 * c12 - n12 * c11)
    pol_component = n11 * c21 - n12 * c22 + sp.I * q0 / a * (n11 * c22 - n12 * c21)
    eta_x_printed = 2 * sp.I * b * (
        -px11 * c11 + n11 * pc11 + px12 * c12 - n12 * pc12
        + sp.I * q0 / a * (-px11 * c12 + n11 * pc12 + px12 * c11 - n12 * pc11)
    )
    pol_x_printed = 2 * sp.I * b * (
        -px11 * c21 + n11 * pc21 + px12 * c22 - n12 * pc22
        + sp.I * q0 / a * (-px11 * c22 + n11 * pc22 + px12 * c21 - n12 * pc21)
    )
    check("eta_x_expanded_component_substitution",
          differentiate_entries(eta_component) - eta_x_printed)
    check("ptilde_x_expanded_component_substitution",
          differentiate_entries(pol_component) - pol_x_printed)

    Pminus = sp.Matrix([[h**2, -h], [-h, 1]]) / (1 + h**2)
    Pplus = sp.Matrix([[1, -h], [-h, h**2]]) / (1 + h**2)
    difference_minus = sp.Matrix([
        [(1 - h**2) * (T * ((1 + h**2) * Tc - h) - h * Tc),
         (1 - h**2) * (T - h**2 * Tc)],
        [(1 - h**2) * (Tc - h**2 * T),
         -(1 - h**2) * (T * ((1 + h**2) * Tc - h) - h * Tc)],
    ]) / ((1 + h**2) * delta)
    difference_plus = sp.Matrix([
        [(1 - h**2) * (h * (T + Tc) - (1 + h**2)),
         (1 - h**2) * (T - h**2 * Tc)],
        [(1 - h**2) * (Tc - h**2 * T),
         -(1 - h**2) * (h * (T + Tc) - (1 + h**2))],
    ]) / ((1 + h**2) * delta)
    check("P_minus_Pminus_all_entries", P - Pminus - difference_minus)
    check("P_minus_Pplus_all_entries", P - Pplus - difference_plus)
    scaled_delta = (1 + h**2) * (1 + 1 / (T * Tc)) - 2 * h * (1 / T + 1 / Tc)
    check("positive_end_scaled_denominator", delta / (T * Tc) - scaled_delta)
    scaled_plus = sp.Matrix([
        [(1 - h**2) * (h * (1 / T + 1 / Tc) - (1 + h**2) / (T * Tc)),
         (1 - h**2) * (1 / Tc - h**2 / T)],
        [(1 - h**2) * (1 / T - h**2 / Tc),
         -(1 - h**2) * (h * (1 / T + 1 / Tc) - (1 + h**2) / (T * Tc))],
    ]) / ((1 + h**2) * scaled_delta)
    check("positive_end_scaled_differences_all_entries", P - Pplus - scaled_plus)

    # Endpoint reconstruction uses the original b,h_0,q_0 relation.
    relation = {b: q0 * (1 + h**2) / (2 * h)}
    for name, Pend in [("minus", Pminus), ("plus", Pplus)]:
        Nend = (a + sp.I * b) * sp.eye(2) - 2 * sp.I * b * Pend
        Nend_dagger = (a - sp.I * b) * sp.eye(2) + 2 * sp.I * b * Pend
        check("medium_endpoint_" + name,
              (H * Nend * Rbg * Nend_dagger * H.inv()
               / (a**2 + b**2) - Rbg).subs(relation))
        product_difference = ((N - Nend) * Rbg * Nd
                              + Nend * Rbg * (Nd - Nend_dagger))
        check("medium_product_difference_" + name,
              N * Rbg * Nd - Nend * Rbg * Nend_dagger - product_difference)
        error_printed = H * (
            -2 * sp.I * b * (P - Pend) * Rbg * Nd
            + 2 * sp.I * b * Nend * Rbg * (P - Pend)
        ) * H.inv() / (a**2 + b**2)
        check("medium_endpoint_error_formula_" + name,
              (H * N * Rbg * Nd * H.inv() / (a**2 + b**2)
               - Rbg - error_printed).subs(relation))
        check("qtilde_endpoint_error_" + name,
              (qtilde - q0 * phase
               + 4 * b * phase * (P[0, 1] - Pend[0, 1])).subs(relation))

    h_example = (3 - sp.sqrt(5)) / 2
    kappa_example = sp.sqrt(5) / 2
    c1_example = 4 * sp.I * kappa_example / (1 + h_example**2)
    c2_example = -h_example**2 * c1_example
    check("example_h_squared", h_example**2 - (7 - 3 * sp.sqrt(5)) / 2)
    check("example_1_plus_h_squared", 1 + h_example**2 - 3 * h_example)
    check("example_c1_radical", c1_example - sp.I * (5 + 3 * sp.sqrt(5)) / 3)
    check("example_c2_radical", c2_example - sp.I * (5 - 3 * sp.sqrt(5)) / 3)
    for name, exact, printed in [
        ("example_c1_decimal", sp.im(c1_example), "3.902734644166456"),
        ("example_c2_decimal", sp.im(c2_example), "-0.569401310833123"),
    ]:
        error = sp.N(abs(exact - sp.Rational(printed)), 40)
        checks.append({"name": name, "passed": bool(error < sp.Rational(1, 2 * 10**15)),
                       "absolute_error": str(error), "tolerance": "5e-16"})

    passed = sum(item["passed"] for item in checks)
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "checks/verify_derivatives.py",
        "source": "chapters/09-derivative-details.tex",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "python_version": sys.version.split()[0],
        "sympy_version": sp.__version__,
        "method": "Exact symbolic identities with independent T,Tc,F,Fc; no finite differences.",
        "scope": ["projection derivatives and recurrences", "physical optical-field derivatives",
                  "N and medium derivatives", "component substitutions", "both endpoint errors",
                  "norming-constant radicals and printed decimals"],
        "checks_total": len(checks),
        "checks_passed": passed,
        "all_passed": passed == len(checks),
        "checks": checks,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Derivative and endpoint symbolic checks: {passed}/{len(checks)} passed.")
    print(f"Report: {args.report.resolve()}")
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
