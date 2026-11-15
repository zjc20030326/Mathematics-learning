"""Exact matrix checks for the background and the Jost differential identity."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
q0 = s.symbols('q0', positive=True)
z = s.symbols('z', nonzero=True)
a, eta0 = s.symbols('a eta0', real=True, nonzero=True)
alpha, beta, x, t = s.symbols('alpha beta x t', real=True)
sig = s.diag(1, -1)
eye = s.eye(2)
k = (z-q0**2/z)/2
lam = (z+q0**2/z)/2
F = 4*beta*k**2-alpha*k-2*beta*q0**2+eta0/(a*(k-a))
Q = s.Matrix([[0, q0], [-q0, 0]])
U = -s.I*k*sig+Q
V = -F*U
Y = s.Matrix([[1, -s.I*q0/z], [-s.I*q0/z, 1]])
gamma = 1+q0**2/z**2
Yi = s.Matrix([[1, s.I*q0/z], [s.I*q0/z, 1]])/gamma
theta = -lam*(x-F*t)
checks = []


def check(name, expression):
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residuals = [s.factor(s.cancel(entry)) for entry in entries]
    passed = all(entry == 0 for entry in residuals)
    checks.append({'name': name, 'passed': passed,
                   'residuals': [str(entry) for entry in residuals]})
    assert passed, (name, residuals)


check('spectral_curve', lam**2-k**2-q0**2)
check('background_square', U*U+(k**2+q0**2)*eye)
check('first_eigenvector', U*Y[:, 0]+s.I*lam*Y[:, 0])
check('second_eigenvector', U*Y[:, 1]-s.I*lam*Y[:, 1])
check('determinant', Y.det()-gamma)
check('inverse', Y*Yi-eye)
check('space_diagonalization', Yi*U*Y+s.I*lam*sig)
check('time_diagonalization', Yi*V*Y-s.I*lam*F*sig)
check('theta_x', s.diff(theta, x)+lam)
check('theta_t', s.diff(theta, t)-lam*F)

# Arbitrary perturbation matrices and arbitrary values of the modified solution.
# This checks matrix multiplication order without imposing an HMB solution.
mu = s.Matrix(2, 2, s.symbols('m11 m12 m21 m22'))
dQ = s.Matrix([[0, s.Symbol('dq12')], [s.Symbol('dq21'), 0]])
dV = s.Matrix(2, 2, s.symbols('dv11 dv12 dv21 dv22'))
mu_x = (U+dQ)*mu+s.I*lam*mu*sig
mu_t = (V+dV)*mu-s.I*lam*F*mu*sig
comm = sig*Yi*mu-Yi*mu*sig
check('x_coefficient_of_total_differential',
      Yi*mu_x-s.I*s.diff(theta, x)*comm-Yi*dQ*mu)
check('t_coefficient_of_total_differential',
      Yi*mu_t-s.I*s.diff(theta, t)*comm-Yi*dV*mu)
Rbg = s.Matrix([[eta0, s.I*eta0*q0/a], [-s.I*eta0*q0/a, -eta0]])
check('medium_background_identity', Rbg-s.I*eta0/a*(-s.I*a*sig+Q))

sources = [ROOT/'section'/name for name in
           ['02-model.tex', '03-direct-scattering.tex', '04-evolution-resonance.tex']]
report = {'generated_utc': datetime.now(timezone.utc).isoformat(),
          'script': 'scripts/verify_background_jost.py',
          'method': 'Exact symbolic matrix identities; no numerical sampling.',
          'domain': 'z != 0, +/- i*q0; k(z) != a for time expressions.',
          'scope': 'Algebraic identities for conditional full-Lax solutions; no inverse-solvability claim.',
          'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in sources},
          'checks_total': len(checks), 'checks_passed': sum(c['passed'] for c in checks),
          'checks': checks}
(ROOT/'calculations/background-jost-identities.json').write_text(
    json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(f"Background and Jost exact identities: {len(checks)}/{len(checks)} passed.")
