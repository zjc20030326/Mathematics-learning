"""Exact algebraic checks for normalized Jost matrices and scattering phases."""
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
started = perf_counter()
sig = s.diag(1, -1)
ell, F = s.symbols('lambda F')
phase, tau = s.symbols('phase tau', nonzero=True)

# At a fixed spectral point, phase=exp(i*lambda*x) and
# tau=exp(i*lambda*F*t).  Nonzero algebraic symbols keep every residual
# rational; their derivatives are imposed only by the chain rule.
E = s.diag(phase**-1, phase)
Ei = s.diag(phase, phase**-1)
T = s.diag(tau/phase, phase/tau)  # exp(i*theta*sigma_3)
M = s.Matrix(2, 2, s.symbols('m11 m12 m21 m22'))
Mi = M.inv()
U = s.Matrix(2, 2, s.symbols('u11 u12 u21 u22'))
V = s.Matrix(2, 2, s.symbols('v11 v12 v21 v22'))
R = s.Matrix(2, 2, s.symbols('r11 r12 r21 r22'))
S = s.Matrix(2, 2, s.symbols('s11 s12 s21 s22'))
checks = []


def comm(a, b):
    return a*b-b*a


def check(name, expression):
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residuals = [s.factor(s.cancel(entry)) for entry in entries]
    passed = all(entry == 0 for entry in residuals)
    checks.append({'name': name, 'passed': passed,
                   'residuals': [str(entry) for entry in residuals]})
    assert passed, (name, residuals)


# Product-rule restoration of phi=M*exp(i*theta*sigma_3).
Mx = U*M+s.I*ell*M*sig
Mt = V*M-s.I*ell*F*M*sig
Tx = s.diff(T, phase)*s.I*ell*phase
Tt = s.diff(T, tau)*s.I*ell*F*tau
check('modified_space_equation_restores_phi_lax', Mx*T+M*Tx-U*M*T)
check('modified_time_equation_restores_phi_lax', Mt*T+M*Tt-V*M*T)

# The plus/minus relation is M_plus=M_minus*E*S*E^{-1}.
A = E*S*Ei
P = M*A
Ax = s.diff(A, phase)*s.I*ell*phase
check('spatial_scattering_exponential_conjugation',
      A-s.Matrix([[S[0, 0], phase**-2*S[0, 1]],
                  [phase**2*S[1, 0], S[1, 1]]]))
check('scattering_space_constancy',
      U*P+s.I*ell*P*sig-Mx*A-M*Ax)
check('scattering_time_commutator',
      Ei*Mi*((V*P-s.I*ell*F*P*sig)-Mt*A)*E
      -comm(s.I*ell*F*sig, S))
time_conjugation = s.diag(tau, tau**-1)*S*s.diag(tau**-1, tau)
check('time_conjugation_solves_scattering_ode',
      s.diff(time_conjugation, tau)*s.I*ell*F*tau
      -comm(s.I*ell*F*sig, time_conjugation))

# Compute all four determinants from the actual matrix relation.  The
# off-diagonal formulas require opposite spatial phase factors.
gamma = M.det()
check('transmission_s11_wronskian',
      s.Matrix.hstack(P[:, 0], M[:, 1]).det()/gamma-S[0, 0])
check('transmission_s22_wronskian',
      s.Matrix.hstack(M[:, 0], P[:, 1]).det()/gamma-S[1, 1])
check('reflection_s21_wronskian_space_phase',
      phase**-2*s.Matrix.hstack(M[:, 0], P[:, 0]).det()/gamma-S[1, 0])
check('reflection_s12_wronskian_space_phase',
      phase**2*s.Matrix.hstack(P[:, 1], M[:, 1]).det()/gamma-S[0, 1])

# At resonance, R_x=[U,R].  B=M^{-1}*R*M has a commutator derivative,
# so E^{-1}*B*E is constant in x.  Endpoint arguments remain analytical.
B = Mi*R*M
Rx = comm(U, R)
Bx = -Mi*Mx*B+Mi*Rx*M+Mi*R*Mx
check('resonance_medium_conjugation_commutator', Bx+s.I*ell*comm(sig, B))
check('resonance_medium_exponential_conjugation_is_constant',
      s.diff(Ei, phase)*s.I*ell*phase*B*E+Ei*Bx*E
      +Ei*B*s.diff(E, phase)*s.I*ell*phase)

# On the real spectral axis, lambda is real and U is anti-Hermitian.
# The Gram matrix H=M^dagger*M obeys the same commutator equation.
ell_real, u_real = s.symbols('lambda_real u_real', real=True)
q = s.symbols('q')
Ureal = s.Matrix([[s.I*u_real, q], [-s.conjugate(q), -s.I*u_real]])
Mx_real = Ureal*M+s.I*ell_real*M*sig
H = M.conjugate().T*M
Hx = Mx_real.conjugate().T*M+M.conjugate().T*Mx_real
check('real_axis_gram_commutator', Hx+s.I*ell_real*comm(sig, H))
check('real_axis_gram_exponential_conjugation_is_constant',
      s.diff(Ei, phase)*s.I*ell_real*phase*H*E+Ei*Hx*E
      +Ei*H*s.diff(E, phase)*s.I*ell_real*phase)

# For mu_{+,2}=d*mu_{-,1}, d=b(t)*exp(-2*i*lambda*x).
# Derive the time coefficient from the two column equations, rather than
# assuming a space-independent proportionality factor for the mu columns.
b = s.symbols('b')
d = b/phase**2
column1 = s.Matrix(s.symbols('column11 column21'))
column2 = d*column1
dx = s.diff(d, phase)*s.I*ell*phase
column1_x = U*column1+s.I*ell*column1
column2_x = U*column2-s.I*ell*column2
column1_t = V*column1-s.I*ell*F*column1
column2_t = V*column2+s.I*ell*F*column2
check('discrete_proportion_space_sign', dx+2*s.I*ell*d)
check('discrete_proportion_space_column_equation',
      dx*column1+d*column1_x-column2_x)
check('discrete_norming_time_coefficient',
      phase**2*(column2_t-d*column1_t)-2*s.I*ell*F*b*column1)
check('discrete_full_lax_column_ratio',
      d*T[1, 1]-(b/tau**2)*T[0, 0])


# Spatial chi and modified mu use the same spatial scattering matrix.
# Here M and P still denote generic mu_minus and mu_plus values.
chi_minus = M*E
chi_plus = P*E
Ex = s.diff(E, phase)*s.I*ell*phase
check('spatial_chi_equation_from_mu',
      Mx*E+M*Ex-U*chi_minus)
check('spatial_chi_time_equation_from_mu',
      Mt*E-(V*chi_minus-chi_minus*s.I*ell*F*sig))
check('spatial_chi_scattering_equals_mu_scattering',
      chi_plus-chi_minus*S)
full_time_phase = s.diag(tau, tau**-1)
check('full_lax_chi_mu_normalizations_agree',
      chi_minus*full_time_phase-M*T)

# RHP data evolve by rho(t)=rho0/tau**2, rho_tilde(t)=rho_tilde0*tau**2.
# T=exp(i*theta*sigma_3) gives J=T*J0*T^{-1}; hence Psi=M_RHP*T
# has the x,t-independent jump J0.  This is a matrix-algebra identity,
# with time phases defined only where F is finite.
rho0, rho_tilde0, c0 = s.symbols('rho0 rho_tilde0 c0')
rho_t = rho0/tau**2
rho_tilde_t = rho_tilde0*tau**2
J = s.Matrix([[1, rho_tilde_t/phase**2],
              [-rho_t*phase**2, 1-rho_t*rho_tilde_t]])
J0 = s.Matrix([[1, rho_tilde0], [-rho0, 1-rho0*rho_tilde0]])
Ti = s.diag(phase/tau, tau/phase)
check('rhp_evolved_jump_conjugates_to_constant', Ti*J*T-J0)
c_t = c0*tau**2
C = c_t/phase**2
check('rhp_residue_timephase_consistency',
      s.Matrix([C-c0*T[0, 0]/T[1, 1], C*T[1, 1]/T[0, 0]-c0]))

sources = [ROOT/'section'/name for name in
           ['03-direct-scattering.tex', '04-evolution-resonance.tex',
            '05-rhp-reconstruction.tex', '06-reflectionless.tex',
            '07-hmb-verification.tex']]
report = {
    'generated_utc': datetime.now(timezone.utc).isoformat(),
    'script': 'scripts/verify_jost_normalization.py',
    'method': 'Exact rational matrix algebra with nonzero formal exponential phases and their chain-rule derivatives; no numerical sampling.',
    'domain': 'Generic invertible 2x2 mu values; nonzero spatial and time phases; lambda,F fixed at a spectral point. Time phases require finite F, hence k != a in HMB. Gram checks additionally assume real lambda and anti-Hermitian spatial U.',
    'scope': 'Checks spatial chi/modified mu/full phi normalization identities, scattering matrix order and evolution, four Wronskian phase factors, conjugation derivative identities, discrete column proportions and the constant jump/residue conjugation. Does not prove Jost existence, analyticity, endpoint convergence, resonance solvability, RHP solvability or an inverse problem.',
    'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in sources},
    'checks_total': len(checks),
    'checks_passed': sum(c['passed'] for c in checks),
    'elapsed_seconds': round(perf_counter()-started, 6),
    'checks': checks,
}
(ROOT/'calculations/jost-normalization-identities.json').write_text(
    json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(f"Jost normalization exact identities: {len(checks)}/{len(checks)} passed.")
