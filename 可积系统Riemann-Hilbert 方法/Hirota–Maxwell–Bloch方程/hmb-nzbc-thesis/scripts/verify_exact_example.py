"""Exact rational PDE verification from the one-orbit residue linear system.

q0=1, b=5/4, kappa=3/4, a=2, eta0=3/5, alpha=1, beta=1/10.
tau=exp(2 Re(s))>0, C=exp(2 i Im(s)); conjugation sends C to 1/C.
No Darboux or projector input is used.
"""
from pathlib import Path
import json
import sympy as s

tau = s.Symbol('tau', positive=True)
C = s.Symbol('C', nonzero=True)
I, rat = s.I, s.Rational
q0, b, kappa, a, eta0 = 1, rat(5, 4), rat(3, 4), 2, rat(3, 5)
alpha, beta = 1, rat(1, 10)
Omega = -alpha+2*eta0/a
k = I*b
F = s.simplify(s.expand_complex(4*beta*k**2-alpha*k-2*beta+eta0/(a*(k-a))))
reF, imF = s.re(F), s.im(F)

# First solve finite1/finite2. Rows are the transposed regular columns.
zeta = [2*I, -I/2]
zbar = [-zz for zz in zeta]
c1 = rat(12, 5)*I
coeff = [c1*tau*C, -c1*tau/(4*C)]
coeff_bar = [-c1*tau/C, c1*tau*C/4]
A = s.Matrix(2, 2, lambda j, l: coeff_bar[l]/(zeta[j]-zbar[l]))
B = s.Matrix(2, 2, lambda j, l: coeff[l]/(zbar[j]-zeta[l]))
f = s.Matrix([[1, -I/zz] for zz in zeta])
g = s.Matrix([[-I/zz, 1] for zz in zbar])
L = s.eye(2)+A*B
regular_v = (L.inv()*(f-A*g)).applyfunc(s.cancel)
regular_u = (g+B*regular_v).applyfunc(s.cancel)
assert (regular_v+A*regular_u-f).applyfunc(s.cancel) == s.zeros(2)
assert (regular_u-B*regular_v-g).applyfunc(s.cancel) == s.zeros(2)

# The light field comes from the z^{-1} coefficient, the medium from resonance.
q = s.cancel(1+I*sum(coeff[j]*regular_v[j, 0] for j in range(2)))
za = 2+s.sqrt(5)
lambda_a = s.sqrt(5)
M = s.Matrix([[1, -I/za], [-I/za, 1]])
for j in range(2):
    M[:, 1] += coeff[j]*regular_v[j, :].T/(za-zeta[j])
    M[:, 0] -= coeff_bar[j]*regular_u[j, :].T/(za-zbar[j])
M = M.applyfunc(lambda v: s.cancel(s.radsimp(v)))
gamma = 1+1/za**2
R = (M*s.diag(eta0*lambda_a/a, -eta0*lambda_a/a)*M.adjugate()/gamma).applyfunc(lambda v: s.cancel(s.radsimp(v)))
p, eta = -R[0, 1], R[0, 0]
def conj(f):
    return s.cancel(s.conjugate(f).subs(s.conjugate(C), 1/C))
qc, pc = conj(q), conj(p)
assert s.cancel(conj(eta)-eta) == 0
assert s.cancel(R[1, 0]+pc) == 0

def dx(f):
    return s.cancel(2*kappa*tau*s.diff(f, tau))

def dt(f):
    return s.cancel(-2*kappa*reF*tau*s.diff(f, tau)-2*I*kappa*imF*C*s.diff(f, C))

qx, qxx = dx(q), dx(dx(q))
qxxx = dx(qxx)
r1 = s.cancel(dt(q)-I*Omega*q-I*alpha*(qxx/2+q*qc*q)-beta*(qxxx+6*q*qc*qx)-2*p)
r2 = s.cancel(dx(p)-2*eta*q+2*I*a*p)
r3 = s.cancel(dx(eta)+qc*p+q*pc)
rb = s.cancel(p*pc+eta**2-eta0**2*(1+rat(1, 4)))
assert [r1, r2, r3, rb] == [0, 0, 0, 0]
report = dict(method='Exact Schur solve of residue system; infinity/resonance reconstruction',
    parameters={'q0':'1', 'b':'5/4', 'kappa':'3/4', 'a':'2', 'eta0':'3/5', 'alpha':'1', 'beta':'1/10'},
    F=str(F), variables='tau=exp(2 Re(s)), C=exp(2 i Im(s)), T=tau*C, s=kappa*(x-F*t)',
    residue_system_residuals=['0', '0'], medium_Hermitian_residual='0',
    exact_residuals_q_p_eta_bloch=list(map(str, [r1, r2, r3, rb])))
out = Path(__file__).resolve().parents[1]/'calculations'/'exact-example-results.json'
out.write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
