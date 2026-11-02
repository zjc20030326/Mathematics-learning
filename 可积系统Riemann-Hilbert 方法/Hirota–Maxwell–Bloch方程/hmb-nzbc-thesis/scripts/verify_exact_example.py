"""Exact rational PDE verification from the one-orbit residue linear system.

Q=1, b=5/4, kappa=3/4, a=2, eta0=3/5, alpha=1, beta=1/10.
T=exp(2 Re(s))>0, C=exp(2 i Im(s)); conjugation sends C to 1/C.
No Darboux or projector input is used.
"""
from pathlib import Path
import json
import sympy as s

T = s.Symbol('T', positive=True)
C = s.Symbol('C', nonzero=True)
I, rat = s.I, s.Rational
Q, b, k, a, eta0 = 1, rat(5, 4), rat(3, 4), 2, rat(3, 5)
alpha, beta = 1, rat(1, 10)
Omega = -alpha+2*eta0/a
lam = I*b
F = s.simplify(s.expand_complex(4*beta*lam**2-alpha*lam-2*beta+eta0/(a*(lam-a))))
reF, imF = s.re(F), s.im(F)

# First solve finite1/finite2. Rows are the transposed regular columns.
zeta = [2*I, -I/2]
zbar = [-zz for zz in zeta]
c1 = rat(12, 5)*I
coeff = [c1/(T*C), -c1*C/(4*T)]
coeff_bar = [-c1*C/T, c1/(4*T*C)]
A = s.Matrix(2, 2, lambda j, l: coeff_bar[l]/(zeta[j]-zbar[l]))
B = s.Matrix(2, 2, lambda j, l: coeff[l]/(zbar[j]-zeta[l]))
f = s.Matrix([[-I/zz, 1] for zz in zeta])
g = s.Matrix([[1, -I/zz] for zz in zbar])
L = s.eye(2)+A*B
regular_v = (L.inv()*(f-A*g)).applyfunc(s.cancel)
regular_u = (g+B*regular_v).applyfunc(s.cancel)
assert (regular_v+A*regular_u-f).applyfunc(s.cancel) == s.zeros(2)
assert (regular_u-B*regular_v-g).applyfunc(s.cancel) == s.zeros(2)

# The light field comes from the w^{-1} coefficient, the medium from resonance.
u = s.cancel(1-I*sum(coeff_bar[j]*regular_u[j, 0] for j in range(2)))
wa = 2+s.sqrt(5)
rhoa = s.sqrt(5)
M = s.Matrix([[1, -I/wa], [-I/wa, 1]])
for j in range(2):
    M[:, 0] += coeff[j]*regular_v[j, :].T/(wa-zeta[j])
    M[:, 1] -= coeff_bar[j]*regular_u[j, :].T/(wa-zbar[j])
M = M.applyfunc(lambda v: s.cancel(s.radsimp(v)))
gamma = 1+1/wa**2
R = (M*s.diag(eta0*rhoa/a, -eta0*rhoa/a)*M.adjugate()/gamma).applyfunc(lambda v: s.cancel(s.radsimp(v)))
v, e = -R[0, 1], R[0, 0]
def conj(f):
    return s.cancel(s.conjugate(f).subs(s.conjugate(C), 1/C))
ub, vb = conj(u), conj(v)
assert s.cancel(conj(e)-e) == 0
assert s.cancel(R[1, 0]+vb) == 0

def dx(f):
    return s.cancel(2*k*T*s.diff(f, T))

def dz(f):
    return s.cancel(-2*k*reF*T*s.diff(f, T)-2*I*k*imF*C*s.diff(f, C))

ux, uxx = dx(u), dx(dx(u))
uxxx = dx(uxx)
r1 = s.cancel(dz(u)-I*Omega*u-I*alpha*(uxx/2+u*ub*u)-beta*(uxxx+6*u*ub*ux)-2*v)
r2 = s.cancel(dx(v)-2*e*u+2*I*a*v)
r3 = s.cancel(dx(e)+ub*v+u*vb)
rb = s.cancel(v*vb+e**2-eta0**2*(1+rat(1, 4)))
assert [r1, r2, r3, rb] == [0, 0, 0, 0]
report = dict(method='Exact Schur solve of residue system; infinity/resonance reconstruction',
    parameters={'Q':'1', 'b':'5/4', 'kappa':'3/4', 'a':'2', 'eta0':'3/5', 'alpha':'1', 'beta':'1/10'},
    F=str(F), variables='T=exp(2 Re(s)), C=exp(2 i Im(s)), s=kappa*(x-F*z)',
    residue_system_residuals=['0', '0'], medium_Hermitian_residual='0',
    exact_residuals_q_p_eta_bloch=list(map(str, [r1, r2, r3, rb])))
out = Path(__file__).resolve().parents[1]/'calculations'/'exact-example-results.json'
out.write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
