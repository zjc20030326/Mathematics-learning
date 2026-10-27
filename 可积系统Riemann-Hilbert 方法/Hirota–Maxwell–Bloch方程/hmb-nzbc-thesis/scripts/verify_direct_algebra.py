"""Exact, parameter-generic elimination of the one-orbit residue system.

t and r are algebraically independent; on the physical slice r=conjugate(t).
The solution is obtained by matrix inversion, not from a projector ansatz.
"""
from pathlib import Path
import json
import sympy as s

h, Q, t, r, w = s.symbols('h Q t r w', nonzero=True)
i = s.I
b = Q*(1+h*h)/(2*h)
k = Q*(1-h*h)/(2*h)
zeta = [i*Q/h, -i*Q*h]
zbar = [-zz for zz in zeta]
c = 4*i*k/(1+h*h)
C = [c/t, -h*h*c/r]
Cb = [-c/r, h*h*c/t]
A = s.Matrix(2, 2, lambda j, l: Cb[l]/(zeta[j]-zbar[l]))
B = s.Matrix(2, 2, lambda j, l: C[l]/(zbar[j]-zeta[l]))
f = s.Matrix([[-i*Q/zz, 1] for zz in zeta])
g = s.Matrix([[1, -i*Q/zz] for zz in zbar])
L = s.eye(2)+A*B
v = (L.inv()*(f-A*g)).applyfunc(s.factor)
u = (g+B*v).applyfunc(s.factor)
dm = (1+h*h)*(t*r+1)-2*h*(t+r)
dp = (1+h*h)*(t*r+1)+2*h*(t+r)
assert s.factor(L.det()-dm*dp/(t*t*r*r*(1+h*h)**2)) == 0
v_expected = s.Matrix([[t*(1+h*h)*(1-h*r), t*(1+h*h)*(r-h)],
                       [r*(1+h*h)*(t-h)/h, r*(1+h*h)*(h*t-1)/h]])/dm
u_expected = s.Matrix([[r*(1+h*h)*(t-h), r*(1+h*h)*(h*t-1)],
                       [t*(1+h*h)*(h*r-1)/h, -t*(1+h*h)*(r-h)/h]])/dm
assert (v-v_expected).applyfunc(s.cancel) == s.zeros(2)
assert (u-u_expected).applyfunc(s.cancel) == s.zeros(2)
assert (v+A*u-f).applyfunc(s.cancel) == s.zeros(2)
assert (u-B*v-g).applyfunc(s.cancel) == s.zeros(2)

sx = s.Matrix([[0, 1], [1, 0]])
H = s.diag(i, -i)
E = s.eye(2)-i*Q/w*sx
M = E.copy()
for j in range(2):
    M[:, 0] += C[j]*v[j, :].T/(w-zeta[j])
    M[:, 1] -= Cb[j]*u[j, :].T/(w-zbar[j])
M = M.applyfunc(s.factor)
# Only after the residue system is solved, factor its resulting matrix.
y = s.Matrix([t-h, 1-h*t])
yb = s.Matrix([[r-h, 1-h*r]])
P = y*yb/dm
lam = (w-Q*Q/w)/2
rho = (w+Q*Q/w)/2
N = (lam+i*b)*s.eye(2)-2*i*b*P
compact = H*N*E*s.diag(1/(rho-i*k), 1/(rho+i*k))*H**-1
assert (M-compact).applyfunc(lambda v: s.simplify(s.factor(v))) == s.zeros(2)
assert (P*P-P).applyfunc(s.factor) == s.zeros(2)
assert s.factor(N.det()-(rho*rho+k*k)) == 0
assert (M.subs(w, -Q*Q/w)-M*(i*w/Q)*sx).applyfunc(lambda v: s.simplify(s.factor(v))) == s.zeros(2)
assert (M.applyfunc(lambda v: s.limit(v, w, s.oo))-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)
assert ((w*M).applyfunc(lambda v: s.limit(v, w, 0))+i*Q*sx).applyfunc(s.simplify) == s.zeros(2)
q_from_infinity = i*s.limit(w*M[0, 1], w, s.oo)
assert s.factor(q_from_infinity+Q+4*b*P[0, 1]) == 0
report = {
    'variables': '0<h<1, Q>0; r=conjugate(t), t!=0 on the physical slice',
    'Schur_determinant': str(s.factor(L.det())),
    'computed_v_rows': str(v), 'computed_u_rows': str(u),
    'checks': {name: 'exact zero' for name in (
        'Schur factorization', 'explicit v', 'explicit u', 'finite1', 'finite2',
        'partial fractions equal compact M', 'rank-one projector identity',
        'numerator determinant', 'sheet symmetry', 'infinity normalization',
        'zero normalization', 'optical reconstruction')},
}
out = Path(__file__).resolve().parents[1]/'calculations'/'direct-algebra-results.json'
out.write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
