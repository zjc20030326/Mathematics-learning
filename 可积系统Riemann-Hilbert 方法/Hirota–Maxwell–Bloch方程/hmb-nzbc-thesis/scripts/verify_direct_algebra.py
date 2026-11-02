"""Exact, parameter-generic elimination of the one-orbit residue system.

T and r are algebraically independent; on the physical slice r=conjugate(T).
The solution is obtained by matrix inversion, not from a projector ansatz.
"""
from pathlib import Path
import json
import sympy as s

h, q0, T, r, z = s.symbols('h q0 T r z', nonzero=True)
i = s.I
b = q0*(1+h*h)/(2*h)
kappa = q0*(1-h*h)/(2*h)
zeta = [i*q0/h, -i*q0*h]
zbar = [-zz for zz in zeta]
c = 4*i*kappa/(1+h*h)
C = [c*T, -h*h*c*r]
Cb = [-c*r, h*h*c*T]
A = s.Matrix(2, 2, lambda j, l: Cb[l]/(zeta[j]-zbar[l]))
B = s.Matrix(2, 2, lambda j, l: C[l]/(zbar[j]-zeta[l]))
f = s.Matrix([[1, -i*q0/zz] for zz in zeta])
g = s.Matrix([[-i*q0/zz, 1] for zz in zbar])
L = s.eye(2)+A*B
v = (L.inv()*(f-A*g)).applyfunc(s.factor)
u = (g+B*v).applyfunc(s.factor)
dm = (1+h*h)*(T*r+1)-2*h*(T+r)
dp = (1+h*h)*(T*r+1)+2*h*(T+r)
assert s.factor(L.det()-dm*dp/((1+h*h)**2)) == 0
v_expected = s.Matrix([[(1+h*h)*(1-h*r), (1+h*h)*(r-h)],
                       [(1+h*h)*(h-T)/h, (1+h*h)*(1-h*T)/h]])/dm
u_expected = s.Matrix([[(1+h*h)*(h-T), (1+h*h)*(1-h*T)],
                       [-(1+h*h)*(1-h*r)/h, -(1+h*h)*(r-h)/h]])/dm
assert (v-v_expected).applyfunc(s.cancel) == s.zeros(2)
assert (u-u_expected).applyfunc(s.cancel) == s.zeros(2)
assert (v+A*u-f).applyfunc(s.cancel) == s.zeros(2)
assert (u-B*v-g).applyfunc(s.cancel) == s.zeros(2)

sx = s.Matrix([[0, 1], [1, 0]])
H = s.diag(i, -i)
E = s.eye(2)-i*q0/z*sx
M = E.copy()
for j in range(2):
    M[:, 1] += C[j]*v[j, :].T/(z-zeta[j])
    M[:, 0] -= Cb[j]*u[j, :].T/(z-zbar[j])
M = M.applyfunc(s.factor)
# Only after the residue system is solved, factor its resulting matrix.
y = s.Matrix([T-h, 1-h*T])
yb = s.Matrix([[r-h, 1-h*r]])
P = y*yb/dm
k = (z-q0*q0/z)/2
lam = (z+q0*q0/z)/2
N = (k+i*b)*s.eye(2)-2*i*b*P
compact = H*N*E*s.diag(1/(lam+i*kappa), 1/(lam-i*kappa))*H**-1
assert (M-compact).applyfunc(lambda v: s.simplify(s.factor(v))) == s.zeros(2)
assert (P*P-P).applyfunc(s.factor) == s.zeros(2)
assert s.factor(N.det()-(lam*lam+kappa*kappa)) == 0
assert (M.subs(z, -q0*q0/z)-M*(i*z/q0)*sx).applyfunc(lambda v: s.simplify(s.factor(v))) == s.zeros(2)
assert (M.applyfunc(lambda v: s.limit(v, z, s.oo))-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)
assert ((z*M).applyfunc(lambda v: s.limit(v, z, 0))+i*q0*sx).applyfunc(s.simplify) == s.zeros(2)
q_from_infinity = i*s.limit(z*M[0, 1], z, s.oo)
assert s.factor(q_from_infinity+q0+4*b*P[0, 1]) == 0
# Independently recover the jump from the new scattering direction phi+=phi- S.
aa,bb,cc,phase = s.symbols('aa bb cc phase', nonzero=True)
dd=(1+bb*cc)/aa
Sm=s.Matrix([[aa,bb],[cc,dd]])
phi_minus=s.Matrix([[1,2],[3,7]])
phi_plus=phi_minus*Sm
def detcols(c1,c2): return s.det(s.Matrix.hstack(c1,c2))
wronskians=s.Matrix([[detcols(phi_plus[:,0],phi_minus[:,1]),detcols(phi_plus[:,1],phi_minus[:,1])],
                      [detcols(phi_minus[:,0],phi_plus[:,0]),detcols(phi_minus[:,0],phi_plus[:,1])]])
assert (wronskians-Sm).applyfunc(s.simplify)==s.zeros(2)
mu_minus=phi_minus
mu_plus=mu_minus*s.Matrix([[aa,bb*phase],[cc/phase,dd]])
Mplus=s.Matrix.hstack(mu_minus[:,0],mu_plus[:,1]/dd)
Mminus=s.Matrix.hstack(mu_plus[:,0]/aa,mu_minus[:,1])
jump=s.Matrix([[1,bb*phase/dd],[-cc/(aa*phase),1-bb*cc/(aa*dd)]])
assert (Mminus.inv()*Mplus-jump).applyfunc(s.factor)==s.zeros(2)
assert s.factor(jump.det()-1)==0
# One-orbit Jost limits: M tends Y on the negative side and Y diag(d,1/d)
# on the positive side, with d=s22; phi+=phi- diag(1/d,d).
Pminus=s.Matrix([[h*h,-h],[-h,1]])/(1+h*h)
Pplus=s.Matrix([[1,-h],[-h,h*h]])/(1+h*h)
def compact_at(P0):
    return H*((k+i*b)*s.eye(2)-2*i*b*P0)*E*s.diag(1/(lam+i*kappa),1/(lam-i*kappa))*H**-1
d=(lam-i*kappa)/(lam+i*kappa)
assert (compact_at(Pminus)-E).applyfunc(lambda v:s.simplify(s.factor(v)))==s.zeros(2)
assert (compact_at(Pplus)-E*s.diag(d,1/d)).applyfunc(lambda v:s.simplify(s.factor(v)))==s.zeros(2)

report = {
    'variables': '0<h<1, q0>0; r=conjugate(T), T!=0 on the physical slice',
    'Schur_determinant': str(s.factor(L.det())),
    'computed_v_rows': str(v), 'computed_u_rows': str(u),
    'checks': {name: 'exact zero' for name in (
        'Schur factorization', 'explicit v', 'explicit u', 'finite1', 'finite2',
        'partial fractions equal compact M', 'rank-one projector identity',
        'numerator determinant', 'sheet symmetry', 'infinity normalization',
        'zero normalization', 'optical reconstruction', 'new-direction Wronskians',
        'new-direction jump', 'jump determinant', 'negative Jost endpoint', 'positive Jost endpoint')},
}
out = Path(__file__).resolve().parents[1]/'calculations'/'direct-algebra-results.json'
out.write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
