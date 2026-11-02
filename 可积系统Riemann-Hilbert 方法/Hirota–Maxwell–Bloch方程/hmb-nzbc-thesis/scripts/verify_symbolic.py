"""Independent zero-curvature and background checks; SymPy, no numerical fitting."""
from pathlib import Path
import json
import sympy as s

x, t = s.symbols('x t', real=True)
k, a, al, be, q0, e0 = s.symbols('k a alpha beta q0 eta0', real=True)
q, qb, p, pb, e = [s.Function(n)(x, t) for n in ('q','qb','p','pb','eta')]
I = s.I
sig = s.diag(1,-1)
A = s.Matrix([[0,q],[-qb,0]])
R = s.Matrix([[e,-p],[-pb,-e]])
U = -I*k*sig+A
V3 = 4*I*be*sig
V2 = -I*al*sig-4*be*A
V1 = s.Matrix([[-2*I*be*q*qb, al*q-2*I*be*s.diff(q,x)],
               [-al*qb-2*I*be*s.diff(qb,x),2*I*be*q*qb]])
d = I*al*q*qb/2-be*(q*s.diff(qb,x)-s.diff(q,x)*qb)
b = 2*be*q*qb*q+I*al*s.diff(q,x)/2+be*s.diff(q,x,2)
c = -2*be*q*qb*qb+I*al*s.diff(qb,x)/2-be*s.diff(qb,x,2)
V0 = s.Matrix([[d,b],[c,-d]])
Vp = k**3*V3+k**2*V2+k*V1+V0
V = Vp+I*R/(k-a)
comm = lambda B,C: B*C-C*B
flowq = I*al*(s.diff(q,x,2)/2+q*q*qb)+be*(s.diff(q,x,3)+6*q*qb*s.diff(q,x))+2*p
flowqb = -I*al*(s.diff(qb,x,2)/2+qb*qb*q)+be*(s.diff(qb,x,3)+6*q*qb*s.diff(qb,x))+2*pb
subs = {s.diff(q,t):flowq,s.diff(qb,t):flowqb,s.diff(p,x):2*e*q-2*I*a*p,
        s.diff(pb,x):2*e*qb+2*I*a*pb,s.diff(e,x):-(qb*p+q*pb)}
Z = (s.diff(U,t)-s.diff(V,x)+comm(U,V)).subs(subs).applyfunc(s.simplify)
RX = (s.diff(R,x)-comm(U.subs(k,a),R)).subs(subs).applyfunc(s.simplify)
bloch = s.simplify(s.diff(p*pb+e**2,x).subs(subs))
bg={q:q0,qb:q0,p:-I*e0*q0/a,pb:I*e0*q0/a,e:e0}
Om=-al*q0**2+2*e0/a
Ub=U.subs(bg)
Vb=(V.subs(bg).doit()+I*Om*sig/2).applyfunc(s.simplify)
F=4*be*k**2-al*k-2*be*q0**2+e0/(a*(k-a))
tc=(Vb+F*Ub).applyfunc(s.simplify)
z=s.symbols('z',nonzero=True)
kz=(z-q0**2/z)/2; lam=(z+q0**2/z)/2
E=s.Matrix([[1,-I*q0/z],[-I*q0/z,1]])
diagonal=(Ub.subs(k,kz)*E+I*lam*E*sig).applyfunc(s.simplify)
result={'zero_curvature':str(Z),'resonance_commutator':str(RX),
        'bloch_derivative':str(bloch),'background_time_identity':str(tc),
        'background_diagonalization':str(diagonal),
        'uniformization':str(s.simplify(lam**2-kz**2-q0**2))}
assert Z==s.zeros(2) and RX==s.zeros(2) and bloch==0 and tc==s.zeros(2) and diagonal==s.zeros(2)
out=Path(__file__).resolve().parents[1]/'calculations'/'symbolic-results.json'
out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
