"""Independent zero-curvature and background checks; SymPy, no numerical fitting."""
from pathlib import Path
import json
import sympy as s

x, z = s.symbols('x z', real=True)
l, a, al, be, Q, e0 = s.symbols('lambda a alpha beta Q eta0', real=True)
q, qb, p, pb, e = [s.Function(n)(x, z) for n in ('q','qb','p','pb','eta')]
I = s.I
sig = s.diag(1,-1)
A = s.Matrix([[0,q],[-qb,0]])
R = s.Matrix([[e,-p],[-pb,-e]])
U = -I*l*sig+A
V3 = 4*I*be*sig
V2 = -I*al*sig-4*be*A
V1 = s.Matrix([[-2*I*be*q*qb, al*q-2*I*be*s.diff(q,x)],
               [-al*qb-2*I*be*s.diff(qb,x),2*I*be*q*qb]])
d = I*al*q*qb/2-be*(q*s.diff(qb,x)-s.diff(q,x)*qb)
b = 2*be*q*qb*q+I*al*s.diff(q,x)/2+be*s.diff(q,x,2)
c = -2*be*q*qb*qb+I*al*s.diff(qb,x)/2-be*s.diff(qb,x,2)
V0 = s.Matrix([[d,b],[c,-d]])
Vp = l**3*V3+l**2*V2+l*V1+V0
V = Vp+I*R/(l-a)
comm = lambda B,C: B*C-C*B
flowq = I*al*(s.diff(q,x,2)/2+q*q*qb)+be*(s.diff(q,x,3)+6*q*qb*s.diff(q,x))+2*p
flowqb = -I*al*(s.diff(qb,x,2)/2+qb*qb*q)+be*(s.diff(qb,x,3)+6*q*qb*s.diff(qb,x))+2*pb
subs = {s.diff(q,z):flowq,s.diff(qb,z):flowqb,s.diff(p,x):2*e*q-2*I*a*p,
        s.diff(pb,x):2*e*qb+2*I*a*pb,s.diff(e,x):-(qb*p+q*pb)}
Z = (s.diff(U,z)-s.diff(V,x)+comm(U,V)).subs(subs).applyfunc(s.simplify)
RX = (s.diff(R,x)-comm(U.subs(l,a),R)).subs(subs).applyfunc(s.simplify)
bloch = s.simplify(s.diff(p*pb+e**2,x).subs(subs))
bg={q:Q,qb:Q,p:-I*e0*Q/a,pb:I*e0*Q/a,e:e0}
Om=-al*Q**2+2*e0/a
Ub=U.subs(bg)
Vb=(V.subs(bg).doit()+I*Om*sig/2).applyfunc(s.simplify)
F=4*be*l**2-al*l-2*be*Q**2+e0/(a*(l-a))
tc=(Vb+F*Ub).applyfunc(s.simplify)
w=s.symbols('w',nonzero=True)
lam=(w-Q**2/w)/2; rho=(w+Q**2/w)/2
E=s.Matrix([[1,-I*Q/w],[-I*Q/w,1]])
diagonal=(Ub.subs(l,lam)*E+I*rho*E*sig).applyfunc(s.simplify)
result={'zero_curvature':str(Z),'resonance_commutator':str(RX),
        'bloch_derivative':str(bloch),'background_time_identity':str(tc),
        'background_diagonalization':str(diagonal),
        'uniformization':str(s.simplify(rho**2-lam**2-Q**2))}
assert Z==s.zeros(2) and RX==s.zeros(2) and bloch==0 and tc==s.zeros(2) and diagonal==s.zeros(2)
out=Path(__file__).resolve().parents[1]/'evidence'/'symbolic-results.json'
out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
