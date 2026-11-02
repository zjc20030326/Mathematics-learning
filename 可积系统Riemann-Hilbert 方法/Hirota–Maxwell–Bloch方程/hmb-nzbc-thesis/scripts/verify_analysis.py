"""Exact algebra for the high-energy transformation and resonance regularity."""
from pathlib import Path
import json
import sympy as s

I=s.I
q,r,qx,rx,Q,w=s.symbols("q r qx rx Q w",nonzero=True)
sigma=s.diag(1,-1)
B=s.Matrix([[0,q],[r,0]])
Bx=s.Matrix([[0,qx],[rx,0]])
lam=(w-Q**2/w)/2
rho=(w+Q**2/w)/2
U=-I*lam*sigma+sigma*B
T=s.eye(2)-I*B/w
C1=I*((Q**2-q*r)*sigma+Bx)
gauge=(U*T+I*rho*T*sigma+I*Bx/w-C1/w).applyfunc(s.simplify)
assert gauge==s.zeros(2)

a,b,e,alpha,beta=s.symbols("a b eta0 alpha beta",real=True,nonzero=True)
F=4*beta*(I*b)**2-alpha*I*b-2*beta*Q**2+e/(a*(I*b-a))
Qreal=s.symbols("Qreal",positive=True)
F=F.subs(Q,Qreal)
A=-4*beta*b**2-2*beta*Qreal**2-e/(a*a+b*b)
Bphase=-alpha*b-e*b/(a*(a*a+b*b))
phase=s.simplify(F-A-I*Bphase)
assert phase==0

v,c=s.symbols("v c",real=True,nonzero=True)
oscillatory=[]
for m in [0,1,3,5,9]:
    f=v**m*(1+v+v*v)*s.exp(I*c/v)
    G=1+v+v*v
    for j in range(5):
        residual=s.simplify(s.diff(f,v,j)-s.exp(I*c/v)*v**(m-2*j)*G)
        assert residual==0
        assert s.simplify(G.subs(v,0)-(-I*c)**j)==0
        G=s.expand(-I*c*G+(m-2*j)*v*G+v*v*s.diff(G,v))
    oscillatory.append({"zero_order":m,"derivatives_checked":4,"residual":"0"})

T0,C,h=s.symbols("T C h",nonzero=True)
y=s.Matrix([T0*C-h,1-h*T0*C])
yb=s.Matrix([T0/C-h,1-h*T0/C])
den=s.expand((yb.T*y)[0])
dencheck=s.simplify(den-((1+h*h)*(T0*T0+1)-2*h*T0*(C+1/C)))
assert dencheck==0
dar=s.simplify(2*I*s.Rational(5,4)/(2+I*s.Rational(5,4))-(50+80*I)/89)
assert dar==0

report={"high_energy_gauge_residual":str(gauge),"phase_decomposition_residual":str(phase),
        "rational_denominator_residual":str(dencheck),"exact_Darboux_factor_residual":str(dar),
        "oscillatory_derivative_checks":oscillatory}
out=Path(__file__).resolve().parents[1]/"calculations"/"analysis-results.json"
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
