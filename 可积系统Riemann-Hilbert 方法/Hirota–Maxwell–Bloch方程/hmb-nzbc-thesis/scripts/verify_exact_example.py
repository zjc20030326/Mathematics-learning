"""Exact rational PDE verification after exponential-coordinate substitution.
Q=1, b=5/4, kappa=3/4, a=2, eta0=3/5, alpha=1, beta=1/10.
T=exp(2 Re(s))>0, C=exp(2 i Im(s)); conjugation sends C to 1/C.
"""
from pathlib import Path
import json,time
import sympy as s
start=time.time()
T=s.Symbol('T',positive=True); C=s.Symbol('C',nonzero=True)
I=s.I; rat=s.Rational
Q=1; b=rat(5,4); k=rat(3,4); a=2; eta0=rat(3,5)
alpha=1; beta=rat(1,10); Omega=-alpha+2*eta0/a
lam=I*b
F=4*beta*lam**2-alpha*lam-2*beta+eta0/(a*(lam-a))
F=s.simplify(s.expand_complex(F))
reF=s.re(F); imF=s.im(F)
yb=s.Matrix([[T/C-rat(1,2),1-T/(2*C)]])
y=s.Matrix([T*C-rat(1,2),1-T*C/2])
P=(y*yb/(yb*y)[0]).applyfunc(s.cancel)
D=s.eye(2)-2*I*b/(a+I*b)*P
Di=s.eye(2)+2*I*b/(a-I*b)*P
R0=s.Matrix([[eta0,I*eta0/a],[-I*eta0/a,-eta0]])
R=(D*R0*Di).applyfunc(s.cancel)
u=s.cancel(-(1+4*b*P[0,1])); ub=s.cancel(-(1+4*b*P[1,0]))
v=R[0,1]; vb=R[1,0]; e=R[0,0]
def dx(f):return s.cancel(2*k*T*s.diff(f,T))
def dz(f):return s.cancel(-2*k*reF*T*s.diff(f,T)-2*I*k*imF*C*s.diff(f,C))
ux=dx(u); uxx=dx(ux); uxxx=dx(uxx)
r1=s.cancel(dz(u)-I*Omega*u-I*alpha*(uxx/2+u*ub*u)-beta*(uxxx+6*u*ub*ux)-2*v)
r2=s.cancel(dx(v)-2*e*u+2*I*a*v)
r3=s.cancel(dx(e)+ub*v+u*vb)
rb=s.cancel(v*vb+e**2-eta0**2*(1+rat(1,4)))
assert [r1,r2,r3,rb]==[0,0,0,0]
report={'parameters':{'Q':'1','b':'5/4','kappa':'3/4','a':'2','eta0':'3/5','alpha':'1','beta':'1/10'},
        'F':str(F),'variables':'T=exp(2 Re(s)), C=exp(2 i Im(s)), s=kappa*(x-F*z)',
        'exact_residuals_q_p_eta_bloch':list(map(str,[r1,r2,r3,rb])),
        'elapsed_seconds':round(time.time()-start,3)}
out=Path(__file__).resolve().parents[1]/'calculations'/'exact-example-results.json'
out.write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
