"""Verify the one-quartet RHP independently of the PDE residual calculation."""
from pathlib import Path
import json
import mpmath as m
from verify_solitons import Q,a,eta0,alpha,beta,Omega,F,seed_y,fields
m.mp.dps=60
b=m.mpf('1.5'); k=m.sqrt(b*b-Q*Q)
H=m.matrix([[1j,0],[0,-1j]])
sig=m.matrix([[1,0],[0,-1]])
sx=m.matrix([[0,1],[1,0]])
K=m.matrix([[0,1],[-1,0]])
def geom(w):return (w-Q*Q/w)/2,(w+Q*Q/w)/2
def proj(x,z):
    y=seed_y(1j*b,x,z)
    return y*y.H/(y.H*y)[0]
def Mrhp(w,x,z):
    l,r=geom(w); E=m.eye(2)-1j*Q/w*sx
    return H*((l+1j*b)*m.eye(2)-2j*b*proj(x,z))*E*m.diag([1/(r-1j*k),1/(r+1j*k)])*H**-1
def maxnorm(A):return max(abs(v) for v in A)
def conjmat(A):return A.apply(m.conj)
x=m.mpf('.37'); z=m.mpf('.21')
zeta=1j*(b+k); hat=-Q*Q/m.conj(zeta)
c1=2j*k/((1-Q*Q/zeta**2)/2)
ch=-Q*Q/m.conj(zeta)**2*m.conj(c1)
checks={}
for name,w0,c in [('outer',zeta,c1),('inner',hat,ch)]:
    l,r=geom(w0); E=m.eye(2)-1j*Q/w0*sx
    N=H*((l+1j*b)*m.eye(2)-2j*b*proj(x,z))*E*H**-1
    res=N[:,0]/((1-Q*Q/w0**2)/2)
    col2=N[:,1]/(r+1j*k)
    C=c*m.exp(2j*r*(x-F(l)*z))
    checks['residue_'+name]=maxnorm(res-C*col2)
for idx,w in enumerate([m.mpc('2.1','0.7'),m.mpc('.3','-.2'),m.mpc('1.2','0')]):
    M=Mrhp(w,x,z); l,r=geom(w)
    checks[f'determinant_{idx}']=abs(m.det(M)-(1+Q*Q/w**2))
    checks[f'sheet_{idx}']=maxnorm(Mrhp(-Q*Q/w,x,z)-M*(1j*w/Q)*sx)
    checks[f'conjugation_{idx}']=maxnorm(M-K*conjmat(Mrhp(m.conj(w),x,z))*K**-1)
    q,p,e=fields(x,z,1); phase=m.exp(1j*Omega*z)
    u=phase*q
    A=m.matrix([[0,u],[-m.conj(u),0]])
    Mx=m.matrix(2,2)
    for i in range(2):
        for j in range(2):Mx[i,j]=m.diff(lambda t:Mrhp(w,t,z)[i,j],x)
    checks[f'spatial_lax_{idx}']=maxnorm(Mx-(-1j*l*sig+A)*M-1j*r*M*sig)
for idx,wa in enumerate([a+m.sqrt(a*a+Q*Q),a-m.sqrt(a*a+Q*Q)]):
    M=Mrhp(wa,x,z); _,r=geom(wa)
    R=M*(eta0*r/a*sig)*M**-1
    q,p,e=fields(x,z,1); v=m.exp(1j*Omega*z)*p
    checks[f'medium_reconstruction_{idx}']=maxnorm(R-m.matrix([[e,-v],[-m.conj(v),-e]]))
assert max(checks.values())<m.mpf('1e-45')
result={'precision_digits':60,'point_x_z':['0.37','0.21'],
        'checks':{name:m.nstr(value,12) for name,value in checks.items()}}
(Path(__file__).resolve().parents[1]/'calculations'/'rhp-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
