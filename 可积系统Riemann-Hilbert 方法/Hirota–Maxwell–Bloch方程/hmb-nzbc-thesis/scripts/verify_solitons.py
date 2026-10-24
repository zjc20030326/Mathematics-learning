"""Direct PDE residuals by arbitrary-precision differentiation of explicit fields.
No use of the PDE to produce derivatives. Run from any working directory.
"""
from pathlib import Path
import json
import mpmath as m

m.mp.dps=60
Q=m.mpf(1); a=m.mpf(2); eta0=m.mpf('0.6')
alpha=m.mpf(1); beta=m.mpf('0.1')
Omega=-alpha*Q**2+2*eta0/a
def F(l): return 4*beta*l*l-alpha*l-2*beta*Q**2+eta0/(a*(l-a))
def seed_y(l,x,z):
    rho=m.sqrt(l*l+Q*Q)
    w=l+rho
    theta=rho*(x-F(l)*z)
    E=m.matrix([[1,-1j*Q/w],[-1j*Q/w,1]])
    return E*m.matrix([m.exp(-1j*theta),m.exp(1j*theta)])
def fields(x,z,n=1):
    # Construct decimal inputs at the active precision, without binary64 rounding.
    ls=[m.mpc(0,m.mpf('1.5')),m.mpc(0,m.mpf('1.9'))][:n]
    q=m.mpc(Q)
    R=m.matrix([[eta0,1j*eta0*Q/a],[-1j*eta0*Q/a,-eta0]])
    dress=[]
    for l in ls:
        l=m.mpc(l); y=seed_y(l,x,z)
        for lk,Pk in dress:
            y=(m.eye(2)-(lk-m.conj(lk))/(l-m.conj(lk))*Pk)*y
        P=y*y.H/(y.H*y)[0]
        delta=l-m.conj(l)
        q=q-2j*delta*P[0,1]
        D=m.eye(2)-delta/(a-m.conj(l))*P
        R=D*R*D.H
        dress.append((l,P))
    # Each pure imaginary simple dressing changes the common endpoint sign.
    phase=(-1)**n*m.exp(-1j*Omega*z)
    return phase*q,-phase*R[0,1],m.re(R[0,0])
def residuals(x,z,n):
    q,p,e=fields(x,z,n)
    qx=m.diff(lambda t:fields(t,z,n)[0],x)
    qxx=m.diff(lambda t:fields(t,z,n)[0],x,2)
    qxxx=m.diff(lambda t:fields(t,z,n)[0],x,3)
    qz=m.diff(lambda t:fields(x,t,n)[0],z)
    px=m.diff(lambda t:fields(t,z,n)[1],x)
    ex=m.diff(lambda t:fields(t,z,n)[2],x)
    return [abs(qz-1j*alpha*(qxx/2+abs(q)**2*q)-beta*(qxxx+6*abs(q)**2*qx)-2*p),
            abs(px-2*e*q+2j*a*p),abs(ex+m.conj(q)*p+q*m.conj(p)),
            abs(abs(p)**2+e**2-eta0**2*(1+Q**2/a**2))]
if __name__=='__main__':
    report={'precision_digits':m.mp.dps,'parameters':{'Q':'1','a':'2','eta0':'0.6','alpha':'1','beta':'0.1'},'cases':[]}
    for n in [1,2]:
        points=[('-1.2','-0.35'),('0','0'),('0.37','0.21'),('1.4','-0.1'),('2.3','0.6')]
        values=[residuals(m.mpf(x),m.mpf(z),n) for x,z in points]
        maxima=[max(r[k] for r in values) for k in range(4)]
        # Explicit asymptotic check, distinct from algebraic proof.
        errors=[]
        for xx in [-30,30]:
            z=m.mpf('0.21'); q,p,e=fields(m.mpf(xx),z,n)
            qb=Q*m.exp(-1j*Omega*z); pb=-1j*eta0*qb/a
            errors.append([m.nstr(abs(q-qb),10),m.nstr(abs(p-pb),10),m.nstr(abs(e-eta0),10)])
        report['cases'].append({'N':n,'spectral_points':['1.5i','1.9i'][:n],
            'points':points,'max_residual_q_p_eta_bloch':[m.nstr(v,12) for v in maxima],
            'boundary_errors_x_minus_plus_30':errors})
        assert max(maxima)<m.mpf('1e-45')
    out=Path(__file__).resolve().parents[1]/'calculations'/'soliton-results.json'
    out.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
