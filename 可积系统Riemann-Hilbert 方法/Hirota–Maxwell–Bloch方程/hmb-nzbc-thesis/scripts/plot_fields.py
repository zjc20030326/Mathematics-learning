"""Vectorized plotting of the same explicit dressing; no PDE solver."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
out=root/'figures'; out.mkdir(exist_ok=True)
Q=1.; a=2.; eta0=.6; alpha=1.; beta=.1
Omega=-alpha*Q**2+2*eta0/a
def F(l):return 4*beta*l*l-alpha*l-2*beta*Q**2+eta0/(a*(l-a))
def calc(X,Z,n):
    shape=X.shape
    eye=np.broadcast_to(np.eye(2,dtype=complex),shape+(2,2))
    R=np.broadcast_to(np.array([[eta0,1j*eta0*Q/a],[-1j*eta0*Q/a,-eta0]]),shape+(2,2)).copy()
    q=np.full(shape,Q,dtype=complex); history=[]
    for l in [1.5j,1.9j][:n]:
        rho=np.sqrt(l*l+Q*Q); w=l+rho
        theta=rho*(X-F(l)*Z)
        v=np.stack([np.exp(-1j*theta),np.exp(1j*theta)],axis=-1)
        E=np.array([[1,-1j*Q/w],[-1j*Q/w,1]])
        y=np.einsum('ij,...j->...i',E,v)
        for lk,Pk in history:
            y=np.einsum('...ij,...j->...i',eye-(lk-lk.conjugate())/(l-lk.conjugate())*Pk,y)
        P=y[..., :,None]*y.conj()[...,None,:]/np.sum(abs(y)**2,axis=-1)[...,None,None]
        delta=l-l.conjugate(); q-=2j*delta*P[...,0,1]
        D=eye-delta/(a-l.conjugate())*P
        R=D@R@np.swapaxes(D.conj(),-1,-2)
        history.append((l,P))
    phase=(-1)**n*np.exp(-1j*Omega*Z)
    return phase*q,-phase*R[...,0,1],R[...,0,0].real
x=np.linspace(-6,6,401); z=np.linspace(-1.5,1.5,251); X,Z=np.meshgrid(x,z)
for n,name in [(1,'one'),(2,'two')]:
    q,p,e=calc(X,Z,n)
    fig,axes=plt.subplots(1,3,figsize=(12,3.35),layout='constrained')
    for ax,data,title,cmap in zip(axes,[abs(q),abs(p),e],[r'$|q(x,z)|$',r'$|p(x,z)|$',r'$\eta(x,z)$'],['viridis','magma','coolwarm']):
        im=ax.pcolormesh(x,z,data,shading='auto',cmap=cmap,rasterized=True)
        ax.set(xlabel=r'$x$',ylabel=r'$z$',title=title)
        fig.colorbar(im,ax=ax,shrink=.88)
    fig.savefig(out/f'{name}-soliton-fields.pdf')
    plt.close(fig)
    err=np.max(abs(abs(p)**2+e**2-eta0**2*(1+Q**2/a**2)))
    print(f'{name}: grid {X.shape}, Bloch error {err:.3e}')
