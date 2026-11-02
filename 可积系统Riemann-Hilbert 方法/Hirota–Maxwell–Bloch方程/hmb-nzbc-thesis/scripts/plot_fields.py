"""Evaluate the one-orbit reflectionless RHP on plotting grids.

Each grid point solves the same 4 by 4 residue block system as
verify_direct_rhp.py. Fields are reconstructed from infinity and resonance.
No projector, Darboux recursion, or PDE discretization is used.
Run from any working directory. Requires NumPy, matplotlib and mpmath.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = dict(Q=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5)


def dispersion(lam, Q, a, eta0, alpha, beta):
    return 4*beta*lam*lam-alpha*lam-2*beta*Q*Q+eta0/(a*(lam-a))


def calc(X, Z, *, Q=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5):
    """Return q,p,eta and diagnostics from batched residue-system solves."""
    X, Z = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(Z, dtype=float))
    shape = X.shape
    kappa = np.sqrt(b*b-Q*Q)
    zeta = np.array([1j*(b+kappa), -1j*(b-kappa)])
    c1 = 4j*kappa/(1-Q*Q/zeta[0]**2)
    constants = np.array([c1, -Q*Q/np.conj(zeta[0])**2*np.conj(c1)])
    lam = (zeta-Q*Q/zeta)/2
    rho = (zeta+Q*Q/zeta)/2
    F = dispersion(lam, Q, a, eta0, alpha, beta)
    C = constants*np.exp(2j*rho*(X[..., None]-F*Z[..., None]))
    A = np.conj(C)[..., None, :]/(zeta[:, None]-np.conj(zeta)[None, :])
    B = C[..., None, :]/(np.conj(zeta)[:, None]-zeta[None, :])
    L = np.broadcast_to(np.eye(4, dtype=complex), shape+(4, 4)).copy()
    L[..., :2, 2:] = A
    L[..., 2:, :2] = -B
    f = np.column_stack((-1j*Q/zeta, np.ones(2)))
    g = np.column_stack((np.ones(2), -1j*Q/np.conj(zeta)))
    rhs = np.broadcast_to(np.vstack((f, g)), shape+(4, 2))
    # Row equilibration preserves the system and limits the grid's dynamic range.
    scale = np.max(np.abs(L), axis=-1)
    answer = np.linalg.solve(L/scale[..., :, None], rhs/scale[..., :, None])
    v, u = answer[..., :2, :], answer[..., 2:, :]
    field = Q-1j*np.sum(np.conj(C)*u[..., :, 0], axis=-1)
    wa = a+np.sqrt(a*a+Q*Q)
    rhoa = (wa+Q*Q/wa)/2
    M = np.broadcast_to(np.array([[1, -1j*Q/wa], [-1j*Q/wa, 1]]), shape+(2, 2)).copy()
    M[..., :, 0] += np.sum(C[..., :, None]*v/(wa-zeta)[..., None], axis=-2)
    M[..., :, 1] -= np.sum(np.conj(C)[..., :, None]*u/(wa-np.conj(zeta))[..., None], axis=-2)
    R = M@np.diag([eta0*rhoa/a, -eta0*rhoa/a])@np.linalg.inv(M)
    phase = np.exp(-1j*(-alpha*Q*Q+2*eta0/a)*Z)
    q, p, eta = phase*field, -phase*R[..., 0, 1], R[..., 0, 0]
    linear = np.max(np.abs(L@answer-rhs)/(1+np.abs(L)@np.abs(answer)+np.abs(rhs)))
    bloch = np.max(np.abs(np.abs(p)**2+eta**2-eta0**2*(1+Q*Q/a**2)))
    return q, p, eta, dict(scaled_linear_residual=float(linear),
                          Bloch_error=float(bloch), eta_imaginary=float(np.max(abs(eta.imag))))


def velocity_period(alpha=1., beta=.1):
    Q, a, eta0, b = (DEFAULT[key] for key in ('Q', 'a', 'eta0', 'b'))
    k = np.sqrt(b*b-Q*Q)
    F = dispersion(1j*b, Q, a, eta0, alpha, beta)
    return float(F.real), float(F.imag), None if abs(F.imag)<1e-14 else float(np.pi/(k*abs(F.imag)))


def main():
    out = ROOT/'figures'
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    x, z = np.linspace(-6, 6, 401), np.linspace(-1.5, 1.5, 251)
    X, Z = np.meshgrid(x, z)
    q, p, eta, base = calc(X, Z)
    assert max(base.values()) < 1e-10
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.35), layout='constrained')
    for ax, data, title, cmap in zip(axes, [abs(q), abs(p), eta.real],
            [r'$|q(x,z)|$', r'$|p(x,z)|$', r'$\eta(x,z)$'], ['viridis', 'magma', 'coolwarm']):
        im = ax.pcolormesh(x, z, data, shading='auto', cmap=cmap, rasterized=True)
        ax.set(xlabel=r'$x$', ylabel=r'$z$', title=title)
        fig.colorbar(im, ax=ax, shrink=.88)
    fig.savefig(out/'one-soliton-fields.pdf')
    plt.close(fig)

    speed, _, period = velocity_period()
    profile_x = np.linspace(-4, 4, 401)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), layout='constrained')
    for frac, ls in zip([0., .25, .5], ['-', '--', ':']):
        zz = frac*period
        qq, _, _, _ = calc(profile_x+speed*zz, zz)
        axes[0].plot(profile_x, abs(qq), ls, label=rf'$z/T_z={frac:g}$')
    axes[0].set(xlabel=r'$X=x-A_bz$', ylabel=r'$|q|$', title='Comoving profiles')
    axes[0].legend()
    zs = np.linspace(-2, 2, 501)
    stationary_alpha = -DEFAULT['eta0']/(DEFAULT['a']*(DEFAULT['a']**2+DEFAULT['b']**2))
    periodic_cases = []
    for al in [stationary_alpha, .5, 1.]:
        aa, bb, tt = velocity_period(alpha=al)
        qq, _, _, _ = calc(aa*zs, zs, alpha=al)
        axes[1].plot(zs, abs(qq), label=rf'$\alpha={al:g}$')
        periodic_cases.append(dict(alpha=al, A_b=aa, B_b=bb, T_z=tt))
    axes[1].set(xlabel=r'$z$', ylabel=r'$|q(A_bz,z)|$', title='Internal oscillation')
    axes[1].legend()
    fig.savefig(out/'one-soliton-dynamics.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), layout='constrained')
    velocity_cases = []
    for ax, be in zip(axes, [0., .15]):
        qq, _, _, errs = calc(X, Z, beta=be)
        aa, bb, tt = velocity_period(beta=be)
        im = ax.pcolormesh(x, z, abs(qq), shading='auto', cmap='viridis',
                           rasterized=True, vmin=0, vmax=4)
        ax.plot(aa*z, z, 'w--', lw=1, label=r'$x=A_bz$')
        ax.set(xlabel=r'$x$', ylabel=r'$z$', title=rf'$\beta={be:g},\ A_b={aa:g}$')
        ax.legend(loc='upper right')
        fig.colorbar(im, ax=ax, shrink=.88)
        assert max(errs.values()) < 1e-10
        velocity_cases.append(dict(beta=be, A_b=aa, B_b=bb, T_z=tt, errors=errs))
    fig.savefig(out/'one-soliton-velocity.pdf')
    plt.close(fig)

    # Independent high-precision residue implementation; no Darboux import.
    import mpmath as mp
    from verify_direct_rhp import fields
    mp.mp.dps = 70
    points = [('-1.2', '-.35'), ('0', '0'), ('.37', '.21'), ('1.4', '-.1'), ('2.3', '.6')]
    maxima = np.zeros(3)
    for xx, zz in points:
        values = calc(float(xx), float(zz))[:3]
        high = fields(mp.mpf(xx), mp.mpf(zz))
        maxima = np.maximum(maxima, [abs(complex(v)-complex(hv)) for v, hv in zip(values, high)])
    assert max(maxima) < 1e-12
    report = dict(method='NumPy solution of equilibrated 4x4 residue block system; infinity/resonance reconstruction',
                  parameters=DEFAULT, grid=dict(x_min=-6, x_max=6, x_count=401,
                  z_min=-1.5, z_max=1.5, z_count=251), base_errors=base,
                  high_precision_comparison_points=points,
                  high_precision_comparison_q_p_eta=list(maxima),
                  dynamics=dict(profile_z_over_T=[0, .25, .5],
                                alpha_cases=periodic_cases, beta_cases=velocity_cases),
                  files=['one-soliton-fields.pdf', 'one-soliton-dynamics.pdf', 'one-soliton-velocity.pdf'])
    (ROOT/'calculations'/'plot-results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
