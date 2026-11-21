"""Evaluate reflectionless RHP soliton solutions and generate thesis figures.

Generates spectral diagrams (analytic-regions.pdf, discrete-spectrum.pdf) and
three-panel soliton field figures (one-soliton-fields.pdf, two-soliton-fields.pdf)
by solving the finite-dimensional residue block system directly at each grid point.
Run from any working directory. Requires NumPy, matplotlib and mpmath.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge
from matplotlib.transforms import ScaledTranslation

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = dict(q0=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5)
NUMERIC_WIDTH_IN = 146.4/25.4

QA = Path(os.environ.get('HMB_FIGURE_QA', Path(tempfile.gettempdir())/'hmb-zhang-fan-figure-qa'))
AUDIT = False
SKILL_SCRIPTS = Path(os.environ.get(
    'NATURE_FIGURE_SCRIPTS',
    Path.home()/'.gemini'/'config'/'skills'/'nature-figure'/'scripts'))


def dispersion(k, q0, a, eta0, alpha, beta):
    return 4*beta*k*k-alpha*k-2*beta*q0*q0+eta0/(a*(k-a))


def eval_reflectionless(X, TIME, zeta, constants, *, q0=1., a=2., eta0=.6, alpha=1., beta=.1):
    """Evaluate (q, p, eta) and diagnostics via equilibrated 2m x 2m residue-system solves."""
    X, TIME = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(TIME, dtype=float))
    shape = X.shape
    zeta = np.asarray(zeta, dtype=complex)
    constants = np.asarray(constants, dtype=complex)
    m = len(zeta)
    k = (zeta-q0*q0/zeta)/2
    lam = (zeta+q0*q0/zeta)/2
    F = dispersion(k, q0, a, eta0, alpha, beta)
    C = constants*np.exp(-2j*lam*(X[..., None]-F*TIME[..., None]))
    A = np.conj(C)[..., None, :]/(zeta[:, None]-np.conj(zeta)[None, :])
    B = C[..., None, :]/(np.conj(zeta)[:, None]-zeta[None, :])
    L = np.broadcast_to(np.eye(2*m, dtype=complex), shape+(2*m, 2*m)).copy()
    L[..., :m, m:] = A
    L[..., m:, :m] = -B
    f = np.column_stack((np.ones(m), -1j*q0/zeta))
    g = np.column_stack((-1j*q0/np.conj(zeta), np.ones(m)))
    rhs = np.broadcast_to(np.vstack((f, g)), shape+(2*m, 2))
    scale = np.max(np.abs(L), axis=-1)
    answer = np.linalg.solve(L/scale[..., :, None], rhs/scale[..., :, None])
    v, u = answer[..., :m, :], answer[..., m:, :]
    field = q0+1j*np.sum(C*v[..., :, 0], axis=-1)
    za = a+np.sqrt(a*a+q0*q0)
    lambda_a = (za+q0*q0/za)/2
    M = np.broadcast_to(np.array([[1, -1j*q0/za], [-1j*q0/za, 1]], dtype=complex), shape+(2, 2)).copy()
    M[..., :, 1] += np.sum(C[..., :, None]*v/(za-zeta)[..., None], axis=-2)
    M[..., :, 0] -= np.sum(np.conj(C)[..., :, None]*u/(za-np.conj(zeta))[..., None], axis=-2)
    R = M@np.diag([eta0*lambda_a/a, -eta0*lambda_a/a])@np.linalg.inv(M)
    phase = np.exp(-1j*(-alpha*q0*q0+2*eta0/a)*TIME)
    q, p, eta = phase*field, -phase*R[..., 0, 1], R[..., 0, 0]
    linear = np.max(np.abs(L@answer-rhs)/(1+np.abs(L)@np.abs(answer)+np.abs(rhs)))
    bloch = np.max(np.abs(np.abs(p)**2+eta**2-eta0**2*(1+q0*q0/a**2)))
    return q, p, eta, dict(scaled_linear_residual=float(linear),
                          Bloch_error=float(bloch), eta_imaginary=float(np.max(abs(eta.imag))))


def one_soliton_data(q0=1., b=1.5):
    kappa = np.sqrt(b*b-q0*q0)
    zeta = np.array([1j*(b+kappa), -1j*(b-kappa)], dtype=complex)
    c1 = 4j*kappa/(1-q0*q0/zeta[0]**2)
    return zeta, np.array([c1, -q0*q0/np.conj(zeta[0])**2*np.conj(c1)], dtype=complex)


def two_soliton_data(q0=1.):
    z1 = 2*q0*np.exp(1j*np.pi/3)
    z2 = 3*q0*np.exp(2j*np.pi/3)
    zeta = np.array([z1, z2, -q0*q0/np.conj(z1), -q0*q0/np.conj(z2)], dtype=complex)
    def sp(idx):
        z = zeta[idx]
        return np.prod([z-zj for j, zj in enumerate(zeta) if j != idx])/np.prod([z-np.conj(zj) for zj in zeta])
    c1, c2 = 1/sp(0), 1/sp(1)
    c3 = -q0*q0/np.conj(z1)**2*np.conj(c1)
    c4 = -q0*q0/np.conj(z2)**2*np.conj(c2)
    return zeta, np.array([c1, c2, c3, c4], dtype=complex)


def calc(X, TIME, *, q0=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5):
    zeta, constants = one_soliton_data(q0, b)
    return eval_reflectionless(X, TIME, zeta, constants, q0=q0, a=a, eta0=eta0, alpha=alpha, beta=beta)


def calc_two_soliton(X, TIME, *, q0=1., a=2., eta0=.6, alpha=1., beta=.1):
    zeta, constants = two_soliton_data(q0)
    return eval_reflectionless(X, TIME, zeta, constants, q0=q0, a=a, eta0=eta0, alpha=alpha, beta=beta)


def panel_label(ax, letter):
    transform = ax.transAxes+ScaledTranslation(-8/72, 4/72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, letter, transform=transform, fontsize=11, fontweight='bold', va='bottom')


def export(fig, name, *, axes=None, exclude=None):
    QA.mkdir(parents=True, exist_ok=True)
    kwargs = {'exclude_axes': exclude or []}
    if axes is not None:
        kwargs['axes'] = axes
    fig.canvas.draw()
    if AUDIT:
        require_matplotlib_panel_alignment(
            fig, json_out=QA/(name+'.alignment.json'), tolerance_pt=1.5,
            gutter_tolerance_pt=1.5, strict=True, **kwargs)
    fig.savefig(ROOT/'figures'/(name+'.pdf'), dpi=600)
    if AUDIT:
        fig.savefig(QA/(name+'.svg'))
    fig.savefig(QA/(name+'.png'), dpi=300)
    plt.close(fig)


def spectral_base():
    label_size = 14
    fig, ax = plt.subplots(figsize=(6.18, 4.2))
    fig.subplots_adjust(left=.06, right=.98, bottom=.07, top=.97)
    ax.set_aspect('equal', adjustable='box')
    ax.set(xlim=(-2.75, 2.75), ylim=(-2.10, 2.10))
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    blue = '#E1EDF5'
    ax.add_patch(Rectangle((-2.50, 0), 5., 1.82, facecolor=blue, edgecolor='none', zorder=-3))
    ax.add_patch(Wedge((0, 0), 1, 0, 180, facecolor='white', edgecolor='none', zorder=-2))
    ax.add_patch(Wedge((0, 0), 1, 180, 360, facecolor=blue, edgecolor='none', zorder=-2))
    ax.plot([0, 0], [-1.85, 1.87], color='#8A8A8A', lw=.6, ls=(0, (3, 4)), zorder=-1)
    ax.annotate('', (0, 1.96), (0, 1.76), arrowprops=dict(arrowstyle='->', lw=.7, color='#696969'))
    ax.plot([-2.55, 2.56], [0, 0], color='#242424', lw=.85)
    ang = np.linspace(0, 2*np.pi, 801)
    ax.plot(np.cos(ang), np.sin(ang), color='#242424', lw=.85)
    for x0, x1 in [(-2.25, -1.85), (1.85, 2.25), (-.18, -.57), (.57, .18)]:
        ax.annotate('', (x1, 0), (x0, 0), arrowprops=dict(arrowstyle='-|>', lw=.85, color='#242424'))
    for th0, th1 in [(133, 125), (227, 235)]:
        a0, a1 = np.deg2rad([th0, th1])
        ax.annotate('', (np.cos(a1), np.sin(a1)), (np.cos(a0), np.sin(a0)),
                    arrowprops=dict(arrowstyle='-|>', lw=.85, color='#242424', connectionstyle='arc3,rad=.02'))
    for x, y, text in [(-1.76, 1.53, r'$D_+$'), (-1.76, -1.53, r'$D_-$'),
                        (.43, -.43, r'$D_+$'), (.43, .43, r'$D_-$')]:
        ax.text(x, y, text, ha='center', va='center', fontsize=label_size)
    ax.text(2.40, .14, r'$\operatorname{Re}z$', ha='center', fontsize=label_size)
    ax.text(.14, 1.94, r'$\operatorname{Im}z$', va='center', fontsize=label_size)
    ax.text(.10, -.18, r'$0$', ha='left', va='top', fontsize=label_size)
    ax.plot([-1, 1], [0, 0], ls='', marker='x', ms=5.5, color='#242424')
    ax.text(-1.20, -.13, r'$-q_0$', ha='center', va='top', fontsize=label_size)
    ax.text(1.20, -.13, r'$q_0$', ha='center', va='top', fontsize=label_size)
    ax.plot([0, 0], [-1, 1], ls='', marker='s', ms=4.7, mfc='white', mec='#6D4C64', mew=.85)
    ax.text(.12, 1.04, r'$\mathrm{i}q_0$', va='bottom', fontsize=label_size)
    ax.text(.12, -1.05, r'$-\mathrm{i}q_0$', va='top', fontsize=label_size)
    return fig, ax


def spectral_figures():
    label_size = 14
    fig, ax = spectral_base()
    export(fig, 'analytic-regions')
    fig, ax = spectral_base()
    zeta = 1.54+1.30j
    hzeta = -1/zeta
    plus = [zeta, np.conj(hzeta)]
    minus = [np.conj(zeta), hzeta]
    ax.plot(np.real(plus), np.imag(plus), ls='', marker='o', ms=5.5, mfc='#2B638A', mec='#2B638A', mew=.9)
    ax.plot(np.real(minus), np.imag(minus), ls='', marker='o', ms=5.5, mfc='white', mec='#242424', mew=.95)
    ax.text(zeta.real+.14, zeta.imag+.04, r'$\zeta_j$', fontsize=label_size, ha='left', va='bottom')
    ax.text(zeta.real+.14, -zeta.imag-.04, r'$(\zeta_j)^*$', fontsize=label_size, ha='left', va='top')
    ax.annotate(r'$h(\zeta_j)$', (hzeta.real, hzeta.imag), (-1.72, .50), fontsize=label_size, ha='center', va='bottom',
                arrowprops=dict(arrowstyle='-', lw=.55, color='#5E5E5E'))
    ax.annotate(r'$h(\zeta_j^*)$', (hzeta.real, -hzeta.imag), (-1.72, -.50), fontsize=label_size, ha='center', va='top',
                arrowprops=dict(arrowstyle='-', lw=.55, color='#5E5E5E'))
    export(fig, 'discrete-spectrum')
    return dict(type='symbolic general complex quartet',
                poles_D_plus=['zeta_j', 'h(conj(zeta_j))'],
                poles_D_minus=['conj(zeta_j)', 'h(zeta_j)'],
                h='-q0^2/z')


def three_panel_contour(name, eval_fn, xlim, ylim, xticks, yticks, q_max):
    """Draw a compact 1x3 filled-contour composite figure for (|q|, |p|, eta)."""
    x = np.linspace(xlim[0], xlim[1], 401)
    times = np.linspace(ylim[0], ylim[1], 251)
    X, TIME = np.meshgrid(x, times)
    q, p, eta, errors = eval_fn(X, TIME)
    if max(errors.values()) >= 1e-10:
        raise ValueError(errors)
    values = [abs(q), abs(p), eta.real]
    titles = [r'$|\widetilde q(x,t)|$', r'$|\widetilde p(x,t)|$', r'$\eta(x,t)$']
    bloch_lim = float(np.sqrt(.45))
    limits = [(0, q_max), (0, bloch_lim), (-bloch_lim, bloch_lim)]
    ticks = [[0, 2, 4], [0, .3, .6], [-.6, 0, .6]]
    cmaps = ['viridis', 'magma', 'RdBu_r']

    fig = plt.figure(figsize=(NUMERIC_WIDTH_IN, 68/25.4))
    grid = fig.add_gridspec(
        2, 3, height_ratios=[1, .055],
        left=.11, right=.98, bottom=.11, top=.87,
        wspace=.30, hspace=.58)
    panels, colorbars = [], []
    for j in range(3):
        ax = fig.add_subplot(grid[0, j])
        ax.set_rasterization_zorder(1.1)
        norm = matplotlib.colors.Normalize(*limits[j])
        im = ax.contourf(x, times, values[j], levels=np.linspace(*limits[j], 41),
                         cmap=cmaps[j], norm=norm)
        ax.set(xlabel=r'$x$', xlim=xlim, ylim=ylim)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        if j == 0:
            ax.set_ylabel(r'$t$')
        ax.tick_params(labelsize=9.5, length=3, width=.75, pad=3)
        ax.xaxis.label.set_size(10.5)
        ax.yaxis.label.set_size(10.5)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(.75)
        ax.set_title(titles[j], fontsize=10.5, pad=6)
        panel_label(ax, chr(97+j))
        panels.append(ax)

        cax = fig.add_subplot(grid[1, j])
        cax.set_rasterization_zorder(1.1)
        fig.colorbar(im, cax=cax, orientation='horizontal', ticks=ticks[j])
        cax.tick_params(labelsize=8.5, length=2.5, pad=2)
        colorbars.append(cax)

    export(fig, name, axes=panels, exclude=colorbars)
    return errors


def main(audit=False, only_numerics=False):
    global AUDIT, require_matplotlib_panel_alignment
    AUDIT = audit
    if AUDIT:
        sys.path.insert(0, str(SKILL_SCRIPTS))
        from audit_panel_alignment import require_matplotlib_panel_alignment
    QA.mkdir(parents=True, exist_ok=True)
    child_env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    (ROOT/'figures').mkdir(exist_ok=True)
    matplotlib.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial'],
        'font.size': 10.5, 'axes.labelsize': 10.5, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': .75, 'savefig.dpi': 600})

    previous_report = json.loads((ROOT/'calculations'/'plot-results.json').read_text(encoding='utf-8'))
    geometry_record = previous_report.get('spectral_geometry') if only_numerics else spectral_figures()

    one_errors = three_panel_contour(
        'one-soliton-fields', calc, (-6, 6), (-1.5, 1.5), [-6, 0, 6], [-1.5, 0, 1.5], 4.0)
    two_errors = three_panel_contour(
        'two-soliton-fields', calc_two_soliton, (-8, 8), (-2.5, 2.5), [-8, 0, 8], [-2, 0, 2], 4.2)

    # Cross-check double-precision grids against 70-digit mpmath on 5 points
    import mpmath as mp
    from verify_direct_rhp import fields as one_high_fields, ResidueSolution
    mp.mp.dps = 70
    points = [('-1.2', '-.35'), ('0', '0'), ('.37', '.21'), ('1.4', '-.1'), ('2.3', '.6')]
    one_maxima = np.zeros(3)
    for xx, tt in points:
        values = calc(float(xx), float(tt))[:3]
        high = one_high_fields(mp.mpf(xx), mp.mpf(tt))
        one_maxima = np.maximum(one_maxima, [abs(complex(v)-complex(h)) for v, h in zip(values, high)])
    assert max(one_maxima) < 1e-12, one_maxima

    q0_m, a_m, eta0_m, alpha_m, beta_m = [mp.mpf(s) for s in ('1', '2', '.6', '1', '.1')]
    Fm = lambda k: 4*beta_m*k**2-alpha_m*k-2*beta_m*q0_m**2+eta0_m/(a_m*(k-a_m))
    z1m, z2m = 2*q0_m*mp.exp(1j*mp.pi/3), 3*q0_m*mp.exp(2j*mp.pi/3)
    zm = [z1m, z2m, -q0_m**2/mp.conj(z1m), -q0_m**2/mp.conj(z2m)]
    sp = lambda idx: mp.fprod(zm[idx]-zj for j, zj in enumerate(zm) if j != idx)/mp.fprod(zm[idx]-mp.conj(zj) for zj in zm)
    c1m, c2m = 1/sp(0), 1/sp(1)
    cm = [c1m, c2m, -q0_m**2/mp.conj(z1m)**2*mp.conj(c1m), -q0_m**2/mp.conj(z2m)**2*mp.conj(c2m)]
    Omega_m = -alpha_m*q0_m**2+2*eta0_m/a_m
    two_maxima = np.zeros(3)
    for xx, tt in points:
        values = calc_two_soliton(float(xx), float(tt))[:3]
        sol_m = ResidueSolution(q0_m, zm, cm, mp.mpf(xx), mp.mpf(tt), Fm)
        Rm = sol_m.medium(a_m, eta0_m)
        ph = mp.exp(-1j*Omega_m*mp.mpf(tt))
        high = (ph*1j*sol_m.M1()[0, 1], -ph*Rm[0, 1], Rm[0, 0])
        two_maxima = np.maximum(two_maxima, [abs(complex(v)-complex(h)) for v, h in zip(values, high)])
    assert max(two_maxima) < 1e-12, two_maxima

    names = ['one-soliton-fields', 'two-soliton-fields']
    if not only_numerics:
        names = ['analytic-regions', 'discrete-spectrum']+names
    qa = {}
    if AUDIT:
        for name in names:
            pdf = ROOT/'figures'/(name+'.pdf')
            text = subprocess.run(
                [sys.executable, str(SKILL_SCRIPTS/'audit_pdf_text.py'), str(pdf), '--min-pt', '5', '--json'],
                capture_output=True, text=True, encoding='utf-8', errors='replace', env=child_env)
            text_record = json.loads(text.stdout)
            collision = subprocess.run(
                [sys.executable, str(SKILL_SCRIPTS/'audit_figure_collisions.py'), str(pdf),
                 '--json-out', str(QA/(name+'.collisions.json'))],
                capture_output=True, text=True, encoding='utf-8', errors='replace', env=child_env)
            collision_record = json.loads((QA/(name+'.collisions.json')).read_text(encoding='utf-8'))
            alignment_record = json.loads((QA/(name+'.alignment.json')).read_text(encoding='utf-8'))
            qa[name] = dict(
                alignment_verdict=alignment_record['verdict'],
                alignment_counts=alignment_record['summary'],
                text_exit=text.returncode,
                collision_exit=collision.returncode,
                minimum_glyph_pt=text_record['minimum_found_pt'],
                collision_verdict=collision_record['verdict'],
                collision_counts=collision_record['summary'])
            if text.returncode or collision.returncode:
                raise RuntimeError(f'{name}: PDF typography/collision audit failed')

    report = dict(
        method='Equilibrated direct residue solve (4x4 one-orbit and 8x8 two-orbit); infinity/resonance reconstruction',
        parameters=DEFAULT,
        figure_layout=dict(
            fields='1 row by 3 columns filled contours with per-panel horizontal colorbars: |q|, |p|, eta',
            width_mm=146.4, height_mm=68),
        one_soliton=dict(errors=one_errors, high_precision_comparison_q_p_eta=list(one_maxima)),
        two_soliton=dict(errors=two_errors, high_precision_comparison_q_p_eta=list(two_maxima)),
        spectral_geometry=geometry_record,
        files=['analytic-regions.pdf', 'discrete-spectrum.pdf', 'one-soliton-fields.pdf', 'two-soliton-fields.pdf'],
        figure_qa=qa)
    (ROOT/'calculations'/'plot-results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', action='store_true', help='run nature-figure alignment, text and collision audits')
    parser.add_argument('--only-numerics', action='store_true', help='preserve spectral PDFs and regenerate soliton figures')
    args = parser.parse_args()
    main(args.audit, args.only_numerics)
