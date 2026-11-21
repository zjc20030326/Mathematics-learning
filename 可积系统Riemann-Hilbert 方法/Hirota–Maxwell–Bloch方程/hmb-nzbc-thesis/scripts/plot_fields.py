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
DEFAULT = dict(q0=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5)


def dispersion(k, q0, a, eta0, alpha, beta):
    return 4*beta*k*k-alpha*k-2*beta*q0*q0+eta0/(a*(k-a))


def calc(X, TIME, *, q0=1., a=2., eta0=.6, alpha=1., beta=.1, b=1.5):
    """Return q,p,eta and diagnostics from batched residue-system solves."""
    if q0<=0 or a==0 or eta0==0 or b<=q0:
        raise ValueError('Require q0>0, a!=0, eta0!=0 and b>q0.')
    X, TIME = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(TIME, dtype=float))
    shape = X.shape
    kappa = np.sqrt(b*b-q0*q0)
    zeta = np.array([1j*(b+kappa), -1j*(b-kappa)])
    c1 = 4j*kappa/(1-q0*q0/zeta[0]**2)
    constants = np.array([c1, -q0*q0/np.conj(zeta[0])**2*np.conj(c1)])
    k = (zeta-q0*q0/zeta)/2
    lam = (zeta+q0*q0/zeta)/2
    F = dispersion(k, q0, a, eta0, alpha, beta)
    C = constants*np.exp(-2j*lam*(X[..., None]-F*TIME[..., None]))
    A = np.conj(C)[..., None, :]/(zeta[:, None]-np.conj(zeta)[None, :])
    B = C[..., None, :]/(np.conj(zeta)[:, None]-zeta[None, :])
    L = np.broadcast_to(np.eye(4, dtype=complex), shape+(4, 4)).copy()
    L[..., :2, 2:] = A
    L[..., 2:, :2] = -B
    f = np.column_stack((np.ones(2), -1j*q0/zeta))
    g = np.column_stack((-1j*q0/np.conj(zeta), np.ones(2)))
    rhs = np.broadcast_to(np.vstack((f, g)), shape+(4, 2))
    # Row equilibration preserves the system and limits the grid's dynamic range.
    scale = np.max(np.abs(L), axis=-1)
    answer = np.linalg.solve(L/scale[..., :, None], rhs/scale[..., :, None])
    v, u = answer[..., :2, :], answer[..., 2:, :]
    field = q0+1j*np.sum(C*v[..., :, 0], axis=-1)
    za = a+np.sqrt(a*a+q0*q0)
    lambda_a = (za+q0*q0/za)/2
    M = np.broadcast_to(np.array([[1, -1j*q0/za], [-1j*q0/za, 1]]), shape+(2, 2)).copy()
    M[..., :, 1] += np.sum(C[..., :, None]*v/(za-zeta)[..., None], axis=-2)
    M[..., :, 0] -= np.sum(np.conj(C)[..., :, None]*u/(za-np.conj(zeta))[..., None], axis=-2)
    R = M@np.diag([eta0*lambda_a/a, -eta0*lambda_a/a])@np.linalg.inv(M)
    phase = np.exp(-1j*(-alpha*q0*q0+2*eta0/a)*TIME)
    q, p, eta = phase*field, -phase*R[..., 0, 1], R[..., 0, 0]
    linear = np.max(np.abs(L@answer-rhs)/(1+np.abs(L)@np.abs(answer)+np.abs(rhs)))
    bloch = np.max(np.abs(np.abs(p)**2+eta**2-eta0**2*(1+q0*q0/a**2)))
    return q, p, eta, dict(scaled_linear_residual=float(linear),
                          Bloch_error=float(bloch), eta_imaginary=float(np.max(abs(eta.imag))))


def velocity_period(alpha=1., beta=.1):
    q0, a, eta0, b = (DEFAULT[key] for key in ('q0', 'a', 'eta0', 'b'))
    kappa = np.sqrt(b*b-q0*q0)
    F = dispersion(1j*b, q0, a, eta0, alpha, beta)
    return float(F.real), float(F.imag), None if abs(F.imag)<1e-14 else float(np.pi/(kappa*abs(F.imag)))


# Figure contract: spectral figures show the analyticity domains and one complete
# orbit; the field figure combines height and localization; dynamics separate
# translation from internal oscillation. All values below are deterministic
# evaluations of the residue system, not observations or a PDE discretization.
import os
import argparse
import sys
import tempfile
import subprocess
from matplotlib.patches import Rectangle, Wedge
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

QA = Path(os.environ.get('HMB_FIGURE_QA', Path(tempfile.gettempdir())/'hmb-zhang-fan-figure-qa'))
AUDIT = False
SKILL_SCRIPTS = Path(os.environ.get('NATURE_FIGURE_SCRIPTS',
    Path.home()/'.agents'/'skills'/'nature-figure'/'scripts'))
# Optional QA tools are loaded only with --audit; plotting itself is portable.


def panel_label(ax, letter):
    transform = ax.transAxes+ScaledTranslation(-10/72, 4/72, ax.figure.dpi_scale_trans)
    if hasattr(ax, 'text2D'):
        ax.text2D(0, 1, letter, transform=transform, fontsize=11, fontweight='bold', va='bottom')
    else:
        ax.text(0, 1, letter, transform=transform, fontsize=11, fontweight='bold', va='bottom')


def export(fig, name, *, axes=None, exclude=None):
    QA.mkdir(parents=True,exist_ok=True)
    kwargs = {'exclude_axes': exclude or []}
    if axes is not None:
        kwargs['axes'] = axes
    fig.canvas.draw()
    if AUDIT:
        require_matplotlib_panel_alignment(fig,
            json_out=QA/(name+'.alignment.json'), tolerance_pt=1.5,
            gutter_tolerance_pt=1.5, strict=True, **kwargs)
    fig.savefig(ROOT/'figures'/(name+'.pdf'), dpi=600)
    if AUDIT: fig.savefig(QA/(name+'.svg'))
    # Backend-native previews live outside the manuscript project.
    fig.savefig(QA/(name+'.png'), dpi=300)
    plt.close(fig)


def spectral_base(*, discrete=False):
    """Symbolic HMB z-plane, with D+ on the left of every contour arc.

    Coordinates serve only to position the schematic: no numerical spectral
    parameters are asserted. Unlike GI, HMB analyticity is set by Im(lambda),
    so its continuous spectrum is a real axis and a circle, not two axes.
    """
    # At 0.86 of the 146.4 mm text width, 14 pt prints at about 11.23 pt.
    label_size = 14
    fig, ax = plt.subplots(figsize=(6.18, 4.2))
    fig.subplots_adjust(left=.06, right=.98, bottom=.07, top=.97)
    ax.set_aspect('equal', adjustable='box')
    ax.set(xlim=(-2.75, 2.75), ylim=(-2.10, 2.10))
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    blue = '#E1EDF5'
    ax.add_patch(Rectangle((-2.50, 0), 5., 1.82, facecolor=blue,
                           edgecolor='none', zorder=-3))
    ax.add_patch(Wedge((0, 0), 1, 0, 180, facecolor='white',
                       edgecolor='none', zorder=-2))
    ax.add_patch(Wedge((0, 0), 1, 180, 360, facecolor=blue,
                       edgecolor='none', zorder=-2))
    # The imaginary axis is a coordinate guide, never a HMB jump contour.
    ax.plot([0, 0], [-1.85, 1.87], color='#8A8A8A', lw=.6,
            ls=(0, (3, 4)), zorder=-1)
    ax.annotate('', (0, 1.96), (0, 1.76),
                arrowprops=dict(arrowstyle='->', lw=.7, color='#696969'))
    ax.plot([-2.55, 2.56], [0, 0], color='#242424', lw=.85)
    ang = np.linspace(0, 2*np.pi, 801)
    ax.plot(np.cos(ang), np.sin(ang), color='#242424', lw=.85)
    for x0, x1 in [(-2.25, -1.85), (1.85, 2.25), (-.18, -.57), (.57, .18)]:
        ax.annotate('', (x1, 0), (x0, 0),
                    arrowprops=dict(arrowstyle='-|>', lw=.85, color='#242424'))
    # Upper arc clockwise; lower arc counterclockwise.
    for th0, th1 in [(133, 125), (227, 235)]:
        a0, a1 = np.deg2rad([th0, th1])
        ax.annotate('', (np.cos(a1), np.sin(a1)), (np.cos(a0), np.sin(a0)),
                    arrowprops=dict(arrowstyle='-|>', lw=.85, color='#242424',
                                    connectionstyle='arc3,rad=.02'))
    for x, y, text in [(-1.76, 1.53, r'$D_+$'), (-1.76, -1.53, r'$D_-$'),
                        (.43, -.43, r'$D_+$'), (.43, .43, r'$D_-$')]:
        ax.text(x, y, text, ha='center', va='center', fontsize=label_size)
    ax.text(2.40, .14, r'$\operatorname{Re}z$', ha='center', fontsize=label_size)
    ax.text(.14, 1.94, r'$\operatorname{Im}z$', va='center', fontsize=label_size)
    ax.text(.10, -.18, r'$0$', ha='left', va='top', fontsize=label_size)
    # Intersections of contour arcs and branch points have different symbols.
    ax.plot([-1, 1], [0, 0], ls='', marker='x', ms=5.5, color='#242424')
    ax.text(-1.20, -.13, r'$-q_0$', ha='center', va='top', fontsize=label_size)
    ax.text(1.20, -.13, r'$q_0$', ha='center', va='top', fontsize=label_size)
    ax.plot([0, 0], [-1, 1], ls='', marker='s', ms=4.7,
            mfc='white', mec='#6D4C64', mew=.85)
    ax.text(.12, 1.04, r'$\mathrm{i}q_0$', va='bottom', fontsize=label_size)
    ax.text(.12, -1.05, r'$-\mathrm{i}q_0$', va='top', fontsize=label_size)
    return fig, ax


def spectral_figures():
    label_size = 14
    fig, ax = spectral_base()
    export(fig, 'analytic-regions')
    fig, ax = spectral_base(discrete=True)
    # One general complex quartet, not the later pure-imaginary specialization.
    zeta = 1.54+1.30j
    hzeta = -1/zeta
    plus = [zeta, np.conj(hzeta)]
    minus = [np.conj(zeta), hzeta]
    ax.plot(np.real(plus), np.imag(plus), ls='', marker='o', ms=5.5,
            mfc='#2B638A', mec='#2B638A', mew=.9)
    ax.plot(np.real(minus), np.imag(minus), ls='', marker='o', ms=5.5,
            mfc='white', mec='#242424', mew=.95)
    ax.text(zeta.real+.14, zeta.imag+.04, r'$\zeta_j$', fontsize=label_size,
            ha='left', va='bottom')
    ax.text(zeta.real+.14, -zeta.imag-.04, r'$(\zeta_j)^*$', fontsize=label_size,
            ha='left', va='top')
    ax.annotate(r'$h(\zeta_j)$', (hzeta.real, hzeta.imag), (-1.72, .50),
                fontsize=label_size, ha='center', va='bottom',
                arrowprops=dict(arrowstyle='-', lw=.55, color='#5E5E5E'))
    ax.annotate(r'$h(\zeta_j^*)$', (hzeta.real, -hzeta.imag), (-1.72, -.50),
                fontsize=label_size, ha='center', va='top',
                arrowprops=dict(arrowstyle='-', lw=.55, color='#5E5E5E'))
    export(fig, 'discrete-spectrum')
    return dict(type='symbolic general complex quartet',
                poles_D_plus=['zeta_j', 'h(conj(zeta_j))'],
                poles_D_minus=['conj(zeta_j)', 'h(zeta_j)'],
                h='-q0^2/z',
                pure_imaginary_orbit='specialization used only in the one-orbit solution',
                contour_direction='real exterior right; real interior left; upper circle clockwise; lower circle counterclockwise',
                domain_color='D+ pale blue, D- white')


# Numerical artwork is drawn at the actual thesis text width (146.4 mm).
# Each colorbar compares only the same field; the surfaces and maps in a row
# use identical limits. All grids remain analytic residue-system evaluations.
NUMERIC_WIDTH_IN = 146.4/25.4
NUMERIC_COLORS = ['#276A8D', '#B56A32', '#5B7E53']
NUMERIC_STYLES = ['-', '--', ':']


def numeric_axes(ax):
    ax.tick_params(labelsize=10, length=3, width=.75, pad=3)
    ax.xaxis.label.set_size(11)
    ax.yaxis.label.set_size(11)
    for spine in ax.spines.values():
        spine.set_linewidth(.75)


def field_figure(X, TIME, x, times):
    """Three rows: optical amplitude, polarization amplitude and inversion.

    Left: three-dimensional surfaces on the original 121 by 101 analytic grid.
    Right: contour maps on the original 401 by 251 grid, with per-field scales.
    """
    q, p, eta, errors = calc(X, TIME)
    if max(errors.values()) >= 1e-10:
        raise ValueError(errors)
    xs, ts = np.linspace(-6, 6, 121), np.linspace(-1.5, 1.5, 101)
    XS, TS = np.meshgrid(xs, ts)
    qs, ps, es, _ = calc(XS, TS)
    values = [abs(q), abs(p), eta.real]
    surface_values = [abs(qs), abs(ps), es.real]
    titles = [r'$|\widetilde q(x,t)|$', r'$|\widetilde p(x,t)|$', r'$\eta(x,t)$']
    limits = [(0, 4), (0, np.sqrt(.45)), (-np.sqrt(.45), np.sqrt(.45))]
    ticks = [[0, 2, 4], [0, .3, .6], [-.6, 0, .6]]
    cmaps = ['viridis', 'magma', 'RdBu_r']
    fig = plt.figure(figsize=(NUMERIC_WIDTH_IN, 165/25.4))
    grid = fig.add_gridspec(3, 3, width_ratios=[1, 1, .045],
                           left=.075, right=.925, bottom=.065, top=.945,
                           wspace=.36, hspace=.58)
    panels, colorbars = [], []
    for j in range(3):
        norm = matplotlib.colors.Normalize(*limits[j])
        ax = fig.add_subplot(grid[j, 0], projection='3d')
        ax.plot_surface(XS, TS, surface_values[j], rstride=1, cstride=1,
                        cmap=cmaps[j], norm=norm, linewidth=0,
                        antialiased=False, alpha=1)
        ax.set(xlim=(-6, 6), ylim=(-1.5, 1.5), zlim=limits[j])
        ax.set_xticks([-6, 0, 6]); ax.set_yticks([-1, 0, 1])
        ax.set_zticks(ticks[j])
        ax.grid(False)
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis.pane.fill = False
            axis.pane.set_edgecolor('#E0E0E0')
        ax.set_xlabel(r'$x$', labelpad=1, fontsize=11)
        ax.set_ylabel(r'$t$', labelpad=1, fontsize=11)
        ax.tick_params(pad=1, labelsize=10)
        ax.view_init(elev=27, azim=-58)
        ax.set_box_aspect((1.15, 1, .72))
        # Keep the projected z-ticks clear of the adjacent map's time label.
        position = ax.get_position()
        ax.set_position([position.x0-.035, position.y0, position.width, position.height])
        ax.set_title(titles[j], fontsize=11, pad=8)
        panel_label(ax, chr(97+2*j))
        panels.append(ax)
        contour = fig.add_subplot(grid[j, 1])
        # Rasterize only the filled contour layer to avoid PDF stitching seams.
        contour.set_rasterization_zorder(1.1)
        levels = np.linspace(*limits[j], 41)
        im = contour.contourf(x, times, values[j], levels=levels,
                              cmap=cmaps[j], norm=norm)
        contour.set(xlabel=r'$x$', ylabel=r'$t$', xlim=(-6, 6), ylim=(-1.5, 1.5))
        contour.set_xticks([-6, 0, 6]); contour.set_yticks([-1.5, 0, 1.5])
        numeric_axes(contour)
        panel_label(contour, chr(98+2*j))
        panels.append(contour)
        cax = fig.add_subplot(grid[j, 2])
        cax.set_rasterization_zorder(1.1)
        fig.colorbar(im, cax=cax, ticks=ticks[j])
        cax.tick_params(labelsize=9.5, length=2.5, pad=3)
        colorbars.append(cax)
    export(fig, 'one-soliton-fields', axes=panels, exclude=colorbars)
    return errors


def dynamics_figures(X, TIME, x, times):
    """Separate comoving shape/oscillation from beta-controlled translation."""
    speed, _, period = velocity_period()
    profile_x = np.linspace(-4, 4, 401)
    fig, axes = plt.subplots(1, 2, figsize=(NUMERIC_WIDTH_IN, 88/25.4))
    fig.subplots_adjust(left=.105, right=.97, top=.90, bottom=.34, wspace=.43)
    labels = [r'$0$', r'$0.25$', r'$0.5$']
    for frac, ls, color, label in zip([0., .25, .5], NUMERIC_STYLES, NUMERIC_COLORS, labels):
        tt = frac*period
        qq, _, _, _ = calc(profile_x+speed*tt, tt)
        axes[0].plot(profile_x, abs(qq), ls, color=color, lw=1.4, label=label)
    axes[0].set(xlabel=r'$X=x-A_bt$', ylabel=r'$|\widetilde q|$',
                xlim=(-4, 4), ylim=(-.08, 4.15))
    axes[0].set_xticks([-4, -2, 0, 2, 4]); axes[0].set_yticks([0, 1, 2, 3, 4])
    axes[0].legend(loc='upper center', bbox_to_anchor=(.5, -.25), ncol=3,
                   title=r'$t/T_t$', fontsize=10, title_fontsize=10.5,
                   handlelength=1.8, handletextpad=.35, columnspacing=.9,
                   borderaxespad=0, labelspacing=.5)
    times_line = np.linspace(-2, 2, 501)
    stationary_alpha = -DEFAULT['eta0']/(DEFAULT['a']*(DEFAULT['a']**2+DEFAULT['b']**2))
    periodic_cases = []
    for al, color, ls in zip([stationary_alpha, .5, 1.], NUMERIC_COLORS, NUMERIC_STYLES):
        aa, bb, tt = velocity_period(alpha=al)
        qq, _, _, _ = calc(aa*times_line, times_line, alpha=al)
        axes[1].plot(times_line, abs(qq), ls, color=color, lw=1.4, label=rf'${al:g}$')
        periodic_cases.append(dict(alpha=al, A_b=aa, B_b=bb, T_t=tt))
    axes[1].set(xlabel=r'$t$', ylabel=r'$|\widetilde q(A_bt,t)|$',
                xlim=(-2, 2), ylim=(1.9, 4.15))
    axes[1].set_xticks([-2, -1, 0, 1, 2]); axes[1].set_yticks([2, 3, 4])
    axes[1].legend(loc='upper center', bbox_to_anchor=(.5, -.25), ncol=3,
                   title=r'$\alpha$', fontsize=10, title_fontsize=10.5,
                   handlelength=1.8, handletextpad=.35, columnspacing=.9,
                   borderaxespad=0, labelspacing=.5)
    for j, ax in enumerate(axes):
        numeric_axes(ax)
        panel_label(ax, chr(97+j))
    export(fig, 'one-soliton-dynamics', axes=list(axes))

    fig, axes = plt.subplots(1, 2, figsize=(NUMERIC_WIDTH_IN, 78/25.4))
    fig.subplots_adjust(left=.10, right=.87, top=.85, bottom=.20, wspace=.36)
    velocity_cases = []
    for j, (ax, be) in enumerate(zip(axes, [0., .15])):
        qq, _, _, errs = calc(X, TIME, beta=be)
        aa, bb, tt = velocity_period(beta=be)
        ax.set_rasterization_zorder(1.1)
        im = ax.contourf(x, times, abs(qq), levels=np.linspace(0, 4, 41), cmap='viridis')
        ax.plot(aa*times, times, color='white', ls='--', lw=1.1)
        ax.set(xlabel=r'$x$', ylabel=r'$t$', xlim=(-6, 6), ylim=(-1.5, 1.5))
        ax.set_xticks([-6, -3, 0, 3, 6]); ax.set_yticks([-1.5, 0, 1.5])
        ax.set_title(rf'$\beta={be:g},\ A_b={aa:g}$', fontsize=10.5, pad=9)
        numeric_axes(ax)
        panel_label(ax, chr(97+j))
        if max(errs.values()) >= 1e-10:
            raise ValueError(errs)
        velocity_cases.append(dict(beta=be, A_b=aa, B_b=bb, T_t=tt, errors=errs))
    cax = fig.add_axes([.905, .20, .022, .65])
    cax.set_rasterization_zorder(1.1)
    fig.colorbar(im, cax=cax, ticks=[0, 1, 2, 3, 4])
    cax.tick_params(labelsize=10, length=3, pad=3)
    cax.set_title(r'$|\widetilde q|$', fontsize=11, pad=9)
    export(fig, 'one-soliton-velocity', axes=list(axes), exclude=[cax])
    return periodic_cases, velocity_cases


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


def calc_two_soliton(X, TIME, *, q0=1., a=2., eta0=.6, alpha=1., beta=.1):
    """Evaluate the general complex two-orbit (m=4) reflectionless RHP on plotting grids."""
    X, TIME = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(TIME, dtype=float))
    shape = X.shape
    zeta, constants = two_soliton_data(q0)
    k = (zeta-q0*q0/zeta)/2
    lam = (zeta+q0*q0/zeta)/2
    F = dispersion(k, q0, a, eta0, alpha, beta)
    C = constants*np.exp(-2j*lam*(X[..., None]-F*TIME[..., None]))
    m = 4
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


def two_soliton_figures():
    """Generate two-soliton 3-field surface/contour figure and collision interaction figure."""
    x, times = np.linspace(-8, 8, 401), np.linspace(-2.5, 2.5, 251)
    X, TIME = np.meshgrid(x, times)
    q, p, eta, errors = calc_two_soliton(X, TIME)
    if max(errors.values()) >= 1e-10:
        raise ValueError(errors)
    xs, ts = np.linspace(-8, 8, 121), np.linspace(-2.5, 2.5, 101)
    XS, TS = np.meshgrid(xs, ts)
    qs, ps, es, _ = calc_two_soliton(XS, TS)
    values = [abs(q), abs(p), eta.real]
    surface_values = [abs(qs), abs(ps), es.real]
    titles = [r'$|\widetilde q(x,t)|$', r'$|\widetilde p(x,t)|$', r'$\eta(x,t)$']
    limits = [(0, 4.2), (0, np.sqrt(.45)), (-np.sqrt(.45), np.sqrt(.45))]
    ticks = [[0, 2, 4], [0, .3, .6], [-.6, 0, .6]]
    cmaps = ['viridis', 'magma', 'RdBu_r']
    fig = plt.figure(figsize=(NUMERIC_WIDTH_IN, 165/25.4))
    grid = fig.add_gridspec(3, 3, width_ratios=[1, 1, .045],
                           left=.075, right=.925, bottom=.065, top=.945,
                           wspace=.36, hspace=.58)
    panels, colorbars = [], []
    for j in range(3):
        norm = matplotlib.colors.Normalize(*limits[j])
        ax = fig.add_subplot(grid[j, 0], projection='3d')
        ax.plot_surface(XS, TS, surface_values[j], rstride=1, cstride=1,
                        cmap=cmaps[j], norm=norm, linewidth=0,
                        antialiased=False, alpha=1)
        ax.set(xlim=(-8, 8), ylim=(-2.5, 2.5), zlim=limits[j])
        ax.set_xticks([-8, 0, 8]); ax.set_yticks([-2, 0, 2])
        ax.set_zticks(ticks[j])
        ax.grid(False)
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis.pane.fill = False
            axis.pane.set_edgecolor('#E0E0E0')
        ax.set_xlabel(r'$x$', labelpad=1, fontsize=11)
        ax.set_ylabel(r'$t$', labelpad=1, fontsize=11)
        ax.tick_params(pad=1, labelsize=10)
        ax.view_init(elev=27, azim=-58)
        ax.set_box_aspect((1.15, 1, .72))
        position = ax.get_position()
        ax.set_position([position.x0-.035, position.y0, position.width, position.height])
        ax.set_title(titles[j], fontsize=11, pad=8)
        panel_label(ax, chr(97+2*j))
        panels.append(ax)
        contour = fig.add_subplot(grid[j, 1])
        contour.set_rasterization_zorder(1.1)
        levels = np.linspace(*limits[j], 41)
        im = contour.contourf(x, times, values[j], levels=levels,
                              cmap=cmaps[j], norm=norm)
        contour.set(xlabel=r'$x$', ylabel=r'$t$', xlim=(-8, 8), ylim=(-2.5, 2.5))
        contour.set_xticks([-8, 0, 8]); contour.set_yticks([-2, 0, 2])
        numeric_axes(contour)
        panel_label(contour, chr(98+2*j))
        panels.append(contour)
        cax = fig.add_subplot(grid[j, 2])
        cax.set_rasterization_zorder(1.1)
        fig.colorbar(im, cax=cax, ticks=ticks[j])
        cax.tick_params(labelsize=9.5, length=2.5, pad=3)
        colorbars.append(cax)
    export(fig, 'two-soliton-fields', axes=panels, exclude=colorbars)

    # Second figure: time snapshots (a) and collision trajectories with shift (b)
    q0, a, eta0, alpha, beta = (DEFAULT[key] for key in ('q0', 'a', 'eta0', 'alpha', 'beta'))
    zeta, _ = two_soliton_data(q0)
    k = (zeta-q0*q0/zeta)/2
    lam = (zeta+q0*q0/zeta)/2
    F = dispersion(k, q0, a, eta0, alpha, beta)
    v1 = float((lam[0]*F[0]).imag/lam[0].imag)
    v2 = float((lam[1]*F[1]).imag/lam[1].imag)
    s1_of = lambda z: ((z-zeta[0])/(z-np.conj(zeta[0])))*((z-zeta[2])/(z-np.conj(zeta[2])))
    s2_of = lambda z: ((z-zeta[1])/(z-np.conj(zeta[1])))*((z-zeta[3])/(z-np.conj(zeta[3])))
    dx1 = float(np.log(abs(s2_of(zeta[0])))/lam[0].imag)
    dx2 = float(-np.log(abs(s1_of(zeta[1])))/lam[1].imag)

    fig, axes = plt.subplots(1, 2, figsize=(NUMERIC_WIDTH_IN, 88/25.4))
    fig.subplots_adjust(left=.10, right=.87, top=.90, bottom=.30, wspace=.38)
    prof_x = np.linspace(-8, 8, 601)
    for tt, ls, color, lbl in zip([-2., 0., 2.], ['--', '-', ':'], NUMERIC_COLORS, [r'$-2$', r'$0$', r'$2$']):
        qq, _, _, _ = calc_two_soliton(prof_x, tt)
        axes[0].plot(prof_x, abs(qq), ls, color=color, lw=1.35, label=lbl)
    axes[0].set(xlabel=r'$x$', ylabel=r'$|\widetilde q(x,t)|$', xlim=(-8, 8), ylim=(-.08, 4.35))
    axes[0].set_xticks([-8, -4, 0, 4, 8]); axes[0].set_yticks([0, 1, 2, 3, 4])
    axes[0].legend(loc='upper center', bbox_to_anchor=(.5, -.26), ncol=3,
                   title=r'$t$', fontsize=10, title_fontsize=10.5,
                   handlelength=1.8, handletextpad=.35, columnspacing=.9,
                   borderaxespad=0, labelspacing=.5)
    numeric_axes(axes[0])
    panel_label(axes[0], 'a')

    axes[1].set_rasterization_zorder(1.1)
    im = axes[1].contourf(x, times, abs(q), levels=np.linspace(0, 4.2, 43), cmap='viridis')
    t_neg = np.linspace(-2.5, 0., 100)
    t_pos = np.linspace(0., 2.5, 100)
    axes[1].plot(v1*t_neg-dx1/2, t_neg, color='white', ls='--', lw=1.1, label='pre-collision')
    axes[1].plot(v2*t_neg-dx2/2, t_neg, color='white', ls='--', lw=1.1)
    axes[1].plot(v1*t_pos+dx1/2, t_pos, color='#FFB347', ls='-', lw=1.15, label='post-collision')
    axes[1].plot(v2*t_pos+dx2/2, t_pos, color='#FFB347', ls='-', lw=1.15)
    axes[1].set(xlabel=r'$x$', ylabel=r'$t$', xlim=(-8, 8), ylim=(-2.5, 2.5))
    axes[1].set_xticks([-8, -4, 0, 4, 8]); axes[1].set_yticks([-2, 0, 2])
    numeric_axes(axes[1])
    panel_label(axes[1], 'b')
    # Add compact legend below panel (b) matching panel (a)
    h1 = Line2D([0], [0], color='#555555', ls='--', lw=1.2)
    h2 = Line2D([0], [0], color='#D97706', ls='-', lw=1.3)
    axes[1].legend([h1, h2], [r'$t<0$', r'$t>0$'], loc='upper center', bbox_to_anchor=(.5, -.26),
                   ncol=2, title=r'$x_n^{(\pm)}(t)$', fontsize=9.5, title_fontsize=10.5,
                   handlelength=1.8, handletextpad=.35, columnspacing=1.0, borderaxespad=0)
    cax = fig.add_axes([.90, .30, .022, .60])
    cax.set_rasterization_zorder(1.1)
    fig.colorbar(im, cax=cax, ticks=[0, 1, 2, 3, 4])
    cax.tick_params(labelsize=9.5, length=2.5, pad=3)
    cax.set_title(r'$|\widetilde q|$', fontsize=10.5, pad=8)
    export(fig, 'two-soliton-interaction', axes=list(axes), exclude=[cax])
    return dict(errors=errors, v1=v1, v2=v2, Delta_x1=dx1, Delta_x2=dx2)


def main(audit=False, only_numerics=False, only_two_soliton=False):
    global AUDIT, require_matplotlib_panel_alignment
    AUDIT = audit
    if AUDIT:
        sys.path.insert(0,str(SKILL_SCRIPTS))
        from audit_panel_alignment import require_matplotlib_panel_alignment
    QA.mkdir(parents=True, exist_ok=True)
    child_env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    if AUDIT:
        source_check = subprocess.run(
            [sys.executable, str(SKILL_SCRIPTS/'validate_figure.py'), str(Path(__file__)), '--json'],
            capture_output=True, text=True, encoding='utf-8', errors='replace', env=child_env)
        source_record = json.loads(source_check.stdout)
        source_record['source'] = 'scripts/plot_fields.py'
        source_record['manual_review'] = ('Numerical canvases are drawn at the thesis text width (146.4 mm); spectral canvases retain their separate placement contract. '
            'the two production-column/TIFF warnings do not apply to vector thesis artwork. '
            'SVG sidecars and native PNG previews are retained with the QA.')
        (QA/'source-preflight.json').write_text(json.dumps(source_record, indent=2), encoding='utf-8')
        if source_check.returncode:
            raise RuntimeError(source_check.stdout)
    (ROOT/'figures').mkdir(exist_ok=True)
    matplotlib.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans','Arial'],
        'font.size':11, 'axes.labelsize':11, 'xtick.labelsize':10, 'ytick.labelsize':10,
        'legend.fontsize':10, 'pdf.fonttype':42, 'svg.fonttype':'none',
        'axes.spines.top':False, 'axes.spines.right':False, 'axes.linewidth':.7,
        'legend.frameon':False, 'savefig.dpi':600})
    previous_report = json.loads((ROOT/'calculations'/'plot-results.json').read_text(encoding='utf-8'))
    geometry_record = previous_report.get('spectral_geometry') if (only_numerics or only_two_soliton) else spectral_figures()
    if only_two_soliton:
        base = previous_report.get('base_errors')
        alpha_cases = previous_report.get('dynamics', {}).get('alpha_cases')
        beta_cases = previous_report.get('dynamics', {}).get('beta_cases')
        maxima = np.array(previous_report.get('high_precision_comparison_q_p_eta', [0., 0., 0.]))
        names = ['two-soliton-fields', 'two-soliton-interaction']
    else:
        x, times = np.linspace(-6,6,401), np.linspace(-1.5,1.5,251)
        X, TIME = np.meshgrid(x,times)
        base = field_figure(X,TIME,x,times)
        alpha_cases, beta_cases = dynamics_figures(X,TIME,x,times)
        import mpmath as mp
        from verify_direct_rhp import fields
        mp.mp.dps=70
        points=[('-1.2','-.35'),('0','0'),('.37','.21'),('1.4','-.1'),('2.3','.6')]
        maxima=np.zeros(3)
        for xx,tt in points:
            values=calc(float(xx),float(tt))[:3]
            high=fields(mp.mpf(xx),mp.mpf(tt))
            maxima=np.maximum(maxima,[abs(complex(v)-complex(h)) for v,h in zip(values,high)])
        assert max(maxima)<1e-12
        names=['one-soliton-fields','one-soliton-dynamics','one-soliton-velocity','two-soliton-fields','two-soliton-interaction']
        if not only_numerics:
            names=['analytic-regions','discrete-spectrum']+names
    two_sol = two_soliton_figures()
    # Cross-check double-precision two-soliton evaluator against 70-digit mpmath
    import mpmath as mp
    from verify_direct_rhp import ResidueSolution
    mp.mp.dps = 70
    q0_m, a_m, eta0_m, alpha_m, beta_m = [mp.mpf(s) for s in ('1', '2', '.6', '1', '.1')]
    Fm = lambda k: 4*beta_m*k**2-alpha_m*k-2*beta_m*q0_m**2+eta0_m/(a_m*(k-a_m))
    z1m, z2m = 2*q0_m*mp.exp(1j*mp.pi/3), 3*q0_m*mp.exp(2j*mp.pi/3)
    zm = [z1m, z2m, -q0_m**2/mp.conj(z1m), -q0_m**2/mp.conj(z2m)]
    sp = lambda idx: mp.fprod(zm[idx]-zj for j, zj in enumerate(zm) if j != idx)/mp.fprod(zm[idx]-mp.conj(zj) for zj in zm)
    c1m, c2m = 1/sp(0), 1/sp(1)
    cm = [c1m, c2m, -q0_m**2/mp.conj(z1m)**2*mp.conj(c1m), -q0_m**2/mp.conj(z2m)**2*mp.conj(c2m)]
    Omega_m = -alpha_m*q0_m**2+2*eta0_m/a_m
    points = [('-1.2', '-.35'), ('0', '0'), ('.37', '.21'), ('1.4', '-.1'), ('2.3', '.6')]
    two_maxima = np.zeros(3)
    for xx, tt in points:
        values = calc_two_soliton(float(xx), float(tt))[:3]
        sol_m = ResidueSolution(q0_m, zm, cm, mp.mpf(xx), mp.mpf(tt), Fm)
        Rm = sol_m.medium(a_m, eta0_m)
        ph = mp.exp(-1j*Omega_m*mp.mpf(tt))
        high = (ph*1j*sol_m.M1()[0, 1], -ph*Rm[0, 1], Rm[0, 0])
        two_maxima = np.maximum(two_maxima, [abs(complex(v)-complex(h)) for v, h in zip(values, high)])
    assert max(two_maxima) < 1e-12, two_maxima
    two_sol['high_precision_comparison_q_p_eta'] = list(two_maxima)
    qa={}
    child_env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    if AUDIT:
        for name in names:
            pdf=ROOT/'figures'/(name+'.pdf')
            text=subprocess.run([sys.executable,str(SKILL_SCRIPTS/'audit_pdf_text.py'),str(pdf),'--min-pt','5','--json'],
                                capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env)
            text_record = json.loads(text.stdout)
            text_record['pdf'] = pdf.name
            (QA/(name+'.text.json')).write_text(json.dumps(text_record,indent=2),encoding='utf-8')
            collision=subprocess.run([sys.executable,str(SKILL_SCRIPTS/'audit_figure_collisions.py'),str(pdf),
                                     '--json-out',str(QA/(name+'.collisions.json'))],capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env)
            collision_record=json.loads((QA/(name+'.collisions.json')).read_text(encoding='utf-8'))
            collision_record['pdf']=pdf.name
            if collision_record['summary']['warn']:
                collision_record['manual_review']='Pending explicit review of the final PDF findings; see figure acceptance notes.'
            (QA/(name+'.collisions.json')).write_text(json.dumps(collision_record,indent=2),encoding='utf-8')
            alignment_record=json.loads((QA/(name+'.alignment.json')).read_text(encoding='utf-8'))
            qa[name]=dict(alignment_verdict=alignment_record['verdict'],
                          alignment_counts=alignment_record['summary'],
                          text_exit=text.returncode,collision_exit=collision.returncode,
                          minimum_glyph_pt=text_record['minimum_found_pt'],
                          collision_verdict=collision_record['verdict'],
                          collision_counts=collision_record['summary'],
                          manual_review=collision_record.get('manual_review','No collision warnings.'))
            if text.returncode or collision.returncode:
                raise RuntimeError(f'{name}: PDF typography/collision audit failed; inspect QA records.')
    report=dict(method='Equilibrated direct residue solve (4x4 one-orbit and 8x8 two-orbit); infinity/resonance reconstruction',
        variables='x,t; k original spectral parameter; lambda background root; z uniformization; output q,p are original tilded fields',
        parameters=DEFAULT,grid=dict(x_min=-6,x_max=6,x_count=401,t_min=-1.5,t_max=1.5,t_count=251),
        surface_grid=dict(x_count=121,t_count=101,rule='separate analytic sampling grid for vector surface rendering'),
        figure_layout=dict(fields='3 rows by 2 columns; optical, polarization, inversion; surfaces left and contour maps right',
            width_mm=146.4,heights_mm=dict(fields=165,dynamics=88,velocity=78,two_fields=165,two_interaction=88),
            field_color_limits=dict(q=[0,4],p=[0,float(np.sqrt(.45))],eta=[-float(np.sqrt(.45)),float(np.sqrt(.45))]),
            rasterized_layers='filled contour maps and colorbar fills at 600 dpi; all surfaces, axes, text and trajectories remain vector'),
        base_errors=base,high_precision_comparison_points=points,
        high_precision_comparison_q_p_eta=list(maxima),dynamics=dict(profile_t_over_T=[0,.25,.5],
        alpha_cases=alpha_cases,beta_cases=beta_cases),two_soliton=two_sol,spectral_geometry=geometry_record,
        files=['analytic-regions.pdf','discrete-spectrum.pdf','one-soliton-fields.pdf','one-soliton-dynamics.pdf','one-soliton-velocity.pdf','two-soliton-fields.pdf','two-soliton-interaction.pdf'],figure_qa=qa)
    (ROOT/'calculations'/'plot-results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',action='store_true',help='also run installed nature-figure source, layout and PDF audits')
    parser.add_argument('--only-numerics',action='store_true',help='preserve spectral PDFs and regenerate only the numerical figures')
    parser.add_argument('--only-two-soliton',action='store_true',help='preserve existing figures and generate only the two-soliton figures')
    args=parser.parse_args()
    main(args.audit, args.only_numerics, args.only_two_soliton)
