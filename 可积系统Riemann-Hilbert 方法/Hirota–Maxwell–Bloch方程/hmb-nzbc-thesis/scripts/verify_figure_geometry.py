"""Check the actual symbolic HMB figure geometry against the spectral domain.

This checks Matplotlib objects from plot_fields.py without exporting figures.
No numerical HMB propagation or inverse-scattering solvability is asserted.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib.colors as colors
import plot_fields as plots

ROOT = Path(__file__).resolve().parents[1]


def domain_value(z):
    return float(np.imag(z)*(1-1/abs(z)**2))


def main():
    captured = {}
    plots.export = lambda fig, name, **kwargs: captured.update({name: fig})
    plots.spectral_figures()
    fig = captured['analytic-regions']
    ax = fig.axes[0]
    fig.canvas.draw()
    arrows = []
    for annotation in ax.texts:
        if not hasattr(annotation, 'arrow_patch') or annotation.arrow_patch is None:
            continue
        head = complex(*annotation.xy)
        tail = complex(*annotation.xyann)
        if abs(head.imag) < 1e-12 and abs(tail.imag) < 1e-12:
            mid = (head+tail)/2
            kind = 'real axis'
        elif abs(abs(head)-1) < 1e-12 and abs(abs(tail)-1) < 1e-12:
            mid = (head+tail)/abs(head+tail)
            kind = 'upper circle' if mid.imag > 0 else 'lower circle'
        else:
            continue  # imaginary coordinate guide, not a HMB contour arc
        left_normal = 1j*(head-tail)/abs(head-tail)
        left = domain_value(mid+.02*left_normal)
        right = domain_value(mid-.02*left_normal)
        assert left > 0 and right < 0, (kind, tail, head, left, right)
        arrows.append(dict(arc=kind, tail=[tail.real,tail.imag], head=[head.real,head.imag],
                           left_domain_value=left, right_domain_value=right))
    assert len(arrows) == 6
    blue = colors.to_rgba('#E1EDF5')
    samples = [-1.6+1.2j, 1.7+.9j, .2-.4j, -.4-.2j,
               -1.6-1.2j, 1.7-.9j, .2+.4j, -.4+.2j]
    shade_checks = []
    for z in samples:
        screen = ax.transData.transform((z.real,z.imag))
        color = colors.to_rgba('white')
        for patch in sorted(ax.patches,key=lambda p:p.get_zorder()):
            if patch.contains_point(screen):
                color = tuple(patch.get_facecolor())
        is_blue = np.allclose(color, blue)
        assert is_blue == (domain_value(z) > 0), (z,color)
        shade_checks.append(dict(point=[z.real,z.imag],domain='D+' if is_blue else 'D-'))
    ax = captured['discrete-spectrum'].axes[0]
    quartet = []
    for line in ax.lines:
        if line.get_marker() != 'o':
            continue
        points = np.asarray(line.get_xdata())+1j*np.asarray(line.get_ydata())
        plus = line.get_markerfacecolor() == '#2B638A'
        assert all(domain_value(z)*(1 if plus else -1)>0 for z in points)
        quartet.extend(points.tolist())
    assert len(quartet) == 4 and len(set(quartet)) == 4
    for z in quartet:
        assert min(abs(np.asarray(quartet)-np.conj(z))) < 1e-12
        assert min(abs(np.asarray(quartet)+1/z)) < 1e-12
    result = dict(verdict='PASS',q0='schematic unit radius only',
                  contour_arrows=arrows,domain_shading_samples=shade_checks,
                  quartet='four distinct off-contour points closed under conjugation and h(z)=-q0^2/z',
                  plus_points='zeta_j and h(conj(zeta_j))',minus_points='conj(zeta_j) and h(zeta_j)')
    (ROOT/'calculations'/'figure-geometry-check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    plots.plt.close('all')


if __name__ == '__main__':
    main()
