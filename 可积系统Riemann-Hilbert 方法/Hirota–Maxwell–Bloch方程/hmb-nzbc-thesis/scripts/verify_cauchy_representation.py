"""Check the HMB jump convention and oriented Cauchy representation.

The nonzero density example is an independently prescribed scalar additive
RHP, not an existence test for the nonlinear reflected HMB inverse problem.
The contour has D+ on its left: outer real rays rightward, inner segment
leftward, upper semicircle clockwise, lower semicircle counterclockwise.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp


def symbolic_jump():
    rho, rt, e = sp.symbols('rho rt e', nonzero=True)
    J = sp.Matrix([[1, rt*e], [-rho/e, 1-rho*rt]])
    G = sp.Matrix([[rho*rt, rt*e], [-rho/e, 0]])
    assert sp.simplify(J.det()) == 1
    assert (sp.eye(2)-J.inv()-G).applyfunc(sp.simplify) == sp.zeros(2)
    assert (J*(sp.eye(2)-G)-sp.eye(2)).applyfunc(sp.simplify) == sp.zeros(2)
    return 'det J=1; G=I-J^-1; J(I-G)=I: exact'


def integral(density):
    real = mp.quad(density, [-mp.inf, -1])
    real -= mp.quad(density, [-1, 0, 1])
    real += mp.quad(density, [1, mp.inf])
    circle = lambda v: density(mp.exp(1j*v))*1j*mp.exp(1j*v)
    return real-mp.quad(circle, [0, mp.pi/2, mp.pi])+mp.quad(circle, [mp.pi, 3*mp.pi/2, 2*mp.pi])


def main():
    mp.mp.dps = 60
    a, b = 2j, -2j
    fp = lambda z: 1/(z-b)**2
    fm = lambda z: 1/(z-a)**2
    density = lambda s: fp(s)-fm(s)
    points = [mp.mpc('2', '.7'), mp.mpc('.3', '-.2'),
              mp.mpc('2', '-.7'), mp.mpc('.3', '.2')]
    errors = []
    for z in points:
        plus = mp.im(z)*(1-1/abs(z)**2) > 0
        value = integral(lambda s: density(s)/(s-z))/(2j*mp.pi)
        errors.append(abs(value-(fp(z) if plus else fm(z))))
    # A second density has a nonzero 1/z coefficient, checking its minus sign.
    w1 = lambda s: 1/(s-b)-1/(s-a)
    coefficient = -integral(w1)/(2j*mp.pi)
    errors.append(abs(coefficient-1))
    assert max(errors) < mp.mpf('1e-50'), errors
    report = {
        'precision_digits': 60,
        'symbolic_jump': symbolic_jump(),
        'contour': 'real rays right; inner segment left; upper circle clockwise; lower circle counterclockwise',
        'test': 'independent additive RHP with W=F_plus-F_minus, F_plus=(z+2i)^-2, F_minus=(z-2i)^-2',
        'points': [str(z) for z in points],
        'cauchy_value_errors': [mp.nstr(e, 16) for e in errors[:-1]],
        'infinity_coefficient_sign_error': mp.nstr(errors[-1], 16),
        'interpretation': 'checks signs and orientation, not reflected HMB inverse solvability',
    }
    out = Path(__file__).resolve().parents[1]/'calculations'/'cauchy-representation-results.json'
    out.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
