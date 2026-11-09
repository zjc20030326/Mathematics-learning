"""Re-run the owning thesis checks without writing to the source project."""
from pathlib import Path
import contextlib
import io
import json
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'hmb-nzbc-thesis'
sys.path.insert(0, str(SOURCE / 'scripts'))
REPORTS = ROOT / 'calculations'
REPORTS.mkdir(exist_ok=True)


def run(name, main=False):
    original = SOURCE / 'scripts' / name
    namespace = {'__file__': str(ROOT / 'checks' / name), '__name__': '__explanation_verification__'}
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        exec(compile(original.read_text(encoding='utf-8-sig'), str(original), 'exec'), namespace)
        if main:
            namespace['main']()
    print(f'{name}: all assertions passed')
    return namespace


def main():
    run('verify_symbolic.py')
    run('verify_direct_algebra.py')
    numerical = run('verify_direct_rhp.py', main=True)
    run('verify_n_orbits.py', main=True)
    mp = numerical['mp']
    mp.mp.dps = 70
    q, p, eta = numerical['fields'](mp.mpf(0), mp.mpf(0))
    assert abs(q + 4) < mp.mpf('1e-60')
    assert abs(p - mp.mpc(0, '.66')) < mp.mpf('1e-60')
    assert abs(eta + mp.mpf('.12')) < mp.mpf('1e-60')
    q0, a, eta0, alpha, beta, b = numerical['parameters']()
    kap = mp.sqrt(b*b-q0*q0)
    F = numerical['dispersion'](1j*b)
    report = {'origin_q': str(q), 'origin_p': str(p), 'origin_eta': str(eta),
              'Omega': str(-alpha*q0*q0+2*eta0/a), 'Bloch_squared': str(abs(p)**2+eta**2),
              'A_b': str(mp.re(F)), 'B_b': str(mp.im(F)),
              'T_t': str(mp.pi/(kap*abs(mp.im(F)))),
              'exp_decay_factor_x_30': str(mp.exp(-60*kap)),
              'source_scripts_unchanged': True}
    (REPORTS / 'expanded-numerics-check.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
