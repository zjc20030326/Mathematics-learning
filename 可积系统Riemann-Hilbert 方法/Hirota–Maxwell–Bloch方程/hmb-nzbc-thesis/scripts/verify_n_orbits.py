"""Check the N-orbit Schur and bordered-determinant formulas independently.

The full 4N by 4N LU solve in verify_direct_rhp is compared with the new
2N by 2N Schur solve. Samples test identities, not global invertibility.
"""
import json
from pathlib import Path

import mpmath as mp

from verify_direct_rhp import ResidueSolution, norm


def orbit_data(N, q0):
    examples = {
        2: [('2', '1'), ('-3', '1.5')],
        3: [('2', '1'), ('1', '2'), ('0', '3')],
        4: [('2', '1'), ('-3', '1.5'), ('1', '3'), ('-1.2', '3.6')],
    }
    representatives = [q0*mp.mpc(x, y) for x, y in examples[N]]
    constants = [mp.mpc(mp.mpf('.2')*(j+1), mp.mpf('.3')+mp.mpf('.1')*j)
                 for j in range(N)]
    zeta = representatives + [-q0**2/mp.conj(z) for z in representatives]
    constants += [-q0**2*mp.conj(c)/mp.conj(z)**2
                  for z, c in zip(representatives, constants)]
    return zeta, constants


def check(N, q0, x, t):
    zeta, constants = orbit_data(N, q0)
    m = 2*N
    def F(k):
        return mp.mpf('.4')*k**2-k-mp.mpf('.2')*q0**2+mp.mpf('.6')/(2*(k-2))
    sol = ResidueSolution(q0, zeta, constants, x, t, F)
    A, B = mp.matrix(m), mp.matrix(m)
    f, g = mp.matrix(m, 2), mp.matrix(m, 2)
    K = mp.eye(2*m)
    for j in range(m):
        f[j, 0], f[j, 1] = 1, -1j*q0/zeta[j]
        g[j, 0], g[j, 1] = -1j*q0/mp.conj(zeta[j]), 1
        for k in range(m):
            A[j, k] = mp.conj(sol.C[k])/(zeta[j]-mp.conj(zeta[k]))
            B[j, k] = sol.C[k]/(mp.conj(zeta[j])-zeta[k])
            K[j, m+k], K[m+j, k] = A[j, k], -B[j, k]
    L = mp.eye(m)+A*B
    rhs = f-A*g
    v = mp.matrix(m, 2)
    for col in range(2):
        v[:, col] = mp.lu_solve(L, rhs[:, col])
    u = g+B*v
    errors = {
        'full_system': sol.linear_error,
        'Schur_v': max(norm(v[j, :].T-sol.v[j]) for j in range(m)),
        'Schur_u': max(norm(u[j, :].T-sol.u[j]) for j in range(m)),
        'block_determinant_relative': abs(mp.det(K)/mp.det(L)-1),
        'same_phase': abs(mp.fprod(z/mp.conj(z) for z in zeta)-1),
        'constant_pair': max(abs(constants[N+j]+q0**2*mp.conj(constants[j])/mp.conj(zeta[j])**2)
                             for j in range(N)),
        'evolving_constant_pair': max(abs(sol.C[N+j]+q0**2*mp.conj(sol.C[j])/mp.conj(zeta[j])**2)
                                      for j in range(N)),
    }
    C = mp.matrix(sol.C)
    r = mp.matrix([1+1j*q0*mp.fsum(mp.conj(sol.C[k])/((zeta[j]-mp.conj(zeta[k]))*mp.conj(zeta[k]))
                              for k in range(m)) for j in range(m)])
    errors['rhs_first_component'] = norm(r-rhs[:, 0])
    border = mp.matrix(m+1)
    for j in range(m):
        border[j, m] = r[j]
        border[m, j] = -1j*C[j]
        for k in range(m):
            border[j, k] = L[j, k]
    border[m, m] = q0
    q = q0+1j*(C.T*v[:, 0])[0]
    errors['bordered_determinant_q'] = abs(mp.det(border)/mp.det(L)-q)
    errors['infinity_reconstruction_q'] = abs(q-1j*sol.M1()[0, 1])
    errors['residue_plus'] = max(norm(sol.v[j]-sol.regular_column(zeta[j], 0)) for j in range(m))
    errors['residue_minus'] = max(norm(sol.u[j]-sol.regular_column(mp.conj(zeta[j]), 1)) for j in range(m))
    sig1 = mp.matrix([[0, 1], [1, 0]])
    sig2 = mp.matrix([[0, -1j], [1j, 0]])
    for name in ('determinant', 'sheet', 'conjugation', 'branch_kernel'):
        errors[name] = mp.mpf(0)
    for z in [q0*w for w in [mp.mpc('2.1', '.7'), mp.mpc('.3', '-.2'), mp.mpf('1.2'), 1j, -1j]]:
        M = sol.M(z)
        errors['determinant'] = max(errors['determinant'], abs(mp.det(M)-(1+q0**2/z**2)))
        errors['sheet'] = max(errors['sheet'], norm(sol.M(-q0**2/z)-M*(1j*z/q0)*sig1))
        errors['conjugation'] = max(errors['conjugation'], norm(M-sig2*sol.M(mp.conj(z)).apply(mp.conj)*sig2))
        if z in (1j*q0, -1j*q0):
            errors['branch_kernel'] = max(errors['branch_kernel'], norm(M-M*(1j*z/q0)*sig1))
    R = sol.medium(mp.mpf(2), mp.mpf('.6'))
    errors['medium_two_sheets'] = norm(R-sol.medium(mp.mpf(2), mp.mpf('.6'), -1))
    errors['medium_Hermitian'] = norm(R-R.H)
    for name, error in errors.items():
        assert error < mp.mpf('1e-65'), (N, q0, x, t, name, mp.nstr(error, 10))
    return errors, abs(mp.det(L))


def main():
    mp.mp.dps = 80
    points = [('-0.4', '0.1'), ('0', '0'), ('0.6', '-0.2')]
    maxima, determinants = {}, {}
    for N in (2, 3, 4):
        for q0 in ('1', '1.4'):
            for x, t in points:
                errors, determinant = check(N, mp.mpf(q0), mp.mpf(x), mp.mpf(t))
                for name, error in errors.items():
                    maxima[name] = max(maxima.get(name, mp.mpf(0)), error)
                determinants[f'N={N},q0={q0},x={x},t={t}'] = mp.nstr(determinant, 16)
    report = {
        'precision_digits': 80,
        'orbit_counts': [2, 3, 4],
        'q0_values': ['1', '1.4'],
        'points_x_t': points,
        'method': 'Independent full block LU versus Schur solve and bordered determinant',
        'representatives_divided_by_q0': {str(N): [[mp.nstr(mp.re(z), 8), mp.nstr(mp.im(z), 8)]
                                                 for z in orbit_data(N, mp.mpf(1))[0][:N]] for N in (2, 3, 4)},
        'max_errors': {name: mp.nstr(error, 16) for name, error in maxima.items()},
        'absolute_Schur_determinants': determinants,
        'scope': 'Sampled identities only. No claim of arbitrary-N global invertibility or PDE verification.',
    }
    target = Path(__file__).resolve().parents[1]/'calculations'/'n-orbit-results.json'
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
