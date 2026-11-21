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


def verify_nondegeneracy_identities(N, q0, x, t):
    """Verify the Cauchy-Hermitian congruence and direct positivity identity in Lemma 6.1."""
    zeta, constants = orbit_data(N, q0)
    m = 2*N
    reps = zeta[:N]
    def F(k):
        return mp.mpf('.4')*k**2-k-mp.mpf('.2')*q0**2+mp.mpf('.6')/(2*(k-2))
    sol = ResidueSolution(q0, zeta, constants, x, t, F)
    C_mat = mp.matrix([[1/(zeta[j]-mp.conj(zeta[k])) for k in range(m)] for j in range(m)])
    D = mp.diag(sol.C)
    A = C_mat*D.apply(mp.conj)
    B = C_mat.apply(mp.conj)*D
    L_N = mp.eye(m)+A*B
    K_blk = mp.eye(2*m)
    for j in range(m):
        for k in range(m):
            K_blk[j, m+k] = A[j, k]
            K_blk[m+j, k] = -B[j, k]
    g_diag = [mp.mpf(1)]*N+[q0/z for z in reps]
    Gamma = mp.diag(g_diag)
    S = 1j*Gamma.H*C_mat*Gamma
    K = mp.matrix([[1j/(reps[j]-mp.conj(reps[k])) for k in range(N)] for j in range(N)])
    T = mp.matrix([[1j*q0/(reps[j]*reps[k]+q0**2) for k in range(N)] for j in range(N)])
    K_conj = K.apply(mp.conj)
    K0 = (K+T*K_conj**-1*T.H)**-1
    T0 = K0*T*K_conj**-1
    S_inv = mp.matrix(m, m)
    for j in range(N):
        for k in range(N):
            S_inv[j, k] = K0[j, k]
            S_inv[j, N+k] = T0[j, k]
            S_inv[N+j, k] = T0.H[j, k]
            S_inv[N+j, N+k] = -mp.conj(K0[j, k])
    D0 = mp.diag(sol.C[:N]+[-mp.conj(v) for v in sol.C[:N]])
    E = mp.diag([mp.conj(v) for v in g_diag]+g_diag)
    EKE = E*K_blk*E**-1
    S_conj = S.apply(mp.conj)
    H = mp.matrix(2*m, 2*m)
    S_conj_inv = S_conj**-1
    for j in range(m):
        for k in range(m):
            H[j, k] = S_inv[j, k]
            H[m+j, m+k] = -S_conj_inv[j, k]
    H = H*EKE
    perm = list(range(0, N))+list(range(3*N, 4*N))+list(range(N, 2*N))+list(range(2*N, 3*N))
    Ht = mp.matrix([[H[r, c] for c in perm] for r in perm])
    Pblk, Xblk, Qblk = Ht[:m, :m], Ht[:m, m:], -Ht[m:, m:]
    det_L = mp.det(L_N)
    det_direct = abs(mp.det(S))**2*mp.det(Pblk)*mp.det(Qblk+Xblk.H*Pblk**-1*Xblk)
    errs = {
        'S_inverse_identity': norm(S*S_inv-mp.eye(m)),
        'H_Hermitian': norm(H-H.H),
        'det_L_positive_identity_rel': abs(det_L-det_direct)/abs(det_L),
        'det_L_imag': abs(mp.im(det_L)),
    }
    assert mp.re(det_L) > 0, det_L
    for name, err in errs.items():
        assert err < mp.mpf('1e-65'), (N, name, err)
    return errs


def verify_two_soliton_hmb():
    """70-digit independent verification of the Section 7.2 / 8.2 two-soliton solution."""
    mp.mp.dps = 70
    q0, a, eta0, alpha, beta = [mp.mpf(s) for s in ('1', '2', '.6', '1', '.1')]
    F = lambda k: 4*beta*k**2-alpha*k-2*beta*q0**2+eta0/(a*(k-a))
    geom = lambda z: ((z-q0**2/z)/2, (z+q0**2/z)/2)
    z1 = 2*q0*mp.exp(1j*mp.pi/3)
    z2 = 3*q0*mp.exp(2j*mp.pi/3)
    zeta = [z1, z2, -q0**2/mp.conj(z1), -q0**2/mp.conj(z2)]

    def s_prime(idx):
        z = zeta[idx]
        return (mp.fprod(z-zj for j, zj in enumerate(zeta) if j != idx)
                / mp.fprod(z-mp.conj(zj) for zj in zeta))

    c1, c2 = 1/s_prime(0), 1/s_prime(1)
    c3 = -q0**2/mp.conj(z1)**2*mp.conj(c1)
    c4 = -q0**2/mp.conj(z2)**2*mp.conj(c2)
    consts = [c1, c2, c3, c4]
    Omega = -alpha*q0**2+2*eta0/a

    def solve2(x, t):
        return ResidueSolution(q0, zeta, consts, x, t, F)

    def fields2(x, t):
        sol = solve2(x, t)
        R = sol.medium(a, eta0)
        phase = mp.exp(-1j*Omega*t)
        return phase*1j*sol.M1()[0, 1], -phase*R[0, 1], R[0, 0]

    points = [('-1.2', '-.35'), ('0', '0'), ('.37', '.21'), ('1.4', '-.1'), ('2.3', '.6')]
    sig = mp.diag([1, -1])
    sx = mp.matrix([[0, 1], [1, 0]])
    K_mat = mp.matrix([[0, 1], [-1, 0]])
    maxima = {}
    min_det = mp.mpf('1e99')
    for xs, ts in points:
        x, t = mp.mpf(xs), mp.mpf(ts)
        sol = solve2(x, t)
        min_det = min(min_det, mp.re(sol.schur_det))
        q, p, eta = fields2(x, t)
        R = sol.medium(a, eta0)
        qx = mp.diff(lambda xx: fields2(xx, t)[0], x)
        qxx = mp.diff(lambda xx: fields2(xx, t)[0], x, 2)
        qxxx = mp.diff(lambda xx: fields2(xx, t)[0], x, 3)
        qt = mp.diff(lambda tt: fields2(x, tt)[0], t)
        px = mp.diff(lambda xx: fields2(xx, t)[1], x)
        etax = mp.diff(lambda xx: fields2(xx, t)[2], x)
        u = mp.exp(1j*Omega*t)*q
        ux = mp.exp(1j*Omega*t)*qx
        uxx = mp.exp(1j*Omega*t)*qxx
        errs = {
            'linear_system': sol.linear_error,
            'same_phase': abs(mp.fprod(z/mp.conj(z) for z in zeta)-1),
            'PDE_q': abs(qt-1j*alpha*(qxx/2+abs(q)**2*q)-beta*(qxxx+6*abs(q)**2*qx)-2*p),
            'PDE_p': abs(px-2*eta*q+2j*a*p),
            'PDE_eta': abs(etax+mp.conj(q)*p+q*mp.conj(p)),
            'medium_two_sheets': norm(R-sol.medium(a, eta0, -1)),
            'medium_Hermitian': norm(R-R.H),
            'medium_trace': abs(R[0, 0]+R[1, 1]),
            'eta_real': abs(mp.im(eta)),
            'Bloch': abs(abs(p)**2+eta**2-eta0**2*(1+q0**2/a**2)),
            'residue_plus': max(norm(sol.C[j]*(sol.v[j]-sol.regular_column(zeta[j], 0))) for j in range(4)),
            'residue_minus': max(norm(mp.conj(sol.C[j])*(sol.u[j]-sol.regular_column(mp.conj(zeta[j]), 1))) for j in range(4)),
        }
        for idx, z in enumerate([mp.mpc('2.1', '.7'), mp.mpc('.3', '-.2'), mp.mpf('1.2')]):
            M = sol.M(z)
            k, lam = geom(z)
            errs[f'determinant_{idx}'] = abs(mp.det(M)-(1+q0**2/z**2))
            errs[f'sheet_{idx}'] = norm(sol.M(-q0**2/z)-M*(1j*z/q0)*sx)
            errs[f'conjugation_{idx}'] = norm(M-K_mat*sol.M(mp.conj(z)).apply(mp.conj)*K_mat**-1)
            Qmat = mp.matrix([[0, u], [-mp.conj(u), 0]])
            U = -1j*k*sig+Qmat
            d = 1j*alpha*abs(u)**2/2-beta*(u*mp.conj(ux)-ux*mp.conj(u))
            off = 2*beta*abs(u)**2*u+1j*alpha*ux/2+beta*uxx
            V1 = mp.matrix([[-2j*beta*abs(u)**2, alpha*u-2j*beta*ux],
                            [-alpha*mp.conj(u)-2j*beta*mp.conj(ux), 2j*beta*abs(u)**2]])
            V0 = mp.matrix([[d, off], [-mp.conj(off), -d]])
            V = 4j*beta*k**3*sig+k**2*(-1j*alpha*sig-4*beta*Qmat)+k*V1+V0+1j*Omega*sig/2+1j*R/(k-a)
            Mx = mp.matrix([[mp.diff(lambda xx: solve2(xx, t).M(z)[r, c], x) for c in range(2)] for r in range(2)])
            Mt = mp.matrix([[mp.diff(lambda tt: solve2(x, tt).M(z)[r, c], t) for c in range(2)] for r in range(2)])
            errs[f'Lax_x_{idx}'] = norm(Mx-U*M-1j*lam*M*sig)
            errs[f'Lax_t_{idx}'] = norm(Mt-V*M+1j*lam*F(k)*M*sig)
        for name, val in errs.items():
            maxima[name] = max(maxima.get(name, mp.mpf(0)), val)
    for name, val in maxima.items():
        assert val < mp.mpf('1e-55'), (name, val)

    # Verify asymptotic velocities, breathing periods, and elastic collision shifts
    k1, lam1 = geom(z1)
    k2, lam2 = geom(z2)
    v1 = mp.im(lam1*F(k1))/mp.im(lam1)
    v2 = mp.im(lam2*F(k2))/mp.im(lam2)
    w1 = -2*abs(lam1)**2*mp.im(F(k1))/mp.im(lam1)
    w2 = -2*abs(lam2)**2*mp.im(F(k2))/mp.im(lam2)
    s1_of = lambda z: ((z-zeta[0])/(z-mp.conj(zeta[0])))*((z-zeta[2])/(z-mp.conj(zeta[2])))
    s2_of = lambda z: ((z-zeta[1])/(z-mp.conj(zeta[1])))*((z-zeta[3])/(z-mp.conj(zeta[3])))
    dx1 = mp.log(abs(s2_of(z1)))/mp.im(lam1)
    dx2 = -mp.log(abs(s1_of(z2)))/mp.im(lam2)
    # Check dressed single-orbit agreement at t = -10 (orbit 1) and t = +10 (orbit 2)
    c1_eff = c1*s2_of(z1)**2/s2_of(0)
    c3_eff = -q0**2/mp.conj(z1)**2*mp.conj(c1_eff)
    q_full_1m = 1j*solve2(v1*(-10)+mp.mpf('.3'), mp.mpf(-10)).M1()[0, 1]
    q_eff_1m = s2_of(0)*1j*ResidueSolution(q0, [zeta[0], zeta[2]], [c1_eff, c3_eff], v1*(-10)+mp.mpf('.3'), mp.mpf(-10), F).M1()[0, 1]
    c2_eff = c2*s1_of(z2)**2/s1_of(0)
    c4_eff = -q0**2/mp.conj(z2)**2*mp.conj(c2_eff)
    q_full_2p = 1j*solve2(v2*10+mp.mpf('.3'), mp.mpf(10)).M1()[0, 1]
    q_eff_2p = s1_of(0)*1j*ResidueSolution(q0, [zeta[1], zeta[3]], [c2_eff, c4_eff], v2*10+mp.mpf('.3'), mp.mpf(10), F).M1()[0, 1]
    assert abs(q_full_1m-q_eff_1m) < mp.mpf('1e-30')
    assert abs(q_full_2p-q_eff_2p) < mp.mpf('1e-16')

    report = {
        'precision_digits': 70,
        'parameters': {'q0': '1', 'a': '2', 'eta0': '0.6', 'alpha': '1', 'beta': '0.1',
                       'zeta1': '2*exp(i*pi/3)', 'zeta2': '3*exp(i*2*pi/3)', 'b1_0': '1', 'b2_0': '1'},
        'points': points,
        'min_det_L2': mp.nstr(min_det, 16),
        'kinematics': {
            'Im_lam1': mp.nstr(mp.im(lam1), 12), 'Im_lam2': mp.nstr(mp.im(lam2), 12),
            'v1': mp.nstr(v1, 12), 'v2': mp.nstr(v2, 12),
            'T_t1': mp.nstr(2*mp.pi/abs(w1), 12), 'T_t2': mp.nstr(2*mp.pi/abs(w2), 12),
            'abs_s2_z1': mp.nstr(abs(s2_of(z1)), 12),
            'Delta_x1': mp.nstr(dx1, 12), 'Delta_x2': mp.nstr(dx2, 12),
            'dressed_orbit1_t_minus10_err': mp.nstr(abs(q_full_1m-q_eff_1m), 10),
            'dressed_orbit2_t_plus10_err': mp.nstr(abs(q_full_2p-q_eff_2p), 10),
        },
        'max_errors': {k: mp.nstr(v, 14) for k, v in maxima.items()},
    }
    target = Path(__file__).resolve().parents[1]/'calculations'/'two-soliton-results.json'
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


def main():
    mp.mp.dps = 80
    points = [('-0.4', '0.1'), ('0', '0'), ('0.6', '-0.2')]
    maxima, determinants, lemma_maxima = {}, {}, {}
    for N in (2, 3, 4):
        for q0 in ('1', '1.4'):
            for x, t in points:
                errors, determinant = check(N, mp.mpf(q0), mp.mpf(x), mp.mpf(t))
                for name, error in errors.items():
                    maxima[name] = max(maxima.get(name, mp.mpf(0)), error)
                determinants[f'N={N},q0={q0},x={x},t={t}'] = mp.nstr(determinant, 16)
        l_errs = verify_nondegeneracy_identities(N, mp.mpf('1.2'), mp.mpf('0.37'), mp.mpf('-0.21'))
        for name, error in l_errs.items():
            lemma_maxima[name] = max(lemma_maxima.get(name, mp.mpf(0)), error)
    report = {
        'precision_digits': 80,
        'orbit_counts': [2, 3, 4],
        'q0_values': ['1', '1.4'],
        'points_x_t': points,
        'method': 'Independent full block LU versus Schur solve, bordered determinant, and Lemma 6.1 congruence identities',
        'representatives_divided_by_q0': {str(N): [[mp.nstr(mp.re(z), 8), mp.nstr(mp.im(z), 8)]
                                                 for z in orbit_data(N, mp.mpf(1))[0][:N]] for N in (2, 3, 4)},
        'max_errors': {name: mp.nstr(error, 16) for name, error in maxima.items()},
        'lemma_6_1_max_errors': {name: mp.nstr(error, 16) for name, error in lemma_maxima.items()},
        'absolute_Schur_determinants': determinants,
        'scope': 'Numerical cross-checks of general-N identities; analytic proof of global det(L_N)>0 is given in Section 6.2.',
    }
    target = Path(__file__).resolve().parents[1]/'calculations'/'n-orbit-results.json'
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    two_report = verify_two_soliton_hmb()
    print(json.dumps({'n_orbit': report['lemma_6_1_max_errors'], 'two_soliton': two_report}, indent=2))


if __name__ == '__main__':
    main()
