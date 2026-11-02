"""Solve the residue system (finite1/finite2) and verify reconstruction.

The solver accepts spectral data, builds a 2m by 2m scalar matrix and solves
two right-hand sides. No projector or preconstructed field enters the solver.
Run from any directory; mpmath is the only third-party dependency.
"""
from pathlib import Path
import argparse
import json
import mpmath as mp


def parameters():
    # Recreate decimals at the current precision (also inside mp.diff).
    return tuple(mp.mpf(s) for s in ('1', '2', '.6', '1', '.1', '1.5'))


def dispersion(lam):
    Q, a, eta0, alpha, beta, _ = parameters()
    return 4*beta*lam**2-alpha*lam-2*beta*Q**2+eta0/(a*(lam-a))


def geometry(w, Q):
    return (w-Q**2/w)/2, (w+Q**2/w)/2


def data():
    Q, _, _, _, _, b = parameters()
    k = mp.sqrt(b*b-Q*Q)
    zeta = [1j*(b+k), -1j*(b-k)]
    c1 = 4j*k/(1-Q**2/zeta[0]**2)
    return zeta, [c1, -Q**2/mp.conj(zeta[0])**2*mp.conj(c1)]


class ResidueSolution:
    def __init__(self, Q, zeta, constants, x, z, F):
        self.Q, self.zeta = Q, zeta
        n = len(zeta)
        self.C = []
        for w, c in zip(zeta, constants):
            lam, rho = geometry(w, Q)
            self.C.append(c*mp.exp(2j*rho*(x-F(lam)*z)))
        A = mp.matrix(n, n)
        B = mp.matrix(n, n)
        L = mp.eye(2*n)
        rhs = mp.matrix(2*n, 2)
        for j in range(n):
            rhs[j, 0], rhs[j, 1] = -1j*Q/zeta[j], 1
            rhs[n+j, 0], rhs[n+j, 1] = 1, -1j*Q/mp.conj(zeta[j])
            for k in range(n):
                A[j, k] = mp.conj(self.C[k])/(zeta[j]-mp.conj(zeta[k]))
                B[j, k] = self.C[k]/(mp.conj(zeta[j])-zeta[k])
                L[j, n+k], L[n+j, k] = A[j, k], -B[j, k]
        answer = mp.matrix(2*n, 2)
        for col in range(2):
            answer[:, col] = mp.lu_solve(L, rhs[:, col])
        self.v = [answer[j, :].T for j in range(n)]
        self.u = [answer[n+j, :].T for j in range(n)]
        self.linear_error = norm(L*answer-rhs)
        self.schur_det = mp.det(mp.eye(n)+A*B)

    def M(self, w):
        R = mp.matrix([[1, -1j*self.Q/w], [-1j*self.Q/w, 1]])
        for j, zz in enumerate(self.zeta):
            R[:, 0] += self.C[j]*self.v[j]/(w-zz)
            R[:, 1] -= mp.conj(self.C[j])*self.u[j]/(w-mp.conj(zz))
        return R

    def regular_column(self, w, col):
        # Evaluate only the regular column, including at the other column's pole.
        if col == 1:
            r = mp.matrix([-1j*self.Q/w, 1])
            for j, zz in enumerate(self.zeta):
                r -= mp.conj(self.C[j])*self.u[j]/(w-mp.conj(zz))
        else:
            r = mp.matrix([1, -1j*self.Q/w])
            for j, zz in enumerate(self.zeta):
                r += self.C[j]*self.v[j]/(w-zz)
        return r

    def M1(self):
        R = mp.matrix([[0, -1j*self.Q], [-1j*self.Q, 0]])
        for j in range(len(self.zeta)):
            R[:, 0] += self.C[j]*self.v[j]
            R[:, 1] -= mp.conj(self.C[j])*self.u[j]
        return R

    def medium(self, a, eta0, sheet=1):
        w = a+sheet*mp.sqrt(a*a+self.Q*self.Q)
        _, rho = geometry(w, self.Q)
        M = self.M(w)
        return M*mp.diag([eta0*rho/a, -eta0*rho/a])*M**-1


def solve(x, z):
    Q, *_ = parameters()
    zz, cc = data()
    return ResidueSolution(Q, zz, cc, x, z, dispersion)


def fields(x, z):
    Q, a, eta0, alpha, _, _ = parameters()
    sol = solve(x, z)
    R = sol.medium(a, eta0)
    phase = mp.exp(-1j*(-alpha*Q*Q+2*eta0/a)*z)
    # Do not discard the imaginary part of eta; test it independently.
    return phase*1j*sol.M1()[0, 1], -phase*R[0, 1], R[0, 0]


def norm(M):
    return max(abs(v) for v in M)


def derivative_matrix(fun, value):
    return mp.matrix([[mp.diff(lambda t: fun(t)[j, k], value)
                       for k in range(2)] for j in range(2)])


def checks_at(x, z):
    Q, a, eta0, alpha, beta, b = parameters()
    sig = mp.diag([1, -1])
    sx = mp.matrix([[0, 1], [1, 0]])
    K = mp.matrix([[0, 1], [-1, 0]])
    sol = solve(x, z)
    zz, cc = data()
    errors = {'linear_system': sol.linear_error}
    errors['theta_condition'] = abs(mp.fprod(w/mp.conj(w) for w in zz)-1)
    errors['constant_pair'] = abs(cc[1]+Q*Q/mp.conj(zz[0])**2*mp.conj(cc[0]))
    for j, w in enumerate(zz):
        errors[f'residue_plus_{j}'] = norm(sol.C[j]*(sol.v[j]-sol.regular_column(w, 1)))
        errors[f'residue_minus_{j}'] = norm(mp.conj(sol.C[j])*(sol.u[j]-sol.regular_column(mp.conj(w), 0)))
    h = Q/(b+mp.sqrt(b*b-Q*Q))
    t = mp.exp(2*mp.sqrt(b*b-Q*Q)*(x-dispersion(1j*b)*z))
    dm = abs(t-h)**2+abs(1-h*t)**2
    dp = abs(t+h)**2+abs(1+h*t)**2
    errors['schur_factorization'] = abs(sol.schur_det-dm*dp/(abs(t)**4*(1+h*h)**2))
    q, p, eta = fields(x, z)
    Omega = -alpha*Q*Q+2*eta0/a
    u = mp.exp(1j*Omega*z)*q
    ux = mp.diff(lambda xx: mp.exp(1j*Omega*z)*fields(xx, z)[0], x)
    uxx = mp.diff(lambda xx: mp.exp(1j*Omega*z)*fields(xx, z)[0], x, 2)
    R = sol.medium(a, eta0)
    errors['medium_two_sheets'] = norm(R-sol.medium(a, eta0, -1))
    errors['medium_Hermitian'] = norm(R-R.H)
    errors['medium_trace'] = abs(mp.trace(R)) if hasattr(mp, 'trace') else abs(R[0, 0]+R[1, 1])
    errors['eta_real'] = abs(mp.im(eta))
    errors['Bloch'] = abs(abs(p)**2+eta**2-eta0**2*(1+Q*Q/(a*a)))
    for index, w in enumerate([mp.mpc('2.1', '.7'), mp.mpc('.3', '-.2'), mp.mpf('1.2')]):
        M = sol.M(w)
        lam, rho = geometry(w, Q)
        errors[f'determinant_{index}'] = abs(mp.det(M)-(1+Q*Q/w**2))
        errors[f'sheet_{index}'] = norm(sol.M(-Q*Q/w)-M*(1j*w/Q)*sx)
        errors[f'conjugation_{index}'] = norm(M-K*sol.M(mp.conj(w)).apply(mp.conj)*K**-1)
        A = mp.matrix([[0, u], [-mp.conj(u), 0]])
        U = -1j*lam*sig+A
        d = 1j*alpha*abs(u)**2/2-beta*(u*mp.conj(ux)-ux*mp.conj(u))
        off = 2*beta*abs(u)**2*u+1j*alpha*ux/2+beta*uxx
        V1 = mp.matrix([[-2j*beta*abs(u)**2, alpha*u-2j*beta*ux],
                        [-alpha*mp.conj(u)-2j*beta*mp.conj(ux), 2j*beta*abs(u)**2]])
        V0 = mp.matrix([[d, off], [-mp.conj(off), -d]])
        V = 4j*beta*lam**3*sig+lam**2*(-1j*alpha*sig-4*beta*A)+lam*V1+V0+1j*Omega*sig/2+1j*R/(lam-a)
        Mx = derivative_matrix(lambda xx: solve(xx, z).M(w), x)
        Mz = derivative_matrix(lambda zz0: solve(x, zz0).M(w), z)
        errors[f'Lax_x_{index}'] = norm(Mx-U*M-1j*rho*M*sig)
        errors[f'Lax_z_{index}'] = norm(Mz-V*M+1j*rho*dispersion(lam)*M*sig)
    for index, w in enumerate([Q, -Q, 1j*Q, -1j*Q]):
        M = sol.M(w)
        errors[f'special_point_det_{index}'] = abs(mp.det(M)-(1+Q*Q/w**2))
        if index >= 2:
            errors[f'branch_kernel_{index}'] = norm(M-M*(1j*w/Q)*sx)
    qx = mp.diff(lambda xx: fields(xx, z)[0], x)
    qxx = mp.diff(lambda xx: fields(xx, z)[0], x, 2)
    qxxx = mp.diff(lambda xx: fields(xx, z)[0], x, 3)
    qz = mp.diff(lambda zz0: fields(x, zz0)[0], z)
    px = mp.diff(lambda xx: fields(xx, z)[1], x)
    etax = mp.diff(lambda xx: fields(xx, z)[2], x)
    errors['PDE_q'] = abs(qz-1j*alpha*(qxx/2+abs(q)**2*q)-beta*(qxxx+6*abs(q)**2*qx)-2*p)
    errors['PDE_p'] = abs(px-2*eta*q+2j*a*p)
    errors['PDE_eta'] = abs(etax+mp.conj(q)*p+q*mp.conj(p))
    # Scaled errors distinguish the specified orders from a mere finite value.
    large = mp.mpf('1e15')
    errors['normalization_infty_scaled'] = norm(large*(sol.M(large)-mp.eye(2))-sol.M1())/large
    small = 1/large
    errors['normalization_zero_scaled'] = norm(small*sol.M(small)+1j*Q*sx)*small
    return errors


def main(compare_legacy=False):
    mp.mp.dps = 70
    points = [('-1.2', '-.35'), ('0', '0'), ('.37', '.21'), ('1.4', '-.1'), ('2.3', '.6')]
    results = [checks_at(mp.mpf(x), mp.mpf(z)) for x, z in points]
    maxima = {name: max(r[name] for r in results) for name in results[0]}
    for name, err in maxima.items():
        limit = mp.mpf('1e-25') if name.startswith('normalization_') else mp.mpf('1e-55')
        assert err < limit, (name, err)
    boundaries = []
    Q, a, eta0, alpha, _, _ = parameters()
    for x in [-30, 30]:
        z = mp.mpf('.21')
        values = fields(mp.mpf(x), z)
        qb = Q*mp.exp(-1j*(-alpha*Q*Q+2*eta0/a)*z)
        row = [abs(v-vb) for v, vb in zip(values, (qb, -1j*eta0*qb/a, eta0))]
        assert max(row) < mp.mpf('1e-27')
        boundaries.append([mp.nstr(v, 14) for v in row])
    report = {
        'precision_digits': 70, 'method': 'LU solve of residue block matrix; no Darboux input',
        'parameters': dict(zip(('Q', 'a', 'eta0', 'alpha', 'beta', 'b'), map(str, parameters()))),
        'points': points, 'max_errors': {k: mp.nstr(v, 14) for k, v in maxima.items()},
        'boundary_errors_x_minus_plus_30': boundaries,
        'normalization_note': 'scaled remainder divided by large w, or multiplied by small w; 1e15 and 1e-15',
    }
    if compare_legacy:
        # Optional historical comparison, after every independent RHP check.
        from verify_solitons import fields as darboux_fields
        mp.mp.dps = 70  # the legacy module sets its own default precision
        comparison = [mp.mpf(0)]*3
        for x, z in points:
            xx, zz = mp.mpf(x), mp.mpf(z)
            direct, other = fields(xx, zz), darboux_fields(xx, zz, 1)
            comparison = [max(old, abs(v-d)) for old, v, d in zip(comparison, direct, other)]
        assert max(comparison) < mp.mpf('1e-55')
        report['legacy_Darboux_comparison_q_p_eta'] = [mp.nstr(v, 14) for v in comparison]
    path = Path(__file__).resolve().parents[1]/'calculations'/'direct-rhp-results.json'
    path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare-legacy', action='store_true',
                        help='also compare with the retained historical Darboux implementation')
    main(parser.parse_args().compare_legacy)
