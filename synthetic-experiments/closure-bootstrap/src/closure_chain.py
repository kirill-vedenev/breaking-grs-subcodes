"""Key recovery for binary Goppa codes from first-order Frobenius closures.

Chain (all steps after `closures_from_secret` read only the public generator
matrix, the supplied closures and the public parameters m, t):

  A. closures E_i at |S| positions  ->  binary quadrics vanishing on the curve
  B. one position: gradient kernel, branch separation, one tangent vector
  C. deep local jet by continuation on the binary quadrics (v_{2j} = v_j^2)
  D. quadrics over K vanishing on the selected Frobenius conjugate branch
  E. aligned tangent planes at every public column
  F. minor code, Sidelnikov-Shestakov, embedding of the public code

The closures are computed from the secret key, so the experiment is synthetic.
The secret is used again only in the audits, which do not feed the recovery.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
import sys
import time
import traceback
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path

import numpy as np
from sage.all import (GF, PolynomialRing, matrix, vector, set_random_seed,
                      shuffle, random_matrix, ZZ)

HERE = Path(__file__).resolve().parent
for _cand in (HERE, HERE.parents[1] / 'minor-code-fullness/src'):
    if (_cand / 'minor.py').exists():
        sys.path.insert(0, str(_cand))
        break
import minor as MC

F2 = GF(2)
CHUNK = 128


# ---------------------------------------------------------------- field tables

class Field:
    """GF(2^m) with elements encoded as integers (polynomial basis)."""

    def __init__(self, E):
        self.E = E
        self.m = int(E.degree())
        q = self.q = int(E.order())
        g = E.multiplicative_generator()
        pw = np.zeros(q - 1, dtype=np.uint16)
        x = E.one()
        for i in range(q - 1):
            pw[i] = int(x.to_integer())
            x *= g
        self.LZ = 2 * (q - 1)
        self.log = np.full(q, self.LZ, dtype=np.int32)
        self.log[pw] = np.arange(q - 1, dtype=np.int32)
        self.exp = np.zeros(4 * (q - 1) + 1, dtype=np.uint16)
        self.exp[:2 * (q - 1)] = np.tile(pw, 2)
        self.pw = pw
        self.elts = [E.from_integer(i) for i in range(q)]

    def mul(self, a, b):
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        assert np.all(a != 0)
        return self.pw[(-self.log[a]) % (self.q - 1)]

    def frob(self, a, s=1):
        """a -> a^(2^s), entrywise."""
        a = np.asarray(a)
        out = self.pw[(self.log[a].astype(np.int64) * (2 ** s)) % (self.q - 1)]
        return np.where(a == 0, 0, out).astype(np.uint16)

    def sqrt(self, a):
        return self.frob(a, self.m - 1)

    def to_sage(self, arr):
        arr = np.asarray(arr)
        elts = self.elts
        return matrix(self.E, arr.shape[0], arr.shape[1],
                      [elts[x] for x in arr.ravel().tolist()])

    def from_sage(self, M):
        return np.array([[int(x.to_integer()) for x in row] for row in M.rows()],
                        dtype=np.uint16).reshape(M.nrows(), M.ncols())

    def bits(self, arr):
        """(..., m) array of the binary coordinates."""
        return ((np.asarray(arr)[..., None] >> np.arange(self.m)) & 1).astype(np.uint8)


def to_f2(a):
    return matrix(F2, np.ascontiguousarray(a, dtype=np.uint8))


def f2_array(M):
    return np.array(M, dtype=np.uint8).reshape(M.nrows(), M.ncols())


def binary_times_field(fld, Q, X):
    """Q (Sage GF(2), r x N) times X (N x c array of field elements)."""
    out = np.zeros((Q.nrows(), X.shape[1]), dtype=np.uint16)
    for bit in range(fld.m):
        plane = to_f2((X >> bit) & 1)
        out |= f2_array(Q * plane).astype(np.uint16) << bit
    return out


def rows_times_field(fld, rows, x):
    """Binary rows (r x N float32 array) times a field vector x (N,)."""
    planes = fld.bits(x).astype(np.float32)
    prod = (rows @ planes).astype(np.int64) & 1
    return (prod << np.arange(fld.m)).sum(axis=1).astype(np.uint16)


# ---------------------------------------------------------------- the instance

def instance(family, m, n, t, seed):
    """goppa: Gamma(alpha, g) with g irreducible of degree t, inside
    GRS_{n-2t}; alternant: binary alternant code of degree t with random
    multipliers, inside GRS_{n-t}. The secret holds the dual multiplier nu of
    the ambient GRS code."""
    set_random_seed(seed)
    E = GF(2 ** m, 'z')
    R = PolynomialRing(E, 'X')
    X = R.gen()
    if family == 'goppa':
        while True:
            g = X ** t + R.random_element(degree=t - 1)
            if g.is_irreducible():
                break
    alpha = list(E)
    shuffle(alpha)
    alpha = alpha[:n]
    fld = Field(E)
    al = np.array([int(a.to_integer()) for a in alpha], dtype=np.uint16)
    if family == 'goppa':
        gv = np.array([int(g(a).to_integer()) for a in alpha], dtype=np.uint16)
        row = fld.inv(gv)
        nu = fld.mul(row, row)
    else:
        rng = np.random.default_rng(seed + 7 * 10 ** 6)
        row = rng.integers(1, fld.q, size=n).astype(np.uint16)
        nu = row
    H = [row]
    for _ in range(t - 1):
        H.append(fld.mul(H[-1], al))
    H = np.array(H, dtype=np.uint16)
    Hb = fld.bits(H).transpose(0, 2, 1).reshape(t * m, n)
    P = to_f2(Hb).right_kernel_matrix()
    secret = dict(alpha=al, nu=nu)
    return E, fld, P, secret


def true_tangents(fld, Pn, secret):
    """Column j: sum_{l != j} nu_l y_l / (alpha_j - alpha_l).

    It lies in F'(alpha_j)/Pi'(alpha_j) + K F(alpha_j)."""
    al, nu = secret['alpha'], secret['nu']
    n = len(al)
    diff = al[:, None] ^ al[None, :]                     # [l, j]
    np.fill_diagonal(diff, 1)
    W = fld.mul(nu[:, None], fld.inv(diff))              # nu_l / (alpha_j - alpha_l)
    np.fill_diagonal(W, 0)
    Yhat = np.zeros(Pn.shape, dtype=np.uint16)
    for a in range(Pn.shape[0]):
        Yhat[a] = np.bitwise_xor.reduce(W[Pn[a].astype(bool)], axis=0)
    return Yhat


def closures_from_secret(fld, Pn, Yhat, positions):
    """Frobenius closure of the tangent plane at each position, as an echelon
    basis over GF(2) of the (m+1)-dimensional space."""
    out = {}
    for i in positions:
        rows = np.vstack([Pn[:, i][None, :], fld.bits(Yhat[:, i]).T])
        out[i] = f2_array(to_f2(rows).row_space().basis_matrix())
    return out


# ---------------------------------------------------------------- the recovery

def polar_rows(p, U, pa, pb):
    """Rows B(p, u) over the pairs, for binary p and binary rows u of U."""
    return (p[pa][None, :] & U[:, pb]) ^ (p[pb][None, :] & U[:, pa])


def polar_field(fld, v, pa, pb, k):
    """Matrix (pairs x k) of the linear map x -> (B(v, x) per pair)."""
    G = np.zeros((len(pa), k), dtype=np.uint16)
    off = pa != pb
    idx = np.nonzero(off)[0]
    G[idx, pa[idx]] = v[pb[idx]]
    G[idx, pb[idx]] = v[pa[idx]]
    return G


def full_coeff(fld, logV, e, pa, pb):
    """[T^e] of V_a(T) V_b(T) for every pair a <= b (diagonal: V_a(T)^2)."""
    out = np.zeros(len(pa), dtype=np.uint16)
    for s in range(0, e + 1, CHUNK):
        i = np.arange(s, min(e + 1, s + CHUNK))
        out ^= np.bitwise_xor.reduce(fld.exp[logV[i][:, pa] + logV[e - i][:, pb]], axis=0)
    return out


def recover(E, fld, Pn, closures, m, D, goppa_t, seed, log=None):
    """Public inputs only: the generator matrix, the closures, m, the degree D
    of the curve and, for a binary Goppa code, the degree t (else 0)."""
    rng = np.random.default_rng(seed)
    k, n = Pn.shape
    t = goppa_t
    pairs = list(combinations_with_replacement(range(k), 2))
    pa, pb = np.array(pairs, dtype=np.int64).T
    npairs = len(pairs)
    rep = dict(k=k, n=n, D=D, closure_positions=len(closures), pairs=npairs)
    clock = time.monotonic()
    timings = {}

    def lap(name):
        nonlocal clock
        now = time.monotonic()
        timings[name] = round(now - clock, 2)
        clock = now
        if log:
            log(f'  {name}: {timings[name]} s')

    # A. binary quadrics from the closures
    blocks = [Pn.T[:, pa] & Pn.T[:, pb]]
    for i, U in closures.items():
        blocks.append(polar_rows(Pn[:, i], U, pa, pb))
    Msys = to_f2(np.vstack(blocks))
    del blocks
    Q = Msys.right_kernel_matrix()
    rep['binary_equations'] = Msys.nrows()
    rep['binary_quadrics'] = Q.nrows()
    del Msys
    lap('A_binary_quadrics')
    if Q.nrows() == 0:
        rep['status'] = 'no_binary_quadrics'
        return rep, None, timings

    # B. branch separation at one position
    order = [int(x) for x in rng.permutation(n)]
    chosen = None
    attempts = []
    for i0 in order[:10]:
        p = Pn[:, i0]
        J = Q * to_f2(polar_field(fld, p.astype(np.uint16), pa, pb, k))
        W = J.right_kernel_matrix()
        basis = [vector(F2, p.tolist())]
        Brows = []
        for w in W.rows():
            if matrix(F2, basis + [w]).rank() > len(basis):
                basis.append(w)
                Brows.append(w)
        Brows = Brows[:max(0, W.nrows() - 1)]
        B = np.array([list(map(int, w)) for w in Brows], dtype=np.uint8).reshape(len(Brows), k)
        b = B.shape[0]
        info = dict(position=i0, gradient_nullity=W.nrows())
        if W.nrows() != m + 1:
            attempts.append(info)
            continue
        wp = list(combinations_with_replacement(range(b), 2))
        mono = np.zeros((len(wp), npairs), dtype=np.uint8)
        for r_, (u, v) in enumerate(wp):
            if u == v:
                mono[r_] = B[u][pa] & B[u][pb]
            else:
                mono[r_] = (B[u][pa] & B[v][pb]) ^ (B[v][pa] & B[u][pb])
        if goppa_t:
            Q2 = Q
        else:
            Q2 = J.left_kernel_matrix() * Q            # quadrics singular at p
        info['second_family'] = int(Q2.nrows())
        Rm = Q2 * to_f2(mono.T)
        info['quadratic_rank'] = int(Rm.rank())
        attempts.append(info)
        if info['quadratic_rank'] != m * (m - 1) // 2:
            continue
        ann = Rm.right_kernel_matrix()
        if ann.nrows() != m:
            continue

        def sym(vec):
            S = matrix(F2, b, b)
            for r_, (u, v) in enumerate(wp):
                S[u, v] = vec[r_]
                S[v, u] = vec[r_]
            return S
        z = None
        for _ in range(200):
            ca = vector(F2, rng.integers(0, 2, m).tolist())
            cb = vector(F2, rng.integers(0, 2, m).tolist())
            SA, SB = sym(ca * ann), sym(cb * ann)
            if not SB.is_invertible():
                continue
            N = SA * SB.inverse()
            cp = N.charpoly()
            if not cp.is_irreducible():
                continue
            theta = cp.change_ring(E).roots()[0][0]
            ker = (N.change_ring(E) - theta).right_kernel()
            if ker.dimension() == 1:
                z = ker.basis()[0]
                break
        if z is None:
            continue
        v1 = np.zeros(k, dtype=np.uint16)
        for u in range(b):
            zu = int(z[u].to_integer())
            if zu:
                v1 ^= np.where(B[u] == 1, zu, 0).astype(np.uint16)
        chosen = (i0, p, J, v1, Q2)
        break
    rep['separation_attempts'] = attempts
    lap('B_branch_separation')
    if chosen is None:
        rep['status'] = 'branch_separation_failed'
        return rep, None, timings
    i0, p, J, v1, Q2 = chosen
    rep['jet_position'] = i0
    prodv = fld.mul(v1[pa], v1[pb])
    assert not binary_times_field(fld, Q2, prodv[:, None]).any(), 'h(v1) != 0'

    # C. deep jet by continuation, v_{2j} = v_j^2
    Jn = f2_array(J).astype(np.uint16)
    Ln = binary_times_field(fld, Q2, polar_field(fld, v1, pa, pb, k))
    A = fld.to_sage(np.vstack([Jn, Ln]))
    rankA = A.rank()
    rep['continuation_rank'] = int(rankA)
    if rankA != k - 2:
        rep['status'] = 'continuation_rank_deficient'
        lap('C_continuation')
        return rep, None, timings
    prow = list(A.pivot_rows())
    pcol = list(A.pivots())
    inv = fld.from_sage(A.matrix_from_rows_and_columns(prow, pcol).inverse())
    loginv = fld.log[inv]
    dq = Q.nrows()
    total_rows = dq + Q2.nrows()
    extra = [int(x) for x in rng.choice(total_rows, size=min(16, total_rows), replace=False)
             if int(x) not in set(prow)]
    used = prow + extra
    Qn = f2_array(Q)
    Q2n = Qn if goppa_t else f2_array(Q2)
    rowsJ = [r_ for r_ in used if r_ < dq]
    rowsL = [r_ - dq for r_ in used if r_ >= dq]
    QJ = Qn[rowsJ].astype(np.float32)
    QL = Q2n[rowsL].astype(np.float32)
    posJ = {r_: c for c, r_ in enumerate(rowsJ)}
    posL = {x: c for c, x in enumerate(rowsL)}
    Aused = np.vstack([Jn, Ln])[used]
    logAused = fld.log[Aused]
    need = 2 * D - n + 1                                 # n + r > 2D
    top = need + 2 + (need % 2)
    V = np.zeros((top + 2, k), dtype=np.uint16)
    V[0], V[1] = p, v1
    if goppa_t:
        V[2] = fld.mul(v1, v1)
    logV = fld.log[V]
    inconsistent = 0
    r_ = 3 if goppa_t else 2
    while r_ <= need + 1:
        V[r_] = 0
        if goppa_t:
            V[r_ + 1] = fld.mul(V[(r_ + 1) // 2], V[(r_ + 1) // 2])
        else:
            V[r_ + 1] = 0
        logV[r_] = fld.log[V[r_]]
        logV[r_ + 1] = fld.log[V[r_ + 1]]
        cJ = rows_times_field(fld, QJ, full_coeff(fld, logV, r_, pa, pb))
        cL = rows_times_field(fld, QL, full_coeff(fld, logV, r_ + 1, pa, pb))
        rhs = np.array([cJ[posJ[x]] if x < dq else cL[posL[x - dq]] for x in used],
                       dtype=np.uint16)
        sol = np.bitwise_xor.reduce(fld.exp[loginv + fld.log[rhs[:len(prow)]][None, :]], axis=1)
        x = np.zeros(k, dtype=np.uint16)
        x[pcol] = sol
        check = np.bitwise_xor.reduce(fld.exp[logAused + fld.log[x][None, :]], axis=1)
        if np.any(check != rhs):
            inconsistent += 1
        V[r_] = x
        logV[r_] = fld.log[x]
        r_ += 2 if goppa_t else 1
    order_r = r_ - 1                                     # V[0..order_r] is known
    if not goppa_t:
        V[order_r + 1:] = 0
        logV = fld.log[V]
    rep['jet_order'] = int(order_r)
    rep['required_order'] = int(need)
    rep['continuation_inconsistent_orders'] = int(inconsistent)
    del Qn, Q2n, QJ, QL
    lap('C_continuation')
    if inconsistent:
        rep['status'] = 'continuation_inconsistent'
        return rep, None, timings

    # D. quadrics over K on the selected branch (all monomials)
    r_use = need
    rows = [(Pn.T[:, pa] & Pn.T[:, pb]).astype(np.uint16)]
    jet_rows = np.zeros((r_use, npairs), dtype=np.uint16)
    for e in range(1, r_use + 1):
        jet_rows[e - 1] = full_coeff(fld, logV, e, pa, pb)
    rows.append(jet_rows)
    S = fld.to_sage(np.vstack(rows))
    del rows, jet_rows
    Rr = S.echelon_form()
    piv = list(Rr.pivots())
    rank = len(piv)
    pivset = set(piv)
    free = [c for c in range(npairs) if c not in pivset]
    rep['K_equations'] = int(S.nrows())
    rep['K_rank'] = int(rank)
    rep['K_rank_expected'] = int(2 * D + 1 - t)
    rep['K_quadrics'] = int(len(free))
    del S
    if len(free) < k - 2:
        rep['status'] = 'too_few_K_quadrics'
        lap('D_K_quadrics')
        return rep, None, timings
    # random elements of the solution space: free coordinates at random
    Hrows = min(len(free), k + 40)
    Fn = rng.integers(0, fld.q, size=(len(free), Hrows)).astype(np.uint16)
    Pp = Rr.matrix_from_rows_and_columns(list(range(rank)), free) * fld.to_sage(Fn)
    Kn = np.zeros((Hrows, npairs), dtype=np.uint16)
    Kn[:, free] = Fn.T
    Kn[:, piv] = fld.from_sage(Pp).T
    del Rr, Pp, Fn
    lap('D_K_quadrics')
    C = np.zeros((Hrows, k, k), dtype=np.uint16)
    off = pa != pb
    C[:, pa[off], pb[off]] = Kn[:, off]
    C[:, pb[off], pa[off]] = Kn[:, off]

    # E. aligned tangent planes at every public column
    Yhat = np.zeros((k, n), dtype=np.uint16)
    bad = []
    for j in range(n):
        Jj = np.bitwise_xor.reduce(C[:, :, Pn[:, j].astype(bool)], axis=2)
        ker = fld.to_sage(Jj).right_kernel_matrix()
        if ker.nrows() != 2:
            bad.append((j, int(ker.nrows())))
            continue
        kn = fld.from_sage(ker)
        yj = Pn[:, j].astype(np.uint16)
        pick = None
        for row in kn:
            if matrix(E, [fld.to_sage(row[None, :]).row(0), fld.to_sage(yj[None, :]).row(0)]).rank() == 2:
                pick = row
                break
        Yhat[:, j] = pick
    rep['gradient_rank_failures'] = [list(x) for x in bad]
    lap('E_tangent_planes')
    return rep, dict(Yhat=Yhat, bad=[j for j, _ in bad], V=V[:order_r + 1], i0=i0), timings


def minor_and_support(E, fld, Pn, Yhat, bad, D, t, seed, rep):
    """F. minor code on the good columns, support, embedding of the public code.
    Here t is the Goppa degree for a binary Goppa code and 0 otherwise."""
    k, n = Pn.shape
    good = [j for j in range(n) if j not in set(bad)]
    Yg = Pn[:, good]
    Hg = Yhat[:, good]
    ng = len(good)
    pr = MC.sample_pairs(k, MC.row_budget(k, D), seed)
    if pr is None:
        pr = [(a, b) for a in range(k) for b in range(a + 1, k)]
    ia = np.array([a for a, _ in pr])
    ib = np.array([b for _, b in pr])
    Mn = np.where(Yg[ia] == 1, Hg[ib], 0) ^ np.where(Yg[ib] == 1, Hg[ia], 0)
    T = fld.to_sage(fld.sqrt(Mn.astype(np.uint16)))
    rank = T.rank()
    rep['minor_rows'] = len(pr)
    rep['minor_rank'] = int(rank)
    rep['minor_expected'] = int(min(D - t, ng))
    rep['minor_length'] = ng
    if rank != D - t or rank >= ng - 1:
        rep['status'] = 'minor_code_not_full'
        return None
    basis = T.echelon_form()[:rank]
    Gs = basis.right_kernel_matrix() if rank > ng - rank else basis
    support = MC.sidelnikov_shestakov(E, Gs, tries=int(E.order()))
    ok = support is not None and MC.grs_multiplier_systematic(E, Gs, support) is not None
    rep['support_recovered'] = bool(ok)
    if not ok:
        rep['status'] = 'support_recovery_failed'
        return None
    Y = to_f2(Yg).change_ring(E)
    lm = MC.grs_multiplier_fast(E, Y, support, D + 1)
    emb = lm is not None and (Y * MC.grs_parity(E, support, D + 1, lm).transpose()).is_zero()
    rep['embedding_verified'] = bool(emb)
    rep['status'] = 'recovered' if emb else 'embedding_failed'
    return dict(good=good, support=support)


# ---------------------------------------------------------------- the audits

def audit_alignment(E, fld, Pn, Yhat, Yhat_true, bad):
    """Which Frobenius conjugate of the true tangent lies in each recovered plane."""
    k, n = Pn.shape
    counts = {}
    for j in range(n):
        if j in bad:
            continue
        A = matrix(E, 2, k, [fld.elts[int(x)] for x in Pn[:, j]] +
                   [fld.elts[int(x)] for x in Yhat[:, j]])
        hit = None
        for s in range(fld.m):
            w = fld.frob(Yhat_true[:, j], s)
            if A.stack(vector(E, [fld.elts[int(x)] for x in w])).rank() == 2:
                hit = s
                break
        counts[hit] = counts.get(hit, 0) + 1
    return {str(a): b for a, b in counts.items()}


def audit_support(E, fld, support, good, secret):
    """Is the recovered support a Moebius image of a conjugate of the secret one?"""
    al = secret['alpha'][good]
    x = support

    def cross(a0, a1, a2, a3):
        return ((a3 - a0) * (a2 - a1)) / ((a3 - a1) * (a2 - a0))
    for s in range(fld.m):
        a = [fld.elts[int(v)] for v in fld.frob(al, s)]
        if all(cross(x[0], x[1], x[2], x[j]) == cross(a[0], a[1], a[2], a[j])
               for j in range(3, len(x))):
            return s
    return None


# ---------------------------------------------------------------- the driver

def run(family, m, n, t, seed, mode, extra=0, log=None):
    started = time.monotonic()
    E, fld, P, secret = instance(family, m, n, t, seed)
    Pn = f2_array(P)
    k = Pn.shape[0]
    goppa_t = t if family == 'goppa' else 0
    D = n - 2 * t - 1 if family == 'goppa' else n - t - 1
    if family != 'goppa':
        mode = 'proved'
    count = {'saving': 2 * D + 1 - goppa_t - n, 'proved': 2 * D + 1 - n}[mode]
    count = max(0, min(n, count + extra))
    rng = np.random.default_rng(seed + 10 ** 6)
    positions = sorted(int(x) for x in rng.choice(n, size=count, replace=False))
    Yhat_true = true_tangents(fld, Pn, secret)
    closures = closures_from_secret(fld, Pn, Yhat_true, positions)
    rep = dict(family=family, m=m, n=n, t=t, k=k, seed=seed, mode=mode, extra_closures=extra,
               closure_dimensions=sorted(set(int(U.shape[0]) for U in closures.values())),
               predicted_binary_quadrics=(max(0, comb(k + 1, 2) - (n + 1) - m * (n - 5 * t - 2))
                                          if family == 'goppa' else None),
               generation_seconds=round(time.monotonic() - started, 2))
    r, data, timings = recover(E, fld, Pn, closures, m, D, goppa_t, seed, log)
    rep.update(r)
    if data is not None:
        tic = time.monotonic()
        out = minor_and_support(E, fld, Pn, data['Yhat'], data['bad'], D, goppa_t, seed, rep)
        timings['F_minor_code'] = round(time.monotonic() - tic, 2)
        tic = time.monotonic()
        rep['audit_tangent_conjugates'] = audit_alignment(E, fld, Pn, data['Yhat'], Yhat_true,
                                                          set(data['bad']))
        if out is not None:
            s = audit_support(E, fld, out['support'], out['good'], secret)
            rep['audit_support_conjugate'] = s
            rep['audit_support_moebius'] = s is not None
        timings['audits'] = round(time.monotonic() - tic, 2)
    rep['timings'] = timings
    rep['total_seconds'] = round(time.monotonic() - started, 2)
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', default='goppa', choices=['goppa', 'alternant'])
    ap.add_argument('--m', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--t', type=int, required=True)
    ap.add_argument('--seeds', required=True, help='a-b inclusive, or a comma list')
    ap.add_argument('--mode', default='saving', choices=['saving', 'proved'])
    ap.add_argument('--extra', type=int, default=0,
                    help='closure positions beyond the count of the mode')
    ap.add_argument('--output', required=True)
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()
    if '-' in args.seeds:
        lo, hi = map(int, args.seeds.split('-'))
        seeds = list(range(lo, hi + 1))
    else:
        seeds = list(map(int, args.seeds.split(',')))
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if 'error' not in row:
                done.add(row['seed'])
    log = (lambda s: print(s, flush=True)) if args.verbose else None
    for seed in seeds:
        if seed in done:
            continue
        try:
            rep = run(args.family, args.m, args.n, args.t, seed, args.mode, args.extra, log)
        except Exception:
            rep = dict(family=args.family, m=args.m, n=args.n, t=args.t, seed=seed, mode=args.mode,
                       error=traceback.format_exc())
        with path.open('a') as handle:
            handle.write(json.dumps(rep, default=int) + '\n')
        short = {a: rep.get(a) for a in ('seed', 'k', 'status', 'binary_quadrics',
                 'predicted_binary_quadrics', 'jet_order', 'minor_rank',
                 'audit_support_moebius', 'total_seconds', 'error') if a in rep}
        print(json.dumps(short, default=int), flush=True)


if __name__ == '__main__':
    main()
