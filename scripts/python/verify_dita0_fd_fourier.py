"""Exact audit of Bengtsson's D6 third basis against the D(1) representative.

This standalone checker uses only Python's fractions/integers.  It works in
Q(zeta_24, sqrt(5)) with zeta_24^8 = zeta_24^4 - 1 and checks:
  * D(1) and F_D are order-six Hadamard matrices;
  * F_D is unbiased to I and D(1);
  * after exact row/column permutations and dephasing, each transposed edge
    is (or is not) in Jaming's two-parameter Fourier family.

Run normally for the exact audit, or with --self-test to also perturb one
coefficient and confirm that the verifier rejects the damaged matrix.
"""

from fractions import Fraction as Q
from itertools import permutations
import sys

N = 8  # [Q(zeta_24):Q]
ZERO = tuple(Q(0) for _ in range(N))
ONE = (Q(1),) + tuple(Q(0) for _ in range(N - 1))


def qred(v):
    v = list(v) + [Q(0)] * max(0, 15 - len(v))
    for k in range(len(v) - 1, 7, -1):
        a = v[k]
        if a:
            v[k] = 0
            v[k - 4] += a
            v[k - 8] -= a
    return tuple(v[:N])


def qadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def qneg(a):
    return tuple(-x for x in a)


def qmul(a, b):
    v = [Q(0)] * 15
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    v[i + j] += x * y
    return qred(v)


def qscale(a, s):
    return tuple(x * s for x in a)


def zpow(k):
    k %= 24
    v = [Q(0)] * (k + 1)
    v[k] = Q(1)
    return qred(v)


def qconj(a):
    out = ZERO
    for k, coeff in enumerate(a):
        if coeff:
            out = qadd(out, qscale(zpow(-k), coeff))
    return out


# A field element is (a,b), representing a + b*sqrt(5), a,b in Q(zeta_24).
E0 = (ZERO, ZERO)
E1 = (ONE, ZERO)
SQRT5 = (ZERO, ONE)


def eadd(x, y):
    return (qadd(x[0], y[0]), qadd(x[1], y[1]))


def eneg(x):
    return (qneg(x[0]), qneg(x[1]))


def escale(x, s):
    return (qscale(x[0], s), qscale(x[1], s))


def emul(x, y):
    return (
        qadd(qmul(x[0], y[0]), qscale(qmul(x[1], y[1]), Q(5))),
        qadd(qmul(x[0], y[1]), qmul(x[1], y[0])),
    )


def econj(x):
    return (qconj(x[0]), qconj(x[1]))


def epowz(k):
    return (zpow(k), ZERO)


def ematmul(A, B):
    return [[sum_e(emul(A[i][k], B[k][j]) for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def sum_e(items):
    total = E0
    for item in items:
        total = eadd(total, item)
    return total


def transpose(A):
    return [list(row) for row in zip(*A)]


def dagger(A):
    return [[econj(x) for x in row] for row in transpose(A)]


def is_zero(x):
    return x == E0


def hadamard_defect(A):
    G = ematmul(A, dagger(A))
    return all(is_zero(eadd(G[i][j], eneg(escale(E1, Q(6) if i == j else Q(0)))))
               for i in range(6) for j in range(6))


def scalar_norm(x):
    return emul(x, econj(x))


def fourier3():
    w = epowz(8)
    w2 = emul(w, w)
    return [[E1, E1, E1], [E1, w, w2], [E1, w2, w]]


def make_fd():
    F3 = fourier3()
    i = epowz(6)
    defect = escale(emul(epowz(9), eadd(E1, eneg(escale(i, Q(2))))), Q(1, 5))
    defect = emul(defect, SQRT5)  # zeta_24^9 (1-2i)/sqrt(5)
    diag = [defect, E1, E1]
    F3D = [[emul(F3[r][c], diag[c]) for c in range(3)] for r in range(3)]
    return [row + row[:] for row in F3] + [
        F3D[r] + [eneg(x) for x in F3D[r]] for r in range(3)
    ]


def make_d1():
    # User's D(z) from the paper, specialized at z=1.
    o, m, i, mi = epowz(0), epowz(12), epowz(6), epowz(18)
    return [
        [o, o, o, o, o, o],
        [o, m, o, m, i, mi],
        [o, mi, i, i, mi, m],
        [o, i, m, o, m, mi],
        [o, o, mi, m, m, i],
        [o, m, m, mi, o, i],
    ]


def make_dbc():
    # Bengtsson et al.'s block-circulant D6 representative, Eqs. (79)-(80).
    exponents = [
        [0, 6, 6, 15, 3, 3],
        [6, 0, 6, 3, 15, 3],
        [6, 6, 0, 3, 3, 15],
        [15, 3, 3, 18, 12, 12],
        [3, 15, 3, 12, 18, 12],
        [3, 3, 15, 12, 12, 18],
    ]
    return [[epowz(k) for k in row] for row in exponents]


def align_fd_to_d1(FD):
    """Apply the exact monomial map D_bc -> D(1) to both bases."""
    row_perm = (0, 1, 3, 2, 4, 5)
    col_perm = (3, 2, 4, 5, 1, 0)
    row_phase = (9, 21, 6, 21, 12, 12)
    col_phase = (0, 9, 12, 12, 9, 15)
    D, Dbc = make_d1(), make_dbc()
    for r in range(6):
        for c in range(6):
            rhs = emul(emul(epowz(row_phase[r]), Dbc[row_perm[r]][col_perm[c]]),
                       epowz(col_phase[c]))
            assert D[r][c] == rhs, ("exact D_bc -> D(1) equivalence", r, c)
    aligned = [[emul(epowz(row_phase[r]), FD[row_perm[r]][c])
                for c in range(6)] for r in range(6)]
    return aligned, row_perm, col_perm, row_phase, col_phase


def flatten_equal(A, B):
    return len(A) == len(B) and all(
        len(a) == len(b) and all(x == y for x, y in zip(a, b))
        for a, b in zip(A, B)
    )


def dephase(A, norm_sq):
    # C_jk = A_jk A_00 / (A_j0 A_0k), using conjugation and |A_jk|^2=norm_sq.
    return [[escale(emul(emul(emul(A[j][k], A[0][0]),
                                  econj(A[j][0])), econj(A[0][k])),
                      Q(1, norm_sq * norm_sq))
             for k in range(6)] for j in range(6)]


def fourier_family(x, y):
    w = epowz(8)
    w2 = emul(w, w)
    m1 = eneg(E1)
    return [
        [E1, E1, E1, E1, E1, E1],
        [E1, eneg(emul(w2, x)), w, eneg(x), w2, eneg(emul(w, x))],
        [E1, emul(w, y), w2, y, w, emul(w2, y)],
        [E1, m1, E1, m1, E1, m1],
        [E1, emul(w2, x), w, x, w2, emul(w, x)],
        [E1, eneg(emul(w, y)), w2, eneg(y), w, eneg(emul(w2, y))],
    ]


def chart_match(A, norm_sq):
    """Find row/column permutations of A^T yielding a dephased F(x,y)."""
    T = transpose(A)
    perms = list(permutations(range(6)))
    # First row r0 and fixed alternating row r3 determine the eligible columns.
    for r0 in range(6):
        for r3 in range(6):
            if r3 == r0:
                continue
            plus, minus = [], []
            for c in range(6):
                ratio = escale(emul(emul(emul(T[r3][c], T[r0][0]),
                                               econj(T[r3][0])), econj(T[r0][c])),
                                    Q(1, norm_sq * norm_sq))
                if ratio == E1:
                    plus.append(c)
                elif ratio == eneg(E1):
                    minus.append(c)
                else:
                    plus = []
                    break
            if len(plus) != 3 or len(minus) != 3:
                continue
            # Canonical Fourier row 3 has signs + - + - + -.  Its first column
            # must be one of the plus positions; arrange the other columns.
            for c0 in plus:
                other_plus = [c for c in plus if c != c0]
                for pp in permutations(other_plus):
                    for mm in permutations(minus):
                        col_perm = (c0, mm[0], pp[0], mm[1], pp[1], mm[2])
                        remaining = [r for r in range(6) if r not in (r0, r3)]
                        for rr in permutations(remaining):
                            row_perm = (r0, rr[0], rr[1], r3, rr[2], rr[3])
                            Nmat = [[T[row_perm[j]][col_perm[k]] for k in range(6)]
                                    for j in range(6)]
                            C = dephase(Nmat, norm_sq)
                            x = eneg(C[1][3])
                            y = C[2][3]
                            if flatten_equal(C, fourier_family(x, y)):
                                assert scalar_norm(x) == E1
                                assert scalar_norm(y) == E1
                                return row_perm, col_perm, x, y
    return None


def check_exact():
    D = make_d1()
    FD = make_fd()
    assert hadamard_defect(D), "D(1) is not an exact CHM"
    assert hadamard_defect(FD), "F_D is not an exact CHM"
    F, row_perm, col_perm, row_phase, col_phase = align_fd_to_d1(FD)
    assert hadamard_defect(F), "aligned F_D is not an exact CHM"
    assert all(scalar_norm(F[r][c]) == E1 for r in range(6) for c in range(6)), \
        "F_D is not unbiased to the computational basis"
    W = ematmul(dagger(D), F)
    assert all(scalar_norm(W[r][c]) == escale(E1, Q(6))
               for r in range(6) for c in range(6)), \
        "F_D is not unbiased to D(1)"
    first = chart_match(F, 1)
    second = chart_match(W, 6)
    return first, second, (row_perm, col_perm, row_phase, col_phase)


def format_field(x):
    terms = []
    names = [f"z^{k}" for k in range(N)]
    for radical_part, label in ((x[0], ""), (x[1], "sqrt(5)*")):
        for k, coeff in enumerate(radical_part):
            if coeff:
                terms.append(f"({coeff})*{label}{names[k]}")
    return " + ".join(terms) if terms else "0"


def main():
    first, second, equivalence = check_exact()
    print("Exact field: Q(zeta_24, sqrt(5)); z = exp(i*pi/12), z^8-z^4+1=0")
    print("D(1) CHM: PASS")
    print("F_D CHM: PASS")
    print("Exact D_bc -> D(1) row/column permutation and dephasing witness: PASS")
    print("  row permutation:", equivalence[0])
    print("  column permutation:", equivalence[1])
    print("  row phase exponents of zeta_24:", equivalence[2])
    print("  column phase exponents of zeta_24:", equivalence[3])
    print("Aligned F_D unbiased to I and D(1): PASS")
    if first:
        print("I -> F_D: exact transposed Fourier chart FOUND")
        print("  row permutation of edge transpose (zero-based):", first[0])
        print("  column permutation of edge transpose (zero-based):", first[1])
        print("  x =", format_field(first[2]))
        print("  y =", format_field(first[3]))
    else:
        print("I -> F_D: no exact transposed Fourier chart in full row/column permutation search")
    if second:
        print("D(1) -> F_D: exact transposed Fourier chart FOUND")
        print("  row permutation of edge transpose (zero-based):", second[0])
        print("  column permutation of edge transpose (zero-based):", second[1])
        print("  x =", format_field(second[2]))
        print("  y =", format_field(second[3]))
    else:
        print("D(1) -> F_D: no exact transposed Fourier chart in full row/column permutation search")
    if "--self-test" in sys.argv:
        damaged = make_fd()
        damaged[0][0] = epowz(1)  # one exact coefficient changed from 1 to zeta_24
        assert not hadamard_defect(damaged), "mutation test failed to detect perturbed coefficient"
        print("Mutation self-test: perturbed FD[0,0] rejected exactly: PASS")
    if (first is None) != (second is None):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
