import sys
import sympy as sp

print("=== PHASE 1: INDEPENDENT ANCHOR AUDIT ===")

# 1. Algebraic constants
omega3 = sp.exp(2 * sp.pi * sp.I / 3)
omega24 = sp.exp(2 * sp.pi * sp.I / 24)
b2 = (1 - 2*sp.I) / sp.sqrt(5)
D_defect = sp.diag(omega24**9 * b2, 1, 1)

# 2. Reconstruct PUBLISHED D_bc from C3 and C4 blocks (Bengtsson et al. 2007, Eqs. 79-80)
C3 = sp.Matrix([
    [1, sp.I, sp.I],
    [sp.I, 1, sp.I],
    [sp.I, sp.I, 1]
])
C4 = sp.Matrix([
    [sp.exp(5*sp.pi*sp.I/4), sp.exp(sp.pi*sp.I/4), sp.exp(sp.pi*sp.I/4)],
    [sp.exp(sp.pi*sp.I/4), sp.exp(5*sp.pi*sp.I/4), sp.exp(sp.pi*sp.I/4)],
    [sp.exp(sp.pi*sp.I/4), sp.exp(sp.pi*sp.I/4), sp.exp(5*sp.pi*sp.I/4)]
])

D_bc = sp.Matrix.vstack(
    sp.Matrix.hstack(C3, C4),
    sp.Matrix.hstack(C4, -sp.I * C3.H)
)

# 3. Reconstruct PUBLISHED F_D (Bengtsson et al. 2007, Eq. 78)
F3 = sp.Matrix([
    [1, 1, 1],
    [1, omega3, omega3**2],
    [1, omega3**2, omega3]
])
F_D = sp.Matrix.vstack(
    sp.Matrix.hstack(F3, F3),
    sp.Matrix.hstack(F3 * D_defect, -F3 * D_defect)
)

# 4. Construct CANDIDATE H_Dita(i) and B3(i) from the proposed orbit
Delta = sp.diag(1, sp.I, -1)
H_Dita_i = sp.Matrix.vstack(
    sp.Matrix.hstack(F3, F3),
    sp.Matrix.hstack(Delta * F3, -Delta * F3)
)

K = sp.diag(1, 1, omega3) * F3
B3_i = sp.Matrix.vstack(
    sp.Matrix.hstack(K, K),
    sp.Matrix.hstack(sp.I * Delta * K, -sp.I * Delta * K)
)

# 5. Check Unitarity / Hadamard property (unnormalized: M * M.H == 6 * I)
print("\n--- Unitarity Checks (M * M^dag == 6 * I_6) ---")
u_Dbc = sp.simplify(D_bc * D_bc.H - 6 * sp.eye(6)).is_zero_matrix
u_FD = sp.simplify(F_D * F_D.H - 6 * sp.eye(6)).is_zero_matrix
u_HDita = sp.simplify(H_Dita_i * H_Dita_i.H - 6 * sp.eye(6)).is_zero_matrix
u_B3 = sp.simplify(B3_i * B3_i.H - 6 * sp.eye(6)).is_zero_matrix

print(f"D_bc unitary:     {u_Dbc}")
print(f"F_D unitary:      {u_FD}")
print(f"H_Dita(i) unitary:{u_HDita}")
print(f"B3(i) unitary:    {u_B3}")

if not (u_Dbc and u_FD and u_HDita and u_B3):
    print("[FATAL] One or more base matrices are not unitary! Stopping.")
    sys.exit(1)

# 6. Mutual Unbiasedness Checks: |(M1^dag M2)_jk|^2 == 6
print("\n--- MUB Property Checks (|(M1^dag M2)_jk|^2 == 6) ---")
cross_pub = sp.simplify(F_D.H * D_bc)
cross_cand = sp.simplify(B3_i.H * H_Dita_i)

mod_pub = set([sp.simplify(sp.Abs(cross_pub[r, c])**2) for r in range(6) for c in range(6)])
mod_cand = set([sp.simplify(sp.Abs(cross_cand[r, c])**2) for r in range(6) for c in range(6)])

print(f"Published pair {{F_D, D_bc}} MUB moduli: {mod_pub}")
print(f"Candidate pair {{B3(i), H_Dita(i)}} MUB moduli: {mod_cand}")

is_pub_mub = (mod_pub == {6})
is_cand_mub = (mod_cand == {6})
print(f"Published pair forms valid MUB: {is_pub_mub}")
print(f"Candidate pair forms valid MUB: {is_cand_mub}")

# 7. Direct Equality Checks
print("\n--- Direct Equality Checks ---")
diff_H = sp.simplify(H_Dita_i - D_bc).is_zero_matrix
diff_B = sp.simplify(B3_i - F_D).is_zero_matrix
print(f"H_Dita(i) == D_bc directly: {diff_H}")
print(f"B3(i) == F_D directly:      {diff_B}")

# 8. Dephasing and Equivalence Check
def dephase(M):
    # Dephase rows so col 0 is all 1, then cols so row 0 is all 1
    D_r = sp.diag(*[1/M[r, 0] for r in range(6)])
    M1 = D_r * M
    D_c = sp.diag(*[1/M1[0, c] for c in range(6)])
    return sp.simplify(M1 * D_c)

print("\n--- Dephased Normal Form Comparison ---")
d_Dbc = dephase(D_bc)
d_HDita = dephase(H_Dita_i)
diff_dephased_H = sp.simplify(d_HDita - d_Dbc).is_zero_matrix
print(f"Dephased H_Dita(i) == Dephased D_bc: {diff_dephased_H}")

if diff_dephased_H:
    print("[SUCCESS] H_Dita(i) is directly CHM-equivalent to D_bc without permutation!")
else:
    print("[INFO] Testing permutation equivalence...")
    # Permutation search over columns
    found_equiv = False
    import itertools
    for p in itertools.permutations(range(6)):
        P = sp.eye(6)[:, list(p)]
        if sp.simplify(dephase(H_Dita_i * P) - d_Dbc).is_zero_matrix:
            print(f"[SUCCESS] H_Dita(i) equivalent to D_bc under column permutation {p}!")
            found_equiv = True
            break
    if not found_equiv:
        print("[WARNING] H_Dita(i) did not match D_bc under direct dephasing/col permutations.")
