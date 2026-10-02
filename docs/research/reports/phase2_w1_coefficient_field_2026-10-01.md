# Exact W1 coefficient-field validation — 2026-10-01

## Result

The joint phase coefficient algebra for the pi/3 W1 candidate is a field. Let `K = Q(q)`, where `q^12 - q^6 + 1 = Phi_36(q)`. Exact Singular calculations give a degree-2 gcd for each pair `(P(a), L0(q,a))`, `(P(b), L1(q,b))`, and `(P(c), L2(q,c))`. Each gcd is one irreducible quadratic over K, with multiplicity 1. The monic factors are:

```text
f_a(u) = u^2 + (16q^11 - 35q^10 + 170q^9 - 18q^8 + 120q^7 + 32q^6
                - 166q^5 - 157q^4 + 8q^3 - 46q^2 - 48q - 66)/323
f_b(u) = u^2 + (-166q^11 + 192q^10 + 170q^9 - 46q^8 - 72q^7 + 32q^6
                + 150q^5 - 35q^4 + 8q^3 + 64q^2 + 120q - 66)/323
f_c(u) = u^2 + (150q^11 - 157q^10 + 170q^9 + 64q^8 - 48q^7 + 32q^6
                + 16q^5 + 192q^4 + 8q^3 - 18q^2 - 72q - 66)/323
```

Put `d_i = -f_i(0)`, so the three phase coordinates satisfy `u_i^2 = d_i`. Singular factored `y^2 - product(d_i)` for all seven nonempty subsets of `{a,b,c}` over K. Each was a single irreducible quadratic. Thus no nontrivial product of the three square classes is a square in K; the classes are independent in `K*/K*2`. Their compositum has degree 8 over K, so the joint phase quotient is a field. Since `[K:Q] = 12`, its degree over Q is 96.

The inverse phase variables are redundant: `gcd(u,P(u)) = 1` (the constant term of P is 104329), so each phase coordinate is a unit in this field. The earlier 30-minute witness run remains inconclusive; this coefficient-field result does not prove W1-emptiness.

A direct `primdecGTZ` computation did not return component data. The field conclusion instead follows from the exact seven-test square-class criterion above and the independently computed joint vector-space dimension 8.

## Reformulation and execution status

`symbolic_export/w1_pi3_exact_abc.m2` is now formulated over Q[q,a,b,c] modulo Phi_36 and the three exact quadratic factors. It maps q,a,b,c into the resulting field, substitutes `a_i^(-1) = 1/a_i`, then constructs the same 17 witness equations. The question remains local to this fixed exact triple.

The Macaulay2 field smoke test and witness rerun are not yet verified. WSL refused to create a VM with `CreateVm/HCS/0x800705aa`, reporting that the host paging file is too small. No witness-elimination conclusion is claimed.

## Reproducible artifacts

- Exact Singular gcd/factorization: `symbolic_export/w1_pi3_coefficient_field_preflight.sing`
- Exact square-class tests: `symbolic_export/w1_pi3_quadratic_square_classes.sing`
- Gcd transcript: `results/phase2_w1_coefficient_field_preflight.log`
- Seven square-class factorizations: `results/phase2_w1_quadratic_square_classes.log`
- Reformulated witness input: `symbolic_export/w1_pi3_exact_abc.m2`
- Macaulay2 field smoke-test source: `symbolic_export/w1_pi3_m2_field_smoke.m2`
