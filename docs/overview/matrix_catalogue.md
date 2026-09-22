# Matrix Catalogue

All matrices below use the repository's unnormalised convention: entries have modulus 1 and `H H^dagger = 6 I`. Normalised bases are obtained by dividing columns by `sqrt(6)`.

| Matrix | Formula source | Normalization | Field | CHM verified? | MUB relations | Equivalence notes |
|---|---|---|---|---|---|---|
| `F6` | Brierley--Weigert Eq. (15); repository `F6` | Unnormalised here | `Q(zeta6)` | Yes; repo tests and BW/Grassl | `{I,F6}`; 48 MU vectors, 16 third bases; no fourth | Heisenberg Fourier representative |
| `F3` | Bengtsson et al. Eq. (61) block construction | Unnormalised here | `Q(zeta3)` | Yes | Used inside `F_D` | `F_D` is twisted-product, not ordinary `F6` |
| `D(x)` | Bengtsson et al. Eq. (11); BW Appendix A Eq. (34) | Source prints normalized; repo drops `1/sqrt(6)` | `Q(i,z)` with `z=exp(2 pi i x)` | Yes | Diţă family; D0 has exact BW result | `D0=D(0)`; `D(x)` equivalences include `x -> x+1/2` and `x -> -x+1/4` |
| `D0` | Bengtsson Eq. (11) at `x=0` | Unnormalised here | `Q(i)` | Yes | 120 MU vectors, 10 third bases, no pair MU | Exact BW nonextendability for triples containing D0 |
| `D_bc` | Bengtsson Eqs. (79)-(80) | Source normalized; repo unnormalised | `Q(zeta24)` | Yes, repo exact/HP checks | In `{I,F_D,D_bc}` | Row+column equivalent to D0; column-only comparison is insufficient |
| `F_D` | Bengtsson Eqs. (61), (78), with `D=diag(zeta24^9 b2,1,1)` and `tan(2 pi c2)=-2` | Source normalized; repo unnormalised | `Q(zeta24,sqrt(5))` | Yes, repo exact checks | `{I,F_D,D_bc}` is an exact MUB triple | This is not ordinary F6; it is the twisted Fourier basis in the D0 representative |
| `K6^(3)(theta,phi,lambda)` | Karlsson Theorem 11; Eqs. (2.4)-(2.5), (5.1) | Unnormalised here | Depends on parameters; algebraic at selected points | Yes for audited original A-block and valid branches | Pair `{I,K}` is automatically MU; third/fourth status depends on pool | Covers H2-reducible CHMs, not all order-six CHMs |
| Diţă slice in repo | Karlsson Dita specialization plus Bengtsson D(x) cross-check | Unnormalised here | `Q(zeta12)` for selected algebraic Karlsson H; larger field for B3 | H verified at lambda 0 and pi/3; no W1 certificate there | Third-MUB and no-fourth claims are pointwise unless separately certified | `lambda=pi/2,3pi/2` numerically match D0; seven T3 points do not |

## Exact formulas

### `F6`

The repository uses

`F6[j,k] = omega6^(j*k)`, `j,k=0,...,5`, with `omega6=exp(-2 pi i/6)`.

The sign convention is equivalent to the opposite Fourier sign by conjugation/permutation; it does not change MUB status.

### `D0`

The unnormalised matrix is

```text
[1  1   1    1    1   1]
[1 -1   i   -i   -i   i]
[1  i  -1    i   -i  -i]
[1 -i   i   -1    i  -i]
[1 -i  -i    i   -1   i]
[1  i  -i   -i    i  -1]
```

This is Bengtsson et al. Eq. (11) with `z=1`, after removing the source normalization.

### `D_bc` and `F_D`

The source block formulas are preserved independently in `src/brierley_weigert_notes.jl`. The repository's exact verification checks `H H^dagger=6I`, unimodularity, and the cross-MUB modulus condition. They are not defined by assigning them equal to the repository's Karlsson candidate.

### Karlsson family

`K6` has block form

```text
[ F2                 Z1                 Z2              ]
[ Z3          1/2 Z3 A Z1       1/2 Z3 B Z2              ]
[ Z4          1/2 Z4 B Z1       1/2 Z4 A Z2              ]
```

with `B=-F2-A`, `Z1,Z2=[[1,1],[z,-z]]`, `Z3,Z4=[[1,z],[1,-z]]`, and

`A11=-1/2 + i sqrt(3)/2 (cos(theta)+exp(-i phi) sin(theta))`,

`A12=-1/2 + i sqrt(3)/2 (-cos(theta)+exp(i phi) sin(theta))`,

`A=[[A11,A12],[conj(A12),-conj(A11)]]`.

Karlsson's Theorem 11 gives the four Mobius relations for the squared `z_i`. The repository's seam and pole handling is an implementation correction/robustness layer; it is not a new family theorem.
