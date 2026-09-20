# Claim Audit

| Claim | Source/evidence | Executed verification? | Scope | Status | Required correction |
|---|---|---:|---|---|---|
| A pair `{I,H}` is represented by an order-six CHM | BW Sec. 2, Eqs. (1)-(3); repo definitions | Yes | Any fixed pair with computational basis chosen | Supported | Keep normalized/unnormalized convention explicit |
| F6 has 48 physical MU vectors and 16 third bases | Grassl Theorem 2; BW Sec. 3.3; repo benchmark | Yes, benchmark PASS | F6/Heisenberg pair | Supported and independently reproduced | Do not generalize to arbitrary CHM |
| D0 has 120 vectors, 10 third bases, no fourth | BW Sec. 4.1, Table 1, Eq. (20); repo benchmark | Yes, benchmark PASS | D0 and equivalent pairs | Supported and independently reproduced | Attribute the theorem to BW; repo T4 is a re-check |
| D_bc is source block-circulant D0 representative | Bengtsson Eqs. (79)-(80); repo matrix code | Yes, exact/HP CHM and equivalence checks | D0 representative | Supported | State row+column equivalence, not column-only |
| F_D is the ordinary F6 | Bengtsson Eq. (61), Eq. (78); repo code | Yes, source formula and exact checks | D0 twisted-product triple | Incorrect if stated | Call it twisted Fourier `F_D`, not ordinary `F6` |
| Karlsson K6^(3) covers all order-six CHMs | Karlsson Theorem 11 and Sec. 7; Wuttig--Tindall v2 Corollary 25 | Yes, source contradicts broad wording | All CHMs | Incorrect as a Karlsson-only claim | Restrict Karlsson to H2-reducible CHMs; record Wuttig--Tindall's newer exhaustive-classification claim separately and mark it not independently reproduced |
| Karlsson A-block used by repo is valid | Karlsson Theorem 11; repo regression tests | Yes | Original transcription | Supported | Keep original formula distinct from literal review transcription |
| Printed review A-block is not unitary | Repo symbolic/numerical audit; review text inspected | Yes in repo; source print inspected | Literal transcription as printed | Repository correction, not a literature theorem | Label as transcription discrepancy; do not imply review family is mathematically invalid |
| Seven Dita points have no fourth MUB | `fourth_mub_obstruction` proof and logs | Yes, certified numerical all-clique W1 | Seven listed points, recovered complete pools, 40/40 cliques | Supported as local certified-numerical | Do not call it a circle or family theorem |
| T3 is new relative to BW | E0 residual log; BW D0 theorem | Yes, but equivalence test is narrower at T3 points | Seven T3 points vs D0 | Defensible only with weaker qualifier | Say not matched by the tested column-only comparison; complete row+column equivalence remains to be checked at every T3 point |
| T4 proves a new D0 nonextendability theorem | Exact M2 log; BW Sec. 4.1 | Yes | One `(D_bc,F_D)` witness | Overstated if called new mathematical result | Call it exact pipeline re-verification; BW already proves stronger D0 result |
| Dense Dita lambda sweep proves no fourth on full circle | `dita_lambda_fourth_dense.csv` | Yes, finite dataset | 628 numerical samples | Unsupported as universal theorem | Report as numerical evidence only |
| 1863/1865 special-locus rows prove no fourth on R | CSV and claim ledger | Yes, finite rows | Sampled special-locus set | Unsupported as region theorem | State sample size and residual open quantifier |
| 512-point grid proves no third/fourth in K6 | CSV; grid excludes theta=0/Dita loci | Yes, sampled | Interior grid only | Incorrect if globalized | Preserve as sampled negative; it misses known loci |
| `dim I=-1` over `CC` is exact emptiness | M2 exports and ALGEBRAIC_ATTACK | The repo documents warning | Inexact coefficient ring | Incorrect | Require exact named number field or a valid numerical exclusion certificate |
| Candidate B3 at lambda=0.4 is algebraic in Q(i,sqrt2,sqrt3) | LLL result in findings | Yes, partial only | 6/36 trivial entries | Unsupported | Treat as inconclusive; BW phase field needs sqrt(5) and zeta24 at D0 |
| Dita third MUB exists for every lambda | Dense/periodicity samples and HP checks | Yes, finite samples | Dita circle samples | Not a theorem | State numerical/HP evidence unless analytic B3 formula is supplied |
| Liang/Chen papers provide a Karlsson-style family | repo assessment plus review cross-check | Partial source check | Cross-family route | Unsupported/closed gap | Do not implement as a parametric family without a primary formula |
| N(6)=3 | Primary sources state open; repo scope says open | No | Global | Open | Never place in theorem/proved ledger |
