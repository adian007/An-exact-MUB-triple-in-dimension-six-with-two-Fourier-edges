# K3 family-membership and MUB-quartet audit — 2026-09-29

This report records C0–C6 for the three supplied pools. The C2 fits are floating-point equivalence searches, not proofs of family membership. No H2-reducible route is used.

## C0. Source theorem and convention

**Source.** Philippe Jaming, Máté Matolcsi, Péter Móra, Ferenc Szöllősi, and Mihály Weiner, “A generalized Pauli problem and an infinite family of MUB-triplets in dimension 6,” arXiv:0902.0882v2 (26 February 2009), [PDF](https://arxiv.org/pdf/0902.0882), [abstract record](https://arxiv.org/abs/0902.0882). The retrieved PDF has 32 pages; the references below use the article's printed page numbers.

**Theorem 1.4, paraphrased (printed p. 5).** For any member \(F(a,b)\) of the paper's order-six Fourier family, the pair consisting of the standard basis and \(F(a,b)\) cannot be extended to a quartet of mutually unbiased orthonormal bases. This is pair-to-quartet nonextension, not merely pair-to-triplet nonextension. The paper's Section 3 expressly treats all parameter values; periodicity identifies them modulo one, and the equivalences described in its reference [5] reduce the search to the triangle with vertices \((0,0),(1/6,0),(1/6,1/12)\) (printed p. 9).

**Paper family and placement, paraphrased from Eq. (3), printed p. 4.** Its two phase parameters are \(x=e^{2\pi ia}\) and \(y=e^{2\pi ib}\), with \(\omega=e^{2\pi i/3}\); the normalized matrix has the factor \(1/\sqrt6\). In zero-based indexing, \(x\) occupies rows \(\{1,4\}\), columns \(\{1,3,5\}\), and \(y\) rows \(\{2,5\}\), columns \(\{1,3,5\}\). This is the transpose placement of the repository convention: repository \(a\) occupies rows \(\{1,3,5\}\), columns \(\{1,4\}\), and repository \(b\) occupies rows \(\{1,3,5\}\), columns \(\{2,5\}\).

**Section 3 method (printed pp. 9–16).** The paper proves its nonextension statement by discretizing the family parameters and vector phases, then performing an exhaustive computer search with explicit interval error estimates. It reports \(N=180\) parameter intervals and \(N'=19\) phase intervals (printed p. 11). Its bounds for the candidate inner-product modulus are
\[
\left|\ |\langle\widetilde u,\widetilde f_j\rangle|-\frac1{\sqrt6}\right|
\le \frac{5\pi}{6N'}\quad(j=0,2,4),
\qquad
\left|\ |\langle\widetilde u,\widetilde f_j\rangle|-\frac1{\sqrt6}\right|
\le \frac{5\pi}{6N'}+\frac{4\pi}{6N}\quad(j=1,3,5)
\]
(printed p. 12, Eqs. (26)–(27)). An improved endpoint estimate is used when the sum of interval lengths is below \(1/\pi\); the paper checks \(5/N'+4/N<1/\pi\) for its chosen values (printed p. 13). Candidate phases are refined over seven generations (p. 14), and the exhaustive search reports no compatible fourth-basis matrix \(\widetilde D\) (p. 16). The algorithm and documentation are identified by the paper at reference [34]; its p. 10 reports a runtime of about six hours on a 3.2 GHz computer. Reference [5] is I. Bengtsson, W. Bruzda, Å. Ericsson, J.-A. Larsson, W. Tadej, and K. Życzkowski, “Mutually unbiased bases and Hadamard matrices of order six,” *J. Math. Phys.* 48 (2007), no. 5, 052106.

**Assessment.** The requested facts are confirmed by the retrieved source: Theorem 1.4 concerns quartets and all real parameters modulo one; the paper's placement is transpose to the repository's; and the proof is computer-assisted with stated error bounds. The published theorem is cited, not computationally reproduced here.

### C0b. Placement check against the repository fitter

The paper matrix was built from Eq. (3) with unit-modulus entries; the common \(1/\sqrt6\) scale cancels in the dephased fit. The repository `fit_fourier_family` tolerance was \(10^{-8}\). Five pseudorandom parameter pairs used seed 20260929.

| Paper parameters \((a,b)\) | `fit_fourier_family(M)` (repo direct) | `fit_fourier_family(M.T)` (repo transpose call) |
|---|---|---|
| (0.3502483584951709, 0.7847217336600212) | NONE; consistency 5.447e-1, total 6.248e-1 | FIT; consistency 5.551e-16, total 4.041e-16 |
| (0.5143357403269661, 0.08997722942380026) | NONE; consistency 4.100e-1, total 4.708e-1 | FIT; consistency 5.090e-16, total 4.965e-16 |
| (0.000875639493777336, 0.5394009336088772) | NONE; consistency 2.025e-1, total 2.415e-1 | FIT; consistency 5.722e-16, total 5.661e-16 |
| (0.3630734994591833, 0.6561419338727186) | NONE; consistency 1.667e-1, total 1.866e-1 | FIT; consistency 4.710e-16, total 5.467e-16 |
| (0.6857523669661448, 0.287783516965283) | NONE; consistency 2.582e-1, total 2.852e-1 | FIT; consistency 5.118e-16, total 4.965e-16 |
| (0, 0) | FIT; consistency 4.727e-16, total 4.710e-16 | FIT; consistency 4.727e-16, total 4.710e-16 |

Thus, for generic tested parameters, the paper's family is detected through the repository transpose call, not its direct call. The \((0,0)\) Fourier matrix is a special overlap of both orientations.

### C0c. Extension symmetry lemma

Let an extension of \((Id,H)\) be an ordered quartet of orthonormal bases \((Id,H,C,D)\), with all six pairs mutually unbiased.

1. Apply the same ambient unitary \(H^\dagger\) to all four bases. Their matrices become \((H^\dagger,Id,H^\dagger C,H^\dagger D)\). A common unitary preserves every inner product, hence preserves orthonormality and mutual unbiasedness. Reorder the first two bases to obtain \((Id,H^\dagger,H^\dagger C,H^\dagger D)\). The operation is reversible, so \((Id,H)\) extends iff \((Id,H^\dagger)\) extends.
2. Complex-conjugate every basis matrix in an extension of \((Id,H)\), obtaining \((Id,\overline H,\overline C,\overline D)\). Conjugation of all matrices preserves absolute inner products. Apply the same ambient unitary \((\overline H)^\dagger=H^T\) to all four; after reordering the first two bases this gives \((Id,H^T,H^T\overline C,H^T\overline D)\). This is reversible, so \((Id,H)\) extends iff \((Id,H^T)\) extends.

These steps use the same ambient unitary on every basis and, in the second step, the same complex conjugation on every basis. They do not independently conjugate or rephase individual bases. Column phases and column permutations merely change vector representatives/order within a basis and preserve the MUB property. The paper's equivalence convention on printed p. 3 likewise requires the conjugation choice to be the same for all matrices. The lemma therefore extends Theorem 1.4 to either orientation of the family.

## C1. Random valid 2-circulant controls

An independently random triple \((\beta,\gamma,\epsilon)\) almost surely does not admit a unit-modulus \(\phi\). To meet the approved constrained-sampling interpretation, \(\beta,\gamma\) were sampled randomly and \(\epsilon,\phi\) were solved on the unit circle. The polynomial used was
\[
\beta\gamma\epsilon\,\phi^2+
(\beta\gamma+\beta^2\epsilon+\gamma\epsilon+\beta\gamma^2\epsilon)\phi+
\beta\gamma\epsilon^2=0.
\]
Seed 20260930 was used; 12 random \((\beta,\gamma)\) attempts yielded five valid, non-Fourier controls. Each received independent random row/column permutations and row/column phase multiplications, using equivalent-matrix seeds 41000–41004. Both polynomial roots and their constraint residuals are listed for every accepted member.

| # | \(\beta\) | \(\gamma\) | \(\epsilon\) | All \(\phi\) roots; \((|\phi|,\ |\mathrm{constraint}|)\) | Rephasing seed | \(X_6\) matrix / constraint residual | Repo F direct / transpose |
|---:|---|---|---|---|---:|---|---|
| 1 | -0.5502842927031812-0.8349773632884664i | 0.7602032785988252-0.6496852893575450i | -0.4042120144643915-0.9146653198643964i | 0.871259934878+0.490821888139i (1, 4.965e-16); -0.801111492710-0.598515142873i (1, 2.289e-16) | 41000 | FIT; 6.939e-16 / 7.772e-16 | NONE / NONE |
| 2 | 0.0154992923721858+0.9998798787534238i | -0.3859754440204046+0.9225090550304921i | 0.9628513804537560-0.2700318854474349i | -0.842732259765-0.538332925197i (1, 1.343e-15); -0.666058864878+0.745899181202i (1, 2.238e-16) | 41001 | FIT; 6.753e-16 / 8.327e-16 | NONE / NONE |
| 3 | -0.3410304889384379-0.9400522355775821i | 0.9939923452253794-0.1094496123033345i | -0.7359545235576676+0.6770309736304585i | 0.996580240943+0.082630644212i (1, 3.331e-16); -0.677494230907+0.735528087218i (1, 3.377e-16) | 41002 | FIT; 3.775e-16 / 1.110e-16 | NONE / NONE |
| 4 | 0.9925663840851134-0.1217044501413282i | -0.5246209436469070+0.8513359298696542i | -0.4787664605848132-0.8779422966340615i | -0.355663477458-0.934614086563i (1, 4.003e-16); 0.990816981886-0.135209868008i (1, 4.003e-16) | 41003 | FIT; 2.168e-15 / 1.742e-15 | NONE / NONE |
| 5 | 0.9573093341953431-0.2890654574010336i | -0.7520843621803215+0.6590668495409392i | 0.1699929234381168+0.9854452830984389i | -0.335844151143+0.941917568656i (1, 1.066e-15); 0.871117096028-0.491075355733i (1, 1.570e-16) | 41004 | FIT; 5.558e-16 / 4.003e-16 | NONE / NONE |

All five fits used tolerance \(10^{-8}\), the \(X_6\) numerical solver, and three solver starts after exhaustive permutation-chart search. All \(\phi\) roots satisfy the requested unit-modulus bound \(10^{-12}\); all reported constraint residuals are below \(1.8\times10^{-15}\). These are NUMERICAL controls.

## C2. All-clique fits and independent reconstruction

Input pools were the three existing 72-by-6 files under `mub_solution_attack_package_v2/01_results/`. Cliques were independently reconstructed with the repository graph routine. The number of six-vector cliques was 4 at each of \(10^{-9},10^{-7},10^{-5}\) for every pool.

| \(\lambda\) | Pool rows | Cliques at \(10^{-9}\) | at \(10^{-7}\) | at \(10^{-5}\) |
|---|---:|---:|---:|---:|
| 0.400000000000 | 72 | 4 | 4 | 4 |
| 1.047197551197 (\(\pi/3\)) | 72 | 4 | 4 | 4 |
| 2.094395102393 (\(2\pi/3\)) | 72 | 4 | 4 | 4 |

The clique labels used in the fit tables correspond to these pool-row indices:

| \(\lambda\) | c0 | c1 | c2 | c3 |
|---|---|---|---|---|
| 0.400000000000 | (2,24,35,47,66,67) | (3,10,25,29,45,70) | (32,57,60,64,68,69) | (33,39,42,49,52,61) |
| 1.047197551197 | (2,22,23,61,65,71) | (10,24,41,42,51,56) | (11,19,28,48,64,69) | (31,40,46,58,62,67) |
| 2.094395102393 | (1,2,19,44,48,69) | (3,4,14,15,51,61) | (10,27,38,53,57,71) | (24,35,47,59,64,70) |

Every original transition matrix was checked in both orientations for flatness, unitarity, the three-distinct-columns predicate, direct \(F_6^{(2)}(a,b)\), transpose-call \(F_6^{(2)}(a,b)\), direct \(X_6\), and transpose \(X_6\). Fourier fitting tolerance was \(10^{-8}\). All F columns below report consistency error \(e\), total maximum residual \(t\), fixed-entry maximum residual \(f\), zero-based row/column permutation \(p/q\), and fitted unit phases \(a,b\). Permutations are compact digit strings: for example, `201345` means \((2,0,1,3,4,5)\). `F^T-call` means the same direct fitter called on \(M^T\). `X/X^T` give status and (matrix, constraint) maximum residuals. Fits and input matrix checks are NUMERICAL.

### \(\lambda=0.400000000000\)

| Clique | Transition | Flat / unitary error | Theorem predicate \(M/M^T\), min residual | F direct: status; \(e,t,f;p/q;a,b\) | F transpose-call: status; \(e,t,f;p/q;a,b\) | X / X^T: status; matrix / constraint residual |
|---:|---|---|---|---|---|---|
| 0 | B3 | 2.22e-16 / 9.02e-16 | no/yes; 1.931e-1 / 3.539e-16 | FIT; 1.59e-15,1.79e-15,1.79e-15; 201345/530241; .981361+.192171i,.544820+.838553i | NONE; .6941,.7800,.7800; 410253/502143; .577772+.816198i,.417963-.908464i | NONE 1.931e-1/9.434e-2; NONE 1.931e-1/9.434e-2 |
| 0 | \(H_D^\dagger B3\) | 1.33e-15 / 8.49e-16 | no/yes; 4.028e-1 / 8.327e-16 | FIT; 1.65e-15,1.60e-15,1.52e-15; 214503/503214; .809071+.587711i,.801117-.598508i | NONE; .4139,.4159,.4159; 201435/021534; .501823+.864970i,.498174-.867077i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 0 | \(B3^\dagger H_D\) | 1.33e-15 / 8.92e-16 | yes/no; 8.368e-16 / 4.028e-1 | NONE; .4139,.4159,.4159; 201435/120435; -.501823-.864970i,.498174-.867077i | FIT; 1.75e-15,1.73e-15,1.67e-15; 521034/204513; .913508-.406820i,-.918882-.394533i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 1 | B3 | 2.22e-16 / 4.44e-16 | no/yes; 1.931e-1 / 3.123e-16 | FIT; 6.75e-16,7.11e-16,7.11e-16; 530412/043521; .274095-.961703i,.324255-.945970i | NONE; .6941,.7800,.7800; 031425/321054; .577772+.816198i,.417963-.908464i | NONE 1.931e-1/9.434e-2; NONE 1.931e-1/9.434e-2 |
| 1 | \(H_D^\dagger B3\) | 6.66e-16 / 6.66e-16 | no/yes; 4.028e-1 / 3.724e-16 | FIT; 9.99e-16,9.55e-16,9.16e-16; 413250/521043; .296410+.955061i,.801117+.598508i | NONE; .4139,.4159,.4159; 304152/521034; .999998-.002107i,-.999998+.002107i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 1 | \(B3^\dagger H_D\) | 6.66e-16 / 6.66e-16 | yes/no; 4.678e-16 / 4.028e-1 | NONE; .4139,.4159,.4159; 130524/021534; .498174-.867077i,-.501823-.864970i | FIT; 8.95e-16,8.90e-16,8.52e-16; 041325/203451; -.296410+.955061i,.809071-.587711i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 2 | B3 | 2.22e-16 / 1.11e-15 | no/yes; 1.931e-1 / 9.310e-16 | FIT; 2.54e-15,2.45e-15,2.31e-15; 150342/043125; .998618+.052552i,.324255-.945970i | NONE; .6941,.7800,.7800; 301425/410325; .577772+.816198i,-.417963+.908464i | NONE 1.931e-1/9.434e-2; NONE 1.931e-1/9.434e-2 |
| 2 | \(H_D^\dagger B3\) | 2.44e-15 / 1.06e-15 | no/yes; 4.028e-1 / 8.693e-16 | FIT; 2.99e-15,3.11e-15,3.11e-15; 415023/302514; -.801117+.598508i,.296410-.955061i | NONE; .4139,.4159,.4159; 014253/102453; .498174-.867077i,.501823+.864970i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 2 | \(B3^\dagger H_D\) | 2.44e-15 / 1.11e-15 | yes/no; 1.013e-15 / 4.028e-1 | NONE; .4139,.4159,.4159; 102435/512043; -.498174+.867077i,-.501823-.864970i | FIT; 3.39e-15,3.55e-15,3.55e-15; 402153/304512; .918882-.394533i,-.678902+.734229i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 3 | B3 | 2.22e-16 / 1.15e-15 | no/yes; 1.931e-1 / 3.002e-16 | FIT; 1.41e-15,1.76e-15,1.76e-15; 310245/103524; .657106+.753798i,-.969906+.243478i | NONE; .6941,.7800,.7800; 120435/502314; -.577772-.816198i,-.417963+.908464i | NONE 1.931e-1/9.434e-2; NONE 1.931e-1/9.434e-2 |
| 3 | \(H_D^\dagger B3\) | 6.66e-16 / 1.07e-15 | no/yes; 4.028e-1 / 5.150e-16 | FIT; 1.99e-15,2.19e-15,2.19e-15; 105324/235041; -.809071-.587711i,-.801117+.598508i | NONE; .4139,.4159,.4159; 204351/521034; .498174-.867077i,-.501823-.864970i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |
| 3 | \(B3^\dagger H_D\) | 6.66e-16 / 8.89e-16 | yes/no; 5.332e-16 / 4.028e-1 | NONE; .4139,.4159,.4159; 021534/021534; -.498174+.867077i,.501823+.864970i | FIT; 1.64e-15,1.91e-15,1.91e-15; 410532/253014; .801117-.598508i,.809071+.587711i | NONE 4.028e-1/2.165e-1; NONE 4.028e-1/2.165e-1 |

### \(\lambda=1.047197551197\) (\(\pi/3\))

| Clique | Transition | Flat / unitary error | Theorem predicate \(M/M^T\), min residual | F direct: status; \(e,t,f;p/q;a,b\) | F transpose-call: status; \(e,t,f;p/q;a,b\) | X / X^T: status; matrix / constraint residual |
|---:|---|---|---|---|---|---|
| 0 | B3 | 2.22e-16 / 1.76e-15 | no/yes; 3.041e-1 / 1.862e-16 | FIT; 3.60e-15,3.58e-15,3.18e-15; 431025/204315; .521141+.853470i,-.753539-.657403i | NONE; .643,.702,.702; 501324/310254; -.444978-.895541i,-.553073+.833133i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 0 | \(H_D^\dagger B3\) | 2.89e-15 / 1.81e-15 | no/yes; 2.887e-1 / 6.062e-16 | FIT; 3.33e-15,3.69e-15,3.69e-15; 103452/124035; -.959701-.281022i,-.890371+.455235i | NONE; .439,.468,.468; 013245/420135; .525557-.850758i,-.474000-.880525i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 0 | \(B3^\dagger H_D\) | 2.89e-15 / 1.55e-15 | yes/no; 6.713e-16 / 2.887e-1 | NONE; .439,.468,.468; 013245/402153; .474000-.880525i,-.525557-.850758i | FIT; 3.13e-15,3.14e-15,3.04e-15; 305412/421530; -.958329-.285667i,.050941-.998702i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 1 | B3 | 2.22e-16 / 1.51e-15 | no/yes; 3.041e-1 / 5.862e-16 | FIT; 2.19e-15,2.34e-15,2.34e-15; 521403/140352; -.216626+.976255i,-.999698-.024586i | NONE; .643,.702,.702; 520314/102435; .553073+.833133i,.444978-.895541i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 1 | \(H_D^\dagger B3\) | 6.66e-16 / 1.48e-15 | no/yes; 2.887e-1 / 2.355e-16 | FIT; 2.34e-15,2.36e-15,2.36e-15; 241503/015234; .723223+.690615i,.231769+.972771i | NONE; .439,.468,.468; 120453/201354; .525557+.850758i,-.474000+.880525i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 1 | \(B3^\dagger H_D\) | 6.66e-16 / 9.21e-16 | yes/no; 3.331e-16 / 2.887e-1 | NONE; .439,.468,.468; 302541/521034; -.525557-.850758i,.474000-.880525i | FIT; 2.42e-15,2.42e-15,2.40e-15; 250413/350142; -.050941-.998702i,.236479-.971637i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 2 | B3 | 2.22e-16 / 1.43e-15 | no/yes; 3.041e-1 / 7.650e-16 | FIT; 2.88e-15,2.88e-15,1.94e-15; 543120/150243; .737148+.675731i,.478557+.878057i | NONE; .643,.702,.702; 520413/531042; -.998051+.062408i,-.998051+.062408i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 2 | \(H_D^\dagger B3\) | 2.00e-15 / 1.40e-15 | no/yes; 2.887e-1 / 6.179e-16 | FIT; 2.14e-15,2.69e-15,2.69e-15; 402153/410523; .890371+.455235i,-.726559-.687104i | NONE; .439,.468,.468; 203541/420135; .525557+.850758i,.474000-.880525i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 2 | \(B3^\dagger H_D\) | 2.00e-15 / 1.55e-15 | yes/no; 6.661e-16 / 2.887e-1 | NONE; .439,.468,.468; 120354/021534; -.474000+.880525i,-.525557-.850758i | FIT; 2.13e-15,2.11e-15,1.98e-15; 405321/104235; .723223+.690615i,.839431-.543467i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 3 | B3 | 2.22e-16 / 7.80e-16 | no/yes; 3.041e-1 / 7.448e-16 | FIT; 2.04e-15,2.04e-15,1.87e-15; 420135/324015; .478557+.878057i,-.737148-.675731i | NONE; .643,.702,.702; 401532/230451; .998051+.062408i,-.998051-.062408i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 3 | \(H_D^\dagger B3\) | 1.78e-15 / 7.67e-16 | no/yes; 2.887e-1 / 8.834e-16 | FIT; 2.41e-15,2.48e-15,2.48e-15; 201354/041352; .890371+.455235i,-.959701+.281022i | NONE; .439,.468,.468; 301245/201354; -.525557+.850758i,-.474000-.880525i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |
| 3 | \(B3^\dagger H_D\) | 1.55e-15 / 8.98e-16 | yes/no; 8.834e-16 / 2.887e-1 | NONE; .439,.468,.468; 301245/210345; .474000-.880525i,.525557+.850758i | FIT; 2.36e-15,2.38e-15,2.38e-15; 105324/014325; -.236479-.971637i,.050941-.998702i | NONE 2.895e-1/1.908e-1; NONE 2.895e-1/1.908e-1 |

### \(\lambda=2.094395102393\) (\(2\pi/3\))

| Clique | Transition | Flat / unitary error | Theorem predicate \(M/M^T\), min residual | F direct: status; \(e,t,f;p/q;a,b\) | F transpose-call: status; \(e,t,f;p/q;a,b\) | X / X^T: status; matrix / constraint residual |
|---:|---|---|---|---|---|---|
| 0 | B3 | 2.22e-16 / 2.44e-15 | no/yes; 2.887e-1 / 6.713e-16 | FIT; 4.44e-15,4.87e-15,4.87e-15; 102453/204315; .236479+.971637i,-.050941+.998702i | NONE; .439,.468,.468; 210543/510432; .525557+.850758i,.474000-.880525i | NONE 2.887e-1/1.908e-1; NONE 2.887e-1/1.908e-1 |
| 0 | \(H_D^\dagger B3\) | 4.88e-15 / 2.41e-15 | no/yes; 3.041e-1 / 3.724e-16 | FIT; 4.68e-15,4.68e-15,4.07e-15; 024531/241350; .216626+.976255i,.999698-.024586i | NONE; .643,.702,.702; 410325/302514; .444978+.895541i,.553073-.833133i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 0 | \(B3^\dagger H_D\) | 4.88e-15 / 2.00e-15 | yes/no; 7.114e-16 / 3.041e-1 | NONE; .643,.702,.702; 410325/320541; -.553073-.833133i,-.444978+.895541i | FIT; 4.72e-15,4.90e-15,4.90e-15; 013542/240351; .216626-.976255i,-.999698-.024586i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 1 | B3 | 2.22e-16 / 1.12e-15 | no/yes; 2.887e-1 / 6.004e-16 | FIT; 2.23e-15,2.38e-15,2.38e-15; 302451/152043; .231769-.972771i,.723223-.690615i | NONE; .439,.468,.468; 210453/104325; .525557+.850758i,-.474000+.880525i | NONE 2.887e-1/1.908e-1; NONE 2.887e-1/1.908e-1 |
| 1 | \(H_D^\dagger B3\) | 8.88e-16 / 1.08e-15 | no/yes; 3.041e-1 / 6.713e-16 | FIT; 2.16e-15,2.13e-15,2.12e-15; 405321/305214; .521141+.853470i,-.953774-.300524i | NONE; .643,.702,.702; 510324/301425; .553073+.833133i,-.444978+.895541i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 1 | \(B3^\dagger H_D\) | 8.88e-16 / 1.02e-15 | yes/no; 6.866e-16 / 3.041e-1 | NONE; .643,.702,.702; 401235/103524; .444978+.895541i,-.553073+.895541i | FIT; 2.04e-15,2.04e-15,2.00e-15; 415023/124035; .478557-.878057i,-.946097+.323883i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 2 | B3 | 2.22e-16 / 7.47e-16 | no/yes; 2.887e-1 / 2.498e-16 | FIT; 1.12e-15,1.25e-15,1.25e-15; 125304/325104; .959701+.281022i,-.726559-.687104i | NONE; .439,.468,.468; 013542/241053; .525557+.850758i,.474000-.880525i | NONE 2.887e-1/1.908e-1; NONE 2.887e-1/1.908e-1 |
| 2 | \(H_D^\dagger B3\) | 6.66e-16 / 7.63e-16 | no/yes; 3.041e-1 / 3.331e-16 | FIT; 1.18e-15,1.42e-15,1.42e-15; 230514/104325; .521141-.853470i,-.753539+.657403i | NONE; .643,.702,.702; 501423/501243; .998051-.062408i,-.998051+.062408i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 2 | \(B3^\dagger H_D\) | 6.66e-16 / 8.88e-16 | yes/no; 3.724e-16 / 3.041e-1 | NONE; .643,.702,.702; 501423/510234; .998051+.062408i,-.998051-.062408i | FIT; 1.36e-15,1.39e-15,1.29e-15; 130425/051243; -.216626+.976255i,.999698+.024586i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 3 | B3 | 2.22e-16 / 5.29e-16 | no/yes; 2.887e-1 / 4.710e-16 | FIT; 1.02e-15,1.14e-15,1.14e-15; 015234/102354; -.726559+.687104i,.959701-.281022i | NONE; .439,.468,.468; 403251/510432; .999557-.029766i,.999557-.029766i | NONE 2.887e-1/1.908e-1; NONE 2.887e-1/1.908e-1 |
| 3 | \(H_D^\dagger B3\) | 6.66e-16 / 6.66e-16 | no/yes; 3.041e-1 / 6.004e-16 | FIT; 1.19e-15,1.20e-15,1.20e-15; 032145/152304; -.946097+.323883i,.478557-.878057i | NONE; .643,.702,.702; 503142/204135; .444978-.895541i,-.553073-.833133i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |
| 3 | \(B3^\dagger H_D\) | 6.66e-16 / 8.88e-16 | yes/no; 5.266e-16 / 3.041e-1 | NONE; .643,.702,.702; 031425/501432; -.553073-.833133i,-.444978+.895541i | FIT; 1.09e-15,1.09e-15,1.05e-15; 052341/120345; -.521141+.853470i,.753539-.657403i | NONE 3.041e-1/1.439e-1; NONE 3.041e-1/1.439e-1 |

### C3 classification and anomaly check

| Transition | Direct F | F by transpose call | X6 in either orientation | Classification |
|---|---:|---:|---:|---|
| B3 | FIT (all 12) | NONE | NONE (all 12) | F |
| \(H_D^\dagger B3\) | FIT (all 12) | NONE | NONE (all 12) | F |
| \(B3^\dagger H_D\) | NONE (all 12) | FIT (all 12) | NONE (all 12) | \(F^T\) |

No `several` or `none` classifications occurred. Every matrix passed the theorem predicate in one orientation and had an F/F-transposed numerical fit in an orientation; there were zero C3 anomalies. These classifications apply only to the finite matrix samples.

## C4. A4 action on the four cliques

The two generators \(G_2,G_3\) and their monomial automorphism equations are recorded in `mub_fold_a4_computation/source_context/01_symbolic_I3/fixed_symmetry_exact.txt`; the package reports exact symbolic verification. Their generated action was closed projectively to 12 elements at each tested \(\lambda\). Each resulting monomial \(U\) was checked numerically by forming \(Q=(H_D^\dagger U H_D)^\dagger\) and verifying \(U H_D Q=H_D\). Maximum equation residuals were \(2.777\times10^{-16}\), \(2.227\times10^{-16}\), and \(1.755\times10^{-16}\), respectively.

Applying all 12 maps to clique 0 and matching six vectors up to phase gave these image counts:

| \(\lambda\) | Images of clique 0 among \((c0,c1,c2,c3)\) | Orbit size | Stabilizer size |
|---|---|---:|---:|
| 0.400000000000 | (3,3,3,3) | 4 | 3 |
| 1.047197551197 | (3,3,3,3) | 4 | 3 |
| 2.094395102393 | (3,3,3,3) | 4 | 3 |

All matched vector sets had projective distance zero to floating precision (maximum overlap error at most \(2.3\times10^{-16}\)). In the fixed closure order \(e,e2,e3,e23,e32,e33,e232,e233,e323,e332,e3233,e3323\), images were:

| \(\lambda\) | Image labels in that order |
|---|---|
| 0.400000000000 | c0,c2,c1,c3,c2,c3,c0,c1,c0,c2,c3,c1 |
| 1.047197551197 | c0,c2,c1,c3,c2,c3,c0,c1,c0,c2,c3,c1 |
| 2.094395102393 | c0,c2,c3,c1,c0,c2,c2,c0,c1,c3,c3,c1 |

The B3 direct-fit parameter representatives \((a,b)\), one per clique, are:

| \(\lambda\) | c0 | c1 | c2 | c3 |
|---|---|---|---|---|
| 0.400000000000 | (.981361+.192171i, .544820+.838553i) | (.274095-.961703i, .324255-.945970i) | (.998618+.052552i, .324255-.945970i) | (.657106+.753798i, -.969906+.243478i) |
| 1.047197551197 | (.521141+.853470i, -.753539-.657403i) | (-.216626+.976255i, -.999698-.024586i) | (.737148+.675731i, .478557+.878057i) | (.478557+.878057i, -.737148-.675731i) |
| 2.094395102393 | (.236479+.971637i, -.050941+.998702i) | (.231769-.972771i, .723223-.690615i) | (.959701+.281022i, -.726559-.687104i) | (-.726559+.687104i, .959701-.281022i) |

The symmetry predicts monomial equivalence of the four transition matrices, and that is observed: all images are F fits and the vector sets map up to phase. The numerical \((a,b)\) outputs are chart-dependent representatives, not invariants under row/column permutation and dephasing. Applying the 12 group elements and refitting confirms F membership, but the returned coordinates are not consistently equal to the target clique's canonical fit coordinates. No explicit parameter-coordinate action law is asserted; only the matrix-level orbit/equivalence statement is NUMERICAL.

## C5. Bridge audit

1. **Pair \((Id,B3)\).** For each of the four cliques at each \(\lambda\), C2's `B3` row gives a direct repository-F fit, including its row permutation \(p\), column permutation \(q\), and fitted \((a,b)\). The paper's convention is the transpose placement (C0b), and C0c extends Theorem 1.4 to that orientation. If the fitted equivalence were exact, the published theorem would exclude a quartet containing \((Id,B3)\). Any quartet \(\{I,H_D,B3,B4\}\) contains that pair, so the same conditional exclusion would apply.
2. **Separate transformed-pair corroboration.** Left-multiply every basis in a hypothetical quartet \(\{I,H_D,B3,B4\}\) by \(H_D^\dagger\), then reorder bases. This gives a quartet containing \((I,T)\), \(T=H_D^\dagger B3\). All 12 such \(T\) matrices have direct repository-F fits in C2; the exact \(p,q,(a,b)\) values are in the corresponding C2 rows. C0b/C0c again match the paper's orientation and cover its transpose. If exact family membership held, Theorem 1.4 would independently exclude the transformed quartet.
3. **The pair \((I,H_D)\) is not the exclusion pair.** A separate direct \(X_6\) fit of \(H_D\) returned FIT at all three samples, with matrix/constraint residuals \(5.551\times10^{-17}/5.551\times10^{-17}\), \(2.483\times10^{-16}/2.289\times10^{-16}\), and \(2.776\times10^{-16}/4.965\times10^{-16}\), respectively. This is a NUMERICAL check of the 2-circulant description. The argument above uses the Fourier pair \((I,B3)\), and separately \((I,H_D^\dagger B3)\), not \((I,H_D)\).

**Gaps and scope.** The transition matrices are normalized basis matrices, while the paper displays the corresponding Hadamard matrices with a \(1/\sqrt6\) scale; the scalar normalization cancels in the dephased fit. The reported \(p,q\) encode row/column permutations, and dephasing removes row/column phases; the fit residual is numerical, so the equivalence is not an exact identity certificate. Column ordering and phases are basis relabellings. Theorem 1.4 covers all real parameters modulo one, and C0c covers both transpose orientations, so there is no parameter-range/orientation gap once exact equivalence is known. The decisive remaining gap is exact membership: residuals around \(10^{-15}\) do not prove exact membership. This is a NUMERICAL application of a published computer-assisted theorem, not an exact exclusion of a quartet.

## C6. Addenda, tiers, and hashes

The package's older pi/3 sample record is superseded by the new numerical fit; see the new [package addendum](../../mub_solution_attack_package_v2/ADDENDUM_pi3_fit.md). No pre-existing package file was changed.

| Evidence input | SHA-256 |
|---|---|
| `mub_solution_attack_package_v2/01_results/pool_0.400000000000.npz` | `30411a5e1e5f8a32bbf354d137bdfe01b9ea2f066b80aa818b42f7685f376da1` |
| `mub_solution_attack_package_v2/01_results/pool_1.047197551197.npz` | `20adf0957e87afdf50259fb29ed74c464fd5bfe851612fbc4f42bdeffbf43262` |
| `mub_solution_attack_package_v2/01_results/pool_2.094395102393.npz` | `ae91d70ada404da588d73843b2e22ef54092b981f1716c7f6766a49000d0613` |
| `scripts/python/k3_family_membership.py` | `df433144a81045b2628a2e5edd0d2f734be9db450fbc7e5e5be7d080caec9a5a` |
| `scripts/python/analyze_k3_fourier_structure.py` | `1bdb7aa3e68e5cceab5f1b2f2ce3f5557dcad3b6f9052ae85a2accaa4baa7827` |
| `mub_fold_a4_computation/source_context/01_symbolic_I3/fixed_symmetry_exact.txt` | `002723b31feab92db3c898249e52ad26fba4df28cfda0eff9fc9356f1939a221` |
| untouched `mub_solution_attack_package_v2/01_results/Fourier_transition_verified_samples.csv` | `6d820ea5773a01797a2d499d3925604c07c454b55051f94b36786903dda13e52` |
| arXiv source PDF `0902.0882v2` | `f06106a71ddeb5cb4a4bd96d082334a189e1a4fe15b0d3ea4bdcd672dfa37b7b` |
| untouched legacy `diagnostic.json` | `6e1089ee82cd1e444a28ef05e8e4653a58d2ad60dd259ebe382b6db9a63d4924` |
| untouched legacy `exact_chart_equations.json` | `fb7e7f592c6c2adab9533b66ba8b788d5f4d5df464597d202bc0f003c32e4c32` |

The claim ledger rows for these findings cite the SHA-256 of this report as their evidence-bundle hash. The report does not contain its own hash. Runs used Python 3.13.14, NumPy 2.4.6, SciPy 1.17.1, and psutil 7.2.2. Pool reconstruction used `six_cliques(V,tol)` at `tol in (1e-9,1e-7,1e-5)`; family fits used `fit_fourier_family(M)`, `fit_fourier_family(M.T)`, `fit_two_circulant(M)`, and `fit_two_circulant(M.T)` sequentially at tolerance `1e-8`.

## What is NOT established

- The results cover only the three tested \(\lambda\) values.
- Pool completeness is not addressed.
- Numerical fits are not family-level or component-level membership results.
- Nothing here addresses the fold or any \(\lambda<\lambda_*\).
- Theorem 1.4 is cited from the published paper and summarized, not independently verified or computationally reproduced here.
- Consequently, the bridge does not establish exact nonextendibility for any approximate pool clique, nor does it address global MUB nonexistence.
