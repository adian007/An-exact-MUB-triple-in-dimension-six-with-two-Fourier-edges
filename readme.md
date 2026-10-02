# An exact MUB triple in dimension six with two Fourier edges

Research code, symbolic computations, numerical experiments, and reproducibility
artifacts for the study of mutually unbiased bases (MUBs) in dimension six.

The central research object is an explicit exact MUB triple

\[
\left\{
I,\frac{D(\rho)}{\sqrt{6}},B
\right\},
\qquad
\rho=e^{i\pi/3},
\]

where \(D(\rho)\) belongs to the Diţă slice of Karlsson's three-parameter
family of \(6\times6\) complex Hadamard matrices.

The main structural result is that the two transition matrices involving the
constructed third basis \(B\) are, up to the relevant equivalences and
transposition, members of the two-parameter Fourier family. A known
non-extendability result for Fourier-family pairs then implies that this
particular MUB triple cannot be extended to a fourth basis.

This repository contains the computational and symbolic material supporting
the paper, together with research notes, protocols, provenance records,
reproduction instructions, and intermediate results.

---

## 1. Research problem

A set of orthonormal bases

\[
\mathcal B_1,\mathcal B_2,\ldots,\mathcal B_k
\]

in \(\mathbb C^d\) is called mutually unbiased when, for every pair of
distinct bases,

\[
|\langle u_i,v_j\rangle|^2=\frac{1}{d}.
\]

For \(d=6\), this becomes

\[
|\langle u_i,v_j\rangle|=\frac{1}{\sqrt 6}.
\]

The dimension-six case is particularly interesting because the structure and
existence of complete sets of MUBs remain an important open problem.

This project studies the problem through the connection between MUBs and
\(6\times6\) complex Hadamard matrices.

After fixing the computational basis \(I\), an unbiased basis can be
represented by a complex Hadamard matrix \(H\) satisfying

\[
HH^\dagger=6I
\]

and

\[
|H_{jk}|=1.
\]

The project focuses on Karlsson's three-parameter family and, in particular,
on its one-parameter Diţă slice.

---

## 2. Main result

The main exact construction studied in this repository is

\[
\left\{
I,
\frac{D(\rho)}{\sqrt6},
B
\right\},
\qquad
\rho=e^{i\pi/3}.
\]

Here:

- \(I\) is the computational basis;
- \(D(\rho)\) is a complex Hadamard matrix on the Diţă slice;
- \(B\) is an explicitly constructed third basis;
- all three bases are mutually unbiased.

The key structural observation is that the two transition matrices connecting
the third basis \(B\) to the other two bases belong, up to equivalence and
transposition, to the two-parameter Fourier family.

Consequently, a known theorem concerning the non-extendability of
Fourier-family pairs applies to this triple.

Therefore, the constructed triple is not contained in a set of four MUBs.

### Important scope of the result

This repository does **not** claim to prove that four MUBs cannot exist in
dimension six.

The result concerns the specific exact MUB triple constructed here.

In particular, the following statements should not be conflated:

- the constructed triple is non-extendable;
- every possible MUB triple in dimension six is non-extendable;
- four MUBs do not exist in dimension six.

Only the first statement is established by the exact construction and the
Fourier-family obstruction used in this project.

---

## 3. The Diţă slice

The project uses the following one-parameter slice of the Karlsson family:

\[
D(z)=
\begin{pmatrix}
1&1&1&1&1&1\\
1&-1&z&-z&i&-i\\
1&-i&i&i&-i&-1\\
1&i&-z&z&-1&-i\\
1&\bar z&-i&-1&-z&i\\
1&-\bar z&-1&-i&z&i
\end{pmatrix},
\qquad |z|=1.
\]

The parameter is often written as

\[
z=e^{i\lambda}.
\]

The main exact point studied in the paper is

\[
z=\rho=e^{i\pi/3}.
\]

A second exact construction is also studied at

\[
z=1.
\]

The computational study additionally examines several parameter values across
the Diţă slice.

---

## 4. Why the Fourier family matters

The two transition matrices associated with the constructed third basis have
a special structure: after the appropriate equivalences and transposition,
they belong to the two-parameter Fourier family.

This is important because Fourier-family pairs have a known obstruction to
extension to a fourth mutually unbiased basis.

The exact argument therefore has the following structure:

```text
Exact MUB triple
      |
      v
Two transition matrices
      |
      v
Fourier-family structure
      |
      v
Known Fourier-family obstruction
      |
      v
This particular triple is non-extendable
