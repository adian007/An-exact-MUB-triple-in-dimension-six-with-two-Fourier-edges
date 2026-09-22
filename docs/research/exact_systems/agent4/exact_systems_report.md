# Agent 4: exact algebraic formulations and normalization audit

**Status:** INCOMPLETE. **Computation:** NOT_RUN. This is an equation and source audit only: no solver, Groebner basis, elimination, or sweep was run.

## Blocking normalization finding

For an unnormalized CHM C, C C^dagger = 6 I. With v=z/sqrt(6) and z[0]=1, MU to C[:,k]/sqrt(6) is |C[:,k]^dagger z|^2=6. If U=C/sqrt(6) is unitary, the equivalent equation is |U[:,k]^dagger z|^2=1.

Evidence that numerical third-basis inputs are unit-normalized:

* src/mub_zauner_6d_liang_chen.jl lines 236-238 maps a pool solution to [1;z]/sqrt(6).
* src/Cliques.jl line 61 sets B=hcat(pool[i]...), and lines 62-75 check unit column norms and B'B=I.
* src/Benchmarks.jl line 118 passes hcat([out.pool[i] for i in b.pool_indices]...) to w1_witness_certificate.
* src/Certification.jl lines 92-97 uses that B3 in the inner-product product and subtracts 6.

Thus numerical W1 receives a unitary B3 but uses the unnormalized RHS 6. For a unitary B3 the correct equation is product-1=0. The appropriate correction is either preserve unit-normalized B3 and use RHS 1, or explicitly set T=sqrt(6)*B3 and retain RHS 6. The same mismatch exists in legacy _witness_mu_to_basis at src/mub_zauner_6d_liang_chen.jl lines 661-671.

symbolic_export/w1_D0_exact.m2 lines 1-6 declares an exact D_bc/F_D calculation and uses an unnormalized convention, so its displayed RHS 6 is coherent if its declared matrices are unnormalized CHMs. This report does not validate its matrices or certificate.

### Claim impact and remediation gate

The current numerical W1 system need not equal physical W1(H,B3); its emptiness cannot alone certify physical nonextendability for the unitary third bases produced by the clique code. This requires referee review before reliance on numerical T3-style W1 claims. The audit does not retract an independently checked exact unnormalized ideal certificate.

Do not patch code in this campaign. Before reuse require independent implementation to: (i) report B3 scale; (ii) verify H H^dagger=6I, T T^dagger=6I, and abs2((H^dagger T)[j,k])=6 for unnormalized T; (iii) compare normalized and unnormalized formulas under T=sqrt(6)U; and (iv) pass a positive-control MUB witness in dimensions 2 or 3. A cheap structural test specification: before W1 construction, if max(abs(B3'B3-I)) is small choose RHS 1; if max(abs(B3'B3-6I)) is small choose RHS 6; otherwise reject B3.

## Exact systems

Let K be an explicitly declared exact field containing all entries and coefficient conjugates. Let w_i be formal inverse variables, not complex conjugates, and z_0=w_0=1. For unnormalized C define q(C,k)=(sum_j conjugate(C[j,k]) z_j)(sum_j C[j,k] w_j)-6. The physical locus additionally has w_i=conjugate(z_i). Empty complex variety implies empty physical set; a complex root does not prove a physical MU vector.

1. Pair pool: I_pair(H)=<z_i w_i-1 (i=1..5), q(H,k) (k=1..6)> in K[z1..z5,w1..w5]. Its empty complex variety excludes a physical vector MU to fixed H. Retaining only five H equations needs a separate derivation in the selected algebraic model.

2. Fixed triple: I_W1(H,T)=I_pair(H)+<q(T,k) (k=1..6)>, 17 generators in 10 variables. A unit ideal excludes even one fourth-MU vector, hence a fourth ONB, for that fixed valid triple only. A root is not a basis.

3. Fourth ONB: introduce (z[a,i],w[a,i]), a=1..6 and i=1..5, impose W1 for every a and ordered orthogonality 1+sum_i z[a,i]w[b,i]=0 for a != b. In a complexification include both orientations or use real coordinates. Emptiness excludes a fourth ONB for fixed triple. It is unnecessary when W1 itself is empty.

4. All third bases: use T=(t_jk) and conjugate-copy S=(s_jk), with t_jk s_jk=1, T S^T=6I, and (H^dagger T)[j,k](H^T S)[j,k]=6; impose physical S=conjugate(T). Desired statement: forall T [Third(H,T) => W1(H,T)=empty]. Eliminating vector variables alone cannot establish it. A clique list substitutes only after independent complete-pool proof.

5. Selected orbit: for proved H(p),T(p), use W1 over K(p) only after declaring denominators Delta(p). Generic analysis requires saturation by Delta; every Delta=0 stratum needs direct treatment. A unit-circle physical parameter requires a real constraint. The conclusion concerns selected orbit only.

## Saturation and scope

Record coefficient field and involution, term order, each denominator, saturation factor, exceptional fiber, and exact certificate. Never infer exceptional fibers from a generic calculation. A real-coordinate formulation is larger but represents physical conjugation directly.

Files supplied are systems.json, variable_quantifier_table.csv, and a non-executed Macaulay2 template. They make no mathematical result claim.
