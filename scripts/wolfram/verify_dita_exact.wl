(* Exact Dita-slice verification at lambda=0 and lambda=Pi/3.
   The algebraic branches are the ones documented by
   scripts/julia/exact_karlsson_dita_H.jl and docs/ALGEBRAIC_ATTACK.md.
*)

ClearAll["Global`*"];

zeta12 = RootReduce[Exp[I Pi/6]];
ii = zeta12^3;
sqrt3 = zeta12 + zeta12^11;

conj12[x_] := RootReduce[Conjugate[x]];
mobius[z_, alpha_, beta_] := RootReduce[(alpha z - beta)/(conj12[beta] z - conj12[alpha])];

f2 = {{1, 1}, {1, -1}};
a = {{ii, -1}, {-1, ii}};
b = RootReduce[-f2 - a];

zLeft[z_] := {{1, 1}, {z, -z}};
zRight[z_] := {{1, z}, {1, -z}};

buildH[z1_, z2_, z3_, z4_] := ArrayFlatten[{
    {f2, zLeft[z1], zLeft[z2]},
    {zRight[z3], (1/2) zRight[z3].a.zLeft[z1], (1/2) zRight[z3].b.zLeft[z2]},
    {zRight[z4], (1/2) zRight[z4].b.zLeft[z1], (1/2) zRight[z4].a.zLeft[z2]}
    }];

jobs = {
   {"lambda=0", 1, I, -I, 1},
   {"lambda=pi/3", zeta12^2, I, -I, 1 - zeta12^2}
   };

matrixZeroQ[matrix_] := RootReduce[matrix] === ConstantArray[0, Dimensions[matrix]];
unitaryDefect[matrix_] := RootReduce[matrix . ConjugateTranspose[matrix] - 6 IdentityMatrix[6]];
modulusDefect[matrix_] := RootReduce[Map[RootReduce[# Conjugate[#] - 1] &, matrix, {2}]];

Print["Exact Dita-slice verification"];
Print["Field: Q(zeta_12) = Q(i, Sqrt[3])"];
Print["A = ", MatrixForm[a]];
Print["B = ", MatrixForm[b]];
Print["A identity: ", matrixZeroQ[RootReduce[a.ConjugateTranspose[a] - 2 IdentityMatrix[2]]]];
Print["B identity: ", matrixZeroQ[RootReduce[b.ConjugateTranspose[b] - 2 IdentityMatrix[2]]]];

Do[
  {tag, z1, z2, z3, z4} = job;
  z1sq = RootReduce[z1^2];
  z2sq = RootReduce[z2^2];
  z3sq = RootReduce[z3^2];
  z4sq = RootReduce[z4^2];
  z3FromA = mobius[z1sq, a[[1, 2]]^2, a[[1, 1]]^2];
  z4FromB = mobius[z1sq, b[[1, 2]]^2, b[[1, 1]]^2];
  z2FromB = RootReduce[(b[[1, 1]]^2 - z3sq conj12[b[[1, 2]]^2]) /
      (b[[1, 2]]^2 - z3sq conj12[b[[1, 1]]^2])];
  h = RootReduce[buildH[z1, z2, z3, z4]];
  ud = unitaryDefect[h];
  md = modulusDefect[h];
  Print["--- ", tag, " ---"];
  Print["z2^2 = ", z2sq, "; derived = ", z2FromB,
    "; match: ", RootReduce[z2sq - z2FromB] === 0];
  Print["z3^2 = ", z3sq, "; Mobius A(z1^2) match: ", RootReduce[z3sq - z3FromA] === 0];
  Print["z4^2 = ", z4sq, "; Mobius B(z1^2) match: ", RootReduce[z4sq - z4FromB] === 0];
  Print["H unitary defect zero: ", matrixZeroQ[ud]];
  Print["H modulus defect zero: ", matrixZeroQ[md]];
  Print["numeric max unitary defect: ", N[Max[Abs[ud]], 8]];
  Print["numeric max modulus defect: ", N[Max[Abs[md]], 8]];
  , {job, jobs}];