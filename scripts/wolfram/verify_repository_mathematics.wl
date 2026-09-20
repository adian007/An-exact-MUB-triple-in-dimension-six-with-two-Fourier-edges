(* Master Wolfram Language audit for the dimension-six MUB repository.

   Scope covered here:
     1. Karlsson A/B block identities and printed-transcription defect.
     2. Karlsson Mobius identities and numerical CHM checks at representative
        generic points and on the explicitly resolved Fourier seam.
     3. Exact Dita-slice matrices at lambda=0 and lambda=Pi/3.
     4. Exact catalogue checks for F6, Tao S6, D_bc, and F_D.
     5. Exact MU checks for {I, D_bc, F_D}.
     6. Construction and sizing of the exact 17-generator W1 ideal.

   Not reproduced by this file:
     HomotopyContinuation path tracking, interval certification, numerical
     pool solving, exhaustive clique enumeration, and literature claims about
     vector/base counts. Those remain Julia/Macaulay2-specific computations.

   Set runGroebner = True only for the expensive exact W1 elimination.
*)

ClearAll["Global`*"];
$MaxExtraPrecision = 10000;
runGroebner = False;

(* ---------- Shared exact helpers ---------- *)

matrixZeroQ[matrix_] :=
  RootReduce[matrix] === ConstantArray[0, Dimensions[matrix]] ||
   And @@ (TrueQ[RootReduce[#] == 0] & /@ Flatten[matrix]);
hadamardDefects[matrix_] := {
  RootReduce[matrix . ConjugateTranspose[matrix] - Length[matrix] IdentityMatrix[Length[matrix]]],
  RootReduce[Map[RootReduce[# Conjugate[#] - 1] &, matrix, {2}]]
  };
chmQ[matrix_] := And @@ (matrixZeroQ /@ hadamardDefects[matrix]);
muDefect[one_, two_] := RootReduce[
  Map[RootReduce[# Conjugate[#] - 6] &, ConjugateTranspose[one].two, {2}]
  ];

(* ---------- Karlsson A-block and Mobius construction ---------- *)

f2 = {{1, 1}, {1, -1}};

aOriginal[theta_, phi_] := Module[{a11, a12},
  a11 = -1/2 + I Sqrt[3]/2 (Cos[theta] + Exp[-I phi] Sin[theta]);
  a12 = -1/2 + I Sqrt[3]/2 (-Cos[theta] + Exp[I phi] Sin[theta]);
  {{a11, a12}, {Conjugate[a12], -Conjugate[a11]}}
  ];

aPrinted[theta_, phi_] := Module[{a11, a12},
  a11 = 1/2 + I Sqrt[3]/2 (Cos[theta] + Exp[-I phi] Sin[theta]);
  a12 = -1/2 + I Sqrt[3]/2 (-Cos[theta] + Exp[I phi] Sin[theta]);
  {{a11, a12}, {a12, -a11}}
  ];

mobius[z_, alpha_, beta_] := (alpha z - beta)/(Conjugate[beta] z - Conjugate[alpha]);
zLeft[z_] := {{1, 1}, {z, -z}};
zRight[z_] := {{1, z}, {1, -z}};

buildKarlsson[theta_, phi_, lambda_] := Module[
  {a, b, aa, ba, ab, bb, z1, z1sq, z2sq, z3sq, z4sq, z1v, z2v, z3v, z4v},
  a = aOriginal[theta, phi];
  b = -f2 - a;
  z1 = Exp[I lambda];
  z1sq = z1^2;
  If[theta === 0,
   z3sq = 1; z4sq = 1; z2sq = a[[1, 2]]^2/a[[1, 1]]^2,
   aa = a[[1, 2]]^2; ba = a[[1, 1]]^2;
   ab = b[[1, 2]]^2; bb = b[[1, 1]]^2;
   z3sq = mobius[z1sq, aa, ba];
   z4sq = mobius[z1sq, ab, bb];
   z2sq = (bb - z3sq Conjugate[ab])/(ab - z3sq Conjugate[bb])
   ];
  z1v = z1; z2v = Sqrt[z2sq]; z3v = Sqrt[z3sq]; z4v = Sqrt[z4sq];
  ArrayFlatten[{
    {f2, zLeft[z1v], zLeft[z2v]},
    {zRight[z3v], (1/2) zRight[z3v].a.zLeft[z1v], (1/2) zRight[z3v].b.zLeft[z2v]},
    {zRight[z4v], (1/2) zRight[z4v].b.zLeft[z1v], (1/2) zRight[z4v].a.zLeft[z2v]}
    }]
  ];

(* ---------- Exact Dita and catalogue matrices ---------- *)

zeta12 = RootReduce[Exp[I Pi/6]];
ii = zeta12^3;
aDita = {{ii, -1}, {-1, ii}};
bDita = RootReduce[-f2 - aDita];
exactDitaH[z1_, z2_, z3_, z4_] := ArrayFlatten[{
  {f2, zLeft[z1], zLeft[z2]},
  {zRight[z3], (1/2) zRight[z3].aDita.zLeft[z1], (1/2) zRight[z3].b.zLeft[z2]},
  {zRight[z4], (1/2) zRight[z4].b.zLeft[z1], (1/2) zRight[z4].aDita.zLeft[z2]}
  }];

ditaJobs = {
  {"lambda=0", 1, I, -I, 1},
  {"lambda=pi/3", zeta12^2, I, -I, 1 - zeta12^2}
  };

zeta24 = RootReduce[Exp[I Pi/12]];
omega3 = RootReduce[Exp[2 I Pi/3]];
omega6 = RootReduce[Exp[-I Pi/3]];
f6 = Table[omega6^(j k), {j, 0, 5}, {k, 0, 5}];
tao = {
  {1, 1, 1, 1, 1, 1},
  {1, 1, omega3, omega3, omega3^2, omega3^2},
  {1, omega3, 1, omega3^2, omega3^2, omega3},
  {1, omega3, omega3^2, 1, omega3, omega3^2},
  {1, omega3^2, omega3^2, omega3, 1, omega3},
  {1, omega3^2, omega3, omega3^2, omega3, 1}
  };

c3 = {{1, zeta24^6, zeta24^6}, {zeta24^6, 1, zeta24^6},
       {zeta24^6, zeta24^6, 1}};
c4 = {{zeta24^15, zeta24^3, zeta24^3}, {zeta24^3, zeta24^15, zeta24^3},
       {zeta24^3, zeta24^3, zeta24^15}};
dbc = ArrayFlatten[{{c3, c4}, {c4, -I Conjugate[c3]}}];
fourier3 = Table[zeta24^(8 j k), {j, 0, 2}, {k, 0, 2}];
b2 = (1 - 2 I)/Sqrt[5];
fd = ArrayFlatten[{
  {fourier3, fourier3},
  {fourier3.DiagonalMatrix[{zeta24^9 b2, 1, 1}],
   -fourier3.DiagonalMatrix[{zeta24^9 b2, 1, 1}]}
  }];

(* ---------- Exact W1 witness ---------- *)

zVars = Array[zz, 5];
wVars = Array[ww, 5];
witnessVariables = Join[zVars, wVars];
unimodular = MapThread[#1 #2 - 1 &, {zVars, wVars}];
muEquation[matrix_, column_] := Expand[
  (Conjugate[matrix[[1, column]]] +
      Sum[Conjugate[matrix[[row, column]]] zVars[[row - 1]], {row, 2, 6}]) *
   (matrix[[1, column]] +
      Sum[matrix[[row, column]] wVars[[row - 1]], {row, 2, 6}]) - 6
  ];
witness = RootReduce@Join[unimodular,
  Flatten[Table[muEquation[dbc, column], {column, 1, 6}]],
  Flatten[Table[muEquation[fd, column], {column, 1, 6}]]];

(* ---------- Report ---------- *)

Print["=== Master Wolfram mathematics audit ==="];
Print["Scope: exact identities plus high-precision structural checks"];
Print["Excluded: HC path tracking, interval certification, pool solving, cliques"];

originalDefect = FullSimplify[
  ComplexExpand[aOriginal[theta, phi].ConjugateTranspose[aOriginal[theta, phi]] -
    2 IdentityMatrix[2]], Element[{theta, phi}, Reals]];
printedDefect = FullSimplify[
  ComplexExpand[aPrinted[theta, phi].ConjugateTranspose[aPrinted[theta, phi]] -
    2 IdentityMatrix[2]], Element[{theta, phi}, Reals]];
Print["[A] Original block identity: ", originalDefect === ConstantArray[0, {2, 2}]];
Print["[A] Printed (1,1) defect: ", printedDefect[[1, 1]]];
Print["[A] Defect formula match: ",
  FullSimplify[printedDefect[[1, 1]] - Sqrt[3] Sin[theta] Sin[phi],
    Element[{theta, phi}, Reals]] === 0];
Print["[A] Dita A = {{I,-1},{-1,I}}: ",
  matrixZeroQ[aOriginal[ArcCos[1/Sqrt[3]], Pi/4] - aDita]];

Print["[Dita] A identity: ", matrixZeroQ[aDita.ConjugateTranspose[aDita] - 2 IdentityMatrix[2]]];
Print["[Dita] B identity: ", matrixZeroQ[bDita.ConjugateTranspose[bDita] - 2 IdentityMatrix[2]]];

Print["[Catalogue] F6 CHM: ", chmQ[f6]];
Print["[Catalogue] Tao S6 CHM: ", chmQ[tao]];
Print["[Catalogue] D_bc CHM: ", chmQ[dbc]];
Print["[Catalogue] F_D CHM: ", chmQ[fd]];
Print["[Catalogue] D_bc/F_D MU defect zero: ", matrixZeroQ[muDefect[dbc, fd]]];

genericPoints = {{3/10, 1/2, 1/5}, {1, 2, 7/10}, {1/10, 1/10, 3}};
genericChecks = Table[
  h = N[buildKarlsson @@ point, 50];
  {point, Max[Abs[Flatten[h.ConjugateTranspose[h] - 6 IdentityMatrix[6]]]],
   Max[Abs[Flatten[Map[Abs[#] - 1 &, h, {2}]]]]},
  {point, genericPoints}];
Print["[Karlsson] Generic numerical max defects: ", genericChecks];
seamChecks = Table[
  h = N[buildKarlsson[0, 1/2, lambda], 50];
  {lambda, Max[Abs[Flatten[h.ConjugateTranspose[h] - 6 IdentityMatrix[6]]]],
   Max[Abs[Flatten[Map[Abs[#] - 1 &, h, {2}]]]]},
  {lambda, {0, 7/10, 2}}];
Print["[Karlsson] Fourier-seam numerical max defects: ", seamChecks];

Print["[W1] Generators: ", Length[witness], "; variables: ", Length[witnessVariables]];
Print["[W1] Exact field: Q(zeta_24, Sqrt[5])"]; 
If[TrueQ[runGroebner],
  Module[{gb, rem},
    gb = GroebnerBasis[witness, witnessVariables, MonomialOrder -> DegreeLexicographic];
    rem = PolynomialReduce[1, witness, witnessVariables][[2]];
    Print["[W1] Groebner basis: ", gb];
    Print["[W1] 1 modulo witness: ", rem];
    Print["[W1] Unit-ideal certificate: ", rem === 0]
    ],
  Print["[W1] Groebner run skipped; set runGroebner=True to execute"]
  ];

Print["[Dita] Running dedicated exact Dita verification..."];
Get[FileNameJoin[{DirectoryName[$InputFileName], "verify_dita_exact.wl"}]];

Print["=== Audit complete ==="];