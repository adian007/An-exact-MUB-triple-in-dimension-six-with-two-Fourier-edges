(* Exact n_wit=1 witness at the D0-equivalent pair (D_bc, F_D).
   Companion to symbolic_export/w1_D0_exact.m2 and w1_D0_groebner.m2.

   The system has 17 generators in 10 variables over
   Q(zeta_24, Sqrt[5]). Set runGroebner = True to run the exact
   Groebner calculation in Wolfram Language.
*)

ClearAll["Global`*"];
$MaxExtraPrecision = 10000;

(* Keep the fast construction check as the default. Set True for the full run. *)
runGroebner = False;

zeta24 = RootReduce[Exp[I Pi/12]];
s5 = Sqrt[5];

c3 = {
   {1, zeta24^6, zeta24^6},
   {zeta24^6, 1, zeta24^6},
   {zeta24^6, zeta24^6, 1}
   };
c4 = {
   {zeta24^15, zeta24^3, zeta24^3},
   {zeta24^3, zeta24^15, zeta24^3},
   {zeta24^3, zeta24^3, zeta24^15}
   };
dbc = ArrayFlatten[{{c3, c4}, {c4, -I Conjugate[c3]}}];

fourier3 = Table[zeta24^(8 j k), {j, 0, 2}, {k, 0, 2}];
b2 = (1 - 2 I)/s5;
dDiag = DiagonalMatrix[{zeta24^9 b2, 1, 1}];
fourier3d = fourier3 . dDiag;
fd = ArrayFlatten[{{fourier3, fourier3}, {fourier3d, -fourier3d}}];

zVars = Array[zz, 5];
wVars = Array[ww, 5];
variables = Join[zVars, wVars];

unimodular = MapThread[#1 #2 - 1 &, {zVars, wVars}];
muEquation[matrix_, column_] := Expand[
   (Conjugate[matrix[[1, column]]] +
       Sum[Conjugate[matrix[[row, column]]] zVars[[row - 1]], {row, 2, 6}]) *
    (matrix[[1, column]] +
       Sum[matrix[[row, column]] wVars[[row - 1]], {row, 2, 6}]) - 6
   ];

witness = RootReduce@Join[
    unimodular,
    Flatten[Table[muEquation[dbc, column], {column, 1, 6}]],
    Flatten[Table[muEquation[fd, column], {column, 1, 6}]]
    ];

Print["Exact D0-equivalent W1 witness"];
Print["Field: Q(zeta_24, Sqrt[5])"];
Print["Generators: ", Length[witness], "; variables: ", Length[variables]];
Print["CHM defects: ", RootReduce[dbc . ConjugateTranspose[dbc] - 6 IdentityMatrix[6]] === ConstantArray[0, {6, 6}], ", ",
  RootReduce[fd . ConjugateTranspose[fd] - 6 IdentityMatrix[6]] === ConstantArray[0, {6, 6}]];

If[TrueQ[runGroebner],
  Module[{basis, remainder},
   {basis, remainder} = TimeConstrained[
      {GroebnerBasis[witness, variables, MonomialOrder -> DegreeLexicographic],
       PolynomialReduce[1, witness, variables][[2]]},
      Infinity,
      {"TIMEOUT", "TIMEOUT"}
      ];
   Print["Groebner basis: ", basis];
   Print["1 modulo witness: ", remainder];
   Print["CERTIFICATE: ", remainder === 0]
   ]
  ];