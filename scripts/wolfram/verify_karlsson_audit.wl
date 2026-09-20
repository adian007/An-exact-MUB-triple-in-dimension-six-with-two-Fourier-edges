(* Independent symbolic verification of the Karlsson A-block audit.
   This checks the load-bearing formula used by the Julia implementation:
     A A^dagger = 2 I
   and the exact defect of the printed McNulty-Weigert transcription.
*)

ClearAll["Global`*"];

assumptions = Element[{theta, phi}, Reals];
f2 = {{1, 1}, {1, -1}};

a11 = -1/2 + I Sqrt[3]/2 (Cos[theta] + Exp[-I phi] Sin[theta]);
a12 = -1/2 + I Sqrt[3]/2 (-Cos[theta] + Exp[I phi] Sin[theta]);

aOriginal = {{a11, a12}, {Conjugate[a12], -Conjugate[a11]}};
aPrinted = {{a11 + 1, a12}, {a12, -a11 - 1}};

originalDefect = FullSimplify[
   ComplexExpand[aOriginal . ConjugateTranspose[aOriginal] - 2 IdentityMatrix[2]],
   assumptions
   ];
printedDefect = FullSimplify[
   ComplexExpand[aPrinted . ConjugateTranspose[aPrinted] - 2 IdentityMatrix[2]],
   assumptions
   ];

ditaA = FullSimplify[
   aOriginal /. {theta -> ArcCos[1/Sqrt[3]], phi -> Pi/4},
   assumptions
   ];
ditaB = FullSimplify[-f2 - ditaA];

Print["Karlsson A-block symbolic audit"];
Print["Original AA^dagger - 2 I = ", MatrixForm[originalDefect]];
Print["Printed (1,1) defect = ", printedDefect[[1, 1]]];
Print["Target defect = ", Sqrt[3] Sin[theta] Sin[phi]];
Print["Printed defect identity: ",
  FullSimplify[printedDefect[[1, 1]] - Sqrt[3] Sin[theta] Sin[phi], assumptions] === 0];
Print["Original block identity: ", originalDefect === ConstantArray[0, {2, 2}]];
Print["Dita A = ", MatrixForm[ditaA]];
Print["Dita B = ", MatrixForm[ditaB]];
Print["Dita A target: ", ditaA === {{I, -1}, {-1, I}}];
Print["Dita A identity: ", FullSimplify[ditaA . ConjugateTranspose[ditaA] - 2 IdentityMatrix[2]] === ConstantArray[0, {2, 2}]];
Print["Dita B identity: ", FullSimplify[ditaB . ConjugateTranspose[ditaB] - 2 IdentityMatrix[2]] === ConstantArray[0, {2, 2}]];