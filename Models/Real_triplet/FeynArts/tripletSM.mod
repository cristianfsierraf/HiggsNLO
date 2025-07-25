(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"ReTripletModel",
  Authors -> {"Cristian Sierra"},
  Emails -> {"cfsierraf@gmail.com"},
  Version -> "1.0",
  Date -> "2025-04-29",
  Institutions -> {"Tsung-Dao Lee Institute"},
  References -> {"arXiv:2104.10709"}};

FR$ClassesTranslation={};

FR$InteractionOrderPerturbativeExpansion={{HIG, 0}, {QCD, 0}, {QED, 0}};

FR$GoldstoneList={S[2], S[3]};

(*     Declared indices    *)

IndexRange[ Index[Generation] ] = Range[ 3 ]

IndexRange[ Index[Colour] ] = NoUnfold[ Range[ 3 ] ]

IndexRange[ Index[Gluon] ] = NoUnfold[ Range[ 8 ] ]

IndexRange[ Index[SU2W] ] = Range[ 3 ]

(*     Declared particles    *)

M$ClassesDescription = {
F[1] == {
    SelfConjugate -> False,
    Indices -> {Index[Generation]},
    PropagatorLabel -> "v",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Mnu },

F[2] == {
    SelfConjugate -> False,
    Indices -> {Index[Generation]},
    QuantumNumbers -> {-Q},
    PropagatorLabel -> "l",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Ml },

F[3] == {
    SelfConjugate -> False,
    Indices -> {Index[Generation], Index[Colour]},
    QuantumNumbers -> {(2*Q)/3},
    PropagatorLabel -> "uq",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Mu },

F[4] == {
    SelfConjugate -> False,
    Indices -> {Index[Generation], Index[Colour]},
    QuantumNumbers -> {-1/3*Q},
    PropagatorLabel -> "dq",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Md },

V[1] == {
    SelfConjugate -> True,
    Indices -> {},
    PropagatorLabel -> "\\gamma",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> 0 },

V[2] == {
    SelfConjugate -> True,
    Indices -> {},
    PropagatorLabel -> "Z",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> MZ },

V[3] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {Q},
    PropagatorLabel -> "W",
    PropagatorType -> Sine,
    PropagatorArrow -> Forward,
    Mass -> MW },

V[4] == {
    SelfConjugate -> True,
    Indices -> {Index[Gluon]},
    PropagatorLabel -> "G",
    PropagatorType -> Cycles,
    PropagatorArrow -> None,
    Mass -> 0 },

U[1] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> uA,
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> 0 },

U[2] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> uZ,
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> MZ },

U[3] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {Q, GhostNumber},
    PropagatorLabel -> uWp,
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> MW },

U[4] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {-Q, GhostNumber},
    PropagatorLabel -> uWm,
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> MW },

U[5] == {
    SelfConjugate -> False,
    Indices -> {Index[Gluon]},
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> uG,
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> 0 },

S[1] == {
    SelfConjugate -> True,
    Indices -> {},
    PropagatorLabel -> "h",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> mh },

S[2] == {
    SelfConjugate -> True,
    PropagatorLabel -> "G0",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> MZ,
    Indices -> {} },

S[3] == {
    SelfConjugate -> False,
    PropagatorLabel -> "Gch",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    QuantumNumbers -> {Q},
    Mass -> MW,
    Indices -> {} },

S[4] == {
    SelfConjugate -> True,
    Indices -> {},
    PropagatorLabel -> "H0",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> MH0 },

S[6] == {
    SelfConjugate -> False,
    Indices -> {},
    QuantumNumbers -> {Q},
    PropagatorLabel -> "Hch",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> Forward,
    Mass -> MHch }
}


(*        Definitions       *)


Mnu[ 1 ] := Mnue;
Mnu[ 2 ] := Mnum;
Mnu[ 3 ] := Mnut;
Ml[ 1 ] := ME;
Ml[ 2 ] := MM;
Ml[ 3 ] := MTA;
Mu[ 1, _ ] := MU;
Mu[ 1 ] := MU;
Mu[ 2, _ ] := MC;
Mu[ 2 ] := MC;
Mu[ 3, _ ] := MT;
Mu[ 3 ] := MT;
Md[ 1, _ ] := MD;
Md[ 1 ] := MD;
Md[ 2, _ ] := MS;
Md[ 2 ] := MS;
Md[ 3, _ ] := MB;
Md[ 3 ] := MB;
MZ[ ___ ] := MZ;
MW[ ___ ] := MW;
mh[ ___ ] := mh;
MH0[ ___ ] := MH0;
MHch[ ___ ] := MHch;


TheLabel[ F[1, {1}] ] := "ve";
TheLabel[ F[1, {2}] ] := "vm";
TheLabel[ F[1, {3}] ] := "vt";
TheLabel[ F[2, {1}] ] := "e";
TheLabel[ F[2, {2}] ] := "m";
TheLabel[ F[2, {3}] ] := "tau";
TheLabel[ F[3, {1, _}] ] := "u";
TheLabel[ F[3, {1}] ] := "u";
TheLabel[ F[3, {2, _}] ] := "c";
TheLabel[ F[3, {2}] ] := "c";
TheLabel[ F[3, {3, _}] ] := "t";
TheLabel[ F[3, {3}] ] := "t";
TheLabel[ F[4, {1, _}] ] := "d";
TheLabel[ F[4, {1}] ] := "d";
TheLabel[ F[4, {2, _}] ] := "s";
TheLabel[ F[4, {2}] ] := "s";
TheLabel[ F[4, {3, _}] ] := "b";
TheLabel[ F[4, {3}] ] := "b";
TheLabel[ V[4, {__}] ] := TheLabel[V[4]];
TheLabel[ U[5, {__}] ] := TheLabel[U[5]];


(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[1] , S[1] , S[1] , S[1] ] == {{(-6*I)*lam, 0}},

C[ S[4] , S[4] , S[4] , S[4] ] == {{(-6*I)*lam2, 0}},

C[ S[4] , S[4] , S[6] , -S[6] ] == {{(-2*I)*lam2, 0}},

C[ S[6] , S[6] , -S[6] , -S[6] ] == {{(-4*I)*lam2, 0}},

C[ S[1] , S[1] , S[4] , S[4] ] == {{(-I)*lam3, 0}},

C[ S[1] , S[1] , S[6] , -S[6] ] == {{(-I)*lam3, 0}},

C[ S[1] , S[1] , S[1] ] == {{(-6*I)*lam*v, 0}},

C[ S[1] , S[4] , S[4] ] == {{(-I)*lam3*v, 0}},

C[ S[1] , S[6] , -S[6] ] == {{(-I)*lam3*v, 0}},

C[ S[6] , -S[6] , V[1] , V[1] ] == {{(2*I)*EL^2, 0}},

C[ S[6] , -S[6] , V[1] ] == {{(-I)*gc11, 0}, {I*gc11, 0}},

C[ S[1] , V[4, {e1x2}] , V[4, {e2x2}] ] == {{(-I)*gc12*IndexDelta[e1x2, e2x2], 0}, {I*gc12*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -U[5, {e1x1}] , U[5, {e2x1}] , V[4, {e3x2}] ] == {{gc13*SUNF[e3x2, e1x1, e2x1], 0}, {gc13*SUNF[e3x2, e1x1, e2x1], 0}, {0, 0}},

C[ V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] ] == {{-(gc14*SUNF[e1x2, e2x2, e3x2]), 0}, {gc14*SUNF[e1x2, e2x2, e3x2], 0}, {gc14*SUNF[e1x2, e2x2, e3x2], 0}, {-(gc14*SUNF[e1x2, e2x2, e3x2]), 0}, {-(gc14*SUNF[e1x2, e2x2, e3x2]), 0}, {gc14*SUNF[e1x2, e2x2, e3x2], 0}},

C[ V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] , V[4, {e4x2}] ] == {{(-I)*gc15*(SUNF[e1x2, e2x2, e3x2, e4x2] + SUNF[e1x2, e3x2, e2x2, e4x2]), 0}, {I*gc15*(SUNF[e1x2, e2x2, e3x2, e4x2] - SUNF[e1x2, e4x2, e2x2, e3x2]), 0}, {I*gc15*(SUNF[e1x2, e3x2, e2x2, e4x2] + SUNF[e1x2, e4x2, e2x2, e3x2]), 0}},

C[ S[1] , V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] ] == {{-(gc16*SUNF[e1x2, e2x2, e3x2]), 0}, {gc16*SUNF[e1x2, e2x2, e3x2], 0}, {gc16*SUNF[e1x2, e2x2, e3x2], 0}, {-(gc16*SUNF[e1x2, e2x2, e3x2]), 0}, {-(gc16*SUNF[e1x2, e2x2, e3x2]), 0}, {gc16*SUNF[e1x2, e2x2, e3x2], 0}},

C[ S[1] , V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] , V[4, {e4x2}] ] == {{(-I)*gc17*(SUNF[e1x2, e2x2, e3x2, e4x2] + SUNF[e1x2, e3x2, e2x2, e4x2]), 0}, {I*gc17*(SUNF[e1x2, e2x2, e3x2, e4x2] - SUNF[e1x2, e4x2, e2x2, e3x2]), 0}, {I*gc17*(SUNF[e1x2, e3x2, e2x2, e4x2] + SUNF[e1x2, e4x2, e2x2, e3x2]), 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , S[1] ] == {{I*gc18[e1x2]*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc18[e1x2]*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , S[1] ] == {{I*gc19[e1x2]*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc19[e1x2]*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , S[1] ] == {{0, 0}, {I*gc20R[e1x2]*IndexDelta[e1x2, e2x2], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[4, {e3x2}] ] == {{I*gc21*IndexDelta[e1x2, e2x2]*SUNT[e3x2, e1x3, e2x3], 0}, {I*gc21*IndexDelta[e1x2, e2x2]*SUNT[e3x2, e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[4, {e3x2}] ] == {{I*gc22*IndexDelta[e1x2, e2x2]*SUNT[e3x2, e1x3, e2x3], 0}, {I*gc22*IndexDelta[e1x2, e2x2]*SUNT[e3x2, e1x3, e2x3], 0}},

C[ S[4] , -S[6] , V[1] , V[3] ] == {{((-I)*EL^2)/SW, 0}},

C[ S[4] , -S[6] , V[3] ] == {{(-I)*gc24, 0}, {I*gc24, 0}},

C[ V[1] , V[3] , -V[3] ] == {{(-I)*gc25, 0}, {I*gc25, 0}, {I*gc25, 0}, {(-I)*gc25, 0}, {(-I)*gc25, 0}, {I*gc25, 0}},

C[ -S[6] , -S[6] , V[3] , V[3] ] == {{((-2*I)*EL^2)/SW^2, 0}},

C[ S[4] , S[6] , V[1] , -V[3] ] == {{((-I)*EL^2)/SW, 0}},

C[ S[4] , S[6] , -V[3] ] == {{(-I)*gc28, 0}, {I*gc28, 0}},

C[ S[1] , S[1] , V[3] , -V[3] ] == {{((I/2)*EL^2)/SW^2, 0}},

C[ S[4] , S[4] , V[3] , -V[3] ] == {{((2*I)*EL^2)/SW^2, 0}},

C[ S[6] , -S[6] , V[3] , -V[3] ] == {{(I*EL^2)/SW^2, 0}},

C[ S[1] , V[3] , -V[3] ] == {{0, 0}, {0, 0}, {I*gc32, 0}},

C[ V[1] , V[1] , V[3] , -V[3] ] == {{(-I)*gc33, 0}, {(-I)*gc33, 0}, {(2*I)*gc33, 0}},

C[ V[3] , -V[3] , V[2] ] == {{(-I)*gc34, 0}, {I*gc34, 0}, {I*gc34, 0}, {(-I)*gc34, 0}, {(-I)*gc34, 0}, {I*gc34, 0}},

C[ S[6] , S[6] , -V[3] , -V[3] ] == {{((-2*I)*EL^2)/SW^2, 0}},

C[ V[3] , V[3] , -V[3] , -V[3] ] == {{(-I)*gc36, 0}, {(-I)*gc36, 0}, {(2*I)*gc36, 0}},

C[ S[6] , -S[6] , V[1] , V[2] ] == {{((2*I)*CW*EL^2)/SW, 0}},

C[ S[6] , -S[6] , V[2] ] == {{(-I)*gc38, 0}, {I*gc38, 0}},

C[ S[4] , -S[6] , V[3] , V[2] ] == {{((-I)*CW*EL^2)/SW^2, 0}},

C[ S[4] , S[6] , -V[3] , V[2] ] == {{((-I)*CW*EL^2)/SW^2, 0}},

C[ V[1] , V[3] , -V[3] , V[2] ] == {{(-2*I)*gc41, 0}, {I*gc41, 0}, {I*gc41, 0}},

C[ S[1] , S[1] , V[2] , V[2] ] == {{((I/2)*EL^2*(CW^2 + SW^2)^2)/(CW^2*SW^2), 0}},

C[ S[6] , -S[6] , V[2] , V[2] ] == {{((2*I)*CW^2*EL^2)/SW^2, 0}},

C[ S[1] , V[2] , V[2] ] == {{0, 0}, {0, 0}, {I*gc44, 0}},

C[ V[3] , -V[3] , V[2] , V[2] ] == {{(-I)*gc45, 0}, {(-I)*gc45, 0}, {(2*I)*gc45, 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[1] ] == {{I*gc46*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc46*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , V[1] ] == {{I*gc47*IndexDelta[e1x2, e2x2], 0}, {I*gc47*IndexDelta[e1x2, e2x2], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[1] ] == {{I*gc48*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc48*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[1, {e1x2}] , F[2, {e2x2}] , V[3] ] == {{I*gc49*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -F[3, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[3] ] == {{I*gc50[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {0, 0}},

C[ -F[2, {e1x2}] , F[1, {e2x2}] , -V[3] ] == {{I*gc51*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -F[4, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , -V[3] ] == {{I*gc52[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {0, 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[2] ] == {{I*gc53L*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc53R*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , V[2] ] == {{I*gc54L*IndexDelta[e1x2, e2x2], 0}, {I*gc54R*IndexDelta[e1x2, e2x2], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[2] ] == {{I*gc55L*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc55R*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[1, {e1x2}] , F[1, {e2x2}] , V[2] ] == {{I*gc56*IndexDelta[e1x2, e2x2], 0}, {0, 0}}

}

(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)

(* Parameter replacement lists (These lists were created by FeynRules) *)

(* FA Couplings *)

M$FACouplings = {
     gc11 -> EL,
     gc12 -> GH,
     gc13 -> GS,
     gc14 -> -GS,
     gc15 -> -GS^2,
     gc16 -> -(GH*GS),
     gc17 -> -(GH*GS^2),
     gc18[e1x2_] -> -(Md[e1x2]/v),
     gc19[e1x2_] -> -(Mu[e1x2]/v),
     gc20R[e1x2_] -> -(Ml[e1x2]/v),
     gc21 -> GS,
     gc22 -> GS,
     gc24 -> -(EL/SW),
     gc25 -> EL,
     gc28 -> EL/SW,
     gc32 -> (EL^2*v)/(2*SW^2),
     gc33 -> -EL^2,
     gc34 -> (CW*EL)/SW,
     gc36 -> EL^2/SW^2,
     gc38 -> (CW*EL)/SW,
     gc41 -> (CW*EL^2)/SW,
     gc44 -> (EL^2*(CW^2 + SW^2)^2*v)/(2*CW^2*SW^2),
     gc45 -> -((CW^2*EL^2)/SW^2),
     gc46 -> -1/3*EL,
     gc47 -> -EL,
     gc48 -> (2*EL)/3,
     gc49 -> EL/(Sqrt[2]*SW),
     gc50[e1x2_, e2x2_] -> (EL*CKM[e1x2, e2x2])/(Sqrt[2]*SW),
     gc51 -> EL/(Sqrt[2]*SW),
     gc52[e1x2_, e2x2_] -> (EL*Conjugate[CKM[e2x2, e1x2]])/(Sqrt[2]*SW),
     gc53L -> -1/6*(EL*(3*CW^2 + SW^2))/(CW*SW),
     gc53R -> (EL*SW)/(3*CW),
     gc54L -> -1/2*(EL*(CW^2 - SW^2))/(CW*SW),
     gc54R -> (EL*SW)/CW,
     gc55L -> (CW*EL)/(2*SW) - (EL*SW)/(6*CW),
     gc55R -> (-2*EL*SW)/(3*CW),
     gc56 -> (EL*(CW^2 + SW^2))/(2*CW*SW)};

