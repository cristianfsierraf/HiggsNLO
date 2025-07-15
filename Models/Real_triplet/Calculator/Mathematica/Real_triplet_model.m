(* ::Package:: *)

(* ::Input:: *)
(*(*Code to calculate semi-automatically the NLO Higgs-strahlung cross section in the real triplet model*)*)
(*(*Authors:*)
(*      Main authors: Cristian Sierra, cristian.sierra@cern.ch*)
(*                    Qin Yu,*)
(*       Contributors: Mohamed Sassi,*)
(*                     Anjan Barik,*)
(*                     Manuel Diaz*)*)
(*(*Date: 04-07-2025*)*)


(* ::Section:: *)
(*Load FeynCalc*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Global`$LoadAddOns={"FeynHelpers"};*)
(*$LoadFeynArts=True;*)
(*<<FeynCalc`*)
(*$FAVerbose=0;*)
(*$CKM=True;*)


(* ::Input:: *)
(*FAPatch[PatchModelsOnly->True]*)


(* ::Title:: *)
(*Higgsstrahlung@NLO*)


(* ::Chapter::Closed:: *)
(*Self-energies*)


(* ::Section:: *)
(*H*)


(* ::Input:: *)
(*Topo\[CapitalSigma]=CreateTopologies[1,1->1,ExcludeTopologies->{AllBoxes,Triangles}];*)


(* ::Input:: *)
(*\[CapitalSigma]HHDiag=InsertFields[Topo\[CapitalSigma],{S[1]}->{S[1]},Model->{"tripletSM"},*)
(*GenericModel->{"tripletSM"},InsertionLevel->{Particles},*)
(*ExcludeParticles->{S[1],V[1|2|3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[\[CapitalSigma]HHDiag,ColumnsXRows->{2,2},Numbering->None];*)


(* ::Input:: *)
(*\[CapitalSigma]HHAmp[0]=FCFAConvert[CreateFeynAmp[\[CapitalSigma]HHDiag,PreFactor->-I,Truncated->True],*)
(*IncomingMomenta->{k},OutgoingMomenta->{k},TransversePolarizationVectors->{k},*)
(*Contract->True,LoopMomenta->{q},ChangeDimension->D,*)
(*UndoChiralSplittings->True,List->False,DropSumOver->False,*)
(*LorentzIndexNames->{\[Micro],\[Nu]},FinalSubstitutions->{lam3->\[Lambda]3,MH0->MH,MHch->MH}]//*)
(*DiracSimplify[#,DiracSubstitute67->True]&//FCGVToSymbol//DiracSimplify;*)


(* ::Input:: *)
(*(*When setting to True the following command rewrites A0 in terms of B0*)*)


(* ::Input:: *)
(*SetOptions[A0,A0ToB0->False];*)


(* ::Input:: *)
(*(*Tensor decomposition into Passarino Veltman functions:*)*)


(* ::Input:: *)
(*\[CapitalSigma]HHAmp[1]=  TID[\[CapitalSigma]HHAmp[0],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*\[CapitalSigma]HH[0]=1/(2 Pi)^4PaVeLimitTo4[\[CapitalSigma]HHAmp[1]];*)


(* ::Input:: *)
(*\[CapitalSigma]HH[0]//FullSimplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalization in MSbar scheme*)


(* ::Input:: *)
(*\[CapitalSigma]HH[1]=\[CapitalSigma]HH[0]/.Pair[ Momentum[k], Momentum[k]]->kSquare;*)


(* ::Input:: *)
(*\[CapitalSigma]HHat[0]=FCHideEpsilon[PaXEvaluate[\[CapitalSigma]HH[1],PaXAnalytic->True]];*)


(* ::Input:: *)
(*\[CapitalSigma]HHat[1]=\[CapitalSigma]HHat[0]-SMP["Delta"]*Coefficient[\[CapitalSigma]HHat[0],SMP["Delta"]];*)


(* ::Input:: *)
(*(* Derivative of Overscript[\[CapitalSigma], ^]^h(k^2) *)*)


(* ::Input:: *)
(*D\[CapitalSigma]HHat[0]=D[\[CapitalSigma]HHat[1],kSquare];*)


(* ::Input:: *)
(*(* \[PartialD]Overscript[\[CapitalSigma], ^]^h/\[PartialD]k^2 at Subsuperscript[m, h, 2] *)*)


(* ::Input:: *)
(*D\[CapitalSigma]HHatFormh2=D\[CapitalSigma]HHat[0]/.kSquare->mh^2;*)


(* ::Input:: *)
(*\[CapitalSigma]HHatFors=\[CapitalSigma]HHat[1]/.kSquare->s;*)


(* ::Section:: *)
(*ZZ*)


(* ::Input:: *)
(*\[CapitalSigma]ZZDiag=InsertFields[Topo\[CapitalSigma],{V[2]}->{V[2]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[Except[6]],V[1|3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[\[CapitalSigma]ZZDiag,ColumnsXRows->{2,1},Numbering->Simple];*)


(* ::Input:: *)
(*\[CapitalSigma]ZZAmp[0]=DiracSimplify[FCFAConvert[CreateFeynAmp[\[CapitalSigma]ZZDiag,PreFactor->-I,Truncated->True],IncomingMomenta->{k},OutgoingMomenta->{k},TransversePolarizationVectors->{k},Contract->True,LoopMomenta->{q},ChangeDimension->D,UndoChiralSplittings->True,List->False,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu]},FinalSubstitutions->{MH0->MH,MHch->MH}],DiracSubstitute67->True]/.gc38->(FCGV["CW"]*FCGV["EL"])/FCGV["SW"]//DiracSimplify//FCGVToSymbol;*)


(* ::Input:: *)
(*(*Tensor decomposition into Passarino Veltman functions:*)*)


(* ::Input:: *)
(*\[CapitalSigma]ZZAmp[1]= TID[\[CapitalSigma]ZZAmp[0],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*SetOptions[A0,A0ToB0->True];*)


(* ::Input:: *)
(*\[CapitalSigma]ZZ=1/(2 Pi)^4PaVeLimitTo4[\[CapitalSigma]ZZAmp[1]];*)


(* ::Input:: *)
(*(*Projecting out the transverse part, definition from hep-ph/0612057 Eq.50 \[CapitalSigma]ZZ= (-g\[Micro]\[Nu]+k\[Micro]k\[Nu]/k^2)\[CapitalSigma]ZZT + (k\[Micro]k\[Nu]/k^2) \[CapitalSigma]ZZL*)*)


(* ::Input:: *)
(*TransverseProjector=(1/3)(-Pair[LorentzIndex[\[Nu]],LorentzIndex[\[Micro]]]+(Pair[LorentzIndex[\[Nu]],Momentum[k]] Pair[LorentzIndex[\[Micro]],Momentum[k]]/Pair[Momentum[k],Momentum[k]]));*)
(*\[CapitalSigma]ZZT[0]=Contract[TransverseProjector*\[CapitalSigma]ZZ];*)
(*\[CapitalSigma]ZZT[0]//Simplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalization in MSbar scheme*)


(* ::Input:: *)
(*\[CapitalSigma]ZZT[1]=\[CapitalSigma]ZZT[0]/.Pair[ Momentum[k], Momentum[k]]->kSquare;*)


(* ::Input:: *)
(*\[CapitalSigma]ZZTHat[0]=FCHideEpsilon[PaXEvaluate[\[CapitalSigma]ZZT[1],PaXAnalytic->True]];*)


(* ::Input:: *)
(*(* MS bar scheme *)*)


(* ::Input:: *)
(*\[CapitalSigma]ZZTHat[1]=\[CapitalSigma]ZZTHat[0]-SMP["Delta"]*Coefficient[\[CapitalSigma]ZZTHat[0],SMP["Delta"]];*)


(* ::Input:: *)
(*(* Subsuperscript[Overscript[\[CapitalSigma], ^], T, ZZ] at Subsuperscript[M, Z, 2] and s *)*)


(* ::Input:: *)
(*\[CapitalSigma]ZZTHatForMZ2=\[CapitalSigma]ZZTHat[1]/.kSquare->MZ^2;*)


(* ::Input:: *)
(*\[CapitalSigma]ZZTHatFors=\[CapitalSigma]ZZTHat[1]/.kSquare->s;*)


(* ::Input:: *)
(*(* Derivative of \!\(\**)
(*SubsuperscriptBox[*)
(*OverscriptBox["\[CapitalSigma]", "^"], *)
(*StyleBox["T", "TI"], *)
(*StyleBox[*)
(*RowBox[{"Z", "Z"}], "TI"]](\**)
(*SuperscriptBox[*)
(*StyleBox["k", "TI"], "2"])\) *)*)


(* ::Input:: *)
(*D\[CapitalSigma]ZZTHat[0]=D[\[CapitalSigma]ZZTHat[1],kSquare];*)


(* ::Input:: *)
(*(* \[PartialD]Subsuperscript[Overscript[\[CapitalSigma], ^], T, ZZ]/\[PartialD]k^2 at Subsuperscript[M, Z, 2] *)*)


(* ::Input:: *)
(*D\[CapitalSigma]ZZTHatforMZ2=D\[CapitalSigma]ZZTHat[0]/.kSquare->MZ^2;*)


(* ::Section:: *)
(*\[Gamma]Z*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZDiag=InsertFields[Topo\[CapitalSigma],{V[1]}->{V[2]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[Except[6]],V[1|3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[\[CapitalSigma]\[Gamma]ZDiag,ColumnsXRows->{2,1},Numbering->None];*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZAmp[0]=DiracSimplify[FCFAConvert[CreateFeynAmp[\[CapitalSigma]\[Gamma]ZDiag,PreFactor->-I,Truncated->True],IncomingMomenta->{k},OutgoingMomenta->{k},TransversePolarizationVectors->{k},Contract->True,LoopMomenta->{q},ChangeDimension->D,UndoChiralSplittings->True,List->False,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu]},FinalSubstitutions->{MH0->MH,MHch->MH}],DiracSubstitute67->True]/.gc11->FCGV["EL"]/.gc38->(FCGV["CW"]*FCGV["EL"])/FCGV["SW"]//DiracSimplify//FCGVToSymbol;*)


(* ::Input:: *)
(*(*Tensor decomposition into Passarino Veltman functions:*)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZAmp[1]= TID[\[CapitalSigma]\[Gamma]ZAmp[0],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]Z=1/(2 Pi)^4PaVeLimitTo4[\[CapitalSigma]\[Gamma]ZAmp[1]];*)


(* ::Input:: *)
(*(*Projecting out the transverse part*)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZT[0]=Contract[TransverseProjector*\[CapitalSigma]\[Gamma]Z];*)
(*\[CapitalSigma]\[Gamma]ZT[0]//Simplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalization in MS bar scheme*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZT[1]=\[CapitalSigma]\[Gamma]ZT[0]/.Pair[ Momentum[k], Momentum[k]]->kSquare;*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZTHat[0]=FCHideEpsilon[PaXEvaluate[\[CapitalSigma]\[Gamma]ZT[1],PaXAnalytic->True]];*)


(* ::Input:: *)
(*(* MS bar scheme *)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZTHat[1]=\[CapitalSigma]\[Gamma]ZTHat[0]-SMP["Delta"]*Coefficient[\[CapitalSigma]\[Gamma]ZTHat[0],SMP["Delta"]];*)


(* ::Input:: *)
(*(* Subsuperscript[Overscript[\[CapitalSigma], ^], T, \[Gamma]Z] at 0 and s *)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZTHatFor0=\[CapitalSigma]\[Gamma]ZTHat[1]//Series[#,{kSquare,0,0}]&//Normal;*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]ZTHatFors=\[CapitalSigma]\[Gamma]ZTHat[1]/.kSquare->s;*)


(* ::Section:: *)
(*\[Gamma]\[Gamma]*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]Diag=InsertFields[Topo\[CapitalSigma],{V[1]}->{V[1]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[Except[6]],V[1|3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[\[CapitalSigma]\[Gamma]\[Gamma]Diag,ColumnsXRows->{2,1},Numbering->None];*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]Amp[0]=DiracSimplify[FCFAConvert[CreateFeynAmp[\[CapitalSigma]\[Gamma]\[Gamma]Diag,PreFactor->-I,Truncated->True],IncomingMomenta->{k},OutgoingMomenta->{k},TransversePolarizationVectors->{k},Contract->True,LoopMomenta->{q},ChangeDimension->D,UndoChiralSplittings->True,List->False,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu]},FinalSubstitutions->{MH0->MH,MHch->MH}],DiracSubstitute67->True]/.gc11->FCGV["EL"]//DiracSimplify//FCGVToSymbol;*)


(* ::Input:: *)
(*(*Tensor decomposition into Passarino Veltman functions:*)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]Amp[1]= TID[\[CapitalSigma]\[Gamma]\[Gamma]Amp[0],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]=1/(2 Pi)^4PaVeLimitTo4[\[CapitalSigma]\[Gamma]\[Gamma]Amp[1]];*)


(* ::Input:: *)
(*(*Projecting out the transverse part*)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]T[0]=Contract[TransverseProjector*\[CapitalSigma]\[Gamma]\[Gamma]];*)
(*\[CapitalSigma]\[Gamma]\[Gamma]T[0]//FullSimplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalization in MS bar scheme*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]T[1]=\[CapitalSigma]\[Gamma]\[Gamma]T[0]/.Pair[ Momentum[k], Momentum[k]]->kSquare;*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]THat[0]=FCHideEpsilon[PaXEvaluate[\[CapitalSigma]\[Gamma]\[Gamma]T[1],PaXAnalytic->True]];*)


(* ::Input:: *)
(*(* MS bar scheme *)*)


(* ::Input:: *)
(*\[CapitalSigma]\[Gamma]\[Gamma]THat[1]=\[CapitalSigma]\[Gamma]\[Gamma]THat[0]-SMP["Delta"]*Coefficient[\[CapitalSigma]\[Gamma]\[Gamma]THat[0],SMP["Delta"]];*)


(* ::Input:: *)
(*(* Derivative of \!\(\**)
(*SubsuperscriptBox[*)
(*OverscriptBox["\[CapitalSigma]", "^"], *)
(*StyleBox["T", "TI"], "\[Gamma]\[Gamma]"](\**)
(*SuperscriptBox[*)
(*StyleBox["k", "TI"], "2"])\) *)*)


(* ::Input:: *)
(*D\[CapitalSigma]\[Gamma]\[Gamma]THat[0]=D[\[CapitalSigma]\[Gamma]\[Gamma]THat[1],kSquare];*)


(* ::Input:: *)
(*(* \[PartialD]Subsuperscript[Overscript[\[CapitalSigma], ^], T, \[Gamma]\[Gamma]]/\[PartialD]k^2 at 0 *)*)


(* ::Input:: *)
(*D\[CapitalSigma]\[Gamma]\[Gamma]THatFor0=D\[CapitalSigma]\[Gamma]\[Gamma]THat[0]//Series[#,{kSquare,0,0}]&//Normal;*)


(* ::Section:: *)
(*WW*)


(* ::Input:: *)
(*\[CapitalSigma]WWDiag=InsertFields[Topo\[CapitalSigma],{V[3]}->{V[3]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[1],V[1|2|3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[\[CapitalSigma]WWDiag,ColumnsXRows->{3,1},Numbering->None,ImageSize->{512,256}];*)


(* ::Input:: *)
(*\[CapitalSigma]WWAmp[0]=DiracSimplify[FCFAConvert[CreateFeynAmp[\[CapitalSigma]WWDiag,PreFactor->-I,Truncated->True],IncomingMomenta->{k},OutgoingMomenta->{k},TransversePolarizationVectors->{k},Contract->True,LoopMomenta->{q},ChangeDimension->D,UndoChiralSplittings->True,List->False,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu]},FinalSubstitutions->{MH0->MH,MHch->MH}],DiracSubstitute67->True]/.gc28->(FCGV["EL"])/FCGV["SW"]/.gc24->-((FCGV["EL"])/FCGV["SW"])//DiracSimplify//FCGVToSymbol;*)


(* ::Input:: *)
(*(*Tensor decomposition into Passarino Veltman functions:*)*)


(* ::Input:: *)
(*\[CapitalSigma]WWAmp[1]= TID[\[CapitalSigma]WWAmp[0],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*\[CapitalSigma]WW=1/(2 Pi)^4PaVeLimitTo4[\[CapitalSigma]WWAmp[1]];*)


(* ::Input:: *)
(*\[CapitalSigma]WWT[0]=Contract[TransverseProjector*\[CapitalSigma]WW];*)
(*\[CapitalSigma]WWT[0]//FullSimplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalization in MS bar scheme*)


(* ::Input:: *)
(*\[CapitalSigma]WWT[1]=\[CapitalSigma]WWT[0]/.Pair[ Momentum[k], Momentum[k]]->kSquare;*)


(* ::Input:: *)
(*\[CapitalSigma]WWTHat[0]=FCHideEpsilon[PaXEvaluate[\[CapitalSigma]WWT[1],PaXAnalytic->True]];*)


(* ::Input:: *)
(*(* MS bar scheme *)*)


(* ::Input:: *)
(*\[CapitalSigma]WWTHat[1]=\[CapitalSigma]WWTHat[0]-SMP["Delta"]*Coefficient[\[CapitalSigma]WWTHat[0],SMP["Delta"]];*)


(* ::Input:: *)
(*(* Subsuperscript[Overscript[\[CapitalSigma], ^], T, WW] at Subsuperscript[M, W, 2] *)*)


(* ::Input:: *)
(*\[CapitalSigma]WWTHatForMW2=\[CapitalSigma]WWTHat[1]/.kSquare->MW^2;*)


(* ::Chapter:: *)
(*Cross section*)


(* ::Section:: *)
(*LO*)


(* ::Subsubsection:: *)
(*Amplitude*)


(* ::Input:: *)
(*topoLO=CreateTopologies[0,2->2];*)
(*diagLO=InsertFields[topoLO,{-F[2,{1}],F[2,{1}]}->{V[2],S[1]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[1],S[2],V[3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}]; *)


(* ::Input:: *)
(*Paint[diagLO,ColumnsXRows->{1,1},Numbering->None,ImageSize->{256,256} ];*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)


(* ::Input:: *)
(*(*set the kinematics*)*)


(* ::Input:: *)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,ME,ME,MZ,mh];*)


(* ::Input:: *)
(*(*define the matrix element*)*)
(*(*PreFactor is set as -I so that ampLO[0]=M^tree instead of iM^tree*)*)


(* ::Input:: *)
(*ampLO[0]=FCFAConvert[CreateFeynAmp[diagLO,PreFactor->-I],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},TransversePolarizationVectors->{k1},ChangeDimension->D,Contract->True,UndoChiralSplittings->True,List->False,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu]}]//FCGVToSymbol//DiracSimplify[#,DiracSubstitute67->False]&;*)


(* ::Input:: *)
(*(*Check by run following code commented out: The LO amplitude is the same as M^tree in paper*)*)


(* ::Input:: *)
(*(*ampLO[0]/.gc44->(EL^2*(v)/(2*CW^2*SW^2))/.gc54L->(gv+ga) EL/.gc54R->(gv-ga) EL/.v->(2 SW CW MZ/EL)/.ME->0//DiracSimplify[#,DiracSubstitute67->True]&//FullSimplify*)*)


(* ::Input:: *)
(*(*Matrix element squared. Spin average 1/4 is considered here. D dim -> 4 dim*)*)


(* ::Input:: *)
(*ampLOSquared[0]=ampLO[0]( ComplexConjugate[ampLO[0]])//ChangeDimension[#,4]&//FeynAmpDenominatorExplicit//FermionSpinSum[#,ExtraFactor->1/2^2]&//DoPolarizationSums[#,k1]&//TrickMandelstam[#,{s,t,u,2ME^2+MZ^2+mh^2}]& //DiracSimplify;*)


(* ::Subsubsection:: *)
(*Matrix Element and cross section*)


(* ::Input:: *)
(*(*LO matrix element squared. We set ME=0 here*)*)


(* ::Input:: *)
(*MLO2aux=ampLOSquared[0]/.gc44->(EL^2*(v)/(2*CW^2*SW^2))/.gc54L->(gv+ga) EL/.gc54R->(gv-ga) EL/.v->(2 SW CW MZ/EL)/.ME->0/.u->(MZ^2+mh^2-s-t);*)


(* ::Input:: *)
(*(*The differential cross section is given by (dSubscript[\[Sigma], LO])/(dt)\[LongEqual]1/(16\[Pi]s^2)\[LeftBracketingBar]M\[RightBracketingBar]^2 where M is the matrix element. We integrate over t and replace u from s+t+u\[LongEqual]MZ^2+mh^2*)*)


(* ::Input:: *)
(*MLO2[t_]:=Evaluate[MLO2aux]*)


(* ::Input:: *)
(*\[Sigma]Integrated=1/(16 Pi s^2)*Integrate[MLO2[t],t];*)


(* ::Input:: *)
(*(*K\[ADoubleDot]ll\[EAcute]n function*)*)
(* \[Kappa][x_,y_,z_]:=Sqrt[x^2+y^2+z^2-2x y-2x z-2y z];*)
(*(*Kinematic limits for the integration in terms of \[Kappa]\[LongEqual]Sqrt[\[Lambda]], see definition of the Mandelstam variable t and take cos\[Theta]\[LongEqual]\[PlusMinus]1*)*)
(*tUpper[s_]:=1/2 (MZ^2+mh^2-s+\[Kappa][s,MZ^2,mh^2]);*)
(*tLower[s_]:=1/2 (MZ^2+mh^2-s-\[Kappa][s,MZ^2,mh^2]);*)


(* ::Input:: *)
(*(*LO cross section*)*)
(*\[Sigma]LOaux=(\[Sigma]Integrated/.t->tUpper[s])- (\[Sigma]Integrated/.t->tLower[s])//FullSimplify;*)


(* ::Input:: *)
(*(*Eq.(3.2) in the paper, typo in paper, missing square in the \[Kappa] inside the parenthesis*)*)


(* ::Input:: *)
(*\[Sigma]LOliterature=(EL^4 (ga^2+gv^2) \[Kappa][s,MZ^2,mh^2] (12 MZ^2 s+\[Kappa][s,MZ^2,mh^2]^2))/(192 CW^2 \[Pi] (MZ^2-s)^2 s^2 SW^2);*)


(* ::Input:: *)
(*(*Validation*)*)


(* ::Input:: *)
(*\[Sigma]LOaux==\[Sigma]LOliterature//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]& *)


(* ::Section:: *)
(*NLO*)


(* ::Subsection:: *)
(*Amplitude*)


(* ::Input:: *)
(*LoopTopo=CreateTopologies[1,2->2,ExcludeTopologies->{WFCorrections,AllBoxes,SelfEnergies}];*)


(* ::Input:: *)
(*LoopDiag=InsertFields[LoopTopo,{-F[2,{1}],F[2,{1}]}->{V[2],S[1]},Model->{"tripletSM"},GenericModel->{"tripletSM"},InsertionLevel->{Particles},ExcludeParticles->{S[1],S[2],V[3|4],U[1|2|3|4|11|12|31|32],F[_,{_}]}];*)


(* ::Input:: *)
(*Paint[LoopDiag,ColumnsXRows->{3,2},Numbering->Simple,ImageSize->{512,256}];*)


(* ::Input:: *)
(*(*PreFactor is set as -I so that ampsNLO[0]=Subsuperscript[M, Z^*\[Rule]Zh, vert] instead of iSubsuperscript[M, Z^*\[Rule]Zh, vert]*)*)


(* ::Input:: *)
(*ampsNLO[0]=FCFAConvert[CreateFeynAmp[LoopDiag,PreFactor->-I],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},LoopMomenta->{q},TransversePolarizationVectors->{k1},ChangeDimension->D,UndoChiralSplittings->True,List->True,DropSumOver->False,LorentzIndexNames->{\[Micro],\[Nu],\[Micro]1,\[Nu]1,\[Micro]2,\[Nu]2},FinalSubstitutions->{lam3->\[Lambda]3,MH0->MH,MHch->MH}]/.gc11->FCGV["EL"]/.gc38->(FCGV["CW"]*FCGV["EL"])/FCGV["SW"]/.gc47->-FCGV["EL"]//DiracSimplify[#,DiracSubstitute67->False]&//FCGVToSymbol;*)


(* ::Input:: *)
(*(*Splitting the amplitudes into ZZ and \[Gamma]Z terms to simplify the analysis*)*)


(* ::Input:: *)
(*amps\[Gamma]ZNLO[1]={ ampsNLO[0][[1]],ampsNLO[0][[2]],ampsNLO[0][[5]]};*)


(* ::Input:: *)
(*ampsZZNLO[1]={ ampsNLO[0][[3]],ampsNLO[0][[4]],ampsNLO[0][[6]]};*)


(* ::Input:: *)
(*(*Total of each contribution*)*)


(* ::Input:: *)
(*amps\[Gamma]ZNLO[2] =Total[amps\[Gamma]ZNLO[1]];*)


(* ::Input:: *)
(*ampsZZNLO[2] =Total[ampsZZNLO[1]];*)


(* ::Subsection:: *)
(*e + e - Z vertex contribution*)


(* ::Subsubsection:: *)
(*Definitions*)


(* ::Input:: *)
(*(*The LO amplitude without taking ME->0*)*)


(* ::Input:: *)
(*AmpLO=ampLO[0]//DiracSimplify[#,DiracSubstitute67->False]&;*)


(* ::Input:: *)
(*(*The NLO amplitude without taking ME->0, written in PaVe form*)*)


(* ::Input:: *)
(*AmpNLOZZ=TID[ampsZZNLO[2],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*(*The following command allow us to use the Breitenlohner-Maison-\[CloseCurlyQuote]t Hooft-Veltman scheme when dealing with matrices in D dimensions, in particular to handle gamma5*)*)


(* ::Input:: *)
(*FCSetDiracGammaScheme["BMHV"];*)


(* ::Input:: *)
(*(*Product aiming to get Eq. 3.8 in D dimension, notice 2\[GothicCapitalR][AB^*]\[LongEqual]AB^*+BA^**)*)


(* ::Input:: *)
(*ampsquaredZZ=(AmpNLOZZ ComplexConjugate[AmpLO]+AmpLO ComplexConjugate[AmpNLOZZ])/.gc54L->(gv+ga) EL/.gc54R->(gv-ga) EL/.gc44->(EL^2*v)/(2*CW^2*SW^2)//FeynAmpDenominatorExplicit//FermionSpinSum[#,ExtraFactor->1/2^2]& //DoPolarizationSums[#,k1]&;*)


(* ::Input:: *)
(*(*D dim -> 4 dim, take ME->0*)*)
(*(*Notice: //PaVeLimitTo4[(1/(2Pi)^4)*#]& is equal to //FCReplaceD[(1/(2Pi)^D)*#,D->4-2Epsilon]&//Series[#,{Epsilon,0,0}]&//Normal*)*)
(*(*Notice: ME->0 should be taken after DiracSimplify because it leads to new terms with ME*)*)


(* ::Input:: *)
(*TwoReampZZ=1/(2Pi)^4(PaVeLimitTo4[ampsquaredZZ]//DiracSimplify)/.ME->0//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&;*)


(* ::Input:: *)
(*(*The resultant expression is given in terms of tensor PaVe functions*)*)


(* ::Input:: *)
(*TwoReampZZ//FullSimplify//TraditionalForm*)


(* ::Input:: *)
(*(*We can look at the PaVe part alone by take out the factor of TwoReampZZ*)*)


(* ::Input:: *)
(*factorPavepartZZaux=1/(16 (MZ^2)  (\[Pi]^2) (SW^4) ((MZ^2-s)^2) ) EL^6 \[Lambda]3 v^2 (ga^2+gv^2);*)


(* ::Input:: *)
(*(*Writing down the Pavepart in terms of scalars PaVe functions A0,B0,C0 with the //PaVeReduce command:*)*)


(* ::Input:: *)
(*PavepartZZaux=(factorPavepartZZaux )^-1 TwoReampZZ/.u->MZ^2+mh^2-s-t//PaVeReduce;*)


(* ::Input:: *)
(*(*Use //PaVeLimitTo4 again because //PaVeReduce works in D dim. The 1/(2Pi)^4 is set before*)*)


(* ::Input:: *)
(*PavepartZZ=PaVeLimitTo4[PavepartZZaux]//Contract;*)


(* ::Subsubsection:: *)
(*Literature expression from Eq . (3.8) and C .8 in the paper*)


(* ::Input:: *)
(*(*Let's do the general case*)*)


(* ::Input:: *)
(*C22s=(1/(2*(Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))^2))*(2*Mz^2*(Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))-6*Mz^2*(Mz^2*(Mh^2-Mz^2+s)+(M1^2-M2^2)*(Mh^2+Mz^2-s))*PaVe[0,{MZ^2},{MH^2,MH^2}]-(Mh^6+(Mz^2-s)^2*(3*Mz^2-s)-Mh^4*(5*Mz^2+3*s)+Mh^2*(Mz^4+3*s^2)-2*(M1^2-M2^2)*(Mh^4+Mh^2*(4*Mz^2-s)+(Mz^2-s)^2))*PaVe[0,{mh^2},{MH^2,MH^2}]+((Mh^2-3*Mz^2)*(Mh^2-Mz^2)^2-(3*Mh^4+Mz^4)*s+(3*Mh^2+5*Mz^2)*s^2-s^3+(1/s)*(M1^2-M2^2)*((Mh^2-Mz^2)^3-5*(Mh^4-Mz^4)*s+(7*Mh^2-Mz^2)*s^2-3*s^3))*PaVe[0,{s},{MH^2,MH^2}]+2*(Mz^4*(Mh^4-2*Mh^2*(Mz^2-2*s)+(Mz^2-s)^2(*this square is missing in the paper*))+(M1^4+M2^4)*(Mh^4+2*Mh^2*(2*Mz^2-s)+(Mz^2-s)^2)-2*Mz^2*(Mh^4*(-2*M1^2+M2^2)+(Mz^2-s)^2*(M1^2-2*M2^2)+Mh^2*(Mz^2+s)*(M1^2+M2^2))-2*M1^2*M2^2*(Mh^4+2*Mh^2*(2*Mz^2-s)+(Mz^2-s)^2))*PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]);*)


(* ::Input:: *)
(*C23s=(1/(2*(Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))^2))*((Mh^2+Mz^2-s)*(Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))+(Mz^2*(-Mh^4+5*(Mz^2-s)^2-4*Mh^2*(Mz^2+s))-(M1^2-M2^2)*(Mh^4+(Mz^2-s)^2+2*Mh^2*(5*Mz^2-s)))*PaVe[0,{MZ^2},{MH^2,MH^2}]-((Mh^6+2*(Mz^2-s)^3-2*Mh^4*(3*Mz^2+2*s)+Mh^2*(3*Mz^4-8*Mz^2*s+5*s^2))-6*(M1^2-M2^2)*Mh^2*(Mh^2+Mz^2-s))*PaVe[0,{mh^2},{MH^2,MH^2}]+((Mh^6-(Mz^2-s)^2*(3*Mz^2+2*s)-Mh^4*(5*Mz^2+4*s)+Mh^2*(7*Mz^4-4*Mz^2*s+5*s^2))+(M1^2-M2^2)*((Mh^2-Mz^2)^2+2*(-3*Mh^4+Mh^2*Mz^2+2*Mz^4)*s+(3*Mh^2-5*Mz^2)*s^2+2*s^3))*PaVe[0,{s},{MH^2,MH^2}]+2*((Mz^2*((Mz^2-s)^3+Mh^4*(Mz^2+2*s)-Mh^2*(2*Mz^4-3*Mz^2*s+(*there is a wrong extra factor of 2 here in the paper*)s^2)))+3*(M1^4+M2^4)*Mh^2*(Mh^2+Mz^2-s)+M1^2*(Mh^2-Mz^2+s)*(Mh^4+2*Mh^2*(2*Mz^2-s)+(Mz^2-s)^2)+2*M2^2*((Mz^2-s)^3-Mh^4*(2*Mz^2+s)+Mh^2*(Mz^4-3*Mz^2*s+2*s^2))-6*M1^2*M2^2*Mh^2*(Mh^2+Mz^2-s))*PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]);*)


(* ::Input:: *)
(*C24s=(1/(4*(Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))))*((Mh^4+(Mz^2-s)^2-2*Mh^2*(Mz^2+s))-(Mz^2*(Mh^2-Mz^2+s)+(M1^2-M2^2)*(Mh^2+Mz^2-s))*PaVe[0,{MZ^2},{MH^2,MH^2}]+Mh^2*(Mh^2-Mz^2-s+2*M1^2-2*M2^2)*PaVe[0,{mh^2},{MH^2,MH^2}]+(s*(-Mh^2-Mz^2+s)-(M1^2-M2^2)*(Mh^2-Mz^2+s))*PaVe[0,{s},{MH^2,MH^2}]+(2*(Mh^2*Mz^2*s+(M1^4+M2^4)*Mh^2+M1^2*Mh^2*(Mh^2-Mz^2-s)+M2^2*((Mz^2-s)^2-Mh^2*(Mz^2+s))-2*M1^2*M2^2*Mh^2))*PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]);*)


(* ::Input:: *)
(*(*The factor of part from eq 3.8 with only PaVe functions, including spin average 1/4 and the sum factor 2*)*)


(* ::Input:: *)
(*factorPavepartliterature=-(1/4)*2*1/(16\[Pi]^2) (4EL^3)/(MZ SW CW) (gv^2+ga^2)/(s-MZ^2)^2 (-lam3*v)(- gc38)^2/.lam3->\[Lambda]3/.gc38->(CW*EL)/SW/.v->(2 SW CW MZ/EL)//FullSimplify;*)


(* ::Input:: *)
(*(*Part from eq 3.8 with only PaVe functions*)*)


(* ::Input:: *)
(*Pavepartliterature=2((t u +2s MZ^2 -mh^2MZ^2)(2C24s-1/2 PaVe[0,{mh^2},{MH^2,MH^2}])+(s+MZ^2-mh^2)(mh^2MZ^2-t u)(C22s-C23s))/.Mz->MZ/.Mh->mh/.M1->MH/.M2->MH;*)


(* ::Subsubsection:: *)
(*Validation with the literature*)


(* ::Input:: *)
(*(*Factor matching*)*)


(* ::Input:: *)
(*factorPavepartliterature==factorPavepartZZaux/.v->(2 SW CW MZ/EL)//Simplify*)


(* ::Input:: *)
(*(*B0(s) matching *)*)


(* ::Input:: *)
(*Coefficient[Pavepartliterature, PaVe[0,{s},{MH^2,MH^2}]]==Coefficient[PavepartZZ, PaVe[0,{s},{MH^2,MH^2}]]//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&//Simplify*)


(* ::Input:: *)
(*(*B0(k2) matching*)*)


(* ::Input:: *)
(*Coefficient[Pavepartliterature, PaVe[0,{mh^2},{MH^2,MH^2}]]==Coefficient[PavepartZZ, PaVe[0,{mh^2},{MH^2,MH^2}]]//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&*)


(* ::Input:: *)
(*(*B0(k1) matching*)*)


(* ::Input:: *)
(*Coefficient[Pavepartliterature, PaVe[0,{MZ^2},{MH^2,MH^2}]]==Coefficient[PavepartZZ, PaVe[0,{MZ^2},{MH^2,MH^2}]]//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&*)


(* ::Input:: *)
(*(*C0(k1,k2) matching *)*)


(* ::Input:: *)
(*Coefficient[Pavepartliterature,PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]]==Coefficient[PavepartZZ,PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]]//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&*)


(* ::Input:: *)
(*(*Independent term from the literature*)*)


(* ::Input:: *)
(*Indtermliterature= Pavepartliterature-(Coefficient[Pavepartliterature, PaVe[0,{s},{MH^2,MH^2}]]* PaVe[0,{s},{MH^2,MH^2}]+Coefficient[Pavepartliterature, PaVe[0,{MZ^2},{MH^2,MH^2}]]* PaVe[0,{MZ^2},{MH^2,MH^2}]+Coefficient[Pavepartliterature, PaVe[0,{mh^2},{MH^2,MH^2}]]* PaVe[0,{mh^2},{MH^2,MH^2}]+Coefficient[Pavepartliterature,PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]]*PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]);*)


(* ::Input:: *)
(*(*Independent term from the code*)*)


(* ::Input:: *)
(*Indterm= PavepartZZ-(Coefficient[PavepartZZ,PaVe[0,{s},{MH^2,MH^2}]]* PaVe[0,{s},{MH^2,MH^2}]+Coefficient[PavepartZZ,PaVe[0,{MZ^2},{MH^2,MH^2}]]* PaVe[0,{MZ^2},{MH^2,MH^2}]+Coefficient[PavepartZZ,PaVe[0,{mh^2},{MH^2,MH^2}]]* PaVe[0,{mh^2},{MH^2,MH^2}]+Coefficient[PavepartZZ,PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]]*PaVe[0,{mh^2,MZ^2,s},{MH^2,MH^2,MH^2}]);*)


(* ::Input:: *)
(*(*Independent term matching *)*)


(* ::Input:: *)
(*Indtermliterature==Indterm//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&*)


(* ::Subsection:: *)
(*e + e - \[Gamma] vertex contribution*)


(* ::Input:: *)
(*AmpNLO\[Gamma]Z=TID[amps\[Gamma]ZNLO[2],q,UsePaVeBasis->True,ToPaVe->True]//DiracSimplify;*)


(* ::Input:: *)
(*(*Product aiming to get Eq. 3.8 in D dimension*)*)


(* ::Input:: *)
(*ampsquared\[Gamma]Z=(AmpNLO\[Gamma]Z ComplexConjugate[AmpLO]+AmpLO ComplexConjugate[AmpNLO\[Gamma]Z])/.gc54L->(gv+ga) EL/.gc54R->(gv-ga) EL/.gc44->(EL^2*v)/(2*CW^2*SW^2)//FeynAmpDenominatorExplicit//FermionSpinSum[#,ExtraFactor->1/2^2]& //DoPolarizationSums[#,k1]&;*)


(* ::Input:: *)
(*(*D dim -> 4 dim, take ME->0*)*)


(* ::Input:: *)
(*TwoReamp\[Gamma]Z=(1/(2Pi)^4)*(PaVeLimitTo4[ampsquared\[Gamma]Z]//DiracSimplify)/.ME->0//TrickMandelstam[#,{s,t,u,MZ^2+mh^2}]&;*)


(* ::Input:: *)
(*(*The resultant expression is given in terms of tensor PaVe functions, this validates the earlier identification in the paper of the PaVe part of the expression factorizing from both contributions*)*)


(* ::Input:: *)
(*TwoReamp\[Gamma]Z//FullSimplify//TraditionalForm*)


(* ::Subsection:: *)
(*Renormalized propagator*)


(* ::Input:: *)
(*MWHat2=MW^2+Re[\[CapitalSigma]WWTHatForMW2];*)
(*MZHat2=MZ^2+Re[\[CapitalSigma]ZZTHatForMZ2];*)
(*cHat=Sqrt[MWHat2/MZHat2];*)
(*sHat=Sqrt[1-cHat^2];*)
(*ELHat=EL(1-1/2 \[Delta]Z\[Gamma]\[Gamma]Hat-1/2 sHat/cHat \[Delta]ZZ\[Gamma]Hat);*)
(*IWe3=-1/2;*)
(*Qe=-1;*)
(*gveff=(IWe3-2\[Kappa]Hat sHat^2 Qe)/(2 sHat cHat);*)
(*gaeff=IWe3/(2 sHat cHat);*)
(*\[Kappa]Hat=1-cHat/sHat \[CapitalSigma]\[Gamma]ZTHatFors/s;*)
(*\[CapitalGamma]ZMZ=Im[\[CapitalSigma]ZZTHatFors];*)
(*\[Delta]ZZZHat=-Re[D\[CapitalSigma]ZZTHatforMZ2];*)
(*\[Delta]ZhHat=-Re[D\[CapitalSigma]HHatFormh2];*)
(*\[Delta]Z\[Gamma]\[Gamma]Hat=-Re[D\[CapitalSigma]\[Gamma]\[Gamma]THatFor0];*)
(*\[Delta]ZZ\[Gamma]Hat=2 \[CapitalSigma]\[Gamma]ZTHatFor0/MZ^2;*)


(* ::Input:: *)
(*\[Rho]Hat=1/(s-MZ^2+I \[CapitalGamma]Z*MZ) (1+(Re[\[CapitalSigma]ZZTHatForMZ2]-Re[\[CapitalSigma]ZZTHatFors])/(s-MZ^2)+1/2 \[Delta]ZZZHat-1/2 Re[D\[CapitalSigma]HHatFormh2])/.v->(2 SW CW MZ/EL);*)


(* ::Subsection:: *)
(*Matrix Elements and cross section*)


(* ::Input:: *)
(*(*Including \[CapitalGamma]Z in the LO cross section*)*)


(* ::Input:: *)
(*\[Sigma]LOSMaux=Numerator[\[Sigma]LOaux]/(Denominator[\[Sigma]LOaux]/.(MZ^2-s)^2->((MZ^2-s)^2+MZ^2 \[CapitalGamma]Z^2));*)


(* ::Input:: *)
(*\[Sigma]LOSM[s_]:= Evaluate[\[Sigma]LOSMaux]*)


(* ::Input:: *)
(*(*NLO matrix element squared from self energies, typo in paper esq 3.1 and 3.7, missing factor of 1/2*)*)


(* ::Input:: *)
(*Mself2=(Norm[ELHat]^4MZHat2(Norm[gveff]^2+Norm[gaeff]^2))/(2MZ^2Norm[sHat]^2 Norm[cHat]^2) (t u+2s MZ^2-mh^2MZ^2)Norm[\[Rho]Hat]^2/.ScaleMu->\[Mu]/.u->(MZ^2+mh^2-s-t);*)


(* ::Input:: *)
(*MTree2=N[Mself2];*)


(* ::Input:: *)
(*(*NLO matrix element squared from intereference term with vertex contribution*)*)


(* ::Input:: *)
(*(*Adding the Z boson width*)*)


(* ::Input:: *)
(*TwoReampZZaux=(Numerator[TwoReampZZ])/((Denominator[TwoReampZZ])/.(MZ^2-s)^2->((MZ^2-s)^2+MZ^2 \[CapitalGamma]Z^2));*)


(* ::Input:: *)
(*TwoReamp\[Gamma]Zaux=(Numerator[TwoReamp\[Gamma]Z](MZ^2-s))/((Denominator[TwoReamp\[Gamma]Z]*(MZ^2-s))/.(MZ^2-s)^2->((MZ^2-s)^2+MZ^2 \[CapitalGamma]Z^2));*)


(* ::Input:: *)
(*MvertexMtree=Re[ComplexExpand[N[PaXEvaluate[(TwoReampZZaux+TwoReamp\[Gamma]Zaux)/.v->(2 SW CW MZ/EL)/.u->(MZ^2+mh^2-s-t) ,PaXC0Expand->True]]]];*)


(* ::Input:: *)
(*(*The differential cross section is given by d\[Sigma]NLO/dt=(1/(16 Pi s^2))|Mtot,corr|^2 with Mtot,corr given by Eq.3.6*)*)


(* ::Input:: *)
(*Mtotcorr2[\[Lambda]3_,MH_,\[Mu]_,s_,t_]:=Evaluate[(MTree2+MvertexMtree)-MLO2aux]*)


(* ::Input:: *)
(*\[Sigma]1loop[\[Lambda]3_,MH_,\[Mu]_,s_]:=(1/(16 Pi s^2)*NIntegrate[Mtotcorr2[\[Lambda]3,MH,\[Mu],s,tp],{tp,tLower[s],tUpper[s]}])*)


(* ::Section:: *)
(*Plots*)


(* ::Subsection:: *)
(*Constants from PDG 2025*)


(* ::Input:: *)
(*Clear[mh,MZ,MW,CW,SW,\[Alpha],EL,gv,ga,\[CapitalGamma]Z,convfactor]*)


(* ::Input:: *)
(*mh=Rationalize[125.2];*)
(*MZ=Rationalize[91.1880];*)
(*MW=Rationalize[80.3692];*)
(*CW=Rationalize[MW/MZ];*)
(*SW=Sqrt[1-CW^2];*)
(*\[Alpha]=Rationalize[1/137.036];*)
(*EL=Sqrt[4Pi \[Alpha]];*)
(*gv=Rationalize[(-1/2+2SW^2)/(2CW*SW)];*)
(*ga=Rationalize[(-1/(4CW*SW))];*)
(*\[CapitalGamma]Z=Rationalize[2.4952];*)
(*convfactor=0.3894*10^12;(*from GeV^-2 to fb*)*)


(* ::Subsection:: *)
(*Plots*)


(* ::Input:: *)
(*Plot[convfactor \[Sigma]LOSM[s^2],{s,200,400},PlotStyle->{Red},Frame->True,GridLines->Automatic,FrameLabel->{Style["\!\(\*SqrtBox[\(s\)]\)(GeV)",FontSize->18,FontFamily->"Times"],Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(LO\)]\)(fb)",FontSize->18,FontFamily->"Times"]},PlotRange->All,PlotLegends->Placed[{Style["\!\(\*SuperscriptBox[\(e\), \(+\)]\)\!\(\*SuperscriptBox[\(e\), \(-\)]\)\[Rule] Z h",FontFamily->"Times",FontSize->12]},{Right,Top}]]*)


(* ::Input:: *)
(*Manipulate[Abs[\[Sigma]1loop[\[Lambda]3,MH,\[Mu],s]/\[Sigma]LOSM[s]],{s,240,400},{MH,200,700},{\[Lambda]3,-5,5},{\[Mu],MZ,10^4}]*)


(* ::Input:: *)
(*(*ContourPlot takes around 35 minutes to evaluate this plot on a MacBook Pro 2.2 GHz 6-Core Intel Core i7*)*)


(* ::Input:: *)
(*ContourPlot[Abs[\[Sigma]1loop[\[Lambda]3,MH,MZ,240]]/\[Sigma]LOSM[240],{MH,200,700},{\[Lambda]3,-5,5},Contours->{0.5/100},ContourStyle->{{Blue,Dashed}},ContourShading->None,FrameTicksStyle->Directive[Bold,12],FrameLabel->{Style["\!\(\*SubscriptBox[\(M\), \(\[Sum]\)]\)(GeV)",FontSize->18,FontFamily->"Times",Bold],Style["\!\(\*SubscriptBox[\(\[Lambda]\), \(3\)]\)",FontSize->18,FontFamily->"Times",Bold]},PlotLegends->Placed[{Style["|\!\(\*SubscriptBox[\(\[Delta]\[Sigma]\), \(Zh\)]\)|\[LessEqual]0.5%",FontFamily->"Times",FontSize->12]},{Left,Bottom}]]*)


(* ::Input:: *)
(*(*RegionPlot has a bug that does not allow NIntegrate to be evaluated properly, I already tried some solutions on stackexchange but none worked out, https://mathematica.stackexchange.com/questions/94402/error-with-nintegrate-and-regionplot*)*)


(* ::Input:: *)
(*(*I have not tried the Plot3D option here though, https://mathematica.stackexchange.com/questions/189717/regionplot-evaluating-nintegrate-before-assigning-variable-values-and-resulting*)*)


(* ::Input:: *)
(*(*RegionPlot[ImplicitRegion[(Abs[\[Sigma]1loop[\[Lambda]3,MH,MZ,240]]/\[Sigma]LO[240])\[LessEqual]0.005,{{MH,200,700},{\[Lambda]3,-5,5}}],BoundaryStyle\[Rule]{Black},PlotRange->{{200,700},{-5,5}},Frame->True,FrameTicksStyle->Directive[Bold,12],FrameLabel\[Rule]{Style["Subscript[M, \[Sum]](GeV)",FontSize\[Rule]18,FontFamily\[Rule]"Times",Bold],Style["Subscript[\[Lambda], 3]",FontSize\[Rule]18,FontFamily\[Rule]"Times",Bold]},PlotStyle->Directive[Opacity[0.5,Cyan]],PlotLegends\[Rule]Placed[{Style["|Subscript[\[Delta]\[Sigma], Zh]|\[LessEqual]0.5%",FontFamily\[Rule]"Times",FontSize\[Rule]12]},{Left,Bottom}]]*)*)
