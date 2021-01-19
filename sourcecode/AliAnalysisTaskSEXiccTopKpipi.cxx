/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// \authors: Jinjoo Seo
// \email: jseo@cern.ch
//
////////////////////////////////////////////////////////////

#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TChain.h"
#include "TMath.h"
#include <THnSparse.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include "TObjArray.h"
#include "TArrayI.h"
#include <TClonesArray.h>
#include <TObjectTable.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
#include "AliGenCocktailEventHeader.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenerator.h"
#include "AliVertexingHFUtils.h"
#include "AliVertexerTracks.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>

#include "AliPIDResponse.h"

#include "AliAnalysisTaskSEXiccTopKpipi.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXiccTopKpipi);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEXiccTopKpipi::AliAnalysisTaskSEXiccTopKpipi():
AliAnalysisTaskSE(),
fMcEvent(0x0),
fMcHandler(0x0),
fStack(0x0),
fEvent(0x0),
fBzkG(0.),
fPrimVtx(0x0),
fEtaCut(2.),
fEvtCount(0),
fNentries(0x0),
fNcounters(0x0),
fOutputGen(0x0),
fOutputReco(0x0),
fITS(0x0),
fKaonCuts(0),
fPionCuts(0),
fSoftPionCuts(0),
fProtonTrackArray(0x0),
fKaonTrackArray(0x0),
fPionTrackArray(0x0),
fSoftPionTrackArray(0x0),
fhSparsePx(0x0),
fhSparsePy(0x0),
fhSparsePz(0x0),
fhSparsePT(0x0),
fhSparseY(0x0),
fhSparseM(0x0),
nProton(0),
nKaon(0),
nPion(0),
nSoftPion(0),
fGenXiccTree(0x0),
fGenXiccTreeVariable(0),
fRecoXiccTree(0x0),
fRecoXiccTreeVariable(0),
fRecoTrackTree(0x0),
fRecoTrackTreeVariable(0)
{
    /// Default constructor
}
//________________________________________________________________________
AliAnalysisTaskSEXiccTopKpipi::AliAnalysisTaskSEXiccTopKpipi(const char *name):
AliAnalysisTaskSE(name),
fMcEvent(0x0),
fMcHandler(0x0),
fStack(0x0),
fEvent(0x0),
fBzkG(0.),
fPrimVtx(0x0),
fEtaCut(2.),
fEvtCount(0),
fNentries(0x0),
fNcounters(0x0),
fOutputGen(0x0),
fOutputReco(0x0),
fITS(0x0),
fKaonCuts(0),
fPionCuts(0),
fSoftPionCuts(0),
fProtonTrackArray(0x0),
fKaonTrackArray(0x0),
fPionTrackArray(0x0),
fSoftPionTrackArray(0x0),
fhSparsePx(0x0),
fhSparsePy(0x0),
fhSparsePz(0x0),
fhSparsePT(0x0),
fhSparseY(0x0),
fhSparseM(0x0),
nProton(0),
nKaon(0),
nPion(0),
nSoftPion(0),
fGenXiccTree(0x0),
fGenXiccTreeVariable(0),
fRecoXiccTree(0x0),
fRecoXiccTreeVariable(0),
fRecoTrackTree(0x0),
fRecoTrackTreeVariable(0)
{
    /// Default constructor

    DefineInput(0, TChain::Class());
// Output slot #1 writes into a TH1F container (number of events)
  DefineOutput(1,TH1F::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TTree::Class());
  DefineOutput(4,TTree::Class());
  DefineOutput(5,TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskSEXiccTopKpipi::~AliAnalysisTaskSEXiccTopKpipi()
{
    if (fMcEvent){
        delete fMcEvent;
        fMcEvent = 0x0;
    }
    if (fMcHandler){
        delete fMcHandler;
        fMcHandler = 0x0;
    }
    if (fStack){
        delete fStack;
        fStack = 0x0;
    }
    if (fEvent){
        delete fEvent;
        fEvent = 0x0;
    }
    if (fPrimVtx){
        delete fPrimVtx;
        fPrimVtx = 0x0;
    }
    if (fNentries){
        delete fNentries;
        fNentries = 0x0;
    }
    if(fNcounters){
        delete fNcounters;
        fNcounters = 0x0;
    }
    if (fOutputGen){
        delete fOutputGen;
        fOutputGen = 0x0;
    }
    if (fOutputReco){
        delete fOutputReco;
        fOutputReco = 0x0;
    }
    if (fITS){
        delete fITS;
        fITS = 0x0;
    }
    if (fProtonTrackArray){
        delete fProtonTrackArray;
        fProtonTrackArray = 0x0;
    }
    if (fKaonTrackArray){
        delete fKaonTrackArray;
        fKaonTrackArray = 0x0;
    }
    if (fPionTrackArray){
        delete fPionTrackArray;
        fPionTrackArray = 0x0;
    }
    if (fSoftPionTrackArray){
        delete fSoftPionTrackArray;
        fSoftPionTrackArray = 0x0;
    }
    if (fhSparsePx){
        delete fhSparsePx;
        fhSparsePx = 0x0;
    }
    if (fhSparsePy){
        delete fhSparsePy;
        fhSparsePy = 0x0;
    }
    if (fhSparsePz){
        delete fhSparsePz;
        fhSparsePz = 0x0;
    }
    if (fhSparsePT){
        delete fhSparsePT;
        fhSparsePT = 0x0;
    }
    if (fhSparseM){
        delete fhSparseM;
        fhSparseM = 0x0;
    }
    if (fhSparseY){
        delete fhSparseY;
        fhSparseY = 0x0;
    }
    if (fGenXiccTree){
        delete fGenXiccTree;
        fGenXiccTree = 0x0;
    }
    if (fRecoXiccTree){
        delete fRecoXiccTree;
        fRecoXiccTree = 0x0;
    }
    if (fRecoTrackTree){
        delete fRecoTrackTree;
        fRecoTrackTree = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::Init()
{
    /// Initialization

    //if(fDebug > 1) printf("AliAnalysisTaskSEXiccTopKpipi::Init() \n");

}

//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::UserCreateOutputObjects()
{

    /// Create the output container
    //
    //if(fDebug > 1) printf("AliAnalysisTaskSEXiccTopKpipi::UserCreateOutputObjects() \n");

    fEvent = new AliESDEvent();
    fEvent->CreateStdContent();

    DefineGenXiccTree();
    DefineRecoXiccTree();
    DefineRecoTrackTree();

    // create new ITS detector
    fITS = CreateDetector();

    const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();

    fNentries=new TH1F(nameoutput, "Number of events", 23,-0.5,22.5);
fNentries->GetXaxis()->SetBinLabel(1,"All gen events");
fNentries->GetXaxis()->SetBinLabel(2,"All gen particle");
fNentries->GetXaxis()->SetBinLabel(3,"Physical particle");
fNentries->GetXaxis()->SetBinLabel(4,"Mother label not -1");
fNentries->GetXaxis()->SetBinLabel(5,"Excluding charge=0");
fNentries->GetXaxis()->SetBinLabel(6,"After Kine cuts");
fNentries->GetXaxis()->SetBinLabel(7,"Reco tracks");
fNentries->GetXaxis()->SetBinLabel(8,"Not Reco tracks");
fNentries->GetXaxis()->SetBinLabel(9,"N #pi from #Lambda gen");
fNentries->GetXaxis()->SetBinLabel(10,"N p from #Lambda gen");
fNentries->GetXaxis()->SetBinLabel(11,"N #Lambda gen");
fNentries->GetXaxis()->SetBinLabel(12,"N #Omega gen");
fNentries->GetXaxis()->SetBinLabel(13,"N #Lambda from #Omega gen");
fNentries->GetXaxis()->SetBinLabel(14,"N K from #Omega gen");
fNentries->GetXaxis()->SetBinLabel(15,"N D gen");
fNentries->GetXaxis()->SetBinLabel(16,"N B gen");

fNcounters=new TH1F(nameoutput, "Counters used for normalization", 8,-0.5,7.5);
fNcounters->GetXaxis()->SetBinLabel(1,"Gen events");
fNcounters->GetXaxis()->SetBinLabel(2,"All gen particle");
fNcounters->GetXaxis()->SetBinLabel(3,"Gen #Omega^{ccc}");
fNcounters->GetXaxis()->SetBinLabel(4,"Reco #Omega^{ccc}");
fNcounters->GetXaxis()->SetBinLabel(5,"Gen #Omega^{cc}");
fNcounters->GetXaxis()->SetBinLabel(6,"Reco #Omega^{cc}");
fNcounters->GetXaxis()->SetBinLabel(7,"Gen #Omega^{c}");
fNcounters->GetXaxis()->SetBinLabel(8,"Reco #Omega^{c}");


    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("listOutput");

    fProtonTrackArray = new TArrayI(1000);
    fKaonTrackArray = new TArrayI(1000);
    fPionTrackArray = new TArrayI(2000);
    fSoftPionTrackArray = new TArrayI(2000);

	fProtonCuts = 1;
	fKaonCuts = 2;
	fPionCuts = 3;
	fSoftPionCuts = 4;

    Int_t nbinsP[6]={40,40,40,40,40,40};
    Double_t lowEdgesP[6]={0.,0.,0.,0.,0.,0.};
    Double_t upEdgesP[6]={20.,20.,20.,20.,20.,20.};
    fhSparsePx = new THnSparseF("fhSparsePx","fhSparsePx;px_proton:px_kaon:px_pion:px_spion:px_Xic:px_Xicc",6,nbinsP,lowEdgesP,upEdgesP);
    fhSparsePy = new THnSparseF("fhSparsePy","fhSparsePy;py_proton:py_kaon:py_pion:py_spion:py_Xic:py_Xicc",6,nbinsP,lowEdgesP,upEdgesP);
    fhSparsePz = new THnSparseF("fhSparsePz","fhSparsePz;pz_proton:pz_kaon:pz_pion:pz_spion:pz_Xic:pz_Xicc",6,nbinsP,lowEdgesP,upEdgesP);
    fhSparsePT = new THnSparseF("fhSparsePT","fhSparsePT;pT_proton:pT_kaon:pT_pion:pT_spion:pT_Xic:pT_Xicc",6,nbinsP,lowEdgesP,upEdgesP);

    Int_t nbinsM[6] = {120,120,120,120,120,120};
    Double_t lowEdgesM[6] = {0.89,0.45,0.09,0.09,2.0,3.2};
    Double_t upEdgesM[6] = {0.97,0.53,0.17,0.17,2.8,4.0};
    fhSparseM = new THnSparseF("fhSparseM","fhSparseM;m_proton:m_kaon:m_pion:m_spion:m_Xic:m_Xicc",6,nbinsM,lowEdgesM,upEdgesM);

    Int_t nbinsY[4] = {100,100,100,100};
    Double_t lowEdgesY[4] = {-1.,-1,-1,-1};
    Double_t upEdgesY[4] = {1,1,1,1};
    fhSparseY = new THnSparseF("fhSparseY","fhSparseY;Y_proton:Y_kaon:Y_pion:Y_spion:Y_Xic:Y_Xicc",4,nbinsY,lowEdgesY,upEdgesY);

    fOutput->Add(fhSparsePx);
    fOutput->Add(fhSparsePy);
    fOutput->Add(fhSparsePz);
    fOutput->Add(fhSparsePT);
    fOutput->Add(fhSparseM);
    fOutput->Add(fhSparseY);

    // Post the data
    PostData(1,fNentries);
    PostData(2,fOutput);
    PostData(3,fGenXiccTree);
    PostData(4,fRecoXiccTree);
    PostData(5,fRecoTrackTree);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::UserExec(Option_t */*option*/)

{
    /// Execute analysis for current event:
    fEvtCount++;
 //   cout<< "UserExec " << fEvtCount << endl;
    // Process MC truth

    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
        Printf("ERROR: Could not retrieve MC event handler");
        return;
    }
    fMcEvent = eventHandler->MCEvent();
    if (!fMcEvent) {
        Printf("ERROR: Could not retrieve MC event");
        return;
    }

    // create dummy vertex for beam diamond
    const double kBeamSig = 3e-4;//50e-4;
    double beamPos[3] = {0.,0.,0.};
    double beamSig[3] = {kBeamSig, kBeamSig, kBeamSig };
    AliESDVertex diamond(beamPos,beamSig,"diamond");

    // access the stack
    fStack = ((AliMCEvent*)fMcEvent)->Stack();
    if(!fStack) {Printf("ERROR: Could not retrieve the Stack"); return;}

    // reset ESD event
    fEvent->Reset();
		//cout << "Num Of Tracks Before: " << fEvent->GetNumberOfTracks()<< endl;
    fEvent->SetMagneticField(fITS->GetBField()*10.);
    fEvent->SetDiamond(&diamond);
    // fill ESD event
    RecoEvent();
		//cout << "Num Of Tracks After: " << fEvent->GetNumberOfTracks()<< endl;
    fBzkG = fEvent->GetMagneticField();
    fPrimVtx   = (AliESDVertex*)fEvent->GetPrimaryVertexTracks();
    Double_t sigma[3]={0.,0.,0.};
    Double_t position[3]={0.,0.,0.};
    fPrimVtx->GetSigmaXYZ(sigma);
    fPrimVtx->GetXYZ(position);

    MakeCandidates();

    return;
}
//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::Terminate(Option_t */*option*/)
{
    /// Terminate analysis
    //
    //if(fDebug > 1) printf("AliAnalysisTaskSEXiccTopKpipi: Terminate() \n");

    Printf("Nevents analysed = %d", fEvtCount);
    return;

}
//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::MakeCandidates(){

  PrepareTracks();
  FillRecoTrackTree();

  for(int i=0; i<nProton; i++){ //for(int i=0; i<nProton; i++){
    AliESDtrack* fProtonTrack = (AliESDtrack*)fEvent->GetTrack(fProtonTrackArray->At(i));

    for(int j=0; j<nKaon; j++){ //for(int j=0; j<nKaon; j++){
      AliESDtrack* fKaonTrack = (AliESDtrack*)fEvent->GetTrack(fKaonTrackArray->At(j));

      if(fKaonTrack->Charge()*fProtonTrack->Charge()>=0) continue;
      if(fKaonTrack->GetID()==fProtonTrack->GetID()) continue;

      for(int k=0; k<nPion; k++){ //for(int k=0; k<nPion; k++){
        AliESDtrack* fPionTrack = (AliESDtrack*)fEvent->GetTrack(fPionTrackArray->At(k));

        if(fPionTrack->Charge()*fKaonTrack->Charge()>=0) continue;
        if(fPionTrack->Charge()*fProtonTrack->Charge()<=0) continue;
        if(fPionTrack->GetID()==fKaonTrack->GetID()) continue;
        if(fPionTrack->GetID()==fProtonTrack->GetID()) continue;

		      for(int l=0; l<nPion; l++){ //for(int k=0; k<nPion; k++){
            AliESDtrack* fSoftPionTrack = (AliESDtrack*)fEvent->GetTrack(fSoftPionTrackArray->At(l));
		        if(fSoftPionTrack->Charge()*fProtonTrack->Charge()<=0) continue;
		        if(fSoftPionTrack->Charge()*fKaonTrack->Charge()>=0) continue;
            if(fSoftPionTrack->Charge()*fPionTrack->Charge()<=0) continue;
		        if(fSoftPionTrack->GetID()==fKaonTrack->GetID()) continue;
            if(fSoftPionTrack->GetID()==fProtonTrack->GetID()) continue;
		        if(fSoftPionTrack->GetID()==fPionTrack->GetID()) continue;
		 // cout << "FILL!!----------------------------------------------------------------" << endl;
            FillRecoXiccTree(fSoftPionTrack,fProtonTrack,fKaonTrack,fPionTrack);
            FillXiccHistogram(fProtonTrack,fKaonTrack,fPionTrack,fSoftPionTrack);
          //FillXiccTree(fProtonTrack,fKaonTrack,fPionTrack,fSoftPionTrack);
        }//l
      }//k
    }//j
  }//i
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::PrepareTracks(){
  AliESDtrack* trk;
  Int_t nentr = (Int_t)fEvent->GetNumberOfTracks();
  fProtonTrackArray->Reset();
  fKaonTrackArray->Reset();
  fPionTrackArray->Reset();
  fSoftPionTrackArray->Reset();
    nProton = 0;
    nKaon = 0;
    nPion = 0;
    nSoftPion = 0;

  for(Int_t i=0; i<nentr; i++){
    trk = (AliESDtrack*)fEvent->GetTrack(i);
	//cout << "PID : " <<trk->GetPID() << endl;
	//	trk->Print("");
    if(IsSelected(fProtonCuts,trk)) {
      fProtonTrackArray->AddAt(i,nProton);
      nProton++;
    }
    if(IsSelected(fKaonCuts,trk)) {
      fKaonTrackArray->AddAt(i,nKaon);
      nKaon++;
    }
    if(IsSelected(fPionCuts,trk)) {
      fPionTrackArray->AddAt(i,nPion);
      nPion++;
    }
    if(IsSelected(fSoftPionCuts,trk)) {
      fSoftPionTrackArray->AddAt(i,nSoftPion);
      nSoftPion++;
    }
     //AddAt (Int_t c, Int_t i) : Add Int_t c at position i. Check for out of bounds.
  }
//cout << "proton kaon pion spion : " << nProton << "  " << nKaon << "  " << nPion << "  "<< nSoftPion << endl;

  return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEXiccTopKpipi::IsSelected(Int_t CutFlag, AliESDtrack *trk){

  if(!trk) return kFALSE;
  if(trk->Pt()<0.5) return kFALSE;

  if(CutFlag==1){ //Proton candidate
    //if(fPIDResponse->GetNumberOfSigmasTPC(trk,AliPID::kProton)>4) return kFALSE;
    if(TMath::Abs(trk->GetMass()-TDatabasePDG::Instance()->GetParticle(2212)->Mass())>0.05) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==2){ //Kaon candidate
    //if(fPIDResponse->GetNumberOfSigmasTPC(trk,AliPID::kKaon)>4) return kFALSE;
    if(TMath::Abs(trk->GetMass()-TDatabasePDG::Instance()->GetParticle(321)->Mass())>0.05) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==3){ //Pion candidate
    //if(fPIDResponse->GetNumberOfSigmasTPC(trk,AliPID::kPion)>2) return kFALSE;  //to seperate pion and electron
    if(TMath::Abs(trk->GetMass()-TDatabasePDG::Instance()->GetParticle(211)->Mass())>0.05) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==4){ //Soft pion candidate
    //if(fPIDResponse->GetNumberOfSigmasTPC(trk,AliPID::kPion)>2) return kFALSE;  //to seperate pion and electron
    if(TMath::Abs(trk->GetMass()-TDatabasePDG::Instance()->GetParticle(211)->Mass())>0.05) return kFALSE;
    return kTRUE;
  }
  else{
    cout << "Wrong PID cut flag!" << endl;
    return kFALSE;
  }
}

Double_t AliAnalysisTaskSEXiccTopKpipi::GetCorrectedMass(AliESDtrack *trk){
  Double_t m = trk->M();
  Double_t P[3];
  trk->GetPxPyPz(P);
  Double_t px = P[0];
  Double_t py = P[1];
  Double_t pz = P[2];
  Double_t E = trk->E();
  Double_t CorrectedM = sqrt(E*E-px*px-py*py-pz*pz);

  return CorrectedM;
}

//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::FillXiccHistogram(AliESDtrack *proton, AliESDtrack *kaon, AliESDtrack *pion, AliESDtrack *softpion){

  Double_t PionP[3], KaonP[3], ProtonP[3], SPionP[3];
  proton->GetPxPyPz(ProtonP);
  kaon->GetPxPyPz(KaonP);
  pion->GetPxPyPz(PionP);
  softpion->GetPxPyPz(SPionP);
  Double_t px_proton = ProtonP[0];
  Double_t py_proton = ProtonP[1];
  Double_t pz_proton = ProtonP[2];
  Double_t pT_proton = sqrt(pow(px_proton,2)+pow(py_proton,2));
  Double_t px_kaon = KaonP[0];
  Double_t py_kaon = KaonP[1];
  Double_t pz_kaon = KaonP[2];
  Double_t pT_kaon = sqrt(pow(px_kaon,2)+pow(py_kaon,2));
  Double_t px_pion = PionP[0];
  Double_t py_pion = PionP[1];
  Double_t pz_pion = PionP[2];
  Double_t pT_pion = sqrt(pow(PionP[0],2)+pow(PionP[1],2));
  Double_t px_spion = SPionP[0];
  Double_t py_spion = SPionP[1];
  Double_t pz_spion = SPionP[2];
  Double_t pT_spion = sqrt(pow(SPionP[0],2)+pow(SPionP[1],2));
  Double_t m_proton = GetCorrectedMass(proton);
  Double_t m_kaon = GetCorrectedMass(kaon);
  Double_t m_pion = GetCorrectedMass(pion);
  Double_t m_spion = GetCorrectedMass(softpion);
  Double_t E_proton = proton->E();
  Double_t E_kaon = kaon->E();
  Double_t E_pion = pion->E();
  Double_t E_spion = softpion->E();
  Double_t Y_proton = proton->Y();
  Double_t Y_kaon = kaon->Y();
  Double_t Y_pion = pion->Y();
  Double_t Y_spion = softpion->Y();
  Double_t px_Xic = sqrt(pow(ProtonP[0]+KaonP[0]+PionP[0],2));
  Double_t py_Xic = sqrt(pow(ProtonP[1]+KaonP[1]+PionP[1],2));
  Double_t pz_Xic = sqrt(pow(ProtonP[2]+KaonP[2]+PionP[2],2));
  Double_t pT_Xic = sqrt(pow(px_Xic,2)+pow(py_Xic,2));
  Double_t E_Xic = sqrt(pow(px_Xic,2)+pow(py_Xic,2)+pow(pz_Xic,2)+2.46794*2.46794);
  Double_t m_Xic = sqrt(pow(E_proton+E_kaon+E_pion,2)-pow(px_Xic,2)-pow(py_Xic,2)-pow(pz_Xic,2));
  Double_t px_Xicc = sqrt(pow(px_Xic+SPionP[0],2));
  Double_t py_Xicc = sqrt(pow(py_Xic+SPionP[1],2));
  Double_t pz_Xicc = sqrt(pow(pz_Xic+SPionP[2],2));
  Double_t pT_Xicc = sqrt(pow(px_Xicc,2)+pow(py_Xicc,2));
  Double_t E_Xicc = sqrt(pow(px_Xicc,2)+pow(py_Xicc,2)+pow(pz_Xicc,2)+3.6212*3.6212);
  Double_t m_Xicc = sqrt(pow(E_Xic+E_spion,2)-pow(px_Xicc,2)-pow(py_Xicc,2)-pow(pz_Xicc,2));
  Double_t px_Merge[6] = {px_proton,px_kaon,px_pion,px_spion,px_Xic,px_Xicc}; fhSparsePx->Fill(px_Merge);
  Double_t py_Merge[6] = {py_proton,py_kaon,py_pion,py_spion,py_Xic,py_Xicc}; fhSparsePy->Fill(py_Merge);
  Double_t pz_Merge[6] = {pz_proton,pz_kaon,pz_pion,pz_spion,pz_Xic,pz_Xicc}; fhSparsePz->Fill(pz_Merge);
  Double_t pT_Merge[6] = {pT_proton,pT_kaon,pT_pion,pT_spion,pT_Xic,pT_Xicc}; fhSparsePT->Fill(pT_Merge);
  Double_t m_Merge[6] = {m_proton,m_kaon,m_pion,m_spion,m_Xic,m_Xicc}; fhSparseM->Fill(m_Merge);
  Double_t Y_Merge[4] = {Y_proton,Y_kaon,Y_pion,Y_spion}; fhSparseY->Fill(Y_Merge);
  return;
}

void AliAnalysisTaskSEXiccTopKpipi::DefineGenXiccTree(){
  fGenXiccTree = new TTree("fGenXiccTree","fGenXiccTree");
  vector <TString> fTreeVariableName = {"px_Xicc","py_Xicc","pz_Xicc","mass_Xicc","pT_Xicc","y_Xicc","Phi_Xicc","Rxy_Xicc",
    "px_Xic","py_Xic","pz_Xic","mass_Xic","pT_Xic","y_Xic","Phi_Xic","Rxy_Xic",
    "px_spion","py_spion","pz_spion","mass_spion","pT_spion","y_spion","Phi_spion","Rxy_spion",
    "px_proton","py_proton","pz_proton","mass_proton","pT_proton","y_proton","Phi_proton","Rxy_proton",
    "px_kaon","py_kaon","pz_kaon","mass_kaon","pT_kaon","y_kaon","Phi_kaon","Rxy_kaon",
    "px_pion","py_pion","pz_pion","mass_pion","pT_pion","y_pion","Phi_pion","Rxy_pion",
    "PA_spionXic","SecondaryVertex","DecayLengthXY","Chi2perNDF"};
  fGenXiccTreeVariable = new Float_t [fTreeVariableName.size()];

  for (Int_t ivar=0; ivar<(Float_t)fTreeVariableName.size(); ivar++) fGenXiccTree->Branch(fTreeVariableName[ivar].Data(),&fGenXiccTreeVariable[ivar],Form("%s/f",fTreeVariableName[ivar].Data()));

  return;
}

void AliAnalysisTaskSEXiccTopKpipi::FillGenXiccTree(TParticle* mcXicc){
  if(TMath::Abs(mcXicc->GetPdgCode())!=4422) return;

  TParticle* mcXic;
  TParticle* mcspion;
  TParticle* mcproton;
  TParticle* mckaon;
  TParticle* mcpion;

  for(int i=0; i<mcXicc->GetNDaughters(); i++){
    TParticle* mcpart = fMcEvent->Particle(mcXicc->GetDaughter(i));
    if(mcpart->GetPdgCode()==211) mcspion = mcpart;
    if(mcpart->GetPdgCode()==4232){
      mcXic = mcpart;
      for(int j=0; j<mcXic->GetNDaughters(); j++){
        TParticle* mcpart2 = fMcEvent->Particle(mcXic->GetDaughter(i));
        if(mcpart2->GetPdgCode()==2212) mcproton = mcpart2;
        if(mcpart2->GetPdgCode()==2321) mckaon = mcpart2;
        if(mcpart2->GetPdgCode()==211) mcpion = mcpart2;
      }
    }
  }

  for(int i=0; i<44; i++) fGenXiccTreeVariable[i] = -9999;
  fGenXiccTreeVariable[0] = mcXicc->Px();
  fGenXiccTreeVariable[1] = mcXicc->Py();
  fGenXiccTreeVariable[2] = mcXicc->Pz();
  fGenXiccTreeVariable[3] = mcXicc->GetMass();
  fGenXiccTreeVariable[4] = mcXicc->Pt();
  fGenXiccTreeVariable[5] = mcXicc->Y();
  fGenXiccTreeVariable[6] = mcXicc->Phi();
  fGenXiccTreeVariable[7] = mcXicc->R();
  fGenXiccTreeVariable[8] = mcXic->Px();
  fGenXiccTreeVariable[9] = mcXic->Py();
  fGenXiccTreeVariable[10] = mcXic->Pz();
  fGenXiccTreeVariable[11] = mcXic->GetMass();
  fGenXiccTreeVariable[12] = mcXic->Pt();
  fGenXiccTreeVariable[13] = mcXic->Y();
  fGenXiccTreeVariable[14] = mcXic->Phi();
  fGenXiccTreeVariable[15] = mcXic->R();
  fGenXiccTreeVariable[16] = mcspion->Px();
  fGenXiccTreeVariable[17] = mcspion->Py();
  fGenXiccTreeVariable[18] = mcspion->Pz();
  fGenXiccTreeVariable[19] = mcspion->GetMass();
  fGenXiccTreeVariable[20] = mcspion->Pt();
  fGenXiccTreeVariable[21] = mcspion->Y();
  fGenXiccTreeVariable[22] = mcspion->Phi();
  fGenXiccTreeVariable[23] = mcspion->R();
  fGenXiccTreeVariable[24] = mcproton->Px();
  fGenXiccTreeVariable[25] = mcproton->Py();
  fGenXiccTreeVariable[26] = mcproton->Pz();
  fGenXiccTreeVariable[27] = mcproton->GetMass();
  fGenXiccTreeVariable[28] = mcproton->Pt();
  fGenXiccTreeVariable[29] = mcproton->Y();
  fGenXiccTreeVariable[30] = mcproton->Phi();
  fGenXiccTreeVariable[31] = mcproton->R();
  fGenXiccTreeVariable[32] = mckaon->Px();
  fGenXiccTreeVariable[33] = mckaon->Py();
  fGenXiccTreeVariable[34] = mckaon->Pz();
  fGenXiccTreeVariable[35] = mckaon->GetMass();
  fGenXiccTreeVariable[36] = mckaon->Pt();
  fGenXiccTreeVariable[37] = mckaon->Y();
  fGenXiccTreeVariable[38] = mckaon->Phi();
  fGenXiccTreeVariable[39] = mckaon->R();

  fGenXiccTreeVariable[40] = mckaon->Pt();
  fGenXiccTreeVariable[41] = mckaon->Y();
  fGenXiccTreeVariable[42] = mckaon->Phi();
  fGenXiccTreeVariable[43] = mckaon->R();


  fGenXiccTree->Fill();
}


void AliAnalysisTaskSEXiccTopKpipi::DefineRecoXiccTree(){
  fRecoXiccTree = new TTree("fRecoXiccTree","fRecoXiccTree");
  vector <TString> fTreeVariableName = {"px_Xicc","py_Xicc","pz_Xicc","mass_Xicc","pT_Xicc","y_Xicc","Phi_Xicc","Rxy_Xicc",
    "px_Xic","py_Xic","pz_Xic","mass_Xic","pT_Xic","y_Xic","Phi_Xic","Rxy_Xic",
    "px_spion","py_spion","pz_spion","mass_spion","pT_spion","y_spion","Phi_spion","Rxy_spion",
    "px_proton","py_proton","pz_proton","mass_proton","pT_proton","y_proton","Phi_proton","Rxy_proton",
    "px_kaon","py_kaon","pz_kaon","mass_kaon","pT_kaon","y_kaon","Phi_kaon","Rxy_kaon",
    "px_pion","py_pion","pz_pion","mass_pion","pT_pion","y_pion","Phi_pion","Rxy_pion",
    "PA_spionXic","SecondaryVertex","DecayLengthXY","Chi2perNDF"};
  fRecoXiccTreeVariable = new Float_t [fTreeVariableName.size()];

  for (Int_t ivar=0; ivar<(Float_t)fTreeVariableName.size(); ivar++) fRecoXiccTree->Branch(fTreeVariableName[ivar].Data(),&fRecoXiccTreeVariable[ivar],Form("%s/f",fTreeVariableName[ivar].Data()));

  return;
}

void AliAnalysisTaskSEXiccTopKpipi::FillRecoXiccTree(AliESDtrack* spion, AliESDtrack* proton, AliESDtrack* kaon, AliESDtrack* pion){
  Int_t labspion = spion->GetLabel(); if(labspion==-1) return;
  Int_t labproton = proton->GetLabel(); if(labproton==-1) return;
  Int_t labkaon = kaon->GetLabel(); if(labkaon==-1) return;
  Int_t labpion = pion->GetLabel(); if(labpion==-1) return;

  TParticle* mcspion = fMcEvent->Particle(labspion);
  TParticle* mcproton = fMcEvent->Particle(labproton);
  TParticle* mckaon = fMcEvent->Particle(labkaon);
  TParticle* mcpion = fMcEvent->Particle(labpion);

  if(TMath::Abs(mcspion->GetPdgCode())!=211) return;
  if(TMath::Abs(mcproton->GetPdgCode())!=2212) return;
  if(TMath::Abs(mckaon->GetPdgCode())!=321) return;
  if(TMath::Abs(mcpion->GetPdgCode())!=211) return;

  if(mcproton->GetFirstMother() != mckaon->GetFirstMother()) return;
  if(mcproton->GetFirstMother() != mcpion->GetFirstMother()) return;
  if(mcproton->GetFirstMother()==-1) return;

  TParticle* mcXic = fMcEvent->Particle(mcproton->GetFirstMother());
  if(TMath::Abs(mcXic->GetPdgCode())!=4232) return;
  if(mcXic->GetFirstMother() != mcspion->GetFirstMother()) return;
  if(mcXic->GetFirstMother()==-1) return;

  TParticle* mcXicc = fMcEvent->Particle(mcXic->GetFirstMother());
  if(TMath::Abs(mcXicc->GetPdgCode())!=4422) return;

  for(int i=0; i<44; i++) fRecoXiccTreeVariable[i] = -9999;
  fRecoXiccTreeVariable[0] = mcXicc->Px();
  fRecoXiccTreeVariable[1] = mcXicc->Py();
  fRecoXiccTreeVariable[2] = mcXicc->Pz();
  fRecoXiccTreeVariable[3] = mcXicc->GetMass();
  fRecoXiccTreeVariable[4] = mcXicc->Pt();
  fRecoXiccTreeVariable[5] = mcXicc->Y();
  fRecoXiccTreeVariable[6] = mcXicc->Phi();
  fRecoXiccTreeVariable[7] = mcXicc->R();
  fRecoXiccTreeVariable[8] = mcXic->Px();
  fRecoXiccTreeVariable[9] = mcXic->Py();
  fRecoXiccTreeVariable[10] = mcXic->Pz();
  fRecoXiccTreeVariable[11] = mcXic->GetMass();
  fRecoXiccTreeVariable[12] = mcXic->Pt();
  fRecoXiccTreeVariable[13] = mcXic->Y();
  fRecoXiccTreeVariable[14] = mcXic->Phi();
  fRecoXiccTreeVariable[15] = mcXic->R();
  fRecoXiccTreeVariable[16] = mcspion->Px();
  fRecoXiccTreeVariable[17] = mcspion->Py();
  fRecoXiccTreeVariable[18] = mcspion->Pz();
  fRecoXiccTreeVariable[19] = mcspion->GetMass();
  fRecoXiccTreeVariable[20] = mcspion->Pt();
  fRecoXiccTreeVariable[21] = mcspion->Y();
  fRecoXiccTreeVariable[22] = mcspion->Phi();
  fRecoXiccTreeVariable[23] = mcspion->R();
  fRecoXiccTreeVariable[24] = mcproton->Px();
  fRecoXiccTreeVariable[25] = mcproton->Py();
  fRecoXiccTreeVariable[26] = mcproton->Pz();
  fRecoXiccTreeVariable[27] = mcproton->GetMass();
  fRecoXiccTreeVariable[28] = mcproton->Pt();
  fRecoXiccTreeVariable[29] = mcproton->Y();
  fRecoXiccTreeVariable[30] = mcproton->Phi();
  fRecoXiccTreeVariable[31] = mcproton->R();
  fRecoXiccTreeVariable[32] = mckaon->Px();
  fRecoXiccTreeVariable[33] = mckaon->Py();
  fRecoXiccTreeVariable[34] = mckaon->Pz();
  fRecoXiccTreeVariable[35] = mckaon->GetMass();
  fRecoXiccTreeVariable[36] = mckaon->Pt();
  fRecoXiccTreeVariable[37] = mckaon->Y();
  fRecoXiccTreeVariable[38] = mckaon->Phi();
  fRecoXiccTreeVariable[39] = mckaon->R();

  fRecoXiccTreeVariable[40] = mckaon->Pt();
  fRecoXiccTreeVariable[41] = mckaon->Y();
  fRecoXiccTreeVariable[42] = mckaon->Phi();
  fRecoXiccTreeVariable[43] = mckaon->R();


  fRecoXiccTree->Fill();
}

void AliAnalysisTaskSEXiccTopKpipi::DefineRecoTrackTree(){
  fRecoTrackTree = new TTree("fRecoTrackTree","fRecoTrackTree");
  vector <TString> fTreeVariableName = {"recoPt","recoPhi","recoEta","recoY","recoNhits","genPt","genEta","genRad","recoMClabel","pdgCodePart",
    "pdgCodeMother","recoNReconstructableHits","recoMClabelMother","genPtMoth"};
  fRecoTrackTreeVariable = new Float_t [fTreeVariableName.size()];

  for (Int_t ivar=0; ivar<(Float_t)fTreeVariableName.size(); ivar++) fRecoTrackTree->Branch(fTreeVariableName[ivar].Data(),&fRecoTrackTreeVariable[ivar],Form("%s/f",fTreeVariableName[ivar].Data()));

  return;
}
void AliAnalysisTaskSEXiccTopKpipi::FillRecoTrackTree(){
    int ntr = fEvent->GetNumberOfTracks();
    AliESDtrack* esdtr=0x0;
    TParticle *mcPart=0x0;
    for (int itr=0;itr<ntr;itr++) {
        fNentries->Fill(6);
        esdtr = fEvent->GetTrack(itr);
        Float_t recoPt=esdtr->Pt();
        Float_t recoEta=esdtr->Eta();
        Float_t recoY=esdtr->Y();
        Float_t recoPhi=esdtr->Phi();
        Int_t recoHits = (Int_t)esdtr->GetTRDNclusterdEdx();
        Int_t lbl=esdtr->GetLabel();

        //corresponding generated particle
        mcPart=fMcEvent->Particle(lbl);
        Float_t genPt=mcPart->Pt();
        Float_t genEta=mcPart->Eta();
        Float_t genY=mcPart->Y();
        Float_t genRad=mcPart->R();
        Float_t genPtMoth=-1.;
        int pdgCodePart=mcPart->GetPdgCode();
        int lblMother = mcPart->GetFirstMother();
        int pdgCodeMother = -99999;
        if(lblMother>-1) {
            pdgCodeMother = (fMcEvent->Particle(lblMother))->GetPdgCode();
            genPtMoth = (fMcEvent->Particle(lblMother))->Pt();
        }

        const TBits &hits = esdtr->GetTPCClusterMap();
        //TBits &fakes = esdtr->GetTPCSharedMap(); // here we store fakes, but at the moment they are not set
        Int_t nhits=0;
        //printf(" Hits: ");
        for (int ilr=0;ilr<fITS->GetNActiveLayers();ilr++) {
            //printf("%c", hits.TestBitNumber(ilr) ? '+':'-');
            if(hits.TestBitNumber(ilr))nhits++;
        }

        for(int i=0; i<14; i++) fRecoTrackTreeVariable[i] = -9999;
        //Fill histograms
        fRecoTrackTreeVariable[0] = recoPt;
        fRecoTrackTreeVariable[1] = recoPhi;
        fRecoTrackTreeVariable[2] = recoEta;
        fRecoTrackTreeVariable[3] = recoY;
        fRecoTrackTreeVariable[4] = nhits;
        fRecoTrackTreeVariable[5] = genPt;
        fRecoTrackTreeVariable[6] = genEta;
        fRecoTrackTreeVariable[7] = genRad;
        fRecoTrackTreeVariable[8] = lbl;
        fRecoTrackTreeVariable[9] = pdgCodePart;
        fRecoTrackTreeVariable[10] = pdgCodeMother;
        fRecoTrackTreeVariable[11] = recoHits;
        fRecoTrackTreeVariable[12] = lblMother;
        fRecoTrackTreeVariable[13] = genPtMoth;

        fRecoTrackTree->Fill();
    }
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::RecoEvent(){


    int npart = fMcEvent->GetNumberOfTracks();//fStack->GetNprimary();
    const AliVVertex *vtxMC = fMcEvent->GetPrimaryVertex();
    Double_t covmatrix[6]={0.,0.,0.,0.,0.,0.};
    Double_t positiongen[3]={0.,0.,0.};
    //vtxMC->GetSigmaXYZ(sigmagen);
    vtxMC->GetXYZ(positiongen);


    TParticle *part = 0x0;
    AliESDtrack* esdTr = 0x0;
    fNHitsAssigned.clear();
    fParticlePdg.clear();
    fMotherPdg.clear();

    // loop on particles
    //Printf("N PARTICLES %d",npart);
    for(int i=0; i<npart; i++){

        Int_t pdgPart=-1;
        Int_t lblMoth=-1;
        Int_t pdgMoth=-1;
        Float_t ptPart=0.;
        Float_t pPart=0.;
        Float_t etaPart=0.;
        Float_t rapPart=0.;
        Float_t phiPart=0.;
        Float_t xPart=0.;
        Float_t yPart=0.;
        Float_t zPart=0.;
        Float_t rxyPart=0.;
        Int_t nhistassigned=-1;
        Int_t nhistreco=-1;
		Double_t massPart = 0.;
        bool isReco=false;


        //Printf("______");
        //std::unique_ptr<TParticle>part((TParticle*)fMcEvent->Particle(i)); //smart pointer for TParticle* part
        //part.get()
        part = (TParticle*)fMcEvent->Particle(i);//fStack->Particle(i);
        if(!part) {Printf("no part"); continue;}
        FillGenXiccTree(part);
        fNentries->Fill(1);
        fNcounters->Fill(1);
        pdgPart = part->GetPdgCode();
        lblMoth = part->GetFirstMother();
        if(lblMoth>-1)pdgMoth = ((TParticle*)fMcEvent->Particle(lblMoth))->GetPdgCode();
        TString genname = GetGenerator(i,(AliGenCocktailEventHeader*)fMcEvent->GenEventHeader());
        //if(lblMoth==-1){

            //if(genname.Contains("Hijing")) continue;
        if(genname.Contains("PYTHIA")) continue;
            //printf("lblMoth: %d; lblPart: %d; Gen: %s\n",lblMoth, i,genname.Data());
            //if(genname.Contains("PYTHIA"))nc_py++;
        //}
        ptPart=part->Pt();
        pPart=part->P();
        etaPart=part->Eta();
        rapPart=part->Y();
        phiPart=part->Phi();
        xPart=part->Vx();
        yPart=part->Vy();
        zPart=part->Vz();
        rxyPart=part->R();
		massPart=part->GetMass();
		//Printf("HERE i=%10d pdg=%5d mass=%+5.3e",i,pdgPart,massPart);
       // Printf("HERE i=%10d pdg=%5d lblmoth=%10d pt=%.5f, eta=%.5f, phi=%.5f, x=%.5f y=%.5f z=%.5f, rxy=%.5f",i,pdgPart, lblMoth, ptPart, etaPart, phiPart, xPart, yPart, zPart, rxyPart);
        Bool_t res = kFALSE;
        if(ptPart>0.){
            res = fITS->ProcessTrack(part);
            //res = fITS->ProcessTrack(part.get()); //part.get() to pass the pointer as argument
            //Printf("after processTrack i=%d  res=%d",i,res);
            nhistassigned = fITS->GetNHitsAssigned();
        }
        fNHitsAssigned.push_back(nhistassigned);
        fParticlePdg.push_back(pdgPart);
        fMotherPdg.push_back(lblMoth);
        //genname=GetGenerator(i,(AliGenCocktailEventHeader*)fMcEvent->GenEventHeader());
        //if(genname.Contains("Hijing")) continue;
        //if(genname.Contains("PYTHIA")) continue;
        //printf("lblMoth: %d; lblPart: %d; Gen: %s\n",lblMoth, i,genname.Data());


        //Printf("RECO EVENT: index %d | particle %d particleName=%s | mother %d | n reco hits %d",i,fParticlePdg.at(i),part->GetName(),fMotherPdg.at(i),fNHitsAssigned.at(i));
        /*
        Bool_t isPhysPrimary = fMcEvent->IsPhysicalPrimary(i);
        if(!isPhysPrimary) {Printf("i=%d PDG=%d %s no phys prim",i,part->GetPdgCode(),part->GetName()); continue;}
        if(part->GetStatusCode()<=0) {Printf("i=%d PDG=%d %s status code <0",i,part->GetPdgCode(),part->GetName());  continue;}
        fNentries->Fill(2);
        Int_t motherIndex=part->GetFirstMother();
        if(motherIndex==-1) {Printf("i=%d motherIndex==-1 PDG=%d %s",i,part->GetPdgCode(),part->GetName());  continue;}
        fNentries->Fill(3);
        int charge = TDatabasePDG::Instance()->GetParticle(part->GetPdgCode())->Charge()/3;
        if(charge==0) {Printf("i=%d charge=0 PDG=%d %s",i,part->GetPdgCode(),part->GetName());  continue;}
        fNentries->Fill(4);
        Float_t thisEta = part->Eta();
        //if(TMath::Abs(thisEta) > fEtaCut) continue; //fEtaCut
        fNentries->Fill(5);
        */
        Bool_t isPhysPrimary = fMcEvent->IsPhysicalPrimary(i);
        if(TMath::Abs(pdgPart)==4444)fNcounters->Fill(2);
        if(TMath::Abs(pdgPart)==4432)fNcounters->Fill(4);
        if(TMath::Abs(pdgPart)==4332)fNcounters->Fill(6);

        Int_t statusCode = part->GetStatusCode();
        if(isPhysPrimary && statusCode>0) {
            fNentries->Fill(2);
            if(lblMoth>-1) {
                fNentries->Fill(3);
                if( (TDatabasePDG::Instance()->GetParticle(part->GetPdgCode())->Charge()/3)!=0 ) {
                    fNentries->Fill(4);
                    if(TMath::Abs(part->Eta())<2.){
                        fNentries->Fill(5);
                    }
                }
            }
        }
        if(!isPhysPrimary) continue;
        if(statusCode<=0) continue;
        if(lblMoth==-1) continue;
        if(TDatabasePDG::Instance()->GetParticle(part->GetPdgCode())->Charge()/3==0 ) continue;
        // fast track reco
        //Printf("i track = %d part = %s",i,part->GetName());
        //Bool_t res = fITS->ProcessTrack(part);
        if(res){
            isReco=true;
            esdTr = (AliESDtrack*)fITS->GetProbeTrackInwardAsESDTrack();
            esdTr->SetStatus(AliESDtrack::kTPCin|AliESDtrack::kTPCout|AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
				esdTr->SetTRDntracklets((UChar_t)nhistassigned);
            esdTr->SetLabel(i);
			//	Printf("i track = %d part = %s",i,part->GetName());
		//cout << "Track Info : " << esdTr->GetPID() << endl;
            double dum = fEvent->AddTrack(esdTr);
			//	cout << "Track ID " << dum << endl;
				AliESDtrack *trkdum = (AliESDtrack*) fEvent->GetTrack(dum);
				//cout << "Track PID " << trkdum->GetPID() << endl;
            const TBits &hits = esdTr->GetTPCClusterMap();
            Int_t nhits=0;
            for (int ilr=0;ilr<fITS->GetNActiveLayers();ilr++) {
                //printf("%c", hits.TestBitNumber(ilr) ? '+':'-');
                if(hits.TestBitNumber(ilr))nhits++;
            }
            nhistreco=nhits;
            //Printf("nhits reco %d", nhits);
        }
        else {
            //fNentries->Fill(7);/* continue;*/
            isReco=false;
        }
    }// end loop over particles
    //AliESDUtils::RefitESDVertexTracks(esd);
    double posVtx[3]={0,0,0};
    double sigmaVtx[3]={0,0,0};
    double covMatrix[6]={sigmaVtx[0]*sigmaVtx[0],0.,sigmaVtx[1]*sigmaVtx[1],0.,0.,sigmaVtx[2]*sigmaVtx[2]};
    //AliESDVertex *myVtx = new AliESDVertex(posVtx,sigmaVtx,"Vertex");
    AliESDVertex myVtx(posVtx,sigmaVtx,"Vertex");
    fEvent->SetPrimaryVertexTracks(&myVtx);
}
//________________________________________________________________________
TString AliAnalysisTaskSEXiccTopKpipi::GetGenerator(Int_t label, AliGenCocktailEventHeader* header){

    Int_t nsumpart=0;
    TList *lh = header->GetHeaders();
    Int_t nh = lh->GetEntries();
    for(Int_t i=0; i<nh; i++){
        AliGenEventHeader *gh = (AliGenEventHeader*)lh->At(i);
        TString genname = gh->GetName();
        Int_t npart = gh->NProduced();
        if(label>=nsumpart && label<(nsumpart+npart)) return genname;
        nsumpart+=npart;
    }
    TString empty="";
    return empty;

}
//________________________________________________________________________
R5Detector* AliAnalysisTaskSEXiccTopKpipi::CreateDetector(){

    AliESDtrack::OnlineMode(kTRUE); // to avoid friend track creation
    R5Detector* det = new R5Detector("ALICE","ITS");
	//	R5Detector* det = new R5Detector();

    det->SetPropagateToOrigin(kTRUE); // if we want all tracks to be propagated to DCA to 0/0/0.

    det->SetBField(1.);
    // new ideal Pixel properties?
    //config_0
    Double_t x0IB     = 0.001;
    //config_1
    //Double_t x0IB     = 0.0015;
    //config_2 & config_3
    //Double_t x0IB     = 0.0004;

    Double_t x0OB     = 0.005;
    Double_t xRho     = 0.;
    Double_t resRPhiIB     = 0.0001;
    Double_t resZIB        = 0.0001;
    Double_t resRPhiOB     = 0.0005;
    Double_t resZOB        = 0.0005;
    Double_t eff           = 0.98;
    //
    // select Z span in such a way to have +-1 unit eta coverage for vertex at 2sigmaZ (~12cm) from nominal IP
    // i.e. Zmax >= 12 + R/tan( 2*atan(exp(-1.)) ) = 12 + R/0.851

    det->AddLayer((char*)"vertex",  0.0,  0.1, 0, 0); // dummy vertex for matrix calculation
    //config_0
    det->AddLayer((char*)"bpipe",   1.6,  200., 0.0022);
    //config_1 & config_2 & config_3
    //det->AddLayer((char*)"bpipe",   1.4,  200., 0.0015);
    /* |eta|<1
     det->AddLayer((char*)"ddd1",    1.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff);
     det->AddLayer((char*)"ddd2",    2.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff);
     det->AddLayer((char*)"ddd3",    3.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff);
     det->AddLayer((char*)"ddd3a",   8.0,  21.0, x0IB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"ddd4",   20.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"ddd5",   25.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff);

     //det->AddLayer((char*)"ddd6",  35.0, 80.0, x0OB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"ddd7",   40.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"ddd8",   55.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff);

     //  det->AddLayer((char*)"dddZ",  90., 130., x0OB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"dddY",   80.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff);
     det->AddLayer((char*)"dddX",  100.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff);
     */
    //|eta|<2
    //config_0
    det->AddLayer((char*)"ddd1",    1.8,  63.0, x0IB, xRho, resRPhiIB, resZIB,eff);
    //config_1 & config_2
    //det->AddLayer((char*)"ddd1",    1.5,  63.0, x0IB, xRho, resRPhiIB, resZIB,eff);
    //config_3
    //det->AddLayer((char*)"ddd1",    0.5,  63.0, x0IB, xRho, resRPhiIB, resZIB,eff);

    det->AddLayer((char*)"ddd2",    2.8,  63.0, x0IB, xRho, resRPhiIB, resZIB,eff);
    det->AddLayer((char*)"ddd3",    3.8,  63.0, x0IB, xRho, resRPhiIB, resZIB,eff);
    det->AddLayer((char*)"ddd3a",   8.0,  63.0, x0IB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"ddd4",   20.0,  126.0, x0OB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"ddd5",   25.0,  126.0, x0OB, xRho, resRPhiOB, resZOB,eff);

    //det->AddLayer((char*)"ddd6",  35.0, 80.0, x0OB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"ddd7",   40.0,  240.0, x0OB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"ddd8",   55.0,  240.0, x0OB, xRho, resRPhiOB, resZOB,eff);

    //  det->AddLayer((char*)"dddZ",  90., 130., x0OB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"dddY",   80.0, 390.0, x0OB, xRho, resRPhiOB, resZOB,eff);
    det->AddLayer((char*)"dddX",  100.0, 390.0, x0OB, xRho, resRPhiOB, resZOB,eff);

    //det->Print("opt");

    return det;


    //rescale second param by facto3 --> eta+-2
}
