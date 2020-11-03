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
//
//
// \authors: Andrea Festanti
// \email:   andrea.festanti@cern.ch
//
////////////////////////////////////////////////////////////

#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TChain.h"
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

#include "AliRecoDecayOmegac.h"
#include "AliRecoDecayOmegacc.h"
#include "AliRecoDecayOmegaccc.h"
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
fVertexerTracks(0x0),
fBzkG(0.),
fPrimVtx(0x0),
fEtaCut(2.),
fEvtCount(0),
fNentries(0x0),
fNcounters(0x0),
fOutputGen(0x0),
fOutputReco(0x0),
fITS(0x0),
fPIDResponse(0),
fProtonCandidates(1),
fKaonCuts(2),
fPionCuts(3),
fSoftPionCuts(4),
fProtonTrackArray(0x0),
fKaonTrackArray(0x0),
fPionTrackArray(0x0),
fSoftPionTrackArray(0x0)
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
fVertexerTracks(0x0),
fBzkG(0.),
fPrimVtx(0x0),
fEtaCut(2.),
fEvtCount(0),
fNentries(0x0),
fNcounters(0x0),
fOutputGen(0x0),
fOutputReco(0x0),
fITS(0x0),
fPIDResponse(0),
fProtonCandidates(1),
fKaonCuts(2),
fPionCuts(3),
fSoftPionCuts(4),
fProtonTrackArray(0x0),
fKaonTrackArray(0x0),
fPionTrackArray(0x0),
fSoftPionTrackArray(0x0)
{
    /// Default constructor
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
    if (fVertexerTracks){
        delete fVertexerTracks;
        fVertexerTracks = 0x0;
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

    // create new ITS detector
    fITS = CreateDetector();

    fVertexerTracks = new AliVertexerTracks();

    const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();

    // Post the data
    PostData(1,fNentries);
    PostData(2,fNcounters);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::UserExec(Option_t */*option*/)

{
    /// Execute analysis for current event:
    fEvtCount++;
    cout<< "UserExec " << fEvtCount << endl;
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
    fEvent->SetMagneticField(fITS->GetBField()*10.);
    fEvent->SetDiamond(&diamond);
    // fill ESD event
    RecoEvent();
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
  for(Int_t i=0; i<nentr; i++){
    fEsdTr1 = (AliESDtrack*)fEvent->GetTrack(i);
    fIdxTrack1 = fEsdTr1->GetID();
  }


    return;
}
//________________________________________________________________________
void AliAnalysisTaskSEXiccTopKpipi::PrepareTracks(){
  AliESDTrack* trk;
  AliESDTrack* fProtonCandidates;
  AliESDTrack* fKaonCandidates;
  AliESDTrack* fPionCandidates;
  AliESDTrack* fSoftPionCandidates;
  Int_t nProton = 0; Int_t nKaon = 0; Int_t nPion = 0; Int_t nSoftPion = 0;

  Int_t nentr = (Int_t)fEvent->GetNumberOfTracks();

  for(Int_t i=0; i<nentr; i++){
    trk = (AliESDtrack*)fEvent->GetTrack(i);
    if(IsSelected(fProtonCuts,trk)) {
      fProtonCandidates = trk;
      nProton++;
      fProtonTrackArray->AddAt(i,nProton);
    }
    if(IsSelected(fKaonCuts,trk)) {
      fKaonCandidates = trk;
      nKaon++;
      fKaonTrackArray->AddAt(i,nKaon);
    }
    if(IsSelected(fPionCuts,trk)) {
      fPionCandidates = trk;
      nPion++;
      fPionTrackArray->AddAt(i,nPion);
    }
    if(IsSelected(fSoftPionCuts,trk)) {
      fSoftPionCandidates = trk;
      nSoftPion++;
      fSoftPionTrackArray->AddAt(i,nSoftPion);
    }
     //AddAt (Int_t c, Int_t i) : Add Int_t c at position i. Check for out of bounds.
  }

  delete fProtonCandidates;
  delete fKaonCandidates;
  delete fPionCandidates;
  delete fSoftPionCandidates;

  return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEXiccTopKpipi::IsSelected(Int_t *CutFlag, AliESDTrack *trk){

  if(!trk) return kFALSE;

  if(CutFlag==1){ //Proton candidate
    if(fPIDResponse->GetNumberOfSigmasTPC(trk,kProton)>4) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==2){ //Kaon candidate
    if(fPIDResponse->GetNumberOfSigmasTPC(trk,kKaon)>4) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==3){ //Pion candidate
    if(fPIDResponse->GetNumberOfSigmasTPC(trk,kPion)>4) return kFALSE;
    return kTRUE;
  }
  else if(CutFlag==4){ //Soft pion candidate
    if(fPIDResponse->GetNumberOfSigmasTPC(trk,kPion)>4) return kFALSE;
    return kTRUE;
  }
  else{
    cout << "Wrong PID cut flag!" << endl;
    return kFALSE;
  }
}

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
    Printf("N PARTICLES %d",npart);
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
        bool isReco=false;


        //Printf("______");
        //std::unique_ptr<TParticle>part((TParticle*)fMcEvent->Particle(i)); //smart pointer for TParticle* part
        //part.get()
        part = (TParticle*)fMcEvent->Particle(i);//fStack->Particle(i);
        if(!part) {Printf("no part"); continue;}
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
        //Printf("HERE i=%10d pdg=%5d lblmoth=%10d pt=%.5f, eta=%.5f, phi=%.5f, x=%.5f y=%.5f z=%.5f, rxy=%.5f",i,pdgPart, lblMoth, ptPart, etaPart, phiPart, xPart, yPart, zPart, rxyPart);
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
        //Printf("i track = %d part = %s eta=%f motherindex =%d",i,part->GetName(),thisEta,motherIndex);
        //Bool_t res = fITS->ProcessTrack(part);
        if(res){
            isReco=true;
            esdTr = (AliESDtrack*)fITS->GetProbeTrackInwardAsESDTrack();
            esdTr->SetStatus(AliESDtrack::kTPCin|AliESDtrack::kTPCout|AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
            esdTr->SetTRDntracklets((UChar_t)nhistassigned);
            esdTr->SetLabel(i);
            fEvent->AddTrack(esdTr);
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

        //Printf("FILL Tree REconstruction: index=%d res=%d | %d %d | %f %f %f %f %f %f %f %f %f | %d %d %d",i,res,fParticleVarPdg,fParticleVarPdgMoth,fParticleVarPt,fParticleVarP,fParticleVarEta,fParticleVarRap,fParticleVarPhi,fParticleVarX,fParticleVarY,fParticleVarZ,fParticleVarRxy,fParticleVarNhitsAssigned,fParticleVarNhitsReco,fParticleVarIsReco);
        //if(fIsMCSignalProd && fIsFillReconstrucion) fTreeReconstruction->Fill();
//        esdTr = (AliESDtrack*)fITS->GetProbeTrackInwardAsESDTrack();
//        esdTr->SetStatus(AliESDtrack::kTPCin|AliESDtrack::kTPCout|AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
//        esdTr->SetTRDntracklets((UChar_t)nhistassigned);
//        esdTr->SetLabel(i);
//        fEvent->AddTrack(esdTr);



    }// end loop over particles
    //AliESDUtils::RefitESDVertexTracks(esd);
    double posVtx[3]={0,0,0};
    double sigmaVtx[3]={0,0,0};
    double covMatrix[6]={sigmaVtx[0]*sigmaVtx[0],0.,sigmaVtx[1]*sigmaVtx[1],0.,0.,sigmaVtx[2]*sigmaVtx[2]};
    //AliESDVertex *myVtx = new AliESDVertex(posVtx,sigmaVtx,"Vertex");
    AliESDVertex myVtx(posVtx,sigmaVtx,"Vertex");
    fEvent->SetPrimaryVertexTracks(&myVtx);

    printf("number of reco tracks: %d",fEvent->GetNumberOfTracks());

    //delete esdVtx;
    //delete myVtx;


}

//________________________________________________________________________
R5Detector* AliAnalysisTaskSEXiccTopKpipi::CreateDetector(){

    AliESDtrack::OnlineMode(kTRUE); // to avoid friend track creation
    R5Detector* det = new R5Detector("ALICE","ITS");

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
