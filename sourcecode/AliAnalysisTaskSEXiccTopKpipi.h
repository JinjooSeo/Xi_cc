#ifndef AliAnalysisTaskSEXiccTopKpipi_H
#define AliAnalysisTaskSEXiccTopKpipi_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>

#include "AliGenParam.h"
#include "AliGenerator.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "AliESDcascade.h"
#include "AliESDVertex.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskWeakDecayVertexer_mod.h"
#include "R5Detector.h"
#include "AliPIDResponse.h"

class TH1I;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;


class AliAnalysisTaskSEXiccTopKpipi : public AliAnalysisTaskSE
{
public:

    AliAnalysisTaskSEXiccTopKpipi();
    AliAnalysisTaskSEXiccTopKpipi(const char *name);
    virtual ~AliAnalysisTaskSEXiccTopKpipi();

    /// Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    R5Detector* CreateDetector();
    void RecoEvent();

    TString GetGenerator(Int_t label, AliGenCocktailEventHeader* header);

    void MakeCandidates();
    void PrepareTracks();
    Bool_t IsSelected(Int_t CutFlag, AliESDtrack *trk);
    Bool_t FillHistoXic();
    Bool_t FillHistoXicc();
    AliESDVertex* CallReconstructSecondaryVertexXic(Double_t &dispersion);
    AliESDVertex* CallReconstructSecondaryVertexXicc(Double_t &dispersion);

    AliExternalTrackParam* GetTrackParamForCascade();
    AliESDVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t &dispersion);

    TTree* BuildTreeXiccGen(TString name, TString title);
    TTree* BuildTreeCutVar(TString name, TString title, bool isforomegac, bool isforomegacc, bool isforomegaccc);
    TTree* BuildTreeRecoTracks(TString name, TString title);
    TTree* BuildTreeReconstruction(TString name, TString title);

    void FillGenTree();
    void FillRecoTrackTree();

    void SetIsSignalProd(bool opt) {fIsMCSignalProd=opt;}

    void FillXiccHistogram(AliESDtrack *proton, AliESDtrack *kaon, AliESDtrack *pion, AliESDtrack *softpion);
    void FillXiccTree(AliESDtrack *proton, AliESDtrack *kaon, AliESDtrack *pion, AliESDtrack *softpion);
    void DefineTree();


private:

    AliAnalysisTaskSEXiccTopKpipi(const AliAnalysisTaskSEXiccTopKpipi&);
    AliAnalysisTaskSEXiccTopKpipi& operator=(const AliAnalysisTaskSEXiccTopKpipi&);


    AliMCEvent              *fMcEvent;              //!<! MC event
    AliInputEventHandler    *fMcHandler;            //!<! MCEventHandler
    AliStack                *fStack;                //!<!
    AliESDEvent             *fEvent;                //!<!
    Float_t                 fBzkG;
    AliESDVertex            *fPrimVtx;              //!<!
    Float_t                 fEtaCut;
    Int_t                   fEvtCount;
    TH1F                    *fNentries;             //!<!   histogram with number of events on output slot 1
    TH1F                    *fNcounters;            //!<!   histogram with number of events, gen/filtered charmed baryons
    TList                   *fOutput;            //!<!   list on output slot 2
    TList                   *fOutputGen;            //!<!   list on output slot 2
    TList                   *fOutputReco;           //!<!   list on output slot 2
    R5Detector              *fITS;                  //!<!
  //  AliPIDResponse          *fPIDResponse;

  Bool_t fIsMCSignalProd;
  vector<int> fNHitsAssigned;
vector<int> fParticlePdg;
vector<int> fMotherPdg;

    Int_t                  fProtonCuts;
    Int_t                  fKaonCuts;
    Int_t                  fPionCuts;
    Int_t                  fSoftPionCuts;

    TArrayI                 *fProtonTrackArray;
    TArrayI                 *fKaonTrackArray;
    TArrayI                 *fPionTrackArray;
    TArrayI                 *fSoftPionTrackArray;

    TTree                   *fTree;
    Float_t                 *fTreeVariable = nullptr; //!

    THnSparseF              *fhSparsePx; //!
    THnSparseF              *fhSparsePy; //!
    THnSparseF              *fhSparsePz; //!
    THnSparseF              *fhSparsePT; //!
    THnSparseF              *fhSparseM; //!
    THnSparseF              *fhSparseY; //!


    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEXiccTopKpipi,1);
    /// \endcond
};

#endif
