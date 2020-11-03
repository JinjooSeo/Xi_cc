#if !(defined(__CLING__)  || defined(__CINT__)) || defined(__ROOTCLING__) || defined(__ROOTCINT__)
#include "AliGenPythia.h"
#include "AliDecayerPythia8.h"
#include "AliGenParam.h"
#include "AliGenerator.h"
#include "AliGenCocktail.h"
#include "AliGenHijing.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliStack.h"
#include "AliGenPythiaPlus.h"
#include "TPythia8.h"
#endif

// settings
Bool_t generateLambda=kFALSE;
Bool_t generateOmega=kFALSE;
Bool_t generateKaon=kFALSE;
Bool_t generateOmegaccc=kTRUE;
Bool_t addPYTHIAevent=kTRUE;
Bool_t addHIJINGevent=kFALSE;

// generator configuration
Float_t bminConfig=0.;
Float_t bmaxConfig=5.;
Float_t energyConfig=5441.;
TString systemConfig="p-p";
Int_t seedConfig=0;

AliGenCocktail *Generator();
AliGenCocktail *GeneratorCocktail(TString name);
AliGenerator   *GeneratorHIJING();
AliGenerator*  GeneratorPythia8();
void ForceMultiHFLS4Upgrade(AliTPythia8 *pythia8);

AliRun *gAlice = new AliRun("test","");
/*******************************************************/
void testGeneratorXicc(Int_t nev = 1){

    //`export PYTHIA8DATA="/Users/jinjooseo/alice/AliRoot/PYTHIA8/pythia8210/share/Pythia8/xmldoc"
    //export PYTHIA8="/Users/jinjooseo/alice/AliRoot/PYTHIA8/pythia8210"

    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
   gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
// gSystem->Load("libgeant321");

/*    const char *p8dataenv = gSystem->Getenv("Users/jinjoo/alice/AliRoot/PYTHIA8/pythia8210/share/Pythia8/xmldoc");
		 if (!p8dataenv) {
        const char *p8env = gSystem->Getenv("/Users/jinjoo/alice/AliRoot/PYTHIA8/pythia8210");
        if (!p8env) {
            Error("pythia8.C",
                  "Environment variable PYTHIA8 must contain path to pythia directory!");
            return;
        }
        TString p8d = p8env;
        p8d += "/xmldoc";
        gSystem->Setenv("PYTHIA8DATA", p8d);
    }
    const char* path = gSystem->ExpandPathName("Users/jinjooseo/alice/AliRoot/PYTHIA8/pythia8210/share/Pythia8/xmldoc");
    if (gSystem->AccessPathName(path)) {
        Error("pythia8.C",
              "Environment variable PYTHIA8DATA must contain path to $PYTHIA8/xmldoc directory !");
        return;
    }*/
    // Load libraries

  //  gSystem->Load("libEG");
//    gSystem->Load("libEGPythia8");

//		 gSystem->Load("libeg.so");
//		 gSystem->Load("libEGPythia8.so");


    AliRunLoader* runLoader=0x0;
    cout<<"testGenerator.C: Creating Run Loader ..."<<endl;
    runLoader = AliRunLoader::Open("galice.root",
                            AliConfig::GetDefaultEventFolderName(),
                            "recreate");
    if (!runLoader) {
        gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
        return;
    }
    runLoader->SetCompressionLevel(2);
    runLoader->SetNumberOfEventsPerFile(nev);
    runLoader->CdGAFile();
    gAlice->SetRunLoader(runLoader);

    runLoader->LoadKinematics("RECREATE");
    runLoader->MakeTree("E");

    // Create stack and header
    runLoader->MakeStack();
    AliStack*  stack      = runLoader->Stack();
    AliHeader* header     = runLoader->GetHeader();

    // Intialize generator
    AliGenCocktail *generator = Generator();

    generator->Init();
    generator->SetStack(stack);


    for (Int_t iev=0; iev<nev; iev++)
    {
        // Initialize event
        header->Reset(0,iev);
        runLoader->SetEventNumber(iev);
        stack->Reset();
        runLoader->MakeTree("K");

        // Generate event
        generator->Generate();
        //  Analysis
        Int_t npart = stack->GetNprimary();
      //  printf("Analyse %d Particles\n", npart);
        for (Int_t part=0; part<npart; part++) {
            TParticle *MPart = stack->Particle(part);
            Int_t mpart  = MPart->GetPdgCode();
            //printf("Particle %d\n", mpart);
        }

        // Finish event
        header->SetNprimary(stack->GetNprimary());
        header->SetNtrack(stack->GetNtrack());
        stack->FinishEvent();
        header->SetStack(stack);
        runLoader->TreeE()->Fill();
        runLoader->WriteKinematics("OVERWRITE");
		cout << iev << "/"<<nev << endl;
    }
    //runLoader->WriteKinematics("OVERWRITE");
    // Finalize
    generator->FinishRun();
    // Write file
    runLoader->WriteHeader("OVERWRITE");
    generator->Write();
    runLoader->Write();


}
/*** HIJING ****************************************************/

AliGenerator *GeneratorHIJING()
{
    //
    // Libraries
    gSystem->Load("libHIJING");
    gSystem->Load("libTHijing");

    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy
    gener->SetEnergyCMS(energyConfig);
    gener->SetImpactParameterRange(bminConfig, bmaxConfig);
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectile-target
    if (systemConfig.EqualTo("Pb-Pb")) {
        gener->SetProjectile("A", 208, 82);
        gener->SetTarget    ("A", 208, 82);
        gener->SetSpectators(0);
        gener->SetDataDrivenSpectators();
    }
    else if (systemConfig.EqualTo("Xe-Xe")) {
        gener->SetProjectile("A", 129, 54);
        gener->SetTarget    ("A", 129, 54);
        gener->SetSpectators(0);
        gener->SetDataDrivenSpectators();
    }
    else if (systemConfig.EqualTo("p-Pb")) {
        gener->SetProjectile("P", 1, 1);
        gener->SetTarget    ("A", 208, 82);
        gener->SetBoostLHC(1);
        gener->SetSpectators(1);
    }
    else if (systemConfig.EqualTo("Pb-p")) {
        gener->SetProjectile("A", 208, 82);
        gener->SetTarget    ("P", 1, 1);
        gener->SetBoostLHC(1);
        gener->SetSpectators(1);
    }
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(0);
    // enable shadowing
    gener->SetShadowing(1);
    // neutral pion and heavy particle decays switched off
    //  gener->SetDecaysOff(1);
    // kinematic selection
    gener->SetSelectAll(0);
    gener->SetPtHardMin (2.3);
    return gener;
}
//------------------
/*** COCKTAIL ****************************************************/

AliGenCocktail *GeneratorCocktail(TString name)
{

    // configure projectile/target
    TString projN, targN;
    Int_t projA, projZ, targA, targZ;
    // pp
    if (systemConfig.EqualTo("p-p")) {
        projN = "p"; projA = 1; projZ = 1;
        targN = "p"; targA = 1; targZ = 1;
    }
    // PbPb
    else if (systemConfig.EqualTo("Pb-Pb")) {
        projN = "A"; projA = 208; projZ = 82;
        targN = "A"; targA = 208; targZ = 82;
    }
    // XeXe
    else if (systemConfig.EqualTo("Xe-Xe")) {
        projN = "A"; projA = 129; projZ = 54;
        targN = "A"; targA = 129; targZ = 54;
    }
    else if (systemConfig.EqualTo("Pb-p")) {
        projN = "A"; projA = 208; projZ = 82;
        targN = "p"; targA = 1; targZ = 1;
    }
    else if (systemConfig.EqualTo("p-Pb")) {
        projN = "p"; projA = 1; projZ = 1;
        targN = "A"; targA = 208; targZ = 82;
    }
    // not implemented
    else {
        printf("Cocktail not implemented for %s system\n", systemConfig.Data());
        abort();
    }

    //
    AliGenCocktail *ctl = new AliGenCocktail();
    ctl->SetName(name.Data());
    ctl->SetProjectile(projN, projA, projZ);
    ctl->SetTarget(targN, targA, targZ);
    ctl->SetEnergyCMS(energyConfig);
    ctl->SetSeed(seedConfig);
    return ctl;
}
/*** Generator ****************************************************/
AliGenCocktail *Generator(){
    AliGenCocktail *gen = GeneratorCocktail("HIJING_MULTICHARM");

    AliGenerator *hij = GeneratorHIJING();
    if(addHIJINGevent) gen->AddGenerator(hij,"Hijing",1.);

    AliGenerator *py = GeneratorPythia8();
    if(addPYTHIAevent) gen->AddGenerator(py,"PYTHIA",1.);

    AliGenParam *omegaccc = new AliGenParam(1000, new AliGenITSULib(), 4422, "DIST");//4332, "DIST");//4444, "DIST");
    //AliGenParam *omegaccc = new AliGenParam(1000, 4422);//4332, "DIST");//4444, "DIST");
		AliTPythia8* pythia8 = AliTPythia8::Instance();
   /* pythia8->ReadString("4444:onMode = off");
    pythia8->ReadString("4444:onIfAny = 4432");

    pythia8->ReadString("4432:onMode = off");
    pythia8->ReadString("4432:onIfAny = 4332");

    pythia8->ReadString("4332:onMode = off");
    pythia8->ReadString("4332:onIfAny = 3334");

    pythia8->ReadString("3334:onMode = off");
    pythia8->ReadString("3334:onIfAny = 3122");*/
    //ForceMultiHFLS4Upgrade(pythia8);

    AliDecayerPythia8 *d = new AliDecayerPythia8();
    d->Init();
    d->SetForceDecay(kHFLS4Upgrade);
    //d->HeavyFlavourOff();
    //omegaccc->SetDecayer(d);
    omegaccc->SetPtRange(0, 20.);
    omegaccc->SetYRange(-1., +1.);
    omegaccc->SetPhiRange(0., 360.);
    //omegaccc->SetPreserveFullDecayChain(1);
    omegaccc->SetForceDecay(kHFLS4Upgrade);
    omegaccc->SetSelectAll(kTRUE);
    omegaccc->SetPreserveFullDecayChain(1);
    omegaccc->SetDecayer(d);
    if(generateOmegaccc) {
        omegaccc->Init();
        gen->AddGenerator(omegaccc,"Xicc",1.);
    }

   /* AliGenParam *omega= new AliGenParam(1000, new AliGenSTRANGElib(), AliGenSTRANGElib::kProton, "DIST");
    AliGenParam *lambda= new AliGenParam(1000, new AliGenSTRANGElib(), AliGenSTRANGElib::kPion, "DIST");
    AliGenParam *kaon= new AliGenParam(1000, new AliGenSTRANGElib(), AliGenSTRANGElib::kKaon, "DIST");
   */ //AliTPythia8* pythia8 = AliTPythia8::Instance();
	/*	 pythia8->ReadString("3334:oneChannel = 1 1 0 3122 211");
    pythia8->ReadString("3334:onMode = off");
    pythia8->ReadString("3334:onIfAll = 3122 211");
    //pythia8->ReadString("3334:onIfAny = 3322");
    //
 pythia8->ReadString("3122:oneChannel = 1 1 0 2212 211");
    pythia8->ReadString("3122:onMode = off");
    pythia8->ReadString("3122:onIfAll = 2212 211");

    pythia8->ReadString("310:onMode = off");
    pythia8->ReadString("310:onIfMatch = 211");
*/
    //pythia8->ReadString("3322:onIfAny = off");
    //pythia8->ReadString("3322:onIfAny = 3122");

    //AliDecayerPythia8 *d = new AliDecayerPythia8();
 /*   omega->SetDecayer(d);
    omega->SetPtRange(0, 50.);
    omega->SetYRange(-1., +1.);
    omega->SetPhiRange(0., 360.);
    omega->SetPreserveFullDecayChain(1);

    lambda->SetDecayer(d);
    lambda->SetPtRange(0, 50.);
    lambda->SetYRange(-1., +1.);
    lambda->SetPhiRange(0., 360.);
    lambda->SetPreserveFullDecayChain(1);

    kaon->SetDecayer(d);
    kaon->SetPtRange(0, 50.);
    kaon->SetYRange(-1., +1.);
    kaon->SetPhiRange(0., 360.);
    kaon->SetPreserveFullDecayChain(1);
*/

   /* if(generateOmega) {
        omega->Init();
        gen->AddGenerator(omega,"Proton",1.);
    }

    if(generateLambda) {
        lambda->Init();
        gen->AddGenerator(lambda,"Pion",1.);
    }
    if(generateKaon) {
        kaon->Init();
        gen->AddGenerator(kaon,"Kaon",1.);
    }

*/
    return gen;
    //return omegaccc;

}
/*** PYTHIA ****************************************************/
AliGenerator*  GeneratorPythia8()
{
    AliGenPythiaPlus *pythia = new AliGenPythiaPlus(AliPythia8::Instance());
    pythia->SetMomentumRange(0, 999999.);
    pythia->SetThetaRange(0., 180.);
    pythia->SetYRange(-1.,1.);
    pythia->SetPtRange(0,100.);
    pythia->SetProcess(kPyMbDefault); // pythia->SetProcess(kPyMb);
    pythia->SetEnergyCMS(energyConfig);
    pythia->SetEventListRange(-1, 2);
    //pythia->SetTriggerParticle(3122, 1.);

    return pythia;
}

void ForceMultiHFLS4Upgrade(AliTPythia8 *pythia8){
  //Xi_cc++ -> Xi_c+ pi+
  pythia8->ReadString("4422:oneChannel = 1 1 0 4232 211");
  pythia8->ReadString("4422:onMode = off");
  pythia8->ReadString("4422:onIfAll = 4232 211");

  //Xi_c+ -> p K- pi+
  //pythia8->ReadString("4232:oneChannel = 1 1 0 2212 321 211");
  pythia8->ReadString("4232:onMode = off");
  pythia8->ReadString("4232:onIfAll = 2212 321 211");

  //Omega -> Lambda K
 // pythia8->ReadString("3334:onMode = off");
 // pythia8->ReadString("3334:onIfAll = 3122 321");

  //Lambda -> p pi-
  //pythia8->ReadString("3122:onMode = off");
  //pythia8->ReadString("3122:onIfAll = 2212 211");

}
