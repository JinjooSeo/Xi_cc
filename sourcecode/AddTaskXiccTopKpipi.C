AliAnalysisTaskSEXiccTopKpipi *AddTaskXiccTopKpipi()
{
    AliAnalysisTaskSEXiccTopKpipi *task = new AliAnalysisTaskSEXiccTopKpipi("MCanalysis");
    bool isignalProd = true;
    task->SetIsSignalProd(isignalProd);

    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
        return NULL;
    }


    if(!mgr->GetMCtruthEventHandler()){
        Printf("AliAnalysisHFCorrOnFlySim; This task requires an input MC event handler");
        return NULL;
    }

    mgr->AddTask(task);

    // Create containers for input/output

    TString inname = "cinput";
    TString histoname = "coutputEntries";
    TString listnamegen = "coutputListGen";
    TString listnamereco = "coutputListReco";
    TString treenamereco = "coutputTreeRecoTracks";
    TString treenamecutvarc = "coutputTreeCutVarC";
    TString treenamecutvarcc = "coutputTreeCutVarCC";
    TString treenamecutvarccc = "coutputTreeCutVarCCC";
    TString treenamereconstruction = "coutputTreeReconstruction";
    TString treenamegen = "coutputTreeGenOmegaccc";
    TString treenamecascV0 = "coutputTreeCascV0";
    TString treenameV0 = "coutputTreeV0";
    TString histocounters = "coutputCounters";


    AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":MCanalysis";

    AliAnalysisDataContainer *coutputEntries = mgr->CreateContainer(histoname,TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCounters = mgr->CreateContainer(histocounters,TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputListGen = mgr->CreateContainer(listnamegen,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputListReco = mgr->CreateContainer(listnamereco,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputTreeCutVarGen = mgr->CreateContainer("coutputTreeCutVarGen",TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeCutVarGen->SetSpecialOutput();
    AliAnalysisDataContainer *coutputTreeCutVarReco = mgr->CreateContainer("coutputTreeCutVarReco",TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeCutVarReco->SetSpecialOutput();
    AliAnalysisDataContainer *coutputTreeCutVarRecoTrack = mgr->CreateContainer("coutputTreeCutVarRecoTrack",TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeCutVarRecoTrack->SetSpecialOutput();



    mgr->ConnectInput(task,0,cinput);
    mgr->ConnectOutput(task,1,coutputEntries);
    mgr->ConnectOutput(task,2,coutputListGen);
    mgr->ConnectOutput(task,3,coutputTreeCutVarGen);
    mgr->ConnectOutput(task,4,coutputTreeCutVarReco);
    mgr->ConnectOutput(task,5,coutputTreeCutVarRecoTrack);

    return task;



}
