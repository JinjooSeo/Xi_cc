AliAnalysisTaskSEXiccTopKpipi *AddTaskOmegacccToOmega3Pi()
{
    //
    //

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
    AliAnalysisDataContainer *coutputTreeCutVarCCC = mgr->CreateContainer(treenamecutvarccc,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeCutVarCCC->SetSpecialOutput();
    AliAnalysisDataContainer *coutputTreeReconstruction= 0x0;
    AliAnalysisDataContainer *coutputTreeGen = 0x0;
    AliAnalysisDataContainer *coutputTreeCascV0 = 0x0;
    AliAnalysisDataContainer *coutputTreeV0 = 0x0;
    AliAnalysisDataContainer *coutputTreeCutVarC = mgr->CreateContainer(treenamecutvarc,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeCutVarC->SetSpecialOutput();
    AliAnalysisDataContainer *coutputTreeCutVarCC = 0x0;
    AliAnalysisDataContainer *coutputTreeReco = 0x0;
    if(isignalProd){

        coutputTreeCutVarCC = mgr->CreateContainer(treenamecutvarcc,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeCutVarCC->SetSpecialOutput();
        coutputTreeGen = mgr->CreateContainer(treenamegen,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGen->SetSpecialOutput();
        coutputTreeCascV0 = mgr->CreateContainer(treenamecascV0,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeCascV0->SetSpecialOutput();
        coutputTreeV0 = mgr->CreateContainer(treenameV0,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeV0->SetSpecialOutput();
        coutputTreeReconstruction = mgr->CreateContainer(treenamereconstruction,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeReconstruction->SetSpecialOutput();
        coutputTreeReco = mgr->CreateContainer(treenamereco,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeReco->SetSpecialOutput();
    }



    mgr->ConnectInput(task,0,cinput);
    mgr->ConnectOutput(task,1,coutputEntries);
    //mgr->ConnectOutput(task,2,coutputCounters);
    mgr->ConnectOutput(task,2,coutputListGen);
   // mgr->ConnectOutput(task,4,coutputListReco);
    mgr->ConnectOutput(task,3,coutputTreeCutVarCCC);
   /* mgr->ConnectOutput(task,6,coutputTreeCutVarC);

    if(isignalProd){
        mgr->ConnectOutput(task,7,coutputTreeCutVarCC);
        mgr->ConnectOutput(task,8,coutputTreeGen);
        mgr->ConnectOutput(task,9,coutputTreeCascV0);
        mgr->ConnectOutput(task,10,coutputTreeV0);
        mgr->ConnectOutput(task,11,coutputTreeReconstruction);
        mgr->ConnectOutput(task,12,coutputTreeReco);
    }
*/
    return task;



}
