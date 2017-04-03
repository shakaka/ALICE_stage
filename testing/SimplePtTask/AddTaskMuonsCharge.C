AliAnalysisMuonsCharge *AddTaskMuonsCharge(Bool_t usePhysicsSelection, TString outputFileName = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuonsCharge", "No analysis manager to connect to.");
    return NULL;
  }


  // output file name not defined
  if (outputFileName.IsNull()) outputFileName = mgr->GetCommonFileName();

 //Create and configure task
  AliAnalysisMuonsCharge *task = new AliAnalysisMuonsCharge("AliAnalysisMuonsCharge");
  if (!task) {
    Error("AliAnalysisMuonsCharge","Simple Pt task cannot be created! ");
    return NULL;
  }
  task->GetEventCuts()->SetFilterMask( AliMuonEventCuts::kSelectedTrig | AliMuonEventCuts::kPhysicsSelected );
  task->GetEventCuts()->SetTrigClassPatterns("CMUL7-B-NOPF-MUFAST");

  // Add task to analysis manager
  mgr->AddTask(task);

  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // Create and connect output container
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer( "listOfHisto",   TObjArray::Class(),  AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputEventStat = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutputEventStat);

  return task;
}
