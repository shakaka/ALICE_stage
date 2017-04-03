//To run on the Virtual Analysis Facility with PoD:
// ssh username@alivaf-003.cern.ch -Y
//  vaf-enter
//  vafctl stop
//  vafctl start
//  vafreq <num_of_workers>
//  vafcount
//Wait for a certain number of workers, then:
//  root RunPoD.C
//  exit


void RunPoD(
	    Int_t runNumber = 245346, //246434
	    TString dataset = "Find;"
                    "BasePath=/alice/data/2015/LHC15o/%09d/muon_calo_pass1/AOD175/*/;"
                    "FileName=AliAOD.Muons.root;"
                    "Tree=/aodTree;"
                    "Mode=remote;",  // <-- much faster dataset creation
  Bool_t usePhysicsSelection = kTRUE,
  Int_t numEvents = 99999999,
  Int_t firstEvent = 0
) {

  // Not needed on the VAF
  //gEnv->SetValue("XSec.GSI.DelegProxy","2");

  TList *list = new TList();
  //list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice"));  // normally not needed, only used in special cases
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker

  // Not needed on the VAF
  //TProof::Mgr("alice-caf.cern.ch")->SetROOTVersion("VO_ALICE@ROOT::v5-34-08");

  // Note the difference between CAF and VAF
  //TProof::Open("alice-caf.cern.ch");
  TProof::Open("pod://");

  // Assemble dataset from format
  TString datasetWithRun;
  datasetWithRun.Form(dataset.Data(), runNumber);

  // Check the dataset before running the analysis!
  gProof->ShowDataSet(datasetWithRun.Data());
  //return;  // <-- uncomment this to test search before running the analysis!

  // Not needed on the VAF
  //gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-81-AN", list);

  // A single AliRoot package for *all* AliRoot versions: new on VAF
  TFile::Cp("http://alibrary.web.cern.ch/alibrary/vaf/AliceVaf.par", "AliceVaf.par");
  gProof->UploadPackage("AliceVaf.par");
  gProof->EnablePackage("AliceVaf.par", list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  gProof->Load("AliAnalysisJpsi.cxx+g");  // DON'T use double '+' when running multiple times: it uselessly recompiles everything!
  gROOT->LoadMacro("AddTaskJpsi.C");

  TString outputFileName;
  outputFileName.Form("AnalysisResults_run%09d.root", runNumber);
  AliAnalysisJpsi *simplePtTask = AddTaskJpsi(usePhysicsSelection, outputFileName);

  /*if (usePhysicsSelection) { //not needed since simplePtTask is using AliMuonEventCuts
    simplePtTask->SelectCollisionCandidates(AliVEvent::kAny);
    }*/

  if (!mgr->InitAnalysis()) return;

  mgr->StartAnalysis("proof", datasetWithRun, numEvents, firstEvent);

}
