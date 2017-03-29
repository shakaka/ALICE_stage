// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "AliAnalysisMuonUtility.h"

/// ALIROOT/STEER includes
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

/// ALIROOT/ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"



// ALIPHYSICS/PWG includes
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisTaskSimplePt.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSimplePt);

//__________________________________________________________________________
AliAnalysisTaskSimplePt::AliAnalysisTaskSimplePt() :
  AliAnalysisTaskSE(),
  fAODEvent(0x0),
  fMuonEventCuts(0x0),
  fMuonTrackCuts(0x0),
  fOutput(0x0),
  fEventCounters(0x0)
{
  // Default ctor.
}

//__________________________________________________________________________
AliAnalysisTaskSimplePt::AliAnalysisTaskSimplePt(const char *name) :
  AliAnalysisTaskSE(name),
  fAODEvent(0x0),
  fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
  fMuonTrackCuts(new AliMuonTrackCuts("stdMuonCuts","stdMuonCuts")),
  fOutput(0x0),
  fEventCounters(0x0)
{

  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1, TObjArray::Class());
  // Output slot #2 writes event counters
  DefineOutput(2, AliCounterCollection::Class());

  //  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca );
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuMatchLpt );
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

}

//___________________________________________________________________________
/*AliAnalysisTaskSimplePt& AliAnalysisTaskSimplePt::operator=(const AliAnalysisTaskSimplePt& c)
{
  // Assignment operator
  if (this != &c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
  }*/

//___________________________________________________________________________
/*AliAnalysisTaskSimplePt::AliAnalysisTaskSimplePt(const AliAnalysisTaskSimplePt& c) :
  AliAnalysisTaskSE(c),
  fAODEvent(c.fAODEvent),
  fOutput(c.fOutput)
 {
  // Copy ctor
 }
*/
//___________________________________________________________________________
AliAnalysisTaskSimplePt::~AliAnalysisTaskSimplePt()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) {
    if (fOutput) delete fOutput;
    if (fEventCounters) delete fEventCounters;
  }

  delete fMuonTrackCuts;
  delete fMuonEventCuts;

}

//___________________________________________________________________________
void AliAnalysisTaskSimplePt::NotifyRun()
{
 /// Notify run
  fMuonTrackCuts->SetRun(fInputHandler);

}


//___________________________________________________________________________
void AliAnalysisTaskSimplePt::UserCreateOutputObjects(){

  // Output objects creation
  fOutput = new TObjArray(2000);
  fOutput->SetOwner();

  /*
  //moun number distribution before selection
  TH1F *hNoMuon = new TH1F("hNoMuon", "Moun number Distribution", 300,0,300);
  hNoMuon->Sumw2();
  fOutput->AddAtAndExpand( hNoMuon, kNoMuon );

  //muon number distribution after selection
  TH1F *hNoMuonS = new TH1F("hNoMuonS", "Moun number Distribution Selected", 300, 0, 300);
  hNoMuonS->Sumw2();
  fOutput->AddAtAndExpand( hNoMuonS, kNoMuonS);

  //dimuon number distribution before selection
  TH1F *hNoDiMu = new TH1F("hNoDiMu", "Dimoun number Distribution", 300, 0, 300);
  hNoDiMu->Sumw2();
  fOutput->AddAtAndExpand( hNoDiMu, kNDiMu);

  //dimuon number distribution after selection
  TH1F *hNoDiMuS = new TH1F("hNoDiMuS", "Dimoun number Distribution Selected", 300, 0, 300);
  hNoDiMuS->Sumw2();
  fOutput->AddAtAndExpand( hNoDiMuS, kNDiMuS);*/




  // dimuon Pt distribution
  TH1F *hDimuonPt = new TH1F("hDimuonPt", "Dimuon pt Distribution", 200, 0, 10);
  hDimuonPt->Sumw2();
  fOutput->AddAtAndExpand( hDimuonPt, kDimuonPt );

  // single muon Pt distribution
  TH1F *hMuonPt = new TH1F("hMuonPt", "Single muon pt Distribution", 200, 0, 10);
  hMuonPt->Sumw2();
  fOutput->AddAtAndExpand( hMuonPt, kMuonPt );

  // event z vertex distribution
  TH1F *hEventZv = new TH1F("hEventZv", "Event z vertex Distribution", 200, -20, 20);
  hEventZv->Sumw2();
  fOutput->AddAtAndExpand( hEventZv, kEventZv );

  // single muon Eta distribution
  TH1F *hMuonEta = new TH1F("hMuonEta", "Single muon Eta Distribution", 200, -5, 0);
  hMuonEta->Sumw2();
  fOutput->AddAtAndExpand( hMuonEta, kMuonEta );

  // single muon Rab distribution
  TH1F *hMuonRab = new TH1F("hMuonRab", "Single muon Rab Distribution", 200, 0, 100);
  hMuonRab->Sumw2();
  fOutput->AddAtAndExpand( hMuonRab, kMuonRab );

  // single muon Phi distribution
  TH1F *hMuonPhi = new TH1F("hMuonPhi", "Single muon Phi Distribution", 200, 0, 7);
  hMuonPhi->Sumw2();
  fOutput->AddAtAndExpand( hMuonPhi, kMuonPhi );


  // dimuon Y distribution
  TH1F *hDiMuY = new TH1F("hDiMuY", "dimuon Y Distribution", 200, -5, 0);
  hDiMuY->Sumw2();
  fOutput->AddAtAndExpand( hDiMuY, kDiMuY );

  // dimuon Phi distribution
  TH1F *hDiMuPhi = new TH1F("hDiMuPhi", "dimuon Phi Distribution", 200, 0, 7);
  hDiMuPhi->Sumw2();
  fOutput->AddAtAndExpand( hDiMuPhi, kDiMuPhi );

  // dimuon M distribution
  TH1F *hDiMuM = new TH1F("hDiMuM", "dimuon M Distribution", 200, 0, 5);
  hDiMuM->Sumw2();
  fOutput->AddAtAndExpand( hDiMuM, kDiMuM );

  // dimuon zv distance distribution
  TH1F *hDiMuZv= new TH1F("hDiMuZv", "dimuon zv distance Distribution", 200, -0.5, 0.5);
  hDiMuZv->Sumw2();
  fOutput->AddAtAndExpand( hDiMuZv, kDiMuZv );

  // dimuon Charge distribution
  TH1F *hDiMuCh= new TH1F("hDiMuCh", "dimuon Charge Distribution", 200, -2, 2);
  hDiMuCh->Sumw2();
  fOutput->AddAtAndExpand( hDiMuCh, kDiMuCh );


  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("trigger",1000000);
  fEventCounters->AddRubric("run",1000000);
  fEventCounters->AddRubric("selected","yes/no");
  fEventCounters->Init();


  // Required both here and in UserExec()
  PostData(1, fOutput);
  PostData(2, fEventCounters);
}

//___________________________________________________________________________
void AliAnalysisTaskSimplePt::UserExec(Option_t *)
{
  // Execute analysis for current InputEvent

  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }

  //check event cuts
  Bool_t keepEvent = fMuonEventCuts->IsSelected(fInputHandler);

  //fill event counters
  //  TString firedTriggerClasses = fAODEvent->GetFiredTriggerClasses();
  const TObjArray *trig = fMuonEventCuts->GetSelectedTrigClassesInEvent(InputEvent());
  if (!trig) {
    AliError(Form("ERROR: Could not retrieve list of selected trigger in run %d", fCurrentRunNumber));
    return;
  }

  TString selected = keepEvent ? "yes" : "no";

  TString fillName;

  fillName = Form("trigger:any/run:%d/selected:%s",fCurrentRunNumber,selected.Data());
  fEventCounters->Count(fillName.Data());

  for ( Int_t iTrig = 0; iTrig < trig->GetEntries(); iTrig++ ) {

    TString triggerName = ( (TObjString*) trig->At(iTrig) )->GetString();

    //    cout<<"runNr = "<<fCurrentRunNumber<<" "<<iTrig<<" "<<triggerName<<endl;

    fillName = Form("trigger:%s/run:%d/selected:%s",triggerName.Data(),fCurrentRunNumber,selected.Data());
    fEventCounters->Count(fillName.Data());

  }

  //keep only selected events
  if ( !keepEvent ) return;

  AliAODVertex* vEvent = 0; //new AliAODVertex();
  vEvent = fAODEvent->GetPrimaryVertex();
  Double_t zvEvent = 0;
  if (vEvent) zvEvent = vEvent->GetZ();

  ( (TH1F*)fOutput->UncheckedAt(kEventZv) )->Fill( zvEvent );
  //
  // Loop over Muons
  //

  // muon tracks
  Int_t nTracks = 0;
  nTracks = fAODEvent->GetNumberOfTracks();

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    AliAODTrack *track1 = (AliAODTrack*) fAODEvent->GetTrack(iTrack);

    if (!track1) {
      AliError(Form("ERROR: Could not retrieve track %d", iTrack));
      continue;
    }
    /*    Int_t noMuon = track1->GetNumberOfMuonTracks();
    ( (TH1F*)fOutput->UncheckedAt(kNoMuon) )->Fill( noMuon );
    */


    //some basic checks
    //offline trigger (in LHC15o, Apt = 0.5 GeV/c, Lpt = 1 GeV/c, Hpt = 4 GeV/c)
    // matchTrig>=1: Apt
    // matchTrig>=2: Lpt
    // matchTrig>=3: Hpt
    Int_t matchTrig = track1->GetMatchTrigger();

    Bool_t isMuon = track1->IsMuonTrack();

    Float_t ptMuon = track1->Pt();

    Int_t charge = track1->Charge();


    Double_t etaMuon = track1->Eta();
    Double_t rabMuon = track1->GetRAtAbsorberEnd();
    Double_t phiMuon = track1->Phi();


    Int_t noMuonS = 0;
    //Select muon with AliMuonTrackCuts
    /*UInt_t selectionMask = fMuonTrackCuts->GetSelectionMask( track1 );
    UInt_t cutMask = AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuMatchLpt;
    Bool_t resMask = kTRUE;
    if ( !( ( selectionMask & cutMask ) == cutMask) ) resMask = kFALSE;

    if ( !resMask ) continue;
    */

    // simple way to select muon track
    if ( ! fMuonTrackCuts->IsSelected(track1) ) continue;


    ( (TH1F*)fOutput->UncheckedAt(kMuonPt) )->Fill( ptMuon );
    ( (TH1F*)fOutput->UncheckedAt(kMuonEta) )->Fill( etaMuon );
    ( (TH1F*)fOutput->UncheckedAt(kMuonRab) )->Fill( rabMuon );
    ( (TH1F*)fOutput->UncheckedAt(kMuonPhi) )->Fill( phiMuon );
    //    ( (TH1F*)fOutput->UncheckedAt(kMuonPDCA) )->Fill( pdcaMuon  );


  }
  //Loop to match up Dimouns

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    AliAODTrack *track1 = (AliAODTrack*) fAODEvent->GetTrack(iTrack);

    if (!track1) {
      AliError(Form("ERROR: Could not retrieve track1 %d", iTrack));
      continue;
    }
    if ( ! fMuonTrackCuts->IsSelected(track1) ) continue;

    for (Int_t jTrack = iTrack+1; jTrack < nTracks; jTrack++) {
      AliAODTrack *track2 = (AliAODTrack*) fAODEvent->GetTrack(jTrack);

      if (!track2) {
        AliError(Form("ERROR: Could not retrieve track2 %d", jTrack));
        continue;
      }
      if ( ! fMuonTrackCuts->IsSelected(track2) ) continue;

      AliAODDimuon *dimu = dynamic_cast<AliAODDimuon*>(track1, track2);

      fAODEvent->AddObject(dimu);

    }

  }





  // Loop over Dimuons
  for ( Int_t iDimuon = 0; iDimuon < fAODEvent->GetNumberOfDimuons(); iDimuon++) {

    AliAODDimuon *dimu = dynamic_cast<AliAODDimuon*>( fAODEvent->GetDimuon(iDimuon) );


    Bool_t isMuSelected[2] = {kFALSE,kFALSE};

    for (Int_t iMu = 0; iMu < 2; iMu++) {

      AliAODTrack *track = (AliAODTrack*) dimu->GetMu(iMu);
      if ( fMuonTrackCuts->IsSelected(track) ) isMuSelected[iMu] = kTRUE;

    }

    //    if ( !(isMuSelected[0] && isMuSelected[1]) ) continue;



    // ( (TH1F*)fOutput->UncheckedAt(kDimuonPt) )->Fill( dimu->Pt() );

    Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2();
    TLorentzVector lvMuon1, lvMuon2,lvDimuon;

    AliAODTrack *track2 = (AliAODTrack*) dimu->GetMu(0);
    AliAODTrack *track3 = (AliAODTrack*) dimu->GetMu(1);

    //if(track2->Charge()*track3->Charge()==1) continue;

    Short_t chMu = track2->Charge()*track3->Charge();

    Double_t pMu1 = track2->P();
    Double_t pMu2 = track3->P();

    Float_t energy = TMath::Sqrt(track2->P()*track2->P() + muonMass2);
    lvMuon1.SetPxPyPzE(track2->Px(),track2->Py(),track2->Pz(),energy);


    Float_t energy2 = TMath::Sqrt(track3->P()*track3->P() + muonMass2);
    lvMuon2.SetPxPyPzE(track3->Px(),track3->Py(),track3->Pz(),energy2);

    lvDimuon = lvMuon1 + lvMuon2;




    Double_t zvMu1 = track2->Zv();
    Double_t zvMu2 = track3->Zv();

    Double_t disDiMu = zvMu1-zvMu2;

    Double_t ptDiMu = dimu->Pt();
    Double_t mDiMu = lvDimuon.M();
    Double_t yDiMu = dimu->Y();
    Double_t phiDiMu = dimu->Phi();


    ( (TH1F*)fOutput->UncheckedAt(kDimuonPt) )->Fill(ptDiMu);
    ( (TH1F*)fOutput->UncheckedAt(kDiMuM) )->Fill(mDiMu);
    ( (TH1F*)fOutput->UncheckedAt(kDiMuY) )->Fill(yDiMu);
    ( (TH1F*)fOutput->UncheckedAt(kDiMuPhi) )->Fill(phiDiMu);
    ( (TH1F*)fOutput->UncheckedAt(kDiMuZv) )->Fill(disDiMu);
    ( (TH1F*)fOutput->UncheckedAt(kDiMuCh) )->Fill(chMu);
  }

  // Required both here and in UserCreateOutputObjects()
  PostData(1, fOutput);
  PostData(2, fEventCounters);
}

//___________________________________________________________________________
void AliAnalysisTaskSimplePt::Terminate(Option_t *)
{
  // Display the Pt histogram for dimuon

  fOutput = dynamic_cast<TObjArray*> (GetOutputData(1));
  if (fOutput) {
    /* TH1F *hDimuonPt = dynamic_cast<TH1F *>( fOutput->FindObject("hDimuonPt") );
    TH1F *hEventZv = dynamic_cast<TH1F *>( fOutput->FindObject("hEventZv") );
    TH1F *hMuonEta = dynamic_cast<TH1F *>( fOutput->FindObject("hMuonEta") );
    TH1F *hMuonRab = dynamic_cast<TH1F *>( fOutput->FindObject("hMuonRab") );
    TH1F *hMuonPt = dynamic_cast<TH1F *>( fOutput->FindObject("hMuonPt") );
    TH1F *hMuonPhi = dynamic_cast<TH1F *>( fOutput->FindObject("hMuonPhi") );
    TH1F *hDiMuY = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuY") );
    TH1F *hDiMuPhi = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuPhi") );
    TH1F *hDiMuM = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuM") );
    TH1F *hDiMuZv = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuZv") );
    TH1F *hDiMuCh = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuCh") );

    if (hDimuonPt) {
      hDimuonPt->Draw();
    }
    if (hEventZv) {
      new TCanvas();
      hEventZv->Draw();
    }
    if (hMuonEta) {
      new TCanvas();
      hMuonEta->Draw();
    }
    if (hMuonRab) {
      new TCanvas();
      hMuonRab->Draw();
    }
    if (hMuonPt) {
      new TCanvas();
      hMuonPt->Draw();
    }
    if (hMuonPhi) {
      new TCanvas();
      hMuonPhi->Draw();
    }
    if (hDiMuY) {
      new TCanvas();
      hDiMuY->Draw();
    }
    if (hDiMuPhi) {
      new TCanvas();
      hDiMuPhi->Draw();
    }
    if (hDiMuM) {
      new TCanvas();
      hDiMuM->Draw();
      }
    if (hDiMuZv) {
      new TCanvas();
      hDiMuZv->Draw();
    }
    if (hDiMuCh) {
      new TCanvas();
      hDiMuCh->Draw();
    }*/
  }

  //global statistics
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(2));
  if ( fEventCounters ) {
    if (!gROOT->IsBatch() ) {
      cout<<"Event statistics without any selection "<<endl;
      fEventCounters->Print("trigger/run");
      cout<<"Event statistics with event selection "<<endl;
      fEventCounters->Print("trigger/run","selected:yes");
      new TCanvas();
      fEventCounters->Draw("run","trigger","selected:yes");
    }
  }

}
