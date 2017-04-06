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
#include "AliMultSelection.h"



// ALIPHYSICS/PWG includes
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisJpsi.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisJpsi);

//__________________________________________________________________________
AliAnalysisJpsi::AliAnalysisJpsi() :
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
AliAnalysisJpsi::AliAnalysisJpsi(const char *name) :
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
/*AliAnalysisMuonsCharge& AliAnalysisMuonsCharge::operator=(const AliAnalysisMuonsCharge& c)
{
  // Assignment operator
  if (this != &c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
  }*/

//___________________________________________________________________________
/*AliAnalysisMuonsCharge::AliAnalysisMuonsCharge(const AliAnalysisMuonsCharge& c) :
  AliAnalysisTaskSE(c),
  fAODEvent(c.fAODEvent),
  fOutput(c.fOutput)
 {
  // Copy ctor
 }
*/
//___________________________________________________________________________
AliAnalysisJpsi::~AliAnalysisJpsi()
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
void AliAnalysisJpsi::NotifyRun()
{
 /// Notify run
  fMuonTrackCuts->SetRun(fInputHandler);

}


//___________________________________________________________________________
void AliAnalysisJpsi::UserCreateOutputObjects(){

  // Output objects creation
  fOutput = new TObjArray(2000);
  fOutput->SetOwner();

  // Dimuon Charge distribution
  TH1F *hDiMuCh= new TH1F("hDiMuCh", "Dimuon Charge Distribution", 200, -2, 2);
  hDiMuCh->Sumw2();
  fOutput->AddAtAndExpand( hDiMuCh, kDiMuCh );

  // Dimuon pt Distribution
  TH1F *hDiMuPt = new TH1F("hDiMuPt", "Dimuon pt Distribution", 200, 0, 10);
  hDiMuPt->Sumw2();
  fOutput->AddAtAndExpand( hDiMuPt, kDiMuPt );

  // Dimuon Y distribution
  TH1F *hDiMuY = new TH1F("hDiMuY", "Dimuon Y Distribution", 200, 0, 5);
  hDiMuY->Sumw2();
  fOutput->AddAtAndExpand( hDiMuY, kDiMuY );

  // Dimuon Phi distribution
  TH1F *hDiMuPhi = new TH1F("hDiMuPhi", "Dimuon Phi Distribution", 200, 0, 7);
  hDiMuPhi->Sumw2();
  fOutput->AddAtAndExpand( hDiMuPhi, kDiMuPhi );

  // Dimuon M distribution
  TH1F *hDiMuM = new TH1F("hDiMuM", "Dimuon M Distribution", 200, 2, 5);
  hDiMuM->Sumw2();
  fOutput->AddAtAndExpand( hDiMuM, kDiMuM );

  // Event centrality distribution
  TH1F *hEvCen= new TH1F("hEvCen", "Event centrality distribution", 200, 0, 100);
  hEvCen->Sumw2();
  fOutput->AddAtAndExpand( hEvCen, kEvCen );

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
void AliAnalysisJpsi::UserExec(Option_t *)
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
    fillName = Form("trigger:%s/run:%d/selected:%s",triggerName.Data(),fCurrentRunNumber,selected.Data());
    fEventCounters->Count(fillName.Data());
  }

  //keep only selected events
  if ( !keepEvent ) return;
  AliMultSelection *multSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
  Double_t centralityFromV0 = multSelection->GetMultiplicityPercentile("V0M", false);
  if(centralityFromV0 > 90) return;

  ( (TH1F*)fOutput->UncheckedAt(kEvCen) )->Fill( centralityFromV0 );

  //Loop to match up Dimouns
  Int_t nTracks = 0;
  nTracks = fAODEvent->GetNumberOfTracks();

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


      //Charges selected
      Int_t chDiMu = track1->Charge()*track2->Charge();

      if(!(chDiMu==-1)) continue;
      ( (TH1F*)fOutput->UncheckedAt(kDiMuCh) )->Fill(chDiMu);

      //Building Lorentz vetors
      Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2();
      TLorentzVector lvMuon1, lvMuon2, lvDimuon;

      Float_t energy = TMath::Sqrt(track1->P()*track1->P() + muonMass2);
      Float_t energy2 = TMath::Sqrt(track2->P()*track2->P() + muonMass2);

      lvMuon1.SetPxPyPzE(track1->Px(),track1->Py(),track1->Pz(),energy);
      lvMuon2.SetPxPyPzE(track2->Px(),track2->Py(),track2->Pz(),energy2);
      lvDimuon = lvMuon1 + lvMuon2;

      //Getting Pt of dimuon
      Double_t ptDiMu = lvDimuon.Pt();

      //Getting Pl of dimuon
      Double_t plDiMu = TMath::Sqrt((lvDimuon.P()*lvDimuon.P())-(lvDimuon.Pt()*lvDimuon.Pt()));

      //Getting M of dimuon
      Double_t maDiMu = lvDimuon.M();

      //Getting Phi of dimuon
      Double_t phiDiMu = lvDimuon.Phi();

      //Getting rapidity of dimuon
      Double_t yDiMu = TMath::Log((lvDimuon.E()+plDiMu)/(lvDimuon.E()-plDiMu))/2;

      //Getting the distance between Dimuons
      Double_t zvMu1 = track1->Zv();
      Double_t zvMu2 = track2->Zv();
      Double_t disDiMu = zvMu1-zvMu2;

      ( (TH1F*)fOutput->UncheckedAt(kDiMuPt) )->Fill(ptDiMu);
      ( (TH1F*)fOutput->UncheckedAt(kDiMuM) )->Fill(maDiMu);
      ( (TH1F*)fOutput->UncheckedAt(kDiMuY) )->Fill(yDiMu);
      ( (TH1F*)fOutput->UncheckedAt(kDiMuPhi) )->Fill(phiDiMu);


    }

  }

  // Required both here and in UserCreateOutputObjects()
  PostData(1, fOutput);
  PostData(2, fEventCounters);
}

//___________________________________________________________________________
void AliAnalysisJpsi::Terminate(Option_t *)
{
  /*fOutput = dynamic_cast<TObjArray*> (GetOutputData(1));
  if (fOutput) {
    TH1F *hDiMuCh = dynamic_cast<TH1F *>( fOutput->FindObject("hDiMuCh") );
    if (hDiMuCh) {
      new TCanvas();
      hDiMuCh->Draw();
    }
  }
  */


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
