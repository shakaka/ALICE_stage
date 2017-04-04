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

#include "AliDoubleJpsi.h"

using std::cout;
using std::endl;

ClassImp(AliDoubleJpsi);

//__________________________________________________________________________
AliDoubleJpsi::AliDoubleJpsi() :
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
AliDoubleJpsi::AliDoubleJpsi(const char *name) :
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
AliDoubleJpsi::~AliDoubleJpsi()
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
void AliDoubleJpsi::NotifyRun()
{
 /// Notify run
  fMuonTrackCuts->SetRun(fInputHandler);

}


//___________________________________________________________________________
void AliDoubleJpsi::UserCreateOutputObjects(){

  // Output objects creation
  fOutput = new TObjArray(2000);
  fOutput->SetOwner();

  // Dimuon Charge distribution
  TH2F *hDbJpsi= new TH2F("hDbJpsi", "Double Jpsi invariant mass Distribution", 200, 2, 5, 200, 2, 5);
  hDbJpsi->Sumw2();
  fOutput->AddAtAndExpand( hDbJpsi, kDbJpsi );

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
void AliDoubleJpsi::UserExec(Option_t *)
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
  // AliMultSelection *multSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
  // Double_t centralityFromV0 = multSelection->GetMultiplicityPercentile("V0M", false);
  // if(centralityFromV0 > 90) return;


  //Loop to match up Dimouns
  Int_t nTracks = 0;
  nTracks = fAODEvent->GetNumberOfTracks();

  //Loop to match up 2 dimuons
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



      for (Int_t lTrack = jTrack+1; lTrack < nTracks; lTrack++) {
        AliAODTrack *track3 = (AliAODTrack*) fAODEvent->GetTrack(lTrack);
        if (!track3) {
          AliError(Form("ERROR: Could not retrieve track3 %d", lTrack));
          continue;
        }
        if ( ! fMuonTrackCuts->IsSelected(track3) ) continue;



        for (Int_t oTrack = lTrack+1; oTrack < nTracks; oTrack++) {
          AliAODTrack *track4 = (AliAODTrack*) fAODEvent->GetTrack(oTrack);
          if (!track4) {
            AliError(Form("ERROR: Could not retrieve track4 %d", oTrack));
            continue;
          }
          if ( ! fMuonTrackCuts->IsSelected(track4) ) continue;

          Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2();
          TLorentzVector lvMuon1, lvMuon2, lvDimuon, lvMuon3, lvMuon4, lvDimuon2;

          Float_t energy = TMath::Sqrt(track1->P()*track1->P() + muonMass2);
          Float_t energy2 = TMath::Sqrt(track2->P()*track2->P() + muonMass2);
          Float_t energy3 = TMath::Sqrt(track3->P()*track3->P() + muonMass2);
          Float_t energy4 = TMath::Sqrt(track4->P()*track4->P() + muonMass2);

          lvMuon1.SetPxPyPzE(track1->Px(),track1->Py(),track1->Pz(),energy);
          lvMuon2.SetPxPyPzE(track2->Px(),track2->Py(),track2->Pz(),energy2);
          lvMuon3.SetPxPyPzE(track3->Px(),track3->Py(),track3->Pz(),energy3);
          lvMuon4.SetPxPyPzE(track4->Px(),track4->Py(),track4->Pz(),energy4);
          //Choose 12, 34
          if(track1->Charge()*track2->Charge()==-1){
            if (track3->Charge()*track4->Charge()==-1){
              //Building Lorentz vetors
              lvDimuon = lvMuon1 + lvMuon2;
              lvDimuon2= lvMuon3 + lvMuon4;

              //Getting invariant mass
              Double_t maDiMu1 = lvDimuon.M();
              Double_t maDiMu2 = lvDimuon2.M();

              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu1, maDiMu2);
              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu2, maDiMu1);
            }
          }
          //Choose 13, 24
          if(track1->Charge()*track3->Charge()==-1){
            if (track2->Charge()*track4->Charge()==-1){
              //Building Lorentz vetors
              lvDimuon = lvMuon1 + lvMuon3;
              lvDimuon2= lvMuon2 + lvMuon4;

              //Getting invariant mass
              Double_t maDiMu1 = lvDimuon.M();
              Double_t maDiMu2 = lvDimuon2.M();

              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu1, maDiMu2);
              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu2, maDiMu1);
            }
          }
          //Choose 14,23
          if(track1->Charge()*track4->Charge()==-1){
            if (track2->Charge()*track3->Charge()==-1){
              //Building Lorentz vetors
              lvDimuon = lvMuon1 + lvMuon4;
              lvDimuon2= lvMuon2 + lvMuon3;

              //Getting invariant mass
              Double_t maDiMu1 = lvDimuon.M();
              Double_t maDiMu2 = lvDimuon2.M();

              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu1, maDiMu2);
              ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->Fill(maDiMu2, maDiMu1);
            }
          }


        }
      }
    }
  }
  TH1D *hBin1 = ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->ProjectionX("till 3.06 GeV", 66, 70);
  TH1D *hBin2 = ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->ProjectionX("till 3.12 GeV", 74, 78);
  TH1D *hBin3 = ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->ProjectionX("till 3.18 GeV", 82, 86);
  TH1D *hBin4 = ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->ProjectionX("till 3.24 GeV", 90, 94);
  TH1D *hBin5 = ( (TH2F*)fOutput->UncheckedAt(kDbJpsi) )->ProjectionX("till 3.3 GeV", 98, 102);


  fOutput->AddAtAndExpand( hBin1, kBin1 );
  fOutput->AddAtAndExpand( hBin2, kBin2 );
  fOutput->AddAtAndExpand( hBin3, kBin3 );
  fOutput->AddAtAndExpand( hBin4, kBin4 );
  fOutput->AddAtAndExpand( hBin5, kBin5 );

  // Required both here and in UserCreateOutputObjects()
  PostData(1, fOutput);
  PostData(2, fEventCounters);
}

//___________________________________________________________________________
void AliDoubleJpsi::Terminate(Option_t *)
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
