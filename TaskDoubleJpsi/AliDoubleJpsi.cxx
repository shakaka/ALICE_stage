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
  TH2F *hMaDbJpsi= new TH2F("hMaDbJpsi", "Double Jpsi invariant mass Distribution", 200, 2, 5, 200, 2, 5);
  hMaDbJpsi->Sumw2();
  fOutput->AddAtAndExpand( hMaDbJpsi, kMaDbJpsi );

  TH2F *hPtDbJpsi= new TH2F("hPtDbJpsi", "Double Jpsi Pt Distribution", 200, 0, 15, 200, 0, 15);
  hPtDbJpsi->Sumw2();
  fOutput->AddAtAndExpand( hPtDbJpsi, kPtDbJpsi );

  TH2F *hYDbJpsi= new TH2F("hYDbJpsi", "Double Jpsi Y Distribution", 200, 2, 5, 200, 2, 5);
  hYDbJpsi->Sumw2();
  fOutput->AddAtAndExpand( hYDbJpsi, kYDbJpsi );

  TH2F *hPhiDbJpsi= new TH2F("hPhiDbJpsi", "Double Jpsi Phi Distribution", 200, -4, 4, 200, -4, 4);
  hPhiDbJpsi->Sumw2();
  fOutput->AddAtAndExpand( hPhiDbJpsi, kPhiDbJpsi );

  TH1F *hMaJpsi= new TH1F("hMaJpsi", "Inclusive Jpsi Ma Distribution", 200, 2, 5);
  hMaJpsi->Sumw2();
  fOutput->AddAtAndExpand( hMaJpsi, kMaJpsi );


  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("trigger",1000000);
  fEventCounters->AddRubric("run",1000000);
  fEventCounters->AddRubric("selected","yes/no");
  fEventCounters->AddRubric("numMu",1000000);
  fEventCounters->AddRubric("numDiMu",1000000);


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

  fillName = Form("trigger:any/run:%d/selected:%s/numMu:Sum/numDiMu:Non",fCurrentRunNumber,selected.Data());
  fEventCounters->Count(fillName.Data());

  for ( Int_t iTrig = 0; iTrig < trig->GetEntries(); iTrig++ ) {
    TString triggerName = ( (TObjString*) trig->At(iTrig) )->GetString();
    fillName = Form("trigger:%s/run:%d/selected:%s/numMu:Sum/numDiMu:Non",triggerName.Data(),fCurrentRunNumber,selected.Data());
    fEventCounters->Count(fillName.Data());
  }

  //keep only selected events
  if ( !keepEvent ) return;
  AliMultSelection *multSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
  Double_t centralityFromV0 = multSelection->GetMultiplicityPercentile("V0M", false);
  if(centralityFromV0 > 90) return;


  //Loop to match up Dimouns
  Int_t nTracks = 0;
  nTracks = fAODEvent->GetNumberOfTracks();


  Int_t counterDiMu = 0;
  Int_t counterMu = 0;
  //Loop to match up 2 dimuons
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliAODTrack *track1 = (AliAODTrack*) fAODEvent->GetTrack(iTrack);
    if (!track1) {
      AliError(Form("ERROR: Could not retrieve track1 %d", iTrack));
      continue;
    }
    if ( ! fMuonTrackCuts->IsSelected(track1) ) continue;
    counterMu++;

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
          TLorentzVector lvMuon1, lvMuon2, lvMuon3, lvMuon4, lvDimuon, lvDimuon2, lvDimuon3, lvDimuon4, lvTemp;

          Bool_t getDiMu = 0;

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

              lvDimuon = lvMuon1 + lvMuon2;
              lvDimuon2= lvMuon3 + lvMuon4;

              if(track1->Charge()*track3->Charge()==-1){
                if (track2->Charge()*track4->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon3;
                  lvDimuon4 = lvMuon2 + lvMuon4;
                }
              }

              if(track1->Charge()*track4->Charge()==-1){
                if (track2->Charge()*track3->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon4;
                  lvDimuon4 = lvMuon2 + lvMuon3;
                }
              }
            }
          }
          //Choose 13, 24
          if(track1->Charge()*track3->Charge()==-1){
            if (track2->Charge()*track4->Charge()==-1){

              lvDimuon = lvMuon1 + lvMuon3;
              lvDimuon2= lvMuon2 + lvMuon4;

              if(track1->Charge()*track2->Charge()==-1){
                if (track3->Charge()*track4->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon2;
                  lvDimuon4 = lvMuon3 + lvMuon4;
                }
              }

              if(track1->Charge()*track4->Charge()==-1){
                if (track2->Charge()*track3->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon4;
                  lvDimuon4 = lvMuon2 + lvMuon3;
                }
              }
            }
          }
          //Choose 14,23
          if(track1->Charge()*track4->Charge()==-1){
            if (track2->Charge()*track3->Charge()==-1){

              lvDimuon = lvMuon1 + lvMuon4;
              lvDimuon2= lvMuon2 + lvMuon3;

              if(track1->Charge()*track2->Charge()==-1){
                if (track3->Charge()*track4->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon2;
                  lvDimuon4 = lvMuon3 + lvMuon4;
                }
              }

              if(track1->Charge()*track3->Charge()==-1){
                if (track2->Charge()*track4->Charge()==-1){
                  getDiMu = kTRUE;
                  lvDimuon3 = lvMuon1 + lvMuon3;
                  lvDimuon4 = lvMuon2 + lvMuon4;
                }
              }
            }
          }

          //Getting data wanted


          if(getDiMu == kTRUE){

            counterDiMu+=2;

            Double_t ptDiMu1 = lvDimuon.Pt();
            Double_t ptDiMu2 = lvDimuon2.Pt();
            Double_t ptDiMu3 = lvDimuon3.Pt();
            Double_t ptDiMu4 = lvDimuon4.Pt();

            if(ptDiMu1>ptDiMu2){
              Double_t temp = ptDiMu1;
              ptDiMu1 = ptDiMu2;
              ptDiMu2 = temp;

              lvTemp = lvDimuon;
              lvDimuon = lvDimuon2;
              lvDimuon2 = lvTemp;
            }

            if(ptDiMu3>ptDiMu4){
              Double_t temp = ptDiMu3;
              ptDiMu3 = ptDiMu4;
              ptDiMu4 = temp;

              lvTemp = lvDimuon3;
              lvDimuon3 = lvDimuon4;
              lvDimuon4 = lvTemp;
            }

            Double_t maDiMu1 = lvDimuon.M();
            Double_t maDiMu2 = lvDimuon2.M();
            Double_t maDiMu3 = lvDimuon3.M();
            Double_t maDiMu4 = lvDimuon4.M();


            Double_t phiDiMu1 = lvDimuon.Phi();
            Double_t phiDiMu2 = lvDimuon2.Phi();
            Double_t phiDiMu3 = lvDimuon3.Phi();
            Double_t phiDiMu4 = lvDimuon4.Phi();

            Double_t plDiMu1 = TMath::Sqrt((lvDimuon.P()*lvDimuon.P())-(lvDimuon.Pt()*lvDimuon.Pt()));
            Double_t plDiMu2 = TMath::Sqrt((lvDimuon2.P()*lvDimuon2.P())-(lvDimuon2.Pt()*lvDimuon2.Pt()));
            Double_t yDiMu1 = TMath::Log((lvDimuon.E()+plDiMu1)/(lvDimuon.E()-plDiMu1))/2;
            Double_t yDiMu2 = TMath::Log((lvDimuon2.E()+plDiMu2)/(lvDimuon2.E()-plDiMu2))/2;

            Double_t plDiMu3 = TMath::Sqrt((lvDimuon3.P()*lvDimuon3.P())-(lvDimuon3.Pt()*lvDimuon3.Pt()));
            Double_t plDiMu4 = TMath::Sqrt((lvDimuon4.P()*lvDimuon4.P())-(lvDimuon4.Pt()*lvDimuon4.Pt()));
            Double_t yDiMu3 = TMath::Log((lvDimuon3.E()+plDiMu3)/(lvDimuon3.E()-plDiMu3))/2;
            Double_t yDiMu4 = TMath::Log((lvDimuon4.E()+plDiMu4)/(lvDimuon4.E()-plDiMu4))/2;




            ( (TH2F*)fOutput->UncheckedAt(kMaDbJpsi) )->Fill(maDiMu1, maDiMu2);
            // ( (TH2F*)fOutput->UncheckedAt(kMaDbJpsi) )->Fill(maDiMu2, maDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kPtDbJpsi) )->Fill(ptDiMu1, ptDiMu2);
            // ( (TH2F*)fOutput->UncheckedAt(kPtDbJpsi) )->Fill(ptDiMu2, ptDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kYDbJpsi) )->Fill(yDiMu1, yDiMu2);
            // ( (TH2F*)fOutput->UncheckedAt(kYDbJpsi) )->Fill(yDiMu2, yDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kPhiDbJpsi) )->Fill(phiDiMu1, phiDiMu2);
            // ( (TH2F*)fOutput->UncheckedAt(kPhiDbJpsi) )->Fill(phiDiMu2, phiDiMu1);

            ( (TH1F*)fOutput->UncheckedAt(kMaJpsi) )->Fill(maDiMu1);
            ( (TH1F*)fOutput->UncheckedAt(kMaJpsi) )->Fill(maDiMu2);



            ( (TH2F*)fOutput->UncheckedAt(kMaDbJpsi) )->Fill(maDiMu3, maDiMu4);
            // ( (TH2F*)fOutput->UncheckedAt(kMaDbJpsi) )->Fill(maDiMu2, maDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kPtDbJpsi) )->Fill(ptDiMu3, ptDiMu4);
            // ( (TH2F*)fOutput->UncheckedAt(kPtDbJpsi) )->Fill(ptDiMu2, ptDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kYDbJpsi) )->Fill(yDiMu3, yDiMu4);
            // ( (TH2F*)fOutput->UncheckedAt(kYDbJpsi) )->Fill(yDiMu2, yDiMu1);

            ( (TH2F*)fOutput->UncheckedAt(kPhiDbJpsi) )->Fill(phiDiMu3, phiDiMu4);
            // ( (TH2F*)fOutput->UncheckedAt(kPhiDbJpsi) )->Fill(phiDiMu2, phiDiMu1);

            ( (TH1F*)fOutput->UncheckedAt(kMaJpsi) )->Fill(maDiMu3);
            ( (TH1F*)fOutput->UncheckedAt(kMaJpsi) )->Fill(maDiMu4);
          }



        }
      }
    }
  }
  switch (counterDiMu){
    case 0:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:0",fCurrentRunNumber,selected.Data()));
      break;
    case 2:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:2",fCurrentRunNumber,selected.Data()));
      break;
    case 4:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:4",fCurrentRunNumber,selected.Data()));
      break;
    case 6:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:6",fCurrentRunNumber,selected.Data()));
      break;
    case 8:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:8",fCurrentRunNumber,selected.Data()));
      break;
    case 10:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:10",fCurrentRunNumber,selected.Data()));
      break;
    default:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Non/numDiMu:Above10",fCurrentRunNumber,selected.Data()));
      break;
  }
  switch (counterMu){
    case 0:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:0/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 1:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:1/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 2:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:2/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 3:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:3/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 4:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:4/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 5:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:5/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 6:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:6/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 7:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:7/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 8:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:8/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 9:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:9/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    case 10:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:10/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
    default:
      fEventCounters->Count(Form("trigger:any/run:%d/selected:%s/numMu:Above10/numDiMu:Non",fCurrentRunNumber,selected.Data()));
      break;
  }
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
      cout<<"Muons number in event "<<endl;
      fEventCounters->Print("numMu/run");
      new TCanvas();
      fEventCounters->Draw("run","numDiMu","selected:yes");
      new TCanvas();
      fEventCounters->Draw("run","numMu","selected:yes");
    }
  }

}
