
//--------------------------------------------------------------------------
// Base macro for analyzing the analysis task Simple Pt output file.
// Usage: type root (or aliroot) / .L DrawSimplePt.C++ / DrawSimplePt()
// Input: test.root
// Output: no output
//--------------------------------------------------------------------------


// ROOT includes
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>

//PWG includes
#include <AliCounterCollection.h>

void DrawSimplePt( TString fileName ="testMuonCharge.root" ){
  //testSimplePtMuonTask.root
  //AnalysisResults_run000245346.root
  TFile *file = TFile::Open( fileName.Data() );

  //get counters
  AliCounterCollection *eventCounters = static_cast<AliCounterCollection*>(file->FindObjectAny("eventCounters"));
  if (!eventCounters) {
    printf("error with eventCounters, exit!\n");
    return;
  }

  new TCanvas();
  eventCounters->Draw("run","trigger","selected:yes");

  //get dimuon pt histogram
  TObjArray *listOfHisto = dynamic_cast<TObjArray*> (file->FindObjectAny("listOfHisto"));
  if (!listOfHisto) {
    printf("error with listOfHisto, exit!\n");
    return;
  }

  TH1F *hDimuonPt = dynamic_cast<TH1F *>( listOfHisto->FindObject("hDimuonPt") );
  TH1F *hEventZv = dynamic_cast<TH1F *>( listOfHisto->FindObject("hEventZv") );
  TH1F *hMuonEta = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hMuonEta") );
  TH1F *hMuonRab = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hMuonRab") );
  TH1F *hMuonPt = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hMuonPt") );
  TH1F *hMuonPhi = dynamic_cast<TH1F *>( listOfHisto->FindObject("hMuonPhi") );
  TH1F *hDiMuY = dynamic_cast<TH1F *>( listOfHisto->FindObject("hDiMuY") );
  TH1F *hDiMuPhi = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuPhi") );
  TH1F *hDiMuM = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuM") );
  TH1F *hDiMuZv = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuZv") );
  TH1F *hDiMuCh = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuCh") );


  if (hDimuonPt) {
    new TCanvas();
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
  }

}
