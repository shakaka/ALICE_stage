
//--------------------------------------------------------------------------
// Base macro for analyzing the analysis task Simple Pt output file.
// Usage: type root (or aliroot) / .L DrawSimplePt.C++ / DrawSimplePt()
// Input: test.root
// Output: no output
//--------------------------------------------------------------------------


// ROOT includes
#include <TFile.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TPaveText.h"
//PWG includes
#include <AliCounterCollection.h>


//For background
//
// Exeponential background function
Double_t expBg(Double_t *x, Double_t *par) {
   return exp(par[0]+par[1]*x[0]) ;
}

// Pol2 distribution function
Double_t myPol23(Double_t *x, Double_t *par) {
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])/(par[3]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]);
}
// varWGaus distribution function
Double_t varWGaus(Double_t *x, Double_t *par) {
  return par[0]*exp(-TMath::Power((x[0]-par[1]),2)/(2*TMath::Power((par[2]+par[3]*(x[0]-par[1])/par[1]),2)));
}



//For signal
//
// Gaussian distribution function(n, miu, sigma)
Double_t gausDis(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*(TMath::Power(((x[0]-par[1])/par[2]),2)));
}

// extended crystall ball - tails on both sides
Double_t CrystalBallExtended(Double_t *x,Double_t *par)
// 0:N,1:mu,2:sigma,3:alphaL,4:nL,5:alphaR,6:nR
// 7 parameters
{


  Double_t t = TMath::Sign(1.,par[3])*TMath::Sign(1.,par[5])*(x[0]-par[1])/par[2];


  Double_t absAlpha_L = fabs((Double_t)par[3]);
  Double_t absAlpha_R = fabs((Double_t)par[5]);

  if ( (t > -absAlpha_L) && (t < absAlpha_R) ) {
    return par[0]*exp(-0.5*t*t);
  }
  else if (t <= -absAlpha_L) {
    Double_t a =  TMath::Power(par[4]/absAlpha_L,par[4])*exp(-0.5*absAlpha_L*absAlpha_L);
    Double_t b= par[4]/absAlpha_L - absAlpha_L;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  else/*if (t >= absAlpha_R) */{
    Double_t c =  TMath::Power(par[6]/absAlpha_R,par[6])*exp(-0.5*absAlpha_R*absAlpha_R);
    Double_t d= par[6]/absAlpha_R - absAlpha_R;

    return par[0]*(c/TMath::Power(d + t, par[6]));
  }

}



// Sum of background and peak function
//
Double_t expGaus(Double_t *x, Double_t *par) {
  return expBg(x,par) + gausDis(x,&par[2]);
}

Double_t mpolGaus(Double_t *x, Double_t *par) {
  return myPol23(x, par) + gausDis(x,&par[7]);
}

Double_t varGausGaus(Double_t *x, Double_t *par) {
  return varWGaus(x,par) + gausDis(x,&par[4]);
}




Double_t expCB(Double_t *x, Double_t *par) {
  return expBg(x,par) + CrystalBallExtended(x,&par[2]);
}

Double_t polCB(Double_t *x, Double_t *par) {
  return myPol23(x, par) + CrystalBallExtended(x,&par[7]);
}

Double_t varGausCB(Double_t *x, Double_t *par) {
  return varWGaus(x,par) + CrystalBallExtended(x,&par[4]);
}




void DrawAndFit( TString fileName ="AliConbined.root" ){

  TFile *file = TFile::Open( fileName.Data() );

  //get counters
  // AliCounterCollection *eventCounters = static_cast<AliCounterCollection*>(file->FindObjectAny("eventCounters"));
  // if (!eventCounters) {
  //   printf("error with eventCounters, exit!\n");
  //   return;
  // }
  //
  // new TCanvas();
  // eventCounters->Draw("run","trigger","selected:yes");

  //get dimuon pt histogram
  TObjArray *listOfHisto = dynamic_cast<TObjArray*> (file->FindObjectAny("listOfHisto"));
  if (!listOfHisto) {
    printf("error with listOfHisto, exit!\n");
    return;
  }

  TH1F *hDiMuPt = dynamic_cast<TH1F *>( listOfHisto->FindObject("hDiMuPt") );
  TH1F *hDiMuY = dynamic_cast<TH1F *>( listOfHisto->FindObject("hDiMuY") );
  TH1F *hDiMuPhi = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuPhi") );
  TH1F *hDiMuM = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuM") );
  TH1F *hEvCen = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hEvCen") );
  TH1F *hDiMuCh = dynamic_cast<TH1F *>(listOfHisto ->FindObject("hDiMuCh") );

  TF1 *fitFcnPC = new TF1("fitFcnPC",polCB,2,5,14);
  fitFcnPC->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
  fitFcnPC->SetParName(11,  "nL");
  fitFcnPC->SetParName(12, "alphaR");
  fitFcnPC->SetParName(13, "nR");
  fitFcnPC->SetParameters(2700000,-950000,95000,-120,170,-71.2,10.4,25000,3.09,0.07,0.97);
  fitFcnPC->SetParameter(11, 3.98);
  fitFcnPC->SetParameter(12, 2.3);
  fitFcnPC->SetParameter(13, 3.03);

  //GEANT4 1.06,3.23,2.55,1.56
  // fitFcnPC->FixParameter(7, 1.06);
  // fitFcnPC->FixParameter(8, 3.23);
  // fitFcnPC->FixParameter(9, 2.55);
  // fitFcnPC->FixParameter(10, 1.56);
  //GEANT3 0.97,3.98,2.3,3.03
  // fitFcnPC->FixParameter(7, 0.97);
  // fitFcnPC->FixParameter(8, 3.98);
  // fitFcnPC->FixParameter(9, 2.3);
  // fitFcnPC->FixParameter(10, 3.03);
  // //pp13TeV0.98,6.97,1.86,14.99
  fitFcnPC->FixParameter(7, 0.98);
  fitFcnPC->FixParameter(8, 6.97);
  fitFcnPC->FixParameter(9, 1.86);
  fitFcnPC->FixParameter(10, 14.99);
  fitFcnPC->SetRange(2, 5);




  //pp13TeV0.98,6.97,1.86,14.99

  TF1 *fitFcnVC = new TF1("fitFcnVC",varGausCB,2,5,11);
  fitFcnVC->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
  fitFcnVC->SetParameters(70000,1.67,-0.7,-0.2, 25000,3.09,0.07,0.97,3.98,2.3,3.03);
  //GEANT4 1.06,3.23,2.55,1.56
  // fitFcnVC->FixParameter(7, 1.06);
  // fitFcnVC->FixParameter(8, 3.23);
  // fitFcnVC->FixParameter(9, 2.55);
  // fitFcnVC->FixParameter(10, 1.56);
  //GEANT3 0.97,3.98,2.3,3.03
  // fitFcnVC->FixParameter(7, 0.97);
  // fitFcnVC->FixParameter(8, 3.98);
  // fitFcnVC->FixParameter(9, 2.3);
  // fitFcnVC->FixParameter(10, 3.03);
  // //pp13TeV0.98,6.97,1.86,14.99
  fitFcnVC->FixParameter(7, 0.98);
  fitFcnVC->FixParameter(8, 6.97);
  fitFcnVC->FixParameter(9, 1.86);
  fitFcnVC->FixParameter(10, 14.99);
  fitFcnVC->SetRange(2, 5);

  TF1 *fitBgE = new TF1("fitBgE",varGausCB,2,5,2);
  fitBgE->SetParameters(1, 1);
  fitBgE->SetRange(2, 5);

  TF1 *fitBgP = new TF1("fitBgP",myPol23,2,5,7);
  fitBgP->SetParameters(-1, 1, -1,-1,-1,1,1);
  fitBgP->SetRange(2, 5);

  TF1 *fitBgVWG = new TF1("fitBgVWG",varGausCB,2,5,4);
  fitBgVWG->SetParameters(-1, 1, 1, 1);
  fitBgVWG->SetRange(2, 5);


  //CB2 for integration
  //GEANT4
  TF1 *CB2G4 = new TF1("CB2G4",CrystalBallExtended,2,5,7);
  CB2G4->SetParameters(6055, 3.098, 0.06891, 1.06,3.23,2.55,1.56);

  Int_t nJpsiG4 = (Int_t)(CB2G4->Integral(2.89127, 3.30473)/(3.0/200));
  printf("nJpsiG4 = %d\n", nJpsiG4);

  //GEANT3
  TF1 *CB2G3 = new TF1("CB2G3",CrystalBallExtended,2,5,7);
  CB2G3->SetParameters(6075, 3.098, 0.06804, 0.97,3.98,2.3,3.03);

  Int_t nJpsiG3 = (Int_t)(CB2G3->Integral(2.89388, 3.30212)/(3.0/200));
  printf("nJpsiG3 = %d\n", nJpsiG3);

  //pp
  TF1 *CB2pp = new TF1("CB2pp",CrystalBallExtended,2,5,7);
  CB2pp->SetParameters(6042, 3.098, 0.06821, 0.98,6.97,1.86,14.99);

  Int_t nJpsipp = (Int_t)(CB2pp->Integral(2.89337, 3.302663)/(3.0/200));
  printf("nJpsipp = %d\n", nJpsipp);


  if (hDiMuPt) {
    new TCanvas();
    hDiMuPt->Draw();
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
    gStyle->SetOptFit();
    // hDiMuM->Fit(fitBgE, "R");
    // hDiMuM->Fit(fitBgP, "R+");
    // hDiMuM->Fit(fitBgVWG, "R+");
    // hDiMuM->Fit(fitFcnEG, "R");
    // hDiMuM->Fit(fitFcnPG, "R+");
    // hDiMuM->Fit(fitFcnVG, "R+");
    // hDiMuM->Fit(fitFcnEC, "R");
    hDiMuM->Fit(fitFcnPC, "R+");
    hDiMuM->Fit(fitFcnVC, "R+");
    hDiMuM->Draw();
    CB2pp->Draw("same");
  }
  if (hEvCen) {
    new TCanvas();
    hEvCen->Draw();
  }
  if (hDiMuCh) {
    new TCanvas();
    hDiMuCh->Draw();
  }

}
