
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
#include "TFitResult.h"
#include "TFitResultPtr.h"
//PWG includes
#include <AliCounterCollection.h>


//For background
//
// Exeponential background function
Double_t expBg(Double_t *x, Double_t *par) {

//For excluding signal parts
  // if (x[0] > 2.8 && x[0] < 3.3) {
  //     TF1::RejectPoint();
  //     return 0;
  // }
   return exp(par[0]+par[1]*x[0]) ;
}

// Pol2 distribution function
Double_t myPol23(Double_t *x, Double_t *par) {
  // if (x[0] > 2.8 && x[0] < 3.3) {
  //     TF1::RejectPoint();
  //     return 0;
  // }
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])/(par[3]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]);
}
// varWGaus distribution function
Double_t varWGaus(Double_t *x, Double_t *par) {
  // if (x[0] > 2.8 && x[0] < 3.3) {
  //     TF1::RejectPoint();
  //     return 0;
  // }
  return par[0]*exp(-TMath::Power((x[0]-par[1]),2)/(2*TMath::Power((par[2]+par[3]*(x[0]-par[1])/par[1]),2)));
}

Double_t linearFit(Double_t *x, Double_t *par){
  return par[0];
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

enum eList {
  k25FcnPC4 = 0,
  k25FcnPC3 = 1,
  k25FcnPCpp = 2,
  k47FcnPC4 = 3,
  k47FcnPC3 = 4,
  k47FcnPCpp = 5,
  k25FcnVC4 = 6,
  k25FcnVC3 = 7,
  k25FcnVCpp = 8,
  k47FcnVC4 = 9,
  k47FcnVC3 = 10,
  k47FcnVCpp = 11,
  k25FcnEC4 = 12,
  k25FcnEC3 = 13,
  k25FcnECpp = 14,
  k47FcnEC4 = 15,
  k47FcnEC3 = 16,
  k47FcnECpp = 17
};





void DrawAndFitForSingalJ( TString fileName ="AliCombined.root" ){

  TFile *file = TFile::Open( fileName.Data() );

  //get counters
  AliCounterCollection *eventCounters = static_cast<AliCounterCollection*>(file->FindObjectAny("eventCounters"));
  if (!eventCounters) {
    printf("error with eventCounters, exit!\n");
    return;
  }

  new TCanvas();
  eventCounters->Draw("run","trigger","selected:yes");


  Double_t nCMUL7 = eventCounters->GetSum("selected:yes/trigger:CMUL7-B-NOPF-MUFAST");
  printf("number selected = %f\n", nCMUL7);

  //get dimuon pt histogram
  TObjArray *listOfHisto = dynamic_cast<TObjArray*> (file->FindObjectAny("listOfHisto"));
  if (!listOfHisto) {
    printf("error with listOfHisto, exit!\n");
    return;
  }

  TH1F *hMaJpsi = dynamic_cast<TH1F *>( listOfHisto->FindObject("hMaJpsi") );
//histogram for compare each fit method
  const Int_t nMethod = 18;
  TH1F *hNoJpsi = new TH1F("hNoJpsi", "Number of Jpsi", nMethod, 0, nMethod);

  TObjArray *araFunc = new TObjArray(20);
  araFunc->SetOwner();

  TH1F *myhist[20];




  //Function for polinominal and CB2

  //2.2~4.5   1.1e6, 2.83e5, 5.86e4,-11.7,51.1,-26, 3.8
    //GEANT4 1.06,3.23,2.55,1.56
    TF1 *fit25FcnPC4 = new TF1("fit25FcnPC4",polCB,2,5,14);
    fit25FcnPC4->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit25FcnPC4->SetParName(11,  "nL");
    fit25FcnPC4->SetParName(12, "alphaR");
    fit25FcnPC4->SetParName(13, "nR");
    fit25FcnPC4->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,1.06);
    fit25FcnPC4->FixParameter(8, 3.09);
    fit25FcnPC4->FixParameter(9, 0.07);
    fit25FcnPC4->FixParameter(10, 1.06);
    fit25FcnPC4->FixParameter(11, 3.23);
    fit25FcnPC4->FixParameter(12, 2.55);
    fit25FcnPC4->FixParameter(13, 1.56);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit25FcnPC4->SetRange(2.2, 4.5);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit25FcnPC3 = new TF1("fit25FcnPC3",polCB,2,5,14);
    fit25FcnPC3->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit25FcnPC3->SetParName(11,  "nL");
    fit25FcnPC3->SetParName(12, "alphaR");
    fit25FcnPC3->SetParName(13, "nR");
    fit25FcnPC3->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,0.97);
    fit25FcnPC3->FixParameter(8, 3.09);
    fit25FcnPC3->FixParameter(9, 0.07);
    fit25FcnPC3->FixParameter(10, 0.97);
    fit25FcnPC3->FixParameter(11, 3.98);
    fit25FcnPC3->FixParameter(12, 2.3);
    fit25FcnPC3->FixParameter(13, 3.03);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit25FcnPC3->SetRange(2.2, 4.5);

    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit25FcnPCpp = new TF1("fit25FcnPCpp",polCB,2,5,14);
    fit25FcnPCpp->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit25FcnPCpp->SetParName(11,  "nL");
    fit25FcnPCpp->SetParName(12, "alphaR");
    fit25FcnPCpp->SetParName(13, "nR");
    fit25FcnPCpp->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,0.97);
    fit25FcnPCpp->FixParameter(8, 3.09);
    fit25FcnPCpp->FixParameter(9, 0.07);
    fit25FcnPCpp->FixParameter(10, 0.98);
    fit25FcnPCpp->FixParameter(11, 6.97);
    fit25FcnPCpp->FixParameter(12, 1.86);
    fit25FcnPCpp->FixParameter(13, 14.99);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit25FcnPCpp->SetRange(2.2, 4.5);

  //2.4~4.7
    //GEANT4 1.06,3.23,2.55,1.56
    TF1 *fit47FcnPC4 = new TF1("fit47FcnPC4",polCB,2,5,14);
    fit47FcnPC4->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit47FcnPC4->SetParName(11,  "nL");
    fit47FcnPC4->SetParName(12, "alphaR");
    fit47FcnPC4->SetParName(13, "nR");
    fit47FcnPC4->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,0.97);
    fit47FcnPC4->FixParameter(8, 3.09);
    fit47FcnPC4->FixParameter(9, 0.07);
    fit47FcnPC4->FixParameter(10, 1.06);
    fit47FcnPC4->FixParameter(11, 3.23);
    fit47FcnPC4->FixParameter(12, 2.55);
    fit47FcnPC4->FixParameter(13, 1.56);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit47FcnPC4->SetRange(2.4, 4.7);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit47FcnPC3 = new TF1("fit47FcnPC3",polCB,2,5,14);
    fit47FcnPC3->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit47FcnPC3->SetParName(11,  "nL");
    fit47FcnPC3->SetParName(12, "alphaR");
    fit47FcnPC3->SetParName(13, "nR");
    fit47FcnPC3->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,0.97);
    fit47FcnPC3->FixParameter(8, 3.09);
    fit47FcnPC3->FixParameter(9, 0.07);
    fit47FcnPC3->FixParameter(10, 0.97);
    fit47FcnPC3->FixParameter(11, 3.98);
    fit47FcnPC3->FixParameter(12, 2.3);
    fit47FcnPC3->FixParameter(13, 3.03);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit47FcnPC3->SetRange(2.4, 4.7);

    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit47FcnPCpp = new TF1("fit47FcnPCpp",polCB,2,5,14);
    fit47FcnPCpp->SetParNames("p0", "p1", "p2", "p20", "p21", "p22", "p23", "Ns", "miuS", "sigma", "alphaL");
    fit47FcnPCpp->SetParName(11,  "nL");
    fit47FcnPCpp->SetParName(12, "alphaR");
    fit47FcnPCpp->SetParName(13, "nR");
    fit47FcnPCpp->SetParameters(1.83e4, -5.82e3, 4.83e2, 25.2, 35.2, -15.2, 2.31,700,3.09,0.07,0.97);
    fit47FcnPCpp->FixParameter(8, 3.09);
    fit47FcnPCpp->FixParameter(9, 0.07);
    fit47FcnPCpp->FixParameter(10, 0.98);
    fit47FcnPCpp->FixParameter(11, 6.97);
    fit47FcnPCpp->FixParameter(12, 1.86);
    fit47FcnPCpp->FixParameter(13, 14.99);
    fit25FcnPC4->SetParLimits(7, 0, 1e6);
    fit25FcnPC4->SetParLimits(6, 0, 1e6);
    fit47FcnPCpp->SetRange(2.4, 4.7);



  //Function for VWG and CB2


  //2.2~4.5
    //GEANT4 1.06,3.23,2.55,1.56
    TF1 *fit25FcnVC4 = new TF1("fit25FcnVC4",varGausCB,2,5,11);
    fit25FcnVC4->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnVC4->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700,3.09,0.07,0.97,3.98,2.3,3.03);
    fit25FcnVC4->FixParameter(5, 3.09);
    fit25FcnVC4->FixParameter(6, 0.07);
    fit25FcnVC4->FixParameter(7, 1.06);
    fit25FcnVC4->FixParameter(8, 3.23);
    fit25FcnVC4->FixParameter(9, 2.55);
    fit25FcnVC4->FixParameter(10, 1.56);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit25FcnVC4->SetRange(2.2, 4.5);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit25FcnVC3 = new TF1("fit25FcnVC3",varGausCB,2,5,11);
    fit25FcnVC3->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnVC3->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700,3.09,0.07,0.97,3.98,2.3,3.03);
    fit25FcnVC3->FixParameter(5, 3.09);
    fit25FcnVC3->FixParameter(6, 0.07);
    fit25FcnVC3->FixParameter(7, 0.97);
    fit25FcnVC3->FixParameter(8, 3.98);
    fit25FcnVC3->FixParameter(9, 2.3);
    fit25FcnVC3->FixParameter(10, 3.03);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit25FcnVC3->SetRange(2.2, 4.5);

    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit25FcnVCpp = new TF1("fit25FcnVCpp",varGausCB,2,5,11);
    fit25FcnVCpp->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnVCpp->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700,3.09,0.07,0.97,3.98,2.3,3.03);
    fit25FcnVCpp->FixParameter(5, 3.09);
    fit25FcnVCpp->FixParameter(6, 0.07);
    fit25FcnVCpp->FixParameter(7, 0.98);
    fit25FcnVCpp->FixParameter(8, 6.97);
    fit25FcnVCpp->FixParameter(9, 1.86);
    fit25FcnVCpp->FixParameter(10, 14.99);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit25FcnVCpp->SetRange(2.2, 4.5);

  //2.4~4.7
    //GEANT4 1.06,3.23,2.55,1.56
    TF1 *fit47FcnVC4 = new TF1("fit47FcnVC4",varGausCB,2,5,11);
    fit47FcnVC4->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnVC4->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700, 3.09,0.07,0.97,3.98,2.3,3.03);
    fit47FcnVC4->FixParameter(5, 3.09);
    fit47FcnVC4->FixParameter(6, 0.07);
    fit47FcnVC4->FixParameter(7, 1.06);
    fit47FcnVC4->FixParameter(8, 3.23);
    fit47FcnVC4->FixParameter(9, 2.55);
    fit47FcnVC4->FixParameter(10, 1.56);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit47FcnVC4->SetRange(2.4, 4.7);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit47FcnVC3 = new TF1("fit47FcnVC3",varGausCB,2,5,11);
    fit47FcnVC3->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnVC3->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700,3.09,0.07,0.97,3.98,2.3,3.03);
    fit47FcnVC3->FixParameter(5, 3.09);
    fit47FcnVC3->FixParameter(6, 0.07);
    fit47FcnVC3->FixParameter(7, 0.97);
    fit47FcnVC3->FixParameter(8, 3.98);
    fit47FcnVC3->FixParameter(9, 2.3);
    fit47FcnVC3->FixParameter(10, 3.03);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit47FcnVC3->SetRange(2.4, 4.7);

    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit47FcnVCpp = new TF1("fit47FcnVCpp",varGausCB,2,5,11);
    fit47FcnVCpp->SetParNames("Nb", "miuB", "A", "B", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnVCpp->SetParameters(3.27e3, 1.54, -0.721, -0.16, 700,3.09,0.07,0.97,3.98,2.3,3.03);
    fit47FcnVCpp->FixParameter(5, 3.09);
    fit47FcnVCpp->FixParameter(6, 0.07);
    fit47FcnVCpp->FixParameter(7, 0.98);
    fit47FcnVCpp->FixParameter(8, 6.97);
    fit47FcnVCpp->FixParameter(9, 1.86);
    fit47FcnVCpp->FixParameter(10, 14.99);
    fit25FcnPC4->SetParLimits(4, 0, 1e6);
    fit47FcnVCpp->SetRange(2.4, 4.7);

  // With bg of Exeponential


    TF1 *fit25FcnEC4 = new TF1("fit25FcnEC4",expCB,2,5,10);
    fit25FcnEC4->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnEC4->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit25FcnEC4->FixParameter(4, 3.09);
    fit25FcnEC4->FixParameter(5, 0.07);
    fit25FcnEC4->FixParameter(6, 1.06);
    fit25FcnEC4->FixParameter(7, 3.23);
    fit25FcnEC4->FixParameter(8, 2.55);
    fit25FcnEC4->FixParameter(3, 1.56);
    fit25FcnEC4->SetRange(2.2, 4.5);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit25FcnEC3 = new TF1("fit25FcnEC3",expCB,2,5,10);
    fit25FcnEC3->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnEC3->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit25FcnEC3->FixParameter(4, 3.09);
    fit25FcnEC3->FixParameter(5, 0.07);
    fit25FcnEC3->FixParameter(6, 0.97);
    fit25FcnEC3->FixParameter(7, 3.98);
    fit25FcnEC3->FixParameter(8, 2.3);
    fit25FcnEC3->FixParameter(9, 3.03);
    fit25FcnPC4->SetParLimits(3, 0, 1e6);
    fit25FcnEC3->SetRange(2.2, 4.5);

    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit25FcnECpp = new TF1("fit25FcnECpp",expCB,2,5,10);
    fit25FcnECpp->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit25FcnECpp->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit25FcnECpp->FixParameter(4, 3.09);
    fit25FcnECpp->FixParameter(5, 0.07);
    fit25FcnECpp->FixParameter(6, 0.98);
    fit25FcnECpp->FixParameter(7, 6.97);
    fit25FcnECpp->FixParameter(8, 1.86);
    fit25FcnECpp->FixParameter(9, 14.99);
    fit25FcnPC4->SetParLimits(3, 0, 1e6);
    fit25FcnECpp->SetRange(2.2, 4.5);

    //2.4~4.7
    //GEANT4 1.06,3.23,2.55,1.56
    TF1 *fit47FcnEC4 = new TF1("fit47FcnEC4",expCB,2,5,10);
    fit47FcnEC4->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnEC4->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit47FcnEC4->FixParameter(4, 3.09);
    fit47FcnEC4->FixParameter(5, 0.07);
    fit47FcnEC4->FixParameter(6, 1.06);
    fit47FcnEC4->FixParameter(7, 3.23);
    fit47FcnEC4->FixParameter(8, 2.55);
    fit47FcnEC4->FixParameter(9, 1.56);
    fit25FcnPC4->SetParLimits(3, 0, 1e6);
    fit47FcnEC4->SetRange(2.4, 4.7);

    //GEANT3 0.97,3.98,2.3,3.03
    TF1 *fit47FcnEC3 = new TF1("fit47FcnEC3",expCB,2,5,10);
    fit47FcnEC3->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnEC3->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit47FcnEC3->FixParameter(4, 3.09);
    fit47FcnEC3->FixParameter(5, 0.07);
    fit47FcnEC3->FixParameter(6, 0.97);
    fit47FcnEC3->FixParameter(7, 3.98);
    fit47FcnEC3->FixParameter(8, 2.3);
    fit47FcnEC3->FixParameter(9, 3.03);
    fit25FcnPC4->SetParLimits(3, 0, 1e6);
    fit47FcnEC3->SetRange(2.4, 4.7);


    //pp13TeV0.98,6.97,1.86,14.99
    TF1 *fit47FcnECpp = new TF1("fit47FcnECpp",expCB,2,5,10);
    fit47FcnECpp->SetParNames("p0", "p1", "p2", "Ns", "miuS", "sigma", "alphaL", "nL", "alphaR", "nR");
    fit47FcnECpp->SetParameters(10.83, -1.39, -40.42, 700, 3.09, 0.07, 1.06, 3.23, 2.55, 1.56);
    fit47FcnECpp->FixParameter(4, 3.09);
    fit47FcnECpp->FixParameter(5, 0.07);
    fit47FcnECpp->FixParameter(6, 0.98);
    fit47FcnECpp->FixParameter(7, 6.97);
    fit47FcnECpp->FixParameter(8, 1.86);
    fit47FcnECpp->FixParameter(9, 14.99);
    fit25FcnPC4->SetParLimits(3, 0, 1e6);
    fit47FcnECpp->SetRange(2.4, 4.7);



  araFunc->AddAtAndExpand( fit25FcnPC4, k25FcnPC4 );
  araFunc->AddAtAndExpand( fit25FcnPC3, k25FcnPC3 );
  araFunc->AddAtAndExpand( fit25FcnPCpp, k25FcnPCpp );
  araFunc->AddAtAndExpand( fit47FcnPC4, k47FcnPC4 );
  araFunc->AddAtAndExpand( fit47FcnPC3, k47FcnPC3 );
  araFunc->AddAtAndExpand( fit47FcnPCpp, k47FcnPCpp );
  araFunc->AddAtAndExpand( fit25FcnVC4, k25FcnVC4 );
  araFunc->AddAtAndExpand( fit25FcnVC3, k25FcnVC3 );
  araFunc->AddAtAndExpand( fit25FcnVCpp, k25FcnVCpp );
  araFunc->AddAtAndExpand( fit47FcnVC4, k47FcnVC4 );
  araFunc->AddAtAndExpand( fit47FcnVC3, k47FcnVC3 );
  araFunc->AddAtAndExpand( fit47FcnVCpp, k47FcnVCpp );
  araFunc->AddAtAndExpand( fit25FcnEC4, k25FcnEC4 );
  araFunc->AddAtAndExpand( fit25FcnEC3, k25FcnEC3 );
  araFunc->AddAtAndExpand( fit25FcnECpp, k25FcnECpp );
  araFunc->AddAtAndExpand( fit47FcnEC4, k47FcnEC4 );
  araFunc->AddAtAndExpand( fit47FcnEC3, k47FcnEC3 );
  araFunc->AddAtAndExpand( fit47FcnECpp, k47FcnECpp );




  // for Background testing
  TF1 *fitBgE = new TF1("fitBgE",varGausCB,2,5,2);
  fitBgE->SetParameters(1, 1);
  // fitBgE->SetRange(2, 5);

  TF1 *fitBgP = new TF1("fitBgP",myPol23,2,5,7);
  fitBgP->SetParameters(4.65e6, -2.28e6, 9.49e4,-25.5,49.6,-24.8, 4.5);//1.1e6, 2.83e5, 5.86e4,-11.7,51.1,-26, 3.8(2-5)
  fitBgP->SetParLimits(6, 3, 8);
  // fitBgP->SetRange(2.2, 4.5);

  TF1 *fitBgVWG = new TF1("fitBgVWG",varWGaus,2,5,4);
  fitBgVWG->SetParameters(240000, 1.7, -0.7, -0.2);

  if (hMaJpsi) {

    TF1 *fit = new TF1();
    Double_t nS, miuS, sigma, alphaL, nL, alphaR, nR, p0, p1, p2, p20, p21, p22, p23, nB, miuB, bA, bB, chiDiByNDF;
    TF1 *CB2Fit = new TF1("CB2Fit",CrystalBallExtended,2,5,7);

    Double_t araJpsi[nMethod];
    Double_t araBg[nMethod];
    Double_t araErIntegral[nMethod];


    for (Int_t i = 0; i<nMethod; i++){
      myhist[i] = (TH1F*)hMaJpsi->Clone(((TF1*)araFunc->UncheckedAt(i))->GetName());
      gStyle->SetOptFit();
      TFitResultPtr r = myhist[i]->Fit(((TF1*)araFunc->UncheckedAt(i)), "SR");

      myhist[i]->GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
      myhist[i]->GetYaxis()->SetTitle(Form("N"));
      myhist[i]->GetYaxis()->SetTitleOffset(1.2);


      //Get parameters
          fit = myhist[i]->GetFunction( ((TF1*)araFunc->UncheckedAt(i))->GetName() );

          nS = fit->GetParameter("Ns");
          miuS = fit->GetParameter("miuS");
      	  sigma = fit->GetParameter("sigma");
          alphaL = fit->GetParameter("alphaL");
          nL = fit->GetParameter("nL");
          alphaR = fit->GetParameter("alphaR");
          nR = fit->GetParameter("nR");
          chiDiByNDF = fit->GetChisquare()/fit->GetNDF();
          Int_t fitStatus = r;

          CB2Fit->SetParameters(nS, miuS, sigma, alphaL, nL, alphaR, nR);
          CB2Fit->SetLineColor(2);
          CB2Fit->Draw("same");

          // CB2Fit->Draw();

          Int_t nJpsi = (Int_t)(CB2Fit->Integral(2,5)/(3.0/200));
          if (i<6){
            p0 = fit->GetParameter("p0");
            p1 = fit->GetParameter("p1");
            p2 = fit->GetParameter("p2");
            p20 = fit->GetParameter("p20");
            p21 = fit->GetParameter("p21");
            p22 = fit->GetParameter("p22");
            p23 = fit->GetParameter("p23");

            fitBgP->SetParameters(p0, p1, p2, p20, p21, p22, p23);
            fitBgP->SetLineColor(4);
            fitBgP->SetLineStyle(1);
            fitBgP->SetLineWidth(3);
            fitBgP->Draw("same");
            Int_t nBackground = (Int_t)(fitBgP->Integral(miuS-3*sigma, miuS+3*sigma)/(3.0/200));

          }else if (6<=i && i<12){
            nB = fit->GetParameter("Nb");
            miuB = fit->GetParameter("miuB");
            bA = fit->GetParameter("A");
            bB = fit->GetParameter("B");

            fitBgVWG->SetParameters(nB, miuB, bA, bB);
            fitBgVWG->SetLineColor(4);
            fitBgVWG->SetLineStyle(1);
            fitBgVWG->SetLineWidth(3);
            fitBgVWG->Draw("same");
            Int_t nBackground = (Int_t)(fitBgVWG->Integral(miuS-3*sigma, miuS+3*sigma)/(3.0/200));
          }else{
            p0 = fit->GetParameter("p0");
            p1 = fit->GetParameter("p1");
            p2 = fit->GetParameter("p2");

            fitBgE->SetParameters(p0, p1, p2);
            fitBgE->SetLineColor(4);
            fitBgE->SetLineStyle(1);
            fitBgE->SetLineWidth(3);
            fitBgE->Draw("same");

            Int_t nBackground = (Int_t)(fitBgE->Integral(miuS-3*sigma, miuS+3*sigma)/(3.0/200));
          }

          araBg[i] = nBackground;
          araJpsi[i] = nJpsi;
          // printf("%s = %d\n", ((TF1*)araFunc->UncheckedAt(i))->GetName(), nJpsi);
          // printf("%s = %d\n", ((TF1*)araFunc->UncheckedAt(i))->GetName(), nBackground);

          Double_t sigOverB = (Double_t)( (Double_t)araJpsi[i]/(Double_t)araBg[i] ) ;

          // printf("%s = %f\n", ((TF1*)araFunc->UncheckedAt(i))->GetName(), sigOverB);



//Extracting the needed Matrix
          Int_t nPar  = ((TF1*)araFunc->UncheckedAt(i))->GetNpar();
          Double_t erMatrix[49];
          Int_t iErMatrix = 0;


          for (Int_t iMatrix= ( (nPar-7)*nPar )-1+ (nPar-6); iMatrix<nPar*nPar; iMatrix++){
            // printf("gg%d = %f\n",iMatrix, (r->GetCovarianceMatrix()->GetMatrixArray())[iMatrix]);

            erMatrix[iErMatrix] = (r->GetCovarianceMatrix()->GetMatrixArray())[iMatrix];
            iErMatrix++;

            if((iMatrix+1)%nPar==0){
              iMatrix+=nPar-7;
            }
          }
          // for (Int_t iii= 0; iii<50; iii++){
          //   printf("MM%d = %f\n",iii, erMatrix[iii]);
          // }
          Double_t erIntegral = ( CB2Fit->IntegralError(2,5,CB2Fit->GetParameters(), erMatrix)/(3.0/200) );
          araErIntegral[i] = erIntegral;


          TLegend* leg = new TLegend(0.2, 0.4, 0.5, 0.5);
          leg->SetFillStyle(0);
          leg->SetLineColor(0);
          leg->SetTextColor(kBlack);
          leg->SetMargin(0.1);
          // leg->AddEntry((TObject*)0,"Pb-Pb collisions, #sqrt{s_{NN}}=5 TeV", "");//to write something without numerical values for variables
          // leg->AddEntry((TObject*)0,Form("M_{J/#psi} = %3.3f #pm %4.4f Gev/c^{2}",jpsiMean,jpsiMeanError) , "");//to write something and including numerical variables
          leg->AddEntry((TObject*)0,Form("Nsignal/B = %.2f",sigOverB) , "");
          leg->AddEntry((TObject*)0,Form("N of JPsi = %d #pm %.0f",nJpsi, erIntegral) , "");
          leg->Draw();

          if(i == 0){
            c1->Print("c1j.pdf(");
          }else if(i == 17){
            c1->Print("c1j.pdf)");
          }else{
            c1->Print("c1j.pdf");
          }

          // printf("erI = %f\n", erIntegral);


// make a graph in case to see the differences between methods
      if (hMaJpsi) {
        hNoJpsi->GetXaxis()->SetBinLabel(i+1,((TF1*)araFunc->UncheckedAt(i))->GetName());
        hNoJpsi->SetBinContent(i+1, nJpsi);
        hNoJpsi->SetBinError(i+1, erIntegral);
      }


    }
    if (hMaJpsi) {
      Double_t avgJ, sumJ, sqrtJ, rmsJ, statJ, avgError;

      for(Int_t i = 0; i <nMethod; i++){
        //weightting for pp
        if (i%3==2){
          sumJ += araJpsi[i]*2;
          sqrtJ += araJpsi[i]*araJpsi[i]*2;
          avgError += araErIntegral[i]*2;
        }else{
          sumJ += araJpsi[i];
          sqrtJ += araJpsi[i]*araJpsi[i];
          avgError += araErIntegral[i];
        }

      }
      avgJ = sumJ /(nMethod*4/3);
      rmsJ = TMath::Sqrt( sqrtJ/(nMethod*4/3) );
      avgError = avgError/(nMethod*4/3);

//Another loop cause needing in Average
      for(Int_t i = 0; i <nMethod; i++){
        if (i%3==2){
          statJ +=(araJpsi[i]-avgJ)*(araJpsi[i]-avgJ)*2 ;
        }else{
          statJ +=(araJpsi[i]-avgJ)*(araJpsi[i]-avgJ) ;
        }
      }

      statJ = TMath::Sqrt( statJ/(nMethod*4/3) );

      hNoJpsi->Draw("E");

      TLine *avgline = new TLine(0,avgJ,nMethod,avgJ);
      avgline->SetLineColor(kRed);
      avgline->Draw();

      TLine *statline = new TLine(0,avgJ+statJ,nMethod,avgJ+statJ);
      statline->SetLineColor(kRed);
      statline->SetLineStyle(9);
      statline->Draw();

      TLine *statline2 = new TLine(0,avgJ-statJ,nMethod,avgJ-statJ);
      statline2->SetLineColor(kRed);
      statline2->SetLineStyle(9);
      statline2->Draw();


      TLegend* legRe = new TLegend(0.3, 0.33, 0.89, 0.90);
      legRe->SetFillStyle(0);
      legRe->SetLineColor(0);
      legRe->SetTextColor(kBlack);
      legRe->SetMargin(0.1);
      // leg->AddEntry((TObject*)0,"Pb-Pb collisions, #sqrt{s_{NN}}=5 TeV", "");//to write something without numerical values for variables
      // leg->AddEntry((TObject*)0,Form("M_{J/#psi} = %3.3f #pm %4.4f Gev/c^{2}",jpsiMean,jpsiMeanError) , "");//to write something and including numerical variables
      legRe->AddEntry((TObject*)0,Form("N = %.0f (stat) %.2f (sys) %.2f",avgJ, avgError, statJ) , "");
      legRe->Draw();

    }



  }

}
