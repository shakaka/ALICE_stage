
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


void SOBTrend( TString fileName ="NoJpsi.root" ){

  TFile *file = TFile::Open( fileName.Data() );


  Int_t nx = 200;
  const Int_t nProjBin = 20;
  const Int_t nMethod = 18;
  Double_t araNum[nMethod];
  Double_t araErr[nMethod];
  Int_t tmpCounter;
  Double_t araSumX[nProjBin];
  Double_t araSumY[nProjBin];
  Double_t araErrX[nProjBin];
  Double_t araErrY[nProjBin];
  Double_t araSyErX[nProjBin];
  Double_t araSyErY[nProjBin];


  TH1F *hSOBXSum= new TH1F("hSOBXSum", "SOB Trend of X projection", nProjBin, 2, 5);
  TH1F *hSOBYSum= new TH1F("hSOBYSum", "SOB Trend of Y projection", nProjBin, 2, 5);

  TH1F* hSOBX[nProjBin];
  TH1F* hSOBY[nProjBin];

  for(Int_t runs = 0; runs<nProjBin; runs++){
    file->GetObject(Form("hSOBX%d",runs), hSOBX[runs]);
    file->GetObject(Form("hSOBY%d",runs), hSOBY[runs]);

    for (Int_t i =0; i< nMethod; i++){
      araNum[i] = 0;
    }


    tmpCounter =0;
    //get the mean and error of each histogram X
    for (Int_t binNum = 0; binNum < nMethod; binNum++){
      if(hSOBX[runs]->GetBinContent(binNum+1)<1 && hSOBX[runs]->GetBinContent(binNum+1)>0){
        if (binNum%3==2){
          tmpCounter+=2;
          araNum[binNum] = hSOBX[runs]->GetBinContent(binNum+1)*2;

          // printf("araNumX = %f\n",  araNum[binNum]);
        }else{
          tmpCounter++;
          araNum[binNum] = hSOBX[runs]->GetBinContent(binNum+1);
        }
      }
    }
    for (Int_t i = 0; i<nMethod; i++){
      if(araNum[i] != 0)
        araSumX[runs] += araNum[i];
    }
    if (tmpCounter!=0){
      araSumX[runs] = araSumX[runs]/tmpCounter;
    }
    for (Int_t binNum = 0; binNum<nMethod; binNum++){
      if (araNum[binNum]>0){
        araSyErX[runs] += (araSumX[runs]-araNum[binNum])*(araSumX[runs]-araNum[binNum]);
      }
    }
    if (tmpCounter!=0){
      araSyErX[runs] = TMath::Sqrt(araSyErX[runs]/tmpCounter);
    }
    hSOBXSum->SetBinContent(runs+1, araSumX[runs]);
    hSOBXSum->SetBinError(runs+1, araSyErX[runs]);


    //get the mean and error of each histogram Y
    for (Int_t i =0; i< nMethod; i++){
      araNum[i] = 0;
    }
    tmpCounter = 0; // reset counter
    for (Int_t binNum = 0; binNum < nMethod; binNum++){
      if(hSOBY[runs]->GetBinContent(binNum+1)<1 && hSOBY[runs]->GetBinContent(binNum+1)>0 ){
        if (binNum%3==2){
          tmpCounter+=2;
          araNum[binNum] = hSOBY[runs]->GetBinContent(binNum+1)*2;

          // printf("araNumYD = %f\n",  araNum[binNum]);
        }else{
          tmpCounter++;
          araNum[binNum] = hSOBY[runs]->GetBinContent(binNum+1);

          // printf("araNumY = %f\n",  araNum[binNum]);
        }
      }
    }
    for (Int_t i = 0; i<nMethod; i++){
      if(araNum[i] != 0)
        araSumY[runs] += araNum[i];
      // printf("araNumY = %f\n",  araNum[i]);
    }
    // printf("counter = %d\n",  tmpCounter);
    if (tmpCounter!=0){
      araSumY[runs] = araSumY[runs]/tmpCounter;
    }
    for (Int_t binNum = 0; binNum<nMethod; binNum++){
      if (araNum[binNum]>0){
        araSyErY[runs] += (araSumY[runs]-araNum[binNum])*(araSumY[runs]-araNum[binNum]);
      }
    }
    if (tmpCounter!=0){
      araSyErY[runs] = TMath::Sqrt(araSyErY[runs]/tmpCounter);
    }

    hSOBYSum->SetBinContent(runs+1, araSumY[runs]);
    hSOBYSum->SetBinError(runs+1, araSyErY[runs]);

    }
    hSOBXSum->Draw("L");
    c1->Print("hSOBSum.pdf(");
    hSOBYSum->Draw("L");
    c1->Print("hSOBSum.pdf)");
}
