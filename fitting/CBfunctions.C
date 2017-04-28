//Crystal ball function for signal, parameters are 0:Norm, 1:mean, 2:sigma, 3:alpha, 4:n;
// https://root.cern.ch/phpBB3/viewtopic.php?t=7850

Double_t CrystalBall(Double_t *x,Double_t *par){

  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;        // pokud alpha < 0

  Double_t absAlpha = fabs((Double_t)par[3]);

  if (t >= -absAlpha) {
    return par[0]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
} // CB

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
