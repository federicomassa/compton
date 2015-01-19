#include <TF1.h>
#include <TCanvas.h>
#include <cmath>
#include <TGraphErrors.h>
#include <iostream>
using namespace std;

void piccomedio() {
  double max;
  double dist_da_media[90];
  double distanza[90];
  // TF1* gauss = new TF1("gauss","gaus(0)",-200,200);
  // gauss->SetParameters(4200,0,60);
  //TF1* line = new TF1("line","[0]*x",-200,200);
  // TF1* sum = new TF1("sum","[3]*x+[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",-200,200);
 TF1* sum = new TF1("sum","[0]*TMath::Exp(-(x)*(x)/(2*[1]*[1])) + 0.807519*[0]*TMath::Exp(-(x-[2])*(x-[2])/(2*[1]*[1]))",-200,200);
 sum->SetParName(0,"Ampiezza 1");
 sum->SetParName(1,"Sigma");
 sum->SetParName(2,"Distanza");

  for(int dist = 0; dist <= 89; dist++) {
    sum->SetParameters(500,52,double(dist+30));
    sum->SetNpx(10000);
    max = sum->GetMaximumX(-100,100);
    dist_da_media[dist] = max - double(dist+30)/2 ;
    distanza[dist] = double(dist)+30;
  }
  TGraphErrors* dis = new TGraphErrors(90,distanza,dist_da_media,NULL,NULL);
  //   dis->DrawClone("APE");
  
 
  sum->SetParameters(500,52,120);
   sum->Draw();
}
