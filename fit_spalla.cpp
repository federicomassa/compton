#include <TGraphErrors.h>
#include <iostream>
#include <TF1.h>

using namespace std;

void fit_spalla() {
  //  double angolo[15] = {-60,-55,-50,-45,-40,-35,-30,-25,30,35,40,45,50,55,60};
  //  double errangolo[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  // double picchi[15] = {500,569,607,660,705,805,838,946,911,829,765,685,617,540,515};
  // double spalle[15] = {275,350,370,410,450,510,530,630,650,550,480,425,380,325,300};
  // double diff[15];
  // double errdiff[15] = {28,28,28,28,28,20,20,28,28,20,28,20,20,20,20};

   double angolo[15] = {-60,-55,-50,-45,-40,-35,-30,-25,30,35,40,45,50,55,60};
   double errangolo[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
   double picchi[15] = {500,569,607,660,705,805,838,946,911,829,765,685,617,540,515};
   double corr_picchi[15] = {8,7.5,12,10.2,11.16,-27,-37,-40,-37,-27,11.16,10.2,12,7.5,8};
   double picchi_corretti[15];
   double simul_dist[15] =  {47,52,68,70,83,40.6,40,40,40,40.6,83,70,68,52,47};
   //   double spalle[15] = {275,350,370,410,450,510,530,630,650,550,480,425,380,325,300};
double spalle[15] = {275,350,370,410,450,525,550,650,650,550,480,425,380,325,300};  
   double diff[15];
   double errdiff[15] = {28,28,28,28,28,20,20,28,28,20,28,20,20,20,20};

 for (int i = 0; i < 15; i++) {
   picchi_corretti[i] = picchi[i] + corr_picchi[i];}
 
 for (int i = 0; i < 15; i++) {
    diff[i] = picchi_corretti[i] - spalle[i]-simul_dist[i];}


  double errpicchi[15] = {12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5};
  TGraphErrors* graph = new TGraphErrors(15,angolo,diff,errangolo,errdiff);
  TF1* fitfunc = new TF1("Fitting Function","1252.85/(1+1252.82/[0]*(1-cos((x-[1])*4*atan(1)/180)))/(1+2/[0]*1252.85/(1+1252.82/[0]*(1-cos((x-[1])*4*atan(1)/180))))",-65,65);

  fitfunc->SetParName(0,"Mass");
  fitfunc->SetParName(1,"Angular Shift");
  fitfunc->SetParameters(500,1.5);
  
  graph->Fit(fitfunc,"0RS");
  //  Grafico Fit differenze-angolo
  graph->SetTitle("Distanza spalla-fotopicchi;Angolo (gradi); Distanza spalla-fotopicchi (keV)");
  graph->DrawClone("APE");
  fitfunc->DrawClone("SAME");
  fitfunc->SetParameters(511,-3.5);
  fitfunc->SetLineColor(kGreen);
  fitfunc->DrawClone("SAME");

  TGraphErrors* graph2 = new TGraphErrors(15,picchi_corretti,diff,errpicchi,errdiff);
  TF1* fitfunc2 = new TF1("Fitting Function","x/(1+2*x/[0])",450,1000);
  fitfunc2->SetParameter(0,500);
graph2->Fit(fitfunc2,"0RS");

// Grafico differenze vs picchi
  // graph2->DrawClone("APE");
  // fitfunc2->DrawClone("SAME");
  // fitfunc2->SetParameter(0,511);
  // fitfunc2->SetLineColor(kRed);
  // fitfunc2->DrawClone("SAME");
}
