#include <iostream>
#include <string>
#include <cmath>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TAxis.h>
double me = 511; //keV
double r0 = 2.81794; //fm 

double Ef(double R, double theta) {
  return me*R/(1+R*(1-cos(theta*4*atan(1)/180)));
}

double P(double R, double theta) {
  return Ef(R, theta)/me/R;
}

double Klein_Nishina(double *vars, double *pars) {
  //proporzionale a
  return pars[0]*P(pars[1],vars[0]-pars[2])*P(pars[1],vars[0]-pars[2])*(P(pars[1],vars[0]-pars[2])+1/P(pars[1],vars[0]-pars[2])-1+cos((vars[0]-pars[2])*4*atan(1)/180)*cos((vars[0]-pars[2])*4*atan(1)/180));
}


void KN() {
     const int npoints = 8;
     // double angolo_f[npoints] = {-60,-55,-45,-40,30,40,55,60};
      double angolo_f[npoints] = {-60,-55,-45,-40,30,40,55,60};
     double angolo[npoints];
     
     TCanvas *c1 = new TCanvas(); 
     c1->SetGrid();
     
     double E = (1.17+1.33)/2.;
     
     for (int i = 0; i < npoints; i++) {
         angolo[i] = angolo_f[i]-1.38;}
         
     double energia[npoints] = {519286,542943,577561,595448,550457,514355,395066,376093};
     double thenergia[npoints];
     
     for (int i = 0; i < npoints; i++) {
         thenergia[i] = energia[i];}
     
     double errenergia[npoints], errangolo[npoints];
for (int i = 0; i < npoints; i++) { 
    errenergia[i] = sqrt(thenergia[i]);
    errangolo[i] = 1.3;
}   
// [2]+[0]/((1+1250*cos(x*(4*atan(1)/180.))/[1])*(1+1250*cos(x*(4*atan(1)/180.))/[1]))    *   (1/(1+1250*cos(x*(4*atan(1)/180.))/[1])+(1+1250*cos(x*(4*atan(1)/180.))/[1])-1+cos(x*(4*atan(1)/180.))*cos(x*(4*atan(1)/180.))
    //1252.85
     // TF1* fitfunc = new TF1("Fitting Function", "[0]/((1+1250*(1-cos(x*(4*atan(1)/180.)))/[1])*(1+1250*cos(x*(4*atan(1)/180.))/[1]))*(1/(1+1250*cos(x*(4*atan(1)/180.))/[1])+(1+1250*cos(x*(4*atan(1)/180.))/[1])-1+cos(x*(4*atan(1)/180.))*cos(x*(4*atan(1)/180.)))",-90,90);
 TF1*  fitfunc = new TF1("Fitting Function", Klein_Nishina,-90,90,3); 

 fitfunc->SetParameters(110000,0.5,1.3);
     
     TGraphErrors *graph = new TGraphErrors(npoints,angolo_f,thenergia,errangolo,errenergia);
     graph->SetTitle("Klein Nishina");
     graph->GetXaxis()->SetTitle("Angolo(gradi)");
     graph->GetYaxis()->SetTitle("# Conteggi Inorganico");
     graph->GetYaxis()->SetTitleOffset(1.4);
     graph->DrawClone("APE");
     
     
     TFitResultPtr frp = graph->Fit(fitfunc,"0RS");
     frp->Print();
     fitfunc->SetLineColor(kGreen);
     fitfunc->DrawClone("SAME");
     
    // TPaveText* pt = new TPaveText(.05,.1,.95.,.8);
//     pt->AddText("Chi quadro ridotto: ");
//     pt->AddText("fitfunc->GetChisquare()/fitfunc->GetNDF()");
//     pt->DrawClone("same");

       //TF1* th = new TF1("Theoretical", "1119.91/(1+1119.91/[0]*(1-cos((x-[1])*4*atan(1.)/180.)))",-90,90);
//       th->SetParameters(510.9989,1.38);
//       th->SetLineColor(2);
//       th->DrawClone("same");
//               
     }
