#include <iostream>
#include <string>
#include <math>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>

void KNdx() {
     const int npoints = 8;
     double angolo_f[npoints] = {60,55,50,45,-40,-45,-55,-60};
     double angolo[npoints];
     
     TCanvas *c1 = new TCanvas(); 
     c1->SetGrid();
     
     double E = (1.17+1.33)/2.;
     
     for (int i = 0; i < npoints; i++) {
         angolo[i] = angolo_f[i];}
         
     double energia[npoints] = {16766,18570,21813,23901,25323,23891,19632,18121};
     double thenergia[npoints];
     
     for (i = 0; i < npoints; i++) {
         thenergia[i] = energia[i];}
     
     double errenergia[npoints], errangolo[npoints];
for (int i = 0; i < npoints; i++) { 
    errenergia[i] = sqrt(thenergia[i]);
    errangolo[i] = 1.3;
}   
// [2]+[0]/((1+1250*cos(x*(4*atan(1)/180.))/[1])*(1+1250*cos(x*(4*atan(1)/180.))/[1]))    *   (1/(1+1250*cos(x*(4*atan(1)/180.))/[1])+(1+1250*cos(x*(4*atan(1)/180.))/[1])-1+cos(x*(4*atan(1)/180.))*cos(x*(4*atan(1)/180.))
    //1252.85
     TF1* fitfunc = new TF1("Fitting Function", "[0]/((1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])*(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1]))*(1/(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])+(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])-1+cos((x-[2])*(4*atan(1)/180.))*cos((x-[2])*(4*atan(1)/180.)))",-90,90);
     fitfunc->SetParameters(10000,511,-5);
     
     TGraphErrors *graph = new TGraphErrors(npoints,angolo_f,thenergia,errangolo,errenergia);
     graph->SetTitle("Klein Nishina");
     graph->GetXaxis()->SetTitle("Angolo(°)");
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
