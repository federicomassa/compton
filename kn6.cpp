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


void KN6() {
  //   const int npoints = 6;
     const int npoints = 8;
     double angolo_f[npoints] = {60,55,50,45,-40,-45,-55,-60};
     //       double angolo_f[npoints] = {50,45,-40,-45,-55,-60};
       //  double angolo_f[npoints] = {50,45,-45,-55,-60};
     double angolo[npoints];
     
     TCanvas *c1 = new TCanvas(); 
     c1->SetGrid();
     
     double E = (1.17+1.33)/2.;
     
     for (int i = 0; i < npoints; i++) {
         angolo[i] = angolo_f[i]-1.38;}
     double energia[npoints] = {16766/0.9556/0.3858,18570/0.958/0.3475,21813/0.96/0.3122,23901/0.9464/0.2804,25323/0.8533/0.2519,23891/0.9464/0.2804,19632/0.958/0.3475,18121/0.9558/0.3858};
     //  double energia[npoints] = {21813,23901,25323,23891,19632,18121};
     // double energia[npoints] = {21813/0.96/0.3122,23901/0.9464/0.2804,25323/0.8533/0.2519, 23891/0.9464/0.2804,19632/0.958/0.3475,18121/0.9558/0.3858};
     //     double energia[npoints] = {21813,23901,25323,23891,19632,18121};
     double thenergia[npoints];
     
     for (int i = 0; i < npoints; i++) {
       thenergia[i] = Ef(2.4518,angolo[i]);}
         double integral[npoints] = {297536*pow(thenergia[0],2.23)/0.49,301288*pow(thenergia[1],2.23)/0.532,368041*pow(thenergia[2],2.23)/0.512,185865*pow(thenergia[3],2.23)/0.48,111006*pow(thenergia[4],2.23)/0.47};
     //    double integral[npoints] = {23.8*pow(thenergia[0],3),26.12*pow(thenergia[1],3),25.26*pow(thenergia[2],3),17.39*pow(thenergia[3],3),12.80*pow(thenergia[4],3)};
     double errenergia[npoints], errangolo[npoints];
for (int i = 0; i < npoints; i++) { 
    errenergia[i] = sqrt(energia[i]);
    errangolo[i] = 2;
}   
// [2]+[0]/((1+1250*cos(x*(4*atan(1)/180.))/[1])*(1+1250*cos(x*(4*atan(1)/180.))/[1]))    *   (1/(1+1250*cos(x*(4*atan(1)/180.))/[1])+(1+1250*cos(x*(4*atan(1)/180.))/[1])-1+cos(x*(4*atan(1)/180.))*cos(x*(4*atan(1)/180.))
    //1252.85
//    TF1* fitfunc = new TF1("Fitting Function", "[0]/((1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])*(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1]))*(1/(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])+(1+1250*(1-cos((x-[2])*(4*atan(1)/180.)))/[1])-1+cos((x-[2])*(4*atan(1)/180.))*cos((x-[2])*(4*atan(1)/180.)))",-90,90);
 TF1 *fitfunc = new TF1("Fitting Function", Klein_Nishina, -90,90,3);

     fitfunc->SetParameters(135000,2,-2);
     //     fitfunc->FixParameter(1,500);
     TGraphErrors *graph = new TGraphErrors(npoints,angolo_f,energia,errangolo,errenergia);
     graph->SetTitle("Klein Nishina");
     graph->GetXaxis()->SetTitle("Angolo(°)");
     graph->GetYaxis()->SetTitle("# Conteggi Inorganico");
     graph->GetYaxis()->SetTitleOffset(1.4);
     graph->DrawClone("APE");
     
     TFitResultPtr frp = graph->Fit(fitfunc,"0RS");
     frp->Print();
     fitfunc->SetLineColor(kGreen);
     fitfunc->DrawClone("SAME");
     

      fitfunc->SetParameters(4E4,2,-2);
      fitfunc->SetLineColor(kRed);
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
