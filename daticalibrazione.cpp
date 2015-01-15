#include <string>
#include <cmath>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TLine.h>



TFile* out = new TFile("calib_pseudo.root","RECREATE");
int imax = 10000;
double ry;
 double varerr(double x,TFitResultPtr fit){ 
 double xerr = 0; //approx media errori ch (80*80)
 double varpar[4];
  for (int i = 0; i<4;i++){
    varpar[i] = fit->CovMatrix(i,i);}
  return(varpar[0]+(varpar[1]*xerr+fit->Parameter(1)*xerr+x*x*varpar[1])+(varpar[2]*4*x*x*xerr+varpar[2]*pow((xerr+x*x),2)+4*x*x*xerr*pow(fit->Parameter(2),2))+(varpar[3]*9*x*x*x*x*xerr+varpar[3]*pow(x,6)+pow(fit->Parameter(3),2)*9*x*x*x*x*xerr));}

 double coverr(double x,TFitResultPtr fit){
   double xerr = 0; //approx media errori ch (80*80)
   return( fit->CovMatrix(0,1)*x+fit->CovMatrix(0,2)*(xerr+x*x)+fit->CovMatrix(0,3)*x*x*x+fit->CovMatrix(1,2)*x*x*x+fit->CovMatrix(1,3)*x*x*x*x+fit->CovMatrix(2,3)*x*x*x*x*x);
 }

/*  double err2(double x1,TFitResultPtr fit){
  
       return(varerr(x1,fit)+2*coverr(x1,fit));
       return r;} */

void daticalibrazione(){
  TCanvas *c1 = new TCanvas("c1","c1",1000,700);
  double errE[78];
double chgraph[78];
 double meanE[78];
  //      TCanvas *c1 = new TCanvas();
      TRandom3 rndgen;
      // string sources[7] = {"Am1","Am2","Co1","Co2","Na1","Na2","Cs"};
    // double chval[5] = {7002.736,7314.4,4793.13,7126.22,5580.17};
//     double enval[5] = {1173.2, 1332.5,511.0034,1274.5,661.64}; //keV, dubbio su Am1, 14 keV?
//     double cherr[5] = {63.43949,60.38217,92.03397,59.40552,82.03822}; //FWHM/2
//     double enerr[5] = {0.5,0.5,0.0005,0.5,0.05};
//     TGraphErrors *graph = new TGraphErrors(5,chval,enval,cherr,enerr);
   

   double chval[7] = {525.74,1048.9,7002.736,7314.4,4793.13,7126.22,5580.17};
   double enval[7] = {18, 60, 1173.2, 1332.5,511.0034,1274.5,661.64}; //keV, dubbio su Am1, 14 keV?
   //  double cherr[7] = {81.65605,61.1465,63.43949,60.38217,92.03397,59.40552,82.03822}; //FWHM/2
   double cherr[7] = {69.34696,51.9291,53.8764,51.2800,78.16,50.4505,69.6715}; 
   double rchval[7]; 

 TGraphErrors *meas_graph = new TGraphErrors(7,chval,enval,NULL,NULL);

   for (int ch = 100; ch <= 7800; ch = ch + 100) {
     cout << "Canale: " << ch << endl;
     TH1F *hy = new TH1F("hy","E histogram",100,-33+0.104*ch-1.8E-5*ch*ch+4.03E-9*ch*ch*ch-100,-33+0.104*ch-1.8E-5*ch*ch+4.03E-9*ch*ch*ch+100);
   for (int i = 0; i < imax; i++) {
   
   
   
   for (int n = 0; n < 7; n++) {rchval[n] = rndgen.Gaus(chval[n], cherr[n]);}
  // TF1* noam = new TF1("No_Am", "474.629+0.029255*x-3.66635E-5*x*x+6.70617E-9*x*x*x",0,8192);
//   noam->SetLineStyle(2);
//   noam->SetLineColor(kOrange);
//   noam->DrawClone("SAME");
    TGraphErrors *graph = new TGraphErrors(7,rchval,enval,cherr,NULL);
     graph->SetTitle("Calibrazione;Canale;Energia(keV)");
     graph->GetYaxis()->SetTitleOffset(1.4);
     graph->SetLineColor(kRed);
     graph->SetMarkerColor(kRed);
     graph->SetMarkerStyle(20);
     graph->SetMarkerSize(0.8);
     // graph->DrawClone("APE"); 
     
     //TF1 *fitfunc = new TF1("Fitting Function", "[0]+[1]*x+[2]*x*x", 0,8192);    
         TF1 *fitfunc = new TF1("Fitting Function", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0,8192);
// TF1 *fitfunc = new TF1("Fitting Function", "[0]+[1]*x+[2]*x*x+4.03156E-9*x*x*x", 0,8192);
//     TF1 *fitfunc = new TF1("Cubic Fit", "-33.0933+[0]*x-1.80489E-5*x*x+4.03119E-9*x*x*x", 0,8192);
	 fitfunc->SetParameters(-30,0.1,-1E-5,4E-9);
//     fitfunc->SetParameter(0,0.1);
	 //   fitfunc->SetParameters(-30,0.1,1E-5);
     graph->Fit(fitfunc,"0QRS");
     ry = fitfunc->GetParameter(0)+fitfunc->GetParameter(1)*ch+fitfunc->GetParameter(2)*ch*ch+fitfunc->GetParameter(3)*ch*ch*ch;
     hy->Fill(ry);
     delete graph;
     delete fitfunc;
   }
     errE[int(double((ch-100))/100)] = hy->GetRMS(); //assoluto
     //errE[int(double((ch-100))/100)] = hy->GetRMS()/(hy->GetMean()); //relativo 
   meanE[int(double((ch-100))/100)] = hy->GetMean(); 
   chgraph[int(double((ch-100))/100)] = double(ch);
   delete hy;
   }
   c1->cd();
   TGraphErrors *expected = new TGraphErrors(78,chgraph,meanE,NULL,errE);
   expected->SetFillColor(kYellow);
   expected->SetTitle("Curva di calibrazione;Canali;Energia (keV)");
   expected->DrawClone("E4AL");
   meas_graph->DrawClone("PESame");

   TCanvas *c2 = new TCanvas("c2","graph",1000,700);
   c2->cd();
   TGraphErrors *errgraph = new TGraphErrors(78,chgraph,errE,NULL,NULL);
   errgraph->SetTitle("Incertezza sull'energia;Canale;Incertezza energia (keV)");
   errgraph->SetMarkerStyle(kFullDotSmall);
   errgraph->SetMarkerColor(kRed);
   errgraph->SetMarkerSize(3);
   errgraph->DrawClone("APE");
   TLine* line1 = new TLine(525.74,0,525.74,30);
   TLine* line2 = new TLine(1048.9,0,1048.9,30);
   TLine* line3 = new TLine(7002.736,0,7002.736,30);
   TLine* line4 = new TLine(7314.4,0,7314.4,30);
   TLine* line5 = new TLine(4793.13,0,4793.13,30);
   TLine* line6 = new TLine(7126.22,0,7126.22,30);
   TLine* line7 = new TLine(5580.17,0,5580.17,30);
 line1->Draw();
 line2->Draw();
 line3->Draw();
 line4->Draw();
 line5->Draw();
 line6->Draw();
 line7->Draw();
   errgraph->Write();
   out->Close();
     }             
     
     

