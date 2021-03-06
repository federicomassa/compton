#include <TGraphErrors.h>
#include <string>
#include <sstream>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

void display_spettro() {
  TCanvas* c1 = new TCanvas();
  const char name1[] = {"Egate+50.dat"};
  TGraphErrors* hist = new TGraphErrors(name1,"%lg %lg");
  // string str;
  // double integral;
  // int fill;
  // int i = 0;
  // do {
  //   i++;
  //   getline(in,str, '\t');
  //   getline(in,str); //solo per gli spettri in energia
  //   //    if (i <= 30) cout << "i= " << i << "..." << str << endl;
    
  //   stringstream(str) >> fill;
  //   hist->SetBinContent(i,double(fill));
  //   str = "";
  //   fill = 0;
  // } while(!in.eof());

  double* x = hist->GetX();
  double* y = hist->GetY();
  
  double xv[820];
  double yv[820];
  double errX[820];
  double errY[820];

  for (int i = 0; i < 820; i++) {
    xv[i] = x[i];
    yv[i] = y[i];
    errX[i] = 12.5;
    errY[i] = TMath::Sqrt(y[i]);
  }
  TGraphErrors* hist1 = new TGraphErrors(820,xv,yv,NULL,errY);

   hist1->SetFillStyle(3001);
   hist1->SetFillColor(kGreen);
   hist1->SetNameTitle("gate","Gate");
   
  hist1->SetTitle("Gate, x gradi;Energia (keV); Conteggi normalizzati");
  c1->cd();
  hist1->Draw("AF");


  double Emin = 550;
  double Emax = 750;

  TF1* fitfunc = new TF1("Fitting Function", "gaus(0)",Emin,Emax); //singola gaus
  // TF1* fitfunc = new TF1("Fitting Function", "gaus(0)+gaus(3)",Emin,Emax);
  
  fitfunc->SetParameters(500,(Emin+Emax)/2,50);
  // fitfunc->SetParameters(500,(Emin+Emax)/2-20,30,500,(Emin+Emax)/2+20,30);
  hist1->Fit(fitfunc,"0RS");
  cout << "NDF: " << fitfunc->GetNDF() << endl;
  fitfunc->SetRange(fitfunc->GetParameter(1)-3*fitfunc->GetParameter(2),fitfunc->GetParameter(1)+3*fitfunc->GetParameter(2));
  fitfunc->SetLineColor(kRed);
  fitfunc->SetLineStyle(7);
  fitfunc->DrawClone("same");

  
  // Decommentare per fare l'integrale del fotopicco

  // TF1* fitfunc = new TF1("Fitting Function","gaus(0)",520,600);
  // fitfunc->SetParameters(400,550,50);
  // hist->Fit(fitfunc,"0RS");
  // fitfunc->SetLineColor(kRed);
  // fitfunc->SetRange(0,819);
  // integral = fitfunc->Integral(0,819);
  // cout << integral << endl;
  // fitfunc->DrawClone("SAME");


}
