#include <TGraph.h>
#include <string>
#include <sstream>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TPad.h>
#include <TMultiGraph.h>
#include <TLegend.h>

using namespace std;

void display_spettro() {
  TCanvas* c1 = new TCanvas();
  TMultiGraph* mg = new TMultiGraph();
  const char name1[] = {"Egate-30.dat"};
  const char name2[] = {"Enogate-30.dat"};
  TGraph* hist1 = new TGraph(name1,"%lg %lg");
  TGraph* hist2 = new TGraph(name2,"%lg %lg");
  TGraph* hist2_norm;
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
   hist1->SetFillStyle(3001);
   hist1->SetFillColor(kGreen);
   hist1->GetXaxis()->SetRangeUser(0,1500);
   hist1->GetYaxis()->SetRangeUser(0,5000);
   hist1->SetNameTitle("gate","Gate");
   double int1 = hist1->Integral(-34,1820);
   double int2 = hist2->Integral(-34,1820);
   double norm = int2/int1;
   
   double *x = hist2->GetX();
   double *y = hist2->GetY();
   for (int i = 0; i < 820; i++) y[i] = y[i]/norm;
   
   hist2_norm = new TGraph(820,x,y);

    hist2_norm->SetFillStyle(3002);
    hist2_norm->SetFillColor(kRed);
    hist2_norm->GetXaxis()->SetRangeUser(0,1500);
    hist2_norm->GetYaxis()->SetRangeUser(0,5000);
    hist2_norm->SetNameTitle("nogate","Senza gate");
    mg->Add(hist1);
    mg->Add(hist2_norm);
    c1->cd();
    mg->SetTitle("Confronto spettri con e senza gate: 30 gradi;Energia (keV); Conteggi normalizzati");
    // TLegend* lg = new TLegend();
    // lg->AddEntry("hist1","Gate","f");
    // lg->AddEntry("hist2_norm","Senza gate","f");
    mg->Draw("AF");
    c1->cd();
    // lg->Draw();
    c1->BuildLegend();
    c1->SetGridx();
    c1->Update();    
    

    




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
