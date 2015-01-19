#include <TH1F.h>
#include <fstream>
#include <string>
#include <sstream>
#include <TF1.h>

using namespace std;

void display_spettro() {
  ifstream in("gate-15.dat");
  TH1F* hist = new TH1F("hist","Istogramma",819,-32.523,1818.08);
  string str;
  double integral;
  int fill;
  int i = 0;
  do {
    i++;
    getline(in,str);
    stringstream(str) >> fill;
    hist->SetBinContent(i,double(fill));
    str = "";
    fill = 0;
  } while(!in.eof());
  hist->Draw();


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
