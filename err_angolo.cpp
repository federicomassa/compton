// Script per simulare l'errore sull'angolo, sia dovuto allo spessore
// dello scintillatore (Sara riporta 2.5 cm)
// sia all'estensione del rivelatore (diametro 5.08 cm, distanza dal plastico 30 cm) --> angolo occupato alpha +/- 4.85 °
// calcolato il punto di intersezione
#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TH1F.h>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TGraph.h>

using namespace std;


double E = 1252.85; // keV
double me = 510.9989; //keV

double invE(double Ef) {
  return 2*TMath::ASin(TMath::Sqrt((1-Ef/E)*me/2/Ef))*180/4/atan(1);
}

double Ef(double theta) {
  return E/(1+E/me*(1-cos(theta)));
}

double P(double theta) {
  return Ef(theta)/E;
}

double Klein_Nishina(double theta) {
  //proporzionale a..  questa è la dsigma/dtheta, non dOmega (seno all'inizio)
  return sin(theta)*P(theta)*P(theta)*(P(theta)+1/P(theta)-1+cos(theta)*cos(theta));
}

void disegno() {
  double kn[10000];
  double angolo[10000];
  
  for (int i = 0; i < 10000; i++) {
    angolo[i] = double(i)/10000*3.14;
    kn[i] = Klein_Nishina(angolo[i]);
  }

  TGraphErrors* graph = new TGraphErrors(10000,angolo,kn,NULL,NULL);
  graph->Draw("AP");
}

void err_angolo() {
  double thE[21];
  double meanE[21];
  double difftheta[21];
  double angolo[21];
  double R = 5.08; //larghezza totale rivelatore
  double random1;
  bool count;
  double E_th;
  double ymax = 2.5; //cm, detto da Sara
  double L = 30;
  double xint = 0, yint = 0;  //x di intersezione
    double alphamax = R/2/L;
    double ralpha = 0;
    double rE = 0;
    int imax = 100000;
  TRandom3 rndgen;
  double y = 0;
  string str;
  string strE;
  str = "hist";
  strE = "histE";
  string strnumber;
  strnumber = "";
  string newstring;
  string newstringE;
  double alpha;
     const char* c;
     const char* cE;
  ostringstream convert;
  TFile* out = new TFile("err_angolo.root","RECREATE");
  // cout << "prima del for alpha" << endl;
  for (int n = 0; n < 21; n++) {
    alpha = (10 + double(n)*2.5)/180*4*atan(1);
    cout << "n: " << n << "e alpha: " << alpha*180/3.1415 << endl;
    //  cout << "dentro for alpha" << endl;
    E_th = E/(1+E/me*(1-cos(alpha)));
    thE[n] = E_th;
    convert.str("");
    convert << round(alpha/(4*atan(1))*180*10);
    strnumber = convert.str();
    newstring = str + strnumber;
    newstringE = strE + strnumber;
     c = newstring.c_str();
     cE = newstringE.c_str();
     
     TH1F *h = new TH1F(c,"Incertezza angolo",100,(alpha-5*alphamax)/(4*atan(1))*180,(alpha+5*alphamax)/(4*atan(1))*180);
     TH1F *histE = new TH1F(cE,"Istogramma energia",200,E_th-200,E_th+200);
     //  cout << "Dichiarato h" << endl;
    for(int i = 0; i < imax; i = i) {
      //  cout << "i: " << i << " alpha: " << alpha << endl;
   
   
   
    // h[n-1] = new TH1F(c,c,100,alpha-3*alphamax,alpha+3*alphamax);
    y = rndgen.Uniform(ymax)-ymax/2;
    do{
    ralpha = alpha + rndgen.Uniform(4*alphamax)-2*alphamax; 
    random1 = rndgen.Uniform(1);
    if(Klein_Nishina(ralpha)/0.54 >= random1) { //massimo visto da funzione disegno() 
      count = 1;}
      else count = 0;
    }
    while (count == 0);
    


    rE = rndgen.Gaus(E/(1+E/me*(1-cos(ralpha))),12.5);
    xint = (L*cos(alpha)+L*tan(alpha)*sin(alpha)-y)/(tan(alpha)+1/tan(ralpha));
    yint = y + 1/tan(ralpha)*xint;
    if(pow(xint-L*sin(alpha),2)+pow(yint-L*cos(alpha),2) <= pow((5.08/2),2)) {h->Fill(ralpha/(4*atan(1))*180); histE->Fill(rE); i++;}
    }
    meanE[n] = histE->GetMean();
    cout << "energia media a n = " << n << ": " << meanE[n] << endl;
    angolo[n] = 10 + double(n)*2.5;
    difftheta[n] = invE(meanE[n])-angolo[n];
    cout << "differenza angolo: " << difftheta[n] << endl;
    // cout << "Prima di scrivere h" << endl;
    h->Write();
    histE->Write();
    // cout << "Scritto h" << endl;
    delete h;
    delete histE;
    //  cout << "Cancellato h" << endl;
  }
  
  TGraph* diff = new TGraph(21,angolo,difftheta);
  diff->SetTitle("Errore commesso sull'angolo; Angolo (gradi); Differenza angolo vero e misurato (gradi)");
  diff->Draw("AP");

  out->Close();  
}
