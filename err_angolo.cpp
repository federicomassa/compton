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

using namespace std;

void err_angolo() {
  double R = 5.08;
  double ymax = 2.5; //cm, detto da Sara
  double L = 30;
  double xint = 0, yint = 0;  //x di intersezione
    double alphamax = R/2/L;
    double ralpha = 0;
    int imax = 100000;
  TRandom3 rndgen;
  double y = 0;
  string str;
  str = "hist";
  string strnumber;
  strnumber = "";
  string newstring;
     const char* c;
  ostringstream convert;
  TFile* out = new TFile("err_angolo.root","RECREATE");
  // cout << "prima del for alpha" << endl;
  for (double alpha = 2.5/180*4*atan(1); alpha <= double(60.0/180.0*4.0*atan(1)); alpha = alpha + 2.5/180*4*atan(1)) {
    //  cout << "dentro for alpha" << endl;
    convert.str("");
    convert << round(alpha/(4*atan(1))*180*10);
    strnumber = convert.str();
    newstring = str + strnumber;
     c = newstring.c_str();
     
     TH1F *h = new TH1F(c,"Incertezza angolo",100,(alpha-5*alphamax)/(4*atan(1))*180,(alpha+5*alphamax)/(4*atan(1))*180);
     //  cout << "Dichiarato h" << endl;
    for(int i = 0; i < imax; i = i) {
      //  cout << "i: " << i << " alpha: " << alpha << endl;
   
   
   
    // h[n-1] = new TH1F(c,c,100,alpha-3*alphamax,alpha+3*alphamax);
    y = rndgen.Uniform(ymax)-ymax/2;
    ralpha = alpha + rndgen.Uniform(8*alphamax)-4*alphamax;
    xint = (L*cos(alpha)+L*tan(alpha)*sin(alpha)-y)/(tan(alpha)+1/tan(ralpha));
    yint = y + 1/tan(ralpha)*xint;
    if(pow(xint-L*sin(alpha),2)+pow(yint-L*cos(alpha),2) <= pow((5.08/2),2)) {h->Fill(ralpha/(4*atan(1))*180); i++;}
    }
    // cout << "Prima di scrivere h" << endl;
    h->Write();
    // cout << "Scritto h" << endl;
    delete h;
    //  cout << "Cancellato h" << endl;
  }
  
  out->Close();  
}
