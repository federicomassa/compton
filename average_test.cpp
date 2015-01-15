#include <cmath>
#include <TGraphErrors.h>

using namespace std;

double me = 511; //keV
double E1 = 1172;
double E2 = 1332;

double av_vera(double theta) {
  return 0.5*(E1/(1+E1/me*(1-cos(theta/180*4*atan(1))))+E2/(1+E2/me*(1-cos(theta/180*4*atan(1)))));
    }

double av_falsa(double theta) {
  return 0.5*(E1+E2)/(1+(E1+E2)/(2*me)*(1-cos(theta/180*4*atan(1))));
}

void average() {
  double theta[361];
  double rel[361];
  for (int i = 0; i < 361; i++) {
    theta[i] = -90 + double(i)*0.5;
    rel[i] = (av_vera(theta[i])-av_falsa(theta[i]))/av_vera(theta[i]);
}

  TGraphErrors *graph = new TGraphErrors(361,theta,rel,NULL,NULL);
  graph->SetTitle("Errore relativo; Theta(gradi); Errore relativo");
  //TF1 *func = new TF1("func","Errore relativo;Theta (rad)
  graph->DrawClone("APE");
  
}
