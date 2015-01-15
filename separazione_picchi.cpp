#include <cmath>
#include <TGraphErrors.h>

using namespace std;

double me = 511; //keV
double E1 = 1172;
double E2 = 1332;

double dist(double theta) {
  return (-E1/(1+E1/me*(1-cos(theta/180*4*atan(1))))+E2/(1+E2/me*(1-cos(theta/180*4*atan(1)))));
    }

void separazione() {
  double theta[361];
  double sep[361];
  for (int i = 0; i < 361; i++) {
    theta[i] = -90 + double(i)*0.5;
    sep[i] = dist(theta[i]);
}

  TGraphErrors *graph = new TGraphErrors(361,theta,sep,NULL,NULL);
  graph->SetTitle("Separazione tra i due picchi; Theta(gradi); Separazione (keV)");
  graph->DrawClone("APE");
  
}
