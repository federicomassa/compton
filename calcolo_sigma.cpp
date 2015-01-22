#include <TMath.h>
#include <iostream>
#include <TGraphErrors.h>

using namespace std;


double me = 511; //keV
double r0 = 2.81794E-13; //cm 

double Ef(double Ei, double theta) {
  return Ei/(1+Ei/me*(1-cos(theta*4*atan(1)/180)));
}

double P(double Ei, double theta) {
  return Ef(Ei, theta)/Ei;
}

double KN(double E, double theta) {
  //proporzionale a
  return P(E,theta)*P(E,theta)*(P(E,theta)+1/P(E,theta)-1+cos(theta*4*atan(1)/180)*cos(theta*4*atan(1)/180));
}

double sigma(double E) {
  double g = E/me;
  double c = 8*atan(1)*r0*r0*((1+g)/(g*g)*(2*(1+g)/(1+2*g)-1/g*TMath::Log(1+2*g))+1/(2*g)*TMath::Log(1+2*g)-(1+3*g)/(1+2*g)/(1+2*g)); //fm^2
  return c;
}

void ratio() {
  double E1 = 1172; // keV
  double E2 = 1332; //keV
  double phiratio =  0.8625; //letto da wikipedia - Cobalt-60 = phi2/phi1
  double theta[41];
  double r[41];
  for (int thetai = 20; thetai <= 60; thetai++) {
    theta[thetai-20] = double(thetai);
    r[thetai-20] = KN(E2,double(thetai))/KN(E1,double(thetai))*phiratio;
  }
  TGraphErrors* graph = new TGraphErrors(41,theta,r,NULL,NULL);
  graph->DrawClone("APE");
}

