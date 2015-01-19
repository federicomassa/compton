#include <TF1.h>
#include <cmath>
#include <TGraphErrors.h>
using namespace std;

void sigma_ph() {

  double sigma_ph[7] = {17.46,8.06,4.55,2.92,1.51,0.945,0.436};
  double energia[7] = {3,4,5,6,8,10,15};
  double err_sigma[7] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double err_en[7] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  TF1* fitfunc = new TF1("Fitting Function","[0]/(pow(x,[1]))",2,15);
  fitfunc->SetParameters(150,2);
  TGraphErrors* graph = new TGraphErrors(7,energia,sigma_ph,err_en,err_sigma);
  graph->Fit(fitfunc,"0RS");
  graph->DrawClone("APE");
  fitfunc->Draw("SAME");
}
