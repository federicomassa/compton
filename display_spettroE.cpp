#include <TGraphErrors.h>
#include <TAxis.h>
#include <TCanvas.h>

using namespace std;

void display_spettroE() {
  TCanvas *canv = new TCanvas("canv","Canvas",700,700);
  TGraphErrors* graph = new TGraphErrors("./Gate/Egate-50.dat","%lg %lg"," \t");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.6);
  graph->SetMarkerColor(kRed);
  graph->GetXaxis()->SetNdivisions(74);
  graph->GetXaxis()->SetLabelSize(0.015);
  canv->SetGrid();
  graph->Draw("AP");
  
}
