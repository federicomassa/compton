#include <string>
#include <cmath>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1.h>


void ultimacalibrazione(){
     
  TCanvas *c1 = new TCanvas("c1","Calibrazione",700,700);
      
     string sources[7] = {"Am1","Am2","Co1","Co2","Na1","Na2","Cs"}; 

   double chval[7] = {459.7,921.25,6291.42,6574.48,4315.13,6432.36,4991.01};
   double enval[7] = {18, 60, 1173.2, 1332.5,511.0034,1274.5,661.64}; //keV, dubbio su Am1, 14 keV?
   double cherr[7] = {81.65605,61.1465,63.43949,60.38217,92.03397,59.40552,82.03822}; //FWHM/2
   TGraphErrors *graph = new TGraphErrors(7,chval,enval,cherr,NULL);
   
   
  // TF1* noam = new TF1("No_Am", "474.629+0.029255*x-3.66635E-5*x*x+6.70617E-9*x*x*x",0,8192);
//   noam->SetLineStyle(2);
//   noam->SetLineColor(kOrange);
//   noam->DrawClone("SAME");
    
     graph->SetTitle("Calibrazione;Canale;Energia(keV)");
     graph->GetYaxis()->SetTitleOffset(1.4);
     graph->SetLineColor(kRed);
     graph->SetMarkerColor(kRed);
     graph->SetMarkerStyle(20);
     graph->SetMarkerSize(0.8);
     graph->DrawClone("APE"); 
     
     
     TF1 *fitfunc = new TF1("Fitting Function", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0,8192);
// TF1 *fitfunc = new TF1("Fitting Function", "[0]+[1]*x+[2]*x*x+4.03156E-9*x*x*x", 0,8192);
//     TF1 *fitfunc = new TF1("Cubic Fit", "-33.0933+[0]*x-1.80489E-5*x*x+4.03119E-9*x*x*x", 0,8192);
     fitfunc->SetParameters(-30,0.1,-1E-5,4E-9);
//     fitfunc->SetParameter(0,0.1);
     TFitResultPtr frp = graph->Fit(fitfunc,"0SR");
     cout << "Cubic Fit Details: " << endl;
     frp->Print();
     TMatrixDSym covMatrix(frp->GetCovarianceMatrix());
     covMatrix.Print();
     fitfunc->SetLineColor(kBlue);
     fitfunc->DrawClone("SAME"); 
     
     
     //Slope = 0.347639
     TF1 *linfit = new TF1("Linear Fit", "[0]+[1]*x", 4500,8192);
  //   TF1 *linfit = new TF1("Parabolic Fit", "[0]+[1]*x+[2]*x*x", 4500,7500);
     linfit->SetParameters(-800,0.3);
//     linfit->SetParameters(-100,+0.1,1E-5);
     linfit->SetLineStyle(1);
     TFitResultPtr lfrp = graph->Fit(linfit,"0SR+E");
     cout << "Linear Fit Details: " << endl;
//cout << "Parabolic Fit Details: " << endl;
     lfrp->Print();
     TMatrixDSym lincovMatrix(lfrp->GetCovarianceMatrix());
     lincovMatrix.Print();
     linfit->SetLineColor(3);
    // linfit->DrawClone("SAME");
   
     
     
  ////Cubica: param errors
     double eparam3[4];
     for (int i = 0; i <= 3; i++) {
         
       eparam3[i] = double(pow(frp->CovMatrix(i,i),0.5));}
//     Usando gli errori di tutti i parametri senza covarianza
//     for (int i = 0; i <= 3; i++) {
//     fitfunc->SetParameter(i,fitfunc->GetParameter(i)+ (eparam3[i]));}
////     fitfunc->SetParameter(i,fitfunc->GetParameter(i)+1.65E-3);}
//     fitfunc->SetLineColor(kOrange);
////     fitfunc->DrawClone("SAME");
//     
//     for (int i = 0; i <= 3; i++) {
//     fitfunc->SetParameter(i,fitfunc->GetParameter(i)- 2*(eparam3[i]));}
////     fitfunc->SetParameter(i,fitfunc->GetParameter(i)-2*1.65E-3);}
//     fitfunc->SetLineColor(6);
////     fitfunc->DrawClone("SAME");

//Calcolo errori cubica con covarianza, tentativo

     

     

//
// //  Linfit: param errors
//     double *eparam = new double;
//     for (int i = 0; i <= 1; i++) {
//         *(eparam+i) = pow(lfrp->CovMatrix(i,i),0.5);}


                      
         
     //for (int i = 0; i <= 1; i++) {
//     linfit->SetParameter(i,linfit->GetParameter(i)+ *(eparam+i));}
//     linfit->SetLineColor(3);
//     linfit->SetLineStyle(2);
//     linfit->DrawClone("SAME");
//     
//     for (int i = 0; i <= 1; i++) {
//     linfit->SetParameter(i,linfit->GetParameter(i)-2*(*(eparam+i)));}
//     linfit->SetLineColor(3);
//     linfit->SetLineStyle(2);
//     linfit->DrawClone("SAME"); 

     
   //  TF1* tom = new TF1("Tommaso", "1.39364E-5*x*x+0.0868444*x+10.4319",0,8192);
   //tom->DrawClone("SAME");
     
    // TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
//     leg->SetHeader("Legend Title");
//     leg->AddEntry(graph, "Prova graph", "le");
//     leg->AddEntry(linfit, "Prova linfit", "l");
//     leg->AddEntry(fitfunc, "Prova fitfunc", "l");
//     leg->SetFillColor(0);
//     leg->Draw();
//     
    

     
     c1->SetGrid();
     TF1 *offset = new TF1("offset","27+0.103793*x-1.80489E-5*x*x+4.03119E-9*x*x*x",0,8192);
     offset->SetLineColor(2);
    // offset->Draw("SAME");
     c1->Update();
     }             
     
     
