#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <iostream>
#include <cmath>
#include "TH1F.h"
#include "TRandom3.h"

double me = 510.9989;
double alpha=2.45; // good E_mean/m_e value


double mu(double th){
  
  double T=1252.85*(alpha*(1-cos(3.1415/180*th))/(1+alpha*(1-cos(3.1415/180*th))));
  return T;
}

double E_el(double th) {
  double E = mu(th) + me;
  return E;
}

double beta(double th) {
  double p = pow(E_el(th)*E_el(th)-me*me,0.5);
  double b = p/E_el(th);
  return b;
}

double beta_gamma(double th) {
double p = pow(E_el(th)*E_el(th)-me*me,0.5);
 double bg = p/me;
 return bg;
}

double K=307; // KeVcm2/mol
double SCgrmol=11.065; // gr/mol, grammo/mole dello scintillatore
double I=0.0647; //mean exitation energy, KeV, plastic scintillator
double J=0.2; // parametro misterioso
double density=1.03; // plastic density in g/cm^3
double L=2.5; //scintillator lenght in cm
double avZsuA=0.5425; // Z/A medio dello scintillatore, usando rapporti di peso trovati nel sito

double csi(double th){// 1/4 landau width, even though the concept of width doesn't make sense. also a fundamental parameter in the calculation of the MVP
  double x=density*L;
  double cs=K/(2*SCgrmol)*avZsuA*x/(beta(th)*beta(th));
  return cs;
}



double MVP(double th){ // most probable value of the landau
  double massimo=csi(th)*(TMath::Log(2*me*beta_gamma(th)*beta_gamma(th)/I)+TMath::Log(csi(th)/I)+J-beta(th)*beta(th));
  return massimo;
}


double tau(double th){
 double t=mu(th)/me;
  return t;
}

double F(double th){
  double factor=1-beta(th)*beta(th)+((tau(th)*tau(th))/8-(2*tau(th)+1)*TMath::Log(2))/((tau(th)+1)*(tau(th)+1));
  return factor;
}

double eBethe(double th){ //bethe formula for electrons
  double loss=K/2*density*avZsuA/(beta(th)*beta(th))*(TMath::Log(tau(th)*tau(th)*(tau(th)+2)/((2*I/me)*(2*I/me)))+F(th));//neglecting density and shell corrections
  return loss;
}

/*double bethe(double th){
double betebl=
}*/

double attenlenght=43;//attenuation lenght in cm
double transeps=0.5; // efficienza di trasporto luce da scint a PMT, mettiamo 0.5 tanto per
double QE=0.25; // efficienza quantistica del fotocatodo, 25% è valore molto comune e probabile
double photXkev=10; // how many photons are emitted per kev of energy loss in the detector
double PMTamp=531441;//3^12
double TotalRef=0.20;//frazione di fotoni che rimangono nello scintillatore per rif tot
double avPhotoElectrons(double th){
  double elec;
  if(eBethe(th)<mu(th)){
    elec=eBethe(th)*photXkev*transeps*QE*density*L*TotalRef*exp(-15/attenlenght);
  }
  if(eBethe(th)>mu(th)){
    elec=mu(th)*photXkev*transeps*QE*density*L*TotalRef*exp(-15/attenlenght);
  }
  return elec;
}

double resistence=50;//standard output resistence, ohm
double plastictime=10E-9;// durata del segnale del plastico, ns
double echarge=1.6E-19;//electron charge in Coulomb

double amplV(double th){
 double amp=avPhotoElectrons(th)*PMTamp*echarge*resistence/plastictime;
  return amp;
}
void draw(){
  TF1 *dis= new TF1("dis","1252*(2.446*(1-cos(3.1415/180*x))/(1+2.446*(1-cos(3.1415/180*x))))", -90,90);
  //dis->Draw();
}



void graphk(){

  TF1* kgraph=new TF1("grafico_k", "amplV(x)", 20, 70);
  kgraph->SetTitle("Recon. average signal amplitude; theta (gradi);Amplitude (V)");

  kgraph->Draw();

}

//TCanvas* c=new TCanvas("c","C",1000,500);
//c->Divide(2,1);

double variance(double th){

  double npoints=1000000;
  TH1F *histo =new TH1F("hist","hist",160,200,1800);
  TH1F *histo2= new TH1F("hist2","hist2",1000,22000,1220000);
  TRandom3 rndgen;
  double s_delta;
  double s_elec;
  double s_firststad, s_firststad2;

  for (int i=0;i<npoints;i++){
    if (i%100000==0) {
      cout << i << endl;
    }

    s_delta = rndgen.Poisson(3);
    // cout << s_delta;
    s_elec = rndgen.Gaus(avPhotoElectrons(th),sqrt(avPhotoElectrons(th)));
    s_firststad = s_delta*s_elec;
    s_firststad2 = s_delta*s_delta*s_elec*s_elec;
  
    histo->Fill(s_firststad);
    histo2->Fill(s_firststad2);
  }
  //c->cd(1);
histo->Draw();
//  c->cd(2);
histo2->Draw();
  double Ex2=histo2->GetMean();
  double Ex=histo->GetMean();
  double var=Ex2-Ex*Ex;
  return var;
    
    }


//double variance(double th){
//double npoints=1000000;
//TH1F *histo =new TH1F("hist","hist",20,200,300);//160,200,1800);
//TH1F *histo2= new TH1F("hist2","hist2",1000,22000,1220000);
  
//for (int i=0;i<npoints;i++){
//  if (i%100000==0) {
//    cout << i << endl;
//  }
//  TRandom3 rndgen;
//  double s_delta = rndgen.Poisson(3);
    // cout << s_delta << endl;
//  double s_elec = rndgen.Gaus(avPhotoElectrons(th),sqrt(avPhotoElectrons(th)));
//   double s_firststad = s_delta*s_elec;
//  double s_firststad2 = s_delta*s_delta*s_elec*s_elec;
    /* if (1%10000==0){
    cout << s_firststad << endl;
    }*/
//  histo->Fill(s_elec);
//  histo2->Fill(s_firststad2);
//}
  //c->cd(1);
//histo->Draw();
//  c->cd(2);
//histo2->Draw();
//double Ex2=histo2->GetMean();
//double Ex=histo->GetMean();
//double var=Ex2-Ex*Ex;
//return var;
    
//   }

void gaussian(){
  //genero i fotoelettroni con gaussiana a media quella che è e varianza idem, poi genero poissoniananamente la delta, faccio il prodotto, metto in un histogramma e faccio get rms

}










void KNweight() {
     
    double mvpangolo[100];
    double mvpgraph[100];
    for (int i = 1; i < 101; i++) {
      mvpangolo[i-1] = 20+double(i)/100*70;
      mvpgraph[i-1] = MVP(mvpangolo[i-1]);
      cout << mvpangolo[i-1] << "..." << mvpgraph[i-1] << endl;
    }
    TGraphErrors* mvp = new TGraphErrors(100,mvpangolo,mvpgraph,NULL,NULL);
    mvp->DrawClone("APE");

  //  double sigma[10];
  double soglia[100];
  double coeff[10][10];
  double vecttheta[5] = {-62,-57,-47,53,43};
  double theta;
  int scelta;
        double coeff_land[100];
	double coeff_ang[5];

	//  TCanvas* c=new TCanvas("c","Grafici", 900,600);
	//  c->Divide(3,2);
	// TCanvas* c2=new TCanvas("c2", "Steffi Graf", 1800,1200);
  //  c2->Divide(3,2);


for (int ang=0;ang<5;ang++){
      theta=vecttheta[ang];

      for (int soglia_i=0;soglia_i<101;soglia_i++){
	soglia[soglia_i]=205+soglia_i*(600-250)/100;


       TF1* land = new TF1("Landau Curve","TMath::Landau(x,[0]+0.22,[1],1)",-10,10000);
       land->SetParameters(MVP(theta),4*csi(theta));
       //       land->Draw();
       land->SetNpx(100000);
       coeff_land[soglia_i]=land->Integral(soglia[soglia_i],10000);

       if (soglia[soglia_i]==530){
	 coeff_ang[ang]=coeff_land[soglia_i]; // riempio vettore di coefficienti per vari angoli ad una determinata soglia per poi fare grafico
     
       cout << coeff_land[soglia_i] << endl;

       }
      }//for soglia
  

 
      TGraphErrors *cff_soglia = new TGraphErrors(100, soglia ,coeff_land ,NULL,NULL);

      switch (ang){
      case 0:
	//	c->cd(1);
	//	cff_soglia->DrawClone("APE");
	break;
      case 1:
	//	c->cd(2);
	//	cff_soglia->DrawClone("APE");
	break;  
    case 2:
      //	c->cd(3);
	//	cff_soglia->DrawClone("APE");
	break;
      case 3:
	//	c->cd(4);
	//	cff_soglia->DrawClone("APE");
	break; 
     case 4:
       //	c->cd(5);
	//	cff_soglia->DrawClone("APE");
	break;
      }
 }//giusto
  
  TGraphErrors *con = new TGraphErrors(5,vecttheta, coeff_ang, NULL, NULL);
  // c2->cd();con->Draw("APE");
}


