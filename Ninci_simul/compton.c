#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <iostream>
#include <time.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <vector>

#include <TROOT.h>
#include <TApplication.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <fstream>
#include <string>
#include "TVector2.h"
#include "TVector3.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TEllipse.h>
#include "TObject.h"
#include "TLegend.h"
#include <TTree.h>

using namespace std;

TH1F *hen= new TH1F("hen","Simulazione Montecarlo: Spettro NaI",2000,0,2000);
TH3F *hap = new TH3F("hap", "punto di assorbimento", 100, -3, 3, 100, -3, 3, 100, -1, 6);
TH1F *hapz= new TH1F("hapz","punto z assorbimento",100,-1,6);

Double_t rgauss();
Double_t random1();
Double_t Sigma_Comp(Double_t gamma);
Double_t Ang_gen(Double_t gamma);
Double_t kleinnishina(Double_t x, Double_t par);
Double_t l_compton(Double_t gamma);
Double_t l_pe(Double_t gamma);
int efficienza(Double_t theta_fot, Double_t en_originale, Double_t energia, int nevent, Double_t smearingenergia);
Double_t energia_persa(Double_t s_max, Double_t energia);

void ninci()
{ 
    int nevent=1000000000;
    Double_t energia;
    Double_t angoli[10]={25,30,35,40,45,50,55,60,80,0}; //si cicla sugli angoli da 25 a 60 e 80
    Double_t energiebasse[10]={963.4, 895.4, 827.4, 761.9, 700.4, 643.6, 592, 545.5, 404.6,1330}; // energie corrispondenti al picco di 1,17
    Double_t energiealte[10]={1069.3, 986.1, 904.3, 826.7, 754.7, 689.2, 630.4, 577.9, 422.1,1330}; // energie corrispondenti al picco di 1,33
    Double_t smearingenergie[10]={65, 68, 68, 66, 62, 57, 52, 47, 29, 0};
    int successi[10];
    int j;

    cout << "codice angolo (0-9): " << flush;
    cin >> j;
    cout << angoli[j] << " gradi" << endl;


//    for(int j=0;j<9;j++) 
 {
  hen->Reset();
  successi[j]=efficienza(angoli[j], 1170., energiebasse[j], nevent,smearingenergie[j]);
  successi[j]+=efficienza(angoli[j], 1330., energiealte[j], nevent,smearingenergie[j]);

  //cout << "energia in keV: " << flush;
  //cin >> energia;

    cout << angoli[j] << " gradi " <<  (Double_t)successi[j]/(2*nevent)*100 << " %"<<endl;

 TCanvas c2("c2","c2");
    c2.cd();
    hen->GetXaxis()->SetTitle("Energia (keV)");
    hen->GetXaxis()->CenterTitle();
    hen->GetYaxis()->SetTitle("Conteggi");
    hen->GetYaxis()->CenterTitle();
    hen->Draw();
    c2.Update();
    
    switch (j){
         case 0:
	   c2.Print("hen25.png");
           c2.Print("hen25.root");
	 break;
         case 1:
	   c2.Print("hen30.png");
           c2.Print("hen30.root");
	 break;
         case 2:
	   c2.Print("hen35.png");
           c2.Print("hen35.root");
	 break;
         case 3:
	   c2.Print("hen40.png");
           c2.Print("hen40.root");
	 break;
         case 4:
	   c2.Print("hen45.png");
           c2.Print("hen45.root");
	 break;
         case 5:
	   c2.Print("hen50.png");
           c2.Print("hen50.root");
	 break;
         case 6:
	   c2.Print("hen55.png");
           c2.Print("hen55.root");
	 break;
         case 7:
	   c2.Print("hen60.png");
           c2.Print("hen60.root");
	 break;
         case 8:
	   c2.Print("hen80.png");
           c2.Print("hen80.root");
	 break;
         case 9:
	   c2.Print("hen0.png");
           c2.Print("hen0.root");
	 break;}

    char* out_file_root = new char[256];
    out_file_root = (char*)"compton.root";
    TFile foutput(out_file_root,"recreate");
    hen->Write();
    hap->Write();
    hapz->Write();
    foutput.Close(); 
  } 
}	

/*****************************************************
		fine main
****************************************************/

int efficienza(Double_t theta_fot, Double_t en_originale, Double_t energia, int nevent, Double_t smearingenergia){

  const Double_t PI=TMath::Pi();
  long seed=time(NULL);
  srand(seed);

  Double_t r = 5.09/2;
  Double_t rho, theta0, x0, y0, z_plastic, p_el, p_el_trasverso, range_max, lost_energy;
  TVector3 posizione, direzione, assex, asserotazione, proiezione;
  Double_t zeta_fot_pe, zeta_fot_compton;
  Double_t ang_comp,ang_phi,ang_el;
  

  Double_t n1g=3.667;//densità NaI (g/cm^3)

  Double_t n1=n1g*TMath::Na()*.42697;//num e- in NaI

  Double_t en, en_compton, en_compton_provv, gamma, sigma, lambda_nai, lambda_compton, lambda_pe;
  Double_t me=511;
  int i;
  int successi=0;
// parti relative all'efficienza del plastico
      p_el_trasverso=energia*sin(theta_fot*TMath::DegToRad());
      p_el=sqrt(pow(en_originale-energia+511,2)-511*511);
      ang_el=TMath::ASin(p_el_trasverso/p_el);//in radianti
// fine efficienza plastico

for(i=1;i<=nevent;i++)
    {
//      if (i%(nevent/10000)==0) cout << "\r" << (Double_t)i/nevent << "\t" << (Double_t)successi/i<<"  " << flush;
      rho=sqrt(random1())*r;//centro cerchiolino generatore
      theta0=random1()*2*PI;
      x0=rho*cos(theta0);
      y0=rho*sin(theta0);

      posizione.SetXYZ(x0,y0,0);
      direzione.SetXYZ(0,0,1);
      assex.SetXYZ(1,0,0);

      en=energia+smearingenergia*rgauss();
      en_compton=0;
      gamma=en/me;

// parti relative all'efficienza del plastico
      z_plastic=2*random1();
      range_max=z_plastic/cos(ang_el);
      lost_energy=energia_persa(range_max,en_originale-energia);
      lost_energy=lost_energy+rgauss()*sqrt(lost_energy*0.1);//risoluzione data dal radq(Nfotoni) prodotti (1fotone=0.1keV)
 //      cout << "energia el. "<< en_originale-energia <<" en. persa "<<lost_energy << endl;

      if(lost_energy>360)
// fine efficienza plastico

      while(posizione.Z() < (2*r) && posizione.Z() >= 0 && posizione.X()*posizione.X()+posizione.Y()*posizione.Y() < r*r) 
	{
	  lambda_compton=l_compton(gamma);
	  lambda_pe=l_pe(gamma);

	  zeta_fot_compton= - lambda_compton * log(random1());//spazio percorso dal fotone prima di interagire compton
	  zeta_fot_pe= - lambda_pe * log(random1());//spazio percorso dal fotone prima di interagire fotoelettrico

	  if (zeta_fot_compton<zeta_fot_pe) {

	    if ((posizione+zeta_fot_compton*direzione).Z()<(2*r) && (posizione+zeta_fot_compton*direzione).Z()>=0 && (posizione+zeta_fot_compton*direzione).X()*(posizione+zeta_fot_compton*direzione).X()+(posizione+zeta_fot_compton*direzione).Y()*(posizione+zeta_fot_compton*direzione).Y() < r*r) 
	      {
		posizione = posizione+zeta_fot_compton*direzione;
		ang_comp= Ang_gen(gamma);//genero angolo 
		ang_phi=random1()*2*PI;
		asserotazione=direzione.Cross(assex);
		asserotazione.Rotate(ang_phi,direzione);
		direzione.Rotate(ang_comp,asserotazione);
		en_compton_provv=en;
		en=en/(1.+(en/me)*(1.-cos(ang_comp)));
		en_compton_provv-=en;
		en_compton+=en_compton_provv;
		gamma=en/me;
		sigma=Sigma_Comp(gamma);
	      } 

	    else 
	      {
		posizione.SetXYZ(0,0,1000);
		if(en_compton>300) hen->Fill(en_compton+1000*(0.026*pow(en_compton/1000,0.5))*rgauss());
	      } 
	  }

	  else{

	    if((posizione+zeta_fot_pe*direzione).Z()<(2*r) && (posizione+zeta_fot_pe*direzione).Z()>=0 && (posizione+zeta_fot_pe*direzione).X()*(posizione+zeta_fot_pe*direzione).X()+(posizione+zeta_fot_pe*direzione).Y()*(posizione+zeta_fot_pe*direzione).Y() < r*r)
	        {
		    hapz->Fill((posizione+zeta_fot_pe*direzione).Z());
		    hap->Fill((posizione+zeta_fot_pe*direzione).X(),(posizione+zeta_fot_pe*direzione).Y(),(posizione+zeta_fot_pe*direzione).Z());
	            posizione.SetXYZ(0,0,1000);
		    successi++; 
		    if(en+en_compton>300) hen->Fill(en+en_compton+1000*(0.026*pow((en+en_compton)/1000,0.5))*rgauss());
	        }
	    else 
	      {
		posizione.SetXYZ(0,0,1000);
		if(en_compton>300) hen->Fill(en_compton+1000*(0.026*pow(en_compton/1000,0.5))*rgauss()); //bisogna salvare anche questo
	      } 
	  }
    	}
    }

return successi;
}

/**************************************************************************/

Double_t l_compton(Double_t gamma) {
  gamma*=511;
  if(gamma<5.188) cout << "energia bassa!" << endl;
  Double_t energie[25] = {5.188,6,8,10,15,20,30,33.16,33.17,40,50,60,80,100,150,200,300,400,500,600,800,1000,1022,1250,1500};
  Double_t lambda[25] = {6.56,5.83,4.70,4.05,3.24,2.88,2.55,2.50,2.50,2.41,2.35,2.33,2.35,2.39,2.55,2.73,3.08,3.40,3.71,4.00,4.54,5.04,5.09,5.62,6.18};
  int pointer=0;
  for (int i=0;i<25;i++) { if(gamma>=energie[i]) pointer=i; } //mettendo i+1 si prende il successivo...

  Double_t ret=lambda[pointer]+(gamma-energie[pointer])*(lambda[pointer+1]-lambda[pointer])/(energie[pointer+1]-energie[pointer]);

return ret;
}

Double_t l_pe(Double_t gamma) { //il primo sarebbe 0.36
  gamma*=511;
  Double_t energie[25] = {5.188,6,8,10,15,20,30,33.16,33.17,40,50,60,80,100,150,200,300,400,500,600,800,1000,1022,1250,1500};
  Double_t lambda[25] = {0,0.51,1.10,1.98,5.97,13.24,41.26,54.87,9.15,14.98,27.00,44.63,99.16,186.78,594.12,1350.01,4208.37,9120.48,16232.29,25250.23,48783.99,78138.25,82139.30,121742.18,170439.05};
  int pointer=0;
  for (int i=0;i<25;i++) { if(gamma>=energie[i]) pointer=i; } 

  Double_t ret=lambda[pointer]+(gamma-energie[pointer])*(lambda[pointer+1]-lambda[pointer])/(energie[pointer+1]-energie[pointer]);
  return ret/1000;

}

Double_t energia_persa(Double_t s_max, Double_t energia)
{
  Double_t s=0, e=0, delta_e, delta_x;
  Double_t energie[34] = {0.01,0.0125,0.015,0.0175,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1,1.25};
  Double_t dedx[34] = {22.65,19.03,16.5,14.63,13.19,11.09,9.646,8.581,7.764,7.115,6.588,6.15,5.78,5.191,4.741,4.386,4.099,3.576,3.222,2.969,2.779,2.514,2.342,2.222,2.134,2.066,2.014,1.973,1.941,1.893,1.861,1.839,1.824,1.804};
  int pointer=0,i;
  for (i=0;i<34;i++) { if(energia>=energie[i]*1000) pointer=i+1; } //prendo il più alto, sono pessimista

  while (s<s_max && pointer>0) 
   {
    delta_e = energie[pointer]-energie[pointer-1];
    delta_x = delta_e / dedx[pointer];
    e+=delta_e;
    s+=delta_x;
    pointer-=1;
   }

if(pointer==1) { return energia; cout << "morto!" << endl; }
else return e*1000;

}

Double_t Sigma_Comp(Double_t gamma){

  Double_t PI=TMath::Pi();
  Double_t re=2.8179402e-13;//raggio e- in cm
  Double_t sigma_comp;

  sigma_comp=2*PI*re*re*((1.+gamma)/(gamma*gamma)*((2.*(1.+gamma))/(1.+2.*gamma)-log(1.+2.*gamma)/gamma )+log(1.+2.*gamma)/(2.*gamma)-(1.+3.*gamma)/(pow(1.+2.*gamma,2)));

  return sigma_comp;
}

Double_t rgauss() //genera un numero distribuito gaussiamente
{
  Double_t phi = 2*TMath::Pi()*random1();
  Double_t r = -2*log(random1());
  return sqrt(r)*sin(phi);
}

Double_t random1() //genera un numero random tra 0 e 1
{
  return (rand()%10000001)/10000000.;
}

Double_t Ang_gen(Double_t gamma){
  Double_t ang_gen = random1() * TMath::Pi() * 2;
  Double_t anggeny = random1();

  while(anggeny > kleinnishina(ang_gen,gamma))
    {
      ang_gen = random1() * TMath::Pi() * 2;
      anggeny = random1();
    }
  return ang_gen;
}

Double_t kleinnishina(Double_t x, Double_t gamma) { //normalizzato a 1...
  return

    (pow((1+gamma*(1-cos(x))),-2)
     *(
       (1/(1+gamma*(1-cos(x))))
       +(1+gamma*(1-cos(x)))
       -1+pow(cos(x),2)
       )*0.5);

}


