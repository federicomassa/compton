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

Double_t rgauss();
Double_t random1();
Double_t Sigma_Comp(Double_t gamma);
Double_t Ang_gen(Double_t gamma);
Double_t kleinnishina(Double_t x, Double_t par);
TVector2 intersezione (TVector2 a1, TVector2 a2, TVector2 b1, TVector2 b2);
TVector2 propagazionecompton(TVector2 posizione, TVector2 direzione, Double_t rpaz, Double_t raggio, Double_t nriv, Double_t taglio_en);
TVector2 propagazionesenzacompton(TVector2 posizione, TVector2 direzione, Double_t raggio, Double_t nriv);
//TVector2 propagazionecomptonNOdiscr(TVector2 posizione, TVector2 direzione, Double_t rpaz, Double_t raggio);

/***************************************************************
         DICHIARAZIONE DEGLI HISTOGRAMMI 
****************************************************************/
int n_bin=5000;
TH1F *hintx= new TH1F("hintx","Proziezioni asse x",n_bin,-80,80);
TH1F *hinty= new TH1F("hinty","Proziezioni asse y",n_bin,-80,80);
TH2F *hintxy= new TH2F("hintxy","Intersezioni nel piano con smearing dei fotoni",n_bin,-80,80,n_bin,-80,80);
TH1F *hintx_gen= new TH1F("hintx_gen","Risoluzione asse x (generati VS discreti)",n_bin,-80,80);
TH1F *hinty_gen= new TH1F("hinty_gen","Risoluzione asse y (generati VS discreti)",n_bin,-80,80);
TH2F *hintxy_gen= new TH2F("hintxy_gen","Intersezioni nel piano generate",n_bin,-80,80,n_bin,-80,80);
TH2F *hpuntiriv= new TH2F("hpuntiriv","Punti colpiti sul rivelatore",n_bin,-80,80,n_bin,-80,80);
TH1F *hintx_discr= new TH1F("hintx_discr","Proziezioni asse x discretizzati",n_bin,-80,80);
TH1F *hinty_discr= new TH1F("hinty_discr","Proziezioni asse y discretizzati",n_bin,-80,80);
TH2F *hintxy_discr= new TH2F("hintxy_discr","Intersezioni nel piano discreto con smearing dei fotoni",n_bin,-80,80,n_bin,-80,80);

typedef struct{
        TVector2 x[4];
}  TVector2s;

Int_t main(Int_t argc, char** argv)
{
  const Double_t PI=TMath::Pi();
  int nevent=100000;
  long seed=time(NULL);
  srand(seed);
  int forma;
  Double_t taglio_en=400;
  Double_t raggio=80; //raggio rivelatore in cm
  Double_t rpaz=35; //raggio paziente in cm
  Double_t d=1;//dimensione tumore in cm
  Double_t nriv=floor(2*PI*raggio/0.4);
  cout << "Scelgliere la forma dell'organo tumorale:"<<endl;
  cout << "1 cerchio"<<endl;
  cout << "2 quadrato"<<endl;
  cout << "3 triangolo rettangolo"<<endl;
  cout << "4 triangolo equilatero"<<endl;
  cout << "5 corona circolare"<<endl;
  cout << "6 cuoricino"<<endl;
  cin >> forma;
  switch (forma){
  case 1:
    cout<<"è stato scelto il cerchio"<< endl;
    break; 
  case 2:
    cout<<"è stato scelto il quadrato"<< endl;
    break;
  case 3:
    cout<<"è stato scelto il triangolo rettangolo"<< endl;   
    break;
  case 4:
    cout<<"è stato scelto il triangolo equilatero"<< endl;   
    break;
  case 5:
    cout<<"è stato scelto la corona circolare"<< endl;   
    break;
  case 6:
    cout<<"è stato scelto il cuoricino"<< endl;
    break;
  case 7:
    cout<<"simulazione..." << endl;
    break;
  }
  Double_t rho=sqrt(random1())*rpaz*.9;//centro cerchiolino generatore
  Double_t theta0=random1()*2*PI;
  Double_t x0=rho*cos(theta0);
  Double_t y0=rho*sin(theta0);
    
  /***************************************************************
        generazione degli eventi
  ****************************************************************/
  //TVector2 coppie[nevent][2];

  TVector2 coppie_discr[nevent][2];
  Double_t x,y;
  Double_t prov;
  Double_t r,theta,tria1,tria2; 
  Double_t zeta1=10.;//num e- in acqua
  Double_t zeta2=14.;//num e- in aria supposta di azoto N2
  Double_t n1g=1.;//densità acqua (g/cm^3)
  Double_t n2g=1.205e-3;//densità aria (g/cm^3)
  Double_t nav=TMath::Na();
  Double_t n1=n1g*nav*zeta1/20.;//num e- in acqua
  Double_t n2=n2g*nav*zeta2/28.;//num e- in aria
  Double_t en=511;//energia fotone keV
  Double_t me=511;//massa e- keV
  Double_t gamma=en/me; //energia fotone/massa c^2
  int ncompton=0, npuntiint=0, successi=0;
  
  //dichiarazioni varie
  Double_t quad1, quad2;
  Double_t phi1,phi2;

  TVector2 decpoint, direzione1, direzione2,backup1, backup2, backup1noC, backup2noC;
  int i_ok=0;
  for(int i=0;i<nevent;i++)
    {
      switch (forma){
      case 1:
	{//cout<<i<<endl;
	  r=sqrt(random1())*d*0.5;//generazione punto partenza fotoni
	  theta=random1()*2*PI;
	  x=r*cos(theta)+x0;//coordinate cartesiane punto partenza
	  y=r*sin(theta)+y0;
	}break;
      case 2:
	{
	  quad1=random1()*d;
	  quad2=random1()*d;
	  x=quad1+x0;
	  y=quad2+y0;
	}break;
      case 3://triangolo rettangolo
	{
	  tria1=random1()*d;
	  tria2=random1()*d;
	  if (tria2>tria1)
	    {
	      prov=tria1; //variabile provvisoria
	      tria1=tria2;
	      tria2=prov;
	    }
	  x=tria1+x0;
	  y=tria2+y0;
	}break;
      case 4://triangolo equilatero
	{
	  tria1=random1();
	  tria2=random1()*sqrt(3)/2;
	  if (tria1<0.5 && tria2>tria1*tan(PI/3))
	    {
	      tria1=0.5-tria1;
	      tria2=1-tria2;
	    }
	  if (tria1>0.5 && tria2>(1-tria1)*tan(PI/3))
	    {
	      tria1=1.5-tria1;
	      tria2=1-tria2;
	    }
	  x=tria1*d+x0;
	  y=tria2*d+y0;
	}break;
      case 5://corona
	{//cout<<i<<endl;
	  r=sqrt(0.15+0.7*random1())*d*0.5;//generazione punto partenza fotoni 
	  theta=random1()*2*PI;
	  x=r*cos(theta)+x0;//coordinate cartesiane punto partenza
	  y=r*sin(theta)+y0; 
	}break;
      case 6://cuore
	{
	  r=sqrt(random1());//cerchio circoscritto raggio unitario
	  theta=random1()*2*PI;
	  while (r>0.3*pow(abs(tan(theta)),(1./abs(tan(theta)) ) )) 
	    {
	      r=sqrt(random1());//cerchio circoscritto
	      theta=random1()*2*PI;
	    }
	  x=r*2.*d*cos(theta)+x0;//coordinate cartesiane punto partenza
	  y=abs(r*2.*d*sin(theta))+y0;
	}break;
      }
//ora (x,y) è il punto malato

      Double_t rho1=sqrt(random1())*.2;//i fotoni si annichilano entro un cerchio di raggio 2mm
      Double_t theta1=random1()*2*PI;
      Double_t x1=rho1*cos(theta1);
      Double_t y1=rho1*sin(theta1);
      
      x=x+x1;//ora (x,y) è il punto di annichilazione
      y=y+y1;


      hintx_gen->Fill(x);
      hinty_gen->Fill(y);  
      hintxy_gen->Fill(x, y);
      
      phi1=random1()*PI;//direzione raggi gamma tra 0 e pigreco
      phi2=PI+phi1+rgauss()*TMath::DegToRad();//sporco il fotone
      
      decpoint.Set(x,y); //punto di partenza fotoni
      direzione1.Set(cos(phi1),sin(phi1));
      direzione2.Set(cos(phi2),sin(phi2));


      backup1noC=propagazionesenzacompton(decpoint,direzione1,raggio,nriv);//ha il problema che non dà Mod=80!!!!!!!!!!!
      backup2noC=propagazionesenzacompton(decpoint,direzione2,raggio,nriv);

      backup1=propagazionecompton(decpoint,direzione1,rpaz,raggio,nriv,taglio_en);
      backup2=propagazionecompton(decpoint,direzione2,rpaz,raggio,nriv,taglio_en);

//cout<<backup1noC.Mod()<<" "<<backup2noC.Mod()<<endl;

      if (backup1.Mod()>rpaz && backup2.Mod()>rpaz) 
        {
         coppie_discr[i_ok][0]=backup1; //i_ok è il sottoindice che riempie solo con i fotoni che superano taglio in energia
         coppie_discr[i_ok][1]=backup2;

         hpuntiriv->Fill(coppie_discr[i_ok][0].X(),coppie_discr[i_ok][0].Y());//punti discreti sul rivelatore
         hpuntiriv->Fill(coppie_discr[i_ok][1].X(),coppie_discr[i_ok][1].Y());
         i_ok++;
        }

    }//fine ciclo for nevent 

  cout<<"Soglia in energia usata: "<<taglio_en<<" keV"<<endl;
  cout << "Percentuale coppie di fotoni usati: "<< float(i_ok)*100/nevent << "%" << endl;

  TVector3 dir_retta1, dir_retta2;
  TVector3 dir_retta1_discr, dir_retta2_discr;
  /* string str;
     int lunghezzabarra = 25;
     for(int kk=0;kk<lunghezzabarra;kk++) { str.append("."); }
     int t=0;  */
  
  // altre dichiarazioni per il ciclo seguente

  TVector2 punto1, punto2, punto1_discr, punto2_discr, punto3, punto4, punto3_discr, punto4_discr, ptint_discr;

  for(int j=0;j<i_ok-1;j++)
    {
	  punto1_discr=coppie_discr[j][0];
	  punto2_discr=coppie_discr[j][1];
        
	  /*	  if (j%(nevent/100)==0) {   
		  if (j%(nevent/lunghezzabarra)==0) { str.replace(t,1,"#"); t++; }
		  cout << "\r[" << str << "] " << ((j*100)/nevent)+1 << "%" << flush;
		  }*/
      
	  for (int k=j+1;k<i_ok;k++)
	    {

		  punto3_discr=coppie_discr[k][0];
		  punto4_discr=coppie_discr[k][1];
	  
		  dir_retta1_discr.SetXYZ((punto2_discr-punto1_discr).X(),(punto2_discr-punto1_discr).Y(),0);
		  dir_retta2_discr.SetXYZ((punto4_discr-punto3_discr).X(),(punto4_discr-punto3_discr).Y(),0);
	  
		  /*utile per il futuro (possiamo aumentare le dimensioni degli scintillatori)  
		    if (dir_retta1.Angle(dir_retta2)>1.49 && dir_retta1.Angle(dir_retta2)<1.64) 
		    {
		    TVector2 ptint=intersezione(punto1,punto2,punto3,punto4);

		    hintx->Fill(ptint.X());
		    hinty->Fill(ptint.Y());
		    hintxy->Fill(ptint.X(),ptint.Y());
		    }*/

		  if (dir_retta1_discr.Angle(dir_retta2_discr)>(PI/2.-.1) && dir_retta1_discr.Angle(dir_retta2_discr)<(PI/2.+.1)) 
		    { npuntiint++;
		      ptint_discr=intersezione(punto1_discr,punto2_discr,punto3_discr,punto4_discr);
		      hintx_discr->Fill(ptint_discr.X());
		      hinty_discr->Fill(ptint_discr.Y());
		      hintxy_discr->Fill(ptint_discr.X(),ptint_discr.Y());
		    }   
	    }
    }

  cout << "Scrittura istogrammi in corso..." << endl;

  /***************************************************************
        istogramma 2D scala a colori discreto
  ****************************************************************/
  Double_t w =1000;
  Double_t h = 1000;
  TCanvas c1("c1","c1", w, h);
  c1.cd();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1);
  hintxy_discr->Draw("colz");
  c1.Print("ric_disc.png");
  c1.Print("ric_disc.root");
  
  /***************************************************************
 Sovrapposizione istogramma proiezione sulle x generati VS discreti
  ****************************************************************/
  TCanvas c2("c2","c2", w, h);
  c2.cd();
  int tot_x_gen=hintx_gen->Integral(); //Normalizziamo tutto a 1
  hintx_gen->Scale(1./tot_x_gen);
  hintx_gen->Draw();
  int tot_x_discr=hintx_discr->Integral();
  hintx_discr->Scale(1./tot_x_discr);
  hintx_discr->SetLineColor(2);
  hintx_discr->Draw("same");
  c2.Update();
  c2.Print("genVSdiscrX.png");
  c2.Print("genVSdiscrX.root");
  
  /***************************************************************
Sovrapposizione istogramma proiezione sulle x generati VS discreti
  ****************************************************************/
  TCanvas c3("c3","c3", w, h);
  c3.cd();
  int tot_y_gen=hinty_gen->Integral();
  hinty_gen->Scale(1./tot_y_gen);
  hinty_gen->Draw();
  hinty_discr->SetLineColor(2);
  int tot_y_discr=hinty_discr->Integral();
  hinty_discr->Scale(1./tot_y_discr);
  hinty_discr->Draw("same");
  c3.Update();
  c3.Print("genVSdiscrY.png");
  c3.Print("genVSdiscrY.root");

  /***************************************************************
        istogramma 2D scala a colori eventi generati
  ****************************************************************/
  TCanvas c6("c6","c6", w, h);
  c6.cd();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1);
  hintxy_gen->Draw("colz");
  c6.Print("reco2D_gen.png");
  c6.Print("reco2D_gen.root");

  /***************************************************************
        istogramma 2D scala a colori Punti rivelati 
  ****************************************************************/
  TCanvas c7("c7","c7", w, h);
  c7.cd();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1);
  hpuntiriv->Draw("colz");
  c7.Print("ptiriv.png");
  c7.Print("ptiriv.root");

  /***************************************************************
        scrittura degli histogrammi
  ****************************************************************/

  char* out_file_root = new char[256];
  out_file_root = (char*)"toy.root";
  TFile foutput(out_file_root,"recreate");
  //foutput.Write();
  
  // hintx->Write();
  // hinty->Write();
  hintx_discr->Write();
  hinty_discr->Write();
  // hintxy->Write();
  hintx_gen->Write();
  hinty_gen->Write();
  hintxy_gen->Write();
  hpuntiriv->Write();
  foutput.Close();
      
}

/**********************************************************/

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

/*****************************************************
	intersezione
*****************************************************/

TVector2 intersezione (TVector2 a1, TVector2 a2, TVector2 b1, TVector2 b2)
{
  TVector2 inters;
  Double_t ics;
  Double_t ipsilon;
  
  if ((a2.X()-a1.X())==0)
    {
      ics=a1.X();
      ipsilon=b1.Y()+(b2.Y()-b1.Y())/(b2.X()-b1.X())*(a1.X()-b1.X());
    }
  else if ((b2.Y()-b1.Y())==0)
    {
      ipsilon=b2.Y();
      ics=a1.X()+(a2.X()-a1.X())/(a2.Y()-a1.Y())*(ipsilon-a1.Y());
    }
  else{
    ipsilon=(a1.Y()+((a2.Y()-a1.Y())/(a2.X()-a1.X()))*(b1.X()-a1.X()-b1.Y()*((b2.X()-b1.X())/(b2.Y()-b1.Y()))))/(1.-((a2.Y()-a1.Y())/(a2.X()-a1.X()))*((b2.X()-b1.X())/(b2.Y()-b1.Y())));
    ics=(ipsilon-a1.Y())*(a2.X()-a1.X())/(a2.Y()-a1.Y())+a1.X();
  }
  inters.Set(ics,ipsilon);
  //inters.Print();
  return inters;
}

/**********************************************************/


TVector3 sghembepoint(TVector3 p11, TVector3 p12, TVector3 p21, TVector3 p22 ) {

  TVector3 dir1 = (p12-p11).Unit();
  TVector3 dir2 = (p22-p21).Unit();	
	
	
	
}

/**********************************************************
               Sezione d'urto Compton
**********************************************************/

Double_t Sigma_Comp(Double_t gamma){

  Double_t PI=TMath::Pi();
  Double_t re=2.8179402e-13;//raggio e- in cm
  Double_t sigma_comp;

  sigma_comp=2*PI*re*re*((1.+gamma)/(gamma*gamma)*((2.*(1.+gamma))/(1.+2.*gamma)-log(1.+2.*gamma)/gamma )+log(1.+2.*gamma)/(2.*gamma)-(1.+3.*gamma)/(pow(1.+2.*gamma,2)));

  return sigma_comp;
}

/**********************************************************
               generazione angolo Compton
**********************************************************/

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

/**********************************************************
       funzioni  di propagazione dei fotoni
**********************************************************/

TVector2 propagazionesenzacompton(TVector2 posizione, TVector2 direzione, Double_t raggio, Double_t nriv) {

   TVector2 assex(1,0);
   TVector2 assey(0,1);

   Double_t  bi=2*posizione.X()*(direzione.Unit()*assex)+2*posizione.Y()*(direzione.Unit()*assey); //determinazione dei punti di impatto sul rivelatore
   Double_t ci=pow(posizione.X(),2)+pow(posizione.Y(),2)-raggio*raggio;
   Double_t alpha=(-bi+sqrt(bi*bi-4*ci))/2;
   posizione=posizione+alpha*direzione;

  // Double_t angolodiscretizzato = (floor(posizione.Phi()*nriv/(2*TMath::Pi()))+0.5)/(nriv/(2*TMath::Pi()));//discretizzazione angolo
  // posizione.SetMagPhi(posizione.Mod(),angolodiscretizzato);
     
  return posizione;
}


TVector2 propagazionecompton(TVector2 posizione, TVector2 direzione, Double_t rpaz, Double_t raggio, Double_t nriv, Double_t taglio_en) {

  Double_t bi, ci, alpha;

  Double_t zeta_fot;
  Double_t ang_comp;

  Double_t angolodiscretizzato;

  Double_t n1g=1.;//densità acqua (g/cm^3)
  Double_t n2g=1.205e-3;//densità aria (g/cm^3)

  Double_t n1=n1g*TMath::Na()*.55509;//num e- in acqua
  Double_t n2=n2g*TMath::Na()*.49976;//num e- in aria

  TVector2 assex(1,0);
  TVector2 assey(0,1);
  TVector2 vettorenullo(0,0);

  Double_t en=511;
  Double_t me=511;
  Double_t gamma=1;
  Double_t sigma=Sigma_Comp(gamma);   
  Double_t lambda_H2O=1./(sigma*n1);
cout << "c" << lambda_H2O <<endl;
  Double_t lambda_N=1./(sigma*n2);


  while(posizione.Mod() < rpaz) //controllo se è all'interno del corpo 
    {
      lambda_H2O=1./(sigma*n1);
      zeta_fot= - lambda_H2O * log(random1());//spazio percorso dal fotone prima di interagire
      if ((posizione+zeta_fot*direzione).Mod()<rpaz) 
	{
	  posizione = posizione+zeta_fot*direzione;
	  ang_comp= Ang_gen(gamma);//genero angolo 
	  direzione=direzione.Rotate(ang_comp);
	  en=en/(1.+(en/me)*(1.-cos(ang_comp)));
	  gamma=en/me;
	  sigma=Sigma_Comp(gamma);  	
	} 
      else 
       {
	bi=2*posizione.X()*(direzione.Unit()*assex)+2*posizione.Y()*(direzione.Unit()*assey); //uscita fotoni dal paziente
	ci=pow(posizione.X(),2)+pow(posizione.Y(),2)-rpaz*rpaz;
	alpha=(-bi+sqrt(bi*bi-4*ci))/2;
	posizione=posizione+alpha*direzione;
       }
    }


  while(posizione.Mod() < raggio) //controllo se fa compton in aria
    {
      lambda_N=1./(sigma*n2);
      zeta_fot= - lambda_N * log(random1());

      if ((posizione+zeta_fot*direzione).Mod()<raggio)
	{
	  posizione=posizione+zeta_fot*direzione;
	  ang_comp= Ang_gen(gamma);
	  direzione=direzione.Rotate(ang_comp);//genero angolo 
	  en=en/(1.+(en/me)*(1.-cos(ang_comp)));
	  gamma=en/me;
	  sigma=Sigma_Comp(gamma);  
	}
      else {
	bi=2*posizione.X()*(direzione.Unit()*assex)+2*posizione.Y()*(direzione.Unit()*assey); //uscita fotoni dal paziente
	ci=pow(posizione.X(),2)+pow(posizione.Y(),2)-raggio*raggio;
	alpha=(-bi+sqrt(bi*bi-4*ci))/2;
	posizione=posizione+alpha*direzione;
      }
    }
 
  en=en+.1*en*rgauss();//smearing energia
  if (en>taglio_en) //taglio energia
     {
       angolodiscretizzato = (floor(posizione.Phi()*nriv/(2*TMath::Pi()))+0.5)/(nriv/(2*TMath::Pi()));//discretizzazione angolo
       posizione.SetMagPhi(posizione.Mod(),angolodiscretizzato);
       return posizione;
     }
  else return vettorenullo;
}

