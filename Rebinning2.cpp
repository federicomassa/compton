// Programma per modificare il bin degli istogrammi del Compton. Si tratta di accorpare n canali assieme, a partire da 8192 canali totali.
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;


float calibration(float ch) {
// calibrazione: energia = a + b*ch + c*ch^2 + d*ch^3
  double a = -33.0933;
  double b = 0.103793;
  double c = -1.80489E-5;
  double d = 4.03119E-9;
  return (a+b*ch+c*pow(ch,2)+d*pow(ch,3));
      }

void Rebinning2(){
    int bin;
cout << "Nome file input con estensione esplicita: " << endl;
char finname[50],foutname[50];
cin >> finname;
cout << "Nome file output con estensione esplicita: " << endl;
cin >> foutname;
cout << "Quanti bin accorpare? (1 = non fa nulla)" << endl;
cin >> bin;
string linein;
int lineout = 0, partial = 0;
ifstream fin(finname);
ofstream fout(foutname);
int j = 0;
while (!fin.eof()) {
      int h = 0;
      for (int i = 1; i <= bin && !fin.eof(); i++){
      getline(fin,linein);
      stringstream(linein) >> partial;
      lineout += partial;
      j += 1;
      h += 1;
      }
      //(2j-bin+1)/2 è la media aritmetica tra i canali
      if (h == bin && !fin.eof()) {
      fout << calibration(float((2*j-bin+1))/2.) << '\t' << lineout << endl;}
      //fout << lineout << endl;}
      lineout = 0;
      
      } 
fin.close();
fout.close();
}
