#include <iostream>
using namespace std;
int main() 
{

	float Vinf, Atmos_temp, Wall_temp, X_O2, X_N2, X_O, X_N, X_H, X_He;

	Vinf = 7500.0;        // Space Object velocity neglecting winds, Vinf, m/s
	Atmos_temp = 400.0;   // Atmospheric Gas Temperature, Tinf, K
	Wall_temp = 300.0;    // Space object wall temperature, Tw, K

	X_O2 = 0.0062;      // Mass fraction of Molecular Oxygen	
	X_N2 = 0.1340;      // Mass fraction of Molecular Nitrogen
	X_O = 0.8359;       // Mass fraction of Atomic Oxygen
	X_N = 0.0031;       // Mass fraction of Atomic Nitrogen
	X_H = 0.0003;       // Mass fraction of Hydrogen
	X_He = 0.0205;      // Mass fraction of Helium	
	
	double X[6] = {X_O2, X_N2, X_O, X_N, X_H, X_He};


  	cout << X[0] << endl;
  	cout << X[1] << endl;
  	cout << X[2] << endl;
  	cout << X[3] << endl;
  	cout << X[4] << endl;
  	cout << X[5] << endl;

  	return 0;
}