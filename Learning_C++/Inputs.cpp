#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Inputs.hpp"

#define filename "sphere2.stl"
using namespace std;

// command to compile
// g++ Inputs.cpp -o Inputs -Wall -O -lm -lgsl

int main() 
{

	double Vinf = 7500.0;        // Space Object velocity neglecting winds, Vinf, m/s
	double Tinf = 400.0;   // Atmospheric Gas Temperature, Tinf, K
	double Tw = 300.0;    // Space object wall temperature, Tw, K

	double X_O2 = 0.0062;      // Mass fraction of Molecular Oxygen	
	double X_N2 = 0.1340;      // Mass fraction of Molecular Nitrogen
	double X_O = 0.8359;       // Mass fraction of Atomic Oxygen
	double X_N = 0.0031;       // Mass fraction of Atomic Nitrogen
	double X_H = 0.0003;       // Mass fraction of Hydrogen
	double X_He = 0.0205;      // Mass fraction of Helium	
	
	double m_surf = 27;        // Mass of Surface Material, 27 for Aluminum
	double theta_sc = 0.8;

	// Object data
	double L_char = 5.092;     // Characteristic Length of the object for Moment calculations

	// Attitude data
	double alpha = 25.0 * M_PI / 180;      // Pitch Angle, rad
	double beta = 025.0 * M_PI / 180;       // Yaw Angle, rad

	double X[6] = {X_O2, X_N2, X_O, X_N, X_H, X_He};

	double s = 0.0;
	double mu = 0.0;

	MaxPressSA(Vinf, Tinf, X, m_surf, s, mu);

	 //printf("%1.2f\n%1.2f\n", s, mu);

/* READ and DETERMINE NUMBER OF FACETS FROM MESH FILE */
  	int nfacets = 0;
  	nfacets = read_num_lines();

/* READ IN FACET PROPERTIES FROM MESH FILE */

  	struct facet_struct *pfacet = NULL;
	pfacet = (struct facet_struct *) calloc(nfacets, sizeof(struct facet_struct));
  
  	facet_properties(nfacets, pfacet);

   //  printf("%e\n", (pfacet+nfacets-1)->normal[0]);
  	// printf("%e\n", (pfacet+nfacets-1)->normal[1]);
  	// printf("%e\n", (pfacet+nfacets-1)->normal[2]);

  double CD = 0.0;
  double CS = 0.0;
  double CL = 0.0;
  double Cl = 0.0;
  double Cm = 0.0;
  double Cn = 0.0;

  LFACSA(Vinf, Tinf, Tw, L_char, alpha, beta, theta_sc, s, mu, nfacets, pfacet, CD, CS, CL, Cl, Cm, Cn);

	 printf("%e\n", CD);
   printf("%e\n", CS);
   printf("%e\n", CL);
   printf("%e\n", Cl);
   printf("%e\n", Cm);
   printf("%e\n", Cn);

  	return 0;
}

// This function calculated the speed ratio which is an input to the analytical models in function LFACSA

void MaxPressSA(double Vinf, double Tinf, double X[6], double m_surf, double &s, double &mu)
{

	double k = 1.3806488e-23;  // Boltzmann constant, m2 kg s-2 K-1
	double amu = 1.66053892e-27;       // atomic mass unit
	double mO2 = 32.0;       // Molecular Weight of O2
	double mN2 = 28.0;       // Molecular Weight of N2
	double mO = 16.0;        // Molecular Weight of O
	double mN = 14.0;        // Molecular Weight of N
	double mH = 1.0;         // Molecular Weight of H
	double mHe = 4.0;        // Molecular Weight of He

	double mg =  X[0] * mO2 + X[1] * mN2 + X[2] * mO + X[3] * mN + X[4] * mH + X[5] * mHe;

	double vmp = sqrt(2*k*Tinf/mg/amu);        // Most Probable speed, m/s
	s = Vinf/vmp;             // Speed Ratio

	mu = mg/m_surf;

}

/******************************** MESH READ ***************************/
int read_num_lines()
{

  /* Read in the STL mesh file and determine number of lines*/
  /* Input: STL mesh filename */
  /* Output: Number of Facets */

  int nfacets = 0;
  int num_lines = 0;
  int ch;
  FILE *f = fopen(filename, "r");
  
  /*Check that file exists*/
  if(!f) {
    printf("Mesh File does not exist\n");
    exit(1);
  }
  
  while (EOF != (ch=fgetc(f))) 
    if (ch=='\n')
      num_lines++;

  fclose(f);

  /* DETERMINE NUMBER OF FACETS */
  nfacets = (num_lines-2)/7;

  //cout << nfacets << endl;

  return(nfacets);

}


/******************************** FACET PROPERTIES ***************************/
void facet_properties(int nfacets, struct facet_struct *pfacet)
{

  struct facet_struct *facet;

  /* Input: STL mesh filename */
  /* Output: Facet Properties Stucture Containing: */
  /*         Facet Normal [x, y, z]  */
  /*         Vertex1 [x, y, z] */
  /*         Vertex2 [x, y, z] */
  /*         Vertex3 [z, y, z] */

  FILE *f = fopen(filename, "r");

  /* READ IN STL FILE HEADER */
  char header[1024];
  char line1[1024], line2[1024], line3[1024], line4[1024], line5[1024], line6[1024], line7[1024];
  char *vert1x, *vert1y, *vert1z;
  char *vert2x, *vert2y, *vert2z;
  char *vert3x, *vert3y, *vert3z;
  char *normx, *normy, *normz;
  char *temp;
  int HeaderSize = 1;
  int i, ifacet;
  double temp_area, v1, v2, v3, v4;

  double dist1[3], dist2[3], dist3[3];
  double tc[3];
  double Len[3];

  for(i=0; i<3; i++) {
    dist1[i] = 0.0;
    dist2[i] = 0.0;
    dist3[i] = 0.0;
    tc[i] = 0.0;
  }

  for(i=0; i<HeaderSize; i++) {
    fgets(header, 1024, f);
  }

  for(ifacet=0; ifacet<nfacets; ifacet++) {
    /* Read the first line and assign normal components */
    fgets(line1, 1024, f);
    temp = strtok(line1, " ");
    temp = strtok(NULL, " ");
    normx = strtok(NULL, " ");
    normy = strtok(NULL, " ");
    normz = strtok(NULL, " ");
    /* Throw away second line */
    fgets(line2, 1024, f);
    /* Read the third line and assign vertex #1 position */
    fgets(line3, 1024, f);
    temp = strtok(line3, " ");
    vert1x = strtok(NULL, " ");
    vert1y = strtok(NULL, " ");
    vert1z = strtok(NULL, " ");
    /* Read the fourth line and assign vertex #2 position */
    fgets(line4, 1024, f);
    temp = strtok(line4, " ");
    vert2x = strtok(NULL, " ");
    vert2y = strtok(NULL, " ");
    vert2z = strtok(NULL, " ");
    /* Read the fifth line and assign vertex #3 position */
    fgets(line5, 1024, f);
    temp = strtok(line5, " ");
    vert3x = strtok(NULL, " ");
    vert3y = strtok(NULL, " ");
    vert3z = strtok(NULL, " ");
    /* Throw away the sixth and seventh lines */
    fgets(line6, 1024, f);
    fgets(line7, 1024, f);

   	facet = pfacet + ifacet;

    /* Assign facet normal vector to "facet" structure */
    facet->normal[0] = atof(normx);
    facet->normal[1] = atof(normy);
    facet->normal[2] = atof(normz);

    /* Assign vertex #1 position to "facet" structure */
    facet->vertex1[0] = atof(vert1x);
    facet->vertex1[1] = atof(vert1y);
    facet->vertex1[2] = atof(vert1z);

    /* Assign vertex #2 position to "facet" structure */
    facet->vertex2[0] = atof(vert2x);
    facet->vertex2[1] = atof(vert2y);
    facet->vertex2[2] = atof(vert2z);

    /* Assign vertex #3 position to "facet" structure */
    facet->vertex3[0] = atof(vert3x);
    facet->vertex3[1] = atof(vert3y);
    facet->vertex3[2] = atof(vert3z);

    /* Find 1st side of the triangle */
    dist1[0] = facet->vertex1[0] - facet->vertex2[0];
    dist1[1] = facet->vertex1[1] - facet->vertex2[1];
    dist1[2] = facet->vertex1[2] - facet->vertex2[2];
    
    /* Find 2nd side of the triangle */
    dist2[0] = facet->vertex1[0] - facet->vertex3[0];
    dist2[1] = facet->vertex1[1] - facet->vertex3[1];
    dist2[2] = facet->vertex1[2] - facet->vertex3[2];

    /* Find 3rd side of the triangle */
    dist3[0] = facet->vertex3[0] - facet->vertex2[0];
    dist3[1] = facet->vertex3[1] - facet->vertex2[1];
    dist3[2] = facet->vertex3[2] - facet->vertex2[2];

    /* Compute the length of the 3 sides of the triangle */
    Len[0] = sqrt(dist1[0]*dist1[0] + dist1[1]*dist1[1] + dist1[2]*dist1[2]);
    Len[1] = sqrt(dist2[0]*dist2[0] + dist2[1]*dist2[1] + dist2[2]*dist2[2]);
    Len[2] = sqrt(dist3[0]*dist3[0] + dist3[1]*dist3[1] + dist3[2]*dist3[2]);

    // /* Compute the centroid of the facets and assign to structure */

    // facet->centroid[0] = (facet->vertex1[0] + facet->vertex1[1] + facet->vertex1[2])/3.0;
    // facet->centroid[1] = (facet->vertex2[0] + facet->vertex2[1] + facet->vertex2[2])/3.0;
    // facet->centroid[2] = (facet->vertex3[0] + facet->vertex3[1] + facet->vertex3[2])/3.0;

    /* Compute the incenter of the facets and assign to structure */

    facet->centroid[0] = (Len[0] * facet->vertex3[0] + Len[1] * facet->vertex2[0] + Len[2] * facet->vertex1[0]) / (Len[0] + Len[1] + Len[2]);
    facet->centroid[1] = (Len[0] * facet->vertex3[1] + Len[1] * facet->vertex2[1] + Len[2] * facet->vertex1[1]) / (Len[0] + Len[1] + Len[2]);
    facet->centroid[2] = (Len[0] * facet->vertex3[2] + Len[1] * facet->vertex2[2] + Len[2] * facet->vertex1[2]) / (Len[0] + Len[1] + Len[2]);

    // Compute the area using modified and Stabilized Heron's formula http://http.cs.berkeley.edu/%7Ewkahan/Triangle.pdf
    sort(&Len[0], &Len[0]+3);

    double L1 = Len[2];
    double L2 = Len[1];
    double L3 = Len[0];

    temp_area = L2 + L3;
    v1 = L1 + temp_area;
    temp_area = L1 - L2;
    v2 = L3 - temp_area;
    v3 = L3 + temp_area;
    temp_area = L2 - L3;
    v4 = L1 + temp_area;   

    // Area using Heron Formula
    facet->area = 0.25 * sqrt(v1*v2*v3*v4);   

  }
   fclose(f);
}


// This function reads the stl file and computes the force and moment coefficients

//void LFACSA(double Vinf, double Tinf, double Tw, double L_char, double alpha, double beta, double theta_sc, double s, double mu, double &CD, double &CS, double &CL, double &Cl, double &Cm, double &Cn, int nfacets)
void LFACSA(double Vinf, double Tinf, double Tw, double L_char, double alpha, double beta, double theta_sc, double s, double mu, int nfacets, struct facet_struct *pfacet, double &CD, double &CS, double &CL, double &Cl, double &Cm, double &Cn )

{

double Vinf_vec[3] = {Vinf * cos(alpha) * cos(beta), -1.0 * Vinf * cos(alpha) * sin(beta), Vinf * sin(alpha)};		// Velocity Vector
double Vmag = sqrt(Vinf_vec[0]*Vinf_vec[0] + Vinf_vec[1]*Vinf_vec[1] + Vinf_vec[2]*Vinf_vec[2]);

// Computing the Unit Velocity Vector
double Vinfn[3] = {0,0,0};
double Vinfn_neg[3] = {0,0,0};

Vinfn[0] = Vinf_vec[0] / Vmag;
Vinfn[1] = Vinf_vec[1] / Vmag; 
Vinfn[2] = Vinf_vec[2] / Vmag; 

Vinfn_neg[0] = -1.0 * Vinf_vec[0] / Vmag;
Vinfn_neg[1] = -1.0 * Vinf_vec[1] / Vmag; 
Vinfn_neg[2] = -1.0 * Vinf_vec[2] / Vmag;

// ****************** Accommodation Coefficient Models *****************

double ST_ads = 1.0;   // Tangential Momentum Accommodation Coefficient for adsorbate
double SN_ads = 1.0;   // Normal Momentum Accommodation coefficient for adsorbate

double ST_s = 1.0;       // Tangential Momentum Accommodation Coefficient for clean surface
double AT_s = ST_s * (2 - ST_s);     // Tangental Energy Accomodation Coefficient for clean surface
double AC_s = 2.4 * mu / ( (1+mu) * (1+mu) );     // Goodman's Accommodation Model for clean surface
double AN_s = 2 * AC_s - AT_s;   // Normal Energy Accommodation Coefficient for clean surface
double SN_s = 1 - sqrt(1 - AN_s);  // Normal Momentum Accommodation coefficient for clean surface

// ****************** Parameter Initialization *****************

double Sref = 0;
double theta = 0;
double Cpn_s[3] = {0,0,0}, Ctt_s[3] = {0,0,0},C_s[3] = {0,0,0};
double Cpn_ads[3] = {0,0,0}, Ctt_ads[3] = {0,0,0}, C_ads[3] = {0,0,0};
double CG[3] = {1.3333, 0, 0.2414}, rad_dist[3] = {0,0,0};
double cros_prod_s[3] = {0,0,0}, cros_prod_ads[3] = {0,0,0};
double sum_cp_ct_s[3] = {0,0,0}, sum_cp_ct_ads[3] = {0,0,0};
double M_s[3] = {0,0,0}, M_ads[3] = {0,0,0};

int i;

	for (i=0; i<nfacets; i++){

		theta = asin(dot(Vinfn_neg, (pfacet+i)->normal));

    if(theta > 0){
      Sref += (pfacet+i)->area * sin(theta);
    }

    rad_dist[0] = (pfacet+i)->centroid[0] - CG[0];
    rad_dist[1] = (pfacet+i)->centroid[1] - CG[1];
    rad_dist[2] = (pfacet+i)->centroid[2] - CG[2];

    double t_s[3] = {0,0,0}, t_ads[3] = {0,0,0};
    double Ct_s = 0;
    double Ct_ads = 0;
    double Cp_s = 0;
    double Cp_ads = 0;

    double Vn = dot(Vinfn,(pfacet+i)->normal);    // product of Vinf and normal vector

    if(theta == 0){

      // computing 't' vector for clean surface, http://www.ssdl.gatech.edu/papers/conferencePapers/AIAA-2014-0728.pdf
      t_s[0] = (Vn * ((pfacet+i)->normal[0]) - Vinfn[0]) / sqrt(1-Vn*Vn);
      t_s[1] = (Vn * ((pfacet+i)->normal[1]) - Vinfn[1]) / sqrt(1-Vn*Vn);
      t_s[2] = (Vn * ((pfacet+i)->normal[2]) - Vinfn[2]) / sqrt(1-Vn*Vn);

      Ct_s = -1.0 * (ST_s * cos(theta)/ s / sqrt(M_PI)) * (exp(-1.0*(s*sin(theta))*(s*sin(theta)) ) + sqrt(M_PI) * s *sin(theta)*(1.0+erf(s*sin(theta))));
      
      // computing 't' vector for adsorbate, http://www.ssdl.gatech.edu/papers/conferencePapers/AIAA-2014-0728.pdf
      t_ads[0] = (Vn * ((pfacet+i)->normal[0]) - Vinfn[0]) / sqrt(1-Vn*Vn);
      t_ads[1] = (Vn * ((pfacet+i)->normal[1]) - Vinfn[1]) / sqrt(1-Vn*Vn);
      t_ads[2] = (Vn * ((pfacet+i)->normal[2]) - Vinfn[2]) / sqrt(1-Vn*Vn);

      Ct_ads = -1.0 * (ST_ads * cos(theta)/ s / sqrt(M_PI)) * (exp(-1.0*(s*sin(theta))*(s*sin(theta)) ) + sqrt(M_PI) * s *sin(theta)*(1.0+erf(s*sin(theta))));
    }

    if(theta == 0.5 * M_PI){

       double Cp1_s = exp(-1.0*(s*sin(theta))*(s*sin(theta))) * ((2.0-SN_s)*s*sin(theta)/sqrt(M_PI) + SN_s*sqrt(Tw/Tinf)/2);
       double Cp2_s = (1.0+erf(s*sin(theta))) * ((2.0-SN_s) * ((s*sin(theta))*(s*sin(theta)) + 0.5) + 0.5*SN_s*s*sin(theta)*sqrt(M_PI*Tw/Tinf));
        Cp_s = (Cp1_s+Cp2_s)/(s*s);
        
       double Cp1_ads = exp(-1.0*(s*sin(theta))*(s*sin(theta))) * ((2.0-SN_ads)*s*sin(theta)/sqrt(M_PI) + SN_s*sqrt(Tw/Tinf)/2);
       double Cp2_ads = (1.0+erf(s*sin(theta))) * ((2.0-SN_ads) * ((s*sin(theta))*(s*sin(theta)) + 0.5) + 0.5*SN_ads*s*sin(theta)*sqrt(M_PI*Tw/Tinf));
        Cp_ads = (Cp1_ads+Cp2_ads)/(s*s);
    }

    if(theta != 0 && theta != 0.5*M_PI && theta != - 0.5 * M_PI){

      // computing 't' vector for clean surface, http://www.ssdl.gatech.edu/papers/conferencePapers/AIAA-2014-0728.pdf
      t_s[0] = (Vn * ((pfacet+i)->normal[0]) - Vinfn[0]) / sqrt(1-Vn*Vn);
      t_s[1] = (Vn * ((pfacet+i)->normal[1]) - Vinfn[1]) / sqrt(1-Vn*Vn);
      t_s[2] = (Vn * ((pfacet+i)->normal[2]) - Vinfn[2]) / sqrt(1-Vn*Vn);

      Ct_s = -1.0 * (ST_s * cos(theta)/ s / sqrt(M_PI)) * (exp(-1.0*(s*sin(theta))*(s*sin(theta)) ) + sqrt(M_PI) * s *sin(theta)*(1.0+erf(s*sin(theta))));

      double Cp1_s = exp(-1.0*(s*sin(theta))*(s*sin(theta))) * ((2.0-SN_s)*s*sin(theta)/sqrt(M_PI) + SN_s*sqrt(Tw/Tinf)/2);
      double Cp2_s = (1.0+erf(s*sin(theta))) * ((2.0-SN_s) * ((s*sin(theta))*(s*sin(theta)) + 0.5) + 0.5*SN_s*s*sin(theta)*sqrt(M_PI*Tw/Tinf));
        Cp_s = (Cp1_s+Cp2_s)/(s*s);

      // computing 't' vector for adsorbate, http://www.ssdl.gatech.edu/papers/conferencePapers/AIAA-2014-0728.pdf
      t_ads[0] = (Vn * ((pfacet+i)->normal[0]) - Vinfn[0]) / sqrt(1-Vn*Vn);
      t_ads[1] = (Vn * ((pfacet+i)->normal[1]) - Vinfn[1]) / sqrt(1-Vn*Vn);
      t_ads[2] = (Vn * ((pfacet+i)->normal[2]) - Vinfn[2]) / sqrt(1-Vn*Vn);

      Ct_ads = -1.0 * (ST_ads * cos(theta)/ s / sqrt(M_PI)) * (exp(-1.0*(s*sin(theta))*(s*sin(theta)) ) + sqrt(M_PI) * s *sin(theta)*(1.0+erf(s*sin(theta))));

      double Cp1_ads = exp(-1.0*(s*sin(theta))*(s*sin(theta))) * ((2.0-SN_ads)*s*sin(theta)/sqrt(M_PI) + SN_s*sqrt(Tw/Tinf)/2);
      double Cp2_ads = (1.0+erf(s*sin(theta))) * ((2.0-SN_ads) * ((s*sin(theta))*(s*sin(theta)) + 0.5) + 0.5*SN_ads*s*sin(theta)*sqrt(M_PI*Tw/Tinf));
        Cp_ads = (Cp1_ads+Cp2_ads)/(s*s);
    }

    // Clean Surface Contribution
    Cpn_s[0] = Cp_s * -1.0 * ((pfacet+i)->normal[0]) * ((pfacet+i)->area);
    Cpn_s[1] = Cp_s * -1.0 * ((pfacet+i)->normal[1]) * ((pfacet+i)->area);
    Cpn_s[2] = Cp_s * -1.0 * ((pfacet+i)->normal[2]) * ((pfacet+i)->area);

    Ctt_s[0] = Ct_s * t_s[0] * ((pfacet+i)->area);
    Ctt_s[1] = Ct_s * t_s[1] * ((pfacet+i)->area);
    Ctt_s[2] = Ct_s * t_s[2] * ((pfacet+i)->area);

    C_s[0] += Cpn_s[0] + Ctt_s[0];
    C_s[1] += Cpn_s[1] + Ctt_s[1];
    C_s[2] += Cpn_s[2] + Ctt_s[2];

    sum_cp_ct_s[0] = Cpn_s[0] + Ctt_s[0];
    sum_cp_ct_s[1] = Cpn_s[1] + Ctt_s[1];
    sum_cp_ct_s[2] = Cpn_s[2] + Ctt_s[2];

    cross(rad_dist,sum_cp_ct_s,cros_prod_s);

    M_s[0] += cros_prod_s[0];
    M_s[1] += cros_prod_s[1];
    M_s[2] += cros_prod_s[2];

    // Adsorbate Contribution
    Cpn_ads[0] = Cp_ads * -1.0 * ((pfacet+i)->normal[0]) * ((pfacet+i)->area);
    Cpn_ads[1] = Cp_ads * -1.0 * ((pfacet+i)->normal[1]) * ((pfacet+i)->area);
    Cpn_ads[2] = Cp_ads * -1.0 * ((pfacet+i)->normal[2]) * ((pfacet+i)->area);

    Ctt_ads[0] = Ct_ads * t_s[0] * ((pfacet+i)->area);
    Ctt_ads[1] = Ct_ads * t_s[1] * ((pfacet+i)->area);
    Ctt_ads[2] = Ct_ads * t_s[2] * ((pfacet+i)->area);

    C_ads[0] += Cpn_ads[0] + Ctt_ads[0];
    C_ads[1] += Cpn_ads[1] + Ctt_ads[1];
    C_ads[2] += Cpn_ads[2] + Ctt_ads[2];

    sum_cp_ct_ads[0] = Cpn_ads[0] + Ctt_ads[0];
    sum_cp_ct_ads[1] = Cpn_ads[1] + Ctt_ads[1];
    sum_cp_ct_ads[2] = Cpn_ads[2] + Ctt_ads[2];

    cross(rad_dist,sum_cp_ct_ads,cros_prod_ads);

    M_ads[0] += cros_prod_ads[0];
    M_ads[1] += cros_prod_ads[1];
    M_ads[2] += cros_prod_ads[2];

    //printf("%e\n", t_s[0], "%e\n", t_s[1], "%e\n", t_s[2],);
    //printf("%e\n", (pfacet+i)->centroid[0]);
    // printf("%e\n", rad_dist[0]);
	}

  printf("%e\n", Sref);
  // printf( "%e\n", M_s[1]);
  // printf( "%e\n", M_s[2]);

  C_s[0] = C_s[0]/Sref;
  C_s[1] = C_s[1]/Sref;
  C_s[2] = C_s[2]/Sref;
  
  M_s[0] = M_s[0]/Sref/L_char;
  M_s[1] = M_s[1]/Sref/L_char;
  M_s[2] = M_s[2]/Sref/L_char;

  C_ads[0] = C_ads[0]/Sref;
  C_ads[1] = C_ads[1]/Sref;
  C_ads[2] = C_ads[2]/Sref;
  
  M_ads[0] = M_ads[0]/Sref/L_char;
  M_ads[1] = M_ads[1]/Sref/L_char;
  M_ads[2] = M_ads[2]/Sref/L_char;

  double B2WA_X[3] = {cos(alpha)*cos(beta), -1.0*sin(beta)*cos(alpha), sin(alpha)};
  double B2WA_Y[3] = {sin(beta), cos(beta), 0};
  double B2WA_Z[3] = {-1.0*sin(alpha)*cos(beta), sin(alpha)*sin(beta), cos(alpha)};

  double CD_s, CD_ads, CL_s, CL_ads, CS_s, CS_ads, Cl_s, Cl_ads, Cm_s, Cm_ads, Cn_s, Cn_ads;

  CD_s = dot(C_s, B2WA_X);
  CS_s = dot(C_s, B2WA_Y);
  CL_s = dot(C_s, B2WA_Z);

  CD_ads = dot(C_ads, B2WA_X);
  CS_ads = dot(C_ads, B2WA_Y);
  CL_ads = dot(C_ads, B2WA_Z);

  Cl_s = M_s[0];
  Cm_s = M_s[1];
  Cn_s = M_s[2];

  Cl_ads = M_ads[0];
  Cm_ads = M_ads[1];
  Cn_ads = M_ads[2];

    CD = (1.0-theta_sc) * CD_s + theta_sc * CD_ads;
    CS = (1.0-theta_sc) * CS_s + theta_sc * CS_ads;
    CL = (1.0-theta_sc) * CL_s + theta_sc * CL_ads;

    Cl = (1.0-theta_sc) * Cl_s + theta_sc * Cl_ads;
    Cm = (1.0-theta_sc) * Cm_s + theta_sc * Cm_ads;
    Cn = (1.0-theta_sc) * Cn_s + theta_sc * Cn_ads;
  
}

/*********************************** DOT PRODUCT *************************************/
double dot(double V[], double W[]) {
   
  /* Computes dot product between vectors V and W */
  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */
  /* Output: a = Dot Product of V and W */

  double a = 0.0;

  a = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];
  return(a);

}

/*********************************** CROSS PRODUCT *************************************/
void cross(double V[], double W[], double VEC[]) {
  /* Computes cross product between vectors V and W */

  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */

  /* Output: b = Cross Product of V and W */

  VEC[0] = V[1]*W[2] - V[2]*W[1];
  VEC[1] = V[2]*W[0] - V[0]*W[2];
  VEC[2] = V[0]*W[1] - V[1]*W[0];

}
