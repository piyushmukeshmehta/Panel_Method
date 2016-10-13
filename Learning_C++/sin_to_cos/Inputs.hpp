void MaxPressSA(double Vinf, double Tinf, double X[6], double m_surf, double &s, double &mu);

int read_num_lines();

struct facet_struct {
	double normal[3];
	double vertex1[3];
	double vertex2[3];
	double vertex3[3];
	double centroid[3];
	double area;
};

void facet_properties(int nfacets, struct facet_struct *pfacet);

double dot(double V[], double W[]);

void cross(double V[], double W[], double VEC[]);

void LFACSA(double Vinf, double Tinf, double Tw, double L_char, double alpha, double beta, double theta_sc, double s, double mu, int nfacets, struct facet_struct *pfacet, double &CD, double &CS, double &CL, double &Cl, double &Cm, double &Cn);