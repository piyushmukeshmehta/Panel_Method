#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <hdf5.h>

#include "SpiceUsr.h"
#include "defs.h"
#include "structs.h"
#include "io.h"
#include "dirent.h"
#include "hdf5.h"

//extern double sa_altitude[NVERT];
//extern double sa_density[NVERT];

void read_files(double tsearch[2], char *dirname, double etlist[MAXFILES], char files[MAXFILES][100], int *Nfiles) {

  /* Function to read in the GITM densities from a text file and output a */
  /* matlab structure with the density at different locations. */
  
  /* Usage: */
  /* GITM = importgitm(t, dir, UseHFCd) */
  /* where t is vector of approximate times (in ET, or ephemeris time format) */
  /* around which we will try to find density text files. The time vector t */
  /* should have a time step of roughly 1 day, because it will only look at the */
  /* YYYY/MM/DD part of the file names. Also, the default */
  /* data directory is assumed to be /data0/tmp/NEWdata/GITM/, i.e. the standard */
  /* location on impact01 */
  /* UseHFCd determines whether or not to use varying Cd. Assumed to be off. */
  /* If dir is specified, then searches dir for GITM data files */

  
  DIR *dir = NULL;
  struct dirent *ent = NULL;
  char *filename = NULL;
  char year[5], month[3], day[3], hour[3], minute[3], second[3];

  char utc[20] = "2000-01-01T00:00:00\0";
  //char *utc0 = NULL;
  SpiceDouble et;

  //int i;
  int counter=0;

  if ((dir = opendir(dirname)) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir(dir)) != NULL) {

      /* Assign filename to string */
      filename = ent->d_name;

      if(strcmp(filename, "..")==0) {
	//Do nothing
      } else if(strcmp(filename, ".")==0) {
	//Do nothing
      } else if(strcmp(filename, "temp")==0) {
	//Do nothing
      } else {
	
	/* Assign year, month, day, hour, minute, and second from filename */

	memcpy(year,   filename+3,  4*sizeof(char));
	memcpy(month,  filename+7,  2*sizeof(char));
	memcpy(day,    filename+9,  2*sizeof(char));
	memcpy(hour,   filename+11, 2*sizeof(char));
	memcpy(minute, filename+13, 2*sizeof(char));
	memcpy(second, filename+15, 2*sizeof(char));

	year[4] = '\0';
	month[2] = '\0';
	day[2] = '\0';
	hour[2] = '\0';
	minute[2] = '\0';
	second[2] = '\0';
	
	/* Reformat date to CSPICE compatible version */
	strcpy(utc, year);
	strcat(utc, "-");
	strcat(utc, month);
	strcat(utc, "-");
	strcat(utc, day);
	strcat(utc, "T");
	strcat(utc, hour);
	strcat(utc, ":");
	strcat(utc, minute);
	strcat(utc, ":");
	strcat(utc, second);

	printf("utc = %s\n", utc);
	
	/* Convert UTC to ET */
	str2et_c(utc, &et);
	
	/* Check if data file lies within simulation time */
	if(et >= tsearch[0] && et <= tsearch[1]) {	  
	  strcpy(files[counter], filename);
	  etlist[counter] = et;
	  counter++;

	  if(counter >= MAXFILES) {
	    printf("Maximum number of files exceeded! Exiting!\n");
	    exit(1);
	  }	 
	}

      }

    } /* End while loop over directory */

    closedir(dir);

  } else { /* Could not open directory */
      
    printf("Could not open directory! Missing directory: %s\n", dirname);
    exit(1);

  } /* End "if/else" for directory existence */

  if(files == NULL) {
    printf("Did not find any data files. Is the specified directory correct?\n");
    exit(1);
  }

  *Nfiles = counter;

}

void importgitm(char *dirname, double etlist[MAXFILES], char files[MAXFILES][100], double AltKm[NVERT], double LatDeg[NLAT], double LonDeg[NLON], double Rho[MAXFILES][NLON][NLAT][NVERT], double ndenvec[NSPECIES][MAXFILES][NLON][NLAT][NVERT], double Temperature[MAXFILES][NLON][NLAT][NVERT], int Nfiles, int *natm, int isat) {

  /* GITM data includes 4 ghost cells in each dimension */
  double rho_temp[NLON+4][NLAT+4][NVERT+4];
  double alt_temp[NVERT+4], lat_temp[NLAT+4], lon_temp[NLON+4];

  double temp_temp[NLON+4][NLAT+4][NVERT+4];
  double O1D[NLON][NLAT][NVERT], O4plus[NLON][NLAT][NVERT];
  double N2D[NLON][NLAT][NVERT], N2P[NLON][NLAT][NVERT], N4S[NLON][NLAT][NVERT];
  double O2[NLON][NLAT][NVERT];
  double N2[NLON][NLAT][NVERT];
  double He[NLON][NLAT][NVERT];
  double H[NLON][NLAT][NVERT];
  
  char Temperaturepath[100];
  char O1Dpath[100];
  char O4pluspath[100];
  char O2path[100];
  char N2Dpath[100];
  char N2Ppath[100];
  char N4Spath[100];
  char N2path[100];
  char Hepath[100];
  char Hpath[100];

  int i, j, k;
  int ifile;

  char HDF5file[256];
  char Altitudepath[100];
  char Latitudepath[100];
  char Longitudepath[100];
  char Densitypath[100];

  hid_t file_id, dataset_id;
  herr_t status;

  sprintf(Altitudepath,    "/run_%d/Altitude",       natm[isat]);
  sprintf(Latitudepath,    "/run_%d/Latitude",       natm[isat]);
  sprintf(Longitudepath,   "/run_%d/Longitude",      natm[isat]);
  sprintf(Densitypath,     "/run_%d/Density",        natm[isat]);

  /* Import the data from GITM and return as a 3x3x3 matrix where each element is the total density (kg/km^3). */
  /* Dimension order is latitude (rad), then longitude (rad), and altitude (km). */

  for(ifile=0; ifile<Nfiles; ifile++) {
    
    strcpy(HDF5file, dirname);
    strcat(HDF5file, files[ifile]);
    
    /* Open the HDF5 file */
    file_id = H5Fopen(HDF5file, H5F_ACC_RDWR, H5P_DEFAULT);
    
    /* Open and Read Altitude */
    dataset_id = H5Dopen2(file_id, Altitudepath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, alt_temp);

    /* Open and Read Latitude */
    dataset_id = H5Dopen2(file_id, Latitudepath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lat_temp);

    /* Open and Read Longitude */
    dataset_id = H5Dopen2(file_id, Longitudepath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lon_temp);

    /* Open and Read Density */
    dataset_id = H5Dopen2(file_id, Densitypath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_temp);

    //if (UseHFCd) {

    sprintf(Temperaturepath, "/run_%d/Temperature",    natm[isat]);
    sprintf(O1Dpath,         "/run_%d/O1D",            natm[isat]);
    sprintf(O4pluspath,      "/run_%d/O4plus",         natm[isat]);
    sprintf(O2path,          "/run_%d/O2",             natm[isat]);
    sprintf(N2Dpath,         "/run_%d/N2D",            natm[isat]);
    sprintf(N2Ppath,         "/run_%d/N2P",            natm[isat]);
    sprintf(N4Spath,         "/run_%d/N4S",            natm[isat]);
    sprintf(N2path,          "/run_%d/N2",             natm[isat]);
    sprintf(Hepath,          "/run_%d/He",             natm[isat]);
    sprintf(Hpath,           "/run_%d/H",              natm[isat]);
    
    /* Open and Read Temperature */
    dataset_id = H5Dopen2(file_id, Temperaturepath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_temp);
    
    /* Open and Read Number Density */
    dataset_id = H5Dopen2(file_id, O1Dpath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, O1D);
    
    dataset_id = H5Dopen2(file_id, O4pluspath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, O4plus);
    
    dataset_id = H5Dopen2(file_id, O2path, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, O2);
    
    dataset_id = H5Dopen2(file_id, N2Dpath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, N2D);
    
    dataset_id = H5Dopen2(file_id, N2Ppath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, N2P);
    
    dataset_id = H5Dopen2(file_id, N4Spath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, N4S);
    
    dataset_id = H5Dopen2(file_id, N2path, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, N2);
    
    dataset_id = H5Dopen2(file_id, Hepath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, He);
    
    dataset_id = H5Dopen2(file_id, Hpath, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, H);
    
    /* Sum different electronic states to find total number density for O and N */
    for(i=0; i<NLON; i++) {
      for(j=0; j<NLAT; j++) {
	for(k=0; k<NVERT; k++) {
	  Temperature[ifile][i][j][k] = temp_temp[i+2][j+2][k+2];
	  ndenvec[0][ifile][i][j][k] = O1D[i+2][j+2][k+2] + O4plus[i+2][j+2][k+2];
	  ndenvec[1][ifile][i][j][k] = O2[i+2][j+2][k+2];
	  ndenvec[2][ifile][i][j][k] = N2D[i+2][j+2][k+2] + N2P[i+2][j+2][k+2] + N4S[i+2][j+2][k+2];
	  ndenvec[3][ifile][i][j][k] = N2[i+2][j+2][k+2];
	  ndenvec[4][ifile][i][j][k] = He[i+2][j+2][k+2];
	  ndenvec[5][ifile][i][j][k] = H[i+2][j+2][k+2]; 
	}
      }
    }
    //} /* UseHFCd */


    /* Assign temporary data to full data structures */
    for(i=0; i<NVERT; i++) {
      AltKm[i] = alt_temp[i+2]/1.0e3;
    }

    for(i=0; i<NLAT; i++) {
      LatDeg[i] = lat_temp[i+2]*RAD2DEG;
    }

    for(i=0; i<NLON; i++) {
      LonDeg[i] = lon_temp[i+2]*RAD2DEG;
    }

    for(i=0; i<NLON; i++) {
      for(j=0; j<NLAT; j++) {
	for(k=0; k<NVERT; k++) {
	  Rho[ifile][i][j][k] = rho_temp[i+2][j+2][k+2]*1.0e9;
	}
      }
    }

    /* Close the HDF5 data set */
    status = H5Dclose(dataset_id);

    /* Close the HDF5 file */
    status = H5Fclose(file_id);

  }   

}


void importmsis(char *dirname, double etlist[MAXFILES], char files[MAXFILES][100], double AltKm[NVERT], double LatDeg[NLAT], double LonDeg[NLON], double Rho[MAXFILES][NLON][NLAT][NVERT], double ndenvec[NSPECIES][MAXFILES][NLON][NLAT][NVERT], double Temperature[MAXFILES][NLON][NLAT][NVERT], int Nfiles) {

  char *rho_temp;
  char *temp_temp;
  char *alt_temp, *lat_temp, *lon_temp;
  char *O2, *O, *N2, *N, *He, *H;
  char line1[1024];

  int i, j, k;//, l;
  int ifile, iline;
  int numlines;

  char MSISfile[256];

  /* Define dimensional arrays */
  for(i=0; i<NLON; i++) {
    LonDeg[i] = 5.0*i;
  }

  for(i=0; i<NLAT; i++) {
    LatDeg[i] = -90 + 5.0*i;
  }

  for(i=0; i<NVERT; i++) {
    AltKm[i] = 200.0 + 5.0*i;
  }

  for(ifile=0; ifile<Nfiles; ifile++) {
    
    strcpy(MSISfile, dirname);
    strcat(MSISfile, "/");
    strcat(MSISfile, files[ifile]);

    /* Open MSIS data file */
    FILE *f = fopen(MSISfile, "r");
    if (f == NULL) {
      printf("%s: %s\n", strerror(errno), MSISfile);
      exit(1);
    }

    printf("file = %s\n", MSISfile);

    /* Read Header Line */
    fgets(line1, 1024, f);
    fgets(line1, 1024, f);

    /* Determine number of lines in file */
    numlines = NLAT*NLON*NVERT;

    /* Loop over lines */
    for(iline=0; iline<numlines; iline++) {
      
      fgets(line1, 1024, f);

      alt_temp = strtok(line1, " ");
      lon_temp = strtok(NULL, " ");
      lat_temp = strtok(NULL, " ");
      rho_temp = strtok(NULL, " ");
      temp_temp = strtok(NULL, " ");
      O =  strtok(NULL, " ");
      O2 = strtok(NULL, " ");
      N =  strtok(NULL, " ");
      N2 = strtok(NULL, " ");
      He = strtok(NULL, " ");
      H =  strtok(NULL, " ");
    
      k = iline / (NLON*NLAT);
      j = (iline / NLON) % NLAT;
      i = iline % NLON;

      Rho[ifile][i][j][k] = atof(rho_temp)*1.0e9;

      if (UseHFCd) {
	Temperature[ifile][i][j][k] = atof(temp_temp);
	ndenvec[0][ifile][i][j][k] = atof(O);
	ndenvec[1][ifile][i][j][k] = atof(O2);
	ndenvec[2][ifile][i][j][k] = atof(N);
	ndenvec[3][ifile][i][j][k] = atof(N2);
	ndenvec[4][ifile][i][j][k] = atof(He);
	ndenvec[5][ifile][i][j][k] = atof(H);
      } /* UseHFCd */

    }

    fclose(f);    

  }  

}


void read_earth_grav(const char *data_dir, double (*CMAT)[MAXGDO], double (*SMAT)[MAXGDO]) {

  /* Purpose:  Read in Gravity Coefficients */
  /* Written by T. Kelecy, Boeing-LTS, 5/16/2012 */
  /* Converted to C - A. Walker, 08/02/2013  */

  /* Load the file */
  char filepath[PATH_MAX];
  char *filename = "egm96_to360_ascii.txt";
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, filename);
  FILE *gravfile = fopen(filepath, "r");
  if (gravfile == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  char line1[1024];
  int num_lines = 0;
  int iline, rowinc;
  int i, j;
  
  double *degree = NULL;
  double *order = NULL;
  double *C1 = NULL;
  double *C2 = NULL;
  double *C3 = NULL;
  double *C4 = NULL;
 
  if (EarthGrav == 1) { /* EGM96 */

    num_lines = read_num_lines(data_dir, filename);

    degree = (double *) calloc(num_lines, sizeof(double));
    order = (double *) calloc(num_lines, sizeof(double));
  
    C1 = (double *) calloc(num_lines, sizeof(double));
    C2 = (double *) calloc(num_lines, sizeof(double));
    C3 = (double *) calloc(num_lines, sizeof(double));
    C4 = (double *) calloc(num_lines, sizeof(double));

    for(iline=0; iline<num_lines; iline++) {
    
      /* Read the first line */
      fgets(line1, 1024, gravfile);

      /* Assign the variables with the pieces of the line */
      sscanf(line1,"%lf %lf %lf %lf %lf %lf\n",&degree[iline],&order[iline],&C1[iline],&C2[iline],&C3[iline],&C4[iline]);

    }
  
    if(GDO > 360) {
      printf("Degree & Order of EGM96 Gravity Field is too high! Maximum is 360.\n");
      exit(1);
    }

    if(GDO >= 2) {
      rowinc = 0;
      for(i=2; i<GDO+1; i++) {
	for(j=0; j<i+1; j++) {
	  CMAT[i][j] = C1[rowinc];
	  SMAT[i][j] = C2[rowinc];
	  rowinc += 1;
	}
      }
    }
 
  } else { /* JGM-3 */

    if(GDO > 8) {
      printf("Degree & Order of JGM-3 Gravity Field is too high! Maximum is 8.\n");
      exit(1);
    }

    CMAT[2][0] = -0.00108263602298400;
    SMAT[2][0] = 0.0;
    CMAT[3][0] = 2.53243534575440e-06;
    SMAT[3][0] = 0.0;
    CMAT[4][0] = 1.61933120507190e-06;
    SMAT[4][0] = 0.0;
    CMAT[5][0] = 2.27716101636880e-07;
    SMAT[5][0] = 0.0;
    CMAT[6][0] = -5.39648490498340e-07;
    SMAT[6][0] = 0.0;
    CMAT[7][0] = 3.51368442103180e-07;
    SMAT[7][0] = 0.0;
    CMAT[8][0] = 2.02518715208850e-07;
    SMAT[8][0] = 0.0;
    CMAT[2][1] = -2.41400005222210e-10;
    SMAT[2][1] = 1.54309997378440e-09;
    CMAT[3][1] = 2.19279880189650e-06;
    SMAT[3][1] = 2.68011893797260e-07;
    CMAT[4][1] = -5.08725303650240e-07;
    SMAT[4][1] = -4.49459935081170e-07;
    CMAT[5][1] = -5.37165101876620e-08;
    SMAT[5][1] = -8.06634638285300e-08;
    CMAT[6][1] = -5.98779768563030e-08;
    SMAT[6][1] = 2.11646643543820e-08;
    CMAT[7][1] = 2.05148727976720e-07;
    SMAT[7][1] = 6.93698935259080e-08;
    CMAT[8][1] = 1.60345871413790e-08;
    SMAT[8][1] = 4.01997815995100e-08;
    CMAT[2][2] = 1.57453604276720e-06;
    SMAT[2][2] = -9.03868073018690e-07;
    CMAT[3][2] = 3.09016044555830e-07;
    SMAT[3][2] = -2.11402397859750e-07;
    CMAT[4][2] = 7.84122307523660e-08;
    SMAT[4][2] = 1.48155456947140e-07;
    CMAT[5][2] = 1.05590535386740e-07;
    SMAT[5][2] = -5.23267239876320e-08;
    CMAT[6][2] = 6.01209884373730e-09;
    SMAT[6][2] = -4.65039481322170e-08;
    CMAT[7][2] = 3.28449048364920e-08;
    SMAT[7][2] = 9.28231438850840e-09;
    CMAT[8][2] = 6.57654233167430e-09;
    SMAT[8][2] = 5.38131640550560e-09;
    CMAT[3][3] = 1.00558857414550e-07;
    SMAT[3][3] = 1.97201323898890e-07;
    CMAT[4][3] = 5.92157432140720e-08;
    SMAT[4][3] = -1.20112918313970e-08;
    CMAT[5][3] = -1.49261538673890e-08;
    SMAT[5][3] = -7.10087714069860e-09;
    CMAT[6][3] = 1.18226641159150e-09;
    SMAT[6][3] = 1.84313368806250e-10;
    CMAT[7][3] = 3.52854051915120e-09;
    SMAT[7][3] = -3.06115023827880e-09;
    CMAT[8][3] = -1.94635815553990e-10;
    SMAT[8][3] = -8.72351950476050e-10;
    CMAT[4][4] = -3.98239574041290e-09;
    SMAT[4][4] = 6.52560581133960e-09;
    CMAT[5][4] = -2.29791235026810e-09;
    SMAT[5][4] = 3.87300507708040e-10;
    CMAT[6][4] = -3.26413891178910e-10;
    SMAT[6][4] = -1.78449133488820e-09;
    CMAT[7][4] = -5.85119491486240e-10;
    SMAT[7][4] = -2.63618221578670e-10;
    CMAT[8][4] = -3.18935802118560e-10;
    SMAT[8][4] = 9.11773558872550e-11;
    CMAT[5][5] = 4.30476750450290e-10;
    SMAT[5][5] = -1.64820394686360e-09;
    CMAT[6][5] = -2.15577115139000e-10;
    SMAT[6][5] = -4.32918169895400e-10;
    CMAT[7][5] = 5.81848560308730e-13;
    SMAT[7][5] = 6.39725266392350e-12;
    CMAT[8][5] = -4.61517343066280e-12;
    SMAT[8][5] = 1.61252083467840e-11;
    CMAT[6][6] = 2.21369255567410e-12;
    SMAT[6][6] = -5.52771222059660e-11;
    CMAT[7][6] = -2.49071768205960e-11;
    SMAT[7][6] = 1.05348786292660e-11;
    CMAT[8][6] = -1.83936426976340e-12;
    SMAT[8][6] = 8.62774316741500e-12;
    CMAT[7][7] = 2.55907801498730e-14;
    SMAT[7][7] = 4.47598341447510e-13;
    CMAT[8][7] = 3.42976181846240e-13;
    SMAT[8][7] = 3.81476566866850e-13;
    CMAT[8][7] = -1.58033228917250e-13;
    SMAT[8][8] = 1.53533813971480e-13;
  
  } /* EarthGrav */
   
  free(degree);
  free(order);
  free(C1);
  free(C2);
  free(C3);
  free(C4);

  fclose(gravfile);

}


/******************************** READ NUMBER OF LINES ***************************/
int read_num_lines(const char *data_dir, char *filename) {

  /* Read in the EGM96 coefficient file and determine number of lines*/
  /* Input:  */
  /* Output: Number of Facets */

  int num_lines = 0;
  int ch;

  char filepath[PATH_MAX];
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, filename);
  FILE *f = fopen(filepath, "r");
  
  /*Check that file exists*/
  if(!f) {
    printf("Gravity coefficient file does not exist\n");
    exit(1);
  }
  
  while (EOF != (ch=fgetc(f))) 
    if (ch=='\n')
      num_lines++;

  fclose(f);

  return(num_lines);

}


/******************************** READ INPUT FILE ********************************/
void read_input_file(const char *data_dir, double *a0, double *e0, double *i0, double *Om0, double *om0, double *f0, double q0[4], double omega0[3], char **utc0, int *maxt, double *deltat, int *scn, int *n,  double *sigma, double *mean, double *psig, double *vsig, double *cf107, double *cf107A, double *cAp, struct initCd_struct *initCd, char mesh_filename[100]) {

  /* Read in the propagator input file and assign variables */
  /* Input: prop.inp */
  /* Output: Initial Keplerian elements (a0, e0, i0, Om0, om0, f0)
     Initial time string (utc0)
     Total simulation time (maxt)
     Timestep (deltat)
     Mean and sigma for ballistic coefficient scale factor (mean and sigma)
     Position and velocity 1-sigma variance on satellite distribution (psig and vsig) */

  char filepath[PATH_MAX];
  char *filename = "prop.inp";
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, filename);

  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  char line[1024];
  char *temp;
  char *data;
  int i;

  /* READ IN HEADER AND EMPTY LINES */
  for(i=0; i<4; i++) {
    fgets(line, 1024, f);
  }

  /* READ IN INITIAL STATE VECTOR IN KEPLERIAN ELEMENTS */
  for(i=0; i<6; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* SEMI-MAJOR AXIS (KM) */
      sscanf(data, "%lf\n", a0);
    if(i==1) /* ECCENTRICITY (UNITLESS) */
      sscanf(data, "%lf\n", e0);
    if(i==2) /* INCLINATION (DEGREES) */
      sscanf(data, "%lf\n", i0);
    if(i==3) /* ARGUMENT OF PERIGEE (DEGREES) */
      sscanf(data, "%lf\n", om0);
    if(i==4) /* RIGHT-ASCENSION OF THE ASCENDING NODE (DEGREES) */
      sscanf(data, "%lf\n", Om0);
    if(i==5) /* TRUE ANOMALY (DEGREES) */
      sscanf(data, "%lf\n", f0);

  }

  /* READ IN INITIAL ROTATIONAL STATE */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  for(i=0; i<7; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");
    
    if(i==0) /* QUATERNION (VECTOR 1) */
      sscanf(data, "%lf\n", &q0[0]);
    if(i==1) /* QUATERNION (VECTOR 2) */
      sscanf(data, "%lf\n", &q0[1]);
    if(i==2) /* QUATERNION (VECTOR 3) */
      sscanf(data, "%lf\n", &q0[2]);
    if(i==3) /* QUATERNION (SCALAR 4) */
      sscanf(data, "%lf\n", &q0[3]);
    if(i==4) /* X-AXIS ANGULAR VELOCITY (RAD/S) */
      sscanf(data, "%lf\n", &omega0[0]);
    if(i==5) /* Y-AXIS ANGULAR VELOCITY (RAD/S) */
      sscanf(data, "%lf\n", &omega0[1]);
    if(i==6) /* Z-AXIS ANGULAR VELOCITY (RAD/S) */
      sscanf(data, "%lf\n", &omega0[2]);
    
  }


  /* READ IN TIME CONTROLS */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  fgets(line, 1024, f);
  temp = strtok(line, "#");
  data = strtok(NULL, "#");
  *utc0 = strdup(data);
  (*utc0)[21] = '\0';
  
  for(i=0; i<2; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* TOTAL SIMULATION TIME (SECONDS) */
      sscanf(data, "%d\n", maxt);
    if(i==1) /* TIMESTEP (SECONDS) */
      sscanf(data, "%lf\n", deltat);

  }

  /* READ IN SATELLITE CONTROLS */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  for(i=0; i<6; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* SATELLITE CATALOG NUMBER */
      sscanf(data, "%d\n", scn);
    if(i==1) /* NUMBER OF SATELLITES */
      sscanf(data, "%d\n", n);
    if(i==2) /* BALLISTIC COEFFICIENT SCALE FACTOR STANDARD DEVIATION */
      sscanf(data, "%lf\n", sigma);
    if(i==3) /* BALLISTIC COEFFICIENT SCALE FACTOR MEAN */
      sscanf(data, "%lf\n", mean);
    if(i==4) /* SATELLITE POSITION DISTRIBUTION 1-SIGMA VARIANCE */
      sscanf(data, "%lf\n", psig);
    if(i==5) /* SATELLITE VELOCITY DISTRIBUTION 1-SIGMA VARIANCE */
      sscanf(data, "%lf\n", vsig);

  }

  /* READ IN MSIS INPUT CONTROLS */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* DAILY OBSERVED F10.7 SOLAR FLUX PROXY */
      sscanf(data, "%lf\n", cf107);
    if(i==1) /* 81-DAY CENTERED AVERAGE F10.7 SOLAR FLUX PROXY */
      sscanf(data, "%lf\n", cf107A);
    if(i==2) /* DAILY AVERAGE GEOMAGENTIC INDEX */
      sscanf(data, "%lf\n", cAp);

  }

  /* READ IN DRAG COEFFICIENT CONTROLS */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  for(i=0; i<11; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* SATELLITE GEOMETRY */
      sscanf(data, "%d\n", &initCd->geometry);
    if(i==1) /* GAS-SURFACE INTERACTION MODEL */
      sscanf(data, "%d\n", &initCd->gsi_model);
    if(i==2) /* ADSORPTION MODEL */
      sscanf(data, "%d\n", &initCd->ads_model);
    if(i==3) /* PITCH ANGLE */
      sscanf(data, "%lf\n", &initCd->pitch);
    if(i==4) /* YAW ANGLE */
      sscanf(data, "%lf\n", &initCd->yaw);
    if(i==5) /* CYLINDER RADIUS */
      sscanf(data, "%lf\n", &initCd->radius);
    if(i==6) /* CUBOID OR CYLINDER LENGTH */
      sscanf(data, "%lf\n", &initCd->length);
    if(i==7) /* CUBOID HEIGHT */
      sscanf(data, "%lf\n", &initCd->height);
    if(i==8) /* CUBOID WIDTH */
      sscanf(data, "%lf\n", &initCd->width);
    if(i==9) /* SURFACE MATERIAL PARTICLE MASS */
      sscanf(data, "%lf\n", &initCd->surf_mass);
    if(i==10) /* MESH FILENAME */
      sscanf(data, "%s\n", mesh_filename);

  }

  fclose(f);

}


void importussa76(const char *data_dir, double *altitude, double *density) {

  /* Read in 1976 U.S. Standard Atmosphere Look-up Table*/

  /* Inputs: Standard Atmosphere Altitude Array - sa_altitude */
  /*         Standard Atmosphere Density Array - sa_density */

  /* Outputs: Standard Atmosphere Look-up Table */

  /* Load the file */
  char *filename = "USSA1976.dat";
  char filepath[PATH_MAX];
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, filename);
  FILE *safile = fopen(filepath, "r");
  if (safile == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  char line1[1024];
  int iline;
  //int rowinc;
  //int i, j;

  /* Read header line */
  fgets(line1, 1024, safile);

  for(iline=0; iline<NVERT; iline++) {
    
    /* Read the first line */
    fgets(line1, 1024, safile);

    /* Assign the variables with the pieces of the line */
    sscanf(line1,"%lf %lf\n",&altitude[iline], &density[iline]);

  }

  fclose(safile);

}


void create_constants_structure(struct constants_struct *psatconstant, struct initCd_struct *initCd) {
  
  struct constants_struct *satconstant;

  /* Assign constants to constants structure */
  satconstant = psatconstant + 0;
 
  satconstant->DRAGFLAG = UseDrag;
  satconstant->SPHFLAG = UseSphHarmon;            
  satconstant->DOGF = GDO;
  satconstant->ATMDENMODEL = DensityModel;          
  satconstant->EARTHGRAVMODEL = EarthGrav;
  satconstant->SRPFLAG = UseSRP;
  satconstant->CDFLAG = UseHFCd;     
  satconstant->BCFLAG = UseBCscale;     
  satconstant->TAVGFLAG = UseSTimeAvg;  
  satconstant->MOONFLAG = UseMoon;
  satconstant->SUNFLAG = UseSun;   
  satconstant->RKORDER = RungeKutta;     
  satconstant->JGM3RE = JGM3Re;
  satconstant->JGM3RP = JGM3Rp;   
  satconstant->JGM3GME = JGM3GMe;   
  satconstant->JGM3OME = JGM3ome; 
  satconstant->EGM96GME = EGM96GMe; 
  satconstant->EGM96RE = EGM96Re;
  satconstant->MOONGM = MoonGM; 
  satconstant->SUNGM = SunGM;   
  satconstant->EARTHECC = EECC;     
  satconstant->GEOMETRY = initCd->geometry;   
  satconstant->ADSMODEL = initCd->ads_model;
  satconstant->GSIMODEL = initCd->gsi_model;          

}


/******************************** READ ENSEMBLE FILE ********************************/
void read_state_ensemble(const char *data_dir, int nsat, int scn, double **initstate, int *natm, 
			 double **msismod, double **orientation) {

  /* Read in the initial ensemble of satellite states from file and assign variables */
  /* Input: Satellite name.txt */
  /* Output: Initial cartesian elements (x, y, z, u, v, w) */
  /*         Atmospheric ensemble flag */

  int hs = 0; /* Number of header lines */
  char satfilename[100];
  sprintf(satfilename, "%d.txt", scn);


  char filepath[PATH_MAX];
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, satfilename);
  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }
  char line[1024];
  double ntemp;
  double mtemp[3];
  double attitude[2];
  int i, isat;

  /* READ IN FIRST LINE HEADER */
  for(i=0; i<hs; i++) {
    fgets(line, 1024, f);
  }

  /* READ IN INITIAL STATE VECTOR IN KEPLERIAN ELEMENTS */
  for(isat = 0; isat < nsat; isat++) {
   
    if(DensityModel==1 && DynamicMSIS) {

      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
	     &initstate[isat][0], &initstate[isat][1], &initstate[isat][2],
	     &initstate[isat][3], &initstate[isat][4], &initstate[isat][5],
	     &mtemp[0], &mtemp[1], &mtemp[2], &attitude[0], &attitude[1]);

      for(i=0; i<3; i++) {
	msismod[isat][i] = mtemp[i];
      }

      /* Assign ensemble member attitude to orientation array */
      for(i=0; i<2; i++) {
	orientation[isat][i] = attitude[i];
      }

    } else {

      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf\n", 
	     &initstate[isat][0], &initstate[isat][1], &initstate[isat][2],
	     &initstate[isat][3], &initstate[isat][4], &initstate[isat][5],
	     &ntemp);

      natm[isat] = (int)ntemp;

    }

  }
      
  fclose(f);

}


void RemoveSpaces(char* source)
{
  char* i = source;
  char* j = source;
  while(*j != 0)
  {
    *i = *j++;
    if(*i != ' ')
      i++;
  }
  *i = 0;
}
