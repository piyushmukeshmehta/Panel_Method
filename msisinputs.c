#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>

#include "SpiceUsr.h"
#include "defs.h"
#include "structs.h"
#include "io.h"
#include "msisinputs.h"

/******************************* DETERMINE NUMBER OF SW LINES *****************************/
int read_sw_lines(const char *data_dir) {

  /* Read in the space weather indices for propagation with Dynamic MSIS */
  /* Input: SpaceWeather.dat */
  /* Output: Number of observations (e.g. length of array) */

  int hs = 19; /* Number of header lines */
  int fs = 239; /* Number of footer lines */

  char filepath[PATH_MAX];
  char *SWfilename;
  SWfilename = "SpaceWeather.dat";
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, SWfilename);

  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }
  char line[1024];
  int i;
  int num_lines = 0;
  int nobs;
  
  /* READ IN HEADER LINES */
  for(i=0; i<hs; i++) {
    fgets(line, 1024, f);
  }

  /* Calculate number of lines to read */
  num_lines = read_num_lines(data_dir, SWfilename);
  nobs = num_lines - hs - fs;

  return(nobs);

}


/********************************* READ SPACE WEATHER *****************************/
void read_space_weather(const char *data_dir, struct msis_struct *pmsis, int nobs, double cf107, double cf107A, double cAp) {

  /* Read in the space weather indices for propagation with Dynamic MSIS */
  /* Input: SpaceWeather.dat */
  /* Output: Daily observed F10.7 radio flux proxy (sfu)       */
  /*         81-day centered average of F10.7 (sfu)            */
  /*         Ap index array (3-hour geomagnetic indices + AVG) */

  struct msis_struct *msis_inputs;

  int hs = 19; /* Number of header lines */
  
  char filepath[PATH_MAX];
  char *SWfilename;
  SWfilename = "SpaceWeather.dat";
  strcpy(filepath, data_dir);
  strcat(filepath, "/");
  strcat(filepath, SWfilename);

  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }
  char line[1024];
  int i;

  int itmp;
  double dtmp;
  
  /* READ IN HEADER LINES */
  for(i=0; i<hs; i++) {
    fgets(line, 1024, f);
  }

  /* Loop over observation dates */
  for(i=0; i<nobs; i++) {

    msis_inputs = pmsis + i;

    if(ReadSW) {

    /* Note that MSIS uses the Observed F10.7 with the centered-81 day average */
    fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %d %lf %lf %lf %lf %lf\n", 
	   &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, &itmp, 
	   &msis_inputs->ap_index[0], &msis_inputs->ap_index[1], &msis_inputs->ap_index[2], &msis_inputs->ap_index[3], 
	   &msis_inputs->ap_index[4], &msis_inputs->ap_index[5], &msis_inputs->ap_index[6], &msis_inputs->ap_index[7], 
	   &msis_inputs->ap_index[8], &dtmp, &itmp, &itmp, &dtmp, &itmp, &dtmp, &dtmp, 
	   &msis_inputs->F107, &msis_inputs->F107A, &dtmp);

    } else {

      msis_inputs->F107 = cf107;
      msis_inputs->F107A = cf107A;
      for(i=0; i<9; i++) {
	msis_inputs->ap_index[i] = cAp;
      }

    }

  }

  fclose(f);

}
    
/********************************* FIND CORRECT MSIS TIME *****************************/
int find_msis_time(double year, double doy) {
  
  int index = 0;
  int yr0 = 1957;      /* First year of space weather data */
  int day0 = 274;      /* First day of space weather data */
  double dpy = 365.25; /* Days per year */
  int myear;

  if(ReadSW) {
    myear = (int)year - yr0;
    index =  (int)doy + (int)(dpy*myear+0.25) - day0;
  } else {
    index = 0; 
  }

  return(index);

}

/**************************** COMPUTE AP ARRAY INPUT FOR MSIS *****************************/
void get_msis_ap_array(struct msis_struct *pmsis, int index, double hour, double msis_ap_array[7]) {

  int aphour;
  double ap_temp[32]; /* Temporary ap index array */
  double sum6=0.0;
  double sum7=0.0;
  int i;

  struct msis_struct *day0_inputs, *day1_inputs, *day2_inputs, *day3_inputs; 
  day0_inputs = pmsis + index;        /* Current day structure */
  day1_inputs = pmsis + index - 1;    /* Previous day structure */
  day2_inputs = pmsis + index - 2;    /* 2 days prior structure */
  day3_inputs = pmsis + index - 3;    /* 3 days prior structure */

  /* Define temporary ap index array */
  for(i=0; i<8; i++) {
    ap_temp[i] =    day3_inputs->ap_index[i];
    ap_temp[i+8] =  day2_inputs->ap_index[i];
    ap_temp[i+16] = day1_inputs->ap_index[i];
    ap_temp[i+24] = day0_inputs->ap_index[i];
  }
  
  /* Find the 3-hour ap index for the "current" time */
  aphour = hour/3 + 24;

  msis_ap_array[0] = day0_inputs->ap_index[8];  /* Daily average */
  msis_ap_array[1] = ap_temp[aphour];           /* 3hr index for current time */
  msis_ap_array[2] = ap_temp[aphour-1];         /* 3hr index for 3 hrs before current time */
  msis_ap_array[3] = ap_temp[aphour-2];         /* 3hr index for 6 hrs before current time */
  msis_ap_array[4] = ap_temp[aphour-3];         /* 3hr index for 9 hrs before current time */

  /* 6th array element is the average of 8 3hr AP indices for 12hrs to 33 hrs prior to current time */
  for(i=(aphour-11); i<(aphour-3); i++) {
    sum6 += ap_temp[i];
  }

  msis_ap_array[5] = sum6/8.0;

  /* 7th array element is the average of 8 3hr AP indices for 36 to 57 hrs prior to current time */
  for(i=(aphour-19); i<(aphour-11); i++) {
    sum7 += ap_temp[i];
  }

  msis_ap_array[6] = sum7/8.0;

}
