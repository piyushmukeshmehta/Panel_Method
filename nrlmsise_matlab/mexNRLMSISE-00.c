/* ==========================================================================
 * mex NRLMSISE-00 atmospheric model
 *==========================================================================*/

#include "mex.h"
#include "string.h"
#include "matrix.h"
#include "nrlmsise-00.h"

#define MAXCHARS 80   /* max length of string contained in each field */

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    const char **fnames;      /* pointers to field names */
    const mwSize *dims;
    mxArray    *tmp, *fout;
    double       *pdata=NULL;
    int        ifield, nfields, i;
    mxClassID  *classIDflags;
  /*  mwIndex    jstruct; */
    mwSize     NStructElems;
    mwSize     ndim;
    double *tmp_pr;
    
    struct nrlmsise_input input;
    struct nrlmsise_flags flags;
    struct nrlmsise_output output;
    struct ap_array aph,*ap;
    
    /* ----------------------------------------------------------------- */
    /* mexPrintf("NRHS = %d.\n",  nrhs);
    mexPrintf("NLHS = %d.\n",  nlhs); */
            
//     /* check proper input and output */
//     if(nrhs <1){
//         mexErrMsgTxt("At least one input required.");
//     }else if(nrhs >3){
//         mexErrMsgTxt("At most two inputs required.");
//     }else if(nlhs > 1){
//         mexErrMsgTxt("Too many output arguments.");
//     }else if(!mxIsStruct(prhs[0]) && !mxIsStruct(prhs[1])){
//         mexErrMsgTxt("Input must be a structure.");
//     }
    
    /* ------------------------- INPUTS -------------------------------- */
    /* -------------- get first input arguments -------------- */
    nfields = mxGetNumberOfFields(prhs[0]);
    /* mexPrintf("nfileds = %d.\n",  nfields); */
//     if(nfields>12){
//         mexErrMsgTxt("At most three fields are ammissible.");
//     }
    NStructElems = mxGetNumberOfElements(prhs[0]);
    /* mexPrintf("NstructElems = %d.\n", NStructElems); */
//     if( NStructElems!=1){
//         mexErrMsgTxt("Only one data structure is ammissible.");
//     }
    
    /* check empty field, proper data type, and data type consistency;
     * and get classID for each field. */
//     for(ifield=0; ifield<nfields-1; ifield++) {
//         tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
//         if(tmp == NULL) {
//             mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1);
//             mexErrMsgIdAndTxt( "MATLAB:fieldEmpty",
//                     "Above field is empty!");
//         }
//         if(mxIsChar(tmp) ||
//                 mxIsComplex(tmp) || mxGetNumberOfElements(tmp)!=1){
//             mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1);
//             mexErrMsgIdAndTxt( "MATLAB:fieldNotRealScalar",
//                     "Numeric data in above field must be scalar and noncomplex!");
//         }
//     }
        
    /* -------------- get second input arguments -------------- */
    if(nrhs>1){
        nfields = mxGetNumberOfFields(prhs[1]);
        /* mexPrintf("nfileds = %d.\n",  nfields); */
//         if(nfields>3) {
//             mexErrMsgTxt("At most 3 fields are ammissible.");
//         }
        NStructElems = mxGetNumberOfElements(prhs[1]);
        /* mexPrintf("NstructElems = %d.\n", NStructElems); */
//         if( NStructElems!=1) {
//             mexErrMsgTxt("Only one data structure is ammissible.");
//         }
        
        /* check empty field, proper data type, and data type consistency;
         * and get classID for each field. */
//         for(ifield=0; ifield<nfields; ifield++) {
//             tmp = mxGetFieldByNumber(prhs[1], 0, ifield);
//             if(tmp == NULL) {
//                 mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1);
//                 mexErrMsgIdAndTxt( "MATLAB:fieldEmpty",
//                         "Above field is empty!");
//             }
//             if(mxIsChar(tmp) ||
//                     mxIsComplex(tmp) || mxGetNumberOfElements(tmp)>24){
//                 mexPrintf("%s%d\t%s%d\n", "FIELD: ", ifield+1);
//                 mexErrMsgIdAndTxt( "MATLAB:fieldNotRealScalar",
//                         "Numeric data in above field must be scalar and noncomplex!");
//             }
//         }
        /* assign value to struct nrlmsise_flags flags */
        /* ----------------- SWITCHES ---------------- */
        ifield = mxGetFieldNumber(prhs[1], "switches");
        if(ifield!=-1) {
            tmp = mxGetFieldByNumber(prhs[1], 0, ifield);
//             if ( mxGetN(tmp)*mxGetM(tmp)>24 ) {
//                 mexErrMsgTxt( "MATLAB: NotVector. Input R must be at most 24-vector.");
//             }
            nfields = mxGetNumberOfElements(tmp);
            /* mexPrintf("SWITCHES: nfields = %d.\n",  nfields); */
//             if(nfields>24){
//                 mexErrMsgTxt("At most 24 fields are ammissible.");
//             }
            tmp_pr = mxGetPr(tmp);
            for(i=0; i<nfields; i++) {
                flags.switches[i] = tmp_pr[i];
                /* mexPrintf("FLAGS.SWITCHES: %d\n",flags.switches[i]); */
            }
        }else{
            flags.switches[0] = 0;
            for(ifield=1; ifield<24; ifield++) {
                flags.switches[ifield] = 1;
            }
        }
    }else{
        flags.switches[0] = 0;
        for(ifield=1; ifield<24; ifield++) {
            flags.switches[ifield] = 1;
        }
    }
    /* ----------------- SW ---------------- */
    /* set internally */
    /* ----------------- SWC ---------------- */
    /* set internally */
    
    /* assign value to struct nrlmsise_input input */
    /* ------------- YEAR ------------------------ */
    ifield = mxGetFieldNumber(prhs[0], "year");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.year = mxGetScalar(tmp);
        /* mexPrintf("YEAR : %d\n",input.year); */
    }else{
        input.year = 0;    /* without effects */
    }
    /* ------------- DAY OF YEAR ------------------ */
    ifield = mxGetFieldNumber(prhs[0], "doy");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.doy = mxGetScalar(tmp);
        /* mexPrintf("DOY : %d\n",input.doy); */
    }else{
        input.doy = 172;
    }
    /* ----------------- SECONDS ------------------ */
    ifield = mxGetFieldNumber(prhs[0], "sec");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.sec = mxGetScalar(tmp);
        /* mexPrintf("SEC : %f\n",input.sec); */
    }else{
        input.sec = 29000;
    }
        /* ----------------- ALTITUDE ------------------ */
    ifield = mxGetFieldNumber(prhs[0], "alt");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.alt = mxGetScalar(tmp);
        /* mexPrintf("ALT : %f\n",input.alt); */
    }else{
        input.alt = 400;
    }
    /* ----------------- LATITUDE ------------------ */
    ifield = mxGetFieldNumber(prhs[0], "g_lat");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.g_lat = mxGetScalar(tmp);
        /* mexPrintf("LAT : %f\n",input.g_lat); */
    }else{
        input.g_lat = 60;
    }  
    /* ----------------- LONGITUDE ------------------ */
    ifield = mxGetFieldNumber(prhs[0], "g_long");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.g_long = mxGetScalar(tmp);
        /* mexPrintf("LONG : %f\n",input.g_long); */
    }else{
        input.g_long = -70;
    }
    /* ------------------- LST --------------------- */
    ifield = mxGetFieldNumber(prhs[0], "lst");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.lst = mxGetScalar(tmp);
        /* mexPrintf("LST : %f\n",input.lst); */
    }else{
        input.lst = 16;
    }
    /* ------------------- F107A -------------------- */
    ifield = mxGetFieldNumber(prhs[0], "f107A");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.f107A = mxGetScalar(tmp);
        /* mexPrintf("F107A : %f\n",input.f107A); */
    }else{
        input.f107A = 150;
    }
    /* ------------------- F107 --------------------- */
    ifield = mxGetFieldNumber(prhs[0], "f107");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        input.f107 = mxGetScalar(tmp);
        /* mexPrintf("F107 : %f\n",input.f107); */
    }else{
        input.f107 = 150;
    }
    /* ------------ MAGNETIC INDEX SCALAR ----------- */
    ifield = mxGetFieldNumber(prhs[0], "ap");
    if(ifield!=-1) {
        tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
        if(mxIsChar(tmp) ||
                mxIsComplex(tmp) || mxGetNumberOfElements(tmp)!=1){
            mexPrintf("FIELD: (AP)  %d\n", ifield+1);
            mexErrMsgIdAndTxt( "MATLAB:fieldNotRealScalar",
                    "Numeric data in above field must be scalar and noncomplex!");
        }
        input.ap = mxGetScalar(tmp);
        /* mexPrintf("AP : %f\n",input.ap); */
    }else{
        input.ap = 4;
    }
    /* ------------ MAGNETIC INDEX VECTOR ----------- */ 
    if(flags.switches[9] == -1) {
        ifield = mxGetFieldNumber(prhs[0], "ap_a");
        if(ifield!=-1){
            tmp = mxGetFieldByNumber(prhs[0], 0, ifield);
            nfields = mxGetNumberOfElements(tmp);
            /* mexPrintf("APH fields %d\n",nfields); */
            if(nfields!=7){
                mexErrMsgTxt("Field AP_A must be of lenght 7.");
            }
            tmp_pr = mxGetPr(tmp);
            for(i=0; i<nfields; i++) {
                aph.a[i] = tmp_pr[i];
                 /* mexPrintf("AP_A(%d) =  %f\n",i,aph.a[i]); */
            }
            input.ap_a = &aph;
            /* for(i=0; i<7; i++){
                mexPrintf("AP_A(%d) =  %f\n",i,aph.a[i]);
            } */
        }else{
            for(i=0; i<7; i++) {
                aph.a[i] = 100;
            }
            input.ap_a = &aph;
        }
    }
       
    /* ----------------------  C FUNCTION ------------------------------ */   
    /* call to the C function */
    gtd7(&input,&flags,&output);
    
    /* printf("\n");
    for (i=0;i<9;i++)
        printf("%E ",output.d[i]);
    printf("\n%E ",output.t[0]);
    printf("%E \n",output.t[1]); */
    
  /*  printf("\n\n");
    printf("\nTINF  ");
    printf("     %7.2f",output.t[0]);
    printf("\nTG    ");
    printf("     %7.2f",output.t[1]);
    printf("\nHE    "); 
    printf("   %1.3e",output.d[0]);
    printf("\nO     ");
    printf("   %1.3e",output.d[1]);
    printf("\nN2    ");
    printf("   %1.3e",output.d[2]);
    printf("\nO2    ");
    printf("   %1.3e",output.d[3]);
    printf("\nAR    ");
    printf("   %1.3e",output.d[4]);
    printf("\nH     ");
    printf("   %1.3e",output.d[6]);
    printf("\nN     ");
    printf("   %1.3e",output.d[7]);    
    printf("\nANM 0 ");
    printf("   %1.3e",output.d[8]);
    printf("\nRHO   ");
    printf("   %1.3e",output.d[5]);
    printf("\n");
*/
    
    /* -------------------------- OUTPUT ------------------------------- */
    
    nfields = 2;
    
    /* allocate memory  for storing pointers */
    fnames = mxCalloc(nfields, sizeof(*fnames));
    fnames[0] = "t";
    fnames[1] = "d";
        
    /* create a 1x1 struct matrix for output  */
    plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames);
    mxFree((void *)fnames);
    
    /* Create mxArray data structures to hold the data
    to be assigned for the structure. Then assign the field values */
    tmp  = mxCreateDoubleMatrix(1,2, mxREAL);
    tmp_pr    = mxGetPr(tmp);
    for(i = 0; i<2; i++)
        tmp_pr[i] = (double)output.t[i];
    mxSetFieldByNumber(plhs[0],0,0,tmp);
    tmp  = mxCreateDoubleMatrix(1,2, mxREAL);
    tmp_pr    = mxGetPr(tmp);

    tmp  = mxCreateDoubleMatrix(1,9, mxREAL);
    tmp_pr    = mxGetPr(tmp);
    for(i = 0; i<9; i++)
        tmp_pr[i] = (double)output.d[i];
    mxSetFieldByNumber(plhs[0],0,1,tmp);

    return;
}
