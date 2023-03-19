// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* --------------------------- */
/* .Fortran calls */
/* --------------------------- */
extern void F77_NAME(coefs)(int * p, double * phi, int * q, double * theta, double * d, int * m, double * cks);

static const R_FortranMethodDef FortranEntries[] = {   
    {"coefs",         (DL_FUNC) &F77_NAME(coefs),         7},    
    {NULL, NULL, 0}
};

void R_init_PPMiss(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(info, TRUE);
}
