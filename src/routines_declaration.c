#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP _TruncatedNormal_cholperm(SEXP, SEXP, SEXP);
extern SEXP _TruncatedNormal_lnNpr(SEXP, SEXP, SEXP);
extern SEXP _TruncatedNormal_Phinv(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_TruncatedNormal_cholperm", (DL_FUNC) &_TruncatedNormal_cholperm, 3},
    {"_TruncatedNormal_lnNpr",    (DL_FUNC) &_TruncatedNormal_lnNpr,    3},
    {"_TruncatedNormal_Phinv",    (DL_FUNC) &_TruncatedNormal_Phinv,    3},
    {NULL, NULL, 0}
};

void R_init_TruncatedNormal(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
