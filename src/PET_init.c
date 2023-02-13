// #include <R.h>
// #include <Rinternals.h>
// #include <stdlib.h> // for NULL
// #include <R_ext/Rdynload.h>

// /* FIXME: 
//    Check these declarations against the C/Fortran source code.
// */

// /* .C calls */
// extern void Hough1(void *, void *, void *);
// extern void Hough2(void *, void *, void *);
// extern void Hough3(void *, void *, void *);
// extern void Hough4(void *, void *, void *);
// extern void iradon(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// extern void it(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// extern void radonLI(void *, void *, void *);
// extern void radonNN(void *, void *, void *);
// extern void radonSINC(void *, void *, void *);

// // extern void iradon_smoothed_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// // extern void test_vec_C(double *, int *, double *);

// /* .Call calls */
// extern SEXP loadFile(SEXP, SEXP);
// extern SEXP writeFile(SEXP, SEXP, SEXP, SEXP);

// static const R_CMethodDef CEntries[] = {
//     {"Hough1",    (DL_FUNC) &Hough1,     3},
//     {"Hough2",    (DL_FUNC) &Hough2,     3},
//     {"Hough3",    (DL_FUNC) &Hough3,     3},
//     {"Hough4",    (DL_FUNC) &Hough4,     3},
//     {"iradon",    (DL_FUNC) &iradon,    14},
//     {"it",        (DL_FUNC) &it,        34},
//     {"radonLI",   (DL_FUNC) &radonLI,    3},
//     {"radonNN",   (DL_FUNC) &radonNN,    3},
//     {"radonSINC", (DL_FUNC) &radonSINC,  3},
//     {NULL, NULL, 0}
// };

// static const R_CallMethodDef CallEntries[] = {
//     {"loadFile",  (DL_FUNC) &loadFile,  2},
//     {"writeFile", (DL_FUNC) &writeFile, 4},
//     {NULL, NULL, 0}
// };

// void R_init_PET(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
//     R_useDynamicSymbols(dll, FALSE);
// }

