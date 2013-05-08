
#ifndef __jacobi_h__
#define __jacobi_h__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Eigensolver based on jacobi rotation; this is a rip off implementation
 * from the book Numerical Recipes in C.
 */
extern int jacobi (float **a, int n, float d[], float **v);

#ifdef __cplusplus
}
#endif
#endif /* __jacobi_h__ */
