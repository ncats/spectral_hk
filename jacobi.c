
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "jacobi.h"

#ifndef EPS
# define EPS 1e-20
#endif

#ifndef MAX_ITER
# define MAX_ITER 500
#endif

/**
 * sort eigenvalues and eigenvectors in ascending order
 */
static void eigen_sorter (double d[], double **v, int n)
{
  int i, j, k;
  double p;

  for (i = 0; i < n-1; ++i)
    {
      p = d[k = i];
      for (j = i; j < n; ++j)
        if (d[j] < p)
          p = d[k = j];

      if (k != i)
        {
          d[k] = d[i];
          d[i] = p;
          for (j = 0; j < n; ++j)
            {
              p = v[j][i];
              v[j][i] = v[j][k];
              v[j][k] = p;
            }
        }
    }
}

/**
 * eigensolver from the book Numerical Recipes in C, 1992
 */
int jacobi (double **a, int n, double d[], double **v)
{
  int j, iq, ip, i, err = 0;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z=0;

  b = malloc (sizeof (*b) * n);
  if (b == 0)
    {
      err = -1;
      goto done;
    }

  z = malloc (sizeof (*z) * n);
  if (z == 0)
    {
      err = -1;
      goto done;
    }

  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        v[i][j] = 0.;
      v[i][i] = 1.;
    }

  for (i = 0; i < n; ++i)
    {
      b[i] = d[i] = a[i][i];
      z[i] = 0.;
    }

  for (i = 0; i < MAX_ITER; ++i)
    {
      sm = 0.;
      for (ip = 0; ip < n-1; ++ip)
        for (iq = ip+1; iq < n; ++iq)
          sm += fabs (a[ip][iq]);

      if (sm == 0.)
        goto done;

      if (i < 4)
        tresh = .2 * sm/(n*n);
      else
        tresh = 0.;
      
      for (ip = 0; ip < n-1; ++ip)
        {
          for (iq = ip+1; iq < n; ++iq)
            {
              g = 100.*fabs (a[ip][iq]);
              if (i > 4 && g <= EPS *fabs (d[ip]) && g <= EPS * fabs (d[iq]))
                a[ip][iq] = 0.;
              else if (fabs (a[ip][iq]) > tresh)
                {
                  h = d[iq] - d[ip];
                  if (g <= EPS * fabs (h))
                    t = a[ip][iq] / h;
                  else
                    {
                      theta = .5*h/a[ip][iq];
                      t = 1./(fabs (theta) + sqrt (1.+theta*theta));
                      if (theta < 0.) t = -t;
                    }
                  c = 1./sqrt(1+t*t);
                  s = t*c;
                  tau = s / (1.+c);
                  h = t*a[ip][iq];
                  z[ip] -= h;
                  z[iq] += h;
                  d[ip] -= h;
                  d[iq] += h;
                  a[ip][iq] = 0.;

#define __rotate(a,i,j,k,l) \
  do { \
    double g = a[i][j]; \
    double h = a[k][l]; \
    a[i][j] = g-s*(h+g*tau); \
    a[k][l] = h+s*(g-h*tau); \
  } while (0)

                  for (j = 0; j < ip; ++j)
                    __rotate (a, j, ip, j, iq);
                  for (j = ip+1; j < iq; ++j)
                    __rotate (a, ip, j, j, iq);
                  for (j = iq+1; j < n; ++j)
                    __rotate (a, ip, j, iq, j);
                  for (j = 0; j < n; ++j)
                    __rotate (v, j, ip, j, iq);
                }
            }
        }

      for (ip = 0; ip < n; ++ip)
        {
          b[ip] += z[ip];
          d[ip] = b[ip];
          z[ip] = 0.;
        }
    }
  err = 1; /* didn't converge; too many iterations */

 done:
  if (err >= 0)
    eigen_sorter (d, v, n);

  if (b != 0)
    free (b);
  if (z != 0)
    free (z);

  return err;
}


#ifdef __JACOBI_TEST
int main ()
{
  const int n = 16;
  /** 
   * expected output:

Eigenvalues & Fiedler vector
    0.0000     0.2407
    0.0952    -0.4151
    0.4147     0.2506
    0.6157     0.2079
    1.0555    -0.3954
    1.2361    -0.3954
    1.7745    -0.1910
    1.8350     0.0479
    2.0000    -0.0259
    2.7540     0.2366
    2.9476     0.1553
    3.5116     0.1752
    3.8308     0.2615
    4.2331     0.0679
    4.6543    -0.3380
    5.0418     0.1172
  */

  double m[16][16] = {
    {2, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 2, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-1, 0, 2, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
    {-1, 0, 0, 2, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
    {0, -1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
    {0, -1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
    {0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 0, 0, -1, 0},
    {0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, -1},
    {0, 0, 0, 0, 0, 0, -1, -1, 3, 0, 0, 0, 0, -1, 0, 0},
    {0, 0, -1, 0, 0, 0, 0, 0, 0, 3, 0, -1, -1, 0, 0, 0},
    {0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 3, -1, 0, -1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 3, 0, 0, 0, -1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 2, 0, 0},
    {0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 3, 0},
    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 2}    
  };
  double d[16] = {0};
  double **v, **a;
  int err, i;

  a = malloc (sizeof (double *)*n);
  v = malloc (sizeof (double *)*n);
  for (i = 0; i < n; ++i)
    {
      a[i] = m[i];
      v[i] = malloc (n*sizeof (double));
    }

  err = jacobi (a, 16, d, v);
  if (err == 0)
    {
      printf ("Eigenvalues & Fiedler vector\n");
      for (i = 0; i < 16; ++i)
        printf ("%10.4f %10.4f\n", d[i], v[i][1]);
    }
  else if (err > 0)
    printf ("error: jacobi didn't converge!\n");

  for (i = 0; i < n; ++i)
    free (v[i]);
  free (a);
  free (v);

  return err;
}
#endif

/**
 * Local Variables:
 * compile-command: "gcc -Wall -g -o jacobi jacobi.c -D__JACOBI_TEST -lm"
 * End:
 */
