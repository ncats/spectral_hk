
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "b32.h"
#include "spectral.h"
#include "inchi.h"

/*
 * update as appropriate
 */
#define __SPECTRAL_VERSION "v0.1"

#ifdef HAVE_GSL
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_version.h>
# define SPECTRAL_VERSION __SPECTRAL_VERSION " (GSL-" GSL_VERSION ")"
#elif defined(HAVE_MKL)
# include "mkl.h" /* Intel MKL library */
# define _XSTR(X) _STR(X)
# define _STR(X) #X
# define SPECTRAL_VERSION __SPECTRAL_VERSION \
  " (MKL-" _XSTR(__INTEL_MKL__) "." _XSTR(__INTEL_MKL_MINOR__) \
  "." _XSTR(__INTEL_MKL_UPDATE__) ")"
#else
#warning "**** Please consider using either the GSL or MKL eigensolver. \
They are orders of magnitude faster! The bundled implementation is only \
for completeness sake. ****"
# include "jacobi.h"
# define SPECTRAL_VERSION __SPECTRAL_VERSION " (Built-in Jacobi solver)"
#endif

#ifndef EPS
# define EPS 1e-20
#endif


/*
 * maximum graph size; for large graphs a more specialized eigensolver
 * (e.g., iterative ones) are more appropriate.
 */
#ifndef SPECTRAL_MAXG
# define SPECTRAL_MAXG 5000
#endif

#define __set_edge(i,j) G[(i)*size+(j)] = G[(j)*size+(i)] = 1
#define __get_edge(i,j) G[(i)*size+(j)]

/**
 * internal state of spectral_t
 */
struct __spectral_s 
{
  size_t bsize; /* buffer size of spectrum */
  float *spectrum; /* spectrum buffer */
  float *fiedler; /* fiedler vector */
  sha1_t *sha1; /* sha1 hash */
  inchi_t *inchi;
  unsigned char digest[20]; /* digest buffer */
  char hashkey[31]; /* 9(topology) + 10(connection) + 11(full) */
  char errmsg[BUFSIZ];
};


static void
digest_spectrum (spectral_t *sp, int size)
{
  unsigned char data[2];
  unsigned int uv;
  int i = 0, j;

  /* skip over all disconnected components */
  while (i < size && sp->spectrum[++i] < EPS)
    ;

#ifdef SPECTRAL_DEBUG
  printf ("spectral sequence:");
#endif

  for (j = i; j < size; ++j)
    {
      uv = (int)(sp->spectrum[j] / sp->spectrum[i] + 0.5);
#ifdef SPECTRAL_DEBUG
      printf (" %u", uv);
#endif
      data[0] = uv & 0xff;
      data[1] = (uv & 0xffff) >> 8;
      sha1_update (sp->sha1, data, sizeof (data));
    }

#ifdef SPECTRAL_DEBUG
  { double x;
  printf ("\neigenvalues:\n");
  for (j = i; j < size; ++j)
    {
      uv = (int)(sp->spectrum[j] / sp->spectrum[i] + 0.5);
      x = sp->spectrum[j] < 1. ? 1. - sp->spectrum[j] : sp->spectrum[j] - 1.;
      printf ("%3d: %.10f => %3u\n", j, sp->spectrum[j], uv);
    }
  }
#endif
}

void
spectral_adjacency_graph (double **M, const int *G, int nv, size_t size)
{
  int i, j, v;
  /*
   * adjacency 
   */
  for (i = 0; i < nv; ++i)
    {
      for (j = i+1; j < nv; ++j)
        {
          v = __get_edge (i+1, j+1);
          M[i][j] = v;
        }

      for (j = 0; j < i; ++j)
        M[i][j] = M[j][i];

      M[i][i] = 0.;
    }
}

void
spectral_laplacian_graph (double **M, const int *G, int nv, size_t size)
{
  int i, j, v, d;
  /* 
   * laplacian matrix
   */
  for (i = 0; i < nv; ++i)
    {
      d = 0;
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d += v;
            M[i][j] = v ? -1 : 0.;
          }
      M[i][i] = d;
    }
}

void 
spectral_signless_graph (double **M, const int *G, int nv, size_t size)
{
  int i, j, v, d;
  /* 
   * signless laplacian matrix
   */
  for (i = 0; i < nv; ++i)
    {
      d = 0;
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d += v;
            M[i][j] = v ? 1 : 0.;
          }
      M[i][i] = d;
    }
}

void
spectral_normalized_graph (double **M, const int *G, int nv, size_t size)
{
  int i, j, v, *d;

  /* 
   * normalized laplacian matrix
   */
  d = malloc (nv *sizeof (int));
  for (i = 0; i < nv; ++i)
    {
      d[i] = 0;
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d[i] += v;
          }
      M[i][i] = 1;
      for (j = 0; j < i; ++j)
        M[i][j] = M[j][i] = __get_edge (j+1, i+1) ? -1./sqrt (d[i]*d[j]) : 0.;
    }
  free (d);
}

#ifdef HAVE_GSL
static int
graph_spectrum (float *spectrum, float *fiedler, const inchi_t *g)
{
  int i, j, v, *d, nv = inchi_node_count (g);
  double x;
  const int *G = inchi_matrix_A (g);
  size_t size = inchi_matrix_size (g);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc (nv);
  gsl_matrix *A = gsl_matrix_alloc (nv, nv);
  gsl_matrix *V = gsl_matrix_alloc (nv, nv);
  gsl_vector *L = gsl_vector_alloc (nv);

  /* 
   * normalized laplacian matrix
   */
  d = malloc (nv *sizeof (int));
  for (i = 0; i < nv; ++i)
    {
      d[i] = 0;
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d[i] += v;
          }

      gsl_matrix_set (A, i, i, 1);
      for (j = 0; j < i; ++j)
        {
          x = __get_edge (j+1, i+1) ? -1./sqrt (d[i]*d[j]) : 0.;
          gsl_matrix_set (A, i, j, x);
          gsl_matrix_set (A, j, i, x);
        }
    }
  free (d);

#ifdef SPECTRAL_DEBUG
#if 0
  printf ("G = \n");
  for (i = 0; i < nv; ++i)
    {
      for (j = 0; j < nv; ++j)
        printf (" %-4.5f", gsl_matrix_get (A, i, j));
      printf ("\n");
    }
#else
  printf ("G = [");
  for (i = 0; i < nv; ++i)
    {
      for (j = 0; j < nv; ++j)
        printf (" %-4.5f", gsl_matrix_get (A, i, j));
      printf (";\n");
    }
  printf ("]\n");
#endif
#endif
  
  gsl_eigen_symmv (A, L, V, ws);
  gsl_eigen_symmv_sort (L, V, GSL_EIGEN_SORT_VAL_ASC);

  i = 0;
  while (i < nv && gsl_vector_get (L, ++i) < EPS)
    ;
#ifdef SPECTRAL_DEBUG  
  printf ("Eigenvector of the second smallest eigenvalue (%d: %.5f):\n",
          i, gsl_vector_get (L, i));
#endif
  
  { gsl_vector_view ev = gsl_matrix_column (V, i);
    for (i = 0; i < nv; ++i)
      {
#ifdef SPECTRAL_DEBUG   
        printf ("% 3d: % 11.10f\n", i, gsl_vector_get (&ev.vector, i));
#endif
        fiedler[i] = gsl_vector_get (&ev.vector, i);
      }
  }

  
  for (i = 0; i < nv; ++i)
    spectrum[i] = gsl_vector_get (L, i);
    
  gsl_vector_free (L);
  gsl_matrix_free (V);
  gsl_matrix_free (A);
  gsl_eigen_symmv_free (ws);

  return 0;
}

#elif defined(HAVE_MKL)

static int
graph_spectrum (float *spectrum, float *fiedler, const inchi_t *g)
{
  double *a, *d;
  const int *G = inchi_matrix_A (g);
  size_t size = inchi_matrix_size (g);
  int i, j, v, err = 0, nv = g->nv;

  a = malloc (sizeof (double)*nv*nv);
  d = malloc (nv *sizeof (double));

  /* normalized laplacian */
  for (i = 0; i < nv; ++i)
    {
      d[i] = 0.;
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d[i] += v;
          }

      a[i*nv+i] = 1;

      /* upper triangle */
      for (j = 0; j < i; ++j)
        {
          a[j*nv+i] = __get_edge (j+1, i+1) ? -1./sqrt (d[i]*d[j]) : 0;
          a[i*nv+j] = 0.;
        }
    }

#ifdef SPECTRAL_DEBUG
  printf ("G = \n");
  for (i = 0; i < nv; ++i)
    {
      for (j = 0; j < nv; ++j)
        printf (" %-4.5f", a[i*nv+j]);
      printf ("\n");
    }
#endif

  err = LAPACKE_dsyevd (LAPACK_ROW_MAJOR, 'V', 'U', nv, a, nv, d);
  if (err == 0)
    {
      int k = 0;
      while (k < nv && d[++k] < EPS)
        ;
#ifdef SPECTRAL_DEBUG  
      printf ("Eigenvector of the second smallest eigenvalue (%d: %.5f):\n",
              k, d[k]);
#endif
      
      for (i = 0; i < nv; ++i)
        {
          spectrum[i] = d[i];
          fiedler[i] = a[i*nv+k];
#ifdef SPECTRAL_DEBUG   
          printf ("% 3d: % 11.10f\n", i, fiedler[i]);
#endif
        }
    }

  free (d);
  free (a);

  return err;
}

#else

static int
graph_spectrum (float *spectrum, float *fiedler, const inchi_t *g)
{
  double **evec, **a, *d;
  int i, err = 0, nv = inchi_node_count (g);
  const int *G = inchi_matrix_A (g);
  size_t size = inchi_matrix_size (g);
  
  a = malloc (sizeof (double *)*nv);
  d = malloc (sizeof (double)*nv);
  evec = malloc (sizeof (double *)*nv);
  for (i = 0; i < nv; ++i)
    {
      a[i] = malloc (nv*sizeof (double));
      evec[i] = malloc (nv*sizeof (double));
    }

  spectral_normalized_graph (a, G, nv, size);

#ifdef SPECTRAL_DEBUG
  printf ("G = \n");
  for (i = 0; i < nv; ++i)
    {
      int j;
      for (j = 0; j < nv; ++j)
        printf ("%s%.1f", j==0?"":"&", a[i][j]);
      printf (" \\\\\n");
    }

  printf ("[");
  for (i = 0; i < nv; ++i)
    {
      int j;
      for (j = 0; j < nv; ++j)
        printf (" %.5f", a[i][j]);
      printf (";\n");
    }
  printf ("]\n");
#endif

  err = jacobi (a, nv, d, evec);

  { int k = 0;
    while (k < nv && d[++k] < EPS)
      ;
    
#ifdef SPECTRAL_DEBUG    
    printf ("Eigenvector of the second smallest eigenvalue (%d: %.5f):\n",
            k, d[k]);
#endif
    
    for (i = 0; i < nv; ++i)
      {
        spectrum[i] = d[i];
        fiedler[i] = evec[i][k];
#ifdef SPECTRAL_DEBUG
        printf ("% 3d: % 11.10f\n", i, fiedler[i]);
#endif
        free (a[i]);    
      }
    
    for (i = 0; i < nv; ++i)
      free (evec[i]);
  }
  
  free (d);
  free (a);
  free (evec);

  return err;
}
#endif /* !HAVE_GSL */


static int
spectral_inchi (spectral_t *sp, const char *inchi)
{
  int nv;

  nv = inchi_parse (sp->inchi, inchi);
  if (nv < 0)
    (void) strcpy (sp->errmsg, inchi_error (sp->inchi));
  else if (nv > SPECTRAL_MAXG)
    {
      sprintf (sp->errmsg, "Graph is too large (%d > %d) for eigensolver",
               nv, SPECTRAL_MAXG);
      nv = -1;
    }
  else
    {
      if (sp->bsize < nv)
        {
          sp->spectrum = realloc (sp->spectrum, nv*sizeof (float));
          sp->fiedler = realloc (sp->fiedler, nv*sizeof (float));
          sp->bsize = nv;
        }
#ifdef SPECTRAL_DEBUG
      printf ("## %d /c = %s\n", nv, inchi_layer_c (sp->inchi));
#endif
      
      if (graph_spectrum (sp->spectrum, sp->fiedler, sp->inchi) < 0)
        {
          sprintf (sp->errmsg, "Eigensolver didn't converge within "
                   "specified number of iterations");
          nv = -1;
        }
    }
  
  return nv;
}

#undef __set_edge
#undef __get_edge


spectral_t *
spectral_create ()
{
  spectral_t *sp = malloc (sizeof (struct __spectral_s));
  if (sp != 0)
    {
      sp->bsize = 0;
      sp->spectrum = 0;
      sp->fiedler = 0;
      sp->inchi = inchi_create ();
      (void) memset (sp->hashkey, 0, sizeof (sp->hashkey));
      (void) memset (sp->errmsg, 0, sizeof (sp->errmsg));
      sp->sha1 = sha1_create ();
    }
  return sp;
}

void
spectral_free (spectral_t *sp)
{
  if (sp != 0)
    {
      if (sp->spectrum != 0)
        free (sp->spectrum);
      if (sp->fiedler != 0)
        free (sp->fiedler);
      sha1_free (sp->sha1);
      inchi_free (sp->inchi);
      free (sp);
    }
}

const float *
spectral_spectrum (const spectral_t *sp)
{
  return sp->spectrum;
}

const float *
spectral_fiedler (const spectral_t *sp)
{
  return sp->fiedler;
}

size_t
spectral_size (const spectral_t *sp)
{
  return sp->bsize;
}

const char *
spectral_hashkey (const spectral_t *sp)
{
  return sp->hashkey;
}

const char *
spectral_error (const spectral_t *sp)
{
  return sp->errmsg;
}

const char *
spectral_version ()
{
  return SPECTRAL_VERSION;
}

const char *
spectral_digest (spectral_t *sp, const char *inchi)
{
  char *start;
  int size;

  size = spectral_inchi (sp, inchi);
  if (size < 0)
    return 0;

  sha1_reset (sp->sha1);
  
  /*
   * first block is topology
   */
  digest_spectrum (sp, size);
  sha1_digest (sp->sha1, sp->digest);
  start = sp->hashkey;
  b32_encode45 (&start, sp->digest, 20); /* 9 chars */
  
  /*
   * second block is connection
   */
  sha1_reset (sp->sha1);
  sha1_update (sp->sha1, sp->digest, 20); /* chaining */
  if (size > 0)
    {
      const char *inchi_c = inchi_layer_c (sp->inchi);
      sha1_update (sp->sha1, (const unsigned char *)inchi_c, strlen (inchi_c));
    }
  
  sha1_digest (sp->sha1, sp->digest);
  start = sp->hashkey + 9;
  b32_encode50 (&start, sp->digest, 20); /* 10 chars */
  
  /*
   * final block is the full inchi
   */
  sha1_reset (sp->sha1);
  sha1_update (sp->sha1, sp->digest, 20);
  sha1_update (sp->sha1, (const unsigned char *)inchi, strlen (inchi));
  sha1_digest (sp->sha1, sp->digest);
  start = sp->hashkey + 19;
  b32_encode55 (&start, sp->digest, 20); /* 11 chars */

  return sp->hashkey;
}

int
spectral_ratio (double *ratio, spectral_t *sp, const char *inchi)
{
  int i, size = spectral_inchi (sp, inchi);
  if (size < 0)
    return -1;

  i = 0;
  /* skip over all disconnected components */
  while (i < size && sp->spectrum[++i] < EPS)
    ;

#if 0
  { int j = i;
    fprintf (stderr, "pi: %.5f\n", sp->spectrum[6]/sp->spectrum[3]);
    for (; j < size; ++j)
      {
        fprintf (stderr, "%3d: %.5f %.5f\n",
                 j, sp->spectrum[j]/sp->spectrum[i], sp->spectrum[j]);
      }
  }
#endif
  
  *ratio = sp->spectrum[size-1]/(size*sp->spectrum[i]);
  return 0;
}
