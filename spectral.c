
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "b32.h"
#include "interval.h"
#include "spectral.h"

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
# define SPECTRAL_MAXG 300
#endif

/*
 * predefined precision; careful changing this value will affect
 * the hashkey!
 */
#define ROUND_OFF 0.0005

#define __set_edge(i,j) G[(i)*size+(j)] = G[(j)*size+(i)] = 1
#define __get_edge(i,j) G[(i)*size+(j)]

/**
 * internal state of spectral_t
 */
struct __spectral_s 
{
  size_t bsize; /* buffer size of spectrum */
  float *spectrum; /* spectrum buffer */
  char *inchi_c; /* connection layer */
  size_t isize; /* connection buffer size */
  sha1_t *sha1; /* sha1 hash */
  interval_t *itree; /* interval tree for quantization */
  unsigned char digest[20]; /* digest buffer */
  char hashkey[31]; /* 9(topology) + 10(connection) + 11(full) */
  char errmsg[BUFSIZ];
};


static void
digest_spectrum (spectral_t *sp, int size)
{
  unsigned int uv;
  unsigned char data[4];
  int i = 0;

#ifdef SPECTRAL_DEBUG
  printf ("## Eigen values:\n");
  for (; i < size; ++i)
      printf ("%3d: %.10f\n", i, sp->spectrum[i]);
  printf ("## Encoded values:\n");
  i = 0;
#endif

  /* skip over all disconnected components */
  while (i < size && sp->spectrum[++i] < EPS)
    ;

  for (; i < size; ++i)
    {
      uv = interval_encode32 (sp->itree, sp->spectrum[i], ROUND_OFF);
#ifdef SPECTRAL_DEBUG
      printf ("%3d: %.10f => %u\n", i, sp->spectrum[i], uv);
#endif
      data[0] = uv >> 24;
      data[1] = (uv & 0x00ffffff) >> 16;
      data[2] = (uv & 0x0000ffff) >> 8;
      data[3] = uv & 0xff;
      sha1_update (sp->sha1, data, sizeof (data));
    }
}

void
spectral_adjacency_graph (float **M, const int *G, int nv, size_t size)
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
spectral_laplacian_graph (float **M, const int *G, int nv, size_t size)
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
spectral_signless_graph (float **M, const int *G, int nv, size_t size)
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
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  int i, j, v, *d;
  double x;
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
  printf ("G = \n");
  for (i = 0; i < nv; ++i)
    {
      for (j = 0; j < nv; ++j)
        printf (" %-4.5f", gsl_matrix_get (A, i, j));
      printf ("\n");
    }
#endif
  
  gsl_eigen_symmv (A, L, V, ws);
  gsl_eigen_symmv_sort (L, V, GSL_EIGEN_SORT_VAL_ASC);
  
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
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  double *a, *d;
  int i, j, v, err = 0;

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
  for (i = 0; i < nv; ++i)
    spectrum[i] = d[i];

  free (d);
  free (a);

  return err;
}

#else

static int
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  double **evec, **a, *d;
  int i, err = 0;

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
    printf (" % 4d", i+1);
  printf ("\n");
  for (i = 0; i < nv; ++i)
    {
      int j;
      for (j = 0; j < nv; ++j)
        printf (" % 4.1f", a[i][j]);
      printf (" :%3d\n", i+1);
    }
#endif
  err = jacobi (a, nv, d, evec);

  for (i = 0; i < nv; ++i)
    {
      spectrum[i] = d[i];
      free (a[i]);
      free (evec[i]);
    }
  free (d);
  free (a);
  free (evec);

  return err;
}
#endif /* !HAVE_GSL */

static int
parse_inchi_graph (int **pG, size_t *psize, char *inchi, char errmsg[])
{
  char pc;
  int *G, *ppv, *pv, nv = 0, vv = 0;
  size_t size = strlen (inchi);
  char *ptr = inchi;
  char *end = inchi + size;

  if (size > *psize)
    {
      *pG = realloc (*pG, size*size*sizeof (int));
      *psize = size;
    }
  else
    size = *psize;

  G = *pG;
  { int i = 0, j;
    for (; i < size; ++i)
      for (j = 0; j < size; ++j)
        G[i*size+j] = 0;
  }

  ppv = pv = malloc (sizeof (int)*size);
  for (pc = 0; ptr < end; ++ptr)
    {
      int v = strtol (ptr, &ptr, 10);
      if (v > nv)
        nv = v;

      switch (pc) 
        {
        case '(': /* push */
          *pv++ = vv;
          /* fall through */

        case '-':
          __set_edge (vv, v);
          break;
          
        case ')': /* pop */
          if (pv > ppv)
            {
              --pv;
              __set_edge (*pv, v);
            }
          else
            sprintf (errmsg, "Mismatch ()'s in connection layer");
          break;

        case ',':
          if (pv > ppv)
            {
              __set_edge (pv[-1], v);
            }
          else
            sprintf (errmsg, "Character ',' not within () block");
          break;

        case '/':
        case '\0': /* end */
          break;

        case ';': /* component */
          
          break;

        case '*': /* multiplicity */
          nv = v; /* reset */
          vv = 0;
          break;

        default:
          sprintf (errmsg, "Unknown character '%c' in connection layer", pc);
          free (ppv);
          return -1;
        }
      pc = *ptr;
      vv = v;
    }
  free (ppv);

  return nv;
}

static int
spectral_inchi (spectral_t *sp, const char *inchi)
{
  char *ptr, *start = strstr (inchi, "/c"), *end;
  int nv = 0, *G = 0, *pG = 0;
  size_t size = 0, psize = 0;

  if (strncmp ("InChI=", inchi, 6) != 0)
    {
      sprintf (sp->errmsg, "Inchi string doesn't begins with InChI=");
      return -1;
    }
  else if (start == 0)
    {
      sprintf (sp->errmsg, "InChI string doesn't have connection layer");
      return 0;
    }

  ptr = start + 2; /* skip over /c */
  for (end = ptr; *end != '/' && !isspace (*end) && *end != '\0'; ++end)
    ;

  size = end - ptr;
  start = malloc (size + 1);
  (void) strncpy (start, ptr, size);
  start[size] = '\0';

  size = 0;
  ptr = start;
  for (ptr = strtok_r (ptr, ";", &end); 
       ptr != 0; ptr = strtok_r (0, ";", &end))
    {
      int n = parse_inchi_graph (&pG, &psize, ptr, sp->errmsg);
      if (n > nv)
        {
          if (G != 0)
            free (G);

          G = pG;
          size = psize;
          nv = n;

          { char *p = ptr;
            while (*p != '\0' && *p != '*')
              ++p;

            /* don't retain the multiplicity character from the /c layer */
            if (*p == '*')
              ptr = ++p;
          }
          psize = strlen (ptr);
          if (psize > sp->isize)
            {
              sp->inchi_c = realloc (sp->inchi_c, psize+1);
              sp->isize = psize;
            }
          (void) strcpy (sp->inchi_c, ptr);

          pG = 0;
          psize = 0;
        }
    }
  free (start);

  if (pG != 0)
    free (pG);

  if (nv > SPECTRAL_MAXG)
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
          sp->bsize = nv;
        }
#ifdef SPECTRAL_DEBUG
      printf ("## /c = %s\n", sp->inchi_c);
#endif
      if (graph_spectrum (sp->spectrum, G, nv, size) < 0)
        {
          sprintf (sp->errmsg, "Eigensolver didn't converge within "
                   "specified number of iterations");
          nv = -1;
        }
    }
  free (G);

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
      sp->inchi_c = 0;
      sp->isize = 0;
      (void) memset (sp->hashkey, 0, sizeof (sp->hashkey));
      (void) memset (sp->errmsg, 0, sizeof (sp->errmsg));
      sp->sha1 = sha1_create ();
      /* eigenvalues of a normalized laplacian is bounded by [0,2] */
      sp->itree = interval_create (0., 2., 32);
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
      if (sp->inchi_c != 0)
        free (sp->inchi_c);
      sha1_free (sp->sha1);
      interval_free (sp->itree);
      free (sp);
    }
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
    sha1_update (sp->sha1, (unsigned char *)sp->inchi_c, 
                 strlen (sp->inchi_c));
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
