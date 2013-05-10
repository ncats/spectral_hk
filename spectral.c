
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "b32.h"
#include "interval.h"
#include "spectral.h"

/**
 * unfortunately due to numerical round offs, we can't
 * have hashkeys that are independent of eigensolvers!
 */
#ifdef HAVE_GSL
# include <gsl/gsl_eigen.h>
#else
# include "jacobi.h"
#endif

#ifndef EPS
# define EPS 1e-20
#endif

/*
 * maximum graph size; for large graphs a more specialized eigensolver
 * (e.g., iterative ones) are more appropriate.
 */
#ifndef SPECTRAL_MAXG
# define SPECTRAL_MAXG 200
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
  size_t bufsiz; /* buffer size of spectrum */
  float *spectrum; /* spectrum buffer */
  sha1_t *sha1; /* sha1 hash */
  interval_t *itree; /* interval tree for quantization */
  unsigned char digest[20]; /* digest buffer */
  char hashkey[25]; /* 7(topology) + 8(connection) + 9(full) */
  char errmsg[BUFSIZ];
};


static void
digest_spectrum (spectral_t *sp, int size)
{
  unsigned int uv;
  unsigned char data[4];
  int i = 0;

  /* skip over all disconnected components */
  while (sp->spectrum[++i] < EPS && i < size)
    ;
  for (; i < size; ++i)
    {
      uv = interval_encode32 (sp->itree, sp->spectrum[i], ROUND_OFF);
#ifdef SPECTRAL_DEBUG
      printf ("%3d: %.10f %u\n", i, sp->spectrum[i], uv);
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
spectral_normalized_graph (float **M, const int *G, int nv, size_t size)
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
        M[i][j] = M[j][i] = -1./sqrt (d[i]*d[j]);
    }
  free (d);
}

#ifdef HAVE_GSL
static int
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  int i, j;
  float **a;
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc (nv);
  gsl_matrix *A = gsl_matrix_alloc (nv, nv);
  gsl_matrix *V = gsl_matrix_alloc (nv, nv);
  gsl_vector *L = gsl_vector_alloc (nv);

  a = malloc (sizeof (float *)*nv);
  for (i = 0; i < nv; ++i)
    a[i] = malloc (nv*sizeof (float));
  
  spectral_normalized_graph (a, G, nv, size);

  /* copy over the matrix */
  for (i = 0; i < nv; ++i)
    {
      for (j = i+1; j < nv; ++j)
        {
          gsl_matrix_set (A, i, j, a[i][j]);
          gsl_matrix_set (A, j, i, a[j][i]);
        }
      gsl_matrix_set (A, i, i, a[i][i]);
    }
  
  gsl_eigen_symmv (A, L, V, ws);
  gsl_eigen_symmv_sort (L, V, GSL_EIGEN_SORT_VAL_ASC);
  
  for (i = 0; i < nv; ++i)
    {
      spectrum[i] = gsl_vector_get (L, i);
      free (a[i]);
    }
  free (a);
    
  gsl_vector_free (L);
  gsl_matrix_free (V);
  gsl_matrix_free (A);
  gsl_eigen_symmv_free (ws);

  return 0;
}

#else

#warning **** Please consider using the GSL eigensolver. It's an order \
of magnitude faster! The bundled implementation is only for completeness \
sake. ****

static int
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  float **evec, **a;
  int i, err = 0;

  a = malloc (sizeof (float *)*nv);
  evec = malloc (sizeof (float *)*nv);
  for (i = 0; i < nv; ++i)
    {
      a[i] = malloc (nv*sizeof (float));
      evec[i] = malloc (nv*sizeof (float));
    }

  spectral_normalized_graph (a, G, nv, size);
  err = jacobi (a, nv, spectrum, evec);

  for (i = 0; i < nv; ++i)
    {
      free (a[i]);
      free (evec[i]);
    }
  free (a);
  free (evec);

  return err;
}
#endif /* !HAVE_GSL */

static int
spectral_inchi (spectral_t *sp, const char *inchi)
{
  char *ptr, *start = strstr (inchi, "/c"), *end;
  int nv = 0, vv = 0, *G, *pv, *ppv;
  char pc;
  size_t size;

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

  size = end - start;   /* this is an over estimate of the graph size */
  G = malloc (size * size * sizeof (int));
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
            sprintf (sp->errmsg, "Mismatch ()'s in connection layer");
          break;

        case ',':
          if (pv > ppv)
            {
              __set_edge (pv[-1], v);
            }
          else
            sprintf (sp->errmsg, "Character ',' not within () block");
          break;

        case '/':
        case '\0':
        case ';': /* component */
          /* end */
          break;

        case '*': /* multiplicity */
          nv = v; /* reset */
          vv = 0;
          break;

        default:
          sprintf (sp->errmsg, 
                   "Unknown character '%c' in connection layer", pc);
          return -1;
        }
      pc = *ptr;
      vv = v;
    }

  if (nv > SPECTRAL_MAXG)
    {
      sprintf (sp->errmsg, "Graph is too large (%d > %d) for eigensolver",
               nv, SPECTRAL_MAXG);
      nv = -1;
    }
  else
    {
      if (sp->bufsiz < nv)
        {
          sp->spectrum = realloc (sp->spectrum, nv*sizeof (float));
          sp->bufsiz = nv;
        }

      if (graph_spectrum (sp->spectrum, G, nv, size) < 0)
        {
          sprintf (sp->errmsg, "Eigensolver didn't converge within "
                   "specified number of iterations");
          nv = -1;
        }
    }

  free (ppv);
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
      sp->bufsiz = 0;
      sp->spectrum = 0;
      memset (sp->hashkey, 0, sizeof (sp->hashkey));
      memset (sp->errmsg, 0, sizeof (sp->errmsg));
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
spectral_digest (spectral_t *sp, const char *inchi)
{
  char *start, *end;
  int size;

  size = spectral_inchi (sp, inchi);
  if (size < 0)
    ;
  else
    {
      sha1_reset (sp->sha1);

      /*
       * first block is topology
       */
      digest_spectrum (sp, size);
      sha1_digest (sp->sha1, sp->digest);
      start = sp->hashkey;
      b32_encode35 (&start, sp->digest, 5); /* 7 char */

      /*
       * second block is connection
       */
      sha1_reset (sp->sha1);
      sha1_update (sp->sha1, sp->digest, 20); /* chaining */
      start = strstr (inchi, "/c");
      if (start != 0)
        {
          for (end = start+2; *end != '/' && *end != '\0'; ++end)
            ;
          sha1_update (sp->sha1, (unsigned char *)start, end - start);
        }
      sha1_digest (sp->sha1, sp->digest);
      start = sp->hashkey + 7;
      b32_encode40 (&start, sp->digest, 5); /* 8 char's */

      /*
       * final block is the full inchi
       */
      sha1_reset (sp->sha1);
      sha1_update (sp->sha1, sp->digest, 20);
      sha1_update (sp->sha1, (const unsigned char *)inchi, strlen (inchi));
      sha1_digest (sp->sha1, sp->digest);
      start = sp->hashkey + 15;
      b32_encode45 (&start, sp->digest, 6);
    }

  return size < 0 ? 0 : sp->hashkey;
}
