
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "spectral.h"
#include "b32.h"

/**
 * unfortunately due to numerical round offs, we can't
 * have hashkeys that are independent of eigensolvers!
 */
#ifdef HAVE_GSL
# include <gsl/gsl_eigen.h>
#else
# include "jacobi.h"
#endif

#define MATRIX_LAPLACIAN

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
#define PRECISION 0.0001


#define __set_edge(i,j) G[(i)*size+(j)] = G[(j)*size+(i)] = 1
#define __get_edge(i,j) G[(i)*size+(j)]

/**
 * internal state of spectral_t
 */
struct __spectral_t 
{
  size_t bufsiz; /* buffer size of spectrum */
  float *spectrum; /* spectrum buffer */
  sha1_state *sha1;
  unsigned char digest[20]; /* sha1 digest buffer */
  char hashkey[25]; /* 7(topology) + 8(connection) + 9(full) */
  char errmsg[BUFSIZ];
};

static void
encode_rational32 (unsigned char data[4], float x)
{
  union {
    unsigned int ival;
    float fval;
  } u;

  u.fval = floorf ((x + PRECISION)/PRECISION) * PRECISION; /* rounding */
  u.ival&= 0xffff0000;

#ifdef SPECTRAL_DEBUG
  { int e, m;
    //u.fval = x;
    e = (u.ival >> 23) & 0xff;
    m = (e == 0) ? (u.ival & 0x7fffff)<<1 : (u.ival&0x7fffff) | 0x800000;
    printf ("%.10f: %d %d %.10f %f e=%d m=%d %.10f\n", 
            x, u.ival, u.ival & 0xffff0000, u.fval, rintf (x), e, m, m*powf(2., e-150));
  }
#endif

  data[0] = u.ival >> 24;
  data[1] = (u.ival & 0x00ffffff) >> 16;
  data[2] = (u.ival & 0x0000ffff) >> 8;
  data[3] = u.ival & 0xff;
}

#ifdef HAVE_GSL
static int
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  int i, j, v;
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc (nv);
  gsl_matrix *A = gsl_matrix_alloc (nv, nv);
  gsl_matrix *V = gsl_matrix_alloc (nv, nv);
  gsl_vector *L = gsl_vector_alloc (nv);

#if defined(MATRIX_LAPLACIAN) || defined(MATRIX_SIGNLESS)
  /* 
   * laplacian matrix
   */
  for (i = 0; i < nv; ++i)
    {
      int d = 0;
      for (j = 0; j < nv; ++j)
        {
          if (i != j)
            {
              v = __get_edge (i+1, j+1);
              d += v;
              gsl_matrix_set (A, i, j, v ?
#if defined(MATRIX_SIGNLESS)
                              1.
#else
                              -1.
#endif
                              :0.);
            }
        }
      gsl_matrix_set (A, i, i, d);
    }
#else
  /*
   * adjajency 
   */
  for (i = 0; i < nv; ++i)
    {
      for (j = i+1; j < nv; ++j)
        {
          v = __get_edge (i+1, j+1);
          gsl_matrix_set (A, i, j, v);
          gsl_matrix_set (A, j, i, v);
        }
      gsl_matrix_set (A, i, i, 0.);
    }
#endif
  
  v = gsl_eigen_symmv (A, L, V, ws);
  gsl_eigen_symmv_sort (L, V, GSL_EIGEN_SORT_VAL_ASC);
  
  for (i = 0; i < nv; ++i)
    spectrum[i] = gsl_vector_get (L, i);
    
  gsl_vector_free (L);
  gsl_matrix_free (V);
  gsl_matrix_free (A);
  gsl_eigen_symmv_free (ws);

  return 0;
}

#else

static int
graph_spectrum (float *spectrum, const int *G, int nv, size_t size)
{
  float **evec, **a;
  int v, i, j, err = 0;

  a = malloc (sizeof (float *)*nv);
  evec = malloc (sizeof (float *)*nv);

#if defined(MATRIX_LAPLACIAN) || defined(MATRIX_SIGNLESS)
  /* 
   * laplacian matrix
   */
  for (i = 0; i < nv; ++i)
    {
      int d = 0;
      a[i] = malloc (nv*sizeof (float));
      for (j = 0; j < nv; ++j)
        if (i != j)
          {
            v = __get_edge (i+1, j+1);
            d += v;
#if defined(MATRIX_SIGNLESS)
            a[i][j] = v ? 1 : 0.;
#else
            a[i][j] = v ? -1 : 0.;
#endif
          }
      a[i][i] = d;
      evec[i] = malloc (nv*sizeof (float));
    }

#else 
  /*
   * adjajency 
   */
  for (i = 0; i < nv; ++i)
    {
      a[i] = malloc (nv*sizeof (float));
      for (j = i+1; j < nv; ++j)
        {
          v = __get_edge (i+1, j+1);
          a[i][j] = v;
        }

      for (j = 0; j < i; ++j)
        a[i][j] = a[j][i];

      a[i][i] = 0.;
      evec[i] = malloc (nv*sizeof (float));
    }
#endif

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
#endif

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

  size = end - start;
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
          /* end */
          break;

        case ';': /* component */
        case '*': /* multiplicity */
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
  spectral_t *sp = malloc (sizeof (struct __spectral_t));
  if (sp != 0)
    {
      sp->bufsiz = 0;
      sp->spectrum = 0;
      memset (sp->hashkey, 0, sizeof (sp->hashkey));
      memset (sp->errmsg, 0, sizeof (sp->errmsg));
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
      sha1_free (sp->sha1);
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
  unsigned char data[4];
  int i, size;

  size = spectral_inchi (sp, inchi);
  if (size < 0)
    ;
  else
    {
      sha1_reset (sp->sha1);
      /*
       * first block is topology
       */
      i = 0;
#if defined(MATRIX_LAPLACIAN)
      /* skip over all disconnected components */
      while (sp->spectrum[++i] < EPS)
        ;
#endif
      for (; i < size; ++i)
        {
          encode_rational32 (data, sp->spectrum[i]);
          sha1_update (sp->sha1, data, sizeof (data));
        }
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
