
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "b32.h"
#include "interval.h"
#include "spectral.h"
#include "periodic.h"

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

typedef struct __formula_s {
  int index; /* this is not the same the component index */
  int count;
  const element_t *element;
  struct __formula_s *next;
} formula_t;

typedef struct __hlayer_s {
  int index; /* component index */
  int atom; /* atom index */
  int count; /* number of Hs */
  int group; /* shared Hs if group > 1 */
  const element_t *element;
  struct __hlayer_s *next;
} hlayer_t;

typedef struct __vertex_s {
  int index;
  int degree;
  int hcount;
  int charge;
  const element_t *atom;
  struct __vertex_s **neighbors; /* neighbors[0..degree-1] */
} vertex_t;

typedef struct __graph_s {
  int index; /* component index */
  int multiplier; /* multiplicity */
  int nv; /* number of vertices */
  int *A; /* adjacency matrix */
  vertex_t **G; /* adjacency list G[0..nv-1] */
  size_t size; /* stride size */
  formula_t *formula;
  hlayer_t *hlayer;
} graph_t;

/*
 * maximum graph size; for large graphs a more specialized eigensolver
 * (e.g., iterative ones) are more appropriate.
 */
#ifndef SPECTRAL_MAXG
# define SPECTRAL_MAXG 5000
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
  float *fiedler; /* fiedler vector */
  char *inchi_c; /* connection layer */
  size_t isize; /* connection buffer size */
  sha1_t *sha1; /* sha1 hash */
  interval_t *itree; /* interval tree for quantization */
  unsigned char digest[20]; /* digest buffer */
  char hashkey[31]; /* 9(topology) + 10(connection) + 11(full) */
  char errmsg[BUFSIZ];
};


static void
__digest_spectrum (spectral_t *sp, int size)
{
  unsigned int uv;
  unsigned char data[4];
  int i = 0;

#ifdef SPECTRAL_DEBUG
  printf ("## Eigen values:\n");
  { double p = 0.;
    for (; i < size; ++i)
      {
        printf ("%3d: %.10f %.10f %.10f\n", 
                i, sp->spectrum[i], sp->spectrum[i] - p, 
                sp->spectrum[i]/sp->spectrum[1]);
        p = sp->spectrum[i];
      }
  }
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
      printf ("%3d: %.10f => %3u %10u\n", j, sp->spectrum[j], uv, 
              interval_encode32 (sp->itree, x, ROUND_OFF));
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
graph_spectrum (float *spectrum, float *fiedler, const graph_t *g)
{
  int i, j, v, *d, nv = g->nv;
  double x;
  const int *G = g->A;
  size_t size = g->size;
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
graph_spectrum (float *spectrum, float *fiedler, const graph_t *g)
{
  double *a, *d;
  const int *G = g->A;
  size_t size = g->size;
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
graph_spectrum (float *spectrum, float *fiedler, const graph_t *g)
{
  double **evec, **a, *d;
  int i, err = 0, nv = g->nv;
  const int *G = g->A;
  size_t size = g->size;
  
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

static formula_t *
create_formula (int index, int count, const element_t *el)
{
  formula_t *f = malloc (sizeof (formula_t));
  f->index = index;
  f->count = count;
  f->element = el;
  f->next = 0;
  return f;
}

static void
destroy_formula (formula_t *f)
{
  while (f != 0)
    {
      formula_t *next = f->next;
      free (f);
      f = next;
    }
}

static int
update_formula (formula_t **head, formula_t **current,
                int *index, int multi, int count, const element_t *el)
{
  formula_t *f = create_formula (*index, count, el);
  if (*current != 0)
    (*current)->next = f;
  else
    *head = f;
  *current = f;
  
  /* expand */
  if (multi > 1)
    {
      int k;
      formula_t *h = *head;
      while (h->index != *index)
        h = h->next;
      
      for (k = 1; k < multi; ++k)
        {
          formula_t *p = h;
          ++*index;
          while (p->index == h->index)
            {
              f = create_formula (*index, p->count, p->element);
              (*current)->next = f;
              *current = f;
              p = p->next;
            }
        }
    }
  
  return multi*count;
}

static int
parse_formula (formula_t **formula, char *err, const char *inchi)
{
  char *p = strchr (inchi, '/');
  const element_t *el = 0;
  int count = 0, multi = 1, index = 0, total = 0;
  formula_t *current = 0, *head = 0;

  if (p == 0)
    {
      *formula = 0;
      return 0;
    }
  
  ++p; /* skip over / */
  while (*p != '/')
    {
      if (isalpha (*p))
        {
          if (el != 0)
            {
              formula_t *f = create_formula (index, count, el);
              if (current != 0)
                current->next = f;
              else
                head = f;
              current = f;
              total += multi*count;
            }
          
          el = element_lookup_symbol (p);
          if (el != 0)
            {
              count = 1;
              p += strlen (el->symbol);
            }
          else
            {
              sprintf (err,
                       "** Unknown atom in formula at position %ld: %s **\n",
                       p - inchi, p);
              el = 0;
              count = 0;
              multi = 1;
              ++p;
            }
        }
      else if (isdigit (*p))
        {
          /* number */
          count = (int) strtol (p, &p, 10);
          if (el == 0)
            {
              /* multiplicity.. */
              multi = count;
            }
        }
      else if (*p == '.')
        {
          if (el != 0)
            {
              total += update_formula (&head, &current,
                                       &index, multi, count, el);
            }
          el = 0;
          count = 0;
          multi = 1;
          ++index;
          ++p;
        }
      else
        sprintf (err, "** Uknown character in formula: %c **\n", *p);
    }
  
  if (el != 0)
    {
      total += update_formula (&head, &current, &index, multi, count, el);
    }
  *formula = head;
  
  return total;
}

static hlayer_t *
create_hlayer (int index, int atom)
{
  hlayer_t *h = malloc (sizeof (hlayer_t));
  h->index = index;
  h->atom = atom;
  h->count = 0;
  h->group = 0;
  h->element = 0;
  h->next = 0;
  return h;
}

static void
destroy_hlayer (hlayer_t *h)
{
  while (h != 0)
    {
      hlayer_t *next = h->next;
      free (h);
      h = next;
    }
}

static int
parse_layer_h (hlayer_t **hlayer, char *err, const char *inchi)
{
  char *p = strstr (inchi, "/h");
  int pn = 0, n, count, group = 0, shared = 0, index = 0, total = 0;
  hlayer_t *head = 0, *current = 0;
  
  if (p == 0)
    {
      *hlayer = 0;
      return 0;
    }

  p += 2; /* skip over /h */
  while (*p != '/' && *p != '\0')
    {
      switch (*p)
        {
        case '1': case '2': case '3':
        case '4': case '5': case '6':
        case '7': case '8': case '9':
          n = strtol (p, &p, 10);
          if (n > 0)
            {
              if (*p != '*')
                {
                  hlayer_t *h = create_hlayer (index, n);
                  if (shared)
                    {
                      h->count = count;
                      h->group = group;
                    }
                  
                  if (head == 0)
                    head = h;
                  else
                    current->next = h;
                  current = h;
                  ++total;
                }
              else
                sprintf (err, "** Warning: ignore multiplicity %d "
                         "in /h layer **\n", n);
            }
          break;

        case '-':
          pn = n;
          n = strtol (p+1, &p, 10);
          if (n > 0 && pn > 0)
            { int p = pn+1;
              for (; p <= n; ++p)
                {
                  hlayer_t *h = create_hlayer (index, p);
                  if (head == 0)
                    head = h;
                  else
                    current->next = h;
                  current = h;
                  ++total;
                }
              pn = 0;
            }
          else
            {
              sprintf (err, "** Invalid range in /h layer: %d-%d **\n",
                       pn, n);
            }
          break;

        case ',':
          ++p;
          break;

        case '(':
          shared ^= 1; /* toggel parity */
          ++group;
          ++p;
          break;
          
        case ')':
          shared ^= 1;
          ++p;
          break;

        case 'H':
          count = 1;
          if (isdigit (p[1]))
            count = strtol (p+1, &p, 10);
          else
            ++p;

          if (!shared)
            {
              hlayer_t *h = head;
              while (h != 0)
                {
                  if (h->count == 0)
                    h->count = count;
                  h = h->next;
                }
            }
          break;

        case '*': /* multiplicity (ignore for now) */
          ++p;
          break;
          
        case ';': /* component */
          ++index;
          ++p;
          break;

        default:
          sprintf (err, "** Warning: unknown character in /h layer: '%c'\n",
                   *p);
          ++p;
        }
    } /* while() */

  *hlayer = head;
  return total;
}

static vertex_t *
create_vertex (int index, int degree)
{
  vertex_t *v = malloc (sizeof (vertex_t) + sizeof (vertex_t*)*degree);
  v->index = index;
  v->degree = degree;
  v->hcount = 0;
  v->charge = 0;
  v->atom = 0;
  v->neighbors = (vertex_t **)((char *)v + sizeof (vertex_t));
  (void) memset (v->neighbors, 0, degree*sizeof (vertex_t *));
  return v;
}

static void
destroy_vertex (vertex_t *v)
{
  if (v != 0)
    free (v);
}

static void
destroy_graph (graph_t *g)
{
  if (g->A != 0)
    free (g->A);
  if (g->G != 0)
    {
      int i;
      for (i = 0; i < g->nv; ++i)
        destroy_vertex (g->G[i]);
      free (g->G);
    }
  destroy_formula (g->formula);
  destroy_hlayer (g->hlayer);
}

static void
process_graph (graph_t *g)
{
  int i, j, count, total, *G = g->A;
  size_t size = g->size;
  formula_t *f = g->formula;
  hlayer_t *h;
  
  g->G = malloc (sizeof (vertex_t *) * g->nv);
  /* first pass to allocate the vertices */
  for (i = 0; i < g->nv; ++i)
    {
      int d = 0;
      for (j = 0; j < g->nv; ++j)
        if (i != j && __get_edge (i+1, j+1))
          ++d;
      g->G[i] = create_vertex (i+1, d);
    }
  
  /* align the formula with the component; this doesn't seem to be a clean
   * or right way to do this??  */
  while (f->index < g->index)
    f = f->next;

  count = 0;
  {
    formula_t *current = f;
    int index = f->index;
    while (f != 0 && index == f->index)
      {
        if (f->element->atno != 1)
          count += f->count;
        f = f->next;
      }
    f = current;
  }
  
  if (count != g->nv)
    printf ("** formula misaligned with component: expecting %d but got %d!\n",
            g->nv, count);
  
  printf ("graph G for component %d...%d x %d => formula %d\n",
          g->index, g->multiplier, g->nv, f->index);
  
  /* now do the linking */
  if (f->element->atno == 1)
    f = f->next;
  
  count = f->count;
  for (i = 0; i < g->nv; ++i)
    {
      int k = 0, valence;
      vertex_t *u = g->G[i];
      for (j = 0; j < g->nv; ++j)
        if (i != j && __get_edge (i+1, j+1))
          u->neighbors[k++] = g->G[j];

      u->atom = f->element;
      for (h = g->hlayer; h != 0; h = h->next)
        if (h->atom == u->index)
          {
            if (h->group > 0)
              {
                /* shared; ensure that only the smallest index get the H */
                hlayer_t *hl = g->hlayer;
                while (hl->group != h->group)
                  hl = hl->next;
                
                if (hl == h)
                  u->hcount = h->count;
              }
            else
              u->hcount = h->count;
          }

      printf ("%d %d %s:", u->index, u->hcount, u->atom->symbol);
      for (k = 0; k < u->degree; ++k)
        printf (" %d", u->neighbors[k]->index);
      printf ("\n");

      if (--count == 0)
        {
          f = f->next;
          if (f != 0)
            {
              /* don't include Hs */
              if (f->element->atno == 1)
                f = f->next;
              count = f->count;
            }
        }
    }
}

static int
spectral_inchi (spectral_t *sp, const char *inchi)
{
  char *ptr, *start, *end;
  int *pG = 0, n;
  size_t size = 0, psize = 0;
  graph_t g = {0};

  if (strncmp ("InChI=", inchi, 6) != 0)
    {
      sprintf (sp->errmsg, "Inchi string doesn't begin with InChI=");
      return -1;
    }

  start = strstr (inchi, "/c");
  if (start == 0)
    {
      sprintf (sp->errmsg, "InChI string doesn't have connection layer");
      return 0;
    }

  /* parse formula */
  n = parse_formula (&g.formula, sp->errmsg, inchi);
  { formula_t *formula = g.formula;
    printf ("formula: %d\n", n);
    while (formula != 0)
      {
        printf ("%d: %d %s\n", formula->index,
                formula->count, formula->element->symbol);
        formula = formula->next;
      }
  }
  
  /* parse h layer */
  n = parse_layer_h (&g.hlayer, sp->errmsg, inchi);
  { hlayer_t *h = g.hlayer;
    printf ("/h layer...%d\n", n);
    while (h != 0)
      {
        printf ("%d %d %d %d\n", h->index, h->atom, h->count, h->group);
        h = h->next;
      }
  }

  ptr = start + 2; /* skip over /c */
  for (end = ptr; *end != '/' && !isspace (*end) && *end != '\0'; ++end)
    ;

  size = end - ptr;
  start = malloc (size + 1);
  (void) strncpy (start, ptr, size);
  start[size] = '\0';

  ptr = start;
  for (ptr = strtok_r (ptr, ";", &end); 
       ptr != 0; ptr = strtok_r (0, ";", &end))
    {
      n = parse_inchi_graph (&pG, &psize, ptr, sp->errmsg);
      if (n > g.nv)
        {
          /* keep only the largest component */
          if (g.A != 0)
            free (g.A);

          g.A = pG;
          g.nv = n;
          g.size = psize;
          g.index = 0;

          { char *p = ptr;
            for (; p >= start; --p)
              if (*p == ';' || *p == '\0')
                ++g.index;

            p = ptr;
            g.multiplier = strtol (p, &p, 10);
            if (*p != '*')
              g.multiplier = 1;
            else
              /* don't retain the multiplicity character from the /c layer */
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

  if (g.nv > SPECTRAL_MAXG)
    {
      sprintf (sp->errmsg, "Graph is too large (%d > %d) for eigensolver",
               g.nv, SPECTRAL_MAXG);
      g.nv = -1;
    }
  else
    {
      if (sp->bsize < g.nv)
        {
          sp->spectrum = realloc (sp->spectrum, g.nv*sizeof (float));
          sp->fiedler = realloc (sp->fiedler, g.nv*sizeof (float));
          sp->bsize = g.nv;
        }
#ifdef SPECTRAL_DEBUG
      printf ("## /c = %s\n", sp->inchi_c);
#endif
      process_graph (&g);
      if (graph_spectrum (sp->spectrum, sp->fiedler, &g) < 0)
        {
          sprintf (sp->errmsg, "Eigensolver didn't converge within "
                   "specified number of iterations");
          g.nv = -1;
        }
    }
  destroy_graph (&g);

  return g.nv;
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
      sp->inchi_c = 0;
      sp->isize = 0;
      (void) memset (sp->hashkey, 0, sizeof (sp->hashkey));
      (void) memset (sp->errmsg, 0, sizeof (sp->errmsg));
      sp->sha1 = sha1_create ();
      /* eigenvalues of a normalized laplacian is bounded by [0,2] */
      sp->itree = interval_create (0., 1., 32);
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
      if (sp->inchi_c != 0)
        free (sp->inchi_c);
      sha1_free (sp->sha1);
      interval_free (sp->itree);
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

int
spectral_ratio (double *ratio, spectral_t *sp, const char *inchi)
{
  int i, j, size = spectral_inchi (sp, inchi);
  if (size < 0)
    return -1;

  i = 0;
  /* skip over all disconnected components */
  while (i < size && sp->spectrum[++i] < EPS)
    ;

#if 0
  fprintf (stderr, "pi: %.5f\n", sp->spectrum[6]/sp->spectrum[3]);
  for (j = i; j < size; ++j)
    {
      fprintf (stderr, "%3d: %.5f %.5f\n",
               j, sp->spectrum[j]/sp->spectrum[i], sp->spectrum[j]);
    }
#endif
  
  *ratio = sp->spectrum[size-1]/(size*sp->spectrum[i]);
  return 0;
}
