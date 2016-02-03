
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sha1.h"
#include "b32.h"
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
  struct __graph_s *graph;
  struct __edge_s **edges; /* edges[0..degree-1] */
} vertex_t;

typedef struct __edge_s {
  int index;
  int order;
  vertex_t *u;
  vertex_t *v;
  struct __edge_s *next;
} edge_t;

typedef struct __graph_s {
  int index; /* component index */
  int multiplier; /* multiplicity */
  
#ifdef __A
# undef __A
#endif
#define __A(i,j) (g->A[(i)*g->size+(j)])
  int *A; /* adjacency matrix */

#ifdef __L
# undef __L
#endif
#define __L(i,j) (g->L[(i)*g->size+(j)])
  double *L; /* normalized Laplacian: use macro __L for access */

#ifdef __W
# undef __W
#endif
#define __W(i,j) (g->W[(i)*g->size+(j)])
  double *W; /* weighted normalized Laplacian */
  
  int nv; /* number of vertices */
  vertex_t **V; /* adjacency list V[0..nv-1] */

  int ne; /* number of edges */
  edge_t *E; /* linked list of edges */
  
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
create_vertex (graph_t *g, int index, int degree)
{
  vertex_t *v = malloc (sizeof (vertex_t) + sizeof (edge_t*)*degree);
  v->index = index;
  v->degree = degree;
  v->hcount = 0;
  v->charge = 0;
  v->atom = 0;
  v->graph = g;  
  v->edges = (edge_t **)((char *)v + sizeof (vertex_t));
  (void) memset (v->edges, 0, degree*sizeof (edge_t *));
  return v;
}

static void
destroy_vertex (vertex_t *v)
{
  if (v != 0)
    free (v);
}

static void
destroy_edge (edge_t *e)
{
  edge_t *next;
  while (e != 0)
    {
      next = e->next;
      free (e);
      e = next;
    }
}

static void
destroy_graph (graph_t *g)
{
  static graph_t def = {0};
  
  if (g->A != 0)
    free (g->A);
  if (g->V != 0)
    {
      int i;
      for (i = 0; i < g->nv; ++i)
        destroy_vertex (g->V[i]);
      free (g->V);
    }
  
  if (g->L != 0)
    free (g->L);

  if (g->W != 0)
    free (g->W);

  destroy_edge (g->E);
  destroy_formula (g->formula);
  destroy_hlayer (g->hlayer);

  *g = def;
}

static edge_t *
create_edge (int index, vertex_t *u, vertex_t *v)
{
  edge_t *e = malloc (sizeof (edge_t));
  e->index = index;
  e->order = 1;
  e->u = u;
  e->v = v;
  e->next = 0;
  return e;
}

static vertex_t *
edge_other (const edge_t *e, const vertex_t *v)
{
  return e->u == v ? e->v : e->u;
}

static int
implicit_hcount (const vertex_t *u)
{
  int v = 0, k;
  for (k = 0; k < u->degree; ++k)
    v += u->edges[k]->order;
  
  switch (u->atom->atno)
    {
    case 5: v = 3 - v; break;
    case 6: v = 4 - v; break;
    case 7: case 15:
      if (v <= 3) v = 3 - v;
      else v = 5 - v;
      break;
    case 8: v = 2 - v; break;
    case 9: v = 1 - v; break;
    case 14: v = 4 - v; break;
    case 16:
      if (v <= 2) v = 2 - v;
      /*else if (v <= 4) v = 4 - v;*/
      else v = 6 - v;
      break;
    case 17: v = 1 - v; break;
    case 35: case 53: v = 1 - v; break;
    default: v = -1000;
    }
  
  return v;
}

/*
 * these are from http://cccbdb.nist.gov/
 */
static void
bond_energy_length (const edge_t *e,
                    double *d /* energy */, double *r /* length */)
{
  const vertex_t *u = e->u->atom->atno < e->v->atom->atno
    ? e->u : e->v;
  const vertex_t *v = edge_other (e, u);

  /*
   * d is kJ/mol
   * r is pm
   */
  if (d != 0) *d = 0;
  if (r != 0) *r = 1000;

  switch (u->atom->atno)
    {
    case 1:
      switch (v->atom->atno)
        {
        case 1: /* H-H */
          if (d != 0) *d = 432;
          if (r != 0) *r = 74;
          break;

        case 5: /* H-B */
          if (d != 0) *d = 389;
          if (r != 0) *r = 119;
          break;

        case 6: /* H-C */
          if (d != 0) *d = 411;
          if (r != 0) *r = 109;
          break;

        case 7: /* H-N */
          if (d != 0) *d = 386;
          if (r != 0) *r = 101;
          break;

        case 8: /* H-O */
          if (d != 0) *d = 459;
          if (r != 0) *r = 96;
          break;

        case 9: /* H-F */
          if (d != 0) *d = 565;
          if (r != 0) *r = 92;
          break;

        case 14: /* H-Si */
          if (d != 0) *d = 318;
          if (r != 0) *r = 148;
          break;

        case 15: /* H-P */
          if (d != 0) *d = 322;
          if (r != 0) *r = 144;
          break;

        case 16: /* H-S */
          if (d != 0) *d = 363;
          if (r != 0) *r = 134;
          break;

        case 17: /* H-Cl */
          if (d != 0) *d = 428;
          if (r != 0) *r = 127;
          break;

        case 35: /* H-Br */
          if (d != 0) *d = 362;
          if (r != 0) *r = 141;
          break;
          
        case 53: /* H-I */
          if (d != 0) *d = 295;
          if (r != 0) *r = 161;
        }
      break;

    case 5:
      switch (v->atom->atno)
        {
        case 5: /* B-B */
          if (d != 0) *d = 293;
          if (r != 0) *r = 170.2;
          break;

        case 6: /* B-C */
          if (d != 0) *d = 536;
          if (r != 0) *r = 149.1;
          break;

        case 9:
          if (d != 0) *d = 613;
          if (r != 0) *r = 130.7;
          break;

        case 17:
          if (d != 0) *d = 456;
          if (r != 0) *r = 175;
          break;

        case 35:
          if (d != 0) *d = 377;
          if (r != 0) *r = 188.8;
          break;

        default:
          fprintf (stderr, "** Error: No energy/length entries for %s-%s! **\n",
                   u->atom->symbol, v->atom->symbol);
        }
      break;

    case 6:
      switch (v->atom->atno)
        {
        case 6:
          if (e->order == 1)
            {
              if (d != 0) *d = 346;
              if (r != 0) *r = 154;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 602;
              if (r != 0) *r = 134;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 835;
              if (r != 0) *r = 120;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and C **\n", e->order);
          break;

        case 7:
          if (e->order == 1)
            {
              if (d != 0) *d = 305;
              if (r != 0) *r = 147;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 615;
              if (r != 0) *r = 129;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 887;
              if (r != 0) *r = 116;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and N **\n", e->order);
          break;

        case 8:
          if (e->order == 1)
            {
              if (d != 0) *d = 358;
              if (r != 0) *r = 143;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 799;
              if (r != 0) *r = 120;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 1072;
              if (r != 0) *r = 113;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and O **\n", e->order);
          break;

        case 9: /* C-F */
          if (d != 0) *d = 485;
          if (r != 0) *r = 135;
          break;

        case 14: /* C-Si */
          if (d != 0) *d = 318;
          if (r != 0) *r = 185;
          break;

        case 15: /* C-P */
          if (d != 0) *d = 264;
          if (r != 0) *r = 184;
          break;

        case 16: /* C-S */
          if (e->order == 1)
            {
              if (d != 0) *d = 272;
              if (r != 0) *r = 182;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 573;
              if (r != 0) *r = 160;
            }
          break;
          
        case 17: /* C-Cl */
          if (d != 0) *d = 327;
          if (r != 0) *r = 177;
          break;

        case 35: /* C-Br */
          if (d != 0) *d = 285;
          if (r != 0) *r = 194;

        case 53: /* C-I */
          if (d != 0) *d = 213;
          if (r != 0) *r = 214;
        }
      break;

    case 7:
      switch (v->atom->atno)
        {
        case 7: /* N*N */
          if (e->order == 1)
            {
              if (d != 0) *d = 167;
              if (r != 0) *r = 145;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 418;
              if (r != 0) *r = 125;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 942;
              if (r != 0) *r = 110;
            }
          break;

        case 8: /* N*O */
          if (e->order == 1)
            {
              if (d != 0) *d = 201;
              if (r != 0) *r = 140;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 607;
              if (r != 0) *r = 121;
            }
          break;

        case 9: /* N-F */
          if (d != 0) *d = 283;
          if (r != 0) *r = 136;
          break;

        case 14: /* N-Si */
          if (d != 0) *d = 355;
          if (r != 0) *r = 157.19;
          break;

        case 16: /* N=S */
          if (r != 0) *r = 149.7;
          break;
          
        case 17: /* N-Cl */
          if (d != 0) *d = 313;
          if (r != 0) *r = 175;
          break;
        }
      break;

    case 8:
      switch (v->atom->atno)
        {
        case 8:
          if (e->order == 1)
            {
              if (d != 0) *d = 142;
              if (r != 0) *r = 148;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 494;
              if (r != 0) *r = 121;
            }
          break;

        case 9:
          if (d != 0) *d = 190;
          if (r != 0) *r = 142;
          break;

        case 14:
          if (d != 0) *d = 452;
          if (r != 0) *r = 163;
          break;

        case 15:
          if (e->order == 1)
            {
              if (d != 0) *d = 335;
              if (r != 0) *r = 163;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 544;
              if (r != 0) *r = 150;
            }
          break;

        case 16:
          if (e->order == 1)
            {
              if (r != 0) *r = 157.4;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 522;
              if (r != 0) *r = 143;
            }
          break;

        case 53:
          if (d != 0) *d = 201;
          break;
        }
      break;

    case 9:
      switch (v->atom->atno)
        {
        case 9:
          if (d != 0) *d = 155;
          if (r != 0) *r = 142;
          break;
          
        case 14:
          if (d != 0) *d = 565;
          if (r != 0) *r = 160;

        case 15:
          if (d != 0) *d = 490;
          if (r != 0) *d = 154;
          break;

        case 16:
          if (d != 0) *d = 284;
          if (r != 0) *r = 156;
        }
      break;

    case 14:
      switch (v->atom->atno)
        {
        case 14:
          if (d != 0) *d = 222;
          if (r != 0) *r = 233;
          break;
          
        case 16:
          if (d != 0) *d = 293;
          if (r != 0) *r = 200;
          break;

        case 17:
          if (d != 0) *d = 381;
          if (r != 0) *r = 202;
          break;

        case 35:
          if (d != 0) *d = 310;
          if (r != 0) *r = 215;
          break;

        case 53:
          if (d != 0) *d = 234;
          if (r != 0) *r = 243;
          break;
        }
      break;
      
    case 15:
      switch (v->atom->atno)
        {
        case 15:
          if (d != 0) *d = 201;
          if (r != 0) *r = 221;
          break;
          
        case 16:
          if (e->order == 2)
            {
              if (d != 0) *d = 335;
              if (r != 0) *r = 186;
            }
          break;

        case 17:
          if (d != 0) *d = 326;
          if (r != 0) *r = 203;
          break;
          
        case 35:
          if (d != 0) *d = 264;
          break;

        case 53:
          if (d != 0) *d = 184;
        }
      break;

    case 16:
      switch (v->atom->atno)
        {
        case 16:
          if (e->order == 2)
            {
              if (d != 0) *d = 425;
              if (r != 0) *r = 149;
            }
          break;

        case 17:
          if (d != 0) *d = 255;
          if (r != 0) *r = 207;
          break;
        }
      break;

    case 17:
      if (v->atom->atno == 17)
        {
          if (d != 0) *d = 240;
          if (r != 0) *r = 199;
        }
      else if (v->atom->atno == 53)
        {
          if (d != 0) *d = 208;
          if (r != 0) *r = 232;
        }
      break;

    case 35:
      if (v->atom->atno == 35)
        {
          if (d != 0) *d = 190;
          if (r != 0) *r = 228;
        }
      else if (v->atom->atno == 53)
        {
          if (d != 0) *d = 175;
        }
      break;

    case 53:
      if (v->atom->atno == 53)
        {
          if (d != 0) *d = 148;
          if (r != 0) *r = 267;
        }
      break;
    }
}

static void
_edge_closure (int *visited, vertex_t *u, edge_t **edge,
               vertex_t **const *neighbors, graph_t *g)
{
  if (!visited[u->index])
    {
      int k;

      visited[u->index] = 1;
      for (k = 0; k < u->degree; ++k)
        {
          vertex_t *v = neighbors[u->index][k];
          if (*edge == 0 || (*edge != 0 && (*edge)->u != v))
            {
              edge_t *e = g->E;
              /* not very efficient here.. */
              while (e != 0)
                {
                  /* make sure we don't have an edge in the 
                   * opposite direction */
                  if (e->u == v && e->v == u)
                    break;
                  e = e->next;
                }
              
              if (e == 0)
                {
                  e = create_edge (++g->ne, u, v);
                  if (*edge != 0)
                    (*edge)->next = e;
                  else
                    g->E = e;
                  *edge = e;
                  
                  _edge_closure (visited, v, edge, neighbors, g);
                }
            }
        }
    }
}

static void
vertex_add_edge (vertex_t *v, edge_t *edge)
{
  int k = 0;
  for (; k < v->degree; ++k)
    if (v->edges[k] == 0)
      break;
  
  if (k == v->degree)
    fprintf (stderr, "** Error: vertex %d is already filled! **\n", v->index);
  else
    v->edges[k] = edge;
}

static void
edge_closure (vertex_t **const *neighbors, graph_t *g)
{
  int *visited = malloc ((g->nv+1)*sizeof (int));
  edge_t *edge = 0;
    
  (void) memset (visited, 0, (g->nv+1)*sizeof (int));
  _edge_closure (visited, g->V[0], &edge, neighbors, g);
  free (visited);

  edge = g->E;
  while (edge != 0)
    {
      vertex_t *u = edge->u, *v = edge->v;
      vertex_add_edge (u, edge);
      vertex_add_edge (v, edge);
      edge = edge->next;
    }
}

static void
edge_debug (const graph_t *g)
{
  edge_t *e;
  double r;
  printf ("E[%d] = \n", g->ne);
  for (e = g->E; e != 0; e = e->next)
    { vertex_t *u = e->u, *v = e->v;
      bond_energy_length (e, 0, &r);
      printf ("%d: %d %c %d => r = %.0f\n", e->index, u->index,
              e->order == 1 ? '-' : (e->order == 3 ? '#' : '='), v->index,r);
    }
}

static void
edge_parity (edge_t **edges, int *q, vertex_t *u, short *visited)
{
  int k;
  for (k = 0; k < u->degree; ++k)
    {
      edge_t *e = u->edges[k];
      if (!visited[e->index])
        {
          vertex_t *v = edge_other (e, u);
          if (implicit_hcount (u) > u->hcount
              && implicit_hcount (v) > v->hcount)
            {
              edges[(*q)++] = e;
            }
          
          visited[e->index] = 1;
          edge_parity (edges, q, v, visited);
        }
    }
}

static void
edge_parity_propagate (const edge_t **edges, int size)
{
  int i, j;
  printf ("edge parity propagate...\n");
  for (i = 0; i < size; ++i)
    {
      const edge_t *e = edges[i];
      for (j = 0; j < size; ++j)
        {
          if (i != j)
            {
              if (e->v == edges[j]->u)
                {
                  printf ("%d: %d %d --\n", edges[j]->index,
                          edges[j]->u->index, edges[j]->v->index);
                  
                  break;
                }
            }
        }
      
      if (j == size)
        {
          printf ("**\n");
          /* new path.. reset */
        }
    }
}

static int
_edge_order_assignment (vertex_t *u, short *edges)
{
  int k, h = 0;

  for (k = 0; k < u->degree; ++k)
    {
      edge_t *e = u->edges[k];
      if (!edges[e->index])
        {
          vertex_t *v = edge_other (e, u);

          edges[e->index] = e->order;
          if (implicit_hcount (u) > u->hcount)
            {
              h = implicit_hcount (v) - v->hcount;
              if (h == 0)
                {
                  h = -1;
                  continue;
                }
              
              e->order += h;
            }

          h = _edge_order_assignment (v, edges);
          if (h != 0)
            {
              /* backtracking... */
              h = implicit_hcount (v) - v->hcount;
              e->order += h;
              if (h < 0)
                {
                  /* see where else this h can go */
                  int i;
                  for (i = 0; i < u->degree; ++i)
                    {
                      e = u->edges[i];
                      if (!edges[e->index])
                        h = _edge_order_assignment (u, edges);
                    }
                }
              
              return h;
            }
        }
    }

  if (h == 0 && implicit_hcount (u) != u->hcount)
    {
      h = -1;
      /* handle limited charge */
      switch (u->atom->atno)
        {
        case 8:
          if (u->degree == 1 && u->hcount == 0)
            {
              --u->charge;
              h = 0;
            }
          break;
          
        case 7:
          if (u->degree == 4 && u->hcount == 0)
            {
              ++u->charge;
              h = 0;
            }
          break;
        }
    }
  
  return h;
}

static void
_edge_order_assignment2 (vertex_t *u, short *visited)
{
  int k, h = 0;
  
  if (implicit_hcount (u) > u->hcount)
    {
      printf ("+%d ", u->index);
      for (k = 0; k < u->degree; ++k)
        {
          edge_t *e = u->edges[k];
          if (!visited[e->index])
            {
              vertex_t *v = edge_other (e, u);
              
              visited[e->index] = 1;
              _edge_order_assignment2 (v, visited);
            }
        }
      printf ("\n");
    }
}


static void
edge_order_assignment (graph_t *g)
{
  edge_t *e;
  short *visited = malloc (sizeof (short) * (g->ne+1));
  int hcount, k = 0, i, q;

  for (k = 0; k < g->nv; ++k)
    {
      _edge_order_assignment (g->V[k], visited);
      q = 0;
      for (i = 1; i <= g->ne; ++i)
        if (visited[i]) ++q;
      if (q == g->ne)
        break; /* we're done */
    }
      
  hcount = 0;
  for (i = 0; i < g->nv; ++i)
    hcount += implicit_hcount (g->V[i]) - g->V[i]->hcount;
  printf ("V = %d => hcount = %d\n", g->V[k]->index, hcount);

  edge_debug (g);
  if (hcount > 0)
    {
      fprintf (stderr, "** %d charge(s) detected! **\n", hcount);
    }
  
  free (visited);
}

static void
create_graph_L (vertex_t **const *neighbors, graph_t *g)
{
  /* index of L and W are 1-base */
  int i, k;
  vertex_t *u, *v;
  size_t size = sizeof (double)*g->size * (g->nv+1);

  g->L = malloc (size);
  (void) memset (g->L, 0, size);

  printf ("L = [");  
  for (i = 0; i < g->nv; ++i)
    {
      u = g->V[i];
      __L(u->index, u->index) = 1.;
      for (k = 0; k < u->degree; ++k)
        {
          v = neighbors[u->index][k];
          __L(u->index, v->index) = __L(v->index, u->index)
            =  -1./sqrt ((double)u->degree * v->degree);
        }
      
      for (k = 1; k <= g->nv; ++k)
        printf (" %-4.5f", __L(i+1, k));
      
      printf (";\n");
    }
  printf ("]\n");
}

static void
create_graph_W (vertex_t **const *neighbors, graph_t *g)
{
  size_t size = sizeof (double)*g->size * (g->nv+1);
  g->W = malloc (size);
  (void) memset (g->W, 0, size);
  
}

static void
instrument_graph (graph_t *g)
{
  int i, j, count;
  formula_t *f = g->formula;
  hlayer_t *h;
  vertex_t ***neighbors;
  int *group;

  neighbors = malloc (sizeof (vertex_t **) * (g->nv+1));
  neighbors[0] = 0;
  
  g->V = malloc (sizeof (vertex_t *) * g->nv);
  /* first pass to allocate the vertices */
  for (i = 0; i < g->nv; ++i)
    {
      int d = 0;
      for (j = 0; j < g->nv; ++j)
        if (i != j && __A(i+1, j+1))
          ++d;
      neighbors[i+1] = malloc (d * sizeof (vertex_t *));
      g->V[i] = create_vertex (g, i+1, d);
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
    fprintf (stderr, "** Formula misaligned with component: "
             "expecting %d atoms but got %d! **\n", g->nv, count);
  
  printf ("graph G for component %d/%d => formula %d\n",
            g->index, g->multiplier, f->index);
    
  if (f->element->atno == 1)
    f = f->next;

  /* number of groups should never exceed number of atoms */
  group = malloc (sizeof (int)* g->nv);
  (void) memset (group, 0, sizeof (int)*g->nv);
  
  /*
   * instrument vertices
   */
  count = f->count;
  for (i = 0; i < g->nv; ++i)
    {
      int k = 0;
      vertex_t *u = g->V[i];
      
      /* now do the linking */
      for (j = 0; j < g->nv; ++j)
        if (i != j && __A(i+1, j+1))
          neighbors[u->index][k++] = g->V[j];

      u->atom = f->element;
      for (h = g->hlayer; h != 0; h = h->next)
        if (h->atom == u->index)
          {
            if (h->group > 0)
              {
#if 0
                /* shared; ensure that only the largest index get the H */
                if (h->next == 0 || h->next->group != h->group)
                  u->hcount = h->count;
#else
                if (group[h->group] == 0)
                  {
                    /* smallest index */
                    hlayer_t *n = g->hlayer;
                    while (n->group != h->group)
                      n = n->next;
                    u->hcount = n->count;
                    
                    group[h->group] = u->index;
                  }
#endif
              }
            else
              u->hcount = h->count;
            
            break;
          }
      
      if (--count == 0)
        {
          f = f->next;
          if (f != 0)
            {
              /* don't include Hs */
              if (f->element->atno == 1)
                f = f->next;
              
              if (f != 0)
                count = f->count;
            }
        }
    }

  free (group);

  /* instrument edges */
  edge_closure (neighbors, g);

  /* now adjust the order */
  edge_order_assignment (g);
  
  printf ("V[%d] = \n", g->nv);
  for (i = 0; i < g->nv; ++i)
    {
      int k;
      vertex_t *u = g->V[i];
      
      if (u->hcount > 0)
        printf ("%d %s[H%d]:", u->index, u->atom->symbol,  u->hcount);
      else
        printf ("%d %s:", u->index, u->atom->symbol);
      
      for (k = 0; k < u->degree; ++k)
        printf (" %d", neighbors[u->index][k]->index);
      
      printf ("  -- edges");
      for (k = 0; k < u->degree; ++k)
        { edge_t *e = u->edges[k];
            printf (" %d[%d,%d]", e->index, e->u->index, e->v->index);
        }
        printf ("\n");
    }
  
#if 0
  create_graph_L (neighbors, g);
  create_graph_W (neighbors, g);
#endif
  
  for (i = 1; i <= g->nv; ++i)
    if (neighbors[i] != 0)
      free (neighbors[i]);
  free (neighbors);
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
      n = -1;
    }
  else
    {
      n = g.nv;
      if (sp->bsize < g.nv)
        {
          sp->spectrum = realloc (sp->spectrum, g.nv*sizeof (float));
          sp->fiedler = realloc (sp->fiedler, g.nv*sizeof (float));
          sp->bsize = g.nv;
        }
#ifdef SPECTRAL_DEBUG
      printf ("## /c = %s\n", sp->inchi_c);
#endif
      
      instrument_graph (&g);
      if (graph_spectrum (sp->spectrum, sp->fiedler, &g) < 0)
        {
          sprintf (sp->errmsg, "Eigensolver didn't converge within "
                   "specified number of iterations");
          n = -1;
        }
    }
  destroy_graph (&g);

  return n;
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
