#ifdef __inchi_private_h__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "inchi.h"
#include "periodic.h"

typedef struct __formula_s {
  int index; /* this is not the same the component index */
  short count;
  const element_t *element;
  struct __formula_s *next;
} formula_t;

typedef struct __hlayer_s {
  int index; /* component index */
  int atom; /* atom index */
  short count; /* number of Hs */
  short group; /* shared Hs if group > 1 */
  const element_t *element;
  struct __hlayer_s *next;
} hlayer_t;

typedef struct __vertex_s {
  int index;
  short degree;
  short charge;
  
  short hcount;
  short hshare; /* number of shared h's (if any) for this group */
  short hgroup; /* sharing group for h if != 0 */

#define FLAG_RING    0
#define FLAG_HETERO  1
#define FLAG_CHIRAL  2
#define FLAG_QUERY   4
  unsigned flags;
  
  inchi_t *graph;  
  const element_t *atom;
  struct __edge_s **edges; /* edges[0..degree-1] */
} vertex_t;

typedef struct __edge_s {
  int index;
  short order;
  vertex_t *u;
  vertex_t *v;
  struct __edge_s *next;
} edge_t;

typedef struct __path_s {
  int ne;
  edge_t **edges;
  struct __path_s *prev;  
  struct __path_s *next;
} path_t;

struct __inchi_s {
  int index; /* component index */
  short multiplier; /* multiplicity */
  
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
  
  size_t size; /* stride size for A, L, W */

  formula_t *formula;
  hlayer_t *hlayer;

  short nr;
  path_t *R; /* rings.. R[0..nr-1] */

  char *inchi_c;
  size_t isize;
  
  char errmsg[BUFSIZ];
};

#ifndef MIN
# define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
# define MAX(a,b) ((a)>(b)?(a):(b))
#endif

extern vertex_t *_inchi_edge_other (const edge_t *e, const vertex_t *v);
extern void bond_energy_length (const edge_t *e, double *d, double *r);
extern int _inchi_ring_perception (inchi_t *g);
extern void _inchi_vertex_set (vertex_t *v, unsigned flag);
extern int _inchi_vertex_get (const vertex_t *v, unsigned flag);

#endif /* !__inchi_private_h__ */
