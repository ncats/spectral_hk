
#define __inchi_private_h__
#include "_inchi.h"

static path_t *
path_create (edge_t **edges, int ne)
{
  int i;
  path_t *p = malloc (sizeof (path_t) + sizeof (edge_t *)*ne);
  p->ne = ne;
  p->edges = (edge_t **)((unsigned char *)p + sizeof (path_t));
  for (i = 0; i < ne ; ++i)
    p->edges[i] = edges[i];
  p->prev = 0;
  p->next = 0;
  return p;
}

static void
shortest_path_dfs (path_t **path, edge_t **edges, int *ne,
                   short *visited, const vertex_t *start,
                   const vertex_t *u, const vertex_t *v,
                   const int *cost)
{
  if (u == v)
    {
      /* new path */
      path_t *p = path_create (edges, *ne);
#ifdef INCHI_DEBUG
      { int i;
        printf ("## shortest path...%d", *ne);
        for (i = 0; i < *ne; ++i)
          printf (" [%d,%d]", edges[i]->u->index, edges[i]->v->index);
        printf ("\n");
      }
#endif
      p->prev = *path;
      if (*path != 0)
        (*path)->next = p;
      *path = p;
    }
  else
    {
      int n = u->graph->nv+1;
      int k, s = start->index *n;
      
      visited[u->index] = 1;
      for (k = 0; k < u->degree; ++k)
        {
          edge_t *e = u->edges[k];
          const vertex_t *xu = _inchi_edge_other (e, u);
          if (cost[s+xu->index] + cost[xu->index*n+v->index]
              <= cost[s+v->index] && !visited[xu->index])
            {
#ifdef INCHI_DEBUG
              printf (" + %d: [%d,%d]\n", *ne, e->u->index, e->v->index);
#endif
              
              edges[(*ne)++] = e; /* push */
              shortest_path_dfs (path, edges, ne, visited, start, xu, v, cost);
              
#ifdef INCHI_DEBUG
              e = edges[*ne-1];
              printf (" - %d: [%d,%d]\n", *ne-1, e->u->index, e->v->index);
#endif
              edges[--(*ne)] = 0; /* pop */
            }
        }
      visited[u->index] = 0;
    }
}

static path_t *
shortest_paths (short *visited, const vertex_t *u,
                const vertex_t *v, const int *cost)
{
  path_t *path = 0;
  int ne = u->graph->ne;

  edge_t **edges = malloc (sizeof (edge_t *) * ne);
  (void) memset (edges, 0, sizeof (edge_t *) * ne);

#ifdef INCHI_DEBUG
  printf ("** %d - %d\n", u->index, v->index);
#endif
  ne = 0;
  shortest_path_dfs (&path, edges, &ne, visited, u, u, v, cost);

  if (path != 0)
    {
      path_t *p = path;
      while (p != 0)
        {
          path = p;
          p = p->prev;
        }
    }
  
  free (edges);
  return path;
}

static void
all_path_dfs (path_t **path, edge_t **edges, int *ne,
              short *visited, const vertex_t *start,
              const vertex_t *u, const vertex_t *v)
{
  if (u == v)
    {
      /* new path */
      path_t *p = path_create (edges, *ne);
#ifdef INCHI_DEBUG
      { int i;
        printf ("## all path...%d", *ne);
        for (i = 0; i < *ne; ++i)
          printf (" [%d,%d]", edges[i]->u->index, edges[i]->v->index);
        printf ("\n");
      }
#endif
      
      p->prev = *path;
      if (*path != 0)
        (*path)->next = p;
      *path = p;
    }
  else
    {
      int k;
      visited[u->index] = 1;
      for (k = 0; k < u->degree; ++k)
        {
          edge_t *e = u->edges[k];
          const vertex_t *xu = _inchi_edge_other (e, u);
          if (!visited[xu->index])
            {
              edges[(*ne)++] = e; /* push */
              all_path_dfs (path, edges, ne, visited, start, xu, v);
              edges[--(*ne)] = 0; /* pop */
            }
        }
      visited[u->index] = 0;
    }
}

static path_t *
all_paths (short *visited, const vertex_t *u, const vertex_t *v)
{
  path_t *path = 0;
  int ne = u->graph->ne;
  
  edge_t **edges = malloc (sizeof (edge_t *) * ne);
  (void) memset (edges, 0, sizeof (edge_t *) * ne);

  ne = 0;
  all_path_dfs (&path, edges, &ne, visited, u, u, v);

  if (path != 0)
    {
      path_t *p = path;
      while (p != 0)
        {
          path = p;
          p = p->prev;
        }
    }
  
  free (edges);
  return path;
}

static int
label_vertices (short label, short *vertex, const path_t *path)
{
  int i, ov = 0;
  const edge_t *e;
  for (i = 0; i < path->ne; ++i)
    {
      e = path->edges[i];
      if (vertex[e->u->index] && vertex[e->u->index] != label)
        ++ov;
      vertex[e->u->index] = label;
      if (vertex[e->v->index] && vertex[e->v->index] != label)
        ++ov;
      vertex[e->v->index] = label;
    }
  return ov;
}

/* determine whether q is a subpath of p */
static int
path_contains_path (const path_t *p, const path_t *q)
{
  int i, j;
  const edge_t *e;

  if (q->ne > p->ne) return 0;

  /* sigh.. */
  for (i = 0; i < q->ne; ++i)
    {
      e = q->edges[i];
      for (j = 0; j < p->ne; ++j)
        if (e->index == p->edges[j]->index)
          break;
      
      if (j == p->ne)
        return 0;
    }
  return 1;
}

static int
path_contains_edge (const path_t *p, const edge_t *e)
{
  int i;
  for (i = 0; i < p->ne; ++i)
    if (e->index == p->edges[i]->index)
      return 1;
  return 0;
}

static int
paths_overlap_at_vertex (short *vertices, const path_t *p,
                         const path_t *q, const vertex_t *u)
{
  return (0 == label_vertices (1, vertices, p)
          && 1 == label_vertices (2, vertices, q) /* one vertex overlap */
          && 2 == vertices[u->index]); /* and it's u */
}

int
_inchi_ring_perception (inchi_t *g)
{
  edge_t *e;
  path_t *p, *q, *next;
  
  int n = g->nv+1, i, j, k;
  int *cost = malloc (n*n*sizeof (int));
  short *visited = malloc (n * sizeof (short));

  (void) memset (cost, 0, sizeof (int)*n*n);
  
  for (i = 0; i < n; ++i)
    {
      k = i*n;
      for (j = i+1; j < n; ++j)
        cost[k+j] = cost[j*n+i] = n;
    }
  for (e = g->E; e != 0; e = e->next)
    {
      i = e->u->index;
      j = e->v->index;
      cost[i*n+j] = cost[j*n+i] = 1;
    }

  /* 1. now perform floyd-warshall's all pairs shortest path */
  for (k = 1; k < n; ++k)
    for (i = 1; i < n; ++i)
      for (j = 1; j < n; ++j)
        cost[i*n+j] = MIN (cost[i*n+j], cost[i*n+k] + cost[k*n+j]);

  /* 2. find all basis cycles */
  for (i = 0; i < g->nv; ++i)
    {
      vertex_t *u = g->V[i];
      for (e = g->E; e != 0; e = e->next)
        {
          path_t *p1;
          
          (void) memset (visited, 0, sizeof (short) * n);
          p1 = shortest_paths (visited, u, e->u, cost);
          for (p = p1; p != 0;)
            {
              path_t *p2;
              
              label_vertices (1, visited, p);
              p2 = all_paths (visited, u, e->v);
              for (q = p2; q != 0;)
                {
                  /* see if p & q overlap at u */
                  if (paths_overlap_at_vertex (visited, p, q, u))
                    {
                      int found = 0;
                      path_t *r = g->R;
                      while (r != 0)
                        {
                          if (path_contains_path (r, p)
                              && path_contains_path (r, q)
                              && path_contains_edge (r, e))
                            {
                              found = 1;
                              break;
                            }
                          r = r->prev;
                        }

                      if (!found)
                        {
                          int size = p->ne + q->ne + 1;
                          path_t *ring =
                            malloc (sizeof (path_t) + size*sizeof (edge_t *));

                          ring->ne = size;
                          ring->edges =
                            (edge_t **)((unsigned char *)ring
                                        + sizeof (path_t));
                      
                          for (j = 0; j < p->ne; ++j)
                            ring->edges[j] = p->edges[j];
                          for (k = 0; k < q->ne; ++k, ++j)
                            ring->edges[j] = q->edges[k];
                          
                          ring->edges[j] = e;
                          ring->next = 0;
                          ring->prev = g->R;
                          if (g->R != 0)
                            g->R->next = ring;
                          g->R = ring;
                          
#ifdef INCHI_DEBUG
                          printf (" ** => new ring size %d...", ring->ne);
#endif
                          for (j = 0; j < ring->ne; ++j)
                            {
                              _inchi_vertex_set (ring->edges[j]->u, FLAG_RING);
                              _inchi_vertex_set (ring->edges[j]->v, FLAG_RING);
                              
#ifdef INCHI_DEBUG
                              printf (" [%d,%d]", ring->edges[j]->u->index,
                                      ring->edges[j]->v->index);
#endif
                            }
#ifdef INCHI_DEBUG
                          printf ("\n");
#endif
                          ++g->nr;
                        } /* !found */
                    } /* paths overlap */
                  
                  next = q->next;
                  free (q);
                  q = next;
                }
              
              next = p->next;
              free (p);
              p = next;
            }
        }
    }

  for (q = g->R; q != 0; q = q->prev)
    g->R = q;
  
  free (cost);
  free (visited);
  
  return g->nr;
}
