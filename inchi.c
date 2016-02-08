
#define __inchi_private_h__
#include "_inchi.h"

#define __set_edge(i,j) G[(i)*size+(j)] = G[(j)*size+(i)] = 1
#define __get_edge(i,j) G[(i)*size+(j)]

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

#undef __set_edge
#undef __get_edge

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
create_vertex (inchi_t *g, int index, int degree)
{
  vertex_t *v = malloc (sizeof (vertex_t) + sizeof (edge_t*)*degree);
  v->index = index;
  v->degree = degree;
  v->charge = 0;
  v->atom = 0;
  v->hcount = 0;
  v->hgroup = 0;
  v->hshare = 0;
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

vertex_t *
_inchi_edge_other (const edge_t *e, const vertex_t *v)
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
      if (u->degree <= 2) v = 2 - v;
      else if (u->degree < 4) v = 4 - v;
      else v = 6 - v;
      break;
    case 17: v = 1 - v; break;
    case 35: case 53: v = 1 - v; break;
    default: v = -1000;
    }
  
  return v;
}

static void
_edge_closure (int *visited, vertex_t *u, edge_t **edge,
               vertex_t **const *neighbors, inchi_t *g)
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

static edge_t *
vertex_get_edge (vertex_t *u, vertex_t *v)
{
  int i;
  for (i = 0; i < u->degree; ++i)
    if (v == _inchi_edge_other (u->edges[i], u))
      return u->edges[i];
  
  return 0;
}

static void
edge_closure (vertex_t **const *neighbors, inchi_t *g)
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
edge_debug (const inchi_t *g)
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

static int
sort_vertex (const void *p1, const void *p2)
{
  const vertex_t *v1 = *(vertex_t * const *)p1;
  const vertex_t *v2 = *(vertex_t * const *)p2;
  int d = v2->atom->atno - v1->atom->atno;
  if (d == 0)
    d = v1->index - v2->index; /* order by index */
  return d;
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
          vertex_t *v = _inchi_edge_other (e, u);

          edges[e->index] = e->order;
          /* not shared h in v */
          if (u->hgroup == 0 && v->hgroup == 0
              && implicit_hcount (u) > u->hcount)
            {
              h = implicit_hcount (v) - v->hcount;
              if (h == 0)
                {
                  /* *-[*](-*)(-*)(-*) */
                  if (u->hcount == 0 && u->degree == 4)
                    h = 0;
                  else
                    {
                      h = -1;
                      continue;
                    }
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
      /* see if any of the neighboring atom contains shared h that
       * we can use
       */
      if (u->hgroup == 0)
        {
          vertex_t *nb[10] = {0};

          assert (u->degree < 10);
          for (k = 0; k < u->degree; ++k)
            nb[k] = _inchi_edge_other (u->edges[k], u);
          
          qsort (nb, u->degree, sizeof (nb[0]), sort_vertex);
          
          h = implicit_hcount (u) - u->hcount;
          for (k = 0; k < u->degree; ++k)
            {
              vertex_t *v = nb[k];
              printf ("%d: %d %d %d/%d %d\n", u->index, v->index,
                      v->hgroup, v->hcount, v->hshare, v->atom->atno);
              if (v->hshare >= h)
                {
                  edge_t *e = vertex_get_edge (u, v);
                  if (e != 0)
                    {
                      e->order += h;
                      h = 0;
                    }
                  break;
                }
            }
        }

      /* if thing still isn't khosher, we try to see if it might be due
       * to charge 
       */
      if (h != 0)
        {
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
    }
  
  return h;
}

static void
edge_order_assignment (inchi_t *g)
{
  vertex_t *v;
  hlayer_t *h;
  short *visited = malloc (sizeof (short) * (g->ne+1));
  int k = 0, i, q;

  (void) memset (visited, 0, sizeof (short) * (g->ne+1));
  for (k = 0; k < g->nv; ++k)
    {
      v = g->V[k];
      _edge_order_assignment (v, visited);
      
      q = 0;
      for (i = 1; i <= g->ne; ++i)
        if (visited[i]) ++q;

      if (q == g->ne)
        break; /* we're done */
    }
      
  q = 0;
  for (i = 0; i < g->nv; ++i)
    {
      v = g->V[i];
      q += implicit_hcount (v) - v->hcount;
    }
  printf ("V = %d => hcount = %d\n", g->V[k]->index, q);

  edge_debug (g);
  if (q > 0)
    fprintf (stderr, "** %d charge/tautomer(s) detected! **\n", q);
  
  free (visited);
}

static void
create_graph_L (vertex_t **const *neighbors, inchi_t *g)
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
create_graph_W (vertex_t **const *neighbors, inchi_t *g)
{
  size_t size = sizeof (double)*g->size * (g->nv+1);
  g->W = malloc (size);
  (void) memset (g->W, 0, size);
  
}

static void
instrument_graph (inchi_t *g)
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
              u->hshare = h->count;
            else
              u->hcount = h->count;
            u->hgroup = h->group;
            
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

  _inchi_ring_perception (g);
  { path_t *r = g->R;
    i = j = 0;
    printf ("@@ %d ring(s)...\n", g->nr);
    for (; r != 0; r = r->next)
      {
        printf ("  %d:", i);
        for (j = 0; j < r->ne; ++j)
          printf (" [%d,%d]", r->edges[j]->u->index, r->edges[j]->v->index);
        printf ("\n");
        ++i;
      }
    printf ("\n");
  }

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


inchi_t *
inchi_create ()
{
  static inchi_t def = {0};
  inchi_t *g = malloc (sizeof (inchi_t));
  (void) memcpy (g, &def, sizeof (def));
  
  return g;
}

static void
inchi_destroy (inchi_t *g)
{
  if (g->A != 0)
    {
      free (g->A);
      g->A = 0;
    }
  
  if (g->V != 0)
    {
      int i;
      for (i = 0; i < g->nv; ++i)
        destroy_vertex (g->V[i]);
      free (g->V);
      g->V = 0;
    }
  
  if (g->L != 0)
    {
      free (g->L);
      g->L = 0;
    }

  if (g->W != 0)
    {
      free (g->W);
      g->W = 0;
    }

  if (g->R != 0)
    {
      path_t *p = g->R, *n;
      while (p != 0)
        {
          n = p->next;
          free (p);
          p = n;
        }
      g->R = 0;
    }

  if (g->inchi_c != 0)
    {
      free (g->inchi_c);
      g->inchi_c = 0;
    }
  
  destroy_edge (g->E);
  g->E = 0;
  
  destroy_formula (g->formula);
  g->formula = 0;
  
  destroy_hlayer (g->hlayer);
  g->hlayer = 0;
}

void
inchi_free (inchi_t *g)
{
  if (g != 0)
    {
      inchi_destroy (g);
      free (g);
    }
}

int
inchi_parse (inchi_t *g, const char *inchi)
{
  char *ptr, *start, *end;
  int *pG = 0, n;
  size_t size = 0, psize = 0;

  if (strncmp ("InChI=", inchi, 6) != 0)
    {
      sprintf (g->errmsg, "Inchi string doesn't begin with InChI=");
      return -1;
    }

  start = strstr (inchi, "/c");
  if (start == 0)
    {
      sprintf (g->errmsg, "InChI string doesn't have connection layer");
      return 0;
    }

  /* clean up old */
  inchi_destroy (g);

  /* parse formula */
  n = parse_formula (&g->formula, g->errmsg, inchi);
  { formula_t *formula = g->formula;
    printf ("formula: %d\n", n);
    while (formula != 0)
      {
        printf ("%d: %d %s\n", formula->index,
                formula->count, formula->element->symbol);
        formula = formula->next;
      }
  }
  
  /* parse h layer */
  n = parse_layer_h (&g->hlayer, g->errmsg, inchi);
  { hlayer_t *h = g->hlayer;
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
      n = parse_inchi_graph (&pG, &psize, ptr, g->errmsg);
      if (n > g->nv)
        {
          /* keep only the largest component */
          if (g->A != 0)
            free (g->A);

          g->A = pG;
          g->nv = n;
          g->size = psize;
          g->index = 0;

          { char *p = ptr;
            for (; p >= start; --p)
              if (*p == ';' || *p == '\0')
                ++g->index;

            p = ptr;
            g->multiplier = strtol (p, &p, 10);
            if (*p != '*')
              g->multiplier = 1;
            else
              /* don't retain the multiplicity character from the /c layer */
              ptr = ++p;
          }

          psize = strlen (ptr);
          if (psize > g->isize)
            {
              g->inchi_c = realloc (g->inchi_c, psize+1);
              g->isize = psize;
            }
          (void) strcpy (g->inchi_c, ptr);

          pG = 0;
          psize = 0;
        }
    }
  free (start);

  if (pG != 0)
    free (pG);

  instrument_graph (g);
  
  return n;
}

int
inchi_node_count (const inchi_t *g)
{
  return g->nv;
}

int
inchi_edge_count (const inchi_t *g)
{
  return g->ne;
}

const char *
inchi_layer_c (const inchi_t *g)
{
  return g->inchi_c;
}

const int *
inchi_matrix_A (const inchi_t *g)
{
  return g->A;
}

size_t
inchi_matrix_size (const inchi_t *g)
{
  return g->size;
}

const char *
inchi_error (const inchi_t *g)
{
  return g->errmsg;
}
