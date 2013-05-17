/**
 * basic interval tree to encode numerical values
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "interval.h"

struct __interval_s 
{
  int id;
  double start;
  double end;
  struct __interval_s *left;
  struct __interval_s *right;
};

static int
is_leaf (const interval_t *n)
{
  return n->left == 0 && n->right == 0;
}

static interval_t *
build_interval (double start, double end, int *id, int depth, int bits)
{
  interval_t *node = 0;

  if ((1<<depth) <= bits)
    {
      node = malloc (sizeof (struct __interval_s));
      if (node != 0)
        {
          double split = start + (end - start)/2.;
          node->start = start;
          node->end = end;
          node->id = -1;
          node->left = build_interval (start, split, id, depth+1, bits);
          node->right = build_interval (split, end, id, depth+1, bits);
          if (is_leaf (node))
            node->id = (*id)++;
        }
    }

  return node;
}


interval_t *
interval_create (double start, double end, int bits)
{
  int id = 0;
  return build_interval (start, end, &id, 0, bits);
}

static void
get_range (const interval_t *node, double *range)
{
  if (is_leaf (node))
    *range = node->end - node->start;
  else
    get_range (node->left, range); /* we only need one leaf */
}

double
interval_range (const interval_t *tree)
{
  double range = 0;
  get_range (tree, &range);
  return range;
}

static void
dump_interval_tree (const interval_t *node, FILE *outfp, int depth)
{
  int i;
  if (node != 0)
    {
      for (i = 0; i <= depth; ++i)
        fprintf (outfp, "  ");
      if (is_leaf (node))
        fprintf (outfp, "%5d: [%.5lf, %.5lf)\n", 
                 node->id, node->start, node->end);
      else
        fprintf (outfp, "    [%.5lf, %.5lf)\n", node->start, node->end);
      dump_interval_tree (node->left, outfp, depth+1);
      dump_interval_tree (node->right, outfp, depth+1);
    }
}

void
interval_dump (const interval_t *tree, FILE *outfp)
{
  dump_interval_tree (tree, outfp, 0);
}

void
interval_free (interval_t *inter)
{
  if (inter != 0)
    {
      interval_free (inter->left);
      interval_free (inter->right);
      free (inter);
    }
}

static void
encode32 (const interval_t *node, uint32_t *u32, double value, double slack)
{
  if (node != 0)
    {
      double split = node->start + (node->end-node->start)/2.;
      if (value <= split+slack)
        {
#ifdef INTERVAL_DEBUG
          if (is_leaf (node)) 
            printf ("%d: [%lf,%lf) left\n", 
                    node->id, node->start, split+slack);
          else
            printf ("   [%lf,%lf) left\n", node->start, split+slack);
#endif
          encode32 (node->left, u32, value, slack);
          if (node->id >= 0)
            *u32 |= 1<<node->id;
        }

      if (value >= split-slack)
        {
#ifdef INTERVAL_DEBUG
          if (is_leaf (node)) 
            printf ("%d: [%lf,%lf) left\n", 
                    node->id, node->start, split+slack);
          else
            printf ("   [%lf,%lf) left\n", node->start, split+slack);
#endif
          encode32 (node->right, u32, value, slack);
          if (node->id >= 0)
            *u32 |= 1<<node->id;
        }
    }
}

uint32_t 
interval_encode32 (const interval_t *tree, double value, double slack)
{
  uint32_t u32 = 0;
#ifdef INTERVAL_DEBUG
  printf ("====> %lf +/- %lf [%lf,%lf]\n", 
          value, slack, value-slack,value+slack);
#endif
  encode32 (tree, &u32, value, slack);
  return u32;
}

#ifdef __INTERVAL_TEST
int main (int argc, char *argv[])
{
  int i;
  double v;

  interval_t *tree = interval_create (0., 1., 32);
  printf ("** interval range: %lf\n", interval_range (tree));
  interval_dump (tree, stdout);

  for (i = 1; i < argc; ++i)
    {
      v = atof (argv[i]);
      printf ("%.10lf => %d\n", v, 
              interval_encode32 (tree, v, 0.00005));
    }
  
  interval_free (tree);

  return 0;
}
#endif
/**
 * Local Variables:
 * compile-command: "gcc -Wall -g -o interval interval.c -D__INTERVAL_TEST"
 * End:
 */
