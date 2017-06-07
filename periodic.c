
#include <stdio.h>
#include <string.h>

#include "periodic.h"

/*
 * only defined for the organic subset
 */
static struct __element_s TABLE[] = {
  {"Si", 14, 28.0844, 1.8, 1.39, 8.15},  
  {"Cl", 17, 35.453, 3.0, 3.61, 12.97},
  {"Ca", 20, 40.078, 1.0, -0.19, 6.11},
  {"Fe", 26, 55.845, 1.8, 0.50, 7.90},
  {"Co", 27, 58.933195, 1.8, 0.66, 7.86},
  {"Ge", 32, 72.64, 1.8, 0, 0}, /* FIXME: last two numbers are placeholder */
  {"Se", 34, 78.96, 2.4, 2.02, 9.75},
  {"Br", 35, 79.904, 2.8, 3.45, 11.81},
  {"Te", 54, 131.293, 0.0, 1.97, 9.01},
  {"Ru", 44, 101.07, 2.2, 1.10, 7.37},
  {"Pt", 78, 195.084, 2.2, 1.10, 8.96},  
  {"Pb", 82, 207.2, 1.8, 1.39, 7.42},
  {"Hg", 86, 222.0, 0.0, -0.19, 10.44},
  
  {"H", 1, 1.00794, 2.1, 0.75, 13.6},
  {"B", 5, 10.811, 2.0, 0.28, 8.30},
  {"C", 6, 12.0107, 2.5, 1.26, 11.26},
  {"N", 7, 14.0067, 3.0, 0.44, 14.53},
  {"O", 8, 15.9994, 3.5, 1.46, 13.62},
  {"F", 9, 18.9984, 4.0, 3.45, 17.42},
  {"P", 15, 30.9737, 2.1, 0.75, 10.49},
  {"S", 16, 32.065, 2.5, 2.0, 10.36},
  {"I", 53, 126.90447, 3.23, 10.45}
};

static size_t SIZE = sizeof (TABLE) / sizeof (TABLE[0]);

static int
cmp (const char *p1, const char *p2)
{
  for (; *p1 != '\0' && *p2 != '\0' && *p1 == *p2; ++p1, ++p2)
    ;
  return  (*p1 == '\0' || *p2 == '\0') ? 0 : *p1 - *p2;
}

const element_t *
element_lookup_symbol (const char *symbol)
{
  size_t i;
  for (i = 0; i < SIZE; ++i)
    if (cmp (TABLE[i].symbol, symbol) == 0)
      return &TABLE[i];
  return 0;
}

const element_t *
element_lookup_atno (int atno)
{
  size_t i;
  for (i = 0; i < SIZE; ++i)
    if (TABLE[i].atno == atno)
      return &TABLE[i];
  return 0;
}
