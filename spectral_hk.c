/**
 * Test driver for libspectral.a
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "spectral.h"

typedef struct eigenstats_s {
  double sq_power;
  double mean_sq_power;
  double mean;
  double var;
  double median;
  double skewness;
  double kurtosis;
} eigenstats_t;

/*
 * s is sorted
 */
static void calc_stats (eigenstats_t *stats, const float *s, size_t len)
{
  size_t i;
  double v;
  (void) memset (stats, 0, sizeof (*stats));
  for (i = 0; i < len; ++i)
    {
      stats->mean += s[i];
      stats->sq_power += s[i]* s[i];
    }
  stats->mean /= len; /* always 1 for normalized laplacian */
  stats->mean_sq_power = stats->sq_power / len;
  if (len % 2 == 0)
    stats->median = (s[len/2-1] + s[len/2])/2.;
  else
    stats->median = s[len/2];
  /* second pass */
  for (i = 0; i < len; ++i)
    {
      v = s[i] - stats->mean;
      stats->var += v * v;
      stats->skewness += v*v*v;
    }
  v = stats->skewness;
  stats->skewness = (v/len) / pow(stats->var/(len-1.0), 1.5);
  stats->var /= len;
}

int
main (int argc, char *argv[])
{
  /* decode inchi graph */
  FILE *infp, *outfp;
  spectral_t *spectral = spectral_create ();
  char buffer[1<<15] = {0};
  char inchi[1<<14] = {0};
  char *end = buffer + sizeof (buffer);
  const char *hk;

  fprintf (stderr, "## spectral_hk -- %s\n", spectral_version ());
  if (argc > 1)
    {
      infp = fopen (argv[1], "r");
      if (infp == 0)
        {
          fprintf (stderr, "** error: can't open file '%s' for reading! **",
                   argv[1]);
          return 1;
        }
    }
  else
    infp = stdin;

  if (argc > 2)
    {
      outfp = fopen (argv[2], "w");
      if (outfp == 0)
        {
          fprintf (stderr, "** error: can't open file '%s' for writing! **",
                   argv[2]);
          return 1;
        }
    }
  else
    outfp = stdout;
  
  /*
   * assume line contains INCHI as the first token
   */
  while (fgets (buffer, sizeof (buffer), infp) != 0)
    {
      char *tok = buffer;
      while (!isspace (*tok) && tok < end)
        ++tok;

      if (tok < end && (tok - buffer) < sizeof (inchi))
        {
          (void) strncpy (inchi, buffer, tok - buffer);
          inchi[tok-buffer] = '\0';

          hk = spectral_digest (spectral, inchi);
          if (hk != 0)
            {
              tok = buffer + strlen (buffer);
              while (--tok > buffer && *tok != '\n')
                ;
              *tok = '\0';
                    
              { size_t i, size = spectral_size (spectral);
                const float *v = spectral_spectrum (spectral);
                
                (void) fprintf (outfp, "%s\t%s\t%ld\t", hk, buffer, size);
                for (i = 0; i < size; ++i)
                  {
                    (void) fprintf (outfp, "%.5f", v[i]);
                    if (i+1 < size)
                      (void) fprintf (outfp, ",");
                  }
                
#ifdef FIEDLER_VECTOR
                (void) fprintf (outfp, "\t");
                v = spectral_fiedler (spectral);
                for (i = 0; i < size; ++i)
                  {
                    (void) fprintf (outfp, "%.5e", v[i]);
                    if (i+1 < size)
                      (void) fprintf (outfp, ",");
                  }
#endif
#if 0                
                (void) fprintf (outfp, "\t");           
                { eigenstats_t stats;
                  calc_stats (&stats, spectral_spectrum (spectral), size);
                  (void) fprintf (outfp, "%.5f,", stats.sq_power);
                  (void) fprintf (outfp, "%.5f,", stats.mean_sq_power);
                  (void) fprintf (outfp, "%.5f,", stats.mean);
                  (void) fprintf (outfp, "%.5f,", stats.var);
                  (void) fprintf (outfp, "%.5f,", stats.median);
                  (void) fprintf (outfp, "%.5e", stats.skewness);
                }
#endif
                (void) fprintf (outfp, "\n");
              }
            }
          else
            (void) fprintf (stderr, "error: ** failed to process %s (%s) **\n", 
                            tok, spectral_error (spectral));
        }
      else
        (void) fprintf (stderr,
                        "error: ** InChI string exceeds buffer size! **\n");
    }

  spectral_free (spectral);

  if (infp != stdin)
    (void) fclose (infp);

  if (outfp != stdout)
    (void) fclose (outfp);

  return 0;
}

