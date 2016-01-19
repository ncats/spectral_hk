/**
 * Test driver for libspectral.a
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "spectral.h"


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
                    
              fprintf (outfp, "%s\t%s\t", hk, buffer);
              { size_t i, size = spectral_size (spectral);
                const float *v = spectral_spectrum (spectral);
                for (i = 0; i < size; ++i)
                  {
                    (void) fprintf (outfp, "%.5f", v[i]);
                    if (i+1 < size)
                      (void) fprintf (outfp, ",");
                  }
                fprintf (outfp, "\t");
                v = spectral_fiedler (spectral);
                for (i = 0; i < size; ++i)
                  {
                    (void) fprintf (outfp, "%.5f", v[i]);
                    if (i+1 < size)
                      (void) fprintf (outfp, ",");
                  }
                fprintf (outfp, "\n");
              }
            }
          else
            fprintf (stderr, "error: ** failed to process %s (%s) **\n", 
                     tok, spectral_error (spectral));
        }
      else
        fprintf (stderr, "error: ** InChI string exceeds buffer size! **\n");
    }

  spectral_free (spectral);

  if (infp != stdin)
    fclose (infp);

  if (outfp != stdout)
    fclose (outfp);

  return 0;
}

