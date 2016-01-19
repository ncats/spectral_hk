#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "spectral.h"

int
main (int argc, char *argv[])
{
  FILE *infp, *outfp;
  spectral_t *spectral = spectral_create ();
  char buffer[1<<15];
  char inchi[1<<14];
  char *end = buffer + sizeof (buffer);
  double ratio;

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

          ratio = 0.;
          if (spectral_ratio (&ratio, spectral, inchi) == 0)
            {
              double d = fabs (2*M_PI-ratio);
              //if (d < 0.001)
                fprintf (outfp, "%.5f\t%s", ratio, buffer);
            }
          else
            fprintf (stderr, "error: ** can't calculate spectral ratio! **\n");
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
