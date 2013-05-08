/**
 * Test driver for libspectral.a
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "spectral.h"


int
main (int argc, char *argv[])
{
  /* decode inchi graph */
  FILE *infp, *outfp;
  spectral_t *spectral = spectral_create ();
  char buffer[BUFSIZ];
  const char *hk;

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
      char *tok, *last;
      tok = strtok_r (buffer, " ", &last);
      hk = spectral_digest (spectral, buffer);
      if (hk != 0)
        {
          printf ("%s\t%s", hk, tok);
        }
      else
        fprintf (stderr, "error: ** failed to process %s (%s) **\n", 
                 tok, spectral_error (spectral));
    }

  spectral_free (spectral);

  if (infp != stdin)
    fclose (infp);

  if (outfp != stdout)
    fclose (outfp);

  return 0;
}

