
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sha1.h"

int
main (int argc, char *argv[])
{
  int i;
  unsigned char digest[20];
  sha1_state *sha1 = sha1_create ();
  char *s = "abc";

  sha1_update (sha1, (unsigned char *) "abc", strlen (s));
  sha1_digest (sha1, digest);

  for (i = 0; i < 20; ++i)
    {
      printf ("%02x", digest[i]);
      if ((i+1) % 4 == 0) printf (" ");
    }
  printf ("%s\n", s);

  s = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  sha1_reset (sha1);
  sha1_update (sha1, (unsigned char *) s, strlen (s));
  sha1_digest (sha1, digest);
  for (i = 0; i < 20; ++i)
    {
      printf ("%02x", digest[i]);
      if ((i+1) % 4 == 0) printf (" ");
    }
  printf ("%s\n", s);

  sha1_free (sha1);
  return 0;
}
