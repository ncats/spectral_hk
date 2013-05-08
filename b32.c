
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "b32.h"

/*
 * 2^15 = 32768 = 32^3
 * 
 * A triple in base-32 represents 15 bits.
 */

/*
 * The alphabet is choosen such that common vowels are removed
 */
static const char ALPHA[] = {
  'A','B','C','D','F','G','H','J','K',
  'L','M','N','P','Q','R','S','T','U',
  'V','W','X','Y','Z','1','2','3','4',
  '5','6','7','8','9',
};

static size_t ALPHA_SIZE = sizeof (ALPHA) / sizeof (char);


char *
b32_unrank (char **output, int64_t rank, size_t size)
{
  size_t r, sz;
  char *pbuf;

  /*
   * max decoding size is 60 = 5 x 12
   */
  if (size > 60) size = 60; 
  sz = (size + 4) / 5;

  if (*output == 0)
    {
      *output = pbuf = malloc ((sz+1) * sizeof (char));
      if (pbuf == 0)
        return pbuf;
    }
  else
    pbuf = *output;
  pbuf += sz;
  *pbuf = '\0';

  for (r = 0; r < size; r += 5)
    {
      *--pbuf = ALPHA[rank & 0x1f];
      rank >>= 5;
    }

  return pbuf;
}

int64_t
b32_rank (const char *pbuf)
{
  int64_t rank = 0, i;
  int r = 0;
  const char *p = pbuf + strlen (pbuf);

  for (;--p >= pbuf; r += 5)
    {
      for (i = 0; i < ALPHA_SIZE; ++i)
        if (ALPHA[i] == *p)
          break;
      if (i == ALPHA_SIZE)
        return -1;
      rank |= i << r;
    }
  
  return rank;
}

#define __b32_b3(d) (((d)[0]&0xff)|(((d)[1]&0x7f)<<8)) /*8+7=15[1]*/
#define __b32_b4(d) ((((d)[1]&0x80)>>7)|(((d)[2]&0x0f)<<1)) /*1+4=20[4]*/
#define __b32_b5(d) ((((d)[2]&0xf0)>>4)|(((d)[3]&0x01)<<4)) /*4+1=25[7]*/
#define __b32_b6(d) (((d)[3]&0x3f)>>1) /*0+5=30[2]*/
#define __b32_b7(d) ((((d)[3]&0xc0)>>6)|(((d)[4]&0x07)<<2)) /*2+3=35[5]*/
#define __b32_b8(d) (((d)[4]&0xf8)>>3) /*5+0=40[0]*/
#define __b32_b9(d) ((d)[5]&0x1f) /*0+5=45[3]*/
#define __b32_b10(d) ((((d)[6]&0x03)<<3)|(((d)[5]&0xe0)>>5)) /*2+3=50[6]*/
#define __b32_b11(d) (((d)[6]&0x7f)>>2) /*5+0=55[1]*/
#define __b32_b12(d) ((((d)[7]&0x0f)<<1)|(((d)[6]&0x80)>>7)) /*4+1=60[4]*/

char *
b32_encode15 (char **output, const unsigned char *data, size_t len)
{
  if (len < 2) return 0;
  return b32_unrank (output, __b32_b3 (data), 15);
}

char *
b32_encode20 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;  
  if (len >= 3)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 20);
      pbuf[0] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode25 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 4)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 25);
      pbuf[0] = ALPHA[__b32_b5 (data)];
      pbuf[1] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode30 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 4)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 30);
      pbuf[0] = ALPHA[__b32_b6 (data)];
      pbuf[1] = ALPHA[__b32_b5 (data)];
      pbuf[2] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode35 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 5)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 35);
      pbuf[0] = ALPHA[__b32_b7 (data)];
      pbuf[1] = ALPHA[__b32_b6 (data)];
      pbuf[2] = ALPHA[__b32_b5 (data)];
      pbuf[3] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode40 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 5)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 40);
      pbuf[0] = ALPHA[__b32_b8 (data)];
      pbuf[1] = ALPHA[__b32_b7 (data)];
      pbuf[2] = ALPHA[__b32_b6 (data)];
      pbuf[3] = ALPHA[__b32_b5 (data)];
      pbuf[4] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode45 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 6)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 45);
      pbuf[0] = ALPHA[__b32_b9 (data)];
      pbuf[1] = ALPHA[__b32_b8 (data)];
      pbuf[2] = ALPHA[__b32_b7 (data)];
      pbuf[3] = ALPHA[__b32_b6 (data)];
      pbuf[4] = ALPHA[__b32_b5 (data)];
      pbuf[5] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode50 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 7)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 50);
      pbuf[0] = ALPHA[__b32_b10(data)];
      pbuf[1] = ALPHA[__b32_b9 (data)];
      pbuf[2] = ALPHA[__b32_b8 (data)];
      pbuf[3] = ALPHA[__b32_b7 (data)];
      pbuf[4] = ALPHA[__b32_b6 (data)];
      pbuf[5] = ALPHA[__b32_b5 (data)];
      pbuf[6] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode55 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 7)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 55);
      pbuf[0] = ALPHA[__b32_b11(data)];
      pbuf[1] = ALPHA[__b32_b10(data)];
      pbuf[2] = ALPHA[__b32_b9 (data)];
      pbuf[3] = ALPHA[__b32_b8 (data)];
      pbuf[4] = ALPHA[__b32_b7 (data)];
      pbuf[5] = ALPHA[__b32_b6 (data)];
      pbuf[6] = ALPHA[__b32_b5 (data)];
      pbuf[7] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}

char *
b32_encode60 (char **output, const unsigned char *data, size_t len)
{
  char *pbuf = 0;
  if (len >= 8)
    {
      pbuf = b32_unrank (output, __b32_b3 (data), 60);
      pbuf[0] = ALPHA[__b32_b12(data)];
      pbuf[1] = ALPHA[__b32_b11(data)];
      pbuf[2] = ALPHA[__b32_b10(data)];
      pbuf[3] = ALPHA[__b32_b9 (data)];
      pbuf[4] = ALPHA[__b32_b8 (data)];
      pbuf[5] = ALPHA[__b32_b7 (data)];
      pbuf[6] = ALPHA[__b32_b6 (data)];
      pbuf[7] = ALPHA[__b32_b5 (data)];
      pbuf[8] = ALPHA[__b32_b4 (data)];
    }
  return pbuf;
}


#undef __b32_b3
#undef __b32_b4
#undef __b32_b5
#undef __b32_b6
#undef __b32_b7
#undef __b32_b8
#undef __b32_b9
#undef __b32_b10
#undef __b32_b11
#undef __b32_b12


/***********************************************************************
 ** Tests
 **********************************************************************/
int64_t 
b32_test (int64_t start, int64_t end, int bits, int mask, 
          char *encoder (char **, const unsigned char *, size_t)) 
{
  unsigned length = (bits+7)/8;
  unsigned char* data = malloc (sizeof (char)*length);
  char *pbuf = 0, *pout = 0;
  int64_t m = 0xff, i, k;
  unsigned j, l, pct = 0;

  for (i = start; i < end; ++i)
    {
      m = 0xff;
      for (j = 0, l = 0; j < length; ++j) 
        {
          data[j] = ((i & m) >> l);
          m <<= 8;
          l += 8;
        }
      data[j-1] &= mask;

      encoder (&pbuf, data, length);
      k = b32_rank (pbuf);
      if (k != i) 
        {
          fprintf (stderr, 
                   "Encoded value %s maps to %lld but expecting %lld!\n",
                   pbuf, k, i);
          exit (1);
        }

      b32_unrank (&pout, k, bits);
      if (strcmp (pout, pbuf) != 0) 
        {
          fprintf (stderr, "Unrank fails to map %lld; expecting %s"
                   " but got %s\n", k, pbuf, pout);
          exit (1);
        }

      j = (unsigned)(((double)(i-start)/(end-start))*100.+.5);
      if (j - pct >= 5) 
        {
          printf ("%lld <=> %s...%u%%\n", i, pbuf, j);
          pct = j;
        }
    }
  free (pbuf);
  free (pout);
  free (data);

  return end - start;
}

void b32_test15 ()
{
  char *pbuf = 0;
  unsigned char data[2];
  int i, r;

  for (i = 0; i < 32768; ++i)
    {
      data[1] = (i & 0x7fff) >> 8;
      data[0] = i & 0xff;

      b32_encode15 (&pbuf, data, 2);
      r = b32_rank (pbuf);
      if (r != i)
        {
          printf ("%d <=> %s <=> %d\n", i, pbuf, r);
          fprintf (stderr, "fatal error: bug in b32_encode15!\n");
          exit (1);
        }
    }
  free (pbuf);
}

void b32_test20 ()
{
  unsigned i, r, size = 1<<20;
  unsigned char data[3];
  char *pbuf = 0;

  for (i = 0; i < size; ++i)
    {
      data[2] = (i & 0x3fffff)>>16;
      data[1] = (i & 0x00ffff)>>8;
      data[0] = (i & 0xff);

      b32_encode20 (&pbuf, data, 3);
      r = b32_rank (pbuf);
      if (r != i)
        {
          printf ("%u <=> %s <=> %u\n", i, pbuf, r);
          fprintf (stderr, "fatal error: bug in b32_encode20!\n");
          exit (1);
        }
      printf ("%s %u\n", pbuf, r);
    }
  free (pbuf);
}

void b32_test25 ()
{
  unsigned i, r, size = 1<<25;
  unsigned char data[4];
  char *pbuf = 0;

  for (i = 0; i < size; ++i)
    {
      data[3] = (i & 0x01ffffff)>>24;
      data[2] = (i & 0x00ffffff)>>16;
      data[1] = (i & 0x0000ffff)>>8;
      data[0] = (i & 0xff);

      b32_encode25 (&pbuf, data, 4);
      r = b32_rank (pbuf);
      if (r != i)
        {
          printf ("%u <=> %s <=> %u\n", i, pbuf, r);
          fprintf (stderr, "fatal error: bug in b32_encode25!\n");
          exit (1);
        }
    }
  free (pbuf);
}

void b32_rank_test ()
{
  char *codec = 0;
  int64_t rank[] = {
    4396213991723l,
    3608230658484l,
    3993059728205l,
    3663206239873l,
    3644881046077l,
    2930198488023l
  };
  int64_t r;
  int i;
  
  for (i = 0; i < sizeof (rank)/sizeof(rank[0]); ++i)
    {
      b32_unrank (&codec, rank[i], 45);
      r = b32_rank (codec);
      printf ("%lld <=> %s <=> %lld\n", rank[i], codec, r);
      if (r != rank[i])
        {
          fprintf (stderr, "fatal error: bug in codec!\n");
          exit (1);
        }
    }
  free (codec);
}

void b32_random_test (int64_t size, int bits, int mask,
                      char *encoder (char **, const unsigned char *, size_t))
{
  unsigned char *data;
  char *pbuf = 0, *pout = 0;
  unsigned j, l, pct = 0, length = (bits+7)/8;
  int64_t x, k, i;
  FILE *randfp = fopen ("/dev/urandom", "rb");

  if (randfp == 0)
    {
      fprintf (stderr, "Can't open random device\n");
      exit (1);
    }

  data = malloc (length);
  for (i = 0; i < size; ++i)
    {
      if (length == fread (data, sizeof (char), length, randfp))
        {
          data[length-1] &= mask;

          x = 0;
          for (j = 0, l = 0; j < length; ++j)
            {
              x |= ((int64_t)data[j] & 0xff) << l;
              l += 8;
            }

          encoder (&pbuf, data, length);
          k = b32_rank (pbuf);
          if (k != x) 
            {
              fprintf (stderr, 
                       "Encoded value %s maps to %lld but expecting %lld!\n",
                       pbuf,  k, x);
              exit (1);
            }

          b32_unrank (&pout, k, bits);
          if (strcmp (pout, pbuf) != 0) 
            {
              fprintf (stderr, "Unrank fails to map %lld; expecting %s"
                       " but got %s\n", k, pbuf, pout);
              exit (1);
            }
          
          j = (unsigned)(((double)i/size)*100.+.5);
          if (j - pct >= 5) 
            {
              printf ("%llx <=> %s...%u%%\n", x, pbuf, j);
              pct = j;
            }
        }
    }
  free (pbuf);
  free (pout);
  free (data);

  fclose (randfp);
}

#ifdef B32_TEST
int main (int argc, char *argv[])
{
  /*
            {15, 0x7f},
            {20, 0x0f},
            {25, 0x01},
            {30, 0x3f},
            {35, 0x07},
            {40, 0xff},
            {45, 0x1f},
            {50, 0x03},
            {55, 0x7f},
            {60, 0x0f}
  */
  //b32_test (1l<<44, (1l<<44)+(1<<20), 45, 0x1f, b32_encode45);
  int64_t size = 1l<<30;
  if (argc > 1)
    {
      int s = atoi (argv[1]);
      if (s < 64)
        size = 1l<<s;
      else
        size = s;
    }
  printf ("## random sampling of %lld...\n", size);
  b32_random_test (size, 60, 0x0f, b32_encode60);

  return 0;
}
#endif
