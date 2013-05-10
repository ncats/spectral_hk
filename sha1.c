/* -*- mode: c; -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "sha1.h"

struct __sha1_s 
{
  size_t size; /* bytes processed */
  size_t count; /* overflow blocks */
  uint32_t H[5]; /* hash state */
  unsigned char block[64]; /* message block */
  /* work area */
  uint32_t W[80]; 
};

/*
 * big endian macros
 */
#define __get_uint32(a,i) \
  (((a)[(i)]<<24)|((a)[(i)+1]<<16)|((a)[(i)+2]<<8)|((a)[(i)+3]))

#define __set_uint32(a,i,x)                   \
  do { \
    (a)[(i)] = ((x)>>24) & 0xff; \
    (a)[(i)+1] = ((x)>>16) & 0xff; \
    (a)[(i)+2] = ((x)>>8) & 0xff; \
    (a)[(i)+3] = (x) & 0xff; \
  } while (0)

#define ROTL32(x,n) (((x) << (n)) | ((x) >> (32 - (n))))

/*
 * round functions
 */
#define Ch(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define Parity(x,y,z) ((x) ^ (y) ^ (z))
#define Maj(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))

/*
 * round constants
 */
#define K0 0x5a827999
#define K1 0x6ed9eba1
#define K2 0x8f1bbcdc
#define K3 0xca62c1d6

static const unsigned char sha1pad[64] =
{
 0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


sha1_t *
sha1_create ()
{
  sha1_t *state = malloc (sizeof (struct __sha1_s));
  if (state != 0)
    sha1_reset (state);
  return state;
}

void
sha1_free (sha1_t *state)
{
  if (state != 0)
    free (state);
}

void
sha1_reset (sha1_t *state)
{
  if (state != 0)
    {
      state->size = 0;
      state->count = 0;
      state->H[0] = 0x67452301;
      state->H[1] = 0xefcdab89;
      state->H[2] = 0x98badcfe;
      state->H[3] = 0x10325476;
      state->H[4] = 0xc3d2e1f0;
    }
}

static void
process_block (sha1_t *state, const unsigned char M[64])
{
  int t, i;
  uint32_t T, a, b, c, d, e;
  uint32_t *W = state->W;

  for (t = 0, i = 0; t < 16; ++t, i += 4)
    W[t] = __get_uint32(M, i);

  for (; t < 80; ++t)
    W[t] = ROTL32(W[t-3] ^ W[t-8] ^ W[t-14] ^ W[t-16], 1);

  a = state->H[0];
  b = state->H[1];
  c = state->H[2];
  d = state->H[3];
  e = state->H[4];

  /* [0,20) */
  for (t = 0; t < 20; ++t)
    {
      T = ROTL32(a, 5) + Ch(b,c,d) + e + K0 + W[t];
      e = d;
      d = c;
      c = ROTL32(b, 30);
      b = a;
      a = T;
      /*printf ("t=%2d: %08X %08X %08X %08X %08X\n", t, a, b, c, d, e);*/
    }

  /* [20,40) */
  for (; t < 40; ++t)
    {
      T = ROTL32(a, 5) + Parity(b,c,d) + e + K1 + W[t];
      e = d;
      d = c;
      c = ROTL32(b, 30);
      b = a;
      a = T;      
      /*printf ("t=%2d: %08X %08X %08X %08X %08X\n", t, a, b, c, d, e);*/
    }
  
  /* [40,60) */
  for (; t < 60; ++t)
    {
      T = ROTL32(a, 5) + Maj(b,c,d) + e + K2 + W[t];
      e = d;
      d = c;
      c = ROTL32(b, 30);
      b = a;
      a = T;      
      /*printf ("t=%2d: %08X %08X %08X %08X %08X\n", t, a, b, c, d, e);*/
    }
  
  /* [60,80) */
  for (; t < 80; ++t)
    {
      T = ROTL32(a, 5) + Parity(b,c,d) + e + K3 + W[t];
      e = d;
      d = c;
      c = ROTL32(b, 30);
      b = a;
      a = T;
      /*printf ("t=%2d: %08X %08X %08X %08X %08X\n", t, a, b, c, d, e);*/
    }

  state->H[0] += a;
  state->H[1] += b;
  state->H[2] += c;
  state->H[3] += d;
  state->H[4] += e;
}

void
sha1_update (sha1_t *state, const unsigned char *input, size_t len)
{
  size_t leftover, fill;

  if (len == 0) return;
  
  leftover = state->size & 0x3f;
  fill = 64 - leftover;
  state->size += len;

  /* overflow */
  if (state->size < len)
    ++state->count;
  
  if (leftover > 0 && len >= fill)
    {
      memcpy (state->block + leftover, input, fill);
      process_block (state, state->block);
      input += fill;
      len -= fill;
      leftover = 0;
    }

  while (len >= 64)
    {
      process_block (state, input);
      input += 64;
      len -= 64;
    }

  if (len > 0)
    memcpy (state->block + leftover, input, len);
}

void
sha1_digest (sha1_t *state, unsigned char digest[20])
{
    uint32_t last, padn;
    uint32_t high, low;
    unsigned char pad[8];

    high = (state->size >> 29) | (state->count << 3);
    low  = state->size * 8;

    __set_uint32(pad, 0, high);
    __set_uint32(pad, 4, low);

    last = state->size & 0x3F;
    padn = (last < 56) ? (56 - last) : (120 - last);

    sha1_update (state, sha1pad, padn);
    sha1_update (state, pad, 8);

    __set_uint32(digest, 0, state->H[0]);
    __set_uint32(digest, 4, state->H[1]);
    __set_uint32(digest, 8, state->H[2]);
    __set_uint32(digest, 12, state->H[3]);
    __set_uint32(digest, 16, state->H[4]);
}

