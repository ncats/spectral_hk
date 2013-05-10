
#ifndef __sha1_h__
#define __sha1_h__

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/*
 * opaque sha1 internal state
 */
typedef struct __sha1_s sha1_t;

extern sha1_t *sha1_create ();
extern void sha1_reset (sha1_t *state);
extern void sha1_update (sha1_t *state, 
                         const unsigned char *data, size_t len);
extern void sha1_digest (sha1_t *state, unsigned char digest[20]);
extern void sha1_free (sha1_t *state);

#ifdef __cplusplus
}
#endif
#endif /* __sha1_h__ */
