
#ifndef __sha1_h__
#define __sha1_h__

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/*
 * opaque sha1 internal state
 */
typedef struct _sha1_state sha1_state;

extern sha1_state *sha1_create ();
extern void sha1_reset (sha1_state *state);
extern void sha1_update (sha1_state *state, 
                         const unsigned char *data, size_t len);
extern void sha1_digest (sha1_state *state, unsigned char digest[20]);
extern void sha1_free (sha1_state *state);

#ifdef __cplusplus
}
#endif
#endif /* __sha1_h__ */
