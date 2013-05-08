
#ifndef __spectral_h__
#define __spectral_h__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * opaque spectral internal state
 */
typedef struct __spectral_t spectral_t;

extern spectral_t *spectral_create ();
extern const char * spectral_digest (spectral_t *, const char *inchi);
extern const char * spectral_hashkey (const spectral_t *);
extern const char * spectral_error (const spectral_t *);
extern void spectral_free (spectral_t *);

#ifdef __cplusplus
}
#endif
#endif /* __spectral_h__ */
