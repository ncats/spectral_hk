#ifndef __inchi_h__
#define __inchi_h__

#ifdef __cplusplus
extern "C" {
#endif

/* opaque inchi */
typedef struct __inchi_s inchi_t;

  extern inchi_t *inchi_create ();
  extern void inchi_free (inchi_t *);
  extern int inchi_parse (inchi_t *, const char *inchi);
  extern int inchi_node_count (const inchi_t *);
  extern int inchi_edge_count (const inchi_t *);
  extern const char *inchi_error (const inchi_t *);
  extern const char *inchi_layer_c (const inchi_t *);

  extern size_t inchi_matrix_size (const inchi_t *);  
  extern const int *inchi_matrix_A (const inchi_t *);
  
#ifdef __cplusplus
}
#endif
#endif /* !__inchi_h__ */
