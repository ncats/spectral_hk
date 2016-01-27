
#ifndef __periodic_h__
#define __periodic_h__

#ifdef __cplusplus
extern "C" {
#endif
  
typedef struct __element_s {
  const char *symbol;
  int atno;
  double mass;
  double eneg; /* electronegativity */
  double eaff; /* electron affinity */
  double ione; /* ionization energy */
} element_t;

extern const element_t *element_lookup_symbol (const char *symbol);
extern const element_t *element_lookup_atno (int atno);

#ifdef __cplusplus
}
#endif
#endif /* !__periodic_h__ */
