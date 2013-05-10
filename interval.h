#ifndef __interval_h__
#define __interval_h__

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/*
 * opaque interval internal state
 */
typedef struct __interval_s interval_t;

extern interval_t *interval_create ();
extern void interval_free (interval_t *tree);
extern uint32_t interval_encode32 (const interval_t *tree, 
                                   double value, double slack);

#ifdef __cplusplus
}
#endif
#endif /* __interval_h__ */
