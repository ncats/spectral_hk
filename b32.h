
#ifndef __b32_h__
#define __b32_h__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern char *b32_unrank (char **output, int64_t rank, size_t size);
extern int64_t b32_rank (const char *value);

#define __FUNC(BITS) extern char *b32_encode##BITS \
    (char **output, const unsigned char *data, size_t len)

__FUNC(15);
__FUNC(20);
__FUNC(25);
__FUNC(30);
__FUNC(35);
__FUNC(40);
__FUNC(45);
__FUNC(50);
__FUNC(55);
__FUNC(60);

#undef __FUNC
#ifdef __cplusplus
}
#endif
#endif /* __b32_h__ */
