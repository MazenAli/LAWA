#ifndef LAWA_AUX_COMPILETIME_ASSERT_H
#define LAWA_AUX_COMPILETIME_ASSERT_H 1

#define ct_assert(x) FLENS_DEFAULT_INDEXTYPE __compiletime_assert[(FLENS_DEFAULT_INDEXTYPE)x] __attribute__ ((unused))

#endif // LAWA_AUX_COMPILETIME_ASSERT

