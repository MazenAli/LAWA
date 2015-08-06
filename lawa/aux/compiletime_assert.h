#ifndef LAWA_AUX_COMPILETIME_ASSERT_H
#define LAWA_AUX_COMPILETIME_ASSERT_H 1

#define ct_assert(x) int __compiletime_assert[(int)x] __attribute__ ((unused))

#endif // LAWA_AUX_COMPILETIME_ASSERT

