#include <lawa/math/factorial.h>

namespace lawa {

FLENS_DEFAULT_INDEXTYPE
factorial(FLENS_DEFAULT_INDEXTYPE n)
{
    FLENS_DEFAULT_INDEXTYPE fac = 1;
    for (FLENS_DEFAULT_INDEXTYPE i=2; i<=n; ++i) {
        fac *= i;
    }
    return fac;
}

} // namespace lawa

