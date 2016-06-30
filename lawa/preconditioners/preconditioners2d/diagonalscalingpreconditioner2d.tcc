#include <cassert>

namespace lawa {

template<typename T>
DiagonalScalingPreconditioner2D<T>::DiagonalScalingPreconditioner2D(FLENS_DEFAULT_INDEXTYPE sx, FLENS_DEFAULT_INDEXTYPE sy)
    : _sx(sx), _sy(sy)
{
    assert(sx >= 0);
    assert(sy >= 0);
}

template<typename T>
T
DiagonalScalingPreconditioner2D<T>::operator()(XType XisSpline, FLENS_DEFAULT_INDEXTYPE jx, FLENS_DEFAULT_INDEXTYPE /*kx*/,
                                               XType YisSpline, FLENS_DEFAULT_INDEXTYPE jy, FLENS_DEFAULT_INDEXTYPE /*ky*/) const
{
    assert(XisSpline);
    assert(YisSpline);

    return pow2i<T>(- _sx*jx - _sy*jy);
}

} // namespace lawa

