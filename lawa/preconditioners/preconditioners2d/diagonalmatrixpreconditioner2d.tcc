#include <cmath>

namespace lawa {

template <typename T,  typename BilinearForm>
DiagonalMatrixPreconditioner2D<T, BilinearForm>::DiagonalMatrixPreconditioner2D(BilinearForm &a)
    : _a(a), theta(0.), timestep(0.)
{
}

template <typename T, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,BilinearForm>::operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1,
                                                                   XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2) const
{
    T val = fabs(_a(xtype1,j1,k1,xtype2,j2,k2,xtype1,j1,k1,xtype2,j2,k2));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,BilinearForm>::operator()(const Index2D &index) const
{
    T val = fabs(_a(index,index));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,BilinearForm>::operator[](const Index2D &index)
{
    T val = fabs(_a(index,index));
    if (theta==0) {
        return 1./sqrt(val);
    }
    else {
        return 1./sqrt(1 + theta*timestep*val);
    }
}

template <typename T, typename BilinearForm>
void
DiagonalMatrixPreconditioner2D<T,BilinearForm>
::setThetaTimeStepParameters(T _theta, T _timestep)
{
    theta = _theta;
    timestep = _timestep;
}

/*
template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator[](const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}*/

}   // namespace lawa

