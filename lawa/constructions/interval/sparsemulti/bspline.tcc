#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
BSpline<T,Primal,Interval,SparseMulti>::BSpline(const MRA<T,Primal,Interval,SparseMulti> &_mra)
    : mra(_mra), d(_mra.d)
{
}
    
template <typename T>
BSpline<T,Primal,Interval,SparseMulti>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,Interval,SparseMulti>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    assert(0<=k && k<mra.cardI(j));
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2ih<T>(2*(j)*deriv+j) * mra._leftScalingFactors(0) *
                   mra._leftEvaluator[0](pow2i<T>(j)*x, deriv);
        }
        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % mra._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = (FLENS_DEFAULT_INDEXTYPE)std::ceil(T(k) / T(mra._numInnerParts));
            //std::cerr << "type = " << type << ", shift = " << shift << " " << mra._innerScalingFactors(type) << std::endl;
            return pow2ih<T>(2*j*deriv+j) * mra._innerScalingFactors(type) *
                   mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
        }
        // right boundary
        FLENS_DEFAULT_INDEXTYPE shift = k / mra._numInnerParts + 1;
        return pow2ih<T>(2*j*deriv+j) * mra._rightScalingFactors(0) *
               mra._rightEvaluator[0](pow2i<T>(j)*x-shift,deriv);
    }
    else { // Control may reach end of non-void function
        assert(d==4);
        return (T) 0;
    }
}
    
template <typename T>
Support<T>
BSpline<T,Primal,Interval,SparseMulti>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSupport[0];
        }

        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % mra._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = (FLENS_DEFAULT_INDEXTYPE)std::ceil(T(k) / T(mra._numInnerParts));
            return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
        }

        // right boundary
        FLENS_DEFAULT_INDEXTYPE shift = k / mra._numInnerParts + 1;
        return pow2i<T>(-j) * (mra._rightSupport[0]+shift);
    }
    else { // Control may reach end of non-void function
        assert(d==4);
        return Support<T>(-0.,0.);
    }

}

template <typename T>
flens::DenseVector<flens::Array<T> >
BSpline<T,Primal,Interval,SparseMulti>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    if (d==4) {
        // left boundary
        if (k<mra._numLeftParts) {
            return pow2i<T>(-j) * mra._leftSingularSupport[0];
        }

        // inner part
        if (k!=mra.rangeI(j).lastIndex()) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % mra._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = (FLENS_DEFAULT_INDEXTYPE)std::ceil(T(k) / T(mra._numInnerParts));
            flens::DenseVector<flens::Array<T> > result = mra._innerSingularSupport[type];
            result += shift;
            return pow2i<T>(-j) * result;
        }

        // right part
        FLENS_DEFAULT_INDEXTYPE shift = k / mra._numInnerParts + 1;
        flens::DenseVector<flens::Array<T> > result = mra._rightSingularSupport[0];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    else { // Control may reach end of non-void function
        assert(d==4);
        return flens::DenseVector<flens::Array<T> >();
    }

}

template <typename T>
T
BSpline<T,Primal,Interval,SparseMulti>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    return pow2i<T>(-(j+3));
}
        
} // namespace lawa
