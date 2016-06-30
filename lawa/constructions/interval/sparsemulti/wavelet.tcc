#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Primal,Interval,SparseMulti>::Wavelet(const Basis<T,Primal,Interval,SparseMulti> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
}
    
template <typename T>
Wavelet<T,Primal,Interval,SparseMulti>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Primal,Interval,SparseMulti>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2ih<T>(2*j*deriv+j) * basis._leftScalingFactors(k) *
                   basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
        }
        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % basis._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            //std::cerr << "type = " << type << ", shift = " << shift << " "  << std::endl;
            return pow2ih<T>(2*j*deriv+j) * basis._innerScalingFactors(type) *
                   basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
        }

        // right boundary
        FLENS_DEFAULT_INDEXTYPE shift = basis.cardJ(j)/2;
        return pow2ih<T>(2*j*deriv+j) * basis._rightScalingFactors(0) *
               basis._rightEvaluator[0](pow2i<T>(j)*x-shift,deriv);
    }
    else { // Control may reach end of non-void function
        assert(d==4);
        return (T) 0;
    }

}
    
template <typename T>
Support<T>
Wavelet<T,Primal,Interval,SparseMulti>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSupport[k];
        }

        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % basis._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
        }

        // right boundary
        FLENS_DEFAULT_INDEXTYPE shift = basis.cardJ(j)/2;
        return pow2i<T>(-j) * (basis._rightSupport[0]+shift);
    }
    else { // Control may reach end of non-void function
        assert(d==4);
        return Support<T>(-0.,0.);
    }
}

template <typename T>
flens::DenseVector<flens::Array<T> >
Wavelet<T,Primal,Interval,SparseMulti>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSingularSupport[k];
        }

        k-=(basis._numLeftParts-1);
        // inner part
        if (k!=basis.rangeJ(j).lastIndex()-basis._numLeftParts) {
            FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k % basis._numInnerParts);
            FLENS_DEFAULT_INDEXTYPE shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
            flens::DenseVector<flens::Array<T> > result = basis._innerSingularSupport[type];
            result += shift;
            return pow2i<T>(-j) * result;
        }

        // right boundary
        FLENS_DEFAULT_INDEXTYPE shift = basis.cardJ(j)/2;
        flens::DenseVector<flens::Array<T> > result = basis._rightSingularSupport[0];
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
Wavelet<T,Primal,Interval,SparseMulti>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    return pow2i<T>(-(j+3));
}
    
} // namespace lawa
