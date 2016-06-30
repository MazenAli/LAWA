#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC 1

#include <cassert>
#include <iostream>
#include <lawa/math/ceil.h>

namespace lawa {

template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::Wavelet(const Basis<T,Orthogonal,Interval,Multi> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
    switch (d) {
        case 1:
            initialticsize = pow2i<T>(-1);
            break;
        case 2:
            initialticsize = pow2i<T>(-3);
            break;

        case 3:
            initialticsize = pow2i<T>(-4);
            break;

        case 4:
            initialticsize = pow2i<T>(-4);
            break;

        default: std::cerr << "Wavelet<T,Orthogonal,Interval,Multi> not yet realized"
                    " for d = " << d << ". Stopping..." << std::endl;
                    exit(1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    k -= 1;
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)basis._numLeftParts) {
        //std::cerr << " Left boundary: k = " << k << std::endl;
        return pow2ih<T>(2*j*deriv+j) * basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-basis._numLeftParts) % basis._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //std::cerr << "  k = " << k << " : type = " << type << ", shift = " << shift << std::endl;
        return pow2ih<T>(2*j*deriv+j) * 
        basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    //FLENS_DEFAULT_INDEXTYPE shift = iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * basis._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,Interval,Multi>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    k -= 1;
    
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE) basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type = (FLENS_DEFAULT_INDEXTYPE)((k-basis._numLeftParts) % basis._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
    }
    
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k - (basis.cardJ(j) -1 - basis._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    //FLENS_DEFAULT_INDEXTYPE shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
    return pow2i<T>(-j) * (basis._rightSupport[type]+shift);
}

template <typename T>
flens::DenseVector<flens::Array<T> >
Wavelet<T,Orthogonal,Interval,Multi>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    k -= 1;    
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSingularSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-basis._numLeftParts) % basis._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        flens::DenseVector<flens::Array<T> > result = basis._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k - (basis.cardJ(j)-1 - basis._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    //FLENS_DEFAULT_INDEXTYPE shift = iceil<T>((k+1. - basis._numLeftParts)/basis._numInnerParts);
    flens::DenseVector<flens::Array<T> > result = basis._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}
    
template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    //todo: Critical! Yields totally wrong results if too small, increases run-time when too large!
    //return pow2i<T>(-(j+4));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
flens::DenseVector<flens::Array<long double> > *
Wavelet<T,Orthogonal,Interval,Multi>::getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k,
                                                    FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
                                                    FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const
{
	// No split necessary, so set default values
	refinement_k_restart = 1;

    k -= 1;
    refinement_j = j + basis._addRefinementLevel;
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)basis._numLeftParts) {
        refinement_k_first = basis._leftOffsets[k];
    	split = basis._leftRefCoeffs[k].length() + 1;
        return &(basis._leftRefCoeffs[k]);
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-basis._numLeftParts) % basis._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = (FLENS_DEFAULT_INDEXTYPE)lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //refinement_k_first = pow2i<FLENS_DEFAULT_INDEXTYPE>(basis._addRefinementLevel)*shift+basis._innerOffsets[type];
        refinement_k_first = basis._shiftFactor*shift+basis._innerOffsets[type];
    	split = basis._innerRefCoeffs[type].length() + 1;
        return &(basis._innerRefCoeffs[type]);
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = (FLENS_DEFAULT_INDEXTYPE)pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    //refinement_k_first = pow2i<FLENS_DEFAULT_INDEXTYPE>(basis._addRefinementLevel)*shift+basis._rightOffsets[type];
    refinement_k_first = basis._shiftFactor*shift+basis._rightOffsets[type];
	split = basis._rightRefCoeffs[type].length() + 1;
    return &(basis._rightRefCoeffs[type]);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Wavelet<T,Orthogonal,Interval,Multi>::getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    return j + basis._addRefinementLevel;
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::getL2Norm(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const
{
    return 1.;
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    long double pow2ij = (long double)((FLENS_DEFAULT_INDEXTYPE) 1 << j);
    k -= 1;
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)basis._numLeftParts) {
        return pow2ij * basis._leftH1SemiNorms[k];
    }

    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-basis._numLeftParts) % basis._numInnerParts);
        return pow2ij * basis._innerH1SemiNorms[type];
    }

    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    return pow2ij * basis._rightH1SemiNorms[type];
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
