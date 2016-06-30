#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBSPLINE_TCC 1

#include <cassert>
#include <iostream>
#include <lawa/math/lawa_math.h>

namespace lawa {

template <typename T>
BSpline<T,Orthogonal,Interval,MultiRefinement>
::BSpline(const MRA<T,Orthogonal,Interval,MultiRefinement> &_mra)
    : mra(_mra), d(_mra.d), initialticsize(pow2i<T>(-3))
{
    switch (d) {
        case 1:
            initialticsize = 1.;
            break;
        case 2:
            initialticsize = 1.;
            break;

        case 3:
            initialticsize = pow2i<T>(-3);
            break;

        case 4:
            initialticsize = pow2i<T>(-2);
            break;

        default: std::cerr << "BSpline<T,Orthogonal,Interval,MultiRefinement> not yet realized"
                    " for d = " << d << ". Stopping." << std::endl;
                    exit(-1);
    }

}

template <typename T>
BSpline<T,Orthogonal,Interval,MultiRefinement>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,MultiRefinement>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE) mra._numLeftParts) {
        return pow2ih<T>(2*j*deriv+j) * mra._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }

    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-mra._numLeftParts) % mra._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = (k-mra._numLeftParts)/mra._numInnerParts;
        return pow2ih<T>(2*j*deriv+j) *
               mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }

    // right part
    assert(k<mra.cardIL()+mra.cardII(j)+mra.cardIR());
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-2;
    return pow2ih<T>(2*j*deriv+j) * mra._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}

template <typename T>
Support<T>
BSpline<T,Orthogonal,Interval,MultiRefinement>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSupport[k];
    }

    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-mra._numLeftParts) % mra._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = (k-mra._numLeftParts)/mra._numInnerParts;
        return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
    }

    // right part
    assert(k<mra.cardIL()+mra.cardII(j)+mra.cardIR());
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-2;
    return pow2i<T>(-j) * (mra._rightSupport[type]+shift);
}

template <typename T>
flens::DenseVector<flens::Array<T> >
BSpline<T,Orthogonal,Interval,MultiRefinement>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSingularSupport[k];
    }

    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-mra._numLeftParts) % mra._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = (k-mra._numLeftParts)/mra._numInnerParts;
        flens::DenseVector<flens::Array<T> > result = mra._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }

    // right part
    assert(k<mra.cardIL()+mra.cardII(j)+mra.cardIR());
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-2;
    flens::DenseVector<flens::Array<T> > result = mra._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,MultiRefinement>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    //return pow2i<T>(-(j+3));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
flens::DenseVector<flens::Array<long double> > *
BSpline<T,Orthogonal,Interval,MultiRefinement>::
getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
		FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const
{
	// No split necessary, so set default values
	refinement_k_restart = 1;

    refinement_j = j + 1;
    // left boundary
    if (k<(FLENS_DEFAULT_INDEXTYPE)mra._numLeftParts) {
        refinement_k_first = mra._leftOffsets[k];
        split = mra._leftRefCoeffs[k].length()+1;
        return &(mra._leftRefCoeffs[k]);
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)((k-mra._numLeftParts) % mra._numInnerParts);
        FLENS_DEFAULT_INDEXTYPE shift = (k-mra._numLeftParts)/mra._numInnerParts;
        refinement_k_first = mra._shiftFactor*shift+mra._innerOffsets[type];
        split = mra._innerRefCoeffs[type].length()+1;
        return &(mra._innerRefCoeffs[type]);
    }
    // right part
    assert(k<mra.cardIL()+mra.cardII(j)+mra.cardIR());
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-2;
    refinement_k_first = mra._shiftFactor*shift+mra._rightOffsets[type];
    split = mra._rightRefCoeffs[type].length()+1;
    return &(mra._rightRefCoeffs[type]);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
BSpline<T,Orthogonal,Interval,MultiRefinement>::getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    return j + 1;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_REFINEMENTBSPLINE_TCC
