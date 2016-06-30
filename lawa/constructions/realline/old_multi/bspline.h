#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
        
        BSpline(const FLENS_DEFAULT_INDEXTYPE _d);

        //TODO    BSpline(MRA<T,Orthogonal,R,Multi> &mra); 
        
        virtual
        ~BSpline();
        
        T
        operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const;
        
        Support<T>
        support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;
        
        Support<T>
        max_support() const;

        flens::DenseVector<flens::Array<T> >
        singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        T
        tic(FLENS_DEFAULT_INDEXTYPE j) const;

        const unsigned FLENS_DEFAULT_INDEXTYPE d;
        unsigned FLENS_DEFAULT_INDEXTYPE _numSplines;

    private:
        typedef T (*Evaluator)(T x, unsigned short deriv);

        FLENS_DEFAULT_INDEXTYPE
        _shift(FLENS_DEFAULT_INDEXTYPE k) const;

        FLENS_DEFAULT_INDEXTYPE
        _type(FLENS_DEFAULT_INDEXTYPE k) const;

        Evaluator *_evaluator;
        Support<T> *_support;
        flens::DenseVector<flens::Array<T> > *_singularSupport;

        Support<T> _max_support;
    //        T
    //TODO    tic(FLENS_DEFAULT_INDEXTYPE j) const;
    //    FLENS_DEFAULT_INDEXTYPE polynomialOrder;
};

} // namespace lawa

#include <lawa/constructions/realline/multi/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
