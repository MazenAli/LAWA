/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
 
#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>

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

        flens::DenseVector<flens::Array<long double> > *
        getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first) const;

        FLENS_DEFAULT_INDEXTYPE
        getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const;

        //T
        //getL2Norm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        //T
        //getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        const unsigned FLENS_DEFAULT_INDEXTYPE d;
        unsigned FLENS_DEFAULT_INDEXTYPE _numSplines;

//    private:      // should be private one fine day

        typedef T (*Evaluator)(T x, unsigned short deriv);

        FLENS_DEFAULT_INDEXTYPE
        _shift(FLENS_DEFAULT_INDEXTYPE k) const;

        FLENS_DEFAULT_INDEXTYPE
        _type(FLENS_DEFAULT_INDEXTYPE k) const;

        Evaluator                        *_evaluator;
        Support<T>                       *_support;
        flens::DenseVector<flens::Array<T> >           *_singularSupport;
        Support<T>                       _max_support;

        flens::DenseVector<flens::Array<long double> > *_refCoeffs;
        FLENS_DEFAULT_INDEXTYPE                             *_offsets;
        long double                      *_H1SemiNorms;
        T                                _initialticsize;
        FLENS_DEFAULT_INDEXTYPE                              _addRefinementLevel;    //B-splines for refinement are needed on higher levels
        FLENS_DEFAULT_INDEXTYPE                              _shiftFactor;           //Needed since we have multiple B-spline generators for refinement.
};

} // namespace lawa

#include <lawa/constructions/realline/multi/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_BSPLINE_H
