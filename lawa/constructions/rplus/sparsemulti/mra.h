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
 
#ifndef LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H
#define LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H 1

#include <cassert>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.h>

namespace lawa {
    
template <typename _T>
class MRA<_T,Primal,RPlus,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = RPlus;
        static const Construction Cons = SparseMulti;
        
        typedef BasisFunction<T,Primal,RPlus,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,RPlus,SparseMulti>       BSplineType;
        
        MRA(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE j=1);
        
        ~MRA();
        
        // cardinalities of left index sets.
                FLENS_DEFAULT_INDEXTYPE
        cardIL(FLENS_DEFAULT_INDEXTYPE j=1) const;
        
        // ranges of left index sets.
        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeIL(FLENS_DEFAULT_INDEXTYPE j=-1) const;
        
        FLENS_DEFAULT_INDEXTYPE
        level() const;
        
        void
        setLevel(FLENS_DEFAULT_INDEXTYPE j) const;
        
        template <BoundaryCondition BC>
        void
        enforceBoundaryCondition();
        
        const FLENS_DEFAULT_INDEXTYPE d;     
        const FLENS_DEFAULT_INDEXTYPE j0;          // minimal used(!) level.
        
        BSpline<T,Primal,RPlus,SparseMulti> phi;
        
    private:
        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
        
        mutable FLENS_DEFAULT_INDEXTYPE _j;                // the current level.
    
        friend class BSpline<T,Primal,RPlus,SparseMulti>;
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        unsigned FLENS_DEFAULT_INDEXTYPE _numLeftParts,
                     _numInnerParts;
        Evaluator *_leftEvaluator,
                  *_innerEvaluator;
        Support<T> *_leftSupport,
                   *_innerSupport;
        flens::DenseVector<flens::Array<T> > *_leftSingularSupport,
                               *_innerSingularSupport;
        flens::DenseVector<flens::Array<T> > _leftScalingFactors,
                               _innerScalingFactors;
};
    
} // namespace lawa

#include <lawa/constructions/rplus/sparsemulti/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_MRA_H
