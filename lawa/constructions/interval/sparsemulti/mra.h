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
 
#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_MRA_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>

namespace lawa {
    
template <typename _T>
class MRA<_T,Primal,Interval,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = SparseMulti;
        
        typedef BasisFunction<T,Primal,Interval,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,SparseMulti>       BSplineType;
        
        MRA(int d, int j=1);
        
        ~MRA();
        
        Support<T>
        max_support() const;

        //For uniform schemes, we "int" as index type is sufficient
        int
        cardI(int j) const;

        int
        cardIL(int j=1) const;

        int
        cardII(int j) const;

        int
        cardIR(int j=-1) const;

        flens::Range<int>
        rangeI(int j) const;

        flens::Range<int>
        rangeIL(int j=-1) const;

        flens::Range<int>
        rangeII(int j) const;

        flens::Range<int>
        rangeIR(int j) const;
        
        //For adaptive schemes, we may require "long" as index type for local scaling function repr.
        long
        long_cardI(int j) const;

        long
        long_cardIL(int j=1) const;

        long
        long_cardII(int j) const;

        long
        long_cardIR(int j=-1) const;

        flens::Range<long>
        long_rangeI(int j) const;

        flens::Range<long>
        long_rangeIL(int j=-1) const;

        flens::Range<long>
        long_rangeII(int j) const;

        flens::Range<long>
        long_rangeIR(int j) const;


        int
        level() const;
        
        void
        setLevel(int j) const;
        
        template <BoundaryCondition BC>
        void
        enforceBoundaryCondition();
        
        const int d;     
        const int j0;          // minimal used(!) level.
        
        BSpline<T,Primal,Interval,SparseMulti> phi;
        unsigned int _numSplines;
        
    private:
        flens::DenseVector<flens::Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
                                       // bc(1) = 1 -> Dirichlet BC right.
        
        mutable int _j;                // the current level.
    
        friend class BSpline<T,Primal,Interval,SparseMulti>;
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        unsigned int _numLeftParts,
                     _numInnerParts,
                     _numRightParts;
        Evaluator *_leftEvaluator,
                  *_innerEvaluator,
                  *_rightEvaluator;
        Support<T> *_leftSupport,
                   *_innerSupport,
                   *_rightSupport;
        flens::DenseVector<flens::Array<T> > *_leftSingularSupport,
                               *_innerSingularSupport,
                               *_rightSingularSupport;
        flens::DenseVector<flens::Array<T> > _leftScalingFactors,
                               _innerScalingFactors,
                               _rightScalingFactors;
};
    
} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_MRA_H
