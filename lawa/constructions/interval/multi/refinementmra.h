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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTMRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTMRA_H 1

//#include <lawa/aux/integer.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>

namespace lawa {

template <typename _T>
class MRA<_T,Orthogonal,Interval,MultiRefinement>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = MultiRefinement;

        typedef BasisFunction<T,Orthogonal,Interval,MultiRefinement> BasisFunctionType;
        typedef BSpline<T,Orthogonal,Interval,MultiRefinement> BSplineType;

        MRA(int d, int j=-1);

        ~MRA();

        // cardinalities of whole, left, inner, right index sets.
        int
        cardI(int j) const;

        int
        cardIL(int j=-1) const;

        int
        cardII(int j) const;

        int
        cardIR(int j=-1) const;

        // ranges of whole left, inner, right index sets.
        flens::Range<int>
        rangeI(int j) const;

        flens::Range<int>
        rangeIL(int j=-1) const;

        flens::Range<int>
        rangeII(int j) const;

        flens::Range<int>
        rangeIR(int j) const;

        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
        void
        enforceBoundaryCondition();

        const int d;
        const int j0;          // minimal used(!) level.

        BSpline<T,Orthogonal,Interval,MultiRefinement> phi;

    private:
        flens::DenseVector<flens::Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
                                       // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

        friend class BSpline<T,Orthogonal,Interval,MultiRefinement>;
        friend class Basis<T,Orthogonal,Interval,Multi>;

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

        // Refinement coefficients for representation in B-splines on a higher level
        flens::DenseVector<flens::Array<long double> > *_leftRefCoeffs,
                                         *_innerRefCoeffs,
                                         *_rightRefCoeffs;
        long *_leftOffsets,
             *_innerOffsets,
             *_rightOffsets;
        int _shiftFactor;

};

} // namespace lawa

#include <lawa/constructions/interval/multi/refinementmra.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTMRA_H
