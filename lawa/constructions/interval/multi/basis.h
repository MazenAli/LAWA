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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/interval/multi/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;
        
        typedef Basis<T,Orthogonal,Interval,MultiRefinement> RefinementBasis;
        typedef BasisFunction<T,Orthogonal,Interval,Multi>   BasisFunctionType;
        typedef BSpline<T,Orthogonal,Interval,Multi>         BSplineType;
        typedef Wavelet<T,Orthogonal,Interval,Multi>         WaveletType;

        Basis(const FLENS_DEFAULT_INDEXTYPE d, const FLENS_DEFAULT_INDEXTYPE j=-1);
    
        virtual
        ~Basis();
    
        FLENS_DEFAULT_INDEXTYPE
        level() const;
    
        void
        setLevel(const FLENS_DEFAULT_INDEXTYPE j) const;
    
        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();
    
        const BasisFunctionType &
        generator(XType xtype) const;

        //--- cardinalities of whole, left, inner, right index set.
        FLENS_DEFAULT_INDEXTYPE
        cardJ(const FLENS_DEFAULT_INDEXTYPE j) const;
        
        FLENS_DEFAULT_INDEXTYPE
        cardJL(const FLENS_DEFAULT_INDEXTYPE j=-1) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJI(const FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJR(const FLENS_DEFAULT_INDEXTYPE j=-1) const;
    
        //--- ranges of whole, left, inner, right index set.
        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ(const FLENS_DEFAULT_INDEXTYPE j) const;
    
        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJL(const FLENS_DEFAULT_INDEXTYPE j=-1) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJI(const FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJR(const FLENS_DEFAULT_INDEXTYPE j=-1) const;

        void
        getLowerEnclosingScaling(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                                 FLENS_DEFAULT_INDEXTYPE &j_scaling, FLENS_DEFAULT_INDEXTYPE &k_scaling) const;

        void
        getLowerEnclosingWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                 FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet2) const;

        /// Returns the range of indicated functions from SecondBasis whose supports
        /// intersect the support of a given (multi-)scaling with level j_scaling and translation index
        /// k_scaling from the current Basis. This is required for tree-based algorithms.
        /// The returned level of the scaling is chosen s.t. the corresponding refinements live
        /// on the same scale.
        template <typename SecondBasis>
            void
            getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                          const SecondBasis &secondbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;

        template <typename SecondBasis>
            void
            getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling, FLENS_DEFAULT_INDEXTYPE k_scaling,
                                          const SecondBasis &secondbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;


        /// Returns the range of indicated functions from SecondBasis and SecondRefinementBasis
        /// whose supports intersect the support of a given wavelet with level j_wavelet and
        /// translation index k_wavelet from the current Basis. This is required for tree-based algorithms.
        /// The returned level of the functions is chosen s.t. there is "no scale difference", i.e.,
        /// the corresponding refinements should live on the same scale.
        template <typename SecondRefinementBasis>
            void
            getBSplineNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                                          const SecondRefinementBasis &secondrefinementbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_bspline, FLENS_DEFAULT_INDEXTYPE &k_bspline_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_bspline_last) const;

        template <typename SecondBasis>
            void
            getScalingNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                                          const SecondBasis &secondbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_scaling, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;

        template <typename SecondBasis>
            void
            getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                          const SecondBasis &secondbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        template <typename SecondBasis>
            void
            getLowerWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                               const SecondBasis &secondbasis,
                                               FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                               FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        template <typename SecondBasis>
            void
            getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                               const SecondBasis &secondbasis,
                                               FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                               FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        MRA<T,Orthogonal,Interval,Multi> mra;

        const FLENS_DEFAULT_INDEXTYPE d;
        const FLENS_DEFAULT_INDEXTYPE d_;
        const FLENS_DEFAULT_INDEXTYPE j0;          // minimal used(!) level.

        unsigned FLENS_DEFAULT_INDEXTYPE _numLeftParts,
                     _numInnerParts,
                     _numRightParts;

        FLENS_DEFAULT_INDEXTYPE _addRefinementLevel;    //B-splines for refinement are needed on higher levels
        FLENS_DEFAULT_INDEXTYPE _shiftFactor;           //Needed since we have multiple B-spline generators for refinement.




   // private:      // should be private one fine day
        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
                                       // bc(1) = 1 -> Dirichlet BC right.
        
        mutable FLENS_DEFAULT_INDEXTYPE _j;                // the current level.

        Wavelet<T,Orthogonal,Interval,Multi> psi;

        Basis<T,Orthogonal,Interval,MultiRefinement> refinementbasis;

        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        friend class Wavelet<T,Orthogonal,Interval,Multi>;


        Evaluator *_leftEvaluator, 
                  *_innerEvaluator, 
                  *_rightEvaluator;
        Support<T> *_leftSupport, 
                   *_innerSupport, 
                   *_rightSupport;
        flens::DenseVector<flens::Array<T> > *_leftSingularSupport, 
                               *_innerSingularSupport, 
                               *_rightSingularSupport;
        
        flens::DenseVector<flens::Array<long double> > *_leftRefCoeffs,
                                         *_innerRefCoeffs,
                                         *_rightRefCoeffs;

        long double *_leftH1SemiNorms, *_innerH1SemiNorms, *_rightH1SemiNorms;

        FLENS_DEFAULT_INDEXTYPE *_leftOffsets,
             *_innerOffsets,
             *_rightOffsets;


};

} // namespace lawa

#include <lawa/constructions/interval/multi/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H
