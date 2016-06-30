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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Primal,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef Basis<T,Primal,Interval,Dijkema>         RefinementBasis;
        typedef BasisFunction<T,Primal,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Dijkema> BSplineType;
        typedef Wavelet<T,Primal,Interval,Dijkema> WaveletType;

        Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j=-1);
        
        ~Basis();

        FLENS_DEFAULT_INDEXTYPE
        level() const;

        void
        setLevel(FLENS_DEFAULT_INDEXTYPE j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const BasisFunctionType &
        generator(XType xtype) const;

        // cardinalities of whole, left, inner, right index sets (primal).
        FLENS_DEFAULT_INDEXTYPE
        cardJ(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJL(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJI(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJR(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJL(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJI(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJR(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        /// Returns the range of indicicated functions wavelets from SecondBasis or
        /// SecondRefinemementBasiswhose supports intersect the support
        /// of a refinement B-spline with level j_bspline and translation index k_bspline
        /// from the current RefinementBasis. This is required for tree-based algorithms.
        /// The returned level j_wavelet or j_bspline2 is chosen s.t. there is "no scale difference"
        /// i.e., if we refine both functions the corresponding refinements should
        /// live on the same scale.
        template <typename SecondBasis>
        void
        getWaveletNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline, FLENS_DEFAULT_INDEXTYPE k_bspline,
                                      const SecondBasis &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        void
        getWaveletNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline, FLENS_DEFAULT_INDEXTYPE k_bspline,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        template <typename SecondRefinementBasis>
        void
        getBSplineNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline1, FLENS_DEFAULT_INDEXTYPE k_bspline1,
                                      const SecondRefinementBasis &secondrefinementbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_bspline2,
                                      FLENS_DEFAULT_INDEXTYPE &k_bspline2_first, FLENS_DEFAULT_INDEXTYPE &k_bspline2_last) const;


        /// Returns the range of indicated functions from SecondBasis whose supports
        /// intersect the support of a given (multi-)scaling with level j_scaling and translation index
        /// k_scaling from the current Basis. This is required for tree-based algorithms.
        /// The returned level is chosen s.t. the corresponding refinements live
        /// on the same scale.
                                  
        void
        getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                      const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;
                                      
        void
        getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;

        /*template <typename SecondBasis>
        void
        getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                      const SecondBasis &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
        */
        
        void
        getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                      const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
                                      
        void
        getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
                                      
        /// Returns the range of indicated functions from SecondBasis and SecondRefinementBasis
        /// whose supports intersect the support of a given wavelet with level j_wavelet and
        /// translation index k_wavelet from the current Basis. This is required for tree-based algorithms.
        /// The returned level of the functions is chosen s.t. there is "no scale difference", i.e.,
        /// the corresponding refinements should live on the same scale.
        template <typename SecondBasis>
            void
            getBSplineNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                                          const SecondBasis &secondrefinementbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_bspline, FLENS_DEFAULT_INDEXTYPE &k_bspline_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_bspline_last) const;

        template <typename SecondBasis>
        void
        getScalingNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                                      const SecondBasis &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_scaling, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;
                           
        void
        getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                      const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
        
        void
        getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        template <typename SecondBasis>
        void
        getLowerWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                           const SecondBasis &secondbasis,
                                           FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                           FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        /*template <typename SecondBasis>
        void
        getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                           const SecondBasis &secondbasis,
                                           FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                           FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
        */
                                       
       void
       getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                          const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;
                                                                          
   		void
   		getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
   										   const Basis<T,Primal,Periodic,CDF> &secondbasis,
   										   FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
   										   FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        MRA<T,Primal,Interval,Dijkema> mra;
        MRA<T,Dual,Interval,Dijkema>  mra_;

        const FLENS_DEFAULT_INDEXTYPE d, d_, mu;   // mu = mu(d) = d&1.
        const FLENS_DEFAULT_INDEXTYPE min_j0;      // minimal allowed(!) level;
        const FLENS_DEFAULT_INDEXTYPE j0;          // minimal used(!) level.

    private:
        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable FLENS_DEFAULT_INDEXTYPE _j;                // the current level.

        friend class Wavelet<T,Primal,Interval,Dijkema>;

        flens::DenseVector<flens::Array<long double> > *_leftRefCoeffs,
                                         *_innerRefCoeffs,
                                         *_rightRefCoeffs;

        long double *_leftL2Norms,  *_leftH1SemiNorms,
                    *_innerL2Norms, *_innerH1SemiNorms,
                    *_rightL2Norms, *_rightH1SemiNorms;
        FLENS_DEFAULT_INDEXTYPE *_leftOffsets,
             *_innerOffsets,
             *_rightOffsets;

        Basis(const Basis& secondbasis);


    public:
        Wavelet<T,Primal,Interval,Dijkema> psi;
        Basis<T,Primal,Interval,Dijkema> &refinementbasis;

        flens::RefinementMatrix<T,Interval,Dijkema> M1;

        class LaplaceOperator1D {
            public:
                LaplaceOperator1D(FLENS_DEFAULT_INDEXTYPE _d,
                                  const Basis<_T,Primal,Interval,Dijkema> &_refinementbasis);

                T
                operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2);

            private:
                FLENS_DEFAULT_INDEXTYPE d;
                const Basis<_T,Primal,Interval,Dijkema> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > outer_values;
                flens::DenseVector<flens::Array<long double> > inner_values;
        };

        class IdentityOperator1D {
            public:
                IdentityOperator1D(FLENS_DEFAULT_INDEXTYPE _d,
                                   const Basis<_T,Primal,Interval,Dijkema> &_refinementbasis);

                T
                operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2);

            private:
                FLENS_DEFAULT_INDEXTYPE d;
                const Basis<_T,Primal,Interval,Dijkema> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > outer_values;
                flens::DenseVector<flens::Array<long double> > inner_values;
        };

        LaplaceOperator1D LaplaceOp1D;
        IdentityOperator1D IdentityOp1D;
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H

