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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H 1

#include <lawa/settings/enum.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/refinementmatrix.h>
#include <lawa/constructions/periodic/primal/mra.h>
#include <lawa/constructions/periodic/dual/mra.h>
#include <lawa/constructions/periodic/primal/wavelet.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Primal,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef Basis<T,Primal,Interval,Dijkema>        RefinementBasis;
        typedef BasisFunction<T,Primal,Periodic,CDF> 	BasisFunctionType;
        typedef BSpline<T,Primal,Periodic,CDF> 			BSplineType;
        typedef Wavelet<T,Primal,Periodic,CDF> 			WaveletType;

        Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j=0);

        ~Basis();

        FLENS_DEFAULT_INDEXTYPE
        level() const;

        void
        setLevel(FLENS_DEFAULT_INDEXTYPE j) const;

        const BasisFunctionType &
        generator(XType xtype) const;
        
        FLENS_DEFAULT_INDEXTYPE
        cardJ(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJL(FLENS_DEFAULT_INDEXTYPE j = 0) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJI(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJR(FLENS_DEFAULT_INDEXTYPE j = 0) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJL(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJI(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJR(FLENS_DEFAULT_INDEXTYPE j) const;

        /// Returns the range of indicated functions from SecondBasis whose supports
        /// intersect the support of a given (multi-)scaling with level j_scaling and translation index
        /// k_scaling from the current Basis. This is required for tree-based algorithms.
        /// The returned level is chosen s.t. the corresponding refinements live
        /// on the same scale.
        void
        getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
        							  const Basis<T,Primal,Periodic,CDF> &secondbasis,
        							  FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
        							  FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;

        void
        getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
        							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
        							  FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
        							  FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const;

        void
        getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
        							  const Basis<T,Primal,Periodic,CDF> &secondbasis,
        							  FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
        							  FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        void
        getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
        							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
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

        void
        getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        void
        getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                      const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                      FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                      FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;


        template <typename SecondBasis>
		void
		getLowerWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
										   const SecondBasis &secondbasis,
										   FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
										   FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

		void
		getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
										   const Basis<T,Primal,Periodic,CDF> &secondbasis,
										   FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
										   FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        void
        getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                           const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                           FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                           FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;

        const FLENS_DEFAULT_INDEXTYPE d, d_, j0;
        MRA<T,Primal,Periodic,CDF> mra;
        MRA<T,Dual,Periodic,CDF> mra_;
        Wavelet<T,Primal,Periodic,CDF> psi;
        flens::RefinementMatrix<T,Periodic,CDF> M1;
        
        Basis<T,Primal,Interval,Dijkema> refinementbasis;

    private:
        mutable FLENS_DEFAULT_INDEXTYPE _j;

        friend class Wavelet<T,Primal,Periodic,CDF>;


        flens::DenseVector<flens::Array<long double> > *_periodicRefCoeffs, *_rightRefCoeffs;

        FLENS_DEFAULT_INDEXTYPE 	*_innerOffsets, *_split;
        long double *_periodicL2Norms,  *_periodicH1SemiNorms;

        Basis(const Basis& secondbasis);

};

} // namespace lawa

#include <lawa/constructions/periodic/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H

