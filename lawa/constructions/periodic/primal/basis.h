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

        Basis(int _d, int _d_, int j=0);

        ~Basis();

        int
        level() const;

        void
        setLevel(int j) const;

        const BasisFunctionType &
        generator(XType xtype) const;
        
        int
        cardJ(int j) const;

        int
        cardJL(int j = 0) const;

        int
        cardJI(int j) const;

        int
        cardJR(int j = 0) const;

        const flens::Range<int>
        rangeJ(int j) const;

        const flens::Range<int>
        rangeJL(int j) const;

        const flens::Range<int>
        rangeJI(int j) const;

        const flens::Range<int>
        rangeJR(int j) const;

        /// Returns the range of indicated functions from SecondBasis whose supports
        /// intersect the support of a given (multi-)scaling with level j_scaling and translation index
        /// k_scaling from the current Basis. This is required for tree-based algorithms.
        /// The returned level is chosen s.t. the corresponding refinements live
        /// on the same scale.
        void
        getScalingNeighborsForScaling(int j_scaling1, long k_scaling1,
        							  const Basis<T,Primal,Periodic,CDF> &secondbasis,
        							  int &j_scaling2, long &k_scaling_first,
        							  long &k_scaling_last) const;

        void
        getScalingNeighborsForScaling(int j_scaling1, long k_scaling1,
        							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
        							  int &j_scaling2, long &k_scaling_first,
        							  long &k_scaling_last) const;

        void
        getWaveletNeighborsForScaling(int j_scaling1, long k_scaling1,
        							  const Basis<T,Primal,Periodic,CDF> &secondbasis,
        							  int &j_wavelet, long &k_wavelet_first,
        							  long &k_wavelet_last) const;

        void
        getWaveletNeighborsForScaling(int j_scaling1, long k_scaling1,
        							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
        							  int &j_wavelet, long &k_wavelet_first,
        							  long &k_wavelet_last) const;

        /// Returns the range of indicated functions from SecondBasis and SecondRefinementBasis
        /// whose supports intersect the support of a given wavelet with level j_wavelet and
        /// translation index k_wavelet from the current Basis. This is required for tree-based algorithms.
        /// The returned level of the functions is chosen s.t. there is "no scale difference", i.e.,
        /// the corresponding refinements should live on the same scale.
        template <typename SecondRefinementBasis>
		void
		getBSplineNeighborsForWavelet(int j_wavelet, long k_wavelet,
									  const SecondRefinementBasis &secondrefinementbasis,
									  int &j_bspline, long &k_bspline_first,
									  long &k_bspline_last) const;

        template <typename SecondBasis>
		void
		getScalingNeighborsForWavelet(int j_wavelet, long k_wavelet,
									  const SecondBasis &secondbasis,
									  int &j_scaling, long &k_scaling_first,
									  long &k_scaling_last) const;

        void
        getWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                      const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                      int &j_wavelet2, long &k_wavelet_first,
                                      long &k_wavelet_last) const;

        void
        getWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                      const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                      int &j_wavelet2, long &k_wavelet_first,
                                      long &k_wavelet_last) const;


        template <typename SecondBasis>
		void
		getLowerWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
										   const SecondBasis &secondbasis,
										   int &j_wavelet2, long &k_wavelet_first,
										   long &k_wavelet_last) const;

		void
		getHigherWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
										   const Basis<T,Primal,Periodic,CDF> &secondbasis,
										   int &j_wavelet2, long &k_wavelet_first,
										   long &k_wavelet_last) const;

        void
        getHigherWaveletNeighborsForWavelet(int j_wavelet1, long k_wavelet1,
                                           const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                           int &j_wavelet2, long &k_wavelet_first,
                                           long &k_wavelet_last) const;

        const int d, d_, j0;
        MRA<T,Primal,Periodic,CDF> mra;
        MRA<T,Dual,Periodic,CDF> mra_;
        Wavelet<T,Primal,Periodic,CDF> psi;
        flens::RefinementMatrix<T,Periodic,CDF> M1;
        
        Basis<T,Primal,Interval,Dijkema> refinementbasis;

    private:
        mutable int _j;

        friend class Wavelet<T,Primal,Periodic,CDF>;


        flens::DenseVector<flens::Array<long double> > *_periodicRefCoeffs, *_rightRefCoeffs;

        long 	*_innerOffsets, *_split;
        long double *_periodicL2Norms,  *_periodicH1SemiNorms;

        Basis(const Basis& secondbasis);

};

} // namespace lawa

#include <lawa/constructions/periodic/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_H

