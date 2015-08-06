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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Orthogonal,Interval,MultiRefinement>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = MultiRefinement;

        typedef BasisFunction<T,Orthogonal,Interval,MultiRefinement> BasisFunctionType;
        typedef BSpline<T,Orthogonal,Interval,MultiRefinement> BSplineType;

        Basis(const int d, const int j=-1);

        virtual
        ~Basis();

        int
        level() const;

        void
        setLevel(const int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const BasisFunctionType &
        generator(XType xtype) const;


        /// Returns the range of refinement B-splines from SecondRefinementBasis whose supports
        /// intersect the support of a given B-spline with level j_bspline1 and translation index
        /// k_bspline1 from the current Basis. This is required for tree-based algorithms.
        /// The level j_bspline1 of the B-splines is chosen s.t. there is "no scale difference", i.e.,
        /// if we refine both (possibly different) types of B-splines, the corresponding refinements
        /// should live on the same scale.
        template <typename SecondRefinementBasis>
            void
            getBSplineNeighborsForBSpline(int j_bspline1, long k_bspline1,
                                          const SecondRefinementBasis &secondrefinementbasis,
                                          int &j_wavelet, long &k_bspline2_first,
                                          long &k_bspline2_last) const;

        /// Returns the range of wavelets from SecondBasis whose supports intersect the support
        /// of a refinement B-spline with level j_bspline and translation index k_bspline
        /// from the current RefinementBasis. This is required for tree-based algorithms.
        /// The level j_wavelet of the wavelets is chosen s.t. there is "no scale difference", i.e.,
        /// if we refine both wavelets and refinement B-splines, the corresponding refinements should
        /// live on the same scale.
        template <typename SecondBasis>
            void
            getWaveletNeighborsForBSpline(int j_bspline, long k_bspline,
                                          const SecondBasis &basis,
                                          int &j_wavelet, long &k_wavelet_first,
                                          long &k_wavelet_last) const;

        class LaplaceOperator1D {
            public:
                LaplaceOperator1D(int _d,
                                  const Basis<_T,Orthogonal,Interval,MultiRefinement> &_refinementbasis);

                T
                operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2);

            private:
                int d;
                const Basis<_T,Orthogonal,Interval,MultiRefinement> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > outer_values;
                flens::DenseVector<flens::Array<long double> > inner_values1;
                flens::DenseVector<flens::Array<long double> > inner_values2;
        };

        class IdentityOperator1D {
            public:
                IdentityOperator1D(int _d,
                                   const Basis<_T,Orthogonal,Interval,MultiRefinement> &_refinementbasis);

                T
                operator()(XType xtype1, int j1, long k1, XType xtype2, int j2, long k2);

            private:
                int d;
                const Basis<_T,Orthogonal,Interval,MultiRefinement> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > outer_values;
                flens::DenseVector<flens::Array<long double> > inner_values1;
                flens::DenseVector<flens::Array<long double> > inner_values2;
        };


        MRA<T,Orthogonal,Interval,MultiRefinement> mra;

        const int d;
        const int j0;          // minimal used(!) level.

        mutable int _j;                // the current level.
        
        LaplaceOperator1D  LaplaceOp1D;
        IdentityOperator1D IdentityOp1D;
};

} // namespace lawa

#include <lawa/constructions/interval/multi/refinementbasis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_H
