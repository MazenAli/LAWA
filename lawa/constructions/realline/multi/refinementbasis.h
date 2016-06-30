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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBASIS_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Orthogonal,R,MultiRefinement>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = MultiRefinement;

        typedef BasisFunction<T,Orthogonal,R,MultiRefinement>   BasisFunctionType;
        typedef BSpline<T,Orthogonal,R,MultiRefinement>         BSplineType;

        Basis(const FLENS_DEFAULT_INDEXTYPE d, const FLENS_DEFAULT_INDEXTYPE j);

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


        /// Returns the range of refinement B-splines from SecondRefinementBasis whose supports
        /// intersect the support of a given B-spline with level j_bspline1 and translation index
        /// k_bspline1 from the current Basis. This is required for tree-based algorithms.
        /// The level j_bspline1 of the B-splines is chosen s.t. there is "no scale difference", i.e.,
        /// if we refine both (possibly different) types of B-splines, the corresponding refinements
        /// should live on the same scale.
        template <typename SecondRefinementBasis>
            void
            getBSplineNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline1, FLENS_DEFAULT_INDEXTYPE k_bspline1,
                                          const SecondRefinementBasis &secondrefinementbasis,
                                          FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_bspline2_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_bspline2_last) const;

        /// Returns the range of wavelets from SecondBasis whose supports intersect the support
        /// of a refinement B-spline with level j_bspline and translation index k_bspline
        /// from the current RefinementBasis. This is required for tree-based algorithms.
        /// The level j_wavelet of the wavelets is chosen s.t. there is "no scale difference", i.e.,
        /// if we refine both wavelets and refinement B-splines, the corresponding refinements should
        /// live on the same scale.
        template <typename SecondBasis>
            void
            getWaveletNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline, FLENS_DEFAULT_INDEXTYPE k_bspline,
                                          const SecondBasis &basis,
                                          FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                          FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const;


        class LaplaceOperator1D {
            public:
                LaplaceOperator1D(FLENS_DEFAULT_INDEXTYPE _d,
                                  const Basis<_T,Orthogonal,R,MultiRefinement> &_refinementbasis);

                T
                operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2);

            private:
                FLENS_DEFAULT_INDEXTYPE d;
                const Basis<_T,Orthogonal,R,MultiRefinement> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > values1;
                flens::DenseVector<flens::Array<long double> > values2;
        };

        class IdentityOperator1D {
            public:
                IdentityOperator1D(FLENS_DEFAULT_INDEXTYPE _d,
                                   const Basis<_T,Orthogonal,R,MultiRefinement> &_refinementbasis);

                T
                operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2);

            private:
                FLENS_DEFAULT_INDEXTYPE d;
                const Basis<_T,Orthogonal,R,MultiRefinement> &refinementbasis;

                flens::DenseVector<flens::Array<long double> > values1;
                flens::DenseVector<flens::Array<long double> > values2;
        };

        LaplaceOperator1D  LaplaceOp1D;
        IdentityOperator1D IdentityOp1D;

        MRA<T,Orthogonal,R,MultiRefinement> mra;

        const FLENS_DEFAULT_INDEXTYPE d;
        const FLENS_DEFAULT_INDEXTYPE j0;          // minimal used(!) level.

        mutable FLENS_DEFAULT_INDEXTYPE _j;                // the current level.

};

} // namespace lawa

#include <lawa/constructions/realline/multi/refinementbasis.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_REFINEMENTBASIS_H
