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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Dual,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Dual;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef BasisFunction<T,Dual,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Dual,Interval,Dijkema> BSplineType;
        typedef Wavelet<T,Dual,Interval,Dijkema> WaveletType;

        Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j=-1);

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
        cardJ_(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJ_L(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJ_I(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardJ_R(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ_(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ_L(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ_I(FLENS_DEFAULT_INDEXTYPE j) const;

        const flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeJ_R(FLENS_DEFAULT_INDEXTYPE j=-1) const;

        MRA<T,Primal,Interval,Dijkema> mra;
        MRA<T,Dual,Interval,Dijkema>  mra_;

        flens::RefinementMatrix<T,Interval,Dijkema> M1_;

        const FLENS_DEFAULT_INDEXTYPE d, d_, mu;   // mu = mu(d) = d&1.
        const FLENS_DEFAULT_INDEXTYPE min_j0;      // minimal allowed(!) level;
        const FLENS_DEFAULT_INDEXTYPE j0;          // minimal used(!) level. 

    private:
        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable FLENS_DEFAULT_INDEXTYPE _j;                // the current level.

    public:
        Wavelet<T,Dual,Interval,Dijkema> psi_;
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/dual/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BASIS_H

