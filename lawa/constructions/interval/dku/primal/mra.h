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
#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DKU_PRIMAL_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DKU_PRIMAL_MRA_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/interval/refinementmatrix.h>

namespace lawa {

template <typename _T>
class MRA<_T,Primal,Interval,DKU>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = DKU;

        typedef BasisFunction<T,Primal,Interval,DKU> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,DKU> BSplineType;

        MRA(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE d_, FLENS_DEFAULT_INDEXTYPE j=-1);

        ~MRA();

        // cardinalities of whole, left, inner, right index sets.
        FLENS_DEFAULT_INDEXTYPE
        cardI(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardIL(FLENS_DEFAULT_INDEXTYPE j=0) const;

        FLENS_DEFAULT_INDEXTYPE
        cardII(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        cardIR(FLENS_DEFAULT_INDEXTYPE j=0) const;

        // ranges of whole left, inner, right index sets.
        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeI(FLENS_DEFAULT_INDEXTYPE j) const;

        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeIL(FLENS_DEFAULT_INDEXTYPE j=0) const;

        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeII(FLENS_DEFAULT_INDEXTYPE j) const;

        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rangeIR(FLENS_DEFAULT_INDEXTYPE j) const;

        FLENS_DEFAULT_INDEXTYPE
        level() const;

        void
        setLevel(FLENS_DEFAULT_INDEXTYPE j);

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const FLENS_DEFAULT_INDEXTYPE d, d_, mu;    // mu = mu(d) = d&1.
        const FLENS_DEFAULT_INDEXTYPE l1, l2;     // support phi  = [ l1, l2 ] (real line).
        const FLENS_DEFAULT_INDEXTYPE l1_, l2_;   // support phi_  = [ l1_, l2_ ] (real line).
        const FLENS_DEFAULT_INDEXTYPE min_j0;   // minimal allowed(!) level;
        const FLENS_DEFAULT_INDEXTYPE j0;       // minimal used(!) level.

//        BSpline<T,Primal,R,CDF> phiR;
//        BSpline<T,Primal,R,CDF> phi_R;
//        flens::RefinementMatrix<T,Interval,Primbs> M0;

        BSpline<T,Primal,Interval,DKU> phi;
        flens::RefinementMatrix<T,Interval,DKU> M0;
    private:
        void
        _alpha_();

        void
        _beta_();
    
        void
        _calcM0();

        const FLENS_DEFAULT_INDEXTYPE l, q;        // defining index sets left and right, i.e.
                               // I_L = { 1-l2, ..., l+1 } with l >= -l1 and
                               // I_R = { 2^-q, ..., 2^-l1-1 } with q >= l2.

        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > _bc;    // the boundary conditions
                                         // bc(0) = 1 -> Dirichlet BC left.
                                         // bc(1) = 1 -> Dirichlet BC right.

        FLENS_DEFAULT_INDEXTYPE _j;                // the current level.

                                                    // helper matrices for
                                                    // reproduction of monomials
                                                    // at boundaries.
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _Alpha_, _Beta_;

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _ML, _MR;  // upper left,
                                                      // lower right block of
                                                      // refinement matrix M0.
//        flens::DenseVector<flens::Array<T> > _a; // refinement coefficients of inner
//                                   // scaling functions.
};

} // namespace lawa

#include <lawa/constructions/interval/dku/primal/mra.tcc>

#endif //LAWA_CONSTRUCTIONS_INTERVAL_DKU_PRIMAL_MRA_H

