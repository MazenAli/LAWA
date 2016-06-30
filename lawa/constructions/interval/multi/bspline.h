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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,Interval,Multi>
    : public BasisFunction<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;

        BSpline(const MRA<T,Orthogonal,Interval,Multi> &mra);

        virtual
        ~BSpline();

        T
        operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const;

        Support<T>
        support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        flens::DenseVector<flens::Array<T> >
        singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        T
        tic(FLENS_DEFAULT_INDEXTYPE j) const;

        flens::DenseVector<flens::Array<long double> > *
        getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
        		FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const;

        FLENS_DEFAULT_INDEXTYPE
        getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const;

        T
        getL2Norm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        T
        getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;

        const MRA<T,Orthogonal,Interval,Multi> &mra;
        const unsigned FLENS_DEFAULT_INDEXTYPE d;
        T initialticsize;
};

} // namespace lawa

#include <lawa/constructions/interval/multi/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_H
