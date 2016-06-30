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
 
#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Primal,Interval,SparseMulti>
    : public BasisFunction<_T,Primal,Interval,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = SparseMulti;

        BSpline(const MRA<T,Primal,Interval,SparseMulti> &mra);

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

        const MRA<T,Primal,Interval,SparseMulti> &mra;
        const unsigned FLENS_DEFAULT_INDEXTYPE d;
};

} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_BSPLINE_H
