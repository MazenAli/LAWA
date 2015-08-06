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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBSPLINE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,Interval,MultiRefinement>
    : public BasisFunction<_T,Orthogonal,Interval,MultiRefinement>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = MultiRefinement;

        BSpline(const MRA<T,Orthogonal,Interval,MultiRefinement> &mra);

        virtual
        ~BSpline();

        T
        operator()(T x, int j, long k, unsigned short deriv) const;

        Support<T>
        support(int j, long k) const;

        flens::DenseVector<flens::Array<T> >
        singularSupport(int j, long k) const;

        T
        tic(int j) const;

        flens::DenseVector<flens::Array<long double> > *
        getRefinement(int j, long k, int &refinement_j, long &refinement_k_first,
        		long &split, long &refinement_k_restart) const;

        int
        getRefinementLevel(int j) const;

        const MRA<T,Orthogonal,Interval,MultiRefinement> &mra;
        const unsigned int d;
        T initialticsize;
};

} // namespace lawa

#include <lawa/constructions/interval/multi/refinementbspline.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBSPLINE_H
