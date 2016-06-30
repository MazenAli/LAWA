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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Primal,R,CDF>
    : public BasisFunction<_T,Primal,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    BSpline(FLENS_DEFAULT_INDEXTYPE _d);

    BSpline(const MRA<T,Primal,R,CDF> &mra);

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

    const flens::DenseVector<flens::Array<T> > &
    mask() const;

    const FLENS_DEFAULT_INDEXTYPE d, mu;
    const FLENS_DEFAULT_INDEXTYPE l1, l2;
    const flens::DenseVector<flens::Array<T> > a;
};

//------------------------------------------------------------------------------

template <typename T>
BSpline<T,Primal,R,CDF>
N(FLENS_DEFAULT_INDEXTYPE d);

template <typename T>
flens::DenseVector<flens::Array<T> >
N1();

template <typename T>
flens::DenseVector<flens::Array<T> >
N2();

template <typename T>
flens::DenseVector<flens::Array<T> >
N3();

template <typename T>
flens::DenseVector<flens::Array<T> >
N4();

template <typename T>
flens::DenseVector<flens::Array<T> >
N5();

template <typename T>
flens::DenseVector<flens::Array<T> >
N6();

template <typename T>
flens::DenseVector<flens::Array<T> >
N7();

template <typename T>
flens::DenseVector<flens::Array<T> >
N8();

template <typename T>
flens::DenseVector<flens::Array<T> >
N9();

template <typename T>
flens::DenseVector<flens::Array<T> >
N10();

} // namespace lawa

#include <lawa/constructions/realline/primal/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H

