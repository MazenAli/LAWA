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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H
#define LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/refinementmatrix.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/settings/enum.h>

namespace flens {

template <typename T>
class RefinementMatrix<T,lawa::R,lawa::CDF>
    : public Matrix<RefinementMatrix<T,lawa::R,lawa::CDF> >
{
    public:
        typedef T ElementType;

        template <lawa::FunctionSide Side>
            RefinementMatrix(const lawa::BSpline<T,Side,lawa::R,lawa::CDF> &bw);

        template <lawa::FunctionSide Side>
            RefinementMatrix(const lawa::Wavelet<T,Side,lawa::R,lawa::CDF> &bw);

        flens::DenseVector<flens::Array<T> > band;
};

template <typename T>
struct TypeInfo<RefinementMatrix<T,lawa::R,lawa::CDF> >
{
    typedef RefinementMatrix<T,lawa::R,lawa::CDF> Impl;
};

//------------------------------------------------------------------------------

template <typename X, typename Y>
void
mv(cxxblas::Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,lawa::R,lawa::CDF> &A,
   const flens::DenseVector<X> &x, typename X::ElementType beta, flens::DenseVector<Y> &y);

} // namespace flens

#include <lawa/constructions/realline/refinementmatrix.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H

