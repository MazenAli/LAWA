/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_EXTENSIONS_FLENS_CG_H
#define LAWA_EXTENSIONS_FLENS_CG_H 1

#include <limits>
#include <flens/vectortypes/impl/densevector.h>

namespace lawa {

template <typename A>
struct _cg
{
    typedef A AuxVector;
    typedef typename A::ElementType T;
};

template <typename MA, typename VX, typename VB>
    FLENS_DEFAULT_INDEXTYPE
    cg(const MA &A, VX &x, const VB &b,
       typename _cg<VB>::T tol = std::numeric_limits<typename _cg<VB>::T>::epsilon(),
       FLENS_DEFAULT_INDEXTYPE maxIterations = std::numeric_limits<FLENS_DEFAULT_INDEXTYPE>::max());

template <typename MA, typename VX, typename VB>
    FLENS_DEFAULT_INDEXTYPE
    cgls(const MA &A, VX &x, const VB &b,
       typename _cg<VB>::T tol = std::numeric_limits<typename _cg<VB>::T>::epsilon(),
       FLENS_DEFAULT_INDEXTYPE maxIterations = std::numeric_limits<FLENS_DEFAULT_INDEXTYPE>::max());

template <typename Prec, typename MA, typename VA, typename VB>
    FLENS_DEFAULT_INDEXTYPE
    pcg(const Prec &P, const MA &A, VA &x, const VB &b,
        typename _cg<VB>::T tol = std::numeric_limits<typename _cg<VB>::T>::epsilon(),
        FLENS_DEFAULT_INDEXTYPE maxIterations = std::numeric_limits<FLENS_DEFAULT_INDEXTYPE>::max());

//-- trait specialization for FLENS types --------------------------------------

template <typename I>
struct _cg<flens::DenseVector<I> >
{
    typedef typename flens::DenseVector<I>::NoView AuxVector;
    typedef typename flens::DenseVector<I>::ElementType T;
};

} // namespace lawa

#include <extensions/flens/cg.tcc>

#endif // LAWA_EXTENSIONS_FLENS_CG_H

