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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/periodic/periodicsupport.h>
#include <lawa/constructions/periodic/primal/bspline.h>
#include <lawa/constructions/periodic/dual/bspline.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
struct Wavelet<_T,Primal,Periodic,CDF>
    : public BasisFunction<_T,Primal,Periodic,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = Periodic;
    static const Construction Cons = CDF;

   // Wavelet(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_);

   // Wavelet(const BSpline<T,Primal,Periodic,CDF> &_phi,
   //         const BSpline<T,Dual,Periodic,CDF> &_phi_);

    Wavelet(const Basis<T,Primal,Periodic,CDF> &_basis);

    T
    operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const;

    PeriodicSupport<T>
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

    const flens::DenseVector<flens::Array<T> > &
    mask() const;

    static flens::DenseVector<flens::Array<T> >
    mask(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE d_);

    const FLENS_DEFAULT_INDEXTYPE d, d_, mu;
    const FLENS_DEFAULT_INDEXTYPE vanishingMoments;
    const Wavelet<T, Primal, R, CDF> psiR;

    const Basis<T,Primal,Periodic,CDF> &basis;

private:
    Wavelet(const Wavelet& secondwavelet);
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_WAVELET_H

