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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/periodic/periodicsupport.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Primal,Periodic,CDF>
    : public BasisFunction<_T,Primal,Periodic,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = Periodic;
    static const Construction Cons = CDF;

    //BSpline(FLENS_DEFAULT_INDEXTYPE _d);                        // Periodized Bspline on [0,1]

    BSpline(const MRA<T,Primal,Periodic,CDF> &_mra);

    virtual
    ~BSpline();

    T
    operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const;    // (cumulated) evaluation on [0,1]

    PeriodicSupport<T>
    support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;            // subset of [0,1], may be divided in 2 parts

    flens::DenseVector<flens::Array<T> >
    singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;    // Points in [0,1]

    T
/**/tic(FLENS_DEFAULT_INDEXTYPE j) const;                       // -> ???

    /**
     * Get coefficients in _DIJKEMA_ Basis (not really a refinement here!!)
     */
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
    mask() const;                           // mask of original Spline on R

    const FLENS_DEFAULT_INDEXTYPE d, mu;
    const BSpline<T, Primal, R, CDF> phiR;  // ''original'' Bspline on R

    const MRA<T,Primal,Periodic,CDF> &mra;

private:
    BSpline(const BSpline& secondbspline);
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_H

