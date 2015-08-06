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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/primal/bspline.h>
#include <lawa/constructions/periodic/refinementmatrix.h>

namespace lawa {

template <typename _T>
class MRA<_T,Primal,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Primal,Periodic,CDF> BasisFunctionType;
        typedef BSpline<T,Primal,Periodic,CDF> BSplineType;

        MRA(int _d, int _d_, int j=0);

        ~MRA();

        int
        level() const;

        void
        setLevel(int j) const;

        int
        cardI(int j) const;

        int
        cardIL(int j=0) const;

        int
        cardII(int j) const;

        int
        cardIR(int j=0) const;

        flens::Range<int>
        rangeI(int j) const;


        flens::Range<int>
        rangeIL(int j=0) const;

        flens::Range<int>
        rangeII(int j) const;

        flens::Range<int>
        rangeIR(int j) const;


        const int d, d_, j0, mu; // mu = mu(d) = d&1
        
        // do I need this?
        // const int l1, l2;

    protected:
        mutable int _j;

    private:

        friend class BSpline<T,Primal,Periodic,CDF>;

        flens::DenseVector<flens::Array<long double> > *_periodicRefCoeffs, *_rightRefCoeffs;

        long double *_periodicL2Norms,  *_periodicH1SemiNorms;

        long *_leftOffsets,
             *_periodicOffsets,
             *_rightOffsets,
             *_split;

        MRA(const MRA& secondmra);

    public:
        BSpline<T,Primal,Periodic,CDF> phi;
        flens::RefinementMatrix<T,Periodic,CDF> M0;
};

} // namespace lawa

#include <lawa/constructions/periodic/primal/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_H

