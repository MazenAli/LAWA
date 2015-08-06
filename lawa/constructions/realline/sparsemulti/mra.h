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
 
#ifndef LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H
#define LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>

namespace lawa {

template <typename _T>
class MRA<_T,Primal,R,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = R;
        static const Construction Cons = SparseMulti;

        typedef BasisFunction<T,Primal,R,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,R,SparseMulti>       BSplineType;

        MRA(int _d, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        const int d;
        const int j0;          // minimal used(!) level.

        BSpline<T,Primal,R,SparseMulti> phi;

    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/realline/sparsemulti/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_SPARSEMULTI_MRA_H

