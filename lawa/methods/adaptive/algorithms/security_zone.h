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
 
#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H 1

#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

//Security zone for an 1d-index following Urban:2009, p.235 and KU:2010.
template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Side,Domain,Cons> &basis);

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis);
    
template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Side,Domain,Cons> &basis, const int Jmax);
    
template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis, const int Jmax);

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    void
    index_cone(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis,
               IndexSet<Index1D> &ret);

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
    void
    index_cone(const Index1D &lambda, T c, const Basis<T,Side,Domain,Cons> &basis,
               IndexSet<Index1D> &ret, const int Jmax);


// Special case: no bound for lower levels.
template <typename T>
    IndexSet<Index1D>
    C_WO_XBSpline(const IndexSet<Index1D> &Lambda, T c,
                  const Basis<T,Primal,R,CDF> &basis, bool only_pos=false);

template <typename T>
    IndexSet<Index1D>
    C_WO_XBSpline(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis);


// Computation of a security zone for 2d-tensor basis
template <typename T, typename Basis2D>
    IndexSet<Index2D>
    C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis, bool extralevel=false);
    
    // Security Zone does not include levels higher than J1_max, J2_max
template <typename T, typename Basis2D>
    IndexSet<Index2D>
    C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis, const int J1_max, const int J2_max );

template <typename T, typename Basis2D>
    IndexSet<Index2D>
    C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis);

template <typename T, typename Basis3D>
    IndexSet<Index3D>
    C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/security_zone.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H

