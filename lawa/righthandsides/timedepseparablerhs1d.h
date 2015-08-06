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
 
#ifndef LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H
#define LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/functiontypes/functiontypes.h>
#include <lawa/integrals/integrals.h>

namespace lawa{

template<typename T, typename Basis1D>
class TimedepSeparableRHS1D
{
    public:
        const Basis1D& basis;
        const SeparableFunction2D<T>& F;
        
        IntegralF<Gauss, Basis1D>  integralf;
                
    public:
        TimedepSeparableRHS1D(const Basis1D& _basis, const SeparableFunction2D<T>& _F, int order);
               
        T
        operator()(T t, XType xtype, int j, long k) const;
                           
        T
        operator()(T t, const Index1D &index) const;
};


template <typename T, typename TimedepRHS>
class SumOfTimedepRHS1D
{
private:
    const TimedepRHS &rhs1;
    const TimedepRHS &rhs2;

public:
    SumOfTimedepRHS1D(const TimedepRHS &rhs1, const TimedepRHS &rhs2);
    
    T
    operator()(T t, XType xtype, int j, long k) const;
    
    T
    operator()(T t, const Index1D &index) const;
    
    Coefficients<Lexicographical,T,Index1D>
    operator()(T t, const IndexSet<Index1D> &Lambda) const;
};



} // namespace lawa

#include <lawa/righthandsides/timedepseparablerhs1d.tcc>

#endif // LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H
