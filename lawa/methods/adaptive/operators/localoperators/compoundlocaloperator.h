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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_COMPOUNDLOCALOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_COMPOUNDLOCALOPERATOR_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator=SecondLocalOperator,
          typename FourthLocalOperator=SecondLocalOperator>
class CompoundLocalOperator {

    public:
        typedef typename FirstLocalOperator::T T;
        typedef typename Coefficients<Lexicographical,T,Index>::iterator    coeff_it;
        typedef typename Coefficients<Lexicographical,T,Index>::iterator    const_coeff_it;
        typedef typename IndexSet<Index1D>::const_iterator                  const_set1d_it;

        CompoundLocalOperator(FirstLocalOperator  &_firstLocalOp,
                              SecondLocalOperator &_secondLocalOp);

        CompoundLocalOperator(FirstLocalOperator  &_firstLocalOp,
                              SecondLocalOperator &_secondLocalOp,
                              ThirdLocalOperator  &_thirdLocalOp);
        
        CompoundLocalOperator(FirstLocalOperator  &_firstLocalOp,
                              SecondLocalOperator &_secondLocalOp,
                              ThirdLocalOperator  &_thirdLocalOp,
                              FourthLocalOperator  &_fourthLocalOp);

        void
        eval(const Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av);

        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, const char* evalType);

        template <typename Preconditioner>
        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P);

        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, Coefficients<Lexicographical,T,Index> &P);

        template <typename RightPrec, typename LeftPrec>
        void
        eval(Coefficients<Lexicographical,T,Index> &v,
        	 Coefficients<Lexicographical,T,Index> &Av, RightPrec &rightP, LeftPrec &leftP);

        template <typename Preconditioner>
        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, FLENS_DEFAULT_INDEXTYPE operatornumber);

        template <typename Preconditioner>
        void
        eval(Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, const char* evalType);

        template <typename Preconditioner>
        void
        apply(Coefficients<Lexicographical,T,Index> &v,
              Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P, T eps);

        T
        eval(const Index& ind_row, const Index& ind_col);
        
        void
        clear();

        FLENS_DEFAULT_INDEXTYPE                         numOfLocalOp;
        FirstLocalOperator          &firstLocalOp;
        SecondLocalOperator         &secondLocalOp;
        ThirdLocalOperator          &thirdLocalOp;
        FourthLocalOperator         &fourthLocalOp;
};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/localoperators/compoundlocaloperator.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_COMPOUNDLOCALOPERATOR_H
