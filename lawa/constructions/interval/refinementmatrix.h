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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H
#define LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/refinementmatrix.h>

namespace flens {

template <typename T, lawa::Construction Cons>
class RefinementMatrix<T, lawa::Interval, Cons>
    : public Matrix<RefinementMatrix<T, lawa::Interval, Cons> >
{
    public:
        typedef T ElementType;
                
        RefinementMatrix();
        
        RefinementMatrix(FLENS_DEFAULT_INDEXTYPE nLeft, FLENS_DEFAULT_INDEXTYPE nRight,
                         const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A,
                         FLENS_DEFAULT_INDEXTYPE _min_j0, FLENS_DEFAULT_INDEXTYPE cons_j);

        const typename flens::DenseVector<flens::Array<T> >::ConstView
        operator()(FLENS_DEFAULT_INDEXTYPE j, const Underscore<FLENS_DEFAULT_INDEXTYPE> &u, FLENS_DEFAULT_INDEXTYPE col) const;
        
        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        rows() const;

        flens::Range<FLENS_DEFAULT_INDEXTYPE>
        cols() const;

        FLENS_DEFAULT_INDEXTYPE
        numRows() const;

        FLENS_DEFAULT_INDEXTYPE
        numCols() const;

        FLENS_DEFAULT_INDEXTYPE
        firstRow() const;

        FLENS_DEFAULT_INDEXTYPE
        lastRow() const;

        FLENS_DEFAULT_INDEXTYPE
        firstCol() const;

        FLENS_DEFAULT_INDEXTYPE
        lastCol() const;

        FLENS_DEFAULT_INDEXTYPE
        level() const;

        void
        setLevel(FLENS_DEFAULT_INDEXTYPE j) const;

        flens::DenseVector<flens::Array<flens::DenseVector<flens::Array<T> > > > left, right;
        flens::DenseVector<flens::Array<T> > leftband, rightband;
        flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > lengths;
        FLENS_DEFAULT_INDEXTYPE min_j0;

    private:
        void
        _extractMasks(const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &A);        
        
        FLENS_DEFAULT_INDEXTYPE _cons_j;
        mutable FLENS_DEFAULT_INDEXTYPE _j;
        FLENS_DEFAULT_INDEXTYPE _firstRow, _firstCol, _lastRow, _lastCol;
        mutable FLENS_DEFAULT_INDEXTYPE _additionalRows, _additionalCols;
};

template <typename T, lawa::Construction Cons>
struct TypeInfo<RefinementMatrix<T,lawa::Interval,Cons> >
{
    typedef RefinementMatrix<T,lawa::Interval,Cons> Impl;
    typedef T                                 ElementType;
    
};

template <typename X, lawa::Construction Cons, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,lawa::Interval,Cons> &A,
   const flens::DenseVector<X> &x, typename X::ElementType beta, flens::DenseVector<Y> &y);

} // namespace flens

#include <lawa/constructions/interval/refinementmatrix.tcc>

#endif // LAWA_CONSTRUCTIONS_REFINEMENTMATRIX_H

