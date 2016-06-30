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

#include <cassert>
#include <lawa/math/lawa_math.h>

namespace flens {

template <typename T, lawa::Construction Cons>
RefinementMatrix<T,lawa::Interval,Cons>::RefinementMatrix()
{
}

template <typename T, lawa::Construction Cons>
RefinementMatrix<T,lawa::Interval,Cons>::RefinementMatrix(
                                FLENS_DEFAULT_INDEXTYPE nLeft, FLENS_DEFAULT_INDEXTYPE nRight,
                                const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A,
                                FLENS_DEFAULT_INDEXTYPE _min_j0, FLENS_DEFAULT_INDEXTYPE cons_j)
    : left(nLeft, (nLeft>0) ? A.firstCol() : 1),
      right(nRight, (nRight>0) ? A.lastCol()-nRight+1 : A.lastCol()+2),
      lengths(_(-nRight, nLeft)),
      min_j0(_min_j0), _cons_j(cons_j), _j(cons_j),
      _firstRow(A.firstRow()), _firstCol(A.firstCol()),
      _lastRow(A.lastRow()), _lastCol(A.lastCol()),
      _additionalRows(0), _additionalCols(0)
{
    assert(nLeft>=0);
    assert(nRight>=0);
    assert(_cons_j>=min_j0);

    assert(_firstCol>=1);

    _extractMasks(A);

    for (FLENS_DEFAULT_INDEXTYPE i=left.firstIndex(); i<=left.lastIndex(); ++i) {
        lengths(1+i-left.firstIndex()) = left(i).length()+(A.firstRow()-1);
    }

    lengths(0) = leftband.firstIndex() - 1;//-A.firstRow();

    for (FLENS_DEFAULT_INDEXTYPE i=right.firstIndex(); i<=right.lastIndex(); ++i) {
        lengths(-nRight+i-right.firstIndex()) = right(i).length()+(A.firstRow()-1);
    }
}


template <typename T, lawa::Construction Cons>
const typename flens::DenseVector<flens::Array<T> >::ConstView
RefinementMatrix<T,lawa::Interval,Cons>::operator()(FLENS_DEFAULT_INDEXTYPE j, const Underscore<FLENS_DEFAULT_INDEXTYPE> &/*u*/,
                                              FLENS_DEFAULT_INDEXTYPE col) const
{
    assert(j>=min_j0);

    FLENS_DEFAULT_INDEXTYPE additionalCols = 0, additionalRows = 0;
    if (j>_cons_j) {
        for (FLENS_DEFAULT_INDEXTYPE l=_cons_j; l<j; ++l) {
            additionalCols += lawa::pow2i<T>(l);
        }
        additionalRows = 2*additionalCols;
    } else if (j<_cons_j) {
        for (FLENS_DEFAULT_INDEXTYPE l=_cons_j-1; l>=j; --l) {
            additionalCols -= lawa::pow2i<T>(l);
        }
        additionalRows = 2*additionalCols;
    }

    assert(col>=_firstCol);
    assert(col<=_lastCol+additionalCols);

    if (col<=left.lastIndex()) {
        return left(col);
    }

    if (col>=right.firstIndex()+additionalCols) {
        const flens::DenseVector<flens::Array<T> > &rightCol = right(col-additionalCols);
        return rightCol( _ , rightCol.firstIndex()+additionalRows);
    }

    return (col>(_firstCol+_lastCol+additionalCols)/2) ?
          rightband( _ ,leftband.firstIndex() + 2*(col-left.lastIndex()-1))
        : leftband( _ , leftband.firstIndex() + 2*(col-left.lastIndex()-1));
}

template <typename T, lawa::Construction Cons>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
RefinementMatrix<T,lawa::Interval,Cons>::rows() const
{
    return _(firstRow(), lastRow());
}

template <typename T, lawa::Construction Cons>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
RefinementMatrix<T,lawa::Interval,Cons>::cols() const
{
    return _(firstCol(), lastCol());
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::numRows() const
{
    return lastRow()-firstRow()+1;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::numCols() const
{
    return lastCol()-firstCol()+1;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::firstRow() const
{
    return _firstRow;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::lastRow() const
{
    return _lastRow + _additionalRows;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::firstCol() const
{
    return _firstCol;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::lastCol() const
{
    return _lastCol + _additionalCols;
}

template <typename T, lawa::Construction Cons>
FLENS_DEFAULT_INDEXTYPE
RefinementMatrix<T,lawa::Interval,Cons>::level() const
{
    return _j;
}

// TODO: consider setLevel as private(!) friend method or mra/basis.
template <typename T, lawa::Construction Cons>
void
RefinementMatrix<T,lawa::Interval,Cons>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    if (j<_j) {
        assert(j>=min_j0);
        for (FLENS_DEFAULT_INDEXTYPE l=_j-1; l>=j; --l) {
            _additionalCols -= lawa::pow2i<T>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
    if (j>_j) {
        for (FLENS_DEFAULT_INDEXTYPE l=_j; l<j; ++l) {
            _additionalCols += lawa::pow2i<T>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
}

template <typename T, lawa::Construction Cons>
void
RefinementMatrix<T,lawa::Interval,Cons>::_extractMasks(
                                    const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &A)
{
    // extract left block
    for (FLENS_DEFAULT_INDEXTYPE c=A.firstCol(); c<A.firstCol()+left.length(); ++c) {
        FLENS_DEFAULT_INDEXTYPE r = A.lastRow();
        while (fabs(A(r,c))<=1e-12) {
            --r;
            assert(r>=A.firstRow());
        }
        left(c) = A(_(A.firstRow(),r),c);
        left(c).engine().changeIndexBase(A.firstRow());
    }
    
    // extract right block
    for (FLENS_DEFAULT_INDEXTYPE c=A.lastCol()-right.length()+1; c<=A.lastCol(); ++c) {
        FLENS_DEFAULT_INDEXTYPE r = A.firstRow();
        while (fabs(A(r,c))<=1e-12) {
            ++r;
            assert(r<=A.lastRow());
        }
        right(c) = A(_(r,A.lastRow()), c);
        right(c).engine().changeIndexBase(r);
    }
    
    // extract band (left to middle)
    FLENS_DEFAULT_INDEXTYPE c = A.firstCol()+left.length();
    FLENS_DEFAULT_INDEXTYPE first = A.firstRow();
    while (fabs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }    
    FLENS_DEFAULT_INDEXTYPE last = A.lastRow();
    while (fabs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    leftband = A(_(first,last), c);
    leftband.engine().changeIndexBase(first);
#ifdef CHECK_INTERVAL_CONSTRUCTION
    for (++c; c<=(A.firstCol()+A.lastCol())/2; ++c) {
        FLENS_DEFAULT_INDEXTYPE i=leftband.firstIndex();
        first += 2; last += 2;
        for (FLENS_DEFAULT_INDEXTYPE r=first; r<=last; ++r, ++i) {
            assert(fabs(leftband(i)-A(r,c))<=1e-12);
        }
    }
#endif

    // extract band (middle to right)
    c = A.lastCol()-right.length();
    first = A.firstRow();
    while (fabs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }
    last = A.lastRow();
    while (fabs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    rightband = A(_(first,last), c);
    rightband.engine().changeIndexBase(first);
#ifdef CHECK_INTERVAL_CONSTRUCTION
    for (--c; c>(A.firstCol()+A.lastCol())/2; --c) {
        FLENS_DEFAULT_INDEXTYPE i=rightband.firstIndex();
        first -= 2; last -= 2;
        for (FLENS_DEFAULT_INDEXTYPE r=first; r<=last; ++r, ++i) {
            assert(fabs(rightband(i)-A(r,c))<=1e-12);
        }
    }
    assert(leftband.length()==rightband.length());
//    assert(leftband.firstIndex()-A.firstRow()==rightband.lastIndex()-A.lastRow());
#endif

}

//------------------------------------------------------------------------------

template <typename X, lawa::Construction Cons, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,lawa::Interval,Cons> &A,
   const flens::DenseVector<X> &x, typename X::ElementType beta, flens::DenseVector<Y> &y)
{
    typedef typename X::ElementType T;
    assert(alpha==T(1));
    assert(x.engine().stride()==1);
    assert(y.engine().stride()==1);

    if (transA==cxxblas::NoTrans) {
        assert(A.numCols()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.rows()) || y.engine().fill(T(0));
        } else {
            assert(y.length()==A.numRows());
            y.engine().changeIndexBase(A.firstRow());
        }

        // left upper block
        FLENS_DEFAULT_INDEXTYPE ix = x.firstIndex();
        for (FLENS_DEFAULT_INDEXTYPE c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++ix) {
            FLENS_DEFAULT_INDEXTYPE n = A.left(c).length();
            cxxblas::axpy((int)n,
                          x(ix),
                          A.left(c).engine().data(), 1,
                          y.engine().data(), 1);
        }

        // central band (up to middle)
        FLENS_DEFAULT_INDEXTYPE iy = A.leftband.firstIndex()-A.firstRow();
        FLENS_DEFAULT_INDEXTYPE n = A.leftband.length();
        FLENS_DEFAULT_INDEXTYPE middle = lawa::iceil<FLENS_DEFAULT_INDEXTYPE>(x.length()/2.);
        for (FLENS_DEFAULT_INDEXTYPE c=A.left.lastIndex()+1; c<=middle; ++c, iy+=2, ++ix) {
            cxxblas::axpy((int)n,
                          x(ix),
                          A.leftband.engine().data(), 1,
                          y.engine().data()+iy, 1);
        }
        // central band (right of middle)
        FLENS_DEFAULT_INDEXTYPE end = A.left.firstIndex() + x.length() - A.right.length();
        for (FLENS_DEFAULT_INDEXTYPE c=middle+1; c<end; ++c, iy+=2, ++ix) {
            cxxblas::axpy((int)n,
                          x(ix),
                          A.rightband.engine().data(), 1,
                          y.engine().data()+iy, 1);
        }

        // right lower block
        for (FLENS_DEFAULT_INDEXTYPE c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++ix) {
            FLENS_DEFAULT_INDEXTYPE n = A.right(c).length();
            cxxblas::axpy((int)n, 
                          x(ix),
                          A.right(c).engine().data(), 1,
                          y.engine().data()+y.length()-1-n+1, 1);
        }
    } else { // transA==Trans
        assert(A.numRows()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.cols());
        } else {
            assert(y.length()==A.numCols());
            y.engine().changeIndexBase(A.firstCol());
        }
        FLENS_DEFAULT_INDEXTYPE iy = y.firstIndex();
        // left upper block
        for (FLENS_DEFAULT_INDEXTYPE c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++iy) {
            FLENS_DEFAULT_INDEXTYPE n = A.left(c).length();
            cxxblas::dot((int)n,
                         A.left(c).engine().data(), 1,
                         x.engine().data(), 1, 
                         y(iy));
        }

        // central band (up to middle)
        FLENS_DEFAULT_INDEXTYPE middle = y.length()/2;
        FLENS_DEFAULT_INDEXTYPE ix = A.leftband.firstIndex() - A.firstRow();
        for (FLENS_DEFAULT_INDEXTYPE i=A.left.length()+1; i<=middle; ++i, ix+=2, ++iy) {
            cxxblas::dot((int) A.leftband.length(),
                         A.leftband.engine().data(), 1,
                         x.engine().data()+ix, 1,
                         y(iy));
        }
        // central band (right of middle)
        FLENS_DEFAULT_INDEXTYPE end = y.length() - A.right.length();
        for (FLENS_DEFAULT_INDEXTYPE i=middle+1; i<=end; ++i, ix+=2, ++iy) {
            cxxblas::dot((int) A.rightband.length(),
                         A.rightband.engine().data(), 1,
                         x.engine().data()+ix, 1,
                         y(iy));
        }
        // right lower block
        for (FLENS_DEFAULT_INDEXTYPE c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++iy) {
            FLENS_DEFAULT_INDEXTYPE n = A.right(c).length();
            cxxblas::dot((int)n,
                         A.right(c).engine().data(), 1,
                         x.engine().data() + A.numRows() - n, 1, 
                         y(iy));
        }
    }
}

} // namespace flens

