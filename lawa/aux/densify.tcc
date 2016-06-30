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

namespace lawa {

template <typename I>
void
densify(cxxblas::Transpose trans, const flens::Matrix<I> &A,
        flens::GeMatrix<flens::FullStorage<typename I::ElementType, cxxblas::ColMajor> > &D,
        FLENS_DEFAULT_INDEXTYPE firstRow, FLENS_DEFAULT_INDEXTYPE firstCol)
{
    using flens::_;

    if (trans==cxxblas::NoTrans) {
        FLENS_DEFAULT_INDEXTYPE m = A.impl().numRows();
        FLENS_DEFAULT_INDEXTYPE n = A.impl().numCols();
        D.engine().resize(m,n,firstRow,firstCol);
        flens::DenseVector<flens::Array<typename I::ElementType> > e(n);
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=n; ++i) {
            flens::DenseVector<flens::Array<typename I::ElementType> > y(m);
            e(i) = 1.;
            D(_,i+firstCol-1) = A.impl()*e;
            e(i) = 0.;
        }
    } else {
        assert(trans==cxxblas::Trans);

        FLENS_DEFAULT_INDEXTYPE m = A.impl().numRows();
        FLENS_DEFAULT_INDEXTYPE n = A.impl().numCols();
        D.engine().resize(n,m,firstCol,firstRow);
        flens::DenseVector<flens::Array<typename I::ElementType> > e(m);
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=m; ++i) {
            flens::DenseVector<flens::Array<typename I::ElementType> > y(n);
            e(i) = 1.;
            D(_,i+firstCol-1) = transpose(A.impl())*e;
            e(i) = 0.;
        }
    }
}

} // namespace lawa

