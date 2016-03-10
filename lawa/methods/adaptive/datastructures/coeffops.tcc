#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFOPS_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFOPS_TCC 1

#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <htucker/htucker.h>
#include <flens/flens.cxx>
#include <lawa/methods/adaptive/datastructures/indexops.h>

namespace lawa
{

template <SortingCriterion S, typename T, typename Index>
void
setCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::CoeffVec& _coeffs)
{
    assert(_coeffs.size()==coeffs.rank()*coeffs.dim());
    coeffs.getCoefficients() = _coeffs;
}


template <SortingCriterion S, typename T, typename Index>
void
setCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::size_type i,
                const typename SepCoefficients<S, T, Index>::
                size_type j,
                const typename SepCoefficients<S, T, Index>
                ::Coeff& coeff)
{
    assert(i>=1 && i<=coeffs.rank());
    assert(j>=1 && j<=coeffs.dim());
    coeffs.getCoefficients(i, j) = coeff;
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const SeparableRHSD<T, Basis>& rhs,
                const IndexSet<Index>& indexset)
{
    assert(coeffs.rank()==rhs.rank());
    assert(coeffs.dim()==rhs.dim());
    typedef typename SepCoefficients<S, T, Index>::size_type size_type;

    for (size_type i=1; i<=coeffs.rank(); ++i) {
        for (size_type j=1; j<=coeffs.dim(); ++j) {
            setCoefficients(coeffs, i, j, rhs(i, j, indexset));
        }
    }
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const SeparableRHSD<T, Basis>& rhs,
                const std::vector<IndexSet<Index>>& indexset)
{
    assert(coeffs.rank()==rhs.rank());
    assert(coeffs.dim()==rhs.dim());
    assert(indexset.size()==coeffs.rank()*coeffs.dim() ||
           indexset.size()==coeffs.dim());
    typedef typename SepCoefficients<S, T, Index>::size_type size_type;

    for (size_type i=1; i<=coeffs.rank(); ++i) {
        for (size_type j=1; j<=coeffs.dim(); ++j) {
            if (indexset.size()==coeffs.rank()*coeffs.dim()) {
                setCoefficients(coeffs, i, j, rhs(i, j,
                                indexset[(j-1)*coeffs.rank()+(i-1)]));
            } else {
                setCoefficients(coeffs, i, j, rhs(i, j,
                                indexset[j-1]));
            }
        }
    }
}


template <SortingCriterion S, typename T, typename Index>
std::ostream& operator<<(std::ostream& s,
                         const SepCoefficients<S, T, Index>& coeffs)
{
    typedef typename SepCoefficients<S, T, Index>::size_type size_type;

    for (size_type i=1; i<=coeffs.rank(); ++i) {
        for (size_type j=1; j<=coeffs.dim(); ++j) {
            s << "\n***Rank index " << i << ", dimension " << j << "***\n";
            s << coeffs.getCoefficients(i, j);
        }
    }

    return s;
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
unsigned long
maxintind(const Coefficients<S, T, Index>& coeffs, const Basis& basis)
{
    unsigned long max = 0;
    for (const auto& it : coeffs) {
        unsigned long idx = maptoint(it.first, basis);
        max = (idx>max) ? idx : max;
    }

    return max;
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree,
    const SepCoefficients<S, T, Index>& cp)
{
    assert((unsigned) tree.dim()==cp.dim());

    typedef typename SepCoefficients<S, T, Index>::size_type    size_type;
    typedef typename flens::DenseVector<flens::Array<T> >       DV;
    typedef typename flens::DenseVector
                     <flens::Array<unsigned long> >             IDV;
    typedef typename htucker::DenseVectorList<T>                DVList;

    DVList  list;
    IDV     sizes(tree.dim());

    for (size_type j=1; j<=cp.dim(); ++j) {
        unsigned long sizeij = 0;
        for (size_type i=1; i<=cp.rank(); ++i) {
            sizeij = maxintind(cp.getCoefficients(i, j), tree.basis());
            sizes(j) = (sizeij>sizes(j)) ? sizeij : sizes(j);
        }

        for (size_type i=1; i<=cp.rank(); ++i) {
            DV x(sizes(j));
            for (const auto& it : cp.getCoefficients(i, j)) {
                x(maptoint(it.first, tree.basis())) = it.second;
            }
            list.add(x);
        }
    }

    tree.tree().generateTofElementary(list, cp.rank(), cp.dim());
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const unsigned long col, const Coefficients<S, T, Index>& coeff)
{
    assert(idx.length()>0);
    assert(coeff.size()>0);
    assert(col>0);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >  Matrix;
    using flens::_;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            unsigned long rows = U.numRows();
            unsigned long cols = U.numCols();

            unsigned long max = maxintind(coeff, tree.basis());
            if (!rows || !cols) {
                U.resize(max, col);
            } else if (max>rows || col>cols) {
                Matrix copy(U);
                int sizer = (max>rows) ? max : rows;
                int sizec = (col>cols) ? col : cols;
                U.resize(sizer, sizec);
                U(_(1,rows),_(1,cols)) = copy;
            }

            for (const auto& it : coeff) {
                U(maptoint(it.first, tree.basis()), col) = it.second;
            }

            return;
        }
    }

    std::cerr << "error set: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const SepCoefficients<S, T, Index>& coeff)
{
    assert(idx.length()>0);
    assert(coeff.rank()>=1 && coeff.dim()==1);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;
    typedef typename SepCoefficients<S, T, Index>::size_type   size_type;
    using flens::_;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            size_type rowsnode = U.numRows();
            size_type colsnode = U.numCols();
            size_type colsset  = coeff.rank();
            size_type rowsset  = 0;
            for (size_type i=1; i<=colsset; ++i) {
                auto max = maxintind(coeff(i, 1), tree.basis());
                rowsset = (max>rowsset) ? max : rowsset;
            }

            if (!rowsnode || !colsnode) {
                U.resize(rowsset, colsset);
            } else if (rowsset>rowsnode || colsset>colsnode) {
                Matrix copy(U);
                size_type sizer = (rowsset>rowsnode) ? rowsset : rowsnode;
                size_type sizec = (colsset>colsnode) ? colsset : colsnode;
                U.resize(sizer, sizec);
                U(_(1,rowsnode),_(1,colsnode)) = copy;
            }

            for (size_type i=1; i<=colsset; ++i) {
                for (const auto& it : coeff(i, 1)) {
                    U(maptoint(it.first, tree.basis()), i) = it.second;
                }
            }

            return;
        }
    }

    std::cerr << "error set: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned long col, const T alpha,
     const Coefficients<S, T, Index>& coeff)
{
    assert(idx.length()>0);
    assert(coeff.size()>0);
    assert(col>0);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >  Matrix;
    using flens::_;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            unsigned long rows = U.numRows();
            unsigned long cols = U.numCols();

            unsigned long max = maxintind(coeff, tree.basis());
            if (!rows || !cols) {
                U.resize(max, col);
            } else if (max>rows || col>cols) {
                Matrix copy(U);
                int sizer = (max>rows) ? max : rows;
                int sizec = (col>cols) ? col : cols;
                U.resize(sizer, sizec);
                U(_(1,rows),_(1,cols)) = copy;
            }

            for (const auto& it : coeff) {
                U(maptoint(it.first, tree.basis()), col) += alpha*it.second;
            }

            return;
        }
    }

    std::cerr << "error axpy: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const T alpha, const SepCoefficients<S, T, Index>& coeff)
{
    assert(idx.length()>0);
    assert(coeff.rank()>=1 && coeff.dim()==1);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;
    typedef typename SepCoefficients<S, T, Index>::size_type   size_type;
    using flens::_;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            size_type rowsnode = U.numRows();
            size_type colsnode = U.numCols();
            size_type colsset  = coeff.rank();
            size_type rowsset  = 0;
            for (size_type i=1; i<=colsset; ++i) {
                auto max = maxintind(coeff(i, 1), tree.basis());
                rowsset = (max>rowsset) ? max : rowsset;
            }

            if (!rowsnode || !colsnode) {
                U.resize(rowsset, colsset);
            } else if (rowsset>rowsnode || colsset>colsnode) {
                Matrix copy(U);
                size_type sizer = (rowsset>rowsnode) ? rowsset : rowsnode;
                size_type sizec = (colsset>colsnode) ? colsset : colsnode;
                U.resize(sizer, sizec);
                U(_(1,rowsnode),_(1,colsnode)) = copy;
            }

            for (size_type i=1; i<=colsset; ++i) {
                for (const auto& it : coeff(i, 1)) {
                    U(maptoint(it.first, tree.basis()), i) += alpha*it.second;
                }
            }

            return;
        }
    }

    std::cerr << "error axpy: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx,
        const unsigned long col)
{
    assert(idx.length()>0);
    assert(col>0);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;

    Coefficients<Lexicographical, T, Index1D>   ret;
    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            const Matrix& U = tit.getNode()->getContent()
                              ->getUorB();
            assert(col<=(unsigned) U.numCols());
            for (unsigned long i=1; i<=(unsigned) U.numRows(); ++i) {
                if (U(i, col)!=(T) 0) {
                    ret[maptowav(i, tree.basis())] = U(i, col);
                }
            }

            return ret;
        }
    }

    std::cerr << "error extract: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx)
{
    assert(idx.length()>0);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            const Matrix& U = tit.getNode()->getContent()
                              ->getUorB();
            SepCoefficients<Lexicographical, T, Index1D>   ret(U.numCols(), 1);
            for (unsigned long i=1; i<=(unsigned) U.numRows(); ++i) {
                Index1D index = maptowav(i, tree.basis());
                for (unsigned long j=1; j<=(unsigned) U.numCols(); ++j) {
                    if (U(i, j)!=(T) 0) {
                        ret(j, 1)[index] = U(i, j);
                    }
                }
            }

            return ret;
        }
    }

    std::cerr << "error extract: idx not found\n";
    exit(EXIT_FAILURE);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFOPS_TCC 1
