#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_TCC 1

#define _USE_MATH_DEFINES
#ifndef MAX
    #define MAX(x, y) (x>y) ? x : y
#endif

#ifndef MIN
    #define MIN(x, y) (x<y) ? x : y
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <htucker/htucker.h>
#include <flens/flens.cxx>
#include <lawa/methods/adaptive/algorithms/indexops.h>

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


template <SortingCriterion S, typename T, typename Index>
void
updateCoefficients(SepCoefficients<S, T, Index>& coeffs,
                   const typename SepCoefficients<S, T, Index>
                   ::size_type i,
                   const typename SepCoefficients<S, T, Index>::
                   size_type j,
                   const typename SepCoefficients<S, T, Index>
                   ::Coeff& coeff)
{
    assert(i>=1 && i<=coeffs.rank());
    assert(j>=1 && j<=coeffs.dim());
    coeffs.getCoefficients(i, j).update(coeff);
}


template <SortingCriterion S, typename T, typename Index>
void
addCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::size_type i,
                const typename SepCoefficients<S, T, Index>::
                size_type j,
                const typename SepCoefficients<S, T, Index>
                ::Coeff& coeff)
{
    assert(i>=1 && i<=coeffs.rank());
    assert(j>=1 && j<=coeffs.dim());
    coeffs.getCoefficients(i, j) += coeff;
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
genCoefficients(      SepCoefficients<S, T, Index>& coeffs,
                      SeparableRHSD<T, Basis>& rhs,
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
void
genCoefficientsRnd(SepCoefficients<S, T, Index>&       coeffs,
                   const std::vector<IndexSet<Index>>& indexset,
                   const T                             scale,
                   const long                          seed)
{
    assert(coeffs.dim()==indexset.size());
    typedef typename SepCoefficients<S, T, Index>::size_type size_type;

    /* Random number generator */
    unsigned actual_seed;
    if (seed<0) actual_seed = time(NULL);
    else actual_seed = (unsigned) seed;
    srand(actual_seed);

    for (size_type i=1; i<=coeffs.rank(); ++i) {
        for (size_type j=1; j<=coeffs.dim(); ++j) {
            Coefficients<S, T, Index> x;
            for (auto& lambda : indexset[j-1]) {
                T value         = (T) rand()/(T) RAND_MAX;
                value          *= scale;
                x[lambda]       = value;
            }
            setCoefficients(coeffs, i, j, x);
        }
    }
}


template <SortingCriterion S, typename T, typename Basis>
void
setCoefficientsJ0(      SepCoefficients<S, T, Index1D>& coeffs,
                  const Basis&                        basis)
{
    typedef typename SepCoefficients<S, T, Index1D>::size_type size_type;

    for (size_type i=1; i<=coeffs.rank(); ++i) {
        for (size_type j=1; j<=coeffs.dim(); ++j) {
            Coefficients<S, T, Index1D> x;
            auto range = basis.rangeJ(basis.j0);
            Index1D mu(basis.j0, range.firstIndex(), XBSpline);
            x[mu] = 1.;
            setCoefficients(coeffs, i, j, x);
        }
    }
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genAddCoefficients(      SepCoefficients<S, T, Index>& coeffs,
                         SeparableRHSD<T, Basis>& rhs,
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
                addCoefficients(coeffs, i, j, rhs(i, j,
                                indexset[(j-1)*coeffs.rank()+(i-1)]));
            } else {
                addCoefficients(coeffs, i, j, rhs(i, j,
                                indexset[j-1]));
            }
        }
    }
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genUpdateCoefficients(      SepCoefficients<S, T, Index>& coeffs,
                            SeparableRHSD<T, Basis>& rhs,
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
                updateCoefficients(coeffs, i, j, rhs(i, j,
                                   indexset[(j-1)*coeffs.rank()+(i-1)]));
            } else {
                updateCoefficients(coeffs, i, j, rhs(i, j,
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
unsigned FLENS_DEFAULT_INDEXTYPE
maxintind(const Coefficients<S, T, Index>& coeffs, const Basis& basis)
{
    unsigned FLENS_DEFAULT_INDEXTYPE max = 0;
    for (const auto& it : coeffs) {
        unsigned FLENS_DEFAULT_INDEXTYPE idx = maptoint(it.first, basis);
        max = (idx>max) ? idx : max;
    }

    return max;
}


template <SortingCriterion S, typename T, typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintindhash(const Coefficients<S, T, Index>& coeffs,
              const int dim,
              HTCoefficients<T, Basis>& u)
{
    unsigned FLENS_DEFAULT_INDEXTYPE max = 0;
    for (const auto& it : coeffs) {
        unsigned FLENS_DEFAULT_INDEXTYPE idx = u.map()(it.first, dim);
        max = (idx>max) ? idx : max;
    }

    return max;
}


template <typename T, typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintindhash(const IndexSet<Index>& active,
              HTCoefficients<T, Basis>& u,
              const int dim)
{
    unsigned FLENS_DEFAULT_INDEXTYPE max = 0;
    for (const auto& it : active) {
        unsigned FLENS_DEFAULT_INDEXTYPE idx = u.map()(it, dim);
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
                     <flens::Array<unsigned FLENS_DEFAULT_INDEXTYPE> >             IDV;
    typedef typename htucker::DenseVectorList<T>                DVList;

    DVList  list;
    IDV     sizes(tree.dim());

    for (size_type j=1; j<=cp.dim(); ++j) {
        unsigned FLENS_DEFAULT_INDEXTYPE sizeij = 0;
        for (size_type i=1; i<=cp.rank(); ++i) {
            sizeij   = maxintindhash(cp.getCoefficients(i, j), j, tree);
            sizes(j) = (sizeij>sizes(j)) ? sizeij : sizes(j);
        }

        for (size_type i=1; i<=cp.rank(); ++i) {
            DV x(sizes(j));
            for (const auto& it : cp.getCoefficients(i, j)) {
                x(tree.map()(it.first, j)) = it.second;
            }
            list.add(x);
        }
    }

    tree.tree().generateTofElementary(list, cp.rank(), cp.dim());
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(      HTCoefficients<T, Basis>&        tree,
          SepCoefficients<S, T, Index>&    cp,
    const std::vector<IndexSet<Index1D> >& active)
{
    assert((unsigned) tree.dim()==cp.dim());
    assert(active.size()==cp.dim());

    typedef typename SepCoefficients<S, T, Index>::size_type    size_type;
    typedef typename flens::DenseVector<flens::Array<T> >       DV;
    typedef typename htucker::DenseVectorList<T>                DVList;

    DVList  list;

    for (size_type j=1; j<=cp.dim(); ++j) {
        unsigned FLENS_DEFAULT_INDEXTYPE sizej = 0;
        sizej = maxintindhash(active[j-1], j, tree);

        for (size_type i=1; i<=cp.rank(); ++i) {
            DV x(sizej);
            for (const auto& it : active[j-1]) {
                x(tree.map()(it, j)) = cp(i, j)[it];
            }
            list.add(x);
        }
    }

    tree.tree().generateTofElementary(list, cp.rank(), cp.dim());
}


template <typename T, typename Basis>
void
init(HTCoefficients<T, Basis>&              tree,
     const std::vector<IndexSet<Index1D> >& activex,
     const unsigned                         rank,
     const T                                value)
{
    assert(activex.size()==(unsigned) tree.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;
    typedef typename flens::DenseVector<flens::Array<T> >       DV;
    typedef typename htucker::DenseVectorList<T>                DVList;

    DVList  list;

    for (size_type j=1; j<=(unsigned) tree.dim(); ++j) {
        for (size_type i=1; i<=rank; ++i) {
            unsigned FLENS_DEFAULT_INDEXTYPE size
            = maxintindhash(activex[j-1], tree, j);
            DV x(size);
            if (value!=0.) {
                for (auto i=x.firstIndex(); i<=x.lastIndex(); ++i) {
                    x(i) = value;
                }
            }
            list.add(x);
        }
    }

    tree.tree().generateTofElementary(list, rank, tree.dim());
}


template <typename T, typename Basis>
void
rndinit(HTCoefficients<T, Basis>&              tree,
        const std::vector<IndexSet<Index1D> >& activex,
        const unsigned                         rank,
        const T                                scale,
        const long                             seed)
{
    assert(activex.size()==(unsigned) tree.dim());

    SepCoefficients<Lexicographical, T, Index1D> coeffs(rank, activex.size());
    genCoefficientsRnd(coeffs, activex, scale, seed);

    set(tree, coeffs);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const unsigned FLENS_DEFAULT_INDEXTYPE col, const Coefficients<S, T, Index>& coeff)
{
    assert(idx.length()>0);
    assert(coeff.size()>0);
    assert(col>0);
    assert(idx.min()>=1 && idx.max()<=tree.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >  Matrix;
    typedef flens::DenseVector<flens::Array<T> >                      DV;
    using flens::_;

    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            unsigned FLENS_DEFAULT_INDEXTYPE rows = U.numRows();
            unsigned FLENS_DEFAULT_INDEXTYPE cols = U.numCols();

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintindhash(coeff, idx[0],
                                                                 tree);
            if (!rows || !cols) {
                U.resize(max, col);
            } else if (max>rows || col>cols) {
                FLENS_DEFAULT_INDEXTYPE sizer = (max>rows) ? max : rows;
                FLENS_DEFAULT_INDEXTYPE sizec = (col>cols) ? col : cols;
                U.resize(sizer, sizec);
            } else {
                U(_, col) = DV(rows);
            }

            for (const auto& it : coeff) {
                U(tree.map()(it.first, idx[0]), col) = it.second;
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
                auto max = maxintindhash(coeff(i, 1), idx[0], tree);
                rowsset = (max>rowsset) ? max : rowsset;
            }

            if (!rowsnode || !colsnode) {
                U.resize(rowsset, colsset);
            } else if (rowsset>rowsnode || colsset>colsnode) {
                size_type sizer = (rowsset>rowsnode) ? rowsset : rowsnode;
                size_type sizec = (colsset>colsnode) ? colsset : colsnode;
                U.resize(sizer, sizec);
            } else {
                U.fill((T) 0);
            }

            for (size_type i=1; i<=colsset; ++i) {
                for (const auto& it : coeff(i, 1)) {
                    U(tree.map()(it.first, idx[0]), i) = it.second;
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
set_inplace(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
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
                auto max = maxintindhash(coeff(i, 1), idx[0], tree);
                rowsset = (max>rowsset) ? max : rowsset;
            }

            if (!rowsnode || !colsnode) {
                U.resize(rowsset, colsset);
            } else if (rowsset>rowsnode || colsset>colsnode) {
                size_type sizer = (rowsset>rowsnode) ? rowsset : rowsnode;
                size_type sizec = (colsset>colsnode) ? colsset : colsnode;
                Matrix Ucopy    = U;
                U.resize(sizer, sizec);
                U(_(1, rowsnode), _(1, colsnode)) = Ucopy;
            }

            for (size_type i=1; i<=colsset; ++i) {
                for (const auto& it : coeff(i, 1)) {
                    U(tree.map()(it.first, idx[0]), i) = it.second;
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
     const unsigned FLENS_DEFAULT_INDEXTYPE col, const T alpha,
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
            unsigned FLENS_DEFAULT_INDEXTYPE rows = U.numRows();
            unsigned FLENS_DEFAULT_INDEXTYPE cols = U.numCols();

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintindhash(coeff,
                                                                 idx[0],
                                                                 tree);
            if (!rows || !cols) {
                U.resize(max, col);
            } else if (max>rows || col>cols) {
                Matrix copy(U);
                FLENS_DEFAULT_INDEXTYPE sizer = (max>rows) ? max : rows;
                FLENS_DEFAULT_INDEXTYPE sizec = (col>cols) ? col : cols;
                U.resize(sizer, sizec);
                U(_(1,rows),_(1,cols)) = copy;
            }

            for (const auto& it : coeff) {
                U(tree.map()(it.first, idx[0]), col) += alpha*it.second;
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
                auto max = maxintindhash(coeff(i, 1), idx[0], tree);
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
                    U(tree.map()(it.first, idx[0]), i) += alpha*it.second;
                }
            }

            return;
        }
    }

    std::cerr << "error axpy: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned FLENS_DEFAULT_INDEXTYPE col, const T alpha,
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
            unsigned FLENS_DEFAULT_INDEXTYPE rows = U.numRows();
            unsigned FLENS_DEFAULT_INDEXTYPE cols = U.numCols();

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintindhash(coeff,
                                                                 idx[0],
                                                                 tree);
            if (!rows || !cols) {
                U.resize(max, col);
            } else if (max>rows || col>cols) {
                Matrix copy(U);
                FLENS_DEFAULT_INDEXTYPE sizer = (max>rows) ? max : rows;
                FLENS_DEFAULT_INDEXTYPE sizec = (col>cols) ? col : cols;
                U.resize(sizer, sizec);
                U(_(1,rows),_(1,cols)) = copy;
            }

            for (const auto& it : coeff) {
                FLENS_DEFAULT_INDEXTYPE rowi = tree.map()(it.first, idx[0]);
                U(rowi, col) = alpha*U(rowi, col)+it.second;
            }

            return;
        }
    }

    std::cerr << "error xpay: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
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
                auto max = maxintindhash(coeff(i, 1), idx[0], tree);
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
                    FLENS_DEFAULT_INDEXTYPE rowi = tree.map()(it.first, idx[0]);
                    U(rowi, i) = alpha*U(rowi, i)+it.second;
                }
            }

            return;
        }
    }

    std::cerr << "error xpay: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx,
        const unsigned FLENS_DEFAULT_INDEXTYPE col)
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
            for (unsigned FLENS_DEFAULT_INDEXTYPE i=1; i<=(unsigned) U.numRows(); ++i) {
                if (U(i, col)!=(T) 0) {
                    ret[tree.map()(i, idx[0])] = U(i, col);
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
            for (unsigned FLENS_DEFAULT_INDEXTYPE i=1; i<=(unsigned) U.numRows(); ++i) {
                Index1D index = tree.map()(i, idx[0]);
                for (unsigned FLENS_DEFAULT_INDEXTYPE j=1; j<=(unsigned) U.numCols(); ++j) {
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


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
extract(      HTCoefficients<T, Basis>& tree,
        const IndexSet<Index1D>& Lambda,
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
            for (const auto& index : Lambda) {
                auto i = tree.map()(index, idx[0]);
                for (unsigned FLENS_DEFAULT_INDEXTYPE j=1;
                     j<=(unsigned) U.numCols(); ++j) {
                    ret(j, 1)[index] = U(i, j);
                }
            }

            return ret;
        }
    }

    std::cerr << "error extract: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
extract(      HTCoefficients<T, Basis>& tree,
        const IndexSet<Index1D>&        Lambda,
        const unsigned                  j)
{
    assert(j >= 1 && j <= (unsigned) tree.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;

    htucker::DimensionIndex idx(1);
    idx[0] = j;
    for (auto tit=tree.tree().getGeneralTree().end();
              tit>=tree.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            const Matrix& U = tit.getNode()->getContent()
                              ->getUorB();
            if (U.numCols()>1)
                std::cerr << "warning extract: tensor not rank 1\n";

            Coefficients<Lexicographical, T, Index1D>   ret;
            for (const auto& index : Lambda) {
                auto i = tree.map()(index, idx[0]);
                ret[index] = U(i, 1);
            }

            return ret;
        }
    }

    std::cerr << "error extract: dimension not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, SortingCriterion S, typename Index>
SepCoefficients<S, T, Index>
add(const SepCoefficients<S, T, Index>& left,
    const SepCoefficients<S, T, Index>& right)
{
    if (!(left.rank()*left.dim())) {
        return right;
    } else if (!(right.rank()*right.dim())) {
        return left;
    } else {
        assert(left.dim()==right.dim());
    }

    typedef typename SepCoefficients<S, T, Index>::size_type    size_type;

    SepCoefficients<S, T, Index> sum(left.rank()+right.rank(), left.dim());
    for (size_type i=1; i<=sum.rank(); ++i) {
        for (size_type j=1; j<=sum.dim(); ++j) {
            if (i<=left.rank()) {
                sum(i, j) = left(i, j);
            } else {
                sum(i, j) = right(i-left.rank(), j);
            }
        }
    }

    return sum;
}


template <typename T, SortingCriterion S, typename Index>
SepCoefficients<S, T, Index>
operator+(const SepCoefficients<S, T, Index>& left,
          const SepCoefficients<S, T, Index>& right)
{
    return add(left, right);
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalstandard(Sepop<Optype>& A,
             const SepCoefficients<Lexicographical, T, Index1D>& u,
             const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==standard);

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=A.dim(); ++j) {
                // Too many copies here
                TreeCoefficients1D<T> input(hashtablelength,
                                            A(i, j).getTrialBasis().j0);
                TreeCoefficients1D<T> output(hashtablelength,
                                            A(i, j).getTestBasis().j0);
                input = u(k, j);
                Coefficients<Lexicographical, T, Index1D> temp;
                FillWithZeros(A.getrows(j), temp);
                output = temp;


                A(i, j).eval(input, output, "A");
                fromTreeCoefficientsToCoefficients(output, temp);
                prod(1, j) = temp;
            }
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalsimple(Sepop<Optype>& A,
           const SepCoefficients<Lexicographical, T, Index1D>& u,
           const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==simple);

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=prod.dim(); ++j) {
                if (j!=i) {
                    prod(1, j) = u(k, j);
                }
            }

            // Too many copies here
            TreeCoefficients1D<T> input(hashtablelength,
                                        A(i, 1).getTrialBasis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(i, 1).getTestBasis().j0);
            input = u(k, i);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(A.getrows(i), temp);
            output = temp;

            A(i, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            prod(1, i) = temp;

            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evallaplace(Sepop<Optype>& A,
            const SepCoefficients<Lexicographical, T, Index1D>& u,
            const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==laplace);

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=prod.dim(); ++j) {
                if (j!=i) {
                    prod(1, j) = u(k, j);
                }
            }

            // Too many copies here
            TreeCoefficients1D<T> input(hashtablelength,
                                        A(1, 1).getTrialBasis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(1, 1).getTestBasis().j0);
            input = u(k, i);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(A.getrows(i), temp);
            output = temp;

            A(1, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            prod(1, i) = temp;
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepop<Optype>& A,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());

    if (A.type()==standard) {
        return evalstandard(A, u, hashtablelength);
    } else if (A.type()==simple) {
        return evalsimple(A, u, hashtablelength);
    } else {
        return evallaplace(A, u, hashtablelength);
    }
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
operator*(Sepop<Optype>& A,
          const SepCoefficients<Lexicographical, T, Index1D>& u)
{
    assert(A.dim()==u.dim());
    return eval(A, u);
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalstandard(Sepop<Optype>& A,
             const SepCoefficients<Lexicographical, T, Index1D>& u,
             const std::vector<IndexSet<Index1D> >& rows,
             const std::vector<IndexSet<Index1D> >& cols,
             const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==standard);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=A.dim(); ++j) {
                // Too many copies here
                TreeCoefficients1D<T> input(hashtablelength,
                                            A(i, j).getTrialBasis().j0);
                TreeCoefficients1D<T> output(hashtablelength,
                                            A(i, j).getTestBasis().j0);
                Coefficients<Lexicographical, T, Index1D> rest = u(k, j);
                P(cols[j-1], rest);
                input = rest;
                Coefficients<Lexicographical, T, Index1D> temp;
                FillWithZeros(rows[j-1], temp);
                output = temp;

                A(i, j).eval(input, output, "A");
                fromTreeCoefficientsToCoefficients(output, temp);
                prod(1, j) = temp;
            }
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalsimple(Sepop<Optype>& A,
           const SepCoefficients<Lexicographical, T, Index1D>& u,
           const std::vector<IndexSet<Index1D> >& rows,
           const std::vector<IndexSet<Index1D> >& cols,
           const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==simple);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=prod.dim(); ++j) {
                if (j!=i) {
                    prod(1, j) = u(k, j);
                }
            }

            // Too many copies here
            TreeCoefficients1D<T> input(hashtablelength,
                                        A(i, 1).getTrialBasis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(i, 1).getTestBasis().j0);
            Coefficients<Lexicographical, T, Index1D> rest = u(k, i);
            P(cols[i-1], rest);
            input = rest;
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(rows[i-1], temp);
            output = temp;

            A(i, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            prod(1, i) = temp;
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evallaplace(Sepop<Optype>& A,
            const SepCoefficients<Lexicographical, T, Index1D>& u,
            const std::vector<IndexSet<Index1D> >& rows,
            const std::vector<IndexSet<Index1D> >& cols,
            const std::size_t hashtablelength)
{
    assert(A.dim()==u.dim());
    assert(A.type()==laplace);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    SepCoefficients<Lexicographical, T, Index1D> sum;
    for (size_type i=1; i<=A.rank(); ++i) {
        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, A.dim());
            for (size_type j=1; j<=prod.dim(); ++j) {
                if (j!=i) {
                    prod(1, j) = u(k, j);
                }
            }

            // Too many copies here
            TreeCoefficients1D<T> input(hashtablelength,
                                        A(1, 1).getTrialBasis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(1, 1).getTestBasis().j0);
            Coefficients<Lexicographical, T, Index1D> rest = u(k, i);
            P(cols[i-1], rest);
            input = rest;
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(rows[i-1], temp);
            output = temp;

            A(1, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            prod(1, i) = temp;
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepop<Optype>& A,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::vector<IndexSet<Index1D> >& rows,
     const std::vector<IndexSet<Index1D> >& cols,
     const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    if (A.type()==standard) {
        return evalstandard(A, u, rows, cols, hashtablelength);
    } else if (A.type()==simple) {
        return evalsimple(A, u, rows, cols, hashtablelength);
    } else {
        return evallaplace(A, u, rows, cols, hashtablelength);
    }
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalstandard(Sepop<Optype>& A,
             const HTCoefficients<T, Basis>& u,
             const double eps,
             const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==standard);

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        for (size_type j=1; j<=A.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, idx);

            // Too many copies here
            for (size_type k=1; k<=frame.rank(); ++k) {
                TreeCoefficients1D<T> input(hashtablelength,
                                            u.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength,
                                            A(i, j).getTestBasis().j0);
                input = frame(k, 1);
                Coefficients<Lexicographical, T, Index1D> temp;
                FillWithZeros(A.getrows(j), temp);
                output = temp;

                A(i, j).eval(input, output, "A");
                fromTreeCoefficientsToCoefficients(output, temp);
                frame(k, 1) = temp;
            }
            set(prod, idx, frame);
        }

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evalstandard(A): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple(Sepop<Optype>& A,
           const HTCoefficients<T, Basis>& u,
           const double eps,
           const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==simple);

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        htucker::DimensionIndex idx(1);
        idx[0] = i;

        SepCoefficients<Lexicographical, T, Index1D>
        frame = extract(prod, idx);

        // Too many copies here
        for (size_type k=1; k<=frame.rank(); ++k) {
            TreeCoefficients1D<T> input(hashtablelength,
                                        u.basis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(i, 1).getTestBasis().j0);
            input = frame(k, 1);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(A.getrows(i), temp);
            output = temp;

            A(i, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            frame(k, 1) = temp;
        }
        set(prod, idx, frame);

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evalsimple(A): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace(Sepop<Optype>& A,
            const HTCoefficients<T, Basis>& u,
            const double eps,
            const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        htucker::DimensionIndex idx(1);
        idx[0] = i;

        SepCoefficients<Lexicographical, T, Index1D>
        frame = extract(prod, idx);

        // Too many copies here
        for (size_type k=1; k<=frame.rank(); ++k) {
            TreeCoefficients1D<T> input(hashtablelength,
                                        u.basis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(1, 1).getTestBasis().j0);
            input = frame(k, 1);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(A.getrows(i), temp);
            output = temp;

            A(1, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            frame(k, 1) = temp;
        }
        set(prod, idx, frame);

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evallaplace(A): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(Sepop<Optype>& A,
     const HTCoefficients<T, Basis>& u,
     const double eps,
     const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());

    if (A.type()==standard) {
        return evalstandard(A, u, eps, hashtablelength);
    } else if (A.type()==simple) {
        return evalsimple(A, u, eps, hashtablelength);
    } else {
        return evallaplace(A, u, eps, hashtablelength);
    }
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
operator*(Sepop<Optype>& A,
          const HTCoefficients<T, Basis>& u)
{
    assert(A.dim()==(unsigned) u.dim());
    return eval(A, u);
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalstandard(Sepop<Optype>& A,
             const HTCoefficients<T, Basis>& u,
             const std::vector<IndexSet<Index1D> >& rows,
             const std::vector<IndexSet<Index1D> >& cols,
             const double eps,
             const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==standard);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        for (size_type j=1; j<=A.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, cols[j-1], idx);

            // Too many copies here
            for (size_type k=1; k<=frame.rank(); ++k) {
                TreeCoefficients1D<T> input(hashtablelength,
                                            u.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength,
                                            A(i, j).getTestBasis().j0);
                P(cols[j-1], frame(k, 1));
                input = frame(k, 1);
                Coefficients<Lexicographical, T, Index1D> temp;
                FillWithZeros(rows[j-1], temp);
                output = temp;

                A(i, j).eval(input, output, "A");
                fromTreeCoefficientsToCoefficients(output, temp);
                frame(k, 1) = temp;
            }
            set(prod, idx, frame);
        }

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evalstandard(A, rows, cols): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple_old(Sepop<Optype>& A,
               const HTCoefficients<T, Basis>& u,
               const std::vector<IndexSet<Index1D> >& rows,
               const std::vector<IndexSet<Index1D> >& cols,
               const double eps,
               const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==simple);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        htucker::DimensionIndex idx(1);
        idx[0] = i;

        SepCoefficients<Lexicographical, T, Index1D>
        frame = extract(prod, cols[i-1], idx);

        // Too many copies here
        for (size_type k=1; k<=frame.rank(); ++k) {
            TreeCoefficients1D<T> input(hashtablelength,
                                        u.basis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(i, 1).getTestBasis().j0);
            P(cols[i-1], frame(k, 1));
            input = frame(k, 1);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(rows[i-1], temp);
            output = temp;

            A(i, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            frame(k, 1) = temp;
        }
        set(prod, idx, frame);

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evalsimple(A, rows, cols): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace_old(Sepop<Optype>& A,
                const HTCoefficients<T, Basis>& u,
                const std::vector<IndexSet<Index1D> >& rows,
                const std::vector<IndexSet<Index1D> >& cols,
                const double eps,
                const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef typename Sepop<Optype>::size_type   size_type;

    HTCoefficients<T, Basis> sum;
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (size_type i=1; i<=A.rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        htucker::DimensionIndex idx(1);
        idx[0] = i;

        SepCoefficients<Lexicographical, T, Index1D>
        frame = extract(prod, cols[i-1], idx);

        // Too many copies here
        for (size_type k=1; k<=frame.rank(); ++k) {
            TreeCoefficients1D<T> input(hashtablelength,
                                        u.basis().j0);
            TreeCoefficients1D<T> output(hashtablelength,
                                        A(1, 1).getTestBasis().j0);
            P(cols[i-1], frame(k, 1));
            input = frame(k, 1);
            Coefficients<Lexicographical, T, Index1D> temp;
            FillWithZeros(rows[i-1], temp);
            output = temp;

            A(1, 1).eval(input, output, "A");
            fromTreeCoefficientsToCoefficients(output, temp);
            frame(k, 1) = temp;
        }
        set(prod, idx, frame);

        if (i==1) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            // Too many resizes here
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "evallaplace(A, rows, cols): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple(Sepop<Optype>& A,
                 HTCoefficients<T, Basis>& u,
           const std::vector<IndexSet<Index1D> >& rows,
           const std::vector<IndexSet<Index1D> >& cols,
           const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==simple);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    HTCoefficients<T, Basis> Au(u.dim(), u.basis(), u.map());
    Au.tree().set_tree(u.tree());

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titu = u.tree().getGeneralTree().end();

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titAu = Au.tree().getGeneralTree().end();
    for (;titu>=u.tree().getGeneralTree().begin(); titu--, titAu--) {
        auto nodeu  = titu.getNode();
        auto nodeAu = titAu.getNode();

        /* Apply operator to leafs */
        if (nodeu->isLeaf()) {
            htucker::DimensionIndex idx = nodeu->getContent()->getIndex();
            Matrix& Uu  = nodeu->getContent()->getUorB();
            Matrix& UAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            auto rowsU  = Uu.numRows();
            auto colsU  = Uu.numCols();
            auto max    = maxintindhash(rows[idx[0]-1], idx[0], u);
            assert(max>=(unsigned) rowsU);
            UAu.resize(max, 2*colsU);
            UAu(_(1, rowsU), _(1, colsU)) = Uu;

            /* Apply operator to columns */
            for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
                TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
                typedef typename TreeCoefficients1D<T>::val_type val_type;

                /* Set input */
                for (auto& lambda : cols[idx[0]-1]) {
                    auto i = u.map()(lambda, idx[0]);
                    assert(i<=(unsigned) rowsU);

                    auto j     = lambda.j;
                    auto _k     = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        input.bylevel[j-1-input.offset].
                              map.insert(val_type(_k, Uu(i, k)));
                    } else {
                        input.bylevel[j-input.offset].
                              map.insert(val_type(_k, Uu(i, k)));
                    }
                }

                /* Set output */
                for (auto& lambda : rows[idx[0]-1]) {
                    #ifndef NDEBUG
                        auto i = u.map()(lambda, idx[0]);
                        assert(i<=max);
                    #endif

                    auto j     = lambda.j;
                    auto _k    = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        output.bylevel[j-1-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    } else {
                        output.bylevel[j-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    }
                }

                /* Apply A */
                A(idx[0], 1).eval(input, output, "A");

                /* Save result */
                for (typename CoefficientsByLevel<T>::const_it it=
                     output.bylevel[0].map.begin();
                     it!=output.bylevel[0].map.end(); ++it) {
                        Index1D lambda(output.offset+1, (*it).first,XBSpline);
                        auto i = u.map()(lambda, idx[0]);
                        UAu(i, colsU+k) = (*it).second;
                }
                for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                    if (output.bylevel[i].map.size()==0) break;
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[i].map.begin();
                         it!=output.bylevel[i].map.end(); ++it) {
                            Index1D lambda(output.offset+i,
                                           (*it).first,XWavelet);
                            auto j = u.map()(lambda, idx[0]);
                            UAu(j, colsU+k) = (*it).second;
                    }
                }
            }
        } else if (nodeu->isInner()) {

            /* Set (sparse) transfer tensors */
            auto nodenumel   = nodeu->getContent()->getNumRows();
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(4*nodenumel*nodercnumel, 2*nodelcnumel);

            /* First column B of HT(A) */
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                BAu(_(2*nodercnumel*(i-1)+1, (2*i-1)*nodercnumel)
                   ,_(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Second column B of HT(A) */
            FLENS_DEFAULT_INDEXTYPE offset1 = (2*nodenumel+1)*nodercnumel;
            FLENS_DEFAULT_INDEXTYPE offset2 = 2*nodenumel*nodercnumel;
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                /* Left block column */
                BAu(_(offset1+2*nodercnumel*(i-1)+1,
                      offset1+(2*i-1)*nodercnumel), _(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);

                /* Right block column */
                BAu(_(offset2+2*nodercnumel*(i-1)+1,
                      offset2+(2*i-1)*nodercnumel),
                    _(nodelcnumel+1, 2*nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Set meta data */
            nodeAu->getContent()->setNumRows(2*nodenumel);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        } else {

            /* Set root */
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(2*nodercnumel, 2*nodelcnumel);

            /* Left block column */
            BAu(_(nodercnumel+1, 2*nodercnumel), _(1, nodelcnumel)) = Bu;

            /* Right block column */
            BAu(_(1, nodercnumel), _(nodelcnumel+1, 2*nodelcnumel)) = Bu;

            /* Set meta data */
            nodeAu->getContent()->setNumRows(1);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        }
    }

    return Au;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace(Sepop<Optype>& A,
                  HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index1D> >& rows,
            const std::vector<IndexSet<Index1D> >& cols,
            const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    HTCoefficients<T, Basis> Au(u.dim(), u.basis(), u.map());
    Au.tree().set_tree(u.tree());

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titu = u.tree().getGeneralTree().end();

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titAu = Au.tree().getGeneralTree().end();
    for (;titu>=u.tree().getGeneralTree().begin(); titu--, titAu--) {
        auto nodeu  = titu.getNode();
        auto nodeAu = titAu.getNode();

        /* Apply operator to leafs */
        if (nodeu->isLeaf()) {
            htucker::DimensionIndex idx = nodeu->getContent()->getIndex();
            Matrix& Uu  = nodeu->getContent()->getUorB();
            Matrix& UAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            auto rowsU  = Uu.numRows();
            auto colsU  = Uu.numCols();
            auto max    = maxintindhash(rows[idx[0]-1], idx[0], u);
            assert(max>=(unsigned) rowsU);
            UAu.resize(max, 2*colsU);
            UAu(_(1, rowsU), _(1, colsU))         = Uu;

            /* Apply operator to columns */
            for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
                TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
                typedef typename TreeCoefficients1D<T>::val_type val_type;

                /* Set input */
                for (auto& lambda : cols[idx[0]-1]) {
                    auto i = u.map()(lambda, idx[0]);
                    assert(i<=(unsigned) rowsU);

                    auto j     = lambda.j;
                    auto _k     = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        input.bylevel[j-1-input.offset].
                              map.insert(val_type(_k, Uu(i, k)));
                    } else {
                        input.bylevel[j-input.offset].
                              map.insert(val_type(_k, Uu(i, k)));
                    }
                }

                /* Set output */
                for (auto& lambda : rows[idx[0]-1]) {
                    #ifndef NDEBUG
                        auto i = u.map()(lambda, idx[0]);
                        assert(i<=max);
                    #endif

                    auto j     = lambda.j;
                    auto _k     = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        output.bylevel[j-1-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    } else {
                        output.bylevel[j-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    }
                }

                /* Apply A */
                A(1, 1).eval(input, output, "A");

                /* Save result */
                for (typename CoefficientsByLevel<T>::const_it it=
                     output.bylevel[0].map.begin();
                     it!=output.bylevel[0].map.end(); ++it) {
                        Index1D lambda(output.offset+1, (*it).first,XBSpline);
                        auto i = u.map()(lambda, idx[0]);
                        UAu(i, colsU+k) = (*it).second;
                }
                for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                    if (output.bylevel[i].map.size()==0) break;
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[i].map.begin();
                         it!=output.bylevel[i].map.end(); ++it) {
                            Index1D lambda(output.offset+i,
                                           (*it).first,XWavelet);
                            auto j = u.map()(lambda, idx[0]);
                            UAu(j, colsU+k) = (*it).second;
                    }
                }
            }

            nodeAu->getContent()->setNumRows(2*colsU);
        } else if (nodeu->isInner()) {

            /* Set (sparse) transfer tensors */
            auto nodenumel   = nodeu->getContent()->getNumRows();
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(4*nodenumel*nodercnumel, 2*nodelcnumel);

            /* First column B of HT(A) */
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                BAu(_(2*nodercnumel*(i-1)+1, (2*i-1)*nodercnumel)
                   ,_(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Second column B of HT(A) */
            FLENS_DEFAULT_INDEXTYPE offset1 = (2*nodenumel+1)*nodercnumel;
            FLENS_DEFAULT_INDEXTYPE offset2 = 2*nodenumel*nodercnumel;
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                /* Left block column */
                BAu(_(offset1+2*nodercnumel*(i-1)+1,
                      offset1+(2*i-1)*nodercnumel), _(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);

                /* Right block column */
                BAu(_(offset2+2*nodercnumel*(i-1)+1,
                      offset2+(2*i-1)*nodercnumel),
                    _(nodelcnumel+1, 2*nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Set meta data */
            nodeAu->getContent()->setNumRows(2*nodenumel);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        } else {

            /* Set root */
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(2*nodercnumel, 2*nodelcnumel);

            /* Left block column */
            BAu(_(nodercnumel+1, 2*nodercnumel), _(1, nodelcnumel)) = Bu;

            /* Right block column */
            BAu(_(1, nodercnumel), _(nodelcnumel+1, 2*nodelcnumel)) = Bu;

            /* Set meta data */
            nodeAu->getContent()->setNumRows(1);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        }
    }

    return Au;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalLaplaceD_1Dim(      Sepop<Optype>& A,
                        HTCoefficients<T, Basis>& u,
                  const std::vector<IndexSet<Index1D> >& rows,
                  const std::vector<IndexSet<Index1D> >& cols,
                  const unsigned j,
                  const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);
    assert(j>=1 && j<=A.dim());
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    HTCoefficients<T, Basis> Au(u.dim(), u.basis(), u.map());
    Au.tree().set_tree(u.tree());

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titu = u.tree().getGeneralTree().end();

    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titAu = Au.tree().getGeneralTree().end();
    for (;titu>=u.tree().getGeneralTree().begin(); titu--, titAu--) {
        auto nodeu  = titu.getNode();
        auto nodeAu = titAu.getNode();

        /* Apply operator to leafs */
        if (nodeu->isLeaf()) {
            htucker::DimensionIndex idx = nodeu->getContent()->getIndex();
            Matrix& Uu  = nodeu->getContent()->getUorB();
            Matrix& UAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            auto rowsU  = Uu.numRows();
            auto colsU  = Uu.numCols();
            if ((unsigned) idx[0] == j) {
                UAu.resize(1, 2*colsU);
                nodeAu->getContent()->setNumRows(2*colsU);
            } else {
                auto max    = maxintindhash(rows[idx[0]-1], idx[0], u);
                assert(max>=(unsigned) rowsU);
                UAu.resize(max, 2*colsU);
                UAu(_(1, rowsU), _(1, colsU))         = Uu;

                /* Apply operator to columns */
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
                    TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
                    TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
                    typedef typename TreeCoefficients1D<T>::val_type val_type;

                    /* Set input */
                    for (auto& lambda : cols[idx[0]-1]) {
                        auto i = u.map()(lambda, idx[0]);
                        assert(i<=(unsigned) rowsU);

                        auto j     = lambda.j;
                        auto _k     = lambda.k;
                        auto xtype = lambda.xtype;
                        if (xtype==XBSpline) {
                            input.bylevel[j-1-input.offset].
                                  map.insert(val_type(_k, Uu(i, k)));
                        } else {
                            input.bylevel[j-input.offset].
                                  map.insert(val_type(_k, Uu(i, k)));
                        }
                    }

                    /* Set output */
                    for (auto& lambda : rows[idx[0]-1]) {
                        #ifndef NDEBUG
                            auto i = u.map()(lambda, idx[0]);
                            assert(i<=max);
                        #endif

                        auto j     = lambda.j;
                        auto _k     = lambda.k;
                        auto xtype = lambda.xtype;
                        if (xtype==XBSpline) {
                            output.bylevel[j-1-output.offset].
                                  map.insert(val_type(_k, (T) 0));
                        } else {
                            output.bylevel[j-output.offset].
                                  map.insert(val_type(_k, (T) 0));
                        }
                    }

                    /* Apply A */
                    A(1, 1).eval(input, output, "A");

                    /* Save result */
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[0].map.begin();
                         it!=output.bylevel[0].map.end(); ++it) {
                            Index1D lambda(output.offset+1, (*it).first,XBSpline);
                            auto i = u.map()(lambda, idx[0]);
                            UAu(i, colsU+k) = (*it).second;
                    }
                    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                        if (output.bylevel[i].map.size()==0) break;
                        for (typename CoefficientsByLevel<T>::const_it it=
                             output.bylevel[i].map.begin();
                             it!=output.bylevel[i].map.end(); ++it) {
                                Index1D lambda(output.offset+i,
                                               (*it).first,XWavelet);
                                auto j = u.map()(lambda, idx[0]);
                                UAu(j, colsU+k) = (*it).second;
                        }
                    }
                }

                nodeAu->getContent()->setNumRows(2*colsU);
            }
        } else if (nodeu->isInner()) {

            /* Set (sparse) transfer tensors */
            auto nodenumel   = nodeu->getContent()->getNumRows();
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(4*nodenumel*nodercnumel, 2*nodelcnumel);

            /* First column B of HT(A) */
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                BAu(_(2*nodercnumel*(i-1)+1, (2*i-1)*nodercnumel)
                   ,_(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Second column B of HT(A) */
            FLENS_DEFAULT_INDEXTYPE offset1 = (2*nodenumel+1)*nodercnumel;
            FLENS_DEFAULT_INDEXTYPE offset2 = 2*nodenumel*nodercnumel;
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=nodenumel; ++i) {
                /* Left block column */
                BAu(_(offset1+2*nodercnumel*(i-1)+1,
                      offset1+(2*i-1)*nodercnumel), _(1, nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);

                /* Right block column */
                BAu(_(offset2+2*nodercnumel*(i-1)+1,
                      offset2+(2*i-1)*nodercnumel),
                    _(nodelcnumel+1, 2*nodelcnumel)) =
                Bu(_(nodercnumel*(i-1)+1, nodercnumel*i), _);
            }

            /* Set meta data */
            nodeAu->getContent()->setNumRows(2*nodenumel);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        } else {

            /* Set root */
            auto nodelcnumel = nodeu->getContent()->getLeftChildNumRows();
            auto nodercnumel = nodeu->getContent()->getRightChildNumRows();
            Matrix& Bu  = nodeu->getContent()->getUorB();
            Matrix& BAu = const_cast<Matrix&>(nodeAu->getContent()->getUorB());
            BAu.resize(2*nodercnumel, 2*nodelcnumel);

            /* Left block column */
            BAu(_(nodercnumel+1, 2*nodercnumel), _(1, nodelcnumel)) = Bu;

            /* Right block column */
            BAu(_(1, nodercnumel), _(nodelcnumel+1, 2*nodelcnumel)) = Bu;

            /* Set meta data */
            nodeAu->getContent()->setNumRows(1);
            nodeAu->getContent()->setLeftChildNumRows(2*nodelcnumel);
            nodeAu->getContent()->setRightChildNumRows(2*nodercnumel);
        }
    }

    return Au;
}

template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
evalLaplace1Dim(      Sepop<Optype>& A,
                      HTCoefficients<T, Basis>& u,
                const IndexSet<Index1D>& rows,
                const IndexSet<Index1D>& cols,
                const unsigned j,
                const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);
    assert(j>=1 && j<=A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    Matrix UAu;
    htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
    titu = u.tree().getGeneralTree().end();
    for (;titu>=u.tree().getGeneralTree().begin(); titu--) {
        auto nodeu  = titu.getNode();

        /* Apply operator to leafs */
        if (nodeu->isLeaf()) {
            htucker::DimensionIndex idx = nodeu->getContent()->getIndex();
            if ((unsigned) idx[0] == j) {
                Matrix& Uu  = nodeu->getContent()->getUorB();
                auto rowsU  = Uu.numRows();
                auto colsU  = Uu.numCols();

                auto max    = maxintindhash(rows, idx[0], u);
                assert(max>=(unsigned) rowsU);
                UAu.resize(max, 2*colsU);
                UAu(_(1, rowsU), _(1, colsU))         = Uu;

                /* Apply operator to columns */
                for (FLENS_DEFAULT_INDEXTYPE k=1;
                     k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
                    TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
                    TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
                    typedef typename TreeCoefficients1D<T>::val_type val_type;

                    /* Set input */
                    for (auto& lambda : cols) {
                        auto i = u.map()(lambda, idx[0]);
                        assert(i<=(unsigned) rowsU);

                        auto j     = lambda.j;
                        auto _k     = lambda.k;
                        auto xtype = lambda.xtype;
                        if (xtype==XBSpline) {
                            input.bylevel[j-1-input.offset].
                                  map.insert(val_type(_k, Uu(i, k)));
                        } else {
                            input.bylevel[j-input.offset].
                                  map.insert(val_type(_k, Uu(i, k)));
                        }
                    }

                    /* Set output */
                    for (auto& lambda : rows) {
                        #ifndef NDEBUG
                            auto i = u.map()(lambda, idx[0]);
                            assert(i<=max);
                        #endif

                        auto j     = lambda.j;
                        auto _k     = lambda.k;
                        auto xtype = lambda.xtype;
                        if (xtype==XBSpline) {
                            output.bylevel[j-1-output.offset].
                                  map.insert(val_type(_k, (T) 0));
                        } else {
                            output.bylevel[j-output.offset].
                                  map.insert(val_type(_k, (T) 0));
                        }
                    }

                    /* Apply A */
                    A(1, 1).eval(input, output, "A");

                    /* Save result */
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[0].map.begin();
                         it!=output.bylevel[0].map.end(); ++it) {
                            Index1D lambda(output.offset+1, (*it).first,XBSpline);
                            auto i = u.map()(lambda, idx[0]);
                            UAu(i, colsU+k) = (*it).second;
                    }
                    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                        if (output.bylevel[i].map.size()==0) break;
                        for (typename CoefficientsByLevel<T>::const_it it=
                             output.bylevel[i].map.begin();
                             it!=output.bylevel[i].map.end(); ++it) {
                                Index1D lambda(output.offset+i,
                                               (*it).first,XWavelet);
                                auto j = u.map()(lambda, idx[0]);
                                UAu(j, colsU+k) = (*it).second;
                        }
                    }
                }

                return UAu;
            }
        }
    }

    return UAu;
}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
evallaplace(      Sepop<Optype>&            A,
            const flens::GeMatrix<
                  flens::FullStorage<T,
                  cxxblas::ColMajor> >&     U,
                  HTCoefficients<T, Basis>& u,
            const unsigned                  j,
            const IndexSet<Index1D>&        rows,
            const IndexSet<Index1D>&        cols,
            const std::size_t               hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());
    assert(A.type()==laplace);

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    Matrix Au;

    auto rowsU  = U.numRows();
    auto colsU  = U.numCols();
    auto max    = maxintindhash(rows, j, u);
    assert(max>=(unsigned) rowsU);
    Au.resize(max, 2*colsU);
    Au(_(1, rowsU), _(1, colsU))         = U;

    /* Apply operator to columns */
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
        TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
        TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
        typedef typename TreeCoefficients1D<T>::val_type val_type;

        /* Set input */
        for (auto& lambda : cols) {
            auto i = u.map()(lambda, j);
            assert(i<=(unsigned) rowsU);

            auto _j    = lambda.j;
            auto _k    = lambda.k;
            auto xtype = lambda.xtype;
            if (xtype==XBSpline) {
                input.bylevel[_j-1-input.offset].
                      map.insert(val_type(_k, U(i, k)));
            } else {
                input.bylevel[_j-input.offset].
                      map.insert(val_type(_k, U(i, k)));
            }
        }

        /* Set output */
        for (auto& lambda : rows) {
            #ifndef NDEBUG
                auto i = u.map()(lambda, j);
                assert(i<=max);
            #endif

            auto _j     = lambda.j;
            auto _k     = lambda.k;
            auto xtype  = lambda.xtype;
            if (xtype==XBSpline) {
                output.bylevel[_j-1-output.offset].
                      map.insert(val_type(_k, (T) 0));
            } else {
                output.bylevel[_j-output.offset].
                      map.insert(val_type(_k, (T) 0));
            }
        }

        /* Apply A */
        A(1, 1).eval(input, output, "A");

        /* Save result */
        for (typename CoefficientsByLevel<T>::const_it it=
                output.bylevel[0].map.begin();
                it!=output.bylevel[0].map.end(); ++it) {
            Index1D lambda(output.offset+1, (*it).first,XBSpline);
            auto i          = u.map()(lambda, j);
            Au(i, colsU+k)  = (*it).second;
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
            if (output.bylevel[i].map.size()==0) break;
            for (typename CoefficientsByLevel<T>::const_it it=
                 output.bylevel[i].map.begin();
                 it!=output.bylevel[i].map.end(); ++it) {
                Index1D lambda(output.offset+i, (*it).first,XWavelet);
                auto _j   = u.map()(lambda, j);
                Au(_j, colsU+k) = (*it).second;
            }
        }
    }

    return Au;
}


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
evallaplace(      Sepop<Optype>&            A,
                  Prec&                     P,
            const flens::GeMatrix<
                  flens::FullStorage<T,
                  cxxblas::ColMajor> >&     U,
                  HTCoefficients<T, Basis>& u,
            const unsigned                  j,
            const IndexSet<Index1D>&        rows,
            const IndexSet<Index1D>&        cols,
            const std::size_t               hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());
    assert(A.type()==laplace);

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    Matrix Au;

    auto rowsU  = U.numRows();
    auto colsU  = U.numCols();
    auto max    = maxintindhash(rows, j, u);
    assert(max>=(unsigned) rowsU);
    Au.resize(max, colsU);

    /* Apply operator to columns */
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
        TreeCoefficients1D<T> input(hashtablelength, u.basis().j0);
        TreeCoefficients1D<T> output(hashtablelength, u.basis().j0);
        typedef typename TreeCoefficients1D<T>::val_type val_type;

        /* Set input */
        for (auto& lambda : cols) {
            auto i = u.map()(lambda, j);
            assert(i<=(unsigned) rowsU);

            auto _j    = lambda.j;
            auto _k    = lambda.k;
            auto xtype = lambda.xtype;

            if (xtype==XBSpline) {
                input.bylevel[_j-1-input.offset].
                      map.insert(val_type(_k, P(lambda)*U(i, k)));
            } else {
                input.bylevel[_j-input.offset].
                      map.insert(val_type(_k,P(lambda)*U(i, k)));
            }
        }

        /* Set output */
        for (auto& lambda : rows) {
            #ifndef NDEBUG
                auto i = u.map()(lambda, j);
                assert(i<=max);
            #endif

            auto _j     = lambda.j;
            auto _k     = lambda.k;
            auto xtype  = lambda.xtype;
            if (xtype==XBSpline) {
                output.bylevel[_j-1-output.offset].
                      map.insert(val_type(_k, (T) 0));
            } else {
                output.bylevel[_j-output.offset].
                      map.insert(val_type(_k, (T) 0));
            }
        }

        /* Apply A */
        A(1, 1).eval(input, output, "A");

        /* Save result */
        for (typename CoefficientsByLevel<T>::const_it it=
                output.bylevel[0].map.begin();
                it!=output.bylevel[0].map.end(); ++it) {
            Index1D lambda(output.offset+1, (*it).first,XBSpline);
            auto i    = u.map()(lambda, j);
            Au(i, k)  = P(lambda)*(*it).second;
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
            if (output.bylevel[i].map.size()==0) break;
            for (typename CoefficientsByLevel<T>::const_it it=
                 output.bylevel[i].map.begin();
                 it!=output.bylevel[i].map.end(); ++it) {
                Index1D lambda(output.offset+i, (*it).first,XWavelet);
                auto _j   = u.map()(lambda, j);
                Au(_j, k) = P(lambda)*(*it).second;
            }
        }
    }

    return Au;
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(Sepop<Optype>& A,
           HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >& rows,
     const std::vector<IndexSet<Index1D> >& cols,
     const double eps,
     const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    if (A.type()==standard) {
        return evalstandard(A, u, rows, cols, eps, hashtablelength);
    } else if (A.type()==simple) {
        return evalsimple(A, u, rows, cols, hashtablelength);
    } else {
        return evallaplace(A, u, rows, cols, hashtablelength);
    }
}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >
eval1Dim(      Sepop<Optype>& A,
               HTCoefficients<T, Basis>& u,
         const IndexSet<Index1D>& rows,
         const IndexSet<Index1D>& cols,
         const unsigned j,
         const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(j>=1 && j<=A.dim());

    if (A.type()==laplace) {
        return evalLaplace1Dim(A, u, rows, cols, j, hashtablelength);
    } else {
        std::cerr << "eval1Dim: Not implemented for operator type\n";
        exit(1);
    }
}


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalD_1Dim(      Sepop<Optype>& A,
                 HTCoefficients<T, Basis>& u,
           const std::vector<IndexSet<Index1D> >& rows,
           const std::vector<IndexSet<Index1D> >& cols,
           const unsigned j,
           const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(rows.size()==cols.size());
    assert(rows.size()==A.dim());
    assert(j>=1 && j<=A.dim());

    if (A.type()==laplace) {
        return evalLaplaceD_1Dim(A, u, rows, cols, j, hashtablelength);
    } else {
        std::cerr << "evalD_1Dim: Not implemented for operator type\n";
        exit(1);
    }

}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepop<Optype>&              A,
     const flens::GeMatrix<
           flens::FullStorage<T,
           cxxblas::ColMajor> >&       U,
           HTCoefficients<T, Basis>&   u,
     const unsigned                    j,
     const IndexSet<Index1D>&          rows,
     const IndexSet<Index1D>&          cols,
     const double                      eps,
     const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());
    if (A.type()==laplace) {
        (void) eps;
        return evallaplace(A, U, u, j, rows, cols, hashtablelength);
    } else {
        std::cerr << "Error eval: not implemented for operator type\n";
        exit(1);
    }
}


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepop<Optype>&              A,
           Prec&                       P,
     const flens::GeMatrix<
           flens::FullStorage<T,
           cxxblas::ColMajor> >&       U,
           HTCoefficients<T, Basis>&   u,
     const unsigned                    j,
     const IndexSet<Index1D>&          rows,
     const IndexSet<Index1D>&          cols,
     const double                      eps,
     const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());
    if (A.type()==laplace) {
        (void) eps;
        return evallaplace(A, P, U, u, j, rows, cols, hashtablelength);
    } else {
        std::cerr << "Error eval: not implemented for operator type\n";
        exit(1);
    }
}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps,
        const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > AUj, AjUj;
    AUj = eval(A, Uj, u, j, rows, cols, eps, hashtablelength);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., AUj, Pj, 0., AjUj);

    return AjUj;
}


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval_laplace(      Sepop<Optype>&              A,
                      Prec&                       P,
                const flens::GeMatrix<
                      flens::FullStorage<T,
                      cxxblas::ColMajor> >&       Uj,
                const flens::GeMatrix<
                      flens::FullStorage<T,
                      cxxblas::ColMajor> >&       Pj,
                      HTCoefficients<T, Basis>&   u,
                const unsigned                    j,
                const IndexSet<Index1D>&          rows,
                const IndexSet<Index1D>&          cols,
                const double                      eps,
                const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());
    assert(A.type()==laplace);
    using flens::_;

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > AjUj, VG1, VG2;
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., Uj,
                    Pj(_, _(1, Uj.numCols())), 0., VG1);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., Uj,
                    Pj(_, _(Uj.numCols()+1, Pj.numCols())), 0., VG2);
    AjUj  = eval(A, P, VG2, u, j, rows, cols, eps, hashtablelength);
    VG1   = precsq(P, VG1, u, j, cols);
    AjUj += VG1;

    return AjUj;
}


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
              Prec&                       P,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps,
        const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > AjUj;
    if (A.type()==laplace) {
        AjUj = redeval_laplace(A, P, Uj, Pj, u, j, rows, cols, eps,
                               hashtablelength);
    } else {
        std::cerr << "redeval: not implemented for operator type\n";
    }

    return AjUj;
}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&                    A,
              Sepdiagscal<Basis>&               S,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&             Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&             Pj,
              HTCoefficients<T, Basis>&         u,
        const unsigned                          j,
        const std::vector<IndexSet<Index1D> >&  rows,
        const std::vector<IndexSet<Index1D> >&  cols,
        const double                            eps,
        const std::size_t                       hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > SUj, ASUj,
                                                               SASUj, AjUj;
    SUj   = eval(S, Uj, u, j, cols);
    ASUj  = eval(A, SUj, u, j, rows[j-1], cols[j-1], eps, hashtablelength);
    SASUj = eval(S, ASUj, u, j, rows);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., SASUj, Pj, 0., AjUj);

    return AjUj;
}


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
              Sepdiagscal<Basis>&         S,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps,
        const std::size_t                 hashtablelength)
{
    assert(A.dim()==(unsigned)u.dim());
    assert(j>=1 && j<=A.dim());

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > SUj, ASUj,
                                                               SASUj, AjUj;
    SUj   = fixeval(S, Uj, u, j, cols);
    ASUj  = eval(A, SUj, u, j, rows, cols, eps, hashtablelength);
    SASUj = fixeval(S, ASUj, u, j, rows);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., SASUj, Pj, 0., AjUj);

    return AjUj;
}


FLENS_DEFAULT_INDEXTYPE
maxlevel(const std::vector<IndexSet<Index1D> >& Lambda)
{
    typedef std::vector<IndexSet<Index1D> >::size_type size_type;

    FLENS_DEFAULT_INDEXTYPE jmax = 0;

    for (size_type i=0; i<Lambda.size(); ++i) {
        for (const auto& lambda : Lambda[i]) {
             FLENS_DEFAULT_INDEXTYPE j = lambda.j;
             if (lambda.xtype==XWavelet) ++j;
             jmax = MAX(jmax, j);
        }
    }

    return jmax;
}


flens::DenseVector<
flens::Array<FLENS_DEFAULT_INDEXTYPE> >
maxlevels(const std::vector<IndexSet<Index1D> >& Lambda)
{
    typedef std::vector<IndexSet<Index1D> >::size_type size_type;

    flens::DenseVector<
    flens::Array<FLENS_DEFAULT_INDEXTYPE> > jmaxs(Lambda.size());

    for (size_type i=0; i<Lambda.size(); ++i) {
        FLENS_DEFAULT_INDEXTYPE jmax=0;
        for (const auto& lambda : Lambda[i]) {
             FLENS_DEFAULT_INDEXTYPE j = lambda.j;
             if (lambda.xtype==XWavelet) ++j;
             jmax = MAX(jmax, j);
        }
        jmaxs(i+1) = jmax;
    }

    return jmaxs;
}


FLENS_DEFAULT_INDEXTYPE
maxlevel(const IndexSet<Index1D>& Lambda)
{

    FLENS_DEFAULT_INDEXTYPE jmax = 0;
    for (const auto& lambda : Lambda) {
         FLENS_DEFAULT_INDEXTYPE j = lambda.j;
         if (lambda.xtype==XWavelet) ++j;
         jmax = MAX(jmax, j);
    }

    return jmax;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
compOmegamin2(const Basis& basis,
              const typename Sepdiagscal<Basis>::size_type d,
              const typename Sepdiagscal<Basis>::T order)
{
    assert(d>0);
    return d*std::pow(2., 2.*order*basis.j0);
}


template <typename T>
T
compOmegamax2(const std::vector<IndexSet<Index1D> >& Lambda, const T order)
{
    auto d = Lambda.size();

    FLENS_DEFAULT_INDEXTYPE jmax = maxlevel(Lambda);

    return d*std::pow(2., 2.*order*jmax);
}


template <typename T>
T
compOmegamax4(const std::vector<IndexSet<Index1D> >& Lambda, const T order)
{
    auto d = Lambda.size();

    FLENS_DEFAULT_INDEXTYPE jmax = maxlevel(Lambda);

    return d*std::pow(2., 4.*order*jmax);
}


template <typename Basis>
typename Basis::T
compUnDistFac(const Basis& basis,
              const std::vector<IndexSet<Index1D> >& Lambda,
              const typename Basis::T order)
{
    auto d   = Lambda.size();
    FLENS_DEFAULT_INDEXTYPE j0   = basis.j0;
    FLENS_DEFAULT_INDEXTYPE jmax = maxlevel(Lambda);

    typename Basis::T factor = std::pow(2., (order+1)*j0-jmax-1.);
    for (FLENS_DEFAULT_INDEXTYPE j=j0+1; j<=jmax; ++j) {
        factor += std::pow(2., (order+1.)*j-jmax-2.);
    }
    factor *= std::sqrt(d);

    return 1./factor;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
compOmegamin4(const Basis& basis,
              const typename Sepdiagscal<Basis>::size_type d,
              const typename Sepdiagscal<Basis>::T order)
{
    assert(d>0);
    return d*std::pow(2., 4.*order*basis.j0);
}


template <typename T, typename Basis>
T
compIndexscale(const SepCoefficients<Lexicographical, T, Index1D>& u,
               const Basis& basis, const T order)
{
    typedef typename SepCoefficients<Lexicographical, T, Index1D>::size_type
                                                                   size_type;

    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >    jmax(u.dim());
    for (size_type j=1; j<=u.dim(); ++j) {
        for (size_type k=1; k<=u.rank(); ++k) {
            for (const auto& it : u(k, j)) {
                FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                if (it.first.xtype==XWavelet) ++level;
                jmax(j) = MAX(level, jmax(j));
            }
        }
    }

    T sum = 0;
    for (size_type i=1; i<=u.dim(); ++i) {
        sum += std::pow(2., 2.*order*jmax(i));
    }

    return sum/compOmegamin2(basis, u.dim(), order);
}


template <typename T, typename Basis>
T
compIndexscale(const HTCoefficients<T, Basis>& u, const T order)
{
    typedef typename SepCoefficients<Lexicographical, T, Index1D>::size_type
                     size_type;

    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >    jmax(u.dim());
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=u.dim(); ++j) {
        htucker::DimensionIndex idx(1);
        idx[0] = j;
        SepCoefficients<Lexicographical, T, Index1D> frame = extract(u, idx);

        for (size_type k=1; k<=frame.rank(); ++k) {
            for (const auto& it : frame(k, 1)) {
                FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                if (it.first.xtype==XWavelet) ++level;
                jmax(j) = MAX(level, jmax(j));
            }
        }
    }

    T sum = 0;
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=u.dim(); ++i) {
        sum += std::pow(2., 2.*order*jmax(i));
    }

    return sum/compOmegamin2(u.basis(), u.dim(), order);
}


template <typename T, typename Basis>
T
compIndexscale2(const HTCoefficients<T, Basis>& u, const T order)
{
    typedef typename SepCoefficients<Lexicographical, T, Index1D>::size_type
                     size_type;

    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >    jmax(u.dim());
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=u.dim(); ++j) {
        htucker::DimensionIndex idx(1);
        idx[0] = j;
        SepCoefficients<Lexicographical, T, Index1D> frame = extract(u, idx);

        for (size_type k=1; k<=frame.rank(); ++k) {
            for (const auto& it : frame(k, 1)) {
                FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                if (it.first.xtype==XWavelet) ++level;
                jmax(j) = MAX(level, jmax(j));
            }
        }
    }

    T sum = 0;
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=u.dim(); ++i) {
        sum += std::pow(2., 4.*order*jmax(i));
    }

    return sum/compOmegamin4(u.basis(), u.dim(), order);
}


template <typename Basis>
typename Basis::T
compIndexscale2(const Basis& basis,
                const std::vector<IndexSet<Index1D> >& Lambda,
                const typename Basis::T order)
{
    typedef typename SepCoefficients<Lexicographical, typename Basis::T,
                     Index1D>::size_type
                     size_type;

    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > jmax(Lambda.size());
    for (size_type j=1; j<=Lambda.size(); ++j) {
        for (const auto& it : Lambda[j-1]) {
            FLENS_DEFAULT_INDEXTYPE level = it.j;
            if (it.xtype==XWavelet) ++level;
            jmax(j) = MAX(level, jmax(j));
        }
    }

    typename Basis::T sum = 0;
    for (FLENS_DEFAULT_INDEXTYPE i=1; (unsigned) i<=Lambda.size(); ++i) {
        sum += std::pow(2., 4.*order*jmax(i));
    }

    return sum/compOmegamin4(basis, Lambda.size(), order);
}


template <typename Basis>
typename Basis::T
compIndexscale(const Basis& basis,
               const std::vector<IndexSet<Index1D> >& Lambda,
               const typename Basis::T order)
{
    assert(Lambda.size()>0);

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >    jmax(Lambda.size());
    for (size_type j=1; j<=Lambda.size(); ++j) {
        for (const auto& it : Lambda[j-1]) {
            FLENS_DEFAULT_INDEXTYPE level = it.j;
            if (it.xtype==XWavelet) ++level;
            jmax(j) = MAX(level, jmax(j));
        }
    }

    typename Basis::T sum = 0;
    for (unsigned i=1; i<=Lambda.size(); ++i) {
        sum += std::pow(2., 2.*order*jmax(i));
    }

    return sum/compOmegamin2(basis, Lambda.size(), order);
}


template <typename Basis>
void
setScaling(Sepdiagscal<Basis>& S,
           const typename Sepdiagscal<Basis>::T eps)
{
    S.set_eps(eps);
    S.set_nu(eps/2.);
    S.comp_h();
    S.comp_nplus();
}


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepdiagscal<Basis>& S,
     const SepCoefficients<Lexicographical, T, Index1D>& u)
{
    assert(S.dim()==u.dim());

    typedef typename Sepdiagscal<Basis>::size_type size_type;

    T iscale = compIndexscale(u, S.basis(), S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    SepCoefficients<Lexicographical, T, Index1D> sum;
    T factor1 = S.h()*(1./std::sqrt(omega2));
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i) {
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1.+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor2*factor1, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, S.dim());

            for (size_type j=1; j<=S.dim(); ++j) {
                prod(1, j) = u(k, j);

                for (auto& it : prod(1, j)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
operator*(Sepdiagscal<Basis>& S,
          const SepCoefficients<Lexicographical, T, Index1D>& u)
{
    assert(S.dim()==u.dim());
    return eval(S, u);
}


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepdiagscal<Basis>& S,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::vector<IndexSet<Index1D> >& cols)
{
    assert(S.dim()==u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type size_type;

    T iscale = compIndexscale(u, S.basis(), S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    SepCoefficients<Lexicographical, T, Index1D> sum;
    T factor1 = S.h()*(1./std::sqrt(omega2));
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i) {
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1.+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor2*factor1, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type k=1; k<=u.rank(); ++k) {
            SepCoefficients<Lexicographical, T, Index1D> prod(1, S.dim());

            for (size_type j=1; j<=S.dim(); ++j) {
                prod(1, j) = u(k, j);
                P(cols[j-1], prod(1, j));

                for (auto& it : prod(1, j)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            // Too many resizes here
            sum = sum+prod;
        }
    }

    return sum;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval(Sepdiagscal<Basis>& S,
     const HTCoefficients<T, Basis>& u,
     const double eps)
{
    assert(S.dim()==(unsigned) u.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale(u, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    HTCoefficients<T, Basis> sum;
    T factor1 = S.h()*(1./std::sqrt(omega2));
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor2*factor1, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type j=1; j<=S.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, idx);

            for (size_type k=1; k<=frame.rank(); ++k) {
                for (auto& it : frame(k, 1)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            set_inplace(prod, idx, frame);
        }

        if (i==-1*(signed) S.n()) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            sum.tree() = sum.tree()+prod.tree();
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "eval(S): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif
    return sum;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
operator*(Sepdiagscal<Basis>& S,
          const HTCoefficients<T, Basis>& u)
{
    assert(S.dim()==(unsigned) u.dim());
    return eval(S, u);
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval(Sepdiagscal<Basis>& S,
     HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >& cols,
     const double eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis(), u.map());
    sum.tree().set_tree(u.tree());

    unsigned count=0;
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i,
                                                                    ++count) {
        HTCoefficients<T, Basis> prod(u);
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type j=1; j<=S.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, cols[j-1], idx);

            for (size_type k=1; k<=frame.rank(); ++k) {
                P(cols[j-1], frame(k, 1));
                for (auto& it : frame(k, 1)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            set(prod, idx, frame);
        }

        if (i==-1*(signed) S.n()) {
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            //sum.tree() = sum.tree()+prod.tree();
            //sum.tree().truncate(eps);
            sum.tree() = add_truncate(sum.tree(), prod.tree(), eps);
            #ifdef DEBUG_CANCEL
                sumexact.tree() = sumexact.tree()+prod.tree();
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        }
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "eval(S, cols): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif

    return sum;
}


template <typename T, typename Basis>
void
scale(      Sepdiagscal<Basis>&              S,
            HTCoefficients<T, Basis>&        u,
      const std::vector<IndexSet<Index1D> >& cols,
      const FLENS_DEFAULT_INDEXTYPE          l)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());
    assert(l>=-1*(signed) S.n() && l<=(signed) S.nplus());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T omega2  = compOmegamin2(S.basis(), S.dim(), S.order());
    T factor1 = S.h()*(1./std::sqrt(omega2));
    T factor2 = 2.*std::pow(M_PI, -0.5)*
                (1./(1+std::exp(-1.*S.h()*(T) l)));
    factor2   = std::pow(factor1*factor2, 1./(T) S.dim());
    T alpha   = std::pow(std::log(1.+std::exp((T) l*S.h())), 2.);

    for (size_type j=1; j<=S.dim(); ++j) {
        htucker::DimensionIndex idx(1);
        idx[0] = j;

        SepCoefficients<Lexicographical, T, Index1D>
        frame = extract(u, cols[j-1], idx);

        for (size_type k=1; k<=frame.rank(); ++k) {
            P(cols[j-1], frame(k, 1));
            for (auto& it : frame(k, 1)) {
                FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                if (it.first.xtype==XWavelet) ++level;
                T weight = std::pow(2., 2.*S.order()*level)/omega2;
                it.second *= factor2*std::exp(-alpha*weight);
            }
        }
        set(u, idx, frame);
    }
}


template <typename T, typename Basis, typename Optype>
HTCoefficients<T, Basis>
eval(      Sepop<Optype>&                   A,
           Sepdiagscal<Basis>&              Srows,
           HTCoefficients<T, Basis>&        u,
     const std::vector<IndexSet<Index1D> >& rows,
     const std::vector<IndexSet<Index1D> >& cols,
     const T                                eps)
{
    assert(Srows.dim()==(unsigned) u.dim());
    assert(rows.size()==Srows.dim());
    assert(cols.size()==Srows.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), rows, Srows.order());
    Srows.set_iscale(iscale);
    Srows.comp_n();

    auto Scols = Srows;
    iscale     = compIndexscale(u.basis(), cols, Scols.order());
    Scols.set_iscale(iscale);
    Scols.comp_n();

    /* Precompute summands */
    auto Nrows = Srows.n()+Srows.nplus()+1;
    auto Ncols = Scols.n()+Scols.nplus()+1;
    flens::DenseVector<flens::Array<T>>                        nrms(Nrows*Ncols);
    std::vector<HTCoefficients<T, Basis>>                      prods(Nrows*Ncols);
    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE>>  ids(Nrows*Ncols);
    unsigned count = 0;
    for (FLENS_DEFAULT_INDEXTYPE l1=-1*Scols.n();
                                 l1<=(signed) Scols.nplus();
                                 ++l1) {
        auto v = u;
        scale(Scols, v, cols, l1);
        v = eval(A, v, rows, cols);
        for (FLENS_DEFAULT_INDEXTYPE l0=-1*Srows.n();
                                     l0<=(signed) Srows.nplus();
                                     ++l0) {
            auto tmp = v;
            scale(Srows, tmp, rows, l0);
            nrms(count+1) = nrm2(tmp);
            prods[count]  = tmp;
            ++count;
        }
    }

    /* Sort norms */
    flens::sort(nrms, ids);
    count    = 0;
    T cutoff = 0.;
    for (; count<(unsigned) nrms.length(); ++count) {
        cutoff += nrms(count+1);
        if (cutoff>eps/2) break;
    }

    if (count==(unsigned) nrms.length()) {
        std::cerr << "eval: Warning! Truncation parameter too large!\n";
        count = nrms.length()-1;
    }

    /* Add and truncate */
    HTCoefficients<T, Basis> sum(u.dim(), u.basis(), u.map());
    T refsum = 0.;
    for (unsigned i=count+1; i<=(unsigned) nrms.length(); ++i) {
        refsum += (T) (nrms.length()-i+1)*nrms(i);
    }

    sum    = prods[ids(count+1)-1];
    T eps_ = nrms(count+1);
    sum.truncate(eps/2*eps_/refsum);
    ++count;
    for (; count<(unsigned) nrms.length(); ++count) {
        eps_      += nrms(count+1);
        sum.tree() = add_truncate(sum.tree(), prods[ids(count+1)-1].tree(),
                     eps/2*eps_/refsum);
    }

    return sum;
}


template <typename T, typename Basis, typename Optype>
HTCoefficients<T, Basis>
evaleff(      Sepop<Optype>&                   A,
              Sepdiagscal<Basis>&              Srows,
              HTCoefficients<T, Basis>&        u,
        const std::vector<IndexSet<Index1D> >& rows,
        const std::vector<IndexSet<Index1D> >& cols,
        const T                                eps)
{
    assert(Srows.dim()==(unsigned) u.dim());
    assert(rows.size()==Srows.dim());
    assert(cols.size()==Srows.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), rows, Srows.order());
    Srows.set_iscale(iscale);
    Srows.comp_n();

    auto Scols = Srows;
    iscale     = compIndexscale(u.basis(), cols, Scols.order());
    Scols.set_iscale(iscale);
    Scols.comp_n();

    /* Distribute tolerances */
    T epsR = 0.5*eps;
    T epsL = 0.5*eps;

    /* Precompute summands right */
    auto N = Scols.n()+Scols.nplus()+1;
    flens::DenseVector<flens::Array<T>>                        nrms(N);
    std::vector<HTCoefficients<T, Basis>>                      prods(N);
    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE>>  ids(N);

    unsigned count = 0;
    for (FLENS_DEFAULT_INDEXTYPE l1=-1*Scols.n();
                                 l1<=(signed) Scols.nplus();
                                 ++l1) {
        auto v = u;
        scale(Scols, v, cols, l1);
        nrms(count+1) = nrm2(v);
        prods[count]  = v;
        ++count;
    }

    /* Presort small norms */
    T bound = compOmegamax2(cols, Scols.order());
    flens::sort(nrms, ids);
    count    = 0;
    T cutoff = 0.;
    for (; count<(unsigned) nrms.length(); ++count) {
        cutoff += bound*nrms(count+1);
        if (cutoff>epsR/2) break;
    }

    if (count==(unsigned) nrms.length()) {
        std::cerr << "eval: Warning! Truncation parameter too large!\n";
        count = nrms.length()-1;
    }

    std::vector
    <HTCoefficients<T, Basis>>               presortedprods(N-count);
    flens::DenseVector
    <flens::Array<T>>                        presortednrms(N-count);
    flens::DenseVector
    <flens::Array<FLENS_DEFAULT_INDEXTYPE>>  presortedids(N-count);
    for (unsigned cnt=1; count<(unsigned) nrms.length(); ++count, ++cnt) {
        auto v                = eval(A, prods[ids(count+1)-1], rows, cols);
        presortednrms(cnt)    = nrm2(v);
        presortedprods[cnt-1] = v;
    }

    /* Clean up to save temporary memory usage */
    prods.resize(0);
    nrms.resize(0);
    ids.resize(0);

    /* Sort norms */
    flens::sort(presortednrms, presortedids);
    count  = 0;
    cutoff = 0.;
    for (; count<(unsigned) presortednrms.length(); ++count) {
        cutoff += presortednrms(count+1);
        if (cutoff>epsR/2) break;
    }

    if (count==(unsigned) presortednrms.length()) {
        std::cerr << "eval: Warning! Truncation parameter too large!\n";
        count = presortednrms.length()-1;
    }

    /* Add and truncate */
    HTCoefficients<T, Basis> sum(u.dim(), u.basis(), u.map());
    T refsum = 0.;
    for (unsigned i=count+1; i<=(unsigned) presortednrms.length(); ++i) {
        refsum += (T) (presortednrms.length()-i+1)*presortednrms(i);
    }

    sum    = presortedprods[presortedids(count+1)-1];
    T eps_ = presortednrms(count+1);
    sum.truncate(epsR/2*eps_/refsum);
    ++count;
    for (; count<(unsigned) presortednrms.length(); ++count) {
        eps_      += presortednrms(count+1);
        sum.tree() = add_truncate(sum.tree(),
                     presortedprods[presortedids(count+1)-1].tree(),
                     epsR/2*eps_/refsum);
    }

    /* Scale left */
    sum = applyScale(Srows, sum, rows, epsL);

    return sum;
}


template <typename T, typename Basis, typename Optype>
HTCoefficients<T, Basis>
evaleff2(     Sepop<Optype>&                   A,
              Sepdiagscal<Basis>&              Srows,
              HTCoefficients<T, Basis>&        u,
        const std::vector<IndexSet<Index1D> >& rows,
        const std::vector<IndexSet<Index1D> >& cols,
        const T                                eps)
{
    assert(Srows.dim()==(unsigned) u.dim());
    assert(rows.size()==Srows.dim());
    assert(cols.size()==Srows.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), rows, Srows.order());
    Srows.set_iscale(iscale);
    Srows.comp_n();

    auto Scols = Srows;
    iscale     = compIndexscale(u.basis(), cols, Scols.order());
    Scols.set_iscale(iscale);
    Scols.comp_n();

    /* Distribute tolerances */
    T epsR = 0.95*eps;
    T epsL = 0.05*eps;

    /* Scale right */
    T bound  = 0.001*compOmegamax2(cols, Scols.order());
    bound    = std::max(bound, 1.);
    auto v   = applyScale(Scols, u, cols, epsR/std::sqrt(bound));

    /* Apply operator */
    v = eval(A, v, rows, cols);

    /* Scale left */
    v = applyScale(Srows, v, rows, epsL);

    return v;
}


template <typename T, typename Basis, typename Optype>
HTCoefficients<T, Basis>
evaleff2TT(     Sepop<Optype>&                   A,
                Sepdiagscal<Basis>&              Srows,
                HTCoefficients<T, Basis>&        u,
          const std::vector<IndexSet<Index1D> >& rows,
          const std::vector<IndexSet<Index1D> >& cols,
          const T                                eps)
{
    assert(Srows.dim()==(unsigned) u.dim());
    assert(rows.size()==Srows.dim());
    assert(cols.size()==Srows.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), rows, Srows.order());
    Srows.set_iscale(iscale);
    Srows.comp_n();

    auto Scols = Srows;
    iscale     = compIndexscale(u.basis(), cols, Scols.order());
    Scols.set_iscale(iscale);
    Scols.comp_n();

    /* Distribute tolerances */
    T epsR = 0.95*eps;
    T epsL = 0.05*eps;

    /* Scale right */
    T bound  = 0.001*compOmegamax2(cols, Scols.order());
    bound    = std::max(bound, 1.);
    auto v   = applyScaleTT(Scols, u, cols, epsR/std::sqrt(bound));

    /* Apply operator */
    v = eval(A, v, rows, cols);

    /* Scale left */
    v = applyScaleTT(Srows, v, rows, epsL);

    return v;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
applyScale(      Sepdiagscal<Basis>&              S,
                 HTCoefficients<T, Basis>&        u,
           const std::vector<IndexSet<Index1D> >& cols,
           const T                                eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();

    /* Precompute summands */
    auto N = S.n()+S.nplus()+1;
    flens::DenseVector<flens::Array<T>>                        nrms(N);
    std::vector<HTCoefficients<T, Basis>>                      prods(N);
    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE>>  ids(N);
    unsigned count = 0;
    for (FLENS_DEFAULT_INDEXTYPE l0=-1*S.n();
                                 l0<=(signed) S.nplus();
                                 ++l0) {
        auto v = u;
        scale(S, v, cols, l0);
        nrms(count+1) = nrm2(v);
        prods[count]  = v;
        ++count;
    }

    /* Sort norms */
    flens::sort(nrms, ids);
    count    = 0;
    T cutoff = 0.;
    for (; count<(unsigned) nrms.length(); ++count) {
        cutoff += nrms(count+1);
        if (cutoff>eps/2) break;
    }

    if (count==(unsigned) nrms.length()) {
        std::cerr << "eval: Warning! Truncation parameter too large!\n";
        count = nrms.length()-1;
    }

    /* Add and truncate */
    HTCoefficients<T, Basis> sum(u.dim(), u.basis(), u.map());
    T refsum = 0.;
    for (unsigned i=count+1; i<=(unsigned) nrms.length(); ++i) {
        refsum += (T) (nrms.length()-i+1)*nrms(i);
    }

    sum    = prods[ids(count+1)-1];
    T eps_ = nrms(count+1);
    sum.truncate(eps/2*eps_/refsum);
    ++count;
    for (; count<(unsigned) nrms.length(); ++count) {
        eps_      += nrms(count+1);
        sum.tree() = add_truncate(sum.tree(), prods[ids(count+1)-1].tree(),
                     eps/2*eps_/refsum);
    }

    return sum;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
applyScaleTT(      Sepdiagscal<Basis>&              S,
                   HTCoefficients<T, Basis>&          u,
             const std::vector<IndexSet<Index1D> >& cols,
             const T                                eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    /* Compute scales */
    T iscale = compIndexscale(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    std::cout << "S =>\n" << S << std::endl;

    /* Precompute summands */
    auto N = S.n()+S.nplus()+1;
    flens::DenseVector<flens::Array<T>>                        nrms(N);
    std::vector<HTCoefficients<T, Basis>>                      prods(N);
    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE>>  ids(N);
    unsigned count = 0;
    for (FLENS_DEFAULT_INDEXTYPE l0=-1*S.n();
                                 l0<=(signed) S.nplus();
                                 ++l0) {
        auto v = u;
        scale(S, v, cols, l0);
        nrms(count+1) = nrm2(v);
        prods[count]  = v;
        ++count;
    }

    /* Sort norms */
    flens::sort(nrms, ids);
    count    = 0;
    T cutoff = 0.;
    for (; count<(unsigned) nrms.length(); ++count) {
        cutoff += nrms(count+1);
        if (cutoff>eps/2) break;
    }

    if (count==(unsigned) nrms.length()) {
        std::cerr << "eval: Warning! Truncation parameter too large!\n";
        count = nrms.length()-1;
    }

    /* Add and truncate */
    HTCoefficients<T, Basis> sum(u.dim(), u.basis(), u.map());
    T refsum = 0.;
    for (unsigned i=count+1; i<=(unsigned) nrms.length(); ++i) {
        refsum += (T) (nrms.length()-i+1)*nrms(i);
    }

    sum    = prods[ids(count+1)-1];
    T eps_ = nrms(count+1);
    sum.tree().truncate(eps/2*eps_/refsum);
    ++count;
    for (; count<(unsigned) nrms.length(); ++count) {
        eps_      += nrms(count+1);
        sum.tree() = sum.tree()+prods[ids(count+1)-1].tree();
        sum.tree().truncate(eps/2*eps_/refsum);
    }

    return sum;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval_notrunc(Sepdiagscal<Basis>& S,
             HTCoefficients<T, Basis>& u,
             const std::vector<IndexSet<Index1D> >& cols)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis(), u.map());
    sum.tree().set_tree(u.tree());

    unsigned count=0;
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i,
                                                                    ++count) {
        HTCoefficients<T, Basis> prod(u);
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type j=1; j<=S.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, cols[j-1], idx);

            for (size_type k=1; k<=frame.rank(); ++k) {
                P(cols[j-1], frame(k, 1));
                for (auto& it : frame(k, 1)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            set(prod, idx, frame);
        }

        if (i==-1*(signed) S.n()) {
            sum = prod;
        } else {
            sum.tree() = sum.tree()+prod.tree();
        }
    }

    return sum;
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
fixeval_notrunc(Sepdiagscal<Basis>&                    S,
                HTCoefficients<T, Basis>&              u,
                const std::vector<IndexSet<Index1D> >& cols)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());
    T factor1 = S.h()*(1./std::sqrt(omega2));

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis(), u.map());

    unsigned count=0;
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i,
                                                                    ++count) {
        HTCoefficients<T, Basis> prod(u);
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type j=1; j<=S.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, cols[j-1], idx);

            for (size_type k=1; k<=frame.rank(); ++k) {
                P(cols[j-1], frame(k, 1));
                for (auto& it : frame(k, 1)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 2.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            set(prod, idx, frame);
        }

        if (i==-1*(signed) S.n()) {
            sum = prod;
        } else {
            sum.tree() = sum.tree()+prod.tree();
        }
    }

    return sum;
}


template <typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepdiagscal<Basis>&              S,
     const flens::GeMatrix
           <flens::FullStorage
           <T, cxxblas::ColMajor>>&         Uj,
           HTCoefficients<T, Basis>&        u,
     const unsigned                         j,
     const std::vector<IndexSet<Index1D> >& cols)
{

    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage
                     <T, cxxblas::ColMajor> >       Matrix;

    T iscale = compIndexscale(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));

    auto ncolsU = Uj.numCols();
    auto block  = S.nplus()+S.n()+1;
    auto ncols  = ncolsU*block;
    auto nrows  = Uj.numRows();
    Matrix ret(nrows, ncols);
    for (auto& lambda : cols[j-1]) {
        auto row = u.map()(lambda, j);
        assert(row<=(unsigned) nrows);

        unsigned cnt = 0;
        for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus();
                                     ++i, ++cnt) {
            T factor2 = 2.*std::pow(M_PI, -0.5)*
                        (1./(1+std::exp(-1.*S.h()*(T) i)));
            factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
            T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

            ret(_(1, nrows), _(cnt*ncolsU+1, (cnt+1)*ncolsU)) = Uj;
            for (FLENS_DEFAULT_INDEXTYPE k=1; k<=ncolsU; ++k) {
                FLENS_DEFAULT_INDEXTYPE level = lambda.j;
                if (lambda.xtype==XWavelet) ++level;
                T weight               = std::pow
                                         (2., 2.*S.order()*level)/omega2;
                ret(row, cnt*ncolsU+k) *= factor2*std::exp(-alpha*weight);
            }
        }
    }

    return ret;
}


template <typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
fixeval(      Sepdiagscal<Basis>&          S,
        const flens::GeMatrix
              <flens::FullStorage
              <T, cxxblas::ColMajor>>&     Uj,
              HTCoefficients<T, Basis>&    u,
        const unsigned                     j,
        const IndexSet<Index1D>&           cols)
{

    assert(S.dim()==(unsigned) u.dim());
    using flens::_;

    typedef typename flens::GeMatrix
                     <flens::FullStorage
                     <T, cxxblas::ColMajor> >       Matrix;

    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());
    T factor1 = S.h()*(1./std::sqrt(omega2));

    auto ncolsU = Uj.numCols();
    auto block  = S.nplus()+S.n()+1;
    auto ncols  = ncolsU*block;
    auto nrows  = Uj.numRows();
    Matrix ret(nrows, ncols);
    for (auto& lambda : cols) {
        auto row = u.map()(lambda, j);
        assert(row<=(unsigned) nrows);

        unsigned cnt = 0;
        for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus();
                                     ++i, ++cnt) {
            T factor2 = 2.*std::pow(M_PI, -0.5)*
                        (1./(1+std::exp(-1.*S.h()*(T) i)));
            factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
            T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

            ret(_(1, nrows), _(cnt*ncolsU+1, (cnt+1)*ncolsU)) = Uj;
            for (FLENS_DEFAULT_INDEXTYPE k=1; k<=ncolsU; ++k) {
                FLENS_DEFAULT_INDEXTYPE level = lambda.j;
                if (lambda.xtype==XWavelet) ++level;
                T weight               = std::pow
                                         (2., 2.*S.order()*level)/omega2;
                ret(row, cnt*ncolsU+k) *= factor2*std::exp(-alpha*weight);
            }
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
prec(Precon&                                                     P,
     const
     flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
     HTCoefficients<T, Basis>&                                   u,
     const unsigned                                              j,
     const IndexSet<Index1D>&                                    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = U(i, k)*p;
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
prec(Precon&                                                     P,
     const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& Pj,
     const
     flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
     HTCoefficients<T, Basis>&                                   u,
     const unsigned                                              j,
     const IndexSet<Index1D>&                                    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    T alpha = Pj(1, 2);
    T beta  = Pj(1, 1);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);

        auto l = lambda.j;
        if (lambda.xtype==XWavelet) ++l;
        p      = std::pow(alpha*(1./(p*p)-1.)+beta, -0.5);
        p     *= std::pow(2., -l/2.);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = U(i, k)*p;
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
precsq(Precon&                                                     P,
       const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& Pj,
       const
       flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
       HTCoefficients<T, Basis>&                                   u,
       const unsigned                                              j,
       const IndexSet<Index1D>&                                    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    T alpha = Pj(1, 2);
    T beta  = Pj(1, 1);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);

        auto l = lambda.j;
        if (lambda.xtype==XWavelet) ++l;
        p      = std::pow(alpha*(1./(p*p)-1.)+beta, -1.);
        p     *= std::pow(2., -l);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = U(i, k)*p;
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
remove_prec(Precon&                                                     P,
            const
            flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
            HTCoefficients<T, Basis>&                                   u,
            const unsigned                                              j,
            const IndexSet<Index1D>&                                    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = U(i, k)/p;
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
precsq(Precon&                                                     P,
       const
       flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
       HTCoefficients<T, Basis>&                                   u,
       const unsigned                                              j,
       const IndexSet<Index1D>&                                    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = p*p*U(i, k);
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::DenseVector<flens::Array<T> >
assemble_precsq(Precon&                     P,
                HTCoefficients<T, Basis>&   u,
                const unsigned              j,
                const IndexSet<Index1D>&    active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    htucker::DimensionIndex idx(1);
    idx[0] = j;

    auto& U    = extract(u.tree(), idx);
    auto  numr = U.numRows();

    flens::DenseVector<flens::Array<T> >    ret(numr);
    ret(numr);
    ret.fill((T) 1.);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #endif
        ret(i) = p*p;
    }

    return ret;
}


template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
apply_precsq(const flens::DenseVector<flens::Array<T> >&    prec,
             const flens::GeMatrix
                   <flens::FullStorage
                   <T, cxxblas::ColMajor> >&                U)
{
    assert(prec.length()==U.numRows());
    assert(U.numCols()>0);

    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
    ret(U.numRows(), U.numCols());

    for (long i=1; i<=U.numRows(); ++i) {
        for (long j=1; j<=U.numCols(); ++j) {
            ret(i, j) = U(i, j)*prec(i);
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
remove_precsq(Precon&                                      P,
              const
              flens::GeMatrix
              <flens::FullStorage<T, cxxblas::ColMajor> >& U,
              HTCoefficients<T, Basis>&                    u,
              const unsigned                               j,
              const IndexSet<Index1D>&                     active)
{
    assert(j>=1 && j<=(unsigned) u.dim());

    auto numr = U.numRows();
    auto numc = U.numCols();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>
    ret(numr, numc);
    for (auto& lambda : active) {
        auto i = u.map()(lambda, j);
        auto p = P(lambda);
        #ifndef NDEBUG
            assert(i<=(unsigned) numr);
        #else
            (void) numr;
        #endif
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
            ret (i, k) = U(i, k)/(p*p);
        }
    }

    return ret;
}


template <typename Precon, typename T, typename Basis>
void
rank1prec(Precon&                               P,
          HTCoefficients<T, Basis>&             u,
          const std::vector<IndexSet<Index1D>>& Lambda)
{
    assert(Lambda.size() == (unsigned) u.dim());
    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor>>  Matrix;

    for (auto it=u.tree().getGeneralTree().end();
              it>=u.tree().getGeneralTree().begin(); it--) {
        auto node  = it.getNode();

        if (node->isLeaf()) {
            htucker::DimensionIndex
                idx   = node->getContent()->getIndex();
            Matrix& U = const_cast<Matrix&>(node->getContent()->getUorB());
            auto numr = U.numRows();
            auto numc = U.numCols();
            Matrix Uc = U;
            U.fill((T) 0.);

            for (auto& lambda : Lambda[idx[0]-1]) {
                auto i = u.map()(lambda, idx[0]);
                auto p   = P(lambda);
                #ifndef NDEBUG
                    assert(i<=(unsigned) numr);
                #endif
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
                    U(i, k) = Uc(i, k)*p;
                }
            }
        }
    }
}


template <typename Precon, typename T, typename Basis>
void
rank1prec(Precon&                    P,
          HTCoefficients<T, Basis>&  u,
          const unsigned             j,
          const IndexSet<Index1D>&   Lambda)
{
    assert(j>=1 && j<= (unsigned) u.dim());
    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor>>  Matrix;

    for (auto it=u.tree().getGeneralTree().end();
              it>=u.tree().getGeneralTree().begin(); it--) {
        auto index  = it.getNode()->getContent()->getIndex();

        if ((unsigned) index[0]==j) {
            Matrix& U = const_cast<Matrix&>(it.getNode()
                        ->getContent()->getUorB());
            auto numr = U.numRows();
            auto numc = U.numCols();

            for (auto& lambda : Lambda) {
                auto i = u.map()(lambda, j);
                auto p = P(lambda);
                #ifndef NDEBUG
                    assert(i<=(unsigned) numr);
                #endif
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
                    U(i, k) *= p;
                }
            }

            break;
        }
    }
}


template <typename Precon, typename T, typename Basis>
void
remove_rank1prec(Precon&                               P,
                 HTCoefficients<T, Basis>&             u,
                 const std::vector<IndexSet<Index1D>>& Lambda)
{
    assert(Lambda.size() == (unsigned) u.dim());
    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor>>  Matrix;

    for (auto it=u.tree().getGeneralTree().end();
              it>=u.tree().getGeneralTree().begin(); it--) {
        auto node  = it.getNode();

        if (node->isLeaf()) {
            htucker::DimensionIndex
                idx   = node->getContent()->getIndex();
            Matrix& U = const_cast<Matrix&>(node->getContent()->getUorB());
            auto numr = U.numRows();
            auto numc = U.numCols();
            Matrix Uc = U;
            U.fill((T) 0.);

            for (auto& lambda : Lambda[idx[0]-1]) {
                auto i = u.map()(lambda, idx[0]);
                #ifndef NDEBUG
                    assert(i<=(unsigned) numr);
                #endif
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
                    auto p   = 1./P(lambda);
                    U(i, k)  = Uc(i, k)*p;
                }
            }
        }
    }
}


template <typename Precon, typename T, typename Basis>
void
remove_rank1prec(Precon&                    P,
                 HTCoefficients<T, Basis>&  u,
                 const unsigned             j,
                 const IndexSet<Index1D>&   Lambda)
{
    assert(j>=1 && j<= (unsigned) u.dim());
    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor>>  Matrix;

    for (auto it=u.tree().getGeneralTree().end();
              it>=u.tree().getGeneralTree().begin(); it--) {
        auto index  = it.getNode()->getContent()->getIndex();

        if ((unsigned) index[0]==j) {
            Matrix& U = const_cast<Matrix&>(it.getNode()
                        ->getContent()->getUorB());
            auto numr = U.numRows();
            auto numc = U.numCols();

            for (auto& lambda : Lambda) {
                auto i = u.map()(lambda, j);
                auto p = P(lambda);
                #ifndef NDEBUG
                    assert(i<=(unsigned) numr);
                #endif
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=numc; ++k) {
                    U(i, k) /= p;
                }
            }

            break;
        }
    }
}


template <typename T, typename Basis>
HTCoefficients<T, Basis>
evalS2(Sepdiagscal<Basis>& S,
       HTCoefficients<T, Basis>& u,
       const std::vector<IndexSet<Index1D> >& cols,
       const double eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale2(u.basis(), cols, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin4(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis(), u.map());
    sum.tree().set_tree(u.tree());

    unsigned N = S.n()+S.nplus()+1;
    std::vector<HTCoefficients<T, Basis> >  prods(N);
    flens::DenseVector<flens::Array<T> >    nrms(N);
    flens::DenseVector<flens::Array<int> >  ids(N);

    unsigned count=0;
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i,
                                                                    ++count) {
        HTCoefficients<T, Basis> prod(u);
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1.+std::exp(-1.*S.h()*(T) i)));
        factor2 = std::pow(factor1*factor2, 1./(T) S.dim());
        T alpha = std::pow(std::log(1.+std::exp((T) i*S.h())), 2.);

        for (size_type j=1; j<=S.dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, cols[j-1], idx);

            for (size_type k=1; k<=frame.rank(); ++k) {
                P(cols[j-1], frame(k, 1));
                for (auto& it : frame(k, 1)) {
                    FLENS_DEFAULT_INDEXTYPE level = it.first.j;
                    if (it.first.xtype==XWavelet) ++level;
                    T weight = std::pow(2., 4.*S.order()*level)/omega2;
                    it.second *= factor2*std::exp(-alpha*weight);
                }
            }
            set(prod, idx, frame);
        }

        nrms(count+1) = nrm2(prod);
        prods[count]  = prod;
    }

    #ifdef DEBUG_CANCEL
        sumexact.orthogonalize();
        std::cout << "eval(S, cols): kappa = "
                  << sumnorms/(T) sumexact.tree().L2normorthogonal()
                  << std::endl;
    #endif

    flens::sort(nrms, ids);
    count    = 0;
    T cutoff = 0.;
    for (; count<nrms.length(); ++count) {
        cutoff += nrms(count+1);
        if (cutoff>eps/2) break;
    }

    if (count==nrms.length()) {
        std::cerr << "evalS2: Warning! Truncation parameter too large!\n";
        count = nrms.length()-1;
    }

    T refsum = 0.;
    for (unsigned i=count+1; i<=nrms.length(); ++i) {
        refsum += (T) (nrms.length()-i+1)*nrms(i);
    }

    sum    = prods[ids(count+1)-1];
    T eps_ = nrms(count+1);
    sum.truncate(eps/2*eps_/refsum);
    ++count;
    for (;count<nrms.length(); ++count) {
        eps_ += nrms(count+1);
        sum.tree() = add_truncate(sum.tree(), prods[ids(count+1)-1].tree(),
                     eps/2*eps_/refsum);
    }

    return sum;
}


template <typename T, typename Basis>
void
assemble(Sepdiagscal<Basis>& S,
         HTCoefficients<T, Basis>& Stree,
         const std::vector<IndexSet<Index1D> >& cols)
{
    assert(S.dim()==(unsigned) Stree.dim());
    assert(S.dim()==cols.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;
    typedef typename flens::DenseVector<flens::Array<T> >       DV;
    typedef typename htucker::DenseVectorList<T>                DVList;

    DVList  list;

    /* Assemble HT tree */
    S.set_iscale(compIndexscale(S.basis(), cols, S.order()));
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());
    T factor1 = S.h()*(1./std::sqrt(omega2));

    for (size_type j=1; j<=cols.size(); ++j) {
        auto max = maxintindhash(cols[j-1], j, Stree);
        for (FLENS_DEFAULT_INDEXTYPE k=-1*S.n(); k<=(signed) S.nplus(); ++k) {
            DV vecA(max*max);

            T factor2 = 2.*std::pow(M_PI, -0.5)*
                        (1./(1+std::exp(-1.*S.h()*(T) k)));
            factor2   = std::pow(factor2*factor1, 1./(T) S.dim());
            T alpha   = std::pow(std::log(1.+std::exp((T) k*S.h())), 2.);
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=(signed) max; ++i) {
                auto index        = Stree.map()(i, j);
                FLENS_DEFAULT_INDEXTYPE level         = index.j;
                if (index.xtype==XWavelet) ++level;
                T weight          = std::pow(2., 2.*S.order()*level)/omega2;
                vecA((i-1)*max+i) = factor2*std::exp(-alpha*weight);
            }
            list.add(vecA);
        }
    }

    Stree.tree().generateTofElementary(list, S.n()+S.nplus()+1, S.dim());
}


template <typename Basis>
std::ostream& operator<<(std::ostream& s,
                         const Sepdiagscal<Basis>& S)
{
    s << "dim    = " << S.dim() << std::endl;
    s << "order  = " << S.order() << std::endl;
    s << "eps    = " << S.eps() << std::endl;
    s << "nu     = " << S.nu() << std::endl;
    s << "h      = " << S.h() << std::endl;
    s << "iscale = " << S.iscale() << std::endl;
    s << "nplus  = " << S.nplus() << std::endl;
    s << "n      = " << S.n();

    return s;
}


template <typename T, typename Basis, typename Index>
Coefficients<Lexicographical, T, Index>
contraction(      HTCoefficients<T, Basis>& u,
            const IndexSet<Index>& activex,
            const flens::DenseVector<flens::Array<T> >& sigmas,
            const FLENS_DEFAULT_INDEXTYPE dim)
{
    assert(dim>=1 && dim<=u.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >  Matrix;

    htucker::DimensionIndex idx(1);
    idx[0] = dim;

    Coefficients<Lexicographical, T, Index> ret;
    for(auto tit=u.tree().getGeneralTree().end();
             tit>=u.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            const Matrix& U = tit.getNode()->getContent()->getUorB();
            FLENS_DEFAULT_INDEXTYPE numRows = U.numRows();
            FLENS_DEFAULT_INDEXTYPE numCols = U.numCols();
            assert(sigmas.length()>=numCols);

            for (const auto& mu : activex) {
                FLENS_DEFAULT_INDEXTYPE rowi = u.map()(mu, dim);
                assert(rowi<=numRows);
                T entry = 0.;
                for (FLENS_DEFAULT_INDEXTYPE j=1; j<=numCols; ++j) {
                    entry += sigmas(j)*U(rowi, j)*U(rowi, j);
                }
                ret[mu] = std::sqrt(entry);
            }

            return ret;
        }
    }

    std::cerr << "error contraction: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis, typename Index>
void
contraction(      HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index> >& activex,
            const std::vector<flens::DenseVector<flens::Array<T> > >& sigmas,
            std::vector<Coefficients<Lexicographical, T, Index> >& ret)
{
    assert(activex.size()==(unsigned) u.dim() &&
           sigmas.size()==(unsigned) u.dim());
    if (!ret.size()) ret.resize(u.dim());
    assert(ret.size()==(unsigned) u.dim());

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=u.dim(); ++j) {
        ret[j-1] = contraction(u, activex[j-1], sigmas[j-1], j);
    }
}


template <typename T, typename Basis>
void
scal(const T alpha, HTCoefficients<T, Basis>& x)
{
    x.tree().scal(alpha);
}


template <typename T, typename Basis>
T
dot(const HTCoefficients<T, Basis>& x, const HTCoefficients<T, Basis>& y)
{
    assert(x.dim()==y.dim());
    return x.tree().ScalarProduct(y.tree());
}


template <typename T, typename Basis>
T
nrm2(HTCoefficients<T, Basis>& x, bool isorth)
{
    if (!isorth) x.tree().orthogonalize();
    return x.tree().L2normorthogonal();
}


template <typename T, typename Basis>
void
axpy(const T alpha,
     const HTCoefficients<T, Basis>& x,
           HTCoefficients<T, Basis>& y)
{
    assert(x.dim()==y.dim());

    if (alpha==(T) 0) return;

    y.tree() = alpha*x.tree()+y.tree();
}


template <typename T, typename Basis>
void
restrict(HTCoefficients<T, Basis>& f,
         const IndexSet<Index1D>& activex,
         const unsigned j)
{
    assert(j>0 && j<=(unsigned) f.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    htucker::DimensionIndex idx(1);
    idx[0] = j;

    SepCoefficients<Lexicographical, T, Index1D>
    frame = extract(f, activex, idx);

    for (size_type k=1; k<=frame.rank(); ++k) {
        P(activex, frame(k, 1));
    }

    set(f, idx, frame);
}


template <typename T, typename Basis>
void
restrict(HTCoefficients<T, Basis>& f,
         const std::vector<IndexSet<Index1D> >& activex)
{
    assert(activex.size()==(unsigned) f.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    for (auto tit=f.tree().getGeneralTree().end();
              tit>=f.tree().getGeneralTree().begin();
              tit--) {
        auto node = tit.getNode();

        if (node->isLeaf()) {
            htucker::DimensionIndex idx = node->getContent()->getIndex();
            Matrix& U = const_cast<Matrix&>(node->getContent()->getUorB());
            auto max = maxintindhash(activex[idx[0]-1], idx[0], f);
            assert(max<=(unsigned) U.numRows());

            Matrix copy = U;
            if (max!=(unsigned) U.numRows()) {
                U.resize(max, U.numCols());
            } else {
                U.fill((T) 0);
            }

            for (auto& lambda : activex[idx[0]-1]) {
                auto i = f.map()(lambda, idx[0]);
                for (FLENS_DEFAULT_INDEXTYPE k=1; k<=U.numCols(); ++k) {
                    U(i, k) = copy(i, k);
                }
            }
        }
    }
}


template <typename T, typename Basis>
void
extend(HTCoefficients<T, Basis>& f,
       const IndexSet<Index1D>& activex,
       const unsigned j)
{
    assert(j>0 && j<=(unsigned) f.dim());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;
    using flens::_;

    htucker::DimensionIndex idx(1);
    idx[0] = j;

    for (auto tit=f.tree().getGeneralTree().end();
              tit>=f.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            FLENS_DEFAULT_INDEXTYPE      rowsold = U.numRows();
            FLENS_DEFAULT_INDEXTYPE      cols    = U.numCols();
            unsigned rowsnew = maxintindhash(activex, j, f);

            if (rowsnew>(unsigned) rowsold) {
                Matrix copy = U;
                U.resize(rowsnew, cols);
                U(_(1,rowsold), _) = copy;
            }

            return;
        }
    }

    std::cerr << "extend: idx not found\n";
    exit(EXIT_FAILURE);
}


template <typename T, typename Basis>
void
extend(HTCoefficients<T, Basis>& f,
       const std::vector<IndexSet<Index1D> >& activex)
{
    assert(activex.size()==(unsigned) f.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    for (auto tit=f.tree().getGeneralTree().end();
              tit>=f.tree().getGeneralTree().begin();
              tit--) {
        auto node = tit.getNode();

        if (node->isLeaf()) {
            htucker::DimensionIndex idx = node->getContent()->getIndex();
            Matrix& U = const_cast<Matrix&>(node->getContent()->getUorB());
            auto max = maxintindhash(activex[idx[0]-1], idx[0], f);
            assert(max>=(unsigned) U.numRows());

            if (max!=(unsigned) U.numRows()) {
                Matrix copy = U;
                U.resize(max, U.numCols());
                U(_(1, copy.numRows()), _) = copy;
            }
        }
    }
}


template <typename T, typename Optype, typename Basis>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> >>
reduce_laplace(      Sepop<Optype>&                    A,
                     HTCoefficients<T, Basis>&         U,
               const std::vector<IndexSet<Index1D> >&  rows,
               const std::vector<IndexSet<Index1D> >&  cols,
               const std::size_t                       hashtablelength)
{
    assert(A.dim()==(unsigned) U.dim());
    assert(A.dim()==rows.size());
    assert(A.dim()==cols.size());
    using flens::_;

    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  Matrix;

    std::vector<Matrix>     Bv(A.dim());

    for (auto it  = U.tree().getGeneralTree().end();
              it >= U.tree().getGeneralTree().begin();
              it--) {
        auto node  = it.getNode();

        /* Apply operator to leafs */
        if (node->isLeaf()) {
            htucker::DimensionIndex idx = node->getContent()->getIndex();
            Matrix& Uv  = node->getContent()->getUorB();
            auto rowsU  = Uv.numRows();
            auto colsU  = Uv.numCols();
            auto max    = maxintindhash(rows[idx[0]-1], U, idx[0]);
            assert(max>=(unsigned) rowsU);
            Matrix AU(max, colsU);
            AU(_(1, rowsU), _) = Uv;;

            /* Apply operator to columns */
            for (FLENS_DEFAULT_INDEXTYPE k=1; k<=(FLENS_DEFAULT_INDEXTYPE) colsU; ++k) {
                TreeCoefficients1D<T> input(hashtablelength, U.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength, U.basis().j0);
                typedef typename TreeCoefficients1D<T>::val_type val_type;

                /* Set input */
                for (auto& lambda : cols[idx[0]-1]) {
                    auto i = U.map()(lambda, idx[0]);
                    assert(i<=(unsigned) rowsU);

                    auto j     = lambda.j;
                    auto _k     = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        input.bylevel[j-1-input.offset].
                              map.insert(val_type(_k, Uv(i, k)));
                    } else {
                        input.bylevel[j-input.offset].
                              map.insert(val_type(_k, Uv(i, k)));
                    }
                }

                /* Set output */
                for (auto& lambda : rows[idx[0]-1]) {
                    #ifndef NDEBUG
                        auto i = U.map()(lambda, idx[0]);
                        assert(i<=max);
                    #endif

                    auto j     = lambda.j;
                    auto _k     = lambda.k;
                    auto xtype = lambda.xtype;
                    if (xtype==XBSpline) {
                        output.bylevel[j-1-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    } else {
                        output.bylevel[j-output.offset].
                              map.insert(val_type(_k, (T) 0));
                    }
                }

                /* Apply A */
                A(1, 1).eval(input, output, "A");

                /* Save result */
                for (typename CoefficientsByLevel<T>::const_it it=
                     output.bylevel[0].map.begin();
                     it!=output.bylevel[0].map.end(); ++it) {
                        Index1D lambda(output.offset+1, (*it).first,XBSpline);
                        auto i = U.map()(lambda, idx[0]);
                        AU(i, k) = (*it).second;
                }
                for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                    if (output.bylevel[i].map.size()==0) break;
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[i].map.begin();
                         it!=output.bylevel[i].map.end(); ++it) {
                            Index1D lambda(output.offset+i,
                                           (*it).first,XWavelet);
                            auto j = U.map()(lambda, idx[0]);
                            AU(j, k) = (*it).second;
                    }
                }
            }

            /* U^T*A*U */
            Matrix ret;
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1., Uv,
                            AU, 0., ret);
            Bv[idx[0]-1] = ret;
        }
    }

    return Bv;
}


template <typename T, typename Optype, typename Basis>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> >>
reduce(      Sepop<Optype>&                    A,
             HTCoefficients<T, Basis>&         U,
       const std::vector<IndexSet<Index1D> >&  rows,
       const std::vector<IndexSet<Index1D> >&  cols,
       const std::size_t                       hashtablelength)
{
   assert(A.dim()==(unsigned) U.dim());
   assert(A.dim()==rows.size());
   assert(A.dim()==cols.size());

   if (A.type()==laplace) {
       return reduce_laplace(A, U, rows, cols, hashtablelength);
   } else {
       std::cerr << "reduce: Not implemented for operator type\n";
       exit(1);
   }
}


template <typename T, typename Basis>
std::vector<flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > >
reduce_rhs(const HTCoefficients<T, Basis>&         U,
           const HTCoefficients<T, Basis>&         b)
{
    assert(U.dim()==b.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  Matrix;

    std::vector<Matrix>     c(U.dim());

    for (auto it  = U.tree().getGeneralTree().end(),
              itb = b.tree().getGeneralTree().end();
              it >= U.tree().getGeneralTree().begin();
              it--, itb--) {
        auto node  = it.getNode();
        auto nodeb = itb.getNode();

        if (node->isLeaf()) {
            htucker::DimensionIndex idx = node->getContent()->getIndex();
            Matrix& Uv  = node->getContent()->getUorB();
            Matrix& bv  = nodeb->getContent()->getUorB();
            Matrix ret;
            flens::blas::mm(cxxblas::Trans, cxxblas::NoTrans, 1., Uv,
                            bv, 0., ret);
            c[idx[0]-1] = ret;
        }
    }

    return c;
}


template <typename Optype, typename T, typename Basis>
flens::SyMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
assemble_projected_laplace(      Sepop<Optype>&             A,
                                 HTCoefficients<T, Basis>&  u,
                           const IndexSet<Index1D>&         Lambda,
                           const unsigned                   j)
{
    flens::SyMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
    ret(maxintindhash(Lambda, u, j), flens::Lower);

    auto& a = A(1, 1);
    for (auto& lambdaR : Lambda) {
        auto row = u.map()(lambdaR, j);
        for (auto& lambdaC: Lambda) {
            auto col      = u.map()(lambdaC, j);
            if (col>row) continue;
            ret(row, col) = a(lambdaR, lambdaC);
        }
    }

    return ret;
}


template <typename T, typename Basis>
flens::GeMatrix<
flens::FullStorage<T, flens::ColMajor> >
convert(const SepCoefficients<Lexicographical, T, Index1D>& cp,
              HTCoefficients<T, Basis>&                     tree,
        const unsigned                                      j)
{
    assert(j>=1 && j<=cp.dim());
    assert(cp.dim()==(unsigned) tree.dim());

    int rows = maxintindhash(cp(1, j), j, tree);
    flens::GeMatrix<
    flens::FullStorage<T, flens::ColMajor> > U(rows, cp.rank());

    for (unsigned i=1; i<=cp.rank(); ++i) {
        for (const auto& it : cp(i, j)) {
            U(tree.map()(it.first, j), i) = it.second;
        }
    }

    return U;
}


template <typename T, typename Basis>
flens::GeMatrix<
flens::FullStorage<T, flens::ColMajor> >
convert(      SepCoefficients<Lexicographical, T, Index1D>& cp,
        const IndexSet<Index1D>&                            active,
              HTCoefficients<T, Basis>&                     tree,
        const unsigned                                      j)
{
    assert(j>=1 && j<=cp.dim());
    assert(cp.dim()==(unsigned) tree.dim());

    int rows = maxintindhash(active, j, tree);
    flens::GeMatrix<
    flens::FullStorage<T, flens::ColMajor> > U(rows, cp.rank());

    for (const auto& it : active) {
        unsigned r = tree.map()(it, j);
        for (unsigned i=1; i<=cp.rank(); ++i) {
            U(r, i) = cp(i, j)[it];
        }
    }

    return U;
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
convert(const flens::GeMatrix<
              flens::FullStorage<T, flens::ColMajor> >& U,
              HTCoefficients<T, Basis>&                 tree,
        const IndexSet<Index1D>&                        active,
        const unsigned                                  j)
{
    assert(j>=1 && j<=(unsigned) tree.dim());
    assert(U.numCols()==1);

    Coefficients<Lexicographical, T, Index1D> v;

    for (const auto& it : active) {
        v[it] = U(tree.map()(it, j), 1);
    }

    return v;
}


template <typename Index>
std::ostream& operator<<(std::ostream& s,
                         const Mapwavind<Index>& map)
{
    for (unsigned j=0; j<map.dim(); ++j) {
        s << "*** dimension = " << j+1 << " ***\n";
        for (const auto& it : map.get_active()[j].right) {
            s << it.first << " = " << it.second << std::endl;
        }
        s << "\n";
    }

    return s;
}


template <typename Index>
FLENS_DEFAULT_INDEXTYPE
maxlevel(const IndexSet<Index>& Lambda)
{
    FLENS_DEFAULT_INDEXTYPE jmax = -100;
    for (const auto& it : Lambda) {
        FLENS_DEFAULT_INDEXTYPE level = it.j;
        if (it.xtype==XWavelet) ++level;
        jmax = MAX(level, jmax);
    }
    return jmax;
}


template <typename T, typename Basis>
void
insert(      HTCoefficients<T, Basis>& x,
       const Coefficients<Lexicographical, T, Index1D>& v,
       const IndexSet<Index1D>&                         active,
       const unsigned                                   j)
{
    assert(j>=1 && j<=(unsigned) x.dim());
    assert(v.size()==active.size());

    typedef typename flens::GeMatrix
                     <flens::FullStorage<T, flens::ColMajor> > Matrix;

    htucker::DimensionIndex idx(1);
    idx[0] = j;

    for (auto tit=x.tree().getGeneralTree().end();
              tit>=x.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->getContent()->getIndex()==idx) {
            Matrix& U = const_cast<Matrix&>
                        (tit.getNode()->getContent()
                         ->getUorB());
            assert(U.numCols()==1);
            auto max  = maxintindhash(active, x, j);
            U.resize(max, 1);

            for (const auto& it : v) {
                U(x.map()(it.first, j), 1) = it.second;
            }
        }
    }
}


template <typename Index>
std::vector<IndexSet<Index> >
unify(const std::vector<IndexSet<Index> >& A,
      const std::vector<IndexSet<Index> >& B)
{
    assert(A.size()==B.size());

    std::vector<IndexSet<Index> > C(A.size());

    for (unsigned j=0; j<C.size(); ++j) {
        C[j] = A[j] + B[j];
    }

    return C;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_TCC 1
