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
#include <stdlib.h>
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

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintind(coeff, tree.basis());
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
                size_type sizer = (rowsset>rowsnode) ? rowsset : rowsnode;
                size_type sizec = (colsset>colsnode) ? colsset : colsnode;
                U.resize(sizer, sizec);
            } else {
                U.fill((T) 0);
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

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintind(coeff, tree.basis());
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

            unsigned FLENS_DEFAULT_INDEXTYPE max = maxintind(coeff, tree.basis());
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
                FLENS_DEFAULT_INDEXTYPE rowi = maptoint(it.first, tree.basis());
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
                    FLENS_DEFAULT_INDEXTYPE rowi = maptoint(it.first, tree.basis());
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
            for (unsigned FLENS_DEFAULT_INDEXTYPE i=1; i<=(unsigned) U.numRows(); ++i) {
                Index1D index = maptowav(i, tree.basis());
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
            frame = extract(prod, idx);

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
        frame = extract(prod, idx);

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
        frame = extract(prod, idx);

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
           const HTCoefficients<T, Basis>& u,
           const std::vector<IndexSet<Index1D> >& rows,
           const std::vector<IndexSet<Index1D> >& cols,
           const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==simple);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    HTCoefficients<T, Basis> Au(u.dim(), u.basis());
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
            auto max    = maxintind(rows[idx[0]-1], u.basis());
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
                    auto i = maptoint(lambda, u.basis());
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
                        auto i = maptoint(lambda, u.basis());
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
                A(idx[0], 1).eval(input, output, "A");

                /* Save result */
                for (typename CoefficientsByLevel<T>::const_it it=
                     output.bylevel[0].map.begin();
                     it!=output.bylevel[0].map.end(); ++it) {
                        Index1D lambda(output.offset+1, (*it).first,XBSpline);
                        auto i = maptoint(lambda, u.basis());
                        UAu(i, colsU+k) = (*it).second;
                }
                for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                    if (output.bylevel[i].map.size()==0) break;
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[i].map.begin();
                         it!=output.bylevel[i].map.end(); ++it) {
                            Index1D lambda(output.offset+i,
                                           (*it).first,XWavelet);
                            auto j = maptoint(lambda, u.basis());
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
            const HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index1D> >& rows,
            const std::vector<IndexSet<Index1D> >& cols,
            const std::size_t hashtablelength)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.type()==laplace);
    assert(rows.size()==A.dim() && cols.size()==A.dim());

    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> > Matrix;

    using flens::_;

    HTCoefficients<T, Basis> Au(u.dim(), u.basis());
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
            auto max    = maxintind(rows[idx[0]-1], u.basis());
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
                    auto i = maptoint(lambda, u.basis());
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
                        auto i = maptoint(lambda, u.basis());
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
                        auto i = maptoint(lambda, u.basis());
                        UAu(i, colsU+k) = (*it).second;
                }
                for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
                    if (output.bylevel[i].map.size()==0) break;
                    for (typename CoefficientsByLevel<T>::const_it it=
                         output.bylevel[i].map.begin();
                         it!=output.bylevel[i].map.end(); ++it) {
                            Index1D lambda(output.offset+i,
                                           (*it).first,XWavelet);
                            auto j = maptoint(lambda, u.basis());
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
eval(Sepop<Optype>& A,
     const HTCoefficients<T, Basis>& u,
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


template <typename T, typename Basis>
FLENS_DEFAULT_INDEXTYPE
maxlevel(const HTCoefficients<T, Basis>& u)
{
    FLENS_DEFAULT_INDEXTYPE jmax = u.basis().j0;

    for (auto tit=u.tree().getGeneralTree().end();
              tit>=u.tree().getGeneralTree().begin(); tit--) {
        if (tit.getNode()->isLeaf()) {
            FLENS_DEFAULT_INDEXTYPE rows    = tit.getNode()->getContent()->getUorB().numRows();
            auto lambda = maptowav(rows, u.basis());
            FLENS_DEFAULT_INDEXTYPE j       = lambda.j;
            if (lambda.xtype==XWavelet) ++j;
            jmax = MAX(jmax, j);
        } else {
            break;
        }
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


template <typename T, typename Basis>
T
compOmegamax2(const HTCoefficients<T, Basis>& u, const T order)
{
    auto d = u.dim();

    FLENS_DEFAULT_INDEXTYPE jmax = maxlevel(u);

    return d*std::pow(2., 2.*order*jmax);
}


template <typename T, typename Basis>
T
compUnDistFac(const HTCoefficients<T, Basis>& u, const T order)
{
    auto d   = u.dim();
    FLENS_DEFAULT_INDEXTYPE j0   = u.basis().j0;
    FLENS_DEFAULT_INDEXTYPE jmax = maxlevel(u);

    T factor = std::pow(2., (order+1)*j0-jmax-1.);
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

    return sum/compOmegamin2(u.basis(), u.dim(), order);
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
     const HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >& cols,
     const double eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale(u, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin2(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis());

    unsigned counter = 0;
    for (FLENS_DEFAULT_INDEXTYPE i=-1*S.n(); i<=(signed) S.nplus(); ++i, ++counter) {
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
            sum.truncate(eps);
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            prod.truncate(eps);
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
HTCoefficients<T, Basis>
evalS2(Sepdiagscal<Basis>& S,
       const HTCoefficients<T, Basis>& u,
       const std::vector<IndexSet<Index1D> >& cols,
       const double eps)
{
    assert(S.dim()==(unsigned) u.dim());
    assert(cols.size()==S.dim());

    typedef typename Sepdiagscal<Basis>::size_type  size_type;

    T iscale = compIndexscale2(u, S.order());
    S.set_iscale(iscale);
    S.comp_n();
    T omega2 = compOmegamin4(S.basis(), S.dim(), S.order());

    T factor1 = S.h()*(1./std::sqrt(omega2));
    #ifdef DEBUG_CANCEL
        HTCoefficients<T, Basis> sumexact;
        T                        sumnorms = 0.;
    #endif

    HTCoefficients<T, Basis>    sum(u.dim(), u.basis());

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


        if (i==-1*(signed) S.n()) {
            prod.truncate(eps);
            sum = prod;
            #ifdef DEBUG_CANCEL
                sumexact  = prod;
                prod.orthogonalize();
                sumnorms += prod.tree().L2normorthogonal();
            #endif
        } else {
            prod.truncate(eps);
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
        auto max = maxintind(cols[j-1], S.basis());
        for (FLENS_DEFAULT_INDEXTYPE k=-1*S.n(); k<=(signed) S.nplus(); ++k) {
            DV vecA(max*max);

            T factor2 = 2.*std::pow(M_PI, -0.5)*
                        (1./(1+std::exp(-1.*S.h()*(T) k)));
            factor2   = std::pow(factor2*factor1, 1./(T) S.dim());
            T alpha   = std::pow(std::log(1.+std::exp((T) k*S.h())), 2.);
            for (FLENS_DEFAULT_INDEXTYPE i=1; i<=(signed) max; ++i) {
                auto index        = maptowav(i, S.basis());
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
contraction(const HTCoefficients<T, Basis>& u,
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
                FLENS_DEFAULT_INDEXTYPE rowi = maptoint(mu, u.basis());
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
contraction(const HTCoefficients<T, Basis>& u,
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
    frame = extract(f, idx);

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
            auto max = maxintind(activex[idx[0]-1], f.basis());
            assert(max<=(unsigned) U.numRows());

            Matrix copy = U;
            if (max!=(unsigned) U.numRows()) {
                U.resize(max, U.numCols());
            } else {
                U.fill((T) 0);
            }

            for (auto& lambda : activex[idx[0]-1]) {
                auto i = maptoint(lambda, f.basis());
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
            unsigned rowsnew = maxintind(activex, f.basis());

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
            auto max = maxintind(activex[idx[0]-1], f.basis());
            assert(max>=(unsigned) U.numRows());

            if (max!=(unsigned) U.numRows()) {
                Matrix copy = U;
                U.resize(max, U.numCols());
                U(_(1, copy.numRows()), _) = copy;
            }
        }
    }
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_TCC 1
