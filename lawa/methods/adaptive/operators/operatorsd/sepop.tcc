#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC
#define LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC 1

#include <cassert>
#include <iostream>
#include <htucker/htucker.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

namespace lawa
{

template <typename Optype>
Sepop<Optype>::Sepop(const Opvec& _ops,
                     const size_type _rank, const size_type _dim):
    ops_(_ops),
    rank_(_rank),
    dim_(_dim),
    indexset_()
{
    assert(ops_.size()==rank()*dim() ||
           (ops_.size()==rank() && rank()==dim()));

    if (ops_.size()==rank()*dim()) {
        type_ = standard;
    } else {
        type_ = simple;
    }
}


template <typename Optype>
Sepop<Optype>::Sepop(const Optype& op,
                     const size_type _rank, const size_type _dim):
    ops_(1, &op),
    rank_(_rank),
    dim_(_dim),
    indexset_(),
    type_(laplace){}


template <typename Optype>
typename Sepop<Optype>::size_type
Sepop<Optype>::rank() const
{
    return rank_;
}


template <typename Optype>
typename Sepop<Optype>::size_type
Sepop<Optype>::dim() const
{
    return dim_;
}


template <typename Optype>
SepopType
Sepop<Optype>::type() const
{
    return type_;
}


template <typename Optype>
const IndexSet<Index1D>&
Sepop<Optype>::getIndexset() const
{
    return indexset_;
}


template <typename Optype>
void
Sepop<Optype>::setIndexset(const IndexSet<Index1D>& _indexset)
{
    indexset_ = _indexset;
}


template <typename Optype>
const typename Sepop<Optype>::Opvec&
Sepop<Optype>::ops() const
{
    return ops_;
}


template <typename Optype>
typename Sepop<Optype>::Opvec&
Sepop<Optype>::ops()
{
    return ops_;
}


template <typename Optype>
const Optype&
Sepop<Optype>::ops(const size_type i, const size_type j) const
{
    assert(type()==standard);
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return *ops()[(j-1)*rank()+(i-1)];
}


template <typename Optype>
Optype&
Sepop<Optype>::ops(const size_type i, const size_type j)
{
    assert(type()==standard);
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return *ops()[(j-1)*rank()+(i-1)];
}


template <typename Optype>
const Optype&
Sepop<Optype>::operator()(const size_type i, const size_type j) const
{
    assert(type()==standard);
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return ops(i, j);
}


template <typename Optype>
Optype&
Sepop<Optype>::operator()(const size_type i, const size_type j)
{
    assert(type()==standard);
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return ops(i, j);
}

/*
template <typename Optype>
template <typename T, typename Basis>
HTCoefficients<T, Basis>
Sepop<Optype>::eval(const HTCoefficients<T, Basis>& u,
                    const IndexSet<Index1D>& indexset,
                    const std::size_t hashtablelength)
{
    HTCoefficients<T, Basis> sum;
    if (!indexset.size()) {
        std::cerr << "Sepop<Optype>::eval: empty indexset\n";
        return sum;
    }

    for (size_type i=1; i<=rank(); ++i) {
        HTCoefficients<T, Basis> prod(u);
        std::cout << "Current prod is\n";
        prod.tree().print_w_UorB();
        for (size_type j=1; j<=dim(); ++j) {
            htucker::DimensionIndex idx(1);
            idx[0] = j;

            SepCoefficients<Lexicographical, T, Index1D>
            frame = extract(prod, idx);

            // Too many copies here -> efficiency loss!
            for (size_type k=1; k<=frame.rank(); ++k) {
                TreeCoefficients1D<T> input(hashtablelength,
                                            u.basis().j0);
                TreeCoefficients1D<T> output(hashtablelength,
                                            ops(i, j).getTestBasis().j0);
                input = frame(k, 1);
                Coefficients<Lexicographical, T, Index1D> temp;
                FillWithZeros(indexset, temp);
                output = temp;

                std::cout << "frame(" << k << ", 1)=\n" << frame(k, 1)
                          << "\ninput=\n" << input;
                std::cout << "\noutput=\n" << temp << std::endl;
                ops(i, j).eval(input, output, "A");
                fromTreeCoefficientsToCoefficients(output, temp);
                frame(k, 1) = temp;
            }
            set(prod, idx, frame);
        }

        if (i==1) {
            sum = prod;
        } else {
            sum.tree() = sum.tree()+prod.tree();
        }
    }

    return sum;
}
*/
} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC
