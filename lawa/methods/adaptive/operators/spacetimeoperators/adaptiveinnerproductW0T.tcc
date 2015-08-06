#include <iostream>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>

namespace lawa {

template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::
AdaptiveInnerProductW0T(const Basis2D& _basis, TrialPrec& _p_trial, TestPrec& _p_test)
    : basis(_basis),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression(_basis),
      P_trial_data(), P_test_data(), trialprec(_p_trial), testprec(_p_test), noprec(),
      op_identity_t(_basis.first), op_identity_x(_basis.second),
      op_laplace_t(_basis.first), op_laplace_x(_basis.second),
      data_identity_t(op_identity_t,     noprec, compression_1d_t),
      data_identity_x(op_identity_x,     noprec, compression_1d_x),
      data_laplace_t(op_laplace_t,       noprec, compression_1d_t),
      data_laplace_x(op_laplace_x,       noprec, compression_1d_x)
{
}

template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
T
AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::operator()(const Index2D &row_index, const Index2D &col_index){
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TestPrec>::value) {
        // Left precondioning:
        const_coeff_it it_row_index   = P_test_data.find(row_index);
        //  Entry has already been computed:
        if (it_row_index != P_test_data.end()) {
            prec *= (*it_row_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = testprec(row_index);
            P_test_data[row_index] = tmp;
            prec *= tmp;
        }
    }

    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TrialPrec>::value) {
        // Right precondioning:
        const_coeff_it it_col_index   = P_trial_data.find(col_index);
        //  Entry has already been computed:
        if (it_col_index != P_trial_data.end()) {
            prec *= (*it_col_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = trialprec(col_index);
            P_trial_data[col_index] = tmp;
            prec *= tmp;
        }
    }

    //data_identity_t(row_index.index1,col_index.index1) * data_convection_x(row_index.index2,col_index.index2);

    T id_x = data_identity_x(row_index.index2,col_index.index2);

    return prec * (data_identity_t(row_index.index1, col_index.index1) * (id_x + data_laplace_x(row_index.index2,col_index.index2))
    			+  data_laplace_t(row_index.index1, col_index.index1) * id_x * pow2i<T>(-row_index.index2.j - col_index.index2.j));
}

template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
Coefficients<Lexicographical,T,Index2D>
AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
 return lawa::mv(LambdaRow, *this, x);
}


template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
void
AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
					SparseMatrixT &A, T tol, bool useLinearIndex)
{
	std::cerr << "AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::"
			  << "toFlensSparseMatrix not implemented." << std::endl;
	assert(0);
	exit(1);
}

template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
void
AdaptiveInnerProductW0T<T, Basis2D, TrialPrec, TestPrec>::clear()
{
    data_identity_t.clear();
    data_identity_x.clear();
    data_laplace_t.clear();
    data_laplace_x.clear();
}

} // namespace lawa
