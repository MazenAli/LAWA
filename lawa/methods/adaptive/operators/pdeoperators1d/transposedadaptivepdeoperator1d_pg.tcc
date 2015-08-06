namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
TransposedAdaptivePDEOperator1D_PG<T,TrialBasis,TestBasis>::
TransposedAdaptivePDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
		  	  	  	  	  T _reaction, T _convection, T _diffusion)
: trialbasis1d(_trialbasis), testbasis1d(_testbasis),
  compression1d(_trialbasis),
  pde_op1d(_trialbasis, _testbasis, _reaction, _convection, _diffusion),
  data(pde_op1d, prec1d, compression1d)
{

}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedAdaptivePDEOperator1D_PG<T,TrialBasis,TestBasis>::operator()(const Index1D &row_index,
                                                           const Index1D &col_index)
{
	return data(row_index, col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedAdaptivePDEOperator1D_PG<T,TrialBasis,TestBasis>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    return this->operator()(row_index,col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
void
TransposedAdaptivePDEOperator1D_PG<T,TrialBasis,TestBasis>::toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow,
                                                                    const IndexSet<Index1D>& LambdaCol,
                                                                    SparseMatrixT &A_flens, int /*J*/,
                                                                    bool /*useLinearIndex*/)
{
    lawa::toFlensSparseMatrix(*this,LambdaRow,LambdaCol,A_flens);
}

template <typename T, typename TrialBasis, typename TestBasis>
void
TransposedAdaptivePDEOperator1D_PG<T,TrialBasis,TestBasis>::
clear()
{
    data.clear();
}

}   // namespace lawa
