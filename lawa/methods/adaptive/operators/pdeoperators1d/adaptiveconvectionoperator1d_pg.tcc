namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
AdaptiveConvectionOperator1D_PG<T,TrialBasis,TestBasis>::
AdaptiveConvectionOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis)
: trialbasis1d(_trialbasis), testbasis1d(_testbasis),
  compression1d(_trialbasis), convection_op1d(_trialbasis,_testbasis), prec1d(),
  data(convection_op1d, prec1d, compression1d)
{}

template <typename T, typename TrialBasis, typename TestBasis>
T
AdaptiveConvectionOperator1D_PG<T,TrialBasis,TestBasis>::operator()(const Index1D &row_index,
                                                             const Index1D &col_index)
{
    return data(row_index, col_index);
}


template <typename T, typename TrialBasis, typename TestBasis>
T
AdaptiveConvectionOperator1D_PG<T,TrialBasis,TestBasis>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    return this->operator()(row_index,col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
void
AdaptiveConvectionOperator1D_PG<T,TrialBasis,TestBasis>::
clear()
{
    data.clear();
}

}   // namespace lawa
