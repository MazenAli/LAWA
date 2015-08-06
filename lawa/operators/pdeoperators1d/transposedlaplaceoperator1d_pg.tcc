namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
TransposedLaplaceOperator1D_PG<T,TrialBasis, TestBasis>::TransposedLaplaceOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis)
  : trialbasis(_trialbasis), testbasis(_testbasis), integral(_testbasis, _trialbasis)
{
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedLaplaceOperator1D_PG<T,TrialBasis, TestBasis>::operator()(XType xtype1, int j1, long k1,
                                                          XType xtype2, int j2, long k2) const
{   
    // v_x * u_x, , aber Aufruf ist ()(index_u, index_v)
    return integral(j2, k2, xtype2, 1, j1, k1, xtype1, 1);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedLaplaceOperator1D_PG<T, TrialBasis, TestBasis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                                                                      col_index.xtype, col_index.j, col_index.k);
}

}    //namespace lawa
