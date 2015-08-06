namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
TransposedPDEOperator1D_PG<T, TrialBasis, TestBasis>::
TransposedPDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
		T _reaction, T _convection, T _diffusion)
    : trialbasis(_trialbasis), testbasis(_testbasis), reaction(_reaction), convection(_convection), diffusion(_diffusion),
      integral(_testbasis, _trialbasis)
{
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedPDEOperator1D_PG<T, TrialBasis, TestBasis>::operator()(XType xtype1, int j1, long k1,
                                    XType xtype2, int j2, long k2) const
{   
    // diffusion * v_x *  u_x + convection * v * u_x + reaction * v * u
	// called as ()(index_u, index_v)
    return diffusion * integral(j2, k2, xtype2, 1, j1, k1, xtype1, 1)
            + convection * integral(j2, k2, xtype2, 1, j1, k1, xtype1, 0)
            + reaction * integral(j2, k2, xtype2, 0, j1, k1, xtype1, 0);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedPDEOperator1D_PG<T, TrialBasis, TestBasis>::operator()(const Index1D &row_index, const Index1D &col_index)
const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                            col_index.xtype, col_index.j, col_index.k);
}

}   //lawa

