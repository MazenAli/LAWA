namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
WeightedPDEOperator1D_PG<T, TrialBasis, TestBasis>::
WeightedPDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
                      Function<T>& _reaction_f,
                      Function<T>& _convection_f,
                      Function<T>& _diffusion_f,
                      int order, bool _reactionIsZero,
                      bool _convectionIsZero, bool _diffusionIsZero)
    : trialbasis(_trialbasis), testbasis(_testbasis),
      reaction_f(_reaction_f), convection_f(_convection_f), diffusion_f(_diffusion_f),
      reactionIsZero(_reactionIsZero), convectionIsZero(_convectionIsZero), diffusionIsZero(_diffusionIsZero),
      reaction_integral(reaction_f, testbasis, trialbasis), convection_integral(convection_f, testbasis, trialbasis),
      diffusion_integral(diffusion_f, testbasis, trialbasis)

{
    reaction_integral.quadrature.setOrder(order);
    convection_integral.quadrature.setOrder(order);
    diffusion_integral.quadrature.setOrder(order);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
WeightedPDEOperator1D_PG<T, TrialBasis, TestBasis>::operator()(XType xtype1, int j1, long k1,
                                            		    	   XType xtype2, int j2, long k2) const
{
    // diffusion * v_x *  u_x + convection * v * u_x + reaction * v * u
    T val_reaction    = reactionIsZero   ? 0. : reaction_integral(  j1,k1,xtype1,0, j2,k2,xtype2,0);
    T val_convection  = convectionIsZero ? 0. : convection_integral(j1,k1,xtype1,0, j2,k2,xtype2,1);
    T val_diffusion   = diffusionIsZero  ? 0. : diffusion_integral( j1,k1,xtype1,1, j2,k2,xtype2,1);

    return val_reaction + val_convection + val_diffusion;

/*
    return    diffusion_integral(j1, k1, xtype1, 1, j2, k2, xtype2, 1)
            + convection_integral(j1, k1, xtype1, 0, j2, k2, xtype2, 1)
            + reaction_integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0);
*/
}

template <typename T, typename TrialBasis, typename TestBasis>
T
WeightedPDEOperator1D_PG<T, TrialBasis, TestBasis>::operator()(const Index1D &row_index, const Index1D &col_index)
const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                            col_index.xtype, col_index.j, col_index.k);
}

}   // namespace lawa

