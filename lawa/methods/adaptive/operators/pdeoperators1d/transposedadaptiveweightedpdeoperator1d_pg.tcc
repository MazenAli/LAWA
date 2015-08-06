namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
TransposedAdaptiveWeightedPDEOperator1D_PG(const TrialBasis& _trialbasis1d, const TestBasis& _testbasis1d,
 								Function<T> &_reaction_f,Function<T> &_convection_f, Function<T>& _diffusion_f,
 								int order, bool reactionIsZero, bool convectionIsZero, bool diffusionIsZero)
: trialbasis1d(_trialbasis1d), testbasis1d(_testbasis1d), compression1d(trialbasis1d),
  transposedweightedpdeop1d(trialbasis1d, testbasis1d,_reaction_f,_convection_f,_diffusion_f,order,reactionIsZero,convectionIsZero,diffusionIsZero),
  prec1d(), data(transposedweightedpdeop1d, prec1d, compression1d)
{

}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::operator()(const Index1D &row_index,
                                                              		 	 	   const Index1D &col_index)
{
    return data(row_index, col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
T
TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    return this->operator()(row_index,col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
void
TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
clear()
{
    data.clear();
}

}   //namespace lawa
