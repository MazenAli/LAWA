namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
AdaptiveWeightedPDEOperator1D_PG(const TrialBasis& _trialbasis1d, const TestBasis& _testbasis1d,
 								Function<T> &_reaction_f,Function<T> &_convection_f, Function<T>& _diffusion_f,
 								int order, bool reactionIsZero, bool convectionIsZero, bool diffusionIsZero)
: trialbasis1d(_trialbasis1d), testbasis1d(_testbasis1d), compression1d(trialbasis1d),
  weightedpdeop1d(trialbasis1d, testbasis1d,_reaction_f,_convection_f,_diffusion_f,order,reactionIsZero,convectionIsZero,diffusionIsZero),
  prec1d()//,data(weightedpdeop1d, prec1d, compression1d)
{

}

template <typename T, typename TrialBasis, typename TestBasis>
T
AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::operator()(const Index1D &row_index,
                                                              		 const Index1D &col_index)
{

	//return data(row_index, col_index);

    Entry<Index1D> entry(row_index,col_index);

    typename DataWeightedPDEOp1D::const_iterator it_entry = data.find(entry);
    typename DataWeightedPDEOp1D::const_iterator it_end = data.end();

    typedef typename DataWeightedPDEOp1D::value_type val_type;

    if (it_entry != it_end) {
        return (*it_entry).second;
    }
    else {
        T val = weightedpdeop1d(row_index,col_index);

        if (fabs(val) > 0.) data.insert(val_type(entry,val));
        return val;
    }

}

template <typename T, typename TrialBasis, typename TestBasis>
T
AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col)
{
    Index1D row_index(j_row,k_row,xtype_row);
    Index1D col_index(j_col,k_col,xtype_col);
    return this->operator()(row_index,col_index);
}

template <typename T, typename TrialBasis, typename TestBasis>
void
AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis,TestBasis>::
clear()
{
    data.clear();
}

}   //namespace lawa
