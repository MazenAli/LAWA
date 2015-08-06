namespace lawa {

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::ResidualRhs(LHSType& _lhs_bilform, RHSType& _rhs_fct)
 : lhs_bilform(_lhs_bilform), rhs_fct(_rhs_fct)
{}

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
typename LHSType::T
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::
operator()(const Index &index)
{
	return rhs_fct(index) - lhs_bilform(index);
}

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
Coefficients<Lexicographical,typename LHSType::T,Index>
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;
    FillWithZeros(indexset,r);

    r = rhs_fct(indexset);
    r -= lhs_bilform(indexset);

    return r;
}

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
void
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::
set_param(ParamType& mu)
{
	lhs_bilform.set_param(mu);
	rhs_fct.set_param(mu);

	lhs_bilform.set_active_comp(-1);
	rhs_fct.set_active_comp(-1);
}


template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
void
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::
set_active_u(DataType const* u)
{
	lhs_bilform.set_active_u(u);
}

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
void
ResidualRhs<Index,LHSType, RHSType, ParamType, DataType>::
clear()
{
    lhs_bilform.clear();
    rhs_fct.clear();
}

} // namespace lawa
