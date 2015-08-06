namespace lawa {

template <typename Index, typename LocalOperatorType, typename ParamType>
AffineBilformRhs<Index, LocalOperatorType, ParamType>::
AffineBilformRhs(ThetaStructure<ParamType>& _thetas, std::vector<LocalOperatorType*>& _bilformvec)
 : FlexibleBilformRhs<Index,LocalOperatorType>(_bilformvec), thetas(_thetas)
{
	assert(thetas.size() == this->bilformvec.size());
}

template <typename Index, typename LocalOperatorType, typename ParamType>
void
AffineBilformRhs<Index, LocalOperatorType, ParamType>::
set_param(ParamType& mu)
{
	thetas.set_param(mu);
}


template <typename Index, typename LocalOperatorType, typename ParamType>
typename LocalOperatorType::T
AffineBilformRhs<Index, LocalOperatorType, ParamType>::
operator()(const Index &index)
{
    T val = 0.;

    for(auto& i : this->active_comp){
        Coefficients<Lexicographical,T,Index> r;
        r[index] = 0.;
    	(*(this->bilformvec)[i]).eval(*(this->active_u),r);
    	val += thetas.eval(i) * r[index];
    }

    return val;
}

template <typename Index, typename LocalOperatorType, typename ParamType>
Coefficients<Lexicographical,typename LocalOperatorType::T,Index>
AffineBilformRhs<Index, LocalOperatorType, ParamType>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;
    FillWithZeros(indexset,r);

    for(auto& i : this->active_comp){
        Coefficients<Lexicographical,T,Index> tmp;
        FillWithZeros(indexset,tmp);
    	(*(this->bilformvec)[i]).eval(*(this->active_u),tmp);
    	r += thetas.eval(i)*tmp;
    }

    return r;
}

} // namespace lawa
