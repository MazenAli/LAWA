namespace lawa {

template <typename T, typename Index, typename RHSType, typename ParamType>
AffineRhs<T,Index,RHSType,ParamType>::
AffineRhs(ThetaStructure<ParamType>& _thetas, std::vector<RHSType*>& _rhsvec)
 : FlexibleCompoundRhs<T,Index,RHSType>(_rhsvec), thetas(_thetas)
{
	assert(thetas.size() == this->rhsvec.size());
}

template <typename T, typename Index, typename RHSType, typename ParamType>
void
AffineRhs<T,Index,RHSType,ParamType>::
set_param(ParamType& mu)
{
	thetas.set_param(mu);
}


template <typename T, typename Index, typename RHSType, typename ParamType>
T
AffineRhs<T,Index,RHSType,ParamType>::
operator()(const Index &index)
{
    T val = 0.;

    for(auto& i : this->active_comp){

    //for(size_t i = 0; i < thetas.size(); ++i){
    	val += thetas.eval(i) * (*(this->rhsvec)[i])(index);
    }

    return val;
}

template <typename T, typename Index, typename RHSType, typename ParamType>
Coefficients<Lexicographical,T,Index>
AffineRhs<T,Index,RHSType,ParamType>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;

    T val;
    for(auto& index : indexset){
    	val = 0.;
        for(auto& i : this->active_comp){
        //for(size_t i = 0; i < thetas.size(); ++i){
        	val += thetas.eval(i) * (*(this->rhsvec)[i])(index);
        }
        r.insert(std::make_pair(index,val));
    }

    return r;
}

} // namespace lawa
