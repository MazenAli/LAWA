namespace lawa {

template <typename Index, typename LocalOperatorType, typename ParamType>
AffineLocalOperator<Index,LocalOperatorType,ParamType>::
AffineLocalOperator(ThetaStructure<ParamType>& _thetas, std::vector<LocalOperatorType*>& _localops)
 : FlexibleCompoundLocalOperator<Index,LocalOperatorType>(_localops), thetas(_thetas)
{
	assert(thetas.size() == this->localops.size());
}

template <typename Index, typename LocalOperatorType, typename ParamType>
void
AffineLocalOperator<Index,LocalOperatorType,ParamType>::set_param(ParamType& mu)
{
	thetas.set_param(mu);
}

template <typename Index, typename LocalOperatorType, typename ParamType>
void
AffineLocalOperator<Index,LocalOperatorType,ParamType>::eval(const Coefficients<Lexicographical,T,Index> &v,
		 	 	 	 	 	 	 	 	 	 	 	 	 	  Coefficients<Lexicographical,T,Index> &Av)
{
	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

    for(size_t i = 0; i < thetas.size(); ++i){
    	tmp.setToZero();
		this->localops[i]->eval(v, tmp);
		Av += thetas.eval(i) * tmp;
	}
}

template <typename Index, typename LocalOperatorType, typename ParamType>
void
AffineLocalOperator<Index,LocalOperatorType,ParamType>::eval(size_t i, const Coefficients<Lexicographical,T,Index> &v,
															Coefficients<Lexicographical,T,Index> &Av, bool eval_mu)
{
	assert(i < this->localops.size());

	if(eval_mu){
		Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);
		this->localops[i]->eval(v, tmp);
		Av += thetas.eval(i) * tmp;
	}
	else{
		this->localops[i]->eval(v, Av);
	}

}

template <typename Index, typename LocalOperatorType, typename ParamType>
template <typename Preconditioner>
void
AffineLocalOperator<Index,LocalOperatorType,ParamType>::eval(Coefficients<Lexicographical,T,Index> &v,
														Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P)
{
	for(auto& el : v){
		el.second *= P(el.first);
	}

	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

    for(size_t i = 0; i < thetas.size(); ++i){
    	tmp.setToZero();
    	this->localops[i]->eval(v, tmp);
		Av += thetas.eval(i) * tmp;
	}

	for(auto& el : Av){
		el.second *= P(el.first);
	}
	for(auto& el : v){
		el.second /= P(el.first);
	}
}

template <typename Index, typename LocalOperatorType, typename ParamType>
template <typename RightPrec, typename LeftPrec>
void
AffineLocalOperator<Index,LocalOperatorType,ParamType>::eval(Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP)
 {
	for(auto& el : v){
		el.second *= rightP(el.first);
	}

	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

	for(size_t i = 0; i < thetas.size(); ++i){
		tmp.setToZero();
		this->localops[i]->eval(v, tmp);
		Av += thetas.eval(i) * tmp;
	}

	for(auto& el : Av){
		el.second *= leftP(el.first);
	}
	for(auto& el : v){
		el.second /= rightP(el.first);
	}
 }

} // namespace lawa
