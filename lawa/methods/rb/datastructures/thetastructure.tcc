namespace lawa {

template<typename ParamType>
ThetaStructure<ParamType>::ThetaStructure()
{}

template<typename ParamType>
ThetaStructure<ParamType>::ThetaStructure(const std::vector<ThetaFct>& _thetas)
 : thetas(_thetas)
{}


template<typename ParamType>
size_t
ThetaStructure<ParamType>::size() const
{
	return thetas.size();
}

template<typename ParamType>
void
ThetaStructure<ParamType>::set_param(const ParamType& _param)
{
	current_param = _param;
}


template<typename ParamType>
ParamType&
ThetaStructure<ParamType>::get_param()
{
	return current_param;
}

template<typename ParamType>
typename ParamType::value_type
ThetaStructure<ParamType>::eval(size_t i, const ParamType& mu) const
{
	return (thetas[i])(mu);
}

template<typename ParamType>
typename ParamType::value_type
ThetaStructure<ParamType>::eval(size_t i) const
{
	return eval(i,current_param);
}

} // namespace lawa
