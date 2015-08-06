namespace lawa {

template <typename Index, typename LocalOperatorType>
FlexibleBilformRhs<Index,LocalOperatorType>::FlexibleBilformRhs(std::vector<LocalOperatorType*>& _bilformvec)
 : bilformvec(_bilformvec), active_comp(bilformvec.size()), active_u(nullptr)
{
	for(size_t i=0; i < active_comp.size(); ++i){
		active_comp[i]=i;
	}
}

template <typename Index, typename LocalOperatorType>
typename LocalOperatorType::T
FlexibleBilformRhs<Index,LocalOperatorType>::
operator()(const Index &index)
{
    Coefficients<Lexicographical,T,Index> r;
    r[index] = 0.;

    for(auto& i : active_comp){
    	(*bilformvec[i]).eval(*active_u,r);
    }

    return r[index];
}

template <typename Index, typename LocalOperatorType>
Coefficients<Lexicographical,typename LocalOperatorType::T,Index>
FlexibleBilformRhs<Index,LocalOperatorType>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;
    FillWithZeros(indexset,r);

    for(auto& i : active_comp){
    	(*bilformvec[i]).eval(*active_u,r);
    }

    return r;
}

template <typename Index, typename LocalOperatorType>
void
FlexibleBilformRhs<Index,LocalOperatorType>::set_active_comp(int i)
{
	if(i < 0){
		active_comp.resize(bilformvec.size());
		for(size_t i=0; i < active_comp.size(); ++i){
			active_comp[i]=i;
		}
	}
	else{
		assert((size_t)i < bilformvec.size());
		active_comp.resize(1);
		active_comp[0] = i;
	}
}

template <typename Index, typename LocalOperatorType>
void
FlexibleBilformRhs<Index,LocalOperatorType>::
set_active_u(Coefficients<Lexicographical,T,Index> const* _u)
{
	active_u = _u;
}

template <typename Index, typename LocalOperatorType>
const std::vector<LocalOperatorType*>&
FlexibleBilformRhs<Index,LocalOperatorType>::
get_bilformvec() const
{
    return bilformvec;
}

template <typename Index, typename LocalOperatorType>
void
FlexibleBilformRhs<Index,LocalOperatorType>::
clear()
{
    for(auto& op : bilformvec){
        op->clear();
    }
}

} // namespace lawa
