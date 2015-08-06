namespace lawa {

template <typename T, typename Index, typename RHSType>
FlexibleCompoundRhs<T, Index,RHSType>::FlexibleCompoundRhs(std::vector<RHSType*>& _rhsvec)
 : rhsvec(_rhsvec), active_comp(rhsvec.size())
{
	for(size_t i=0; i < active_comp.size(); ++i){
		active_comp[i]=i;
	}
}

template <typename T, typename Index, typename RHSType>
T
FlexibleCompoundRhs<T, Index,RHSType>::operator()(const Index &index)
{
    T val = 0.;

    for(auto& i : active_comp){
    	val += (*rhsvec[i])(index);
    }

    return val;
}

template <typename T, typename Index, typename RHSType>
Coefficients<Lexicographical,T,Index>
FlexibleCompoundRhs<T,Index,RHSType>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;

    T val;
    for(auto& index : indexset){
    	val = 0.;
        for(auto& i : active_comp){
        	val += (*rhsvec[i])(index);
        }
        r.insert(std::make_pair(index,val));
    }

    return r;
}

template <typename T, typename Index, typename RHSType>
template <typename Prec>
Coefficients<Lexicographical,T,Index>
FlexibleCompoundRhs<T,Index,RHSType>::
operator()(const IndexSet<Index> &indexset, Prec& P)
{
    Coefficients<Lexicographical,T,Index> r;

    T val;
    for(auto& index : indexset){
    	val = 0.;
        for(auto& i : active_comp){
        	val += (*rhsvec[i])(index);
        }
        val *= P(index);
        r.insert(std::make_pair(index,val));
    }

    return r;
}

template <typename T, typename Index, typename RHSType>
void
FlexibleCompoundRhs<T, Index,RHSType>::set_active_comp(int i)
{
	if(i < 0){
		active_comp.resize(rhsvec.size());
		for(size_t i=0; i < active_comp.size(); ++i){
			active_comp[i]=i;
		}
	}
	else{
		assert((size_t)i < rhsvec.size());
		active_comp.resize(1);
		active_comp[0] = i;
	}
}

template <typename T, typename Index, typename RHSType>
void
FlexibleCompoundRhs<T, Index,RHSType>::
clear()
{
    for(auto& rhs : rhsvec){
        rhs->clear();
    }
}

}   // namespace lawa

