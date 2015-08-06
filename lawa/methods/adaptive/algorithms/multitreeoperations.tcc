namespace lawa {

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index1D>  &v,
                Coefficients<Lexicographical,T,Index1D>  &C_v, const char* residualType,
                bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index;
        C_index = C((*it).first, (T)1., basis);
        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it newindex=C_index.begin(); newindex!=C_index.end(); ++newindex) {
                if (C_v.find(*newindex)==C_v.end()) {
                    completeMultiTree(basis,*newindex,C_v, sparsetree);
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }
    return;
}

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v, const char* residualType,
                bool IsMW, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index2D>::iterator                                      set2d_it;


#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif BOOST
    typedef boost::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif CONEONE
    typedef std::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer_x(v.size());
    IndexConeContainer indexConeContainer_y(v.size());

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2;

        //C_index1 = C((*it).first.index1, (T)1., basis.first);
        //C_index2 = C((*it).first.index2, (T)1., basis.second);

        IndexConeContainer::const_iterator indexConeContainer_ptr;
        indexConeContainer_ptr = indexConeContainer_x.find((*it).first.index1);
        if (indexConeContainer_ptr!=indexConeContainer_x.end()) {
            C_index1 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            indexConeContainer_x[(*it).first.index1] = C_index1;
        }
        indexConeContainer_ptr = indexConeContainer_y.find((*it).first.index2);
        if (indexConeContainer_ptr!=indexConeContainer_y.end()) {
            C_index2 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            indexConeContainer_y[(*it).first.index2] = C_index2;
        }

        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index2D newindex((*it_C_index1), (*it).first.index2);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,1,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D newindex((*it).first.index1, (*it_C_index2));
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,2,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
        else if (strcmp(residualType,"large1")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    Index2D newindex((*it_C_index1), (*it_C_index2));
                    if (C_v.find(newindex)==C_v.end()) {
                        completeMultiTree(basis,newindex,C_v,0,sparsetree);
                    }
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }


    return;
}


template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v, IndexSet<Index2D>& Cdiff_v,
                const char* residualType, bool IsMW, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index2D>::iterator                                      set2d_it;


#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif BOOST
    typedef boost::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif CONEONE
    typedef std::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer_x(v.size());
    IndexConeContainer indexConeContainer_y(v.size());

    IndexSet<Index2D>::iterator it_indexset2d;

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2;

        //C_index1 = C((*it).first.index1, (T)1., basis.first);
        //C_index2 = C((*it).first.index2, (T)1., basis.second);

        IndexConeContainer::const_iterator indexConeContainer_ptr;
        indexConeContainer_ptr = indexConeContainer_x.find((*it).first.index1);
        if (indexConeContainer_ptr!=indexConeContainer_x.end()) {
            C_index1 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            indexConeContainer_x[(*it).first.index1] = C_index1;
        }
        indexConeContainer_ptr = indexConeContainer_y.find((*it).first.index2);
        if (indexConeContainer_ptr!=indexConeContainer_y.end()) {
            C_index2 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            indexConeContainer_y[(*it).first.index2] = C_index2;
        }

        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index2D newindex((*it_C_index1), (*it).first.index2);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,Cdiff_v,1,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,Cdiff_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D newindex((*it).first.index1, (*it_C_index2));
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,Cdiff_v,2,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,Cdiff_v,0,sparsetree);
                }
            }
        }
        else if (strcmp(residualType,"large1")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    Index2D newindex((*it_C_index1), (*it_C_index2));
                    if (C_v.find(newindex)==C_v.end()) {
                        completeMultiTree(basis,newindex,C_v,Cdiff_v,0,sparsetree);
                    }
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }

    }

    // Because completeMultiTree only checks if an index is already in Cv, not in v,
    // we have to remove a possible entry v_i from Cdiffv
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
    	it_indexset2d = Cdiff_v.find((*it).first);
    	if(it_indexset2d != Cdiff_v.end()){
    		Cdiff_v.erase(it_indexset2d);
    	}
    }

    return;
}


template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                Coefficients<Lexicographical,T,Index3D>  &C_v, const char* residualType,
                bool IsMW, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index3D>::const_iterator                                set3d_it;

#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif BOOST
    typedef boost::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#elif CONEONE
    typedef std::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer_x(v.size());
    IndexConeContainer indexConeContainer_y(v.size());
    IndexConeContainer indexConeContainer_z(v.size());

    for (const_coeff3d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2, C_index3;
        IndexConeContainer::const_iterator indexConeContainer_ptr;
        indexConeContainer_ptr = indexConeContainer_x.find((*it).first.index1);
        if (indexConeContainer_ptr!=indexConeContainer_x.end()) {
            C_index1 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            indexConeContainer_x[(*it).first.index1] = C_index1;
        }
        indexConeContainer_ptr = indexConeContainer_y.find((*it).first.index2);
        if (indexConeContainer_ptr!=indexConeContainer_y.end()) {
            C_index2 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            indexConeContainer_y[(*it).first.index2] = C_index2;
        }
        indexConeContainer_ptr = indexConeContainer_z.find((*it).first.index3);
        if (indexConeContainer_ptr!=indexConeContainer_z.end()) {
            C_index3 = (*indexConeContainer_ptr).second;
        }
        else {
            C_index3 = C((*it).first.index3, (T)1., basis.second);
            indexConeContainer_z[(*it).first.index3] = C_index3;
        }
        if (strcmp(residualType,"standard")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index3D newindex((*it_C_index1), (*it).first.index2, (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,1,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index3D newindex((*it).first.index1, (*it_C_index2), (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,2,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
            for (const_set1d_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
                Index3D newindex((*it).first.index1, (*it).first.index2, (*it_C_index3));
                if (C_v.find(newindex)==C_v.end()) {
                    if (IsMW) completeMultiTree(basis,newindex,C_v,3,sparsetree);
                    else      completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
        else if (strcmp(residualType,"large1")==0) {
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                    for (const_set1d_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
                        Index3D newindex((*it_C_index1), (*it_C_index2), (*it_C_index3));
                        if (C_v.find(newindex)==C_v.end()) {
                            completeMultiTree(basis,newindex,C_v,0,sparsetree);
                        }
                    }
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: unknown residual type " << residualType << std::endl;
            exit(1);
        }
    }
    return;
}

/*
template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                Coefficients<Lexicographical,T,Index3D>  &C_v, int coordDirec)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;
    typedef IndexSet<Index3D>::const_iterator                                set3d_it;

#ifdef TRONE
    typedef std::tr1::unordered_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#else
    typedef __gnu_cxx::hash_map<Index1D, IndexSet<Index1D>, index_hashfunction<Index1D>, index_eqfunction<Index1D> >
            IndexConeContainer;
#endif

    IndexConeContainer indexConeContainer(v.size());

    for (const_coeff3d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        if (coordDirec == 1) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index1);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index1, (T)1., basis.first);
                indexConeContainer[(*it).first.index1] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it_C_index), (*it).first.index2, (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else if (coordDirec == 2) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index2);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index2, (T)1., basis.second);
                indexConeContainer[(*it).first.index2] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it).first.index1, (*it_C_index), (*it).first.index3);
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else if (coordDirec == 3) {
            IndexSet<Index1D > C_index;
            IndexConeContainer::const_iterator indexConeContainer_ptr;
            indexConeContainer_ptr = indexConeContainer.find((*it).first.index3);
            if (indexConeContainer_ptr!=indexConeContainer.end()) {
                C_index = (*indexConeContainer_ptr).second;
            }
            else {
                C_index = C((*it).first.index3, (T)1., basis.second);
                indexConeContainer[(*it).first.index3] = C_index;
            }
            for (const_set1d_it it_C_index=C_index.begin(); it_C_index!=C_index.end(); ++it_C_index) {
                Index3D newindex((*it).first.index1, (*it).first.index2, (*it_C_index));
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,coordDirec);
                }
            }
        }
        else {
            std::cerr << "extendMultiTree: non-admissible coordinate direction " << coordDirec
                      << std::endl;
            exit(1);
        }
    }
    return;
}
*/

template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                          Coefficients<Lexicographical,T,Index2D>  &C_v, int J, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    IndexSet<Index1D> LambdaB_x1;
    for (int k= basis.first.mra.rangeI(basis.first.j0).firstIndex();
             k<=basis.first.mra.rangeI(basis.first.j0).lastIndex(); ++k) {
        LambdaB_x1.insert(Index1D(basis.first.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.first.rangeJL(j).firstIndex(); k<=basis.first.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.first.rangeJR(j).firstIndex(); k<=basis.first.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x2;
    for (int k= basis.second.mra.rangeI(basis.second.j0).firstIndex();
             k<=basis.second.mra.rangeI(basis.second.j0).lastIndex(); ++k) {
        LambdaB_x2.insert(Index1D(basis.second.j0,k,XBSpline));
    }
    for (int j=basis.second.j0; j<=J; ++j) {
        for (int k=basis.second.rangeJL(j).firstIndex(); k<=basis.second.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.second.rangeJR(j).firstIndex(); k<=basis.second.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
    }

    for (const_set1d_it it_x1=LambdaB_x1.begin(); it_x1!=LambdaB_x1.end(); ++it_x1) {
        for (const_set1d_it it_x2=LambdaB_x2.begin(); it_x2!=LambdaB_x2.end(); ++it_x2) {
            Index2D newindex((*it_x1), (*it_x2));
            if (C_v.find(newindex)==C_v.end()) {
                completeMultiTree(basis,newindex,C_v,0,sparsetree);
            }
        }
    }
}

template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                          Coefficients<Lexicographical,T,Index3D>  &C_v, int J, bool sparsetree)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff3d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    IndexSet<Index1D> LambdaB_x1;
    for (int k= basis.first.mra.rangeI(basis.first.j0).firstIndex();
             k<=basis.first.mra.rangeI(basis.first.j0).lastIndex(); ++k) {
        LambdaB_x1.insert(Index1D(basis.first.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.first.rangeJL(j).firstIndex(); k<=basis.first.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.first.rangeJR(j).firstIndex(); k<=basis.first.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x1.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x2;
    for (int k= basis.second.mra.rangeI(basis.second.j0).firstIndex();
             k<=basis.second.mra.rangeI(basis.second.j0).lastIndex(); ++k) {
        LambdaB_x2.insert(Index1D(basis.second.j0,k,XBSpline));
    }
    for (int j=basis.second.j0; j<=J; ++j) {
        for (int k=basis.second.rangeJL(j).firstIndex(); k<=basis.second.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.second.rangeJR(j).firstIndex(); k<=basis.second.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x2.insert(Index1D(j,k,XWavelet));
        }
    }
    IndexSet<Index1D> LambdaB_x3;
    for (int k= basis.third.mra.rangeI(basis.third.j0).firstIndex();
             k<=basis.third.mra.rangeI(basis.third.j0).lastIndex(); ++k) {
        LambdaB_x3.insert(Index1D(basis.third.j0,k,XBSpline));
    }
    for (int j=basis.first.j0; j<=J; ++j) {
        for (int k=basis.third.rangeJL(j).firstIndex(); k<=basis.third.rangeJL(j).lastIndex(); ++k) {
            LambdaB_x3.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.third.rangeJR(j).firstIndex(); k<=basis.third.rangeJR(j).lastIndex(); ++k) {
            LambdaB_x3.insert(Index1D(j,k,XWavelet));
        }
    }
    for (const_set1d_it it_x1=LambdaB_x1.begin(); it_x1!=LambdaB_x1.end(); ++it_x1) {
        for (const_set1d_it it_x2=LambdaB_x2.begin(); it_x2!=LambdaB_x2.end(); ++it_x2) {
            for (const_set1d_it it_x3=LambdaB_x3.begin(); it_x3!=LambdaB_x3.end(); ++it_x3) {
                Index3D newindex((*it_x1), (*it_x2), (*it_x3));
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v,0,sparsetree);
                }
            }
        }
    }
}

/*
// Non-Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end())  return;
    else                            v[index1d] = 0.;

    int  j = index1d.j;
    long k = index1d.k;

    Support<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                if (overlap(supp,new_supp)>0) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                }
            }
        }
        else {
            bool foundPredecessor = false;
            // First, we check whether there is a predecessor.
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)!=v.end()) {
                        foundPredecessor = true;
                        break;
                    }
                }
            }
            // Second, in case there exists no predecessor in the tree, we add the first candidate
            // to the tree.
            if (!foundPredecessor) {
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                        Index1D new_index1d(Index1D(new_j,new_k,new_type));
                        completeMultiTree(basis,new_index1d,v,sparsetree);
                        break;
                    }
                }
            }
        }
    }
    return;
}
*/

// Non-Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end())  return;
    else                            v[index1d] = 0.;

    int  j = index1d.j;
    long k = index1d.k;

    Support<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                if (overlap(supp,new_supp)>0) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                }
            }
        }
        else {
            bool foundPredecessor = false;
            bool foundCoveredSupp = false;
            long covering_k = new_k_last;
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    Index1D new_index1d(new_j,new_k,new_type);
                    if (v.find(new_index1d)!=v.end()) {
                        foundPredecessor = true;
                        break;
                    }
                	foundCoveredSupp = true;
                	covering_k = std::min(covering_k, new_k);
                }
            }
            if (!foundPredecessor) {
            	if(foundCoveredSupp){
                    Index1D new_index1d(new_j,covering_k,new_type);
                    completeMultiTree(basis,new_index1d,v,sparsetree);
            	}
            	else{
               //     std::cerr << "¤¤¤ ==== > Support to be covered: " << supp << std::endl;
                    for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                        Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                        if (new_supp.l1 > supp.l1) {
                //            std::cerr << "¤¤¤       Covering with: " << basis.generator(new_type).support(new_j,new_k-1) << std::endl;
               //             std::cerr << "¤¤¤       Covering with: " << new_supp << std::endl;
                            completeMultiTree(basis, Index1D(new_j,new_k - 1,new_type),v,sparsetree);
                            completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
                            while(new_supp.l2 < supp.l2){
                            	new_k++;
                            	new_supp = basis.generator(new_type).support(new_j,new_k);
                 //               std::cerr << "¤¤¤       Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
                            }
                 //           std::cerr << "¤¤¤ <==============   " << std::endl;
                            break;
                        }
                    }
            	}
            }
        }
    }
    return;
}

// Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end())  return;
    else                            v[index1d] = 0.;

    int  j = index1d.j;
    long k = index1d.k;

    PeriodicSupport<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
        	if(new_k_first < new_k_last){
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                    }
                }
        	}
        	else{
        		long lastIndex = 0;
        		long firstIndex = 0;
        		if(new_type == XBSpline){
        			firstIndex = basis.mra.rangeI(new_j).firstIndex();
        			lastIndex = basis.mra.rangeI(new_j).lastIndex();
        		}
        		else{
        			firstIndex = basis.rangeJ(new_j).firstIndex();
        			lastIndex = basis.rangeJ(new_j).lastIndex();
        		}
                for (long new_k=new_k_first; new_k <= lastIndex ; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                    }
                }
                for (long new_k=firstIndex; new_k <=new_k_last ; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v);
                    }
                }
        	}
        }
        else {
            bool foundPredecessor = false;
            bool foundCoveredSupp = false;
            long covering_k = new_k_last;

        	if(new_k_first < new_k_last){
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    //if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    if(minimal_overlap(covered_supp, supp) >= supp.length()){
                    	Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                        foundCoveredSupp = true;
                        covering_k = std::min(covering_k, new_k);
                    }
                }
        	}
        	else{
        		long lastIndex = 0;
        		long firstIndex = 0;
        		if(new_type == XBSpline){
        			firstIndex = basis.mra.rangeI(new_j).firstIndex();
        			lastIndex = basis.mra.rangeI(new_j).lastIndex();
        		}
        		else{
        			firstIndex = basis.rangeJ(new_j).firstIndex();
        			lastIndex = basis.rangeJ(new_j).lastIndex();
        		}

            	for (long new_k=firstIndex; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    if(minimal_overlap(covered_supp, supp) >= supp.length()){
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                        foundCoveredSupp = true;
                        covering_k = std::min(covering_k, new_k);
                    }
        		}
                if(!foundPredecessor){
                	if(!foundCoveredSupp){ covering_k = lastIndex;}
            		for(long new_k = new_k_first; new_k <= lastIndex; ++new_k){
                        PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                        if(minimal_overlap(covered_supp, supp) >= supp.length()){
                            Index1D new_index1d(new_j,new_k,new_type);
                            if (v.find(new_index1d)!=v.end()) {
                                foundPredecessor = true;
                                break;
                            }
                            foundCoveredSupp = true;
                            covering_k = std::min(covering_k, new_k);
                        }
                    }
                }
        	}

            if (!foundPredecessor) {
            	if(foundCoveredSupp){
                    Index1D new_index1d(new_j,covering_k,new_type);
                    completeMultiTree(basis,new_index1d,v,sparsetree);
            	}
            	else{ // we have to find more than 1 bf covering the complete support
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type == XBSpline){
						firstIndex = basis.mra.rangeI(new_j).firstIndex();
						lastIndex = basis.mra.rangeI(new_j).lastIndex();
					}
					else{
						firstIndex = basis.rangeJ(new_j).firstIndex();
						lastIndex = basis.rangeJ(new_j).lastIndex();
					}
            		if(new_k_first < new_k_last){
						for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
							PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
							// left boundary of periodic support: max(l1, li2)
	                        if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
	                            T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
	                            bool gap = supp.gaplength() > 0 ? true : false;
	                            bool wrapped = false;
	                            if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
	                        		long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
		                            completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,sparsetree);
		                            if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
		                            	wrapped = true;
		                            }
	                        	}
	                            completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
	                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
	                            	wrapped = true;
	                            }

	                            while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
	                            	new_k++;
									if(new_k > lastIndex) new_k = firstIndex;
	                            	new_supp = basis.generator(new_type).support(new_j,new_k);
	                                completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
		                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
		                            	wrapped = true;
		                            }
	                            }
	                            break;
	                        }
						}
					}
					else{
						bool is_break = false;
						for (long new_k=new_k_first; new_k <= lastIndex ; ++new_k) {
							PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
	                        if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
	                            T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
	                            bool gap = supp.gaplength() > 0 ? true : false;
	                            bool wrapped = false;
	                            if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
	                        		long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
		                            completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,sparsetree);
		                            if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
		                            	wrapped = true;
		                            }
	                        	}
	                            completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
	                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
	                            	wrapped = true;
	                            }

	                            while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
									new_k++;
									if(new_k > lastIndex){
										is_break = false;
										break;
									}
	                            	new_supp = basis.generator(new_type).support(new_j,new_k);
	                                completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
		                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
		                            	wrapped = true;
		                            }
		                            is_break = true;
	                            }
	                            break;
	                        }
						}
						if(is_break == false){
							for (long new_k=firstIndex; new_k <=new_k_last ; ++new_k) {
								PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
		                        if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
		                            T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
		                            bool gap = supp.gaplength() > 0 ? true : false;
		                            bool wrapped = false;
		                            if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
		                        		long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
			                            completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,sparsetree);
			                            if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
			                            	wrapped = true;
			                            }
		                        	}
		                            completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
		                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
		                            	wrapped = true;
		                            }

		                            while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
										new_k++;
										if(new_k > lastIndex) new_k = firstIndex;
		                            	new_supp = basis.generator(new_type).support(new_j,new_k);
		                                completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,sparsetree);
			                            if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
			                            	wrapped = true;
			                            }
		                            }
		                            break;
		                        }
							}
						}
					}
            	}
            }
        }
    }
    return;
}


// ---- Non-Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v,
                  IndexSet<Index1D>& diff_v, bool sparsetree=false)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index1d);
    	v[index1d] = 0.;
    }

    int  j = index1d.j;
    long k = index1d.k;

    Support<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                if (overlap(supp,new_supp)>0) {
                    Index1D new_index1d(Index1D(new_j,new_k,new_type));
                    if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v, diff_v);
                }
            }
        }
        else {
            bool foundPredecessor = false;
            bool foundCoveredSupp = false;
            long covering_k = new_k_last;
            for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                Support<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    Index1D new_index1d(new_j,new_k,new_type);
                    if (v.find(new_index1d)!=v.end()) {
                        foundPredecessor = true;
                        break;
                    }
                	foundCoveredSupp = true;
                	covering_k = std::min(covering_k, new_k);
                }
            }
            if (!foundPredecessor) {
            	if(foundCoveredSupp){
                    Index1D new_index1d(new_j,covering_k,new_type);
                    completeMultiTree(basis,new_index1d,v,diff_v,sparsetree);
            	}
            	else{
                    for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                 //       std::cerr << "¤¤¤ ==== > Support to be covered: " << supp << std::endl;
                        Support<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                        if (new_supp.l1 > supp.l1) {
                //            std::cerr << "¤¤¤       Covering with: " << basis.generator(new_type).support(new_j,new_k-1) << std::endl;
               //             std::cerr << "¤¤¤       Covering with: " << new_supp << std::endl;
                            completeMultiTree(basis, Index1D(new_j,new_k - 1,new_type),v,diff_v,sparsetree);
                            completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
                            while(new_supp.l2 < supp.l2){
                            	new_k++;
                            	new_supp = basis.generator(new_type).support(new_j,new_k);
                //                std::cerr << "¤¤¤       Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
                            }
               //             std::cerr << "¤¤¤ <==============   " << std::endl;
                            break;
                        }
                    }
            	}
            }
        }
    }
    return;
}

// ---- Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v,
                  IndexSet<Index1D>& diff_v, bool sparsetree=false)
{
    int j0 = basis.j0;

    if (v.find(index1d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index1d);
    	v[index1d] = 0.;
    }

    int  j = index1d.j;
    long k = index1d.k;

    PeriodicSupport<typename Basis::T> supp = basis.generator(index1d.xtype).support(j,k);

    int new_j = 0;
    long new_k_first = 0, new_k_last = 0;
    bool checkPredecessors=true;
    XType new_type = XWavelet;
    if (j==j0 && index1d.xtype==XWavelet) {
        basis.getScalingNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XBSpline;
        assert(new_j==j);
    }
    else if (j>j0 && index1d.xtype==XWavelet) {
        basis.getLowerWaveletNeighborsForWavelet(j,k,basis,new_j,new_k_first,new_k_last);
        new_type = XWavelet;
        assert(new_j==j-1);
    }
    else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

    if (checkPredecessors) {
        if (!sparsetree) {
        	if(new_k_first < new_k_last){
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v, diff_v);
                    }
                }
        	}
        	else{
        		long lastIndex = 0;
        		long firstIndex = 0;
        		if(new_type == XBSpline){
        			firstIndex = basis.mra.rangeI(new_j).firstIndex();
        			lastIndex = basis.mra.rangeI(new_j).lastIndex();
        		}
        		else{
        			firstIndex = basis.rangeJ(new_j).firstIndex();
        			lastIndex = basis.rangeJ(new_j).lastIndex();
        		}
                for (long new_k=new_k_first; new_k <= lastIndex ; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v,diff_v);
                    }
                }
                for (long new_k=firstIndex; new_k <=new_k_last ; ++new_k) {
                    PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
                    if (overlap(supp,new_supp)>0) {
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)==v.end()) completeMultiTree(basis,new_index1d,v,diff_v);
                    }
                }
        	}
        }
        else {
            bool foundPredecessor = false;
            bool foundCoveredSupp = false;
            long covering_k = new_k_last;

        	if(new_k_first < new_k_last){
                for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    //if (covered_supp.l1<=supp.l1 && covered_supp.l2>=supp.l2) {
                    if(minimal_overlap(covered_supp, supp) >= supp.length()){
                    	Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                        foundCoveredSupp = true;
                        covering_k = std::min(covering_k, new_k);
                    }
                }
        	}
        	else{
        		long lastIndex = 0;
        		long firstIndex = 0;
        		if(new_type == XBSpline){
        			firstIndex = basis.mra.rangeI(new_j).firstIndex();
        			lastIndex = basis.mra.rangeI(new_j).lastIndex();
        		}
        		else{
        			firstIndex = basis.rangeJ(new_j).firstIndex();
        			lastIndex = basis.rangeJ(new_j).lastIndex();
        		}

            	for (long new_k=firstIndex; new_k<=new_k_last; ++new_k) {
                    PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                    if(minimal_overlap(covered_supp, supp) >= supp.length()){
                        Index1D new_index1d(new_j,new_k,new_type);
                        if (v.find(new_index1d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                        foundCoveredSupp = true;
                        covering_k = std::min(covering_k, new_k);
                    }
        		}
                if(!foundPredecessor){
                	if(!foundCoveredSupp){ covering_k = lastIndex;}
            		for(long new_k = new_k_first; new_k <= lastIndex; ++new_k){
                        PeriodicSupport<typename Basis::T> covered_supp = basis.generator(new_type).support(new_j,new_k);
                        if(minimal_overlap(covered_supp, supp) >= supp.length()){
                            Index1D new_index1d(new_j,new_k,new_type);
                            if (v.find(new_index1d)!=v.end()) {
                                foundPredecessor = true;
                                break;
                            }
                            foundCoveredSupp = true;
                            covering_k = std::min(covering_k, new_k);
                        }
                    }
                }
        	}

            if (!foundPredecessor) {
            	if(foundCoveredSupp){
                    Index1D new_index1d(new_j,covering_k,new_type);
                    completeMultiTree(basis,new_index1d,v,diff_v,sparsetree);
            	}
            	else{ // we have to find more than 1 bf covering the complete support
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type == XBSpline){
						firstIndex = basis.mra.rangeI(new_j).firstIndex();
						lastIndex = basis.mra.rangeI(new_j).lastIndex();
					}
					else{
						firstIndex = basis.rangeJ(new_j).firstIndex();
						lastIndex = basis.rangeJ(new_j).lastIndex();
					}
					if(new_k_first < new_k_last){
						for (long new_k=new_k_first; new_k<=new_k_last; ++new_k) {
							PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
							// left boundary of periodic support: max(l1, li2)
							if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
								T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
								bool gap = supp.gaplength() > 0 ? true : false;
								bool wrapped = false;
								if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
									long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
									completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,diff_v,sparsetree);
									if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
										wrapped = true;
									}
								}
								completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
								if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
									wrapped = true;
								}

								while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
									new_k++;
									if(new_k > lastIndex) new_k = firstIndex;
									new_supp = basis.generator(new_type).support(new_j,new_k);
									completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
									if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
										wrapped = true;
									}
								}
								break;
							}
						}
					}
					else{
						bool is_break = false;
						for (long new_k=new_k_first; new_k <= lastIndex ; ++new_k) {
							PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
							if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
								T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
								bool gap = supp.gaplength() > 0 ? true : false;
								bool wrapped = false;
								if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
									long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
									completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,diff_v,sparsetree);
									if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
										wrapped = true;
									}
								}
								completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
								if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
									wrapped = true;
								}

								while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
									new_k++;
									if(new_k > lastIndex){
										is_break = false;
										break;
									}
									new_supp = basis.generator(new_type).support(new_j,new_k);
									completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
									if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
										wrapped = true;
									}
									is_break = true;
								}
								break;
							}
						}
						if(is_break == false){
							for (long new_k=firstIndex; new_k <=new_k_last ; ++new_k) {
								PeriodicSupport<typename Basis::T> new_supp = basis.generator(new_type).support(new_j,new_k);
								if (std::max(new_supp.l1, new_supp.li2) >= std::max(supp.l1, supp.li2)) {
									T rightbd = supp.li1 > 0 ? supp.li1 : supp.l2;
									bool gap = supp.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp.l1, new_supp.li2) > std::max(supp.l1, supp.li2)) {
										long left_k = new_k-1 >= firstIndex ? new_k - 1 : lastIndex;
										completeMultiTree(basis, Index1D(new_j,left_k,new_type),v,diff_v,sparsetree);
										if(basis.generator(new_type).support(new_j,left_k).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
									if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp.li1 > 0 ? new_supp.li1 : new_supp.l2) < rightbd ||( gap && !wrapped)){
										new_k++;
										if(new_k > lastIndex) new_k = firstIndex;
										new_supp = basis.generator(new_type).support(new_j,new_k);
										completeMultiTree(basis, Index1D(new_j,new_k,new_type),v,diff_v,sparsetree);
										if(basis.generator(new_type).support(new_j,new_k).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
					}
            	}
            }
        }
    }
    return;
}



// Non-Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v, int coordDirec, bool sparsetree,
                  bool isAlreadyMultiTree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (isAlreadyMultiTree) {
        if (v.find(index2d)!=v.end())  return;
        else                           v[index2d] = 0.;
    }
    else {
        //std::cerr << "     Completion to multitree for index = " << index2d << std::endl;
        if (v.find(index2d)==v.end())  v[index2d] = 0.;
        //std::cerr << "     Node inserted for index = " << index2d << std::endl;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
        Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
        //check x-direction
        int new_j_x = 0;
        long new_k_x_first = 0, new_k_x_last = 0;
        bool checkPredecessors=true;
        XType new_type_x = XWavelet;
        if (j_x==j0_x && index_x.xtype==XWavelet) {
            basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XBSpline;
            assert(new_j_x==j_x);
        }
        else if (j_x>j0_x && index_x.xtype==XWavelet) {
            basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XWavelet;
            assert(new_j_x==j_x-1);
        }
        else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (overlap(supp_x,new_supp_x)>0) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                bool foundCoveredSupp = false;
                long covering_k = new_k_x_last;
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> covered_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (covered_supp.l1<=supp_x.l1 && covered_supp.l2>=supp_x.l2) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    	foundCoveredSupp = true;
                    	covering_k = std::min(covering_k, new_k_x);
                    }
                }
                if (!foundPredecessor) {
                	if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,coordDirec,sparsetree);
                	}
                	else{
                        for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                        //    std::cerr << "¤¤¤ ====> X-Direction: Support to be covered: " << supp_x << std::endl;
                            Support<typename Basis::T> new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                            if (new_supp.l1 > supp_x.l1) {
                       //         std::cerr << "¤¤¤       X-Covering with: " << basis.first.generator(new_type_x).support(new_j_x,new_k_x-1) << std::endl;
                       //         std::cerr << "¤¤¤       X-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y),v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_x.l2){
                                	new_k_x++;
                                	new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                       //             std::cerr << "¤¤¤       X-Covering with: " << new_supp << std::endl;
                      //              std::cout << "add " << Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y) << std::endl;
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
                                }
                       //         std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
                        }
                	}
                }
            }
        }
    }
    if (coordDirec==0 || coordDirec==2) {
        Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
        //check y-direction
        bool checkPredecessors=true;
        int new_j_y = 0;
        long new_k_y_first = 0, new_k_y_last = 0;
        XType new_type_y = XWavelet;
        if (j_y==j0_y && index_y.xtype==XWavelet) {
            basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XBSpline;
            assert(new_j_y==j_y);
        }
        else if (j_y>j0_y && index_y.xtype==XWavelet) {
            basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XWavelet;
            assert(new_j_y==j_y-1);
        }
        else {
            checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (overlap(supp_y,new_supp_y)>0) {
                        Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                bool foundCoveredSupp = false;
                long covering_k = new_k_y_last;
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> covered_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (covered_supp.l1<=supp_y.l1 && covered_supp.l2>=supp_y.l2) {
                        Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    	foundCoveredSupp = true;
                    	covering_k = std::min(covering_k, new_k_y);
                    }
                }
                if (!foundPredecessor) {
                	if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,coordDirec,sparsetree);
                	}
                	else{
                        for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                       //     std::cerr << "¤¤¤ ==== > Y-Direction: Support to be covered: " << supp_y << std::endl;
                            Support<typename Basis::T> new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                            if (new_supp.l1 > supp_y.l1) {
                      //          std::cerr << "¤¤¤       Y-Covering with: " << basis.second.generator(new_type_y).support(new_j_y,new_k_y-1) << std::endl;
                      //          std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y-1,new_type_y)),v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_y.l2){
                                	new_k_y++;
                                	new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                       //             std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                    completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,coordDirec,sparsetree);
                                }
                      //          std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
                        }
                	}
                }
            }
        }
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v, int coordDirec, bool sparsetree,
                    bool isAlreadyMultiTree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (isAlreadyMultiTree) {
        if (v.find(index2d)!=v.end())  return;
        else                           v[index2d] = 0.;
    }
    else {
        //std::cerr << "     Completion to multitree for index = " << index2d << std::endl;
        if (v.find(index2d)==v.end())  v[index2d] = 0.;
        //std::cerr << "     Node inserted for index = " << index2d << std::endl;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

	if (coordDirec==0 || coordDirec==1) {
		PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
					for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
	                        foundCoveredSupp = true;
	                        covering_k = std::min(covering_k, new_k_x);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
	                        foundCoveredSupp = true;
	                        covering_k = std::min(covering_k, new_k_x);
						}
					}
					if(!foundPredecessor){
	                	if(!foundCoveredSupp){ covering_k = lastIndex;}
						for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
							PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
								Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
		                        foundCoveredSupp = true;
		                        covering_k = std::min(covering_k, new_k_x);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_x == XBSpline){
							firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
							lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
						}
						else{
							firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
							lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
						}

						if(new_k_x_first < new_k_x_last){
							for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
										if(new_k_x > lastIndex) new_k_x = firstIndex;
										new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}

						else{
							bool is_break = false;
							for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
										completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
			                          	if(new_k_x > lastIndex){
											is_break = false;
											break;
										}
			                          	new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
		                                is_break = true;
									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
									PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
										T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
										bool gap = supp_x.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
											long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
											completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
											new_k_x++;
											if(new_k_x > lastIndex) new_k_x = firstIndex;
											new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
											completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
    }
    if (coordDirec==0 || coordDirec==2) {
		PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
					}
				}

				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
					for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
					if(!foundPredecessor){
	                	if(!foundCoveredSupp){ covering_k = lastIndex;}
						for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
							PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
							if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
								Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_y);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_y == XBSpline){
							firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
							lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
						}
						else{
							firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
							lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
						}
						if(new_k_y_first < new_k_y_last){
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
			                            completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex) new_k_y = firstIndex;
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
			                            completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex){
											is_break = false;
											break;
										}
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
										is_break=true;

									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
									PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
										T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
										bool gap = supp_y.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
											long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
				                            completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec, sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
											new_k_y++;
											if(new_k_y > lastIndex) new_k_y = firstIndex;
											new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// Periodic-NonPeriodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v, int coordDirec, bool sparsetree,
                    bool isAlreadyMultiTree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (isAlreadyMultiTree) {
        if (v.find(index2d)!=v.end())  return;
        else                           v[index2d] = 0.;
    }
    else {
        //std::cerr << "     Completion to multitree for index = " << index2d << std::endl;
        if (v.find(index2d)==v.end())  v[index2d] = 0.;
        //std::cerr << "     Node inserted for index = " << index2d << std::endl;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
		PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
					for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_x);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_x);
						}
					}
					if(!foundPredecessor){
						if(!foundCoveredSupp){covering_k = lastIndex;}
						for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
							PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
								Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_x);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_x == XBSpline){
							firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
							lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
						}
						else{
							firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
							lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
						}
						if(new_k_x_first < new_k_x_last){
							for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
										if(new_k_x > lastIndex) new_k_x = firstIndex;
										new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
										completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
			                          	if(new_k_x > lastIndex){
											is_break = false;
											break;
										}
			                          	new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
		                                is_break = true;
									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
									PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
										T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
										bool gap = supp_x.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
											long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
											completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
											new_k_x++;
											if(new_k_x > lastIndex) new_k_x = firstIndex;
											new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
											completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec==0 || coordDirec==2) {
		Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					Support<typename Basis::T> covered_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (covered_supp.l1<=supp_y.l1 && covered_supp.l2>=supp_y.l2) {
						Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)!=v.end()) {
							foundPredecessor = true;
							break;
						}
						foundCoveredSupp = true;
						covering_k = std::min(covering_k, new_k_y);
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,coordDirec,sparsetree);
					}
					else{
						for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                          //  std::cerr << "¤¤¤ ==== > Y-Direction: Support to be covered: " << supp_y << std::endl;
							Support<typename Basis::T> new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                            if (new_supp.l1 > supp_y.l1) {
                         //       std::cerr << "¤¤¤       Y-Covering with: " << basis.second.generator(new_type_y).support(new_j_y,new_k_y-1) << std::endl;
                         //       std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y-1,new_type_y)),v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_y.l2){
                                	new_k_y++;
                                	new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                         //           std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                    completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,coordDirec,sparsetree);
                                }
                        //        std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
						}
					}
				}
			}
		}
	}
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// NonPeriodic-Periodic Version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v, int coordDirec, bool sparsetree,
                    bool isAlreadyMultiTree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (isAlreadyMultiTree) {
        if (v.find(index2d)!=v.end())  return;
        else                           v[index2d] = 0.;
    }
    else {
        //std::cerr << "     Completion to multitree for index = " << index2d << std::endl;
        if (v.find(index2d)==v.end())  v[index2d] = 0.;
        //std::cerr << "     Node inserted for index = " << index2d << std::endl;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
		Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
					Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
					if (overlap(supp_x,new_supp_x)>0) {
						Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
					Support<typename Basis::T> covered_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
					if (covered_supp.l1<=supp_x.l1 && covered_supp.l2>=supp_x.l2) {
						Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
						if (v.find(new_index2d)!=v.end()) {
							foundPredecessor = true;
							break;
						}
						foundCoveredSupp = true;
						covering_k = std::min(covering_k, new_k_x);
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,coordDirec,sparsetree);
					}
					else{
						for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
							Support<typename Basis::T> new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if (new_supp.l1 > supp_x.l1) {
								completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y),v,coordDirec,sparsetree);
								completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
								while(new_supp.l2 < supp_x.l2){
									new_k_x++;
									new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									std::cout << "add " << Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y) << std::endl;
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,coordDirec,sparsetree);
								}
								break;
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec==0 || coordDirec==2) {
		PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
					}
				}

				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
					for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
					if(!foundPredecessor){
						if(!foundCoveredSupp){covering_k = lastIndex;}
						for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
							PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
							if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
								Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_y);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_y == XBSpline){
							firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
							lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
						}
						else{
							firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
							lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
						}
						if(new_k_y_first < new_k_y_last){
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex) new_k_y = firstIndex;
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec,sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex){
											is_break = false;
											break;
										}
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
										is_break=true;

									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
									PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
										T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
										bool gap = supp_y.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
											long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,coordDirec,sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
											new_k_y++;
											if(new_k_y > lastIndex) new_k_y = firstIndex;
											new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, coordDirec, sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// Non-Periodic Version + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v, int coordDirec, bool sparsetree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (v.find(index2d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index2d);
    	v[index2d] = 0.;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
        Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
        //check x-direction
        int new_j_x = 0;
        long new_k_x_first = 0, new_k_x_last = 0;
        bool checkPredecessors=true;
        XType new_type_x = XWavelet;
        if (j_x==j0_x && index_x.xtype==XWavelet) {
            basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XBSpline;
            assert(new_j_x==j_x);
        }
        else if (j_x>j0_x && index_x.xtype==XWavelet) {
            basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XWavelet;
            assert(new_j_x==j_x-1);
        }
        else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (overlap(supp_x,new_supp_x)>0) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                bool foundCoveredSupp = false;
                long covering_k = new_k_x_last;
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> covered_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (covered_supp.l1<=supp_x.l1 && covered_supp.l2>=supp_x.l2) {
                        Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    	foundCoveredSupp = true;
                    	covering_k = std::min(covering_k, new_k_x);
                    }
                }
                if (!foundPredecessor) {
                	if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
                	}
                	else{
                        //std::cerr << "¤¤¤ ====> X-Direction: Support to be covered: " << supp_x << std::endl;
                        for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                            Support<typename Basis::T> new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                            if (new_supp.l1 > supp_x.l1) {
                                //std::cerr << "¤¤¤       X-Covering with: " << basis.first.generator(new_type_x).support(new_j_x,new_k_x-1) << std::endl;
                                //std::cerr << "¤¤¤       X-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_x.l2){
                                	new_k_x++;
                                	new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                                   // std::cerr << "¤¤¤       X-Covering with: " << new_supp << std::endl;
                                   // std::cout << "add " << Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y) << std::endl;
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
                                }
                               // std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
                        }
                	}
                }
            }
        }
    }
    if (coordDirec==0 || coordDirec==2) {
        Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
        //check y-direction
        bool checkPredecessors=true;
        int new_j_y = 0;
        long new_k_y_first = 0, new_k_y_last = 0;
        XType new_type_y = XWavelet;
        if (j_y==j0_y && index_y.xtype==XWavelet) {
            basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XBSpline;
            assert(new_j_y==j_y);
        }
        else if (j_y>j0_y && index_y.xtype==XWavelet) {
            basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XWavelet;
            assert(new_j_y==j_y-1);
        }
        else {
            checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (overlap(supp_y,new_supp_y)>0) {
                        Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                bool foundCoveredSupp = false;
                long covering_k = new_k_y_last;
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> covered_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (covered_supp.l1<=supp_y.l1 && covered_supp.l2>=supp_y.l2) {
                        Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
                        if (v.find(new_index2d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    	foundCoveredSupp = true;
                    	covering_k = std::min(covering_k, new_k_y);
                    }
                }
                if (!foundPredecessor) {
                	if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,diff_v,coordDirec,sparsetree);
                	}
                	else{
                        for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                           // std::cerr << "¤¤¤ ==== > Y-Direction: Support to be covered: " << supp_y << std::endl;
                            Support<typename Basis::T> new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                            if (new_supp.l1 > supp_y.l1) {
                            //    std::cerr << "¤¤¤       Y-Covering with: " << basis.second.generator(new_type_y).support(new_j_y,new_k_y-1) << std::endl;
                            //    std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y-1,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_y.l2){
                                	new_k_y++;
                                	new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                             //       std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                    completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                }
                             //   std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
                        }
                	}
                }
            }
        }
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// Periodic Version + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v, int coordDirec, bool sparsetree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (v.find(index2d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index2d);
    	v[index2d] = 0.;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

	if (coordDirec==0 || coordDirec==1) {
		PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
					for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
	                        foundCoveredSupp = true;
	                        covering_k = std::min(covering_k, new_k_x);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
	                        foundCoveredSupp = true;
	                        covering_k = std::min(covering_k, new_k_x);
						}
					}
					if(!foundPredecessor){
	                	if(!foundCoveredSupp){ covering_k = lastIndex;}
						for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
							PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
								Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
		                        foundCoveredSupp = true;
		                        covering_k = std::min(covering_k, new_k_x);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
                        completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_x == XBSpline){
							firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
							lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
						}
						else{
							firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
							lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
						}
						if(new_k_x_first < new_k_x_last){
							for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
										if(new_k_x > lastIndex) new_k_x = firstIndex;
										new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
										completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
			                          	if(new_k_x > lastIndex){
											is_break = false;
											break;
										}
			                          	new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
		                                is_break = true;
									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
									PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
										T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
										bool gap = supp_x.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
											long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
											completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
											new_k_x++;
											if(new_k_x > lastIndex) new_k_x = firstIndex;
											new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
											completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
    }
    if (coordDirec==0 || coordDirec==2) {
		PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
					}
				}

				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
					for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
					if(!foundPredecessor){
	                	if(!foundCoveredSupp){ covering_k = lastIndex;}
						for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
							PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
							if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
								Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_y);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,diff_v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_y == XBSpline){
							firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
							lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
						}
						else{
							firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
							lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
						}
						if(new_k_y_first < new_k_y_last){
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex) new_k_y = firstIndex;
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v, coordDirec,sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex){
											is_break = false;
											break;
										}
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
										is_break=true;

									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
									PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
										T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
										bool gap = supp_y.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
											long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v, coordDirec,sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
											new_k_y++;
											if(new_k_y > lastIndex) new_k_y = firstIndex;
											new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// Periodic-NonPeriodic Version + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v, int coordDirec, bool sparsetree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (v.find(index2d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index2d);
    	v[index2d] = 0.;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
		PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
					for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
						PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if (overlap(supp_x,new_supp_x)>0) {
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				if(new_k_x_first < new_k_x_last){
					for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_x);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_x == XBSpline){
						firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
						lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
					}
					else{
						firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
						lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
					}
					for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
						PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
						if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
							Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_x);
						}
					}
					if(!foundPredecessor){
						if(!foundCoveredSupp){covering_k = lastIndex;}
						for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
							PeriodicSupport<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if(minimal_overlap(covered_supp_x, supp_x) >= supp_x.length()){
								Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_x);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_x == XBSpline){
							firstIndex = basis.first.mra.rangeI(new_j_x).firstIndex();
							lastIndex = basis.first.mra.rangeI(new_j_x).lastIndex();
						}
						else{
							firstIndex = basis.first.rangeJ(new_j_x).firstIndex();
							lastIndex = basis.first.rangeJ(new_j_x).lastIndex();
						}
						if(new_k_x_first < new_k_x_last){
							for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
										if(new_k_x > lastIndex) new_k_x = firstIndex;
										new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
	                                    completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_x=new_k_x_first; new_k_x<=lastIndex; ++new_k_x) {
								PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
									T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
									bool gap = supp_x.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
										long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
										completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
									if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
										new_k_x++;
			                          	if(new_k_x > lastIndex){
											is_break = false;
											break;
										}
			                          	new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}
		                                is_break = true;
									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_x=firstIndex; new_k_x<=new_k_x_last; ++new_k_x) {
									PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_x.l1, new_supp_x.li2) >= std::max(supp_x.l1, supp_x.li2)) {
										T rightbd = supp_x.li1 > 0 ? supp_x.li1 : supp_x.l2;
										bool gap = supp_x.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_x.l1, new_supp_x.li2) > std::max(supp_x.l1, supp_x.li2)) {
											long left_k_x = new_k_x-1 >= firstIndex ? new_k_x - 1 : lastIndex;
											completeMultiTree(basis, Index2D(Index1D(new_j_x,left_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,left_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
										if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_x.li1 > 0 ? new_supp_x.li1 : new_supp_x.l2) < rightbd ||( gap && !wrapped)){
											new_k_x++;
											if(new_k_x > lastIndex) new_k_x = firstIndex;
											new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
											completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
											if(basis.first.generator(new_type_x).support(new_j_x,new_k_x).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec==0 || coordDirec==2) {
		Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					Support<typename Basis::T> covered_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (covered_supp.l1<=supp_y.l1 && covered_supp.l2>=supp_y.l2) {
						Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)!=v.end()) {
							foundPredecessor = true;
							break;
						}
						foundCoveredSupp = true;
						covering_k = std::min(covering_k, new_k_y);
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,diff_v,coordDirec,sparsetree);
					}
					else{
						for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                          //  std::cerr << "¤¤¤ ==== > Y-Direction: Support to be covered: " << supp_y << std::endl;
							Support<typename Basis::T> new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                            if (new_supp.l1 > supp_y.l1) {
                           //     std::cerr << "¤¤¤       Y-Covering with: " << basis.second.generator(new_type_y).support(new_j_y,new_k_y-1) << std::endl;
                          //      std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y-1,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                while(new_supp.l2 < supp_y.l2){
                                	new_k_y++;
                                	new_supp = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                          //          std::cerr << "¤¤¤       Y-Covering with: " << new_supp << std::endl;
                                    completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
                                }
                         //       std::cerr << "¤¤¤ <==============   " << std::endl;
                                break;
                            }
						}
					}
				}
			}
		}
	}
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

// NonPeriodic-Periodic Version + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v, int coordDirec, bool sparsetree)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (v.find(index2d)!=v.end()){
    	return;
    }
    else{
    	diff_v.insert(index2d);
    	v[index2d] = 0.;
    }


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  j_x = index_x.j, j_y = index_y.j;
    long k_x = index_x.k, k_y = index_y.k;

    if (coordDirec==0 || coordDirec==1) {
		Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
		//check x-direction
		int new_j_x = 0;
		long new_k_x_first = 0, new_k_x_last = 0;
		bool checkPredecessors=true;
		XType new_type_x = XWavelet;
		if (j_x==j0_x && index_x.xtype==XWavelet) {
			basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XBSpline;
			assert(new_j_x==j_x);
		}
		else if (j_x>j0_x && index_x.xtype==XWavelet) {
			basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
			new_type_x = XWavelet;
			assert(new_j_x==j_x-1);
		}
		else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
					Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
					if (overlap(supp_x,new_supp_x)>0) {
						Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_x_last;
				for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
					Support<typename Basis::T> covered_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
					if (covered_supp.l1<=supp_x.l1 && covered_supp.l2>=supp_x.l2) {
						Index2D new_index2d(Index1D(new_j_x,new_k_x,new_type_x),index_y);
						if (v.find(new_index2d)!=v.end()) {
							foundPredecessor = true;
							break;
						}
						foundCoveredSupp = true;
						covering_k = std::min(covering_k, new_k_x);
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(Index1D(new_j_x,covering_k,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
					}
					else{
						for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
							Support<typename Basis::T> new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
							if (new_supp.l1 > supp_x.l1) {
								completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
								completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
								while(new_supp.l2 < supp_x.l2){
									new_k_x++;
									new_supp = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
									std::cout << "add " << Index2D(Index1D(new_j_x,new_k_x-1,new_type_x),index_y) << std::endl;
									completeMultiTree(basis, Index2D(Index1D(new_j_x,new_k_x,new_type_x),index_y),v,diff_v,coordDirec,sparsetree);
								}
								break;
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec==0 || coordDirec==2) {
		PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
		//check y-direction
		bool checkPredecessors=true;
		int new_j_y = 0;
		long new_k_y_first = 0, new_k_y_last = 0;
		XType new_type_y = XWavelet;
		if (j_y==j0_y && index_y.xtype==XWavelet) {
			basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XBSpline;
			assert(new_j_y==j_y);
		}
		else if (j_y>j0_y && index_y.xtype==XWavelet) {
			basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
			new_type_y = XWavelet;
			assert(new_j_y==j_y-1);
		}
		else {
			checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
		}

		if (checkPredecessors) {
			if (!sparsetree) {
				for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
					PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
					if (overlap(supp_y,new_supp_y)>0) {
						Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
						if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
					}
				}

				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
					for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
						PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if (overlap(supp_y,new_supp_y)>0) {
							Index2D new_index2d(index_x,Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v,diff_v,coordDirec);
						}
					}
				}
			}
			else {
				bool foundPredecessor = false;
				bool foundCoveredSupp = false;
				long covering_k = new_k_y_last;
				if(new_k_y_first < new_k_y_last){
					for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
				}
				else{
					long lastIndex = 0;
					long firstIndex = 0;
					if(new_type_y == XBSpline){
						firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
						lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
					}
					else{
						firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
						lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
					}
					for (long new_k_y=firstIndex; new_k_y<=new_k_y_last; ++new_k_y) {
						PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
						if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
							Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
							if (v.find(new_index2d)!=v.end()) {
								foundPredecessor = true;
								break;
							}
							foundCoveredSupp = true;
							covering_k = std::min(covering_k, new_k_y);
						}
					}
					if(!foundPredecessor){
						if(!foundCoveredSupp){covering_k = lastIndex;}
						for (long new_k_y=new_k_y_first; new_k_y<=lastIndex; ++new_k_y) {
							PeriodicSupport<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
							if(minimal_overlap(covered_supp_y, supp_y) >= supp_y.length()){
								Index2D new_index2d(index_x, Index1D(new_j_y,new_k_y,new_type_y));
								if (v.find(new_index2d)!=v.end()) {
									foundPredecessor = true;
									break;
								}
								foundCoveredSupp = true;
								covering_k = std::min(covering_k, new_k_y);
							}
						}
					}
				}
				if (!foundPredecessor) {
					if(foundCoveredSupp){
						completeMultiTree(basis,Index2D(index_x, Index1D(new_j_y,covering_k,new_type_y)),v,diff_v,coordDirec,sparsetree);
					}
					else{ // we have to find more than 1 bf covering the complete support
						long lastIndex = 0;
						long firstIndex = 0;
						if(new_type_y == XBSpline){
							firstIndex = basis.second.mra.rangeI(new_j_y).firstIndex();
							lastIndex = basis.second.mra.rangeI(new_j_y).lastIndex();
						}
						else{
							firstIndex = basis.second.rangeJ(new_j_y).firstIndex();
							lastIndex = basis.second.rangeJ(new_j_y).lastIndex();
						}
						if(new_k_y_first < new_k_y_last){
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v,coordDirec,sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex) new_k_y = firstIndex;
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									break;
								}
							}
						}
						else{
							bool is_break = false;
							for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
								PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
								// left boundary of periodic support: max(l1, li2)
								if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
									T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
									bool gap = supp_y.gaplength() > 0 ? true : false;
									bool wrapped = false;
									if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
										long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v, coordDirec,sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
											wrapped = true;
										}
									}
									completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
									if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
										wrapped = true;
									}

									while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
										new_k_y++;
										if(new_k_y > lastIndex){
											is_break = false;
											break;
										}
										new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}
										is_break=true;

									}
									break;
								}
							}
							if(is_break==false){
								for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
									PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
									// left boundary of periodic support: max(l1, li2)
									if (std::max(new_supp_y.l1, new_supp_y.li2) >= std::max(supp_y.l1, supp_y.li2)) {
										T rightbd = supp_y.li1 > 0 ? supp_y.li1 : supp_y.l2;
										bool gap = supp_y.gaplength() > 0 ? true : false;
										bool wrapped = false;
										if (std::max(new_supp_y.l1, new_supp_y.li2) > std::max(supp_y.l1, supp_y.li2)) {
											long left_k_y = new_k_y-1 >= firstIndex ? new_k_y - 1 : lastIndex;
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,left_k_y,new_type_y)),v,diff_v, coordDirec,sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,left_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
										if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
											wrapped = true;
										}

										while((new_supp_y.li1 > 0 ? new_supp_y.li1 : new_supp_y.l2) < rightbd ||( gap && !wrapped)){
											new_k_y++;
											if(new_k_y > lastIndex) new_k_y = firstIndex;
											new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
											completeMultiTree(basis, Index2D(index_x, Index1D(new_j_y,new_k_y,new_type_y)), v, diff_v,coordDirec, sparsetree);
											if(basis.second.generator(new_type_y).support(new_j_y,new_k_y).gaplength() > 0){
												wrapped = true;
											}
										}
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}



/*
 * TODO periodic
 */
template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v, int coordDirec, bool sparsetree)
{
	assert(Basis::Domain != Periodic);

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    if (v.find(index3d)!=v.end())  return;
    else                           v[index3d] = 0.;

    Index1D index_x = index3d.index1;
    Index1D index_y = index3d.index2;
    Index1D index_z = index3d.index3;

    int  j_x = index_x.j, j_y = index_y.j, j_z = index_z.j;
    long k_x = index_x.k, k_y = index_y.k, k_z = index_z.k;

    if (coordDirec==0 || coordDirec==1) {
        Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
        //check x-direction
        int new_j_x = 0;
        long new_k_x_first = 0, new_k_x_last = 0;
        bool checkPredecessors=true;
        XType new_type_x = XWavelet;
        if (j_x==j0_x && index_x.xtype==XWavelet) {
            basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XBSpline;
            assert(new_j_x==j_x);
        }
        else if (j_x>j0_x && index_x.xtype==XWavelet) {
            basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,new_j_x,new_k_x_first,new_k_x_last);
            new_type_x = XWavelet;
            assert(new_j_x==j_x-1);
        }
        else checkPredecessors = false;    // Index corresponds to a scaling function -> no predecessor

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> new_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (overlap(supp_x,new_supp_x)>0) {
                        Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                    Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                    if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                        Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }

                if (!foundPredecessor) {
                    for (long new_k_x=new_k_x_first; new_k_x<=new_k_x_last; ++new_k_x) {
                        Support<typename Basis::T> covered_supp_x = basis.first.generator(new_type_x).support(new_j_x,new_k_x);
                        if (covered_supp_x.l1<=supp_x.l1 && covered_supp_x.l2>=supp_x.l2) {
                            Index3D new_index3d(Index1D(new_j_x,new_k_x,new_type_x),index_y,index_z);
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec==0 || coordDirec==2) {
        Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
        //check y-direction
        bool checkPredecessors=true;
        int new_j_y = 0;
        long new_k_y_first = 0, new_k_y_last = 0;
        XType new_type_y = XWavelet;
        if (j_y==j0_y && index_y.xtype==XWavelet) {
            basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XBSpline;
            assert(new_j_y==j_y);
        }
        else if (j_y>j0_y && index_y.xtype==XWavelet) {
            basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,new_j_y,new_k_y_first,new_k_y_last);
            new_type_y = XWavelet;
            assert(new_j_y==j_y-1);
        }
        else {
            checkPredecessors=false; // no "return" here!!! we also need to check the next "if"-clause;
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> new_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (overlap(supp_y,new_supp_y)>0) {
                        Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                    Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                    if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                        Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }

                if (!foundPredecessor) {
                    for (long new_k_y=new_k_y_first; new_k_y<=new_k_y_last; ++new_k_y) {
                        Support<typename Basis::T> covered_supp_y = basis.second.generator(new_type_y).support(new_j_y,new_k_y);
                        if (covered_supp_y.l1<=supp_y.l1 && covered_supp_y.l2>=supp_y.l2) {
                            Index3D new_index3d(index_x,Index1D(new_j_y,new_k_y,new_type_y),index_z);
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec== 0 || coordDirec==3) {
        Support<typename Basis::T> supp_z = basis.third.generator(index_z.xtype).support(j_z,k_z);
        //check z-direction
        bool checkPredecessors=true;
        int new_j_z = 0;
        long new_k_z_first = 0, new_k_z_last = 0;
        XType new_type_z = XWavelet;
        if (j_z==j0_z && index_z.xtype==XWavelet) {
            basis.third.getScalingNeighborsForWavelet(j_z,k_z,basis.third,new_j_z,new_k_z_first,new_k_z_last);
            new_type_z = XBSpline;
            assert(new_j_z==j_z);
        }
        else if (j_z>j0_z && index_z.xtype==XWavelet) {
            basis.third.getLowerWaveletNeighborsForWavelet(j_z,k_z,basis.third,new_j_z,new_k_z_first,new_k_z_last);
            new_type_z = XWavelet;
            assert(new_j_z==j_z-1);
        }
        else {
            checkPredecessors=false; // no "return" here, see above
        }

        if (checkPredecessors) {
            if (!sparsetree) {
                for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                    Support<typename Basis::T> new_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                    if (overlap(supp_z,new_supp_z)>0) {
                        Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                        if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v,coordDirec);
                    }
                }
            }
            else {
                bool foundPredecessor = false;
                for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                    Support<typename Basis::T> covered_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                    if (covered_supp_z.l1<=supp_z.l1 && covered_supp_z.l2>=supp_z.l2) {
                        Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                        if (v.find(new_index3d)!=v.end()) {
                            foundPredecessor = true;
                            break;
                        }
                    }
                }
                if (!foundPredecessor) {
                    for (long new_k_z=new_k_z_first; new_k_z<=new_k_z_last; ++new_k_z) {
                        Support<typename Basis::T> covered_supp_z = basis.third.generator(new_type_z).support(new_j_z,new_k_z);
                        if (covered_supp_z.l1<=supp_z.l1 && covered_supp_z.l2>=supp_z.l2) {
                            Index3D new_index3d(index_x,index_y,Index1D(new_j_z,new_k_z,new_type_z));
                            completeMultiTree(basis,new_index3d,v,coordDirec,sparsetree);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (coordDirec != 0 && coordDirec != 1 && coordDirec != 2 && coordDirec !=3 ) {
        std::cerr << "completeMultiTree: non-admissible coordinate direction " << coordDirec
                  << std::endl;
    }
    return;
}

/*
template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    if (v.find(index3d)!=v.end())  return;
    else                           v[index3d] = 0.;


    Index1D index_x = index3d.index1;
    Index1D index_y = index3d.index2;
    Index1D index_z = index3d.index3;

    int  j_x = index_x.j, j_y = index_y.j, j_z = index_z.j;
    long k_x = index_x.k, k_y = index_y.k, k_z = index_z.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(j_x,k_x);
    //check x-direction
    if (j_x==j0_x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(j_x,k_x,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(j_x==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index3D new_index3d(Index1D(j_scaling,k_scaling,XBSpline),index_y,index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_x>j0_x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(j_x,k_x,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_x-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index3D new_index3d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y,index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(j_y,k_y);
    //check y-direction
    if (j_y==j0_y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(j_y,k_y,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(j_y==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index3D new_index3d(index_x,Index1D(j_scaling,k_scaling,XBSpline),index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_y>j0_y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(j_y,k_y,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_y-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index3D new_index3d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet),index_z);
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }

    Support<typename Basis::T> supp_z = basis.third.generator(index_z.xtype).support(j_z,k_z);
    //check z-direction
    if (j_z==j0_z  && index_z.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.third.getScalingNeighborsForWavelet(j_z,k_z,basis.third,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(j_z==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_z = basis.third.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_z,new_supp_z)>0) {
                Index3D new_index3d(index_x,index_y,Index1D(j_scaling,k_scaling,XBSpline));
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    else if (j_z>j0_z && index_z.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.third.getLowerWaveletNeighborsForWavelet(j_z,k_z,basis.third,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(j_z-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_z = basis.third.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_z,new_supp_z)>0) {
                Index3D new_index3d(index_x,index_y,Index1D(j_wavelet,k_wavelet,XWavelet));
                if (v.find(new_index3d)==v.end()) completeMultiTree(basis,new_index3d,v);
            }
        }
    }
    return;
}
*/

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index2D> &v, int j, T gamma)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    for (long k1=basis.first.mra.rangeI(j0_x).firstIndex(); k1<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
        for (long k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
            Index1D row(j0_x,k1,XBSpline);
            Index1D col(j0_y,k2,XBSpline);
            v[Index2D(row,col)] = 0.;
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0_y+i2-1;
            for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_x,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                v[Index2D(row,col)] = 0.;
                v[Index2D(col,row)] = 0.;
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_x+i1-1;
        for (int i2=1; i2<=j; ++i2) {
            if (T(i1+i2)-gamma*std::max(i1,i2)>(1-gamma)*j) {
                continue;
            }
            int j2=j0_y+i2-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    v[Index2D(row,col)] = 0.;
                }
            }
        }
    }
    return;
}

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index3D> &v, int J, T gamma)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    int j0_z = basis.third.j0;

    for (int i_x=0; i_x<=J; ++i_x) {
        flens::Range<int> range_x(0,0);
        XType xtype_x;
        int j_x;
        if (i_x==0) {   range_x = basis.first.mra.rangeI(j0_x);   j_x = j0_x;       xtype_x = XBSpline; }
        else        {   range_x = basis.first.rangeJ(j0_x+i_x-1); j_x = j0_x+i_x-1; xtype_x = XWavelet; }
        for (int k_x=range_x.firstIndex(); k_x<=range_x.lastIndex(); ++k_x) {
            Index1D index_x(j_x, k_x, xtype_x);

            for (int i_y=0; i_y<=J-i_x; ++i_y) {
                flens::Range<int> range_y(0,0);
                XType xtype_y;
                int j_y;
                if (i_y==0) {   range_y = basis.second.mra.rangeI(j0_y);   j_y = j0_y;       xtype_y = XBSpline; }
                else        {   range_y = basis.second.rangeJ(j0_y+i_y-1); j_y = j0_y+i_y-1; xtype_y = XWavelet; }
                for (int k_y=range_y.firstIndex(); k_y<=range_y.lastIndex(); ++k_y) {
                    Index1D index_y(j_y, k_y, xtype_y);

                    for (int i_z=0; i_z<=J-i_x-i_y; ++i_z) {
                        int sum = i_x+i_y+i_z;
                        int max = std::max(i_x,std::max(i_y,i_z));
                        if (sum - gamma*max > (1-gamma)*J) continue;
                        flens::Range<int> range_z(0,0);
                        XType xtype_z;
                        int j_z;
                        if (i_z==0) {   range_z = basis.second.mra.rangeI(j0_z);   j_z = j0_z;       xtype_z = XBSpline; }
                        else        {   range_z = basis.second.rangeJ(j0_z+i_z-1); j_z = j0_z+i_z-1; xtype_z = XWavelet; }
                        for (int k_z=range_z.firstIndex(); k_z<=range_z.lastIndex(); ++k_z) {
                            Index1D index_z(j_z, k_z, xtype_z);
                            v[Index3D(index_x,index_y,index_z)] = 0.;
                        }
                    }
                }
            }
        }
    }
    return;
}

// Non-Periodic Version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (jx==j0x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(jx==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jx>j0x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jx-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (jy==j0y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(jy==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jy>j0y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jy-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
}

// Periodic Version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

	PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (jx==j0x && index_x.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
												  j_scaling,k_scaling_first,k_scaling_last);
		assert(jx==j_scaling);
		if(k_scaling_first < k_scaling_last){
			for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_scaling=basis.first.mra.rangeI(j_scaling).firstIndex(); k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_scaling=k_scaling_first; k_scaling<=basis.first.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}

	}
	else if (jx>j0x && index_x.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
													   j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jx-1==j_wavelet);
		if(k_wavelet_first < k_wavelet_last){
			for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_wavelet=basis.first.rangeJ(j_wavelet).firstIndex(); k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_wavelet=k_wavelet_first; k_wavelet<=basis.first.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}

    PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (jy==j0y  && index_y.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
												   j_scaling,k_scaling_first,k_scaling_last);
		assert(jy==j_scaling);
		if(k_scaling_first < k_scaling_last){
			for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_scaling=basis.first.mra.rangeI(j_scaling).firstIndex(); k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_scaling=k_scaling_first; k_scaling<=basis.first.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}
	else if (jy>j0y && index_y.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
														j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jy-1==j_wavelet);
		if(k_wavelet_first < k_wavelet_last){
			for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_wavelet=basis.second.rangeJ(j_wavelet).firstIndex(); k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_wavelet=k_wavelet_first; k_wavelet<=basis.second.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}
}

// Periodic-NonPeriodic Version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
				 and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

	PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (jx==j0x && index_x.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
												  j_scaling,k_scaling_first,k_scaling_last);
		assert(jx==j_scaling);
		if(k_scaling_first < k_scaling_last){
			for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_scaling=basis.first.mra.rangeI(j_scaling).firstIndex(); k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_scaling=k_scaling_first; k_scaling<=basis.first.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}

	}
	else if (jx>j0x && index_x.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
													   j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jx-1==j_wavelet);
		if(k_wavelet_first < k_wavelet_last){
			for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_wavelet=basis.first.rangeJ(j_wavelet).firstIndex(); k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_wavelet=k_wavelet_first; k_wavelet<=basis.first.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}

	Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (jy==j0y  && index_y.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
												   j_scaling,k_scaling_first,k_scaling_last);
		assert(jy==j_scaling);
		for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
			Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
			}
		}
	}
	else if (jy>j0y && index_y.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
														j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jy-1==j_wavelet);
		for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
			Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
			}
		}
	}
}

// NonPeriodic-Periodic Version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;


	Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (jx==j0x && index_x.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
												  j_scaling,k_scaling_first,k_scaling_last);
		assert(jx==j_scaling);
		for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
			Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
			}
		}
	}
	else if (jx>j0x && index_x.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
													   j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jx-1==j_wavelet);
		for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
			Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
			}
		}
	}

	PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (jy==j0y  && index_y.xtype==XWavelet) {
		int  j_scaling = 0;
		long k_scaling_first = 0, k_scaling_last = 0;
		basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
												   j_scaling,k_scaling_first,k_scaling_last);
		assert(jy==j_scaling);
		if(k_scaling_first < k_scaling_last){
			for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_scaling=basis.first.mra.rangeI(j_scaling).firstIndex(); k_scaling<=k_scaling_last; ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_scaling=k_scaling_first; k_scaling<=basis.first.mra.rangeI(j_scaling).lastIndex(); ++k_scaling) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}
	else if (jy>j0y && index_y.xtype==XWavelet) {
		int  j_wavelet = 0;
		long k_wavelet_first = 0, k_wavelet_last = 0;
		basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
														j_wavelet,k_wavelet_first,k_wavelet_last);
		assert(jy-1==j_wavelet);
		if(k_wavelet_first < k_wavelet_last){
			for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
		else{
			for (long k_wavelet=basis.second.rangeJ(j_wavelet).firstIndex(); k_wavelet<=k_wavelet_last; ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
			for (long k_wavelet=k_wavelet_first; k_wavelet<=basis.second.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
				}
			}
		}
	}
}

// Non-Periodic version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;


    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (index_x.j==j0_x) {
        for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.first.rangeJ(index_x.j-1).firstIndex(),(long)index_x.k/2 - offset);
        long k_last  = std::min((long)basis.first.rangeJ(index_x.j-1).lastIndex(),(long)index_x.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (index_y.j==j0_y) {
        for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.second.rangeJ(index_y.j-1).firstIndex(),(long)index_y.k/2 - offset);
        long k_last  = std::min((long)basis.second.rangeJ(index_y.j-1).lastIndex(),(long)index_y.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
}

// Periodic version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

	PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (index_x.j==j0_x) {
		for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
			PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = (long)index_x.k/2 - offset;
		long k_last = (long)index_x.k/2 + offset;
		bool is_wrap = false;
		while(k_first < 1){
			k_first = (long)basis.first.rangeJ(index_x.j-1).lastIndex() + k_first;
			is_wrap = true;
		}
		while(k_last > (long)basis.first.rangeJ(index_x.j-1).lastIndex()){
			k_last = k_last - (long)basis.first.rangeJ(index_x.j-1).lastIndex();
		   is_wrap = true;
		}
		if(k_first == k_last || k_first == k_last+1 || (k_first < k_last && is_wrap)){
			k_first = (long)basis.first.rangeJ(index_x.j-1).firstIndex();
			k_last = (long)basis.first.rangeJ(index_x.j-1).lastIndex();
		}

		if(k_first < k_last){
			for (long k=k_first; k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}
		else{
			for (long k=(long)basis.first.rangeJ(index_x.j-1).firstIndex(); k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
			for (long k=k_first; k<=(long)basis.first.rangeJ(index_x.j-1).lastIndex(); ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}

	}

	PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (index_y.j==j0_y) {
		for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
			PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = (long)index_y.k/2 - offset;
		long k_last = (long)index_y.k/2 + offset;
		bool is_wrap = false;
		while(k_first < 1){
			k_first = (long)basis.second.rangeJ(index_y.j-1).lastIndex() + k_first;
			is_wrap = true;
		}
		while(k_last > (long)basis.second.rangeJ(index_y.j-1).lastIndex()){
			k_last = k_last - (long)basis.second.rangeJ(index_y.j-1).lastIndex();
		   is_wrap = true;
		}
		if(k_first == k_last || k_first == k_last+1 || (k_first < k_last && is_wrap)){
			k_first = (long)basis.second.rangeJ(index_y.j-1).firstIndex();
			k_last = (long)basis.second.rangeJ(index_y.j-1).lastIndex();
		}

		if(k_first < k_last){
			for (long k=k_first; k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}
		else{
			for (long k=(long)basis.second.rangeJ(index_y.j-1).firstIndex(); k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
			for (long k=k_first; k<=(long)basis.second.rangeJ(index_y.j-1).lastIndex(); ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}

	}
}

// Periodic-NonPeriodic version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

	PeriodicSupport<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (index_x.j==j0_x) {
		for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
			PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = (long)index_x.k/2 - offset;
		long k_last = (long)index_x.k/2 + offset;
		bool is_wrap = false;
		while(k_first < 1){
			k_first = (long)basis.first.rangeJ(index_x.j-1).lastIndex() + k_first;
			is_wrap = true;
		}
		while(k_last > (long)basis.first.rangeJ(index_x.j-1).lastIndex()){
			k_last = k_last - (long)basis.first.rangeJ(index_x.j-1).lastIndex();
		   is_wrap = true;
		}
		if(k_first == k_last || k_first == k_last+1 || (k_first < k_last && is_wrap)){
			k_first = (long)basis.first.rangeJ(index_x.j-1).firstIndex();
			k_last = (long)basis.first.rangeJ(index_x.j-1).lastIndex();
		}

		if(k_first < k_last){
			for (long k=k_first; k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}
		else{
			for (long k=(long)basis.first.rangeJ(index_x.j-1).firstIndex(); k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
			for (long k=k_first; k<=(long)basis.first.rangeJ(index_x.j-1).lastIndex(); ++k) {
				PeriodicSupport<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
				if (overlap(supp_x,new_supp_x)>0) {
					Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
					//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}

	}

	Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (index_y.j==j0_y) {
		for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
			Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = std::max((long)basis.second.rangeJ(index_y.j-1).firstIndex(),(long)index_y.k/2 - offset);
		long k_last  = std::min((long)basis.second.rangeJ(index_y.j-1).lastIndex(),(long)index_y.k/2 + offset);
		for (long k=k_first; k<=k_last; ++k) {
			Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
}

// NonPeriodic-Periodic version
template <typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;


	Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
	//check x-direction
	if (index_x.j==j0_x) {
		for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
			Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = std::max((long)basis.first.rangeJ(index_x.j-1).firstIndex(),(long)index_x.k/2 - offset);
		long k_last  = std::min((long)basis.first.rangeJ(index_x.j-1).lastIndex(),(long)index_x.k/2 + offset);
		for (long k=k_first; k<=k_last; ++k) {
			Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
			if (overlap(supp_x,new_supp_x)>0) {
				Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
				//std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}

	PeriodicSupport<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
	//check y-direction
	if (index_y.j==j0_y) {
		for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
			PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
			if (overlap(supp_y,new_supp_y)>0) {
				Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
				//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
				if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
			}
		}
	}
	else {
		long k_first = (long)index_y.k/2 - offset;
		long k_last = (long)index_y.k/2 + offset;
		bool is_wrap = false;
		while(k_first < 1){
			k_first = (long)basis.second.rangeJ(index_y.j-1).lastIndex() + k_first;
			is_wrap = true;
		}
		while(k_last > (long)basis.second.rangeJ(index_y.j-1).lastIndex()){
			k_last = k_last - (long)basis.second.rangeJ(index_y.j-1).lastIndex();
		   is_wrap = true;
		}
		if(k_first == k_last || k_first == k_last+1 || (k_first < k_last && is_wrap)){
			k_first = (long)basis.second.rangeJ(index_y.j-1).firstIndex();
			k_last = (long)basis.second.rangeJ(index_y.j-1).lastIndex();
		}

		if(k_first < k_last){
			for (long k=k_first; k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}
		else{
			for (long k=(long)basis.second.rangeJ(index_y.j-1).firstIndex(); k<=k_last; ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
			for (long k=k_first; k<=(long)basis.second.rangeJ(index_y.j-1).lastIndex(); ++k) {
				PeriodicSupport<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
				if (overlap(supp_y,new_supp_y)>0) {
					Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
					//std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
					if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
				}
			}
		}
	}
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){
			basis_origin.first.getScalingNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}
		else{
			basis_origin.first.getWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){
			basis_origin.second.getScalingNeighborsForScaling((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}
		else{
			basis_origin.second.getWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
		}
		else{
			if((*it).index1.xtype==XBSpline){
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.mra.rangeI(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.mra.rangeI(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
			else{
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
		}
    }
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
	IndexSet<Index2D> indexset = get_indexset(coeffs_origin);
	getStableExpansion(basis_origin, basis_target, indexset, coeffs_target);
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{

	// First we take the counterpart cone
	getCounterpart(basis_origin, basis_target, indexset_origin, coeffs_target);

	// Then we insert all HigherWaveletNeighbours
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.first.getWaveletNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_origin.first,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index1.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_x = basis_target.first.rangeJ(j_s+1).lastIndex();
			k_last_x = basis_target.first.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.first.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				for(long k = basis_origin.first.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				assert(j_x == (*it).index1.j+1);
			}
		}
		else{
			basis_origin.first.getHigherWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																	j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j+1);

		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.second.getWaveletNeighborsForScaling((*it).index2.j, (*it).index2.k,basis_origin.second,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index2.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_y = basis_target.second.rangeJ(j_s+1).lastIndex();
			k_last_y = basis_target.second.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.second.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				for(long k = basis_origin.second.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				assert(j_y == (*it).index2.j+1);
			}
		}
		else{
			basis_origin.second.getHigherWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																	j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j+1);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
		else{
			for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
			for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
    }
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_woMixedHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
	IndexSet<Index2D> indexset = get_indexset(coeffs_origin);
	getStableExpansion_woMixedHW(basis_origin, basis_target, indexset, coeffs_target);
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_woMixedHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{

	// First we take the counterpart cone
	getCounterpart(basis_origin, basis_target, indexset_origin, coeffs_target);

	// Then we insert all HigherWaveletNeighbours
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.first.getWaveletNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_origin.first,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index1.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_x = basis_target.first.rangeJ(j_s+1).lastIndex();
			k_last_x = basis_target.first.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.first.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				for(long k = basis_origin.first.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				assert(j_x == (*it).index1.j+1);
			}
		}
		else{
			basis_origin.first.getHigherWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																	j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j+1);

		}

		// Insert all higher neighbours in x-direction into the indexset
		// (with original y-component)
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
		}
		else{
			for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
			for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
		}


		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.second.getWaveletNeighborsForScaling((*it).index2.j, (*it).index2.k,basis_origin.second,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index2.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_y = basis_target.second.rangeJ(j_s+1).lastIndex();
			k_last_y = basis_target.second.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.second.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				for(long k = basis_origin.second.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				assert(j_y == (*it).index2.j+1);
			}
		}
		else{
			basis_origin.second.getHigherWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																	j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j+1);
		}

		// Insert all higher neighbours in y-direction into the indexset
		// (with original x-component)
		if(k_first_y < k_last_y){
			for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
				Index1D ind_y(j_y, k_new_y, XWavelet);
				completeMultiTree(basis_target, Index2D((*it).index1, ind_y), coeffs_target, 0, true);
			}
		}
		else{
			for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
				Index1D ind_y(j_y, k_new_y, XWavelet);
				completeMultiTree(basis_target, Index2D((*it).index1, ind_y), coeffs_target, 0, true);
			}
			for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
				Index1D ind_y(j_y, k_new_y, XWavelet);
				completeMultiTree(basis_target, Index2D((*it).index1, ind_y), coeffs_target, 0, true);
			}
		}
    }
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_onlyTemporalHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
	IndexSet<Index2D> indexset = get_indexset(coeffs_origin);
	getStableExpansion_onlyTemporalHW(basis_origin, basis_target, indexset, coeffs_target);
}

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_onlyTemporalHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{

	// First we take the counterpart cone
	getCounterpart(basis_origin, basis_target, indexset_origin, coeffs_target);

	// Then we insert all HigherWaveletNeighbours in x direction (= time)
    long k_first_x, k_last_x;
    int j_x;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.first.getWaveletNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_origin.first,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index1.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_x = basis_target.first.rangeJ(j_s+1).lastIndex();
			k_last_x = basis_target.first.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.first.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				for(long k = basis_origin.first.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				assert(j_x == (*it).index1.j+1);
			}
		}
		else{
			basis_origin.first.getHigherWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																	j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j+1);

		}

		// Insert all higher neighbours in x-direction into the indexset
		// (with original y-component)
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
		}
		else{
			for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
			for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);
				completeMultiTree(basis_target, Index2D(ind_x, (*it).index2), coeffs_target, 0, true);
			}
		}
    }
}


}   // namespace lawa

