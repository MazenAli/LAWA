#include <fstream>

namespace lawa {

template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename)
{
  typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
  
  std::ofstream data(filename);

	if(!data.good()){
	    std::cerr << "File " << filename << " could not be opened for writing" << std::endl;
	    exit(1);
	}
  data.precision(40);
  data << "# Center_x Center_y Value Xtype1 j1 k1 Xtype2 j2 k2" << std::endl;
  
  const typename Basis::FirstBasisType& basis_x = basis2d.first;
  const typename Basis::SecondBasisType& basis_y = basis2d.second;
  
  for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
      int j1=(*it).first.index1.j, k1=(*it).first.index1.k, j2=(*it).first.index2.j, k2=(*it).first.index2.k;
      XType type1=(*it).first.index1.xtype, type2=(*it).first.index2.xtype;

      //center of the support
      double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
      double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

      data << x << " " << y << " " << std::scientific <<  (*it).second << " " 
           << type1 << " " << j1 << " " << k1 << " "
           << type2 << " " << j2 << " " << k2 << std::endl;
  }
  data.close();
}

template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename)
{
  typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
  
  std::ofstream data(filename);

	if(!data.good()){
	    std::cerr << "File " << filename << " could not be opened for writing" << std::endl;
	    exit(1);
	}
  data.precision(40);
  data << "# Center_x Center_y Value Xtype1 j1 k1 Xtype2 j2 k2" << std::endl;

  const typename Basis::FirstBasisType& basis_x = basis2d.first;
  const typename Basis::SecondBasisType& basis_y = basis2d.second;

  for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
      int j1=(*it).first.index1.j, k1=(*it).first.index1.k, j2=(*it).first.index2.j, k2=(*it).first.index2.k;
      XType type1=(*it).first.index1.xtype, type2=(*it).first.index2.xtype;

      //center of the support
      double x;
      if(basis_x.generator(type1).support(j1,k1).gaplength() > 0){
    	  double h = 0.5*basis_x.generator(type1).support(j1,k1).length();
          x = basis_x.generator(type1).support(j1,k1).li2 + h;
          if(x > 1.){
        	  x -= 1.;
          }
      }
      else{
          x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
      }

      double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);


      data << x << " " << y << " " << std::scientific <<  (*it).second << " "
           << type1 << " " << j1 << " " << k1 << " "
           << type2 << " " << j2 << " " << k2 << std::endl;
  }
  data.close();
}

template <typename T>
void
readCoeffVector2D(Coefficients<Lexicographical,T,Index2D>&coeff, const char* filename, bool append){
  
  if(append == false){
    coeff.clear();
  }
  
  typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
  
  T x,y, coeff_value;
  int type1, type2;
  int j1, k1, j2, k2;
  std::ifstream data(filename);
  if(data.is_open() == false){
    std::cerr << "File " << filename << " could not be opened for reading" << std::endl;
    exit(1);
  }
  
  char comment[256];
  data.getline (comment,256);
   
  while(!data.eof()){
    data >> x >> y >> coeff_value >> type1 >> j1 >> k1 >> type2 >> j2 >> k2;
    Index1D ind1(j1, k1, (XType)type1);
    Index1D ind2(j2, k2, (XType)type2);
    Index2D ind(ind1, ind2);
    coeff.insert(val_type(ind, coeff_value)); 
  }
}

template <typename T, typename Basis>
void
saveIndexSet2D(const IndexSet<Index2D> &indexset, const Basis &basis2d, const char* filename)
{
  typedef typename IndexSet<Index2D>::const_iterator const_ind_it;
  
  std::ofstream data(filename);
  data << "# Center_x Center_y Xtype1 j1 k1 Xtype2 j2 k2" << std::endl;
  
  const typename Basis::FirstBasisType& basis_x = basis2d.first;
  const typename Basis::SecondBasisType& basis_y = basis2d.second;
  
  for (const_ind_it it = indexset.begin(); it != indexset.end(); ++it) {
      int j1=(*it).index1.j, k1=(*it).index1.k, j2=(*it).index2.j, k2=(*it).index2.k;
      XType type1=(*it).index1.xtype, type2=(*it).index2.xtype;

      //center of the support
      double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
      double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

      data << x << " " << y  << " " 
           << type1 << " " << j1 << " " << k1 << " "
           << type2 << " " << j2 << " " << k2 << std::endl;
  }
  data.close();
  
}

template <typename T>
void
readIndexSet2D(IndexSet<Index2D>& indexset, const char* filename, bool append)
{
  
  if(append == false){
    indexset.clear();
  }
    
  T x,y;
  int type1, type2;
  int j1, k1, j2, k2;
  std::ifstream data(filename);
  if(data.is_open() == false){
    std::cerr << "File " << filename << " could not be opened for reading" << std::endl;
    exit(1);
  }
  
  char comment[256];
  data.getline (comment,256);
   
  while(!data.eof()){
    data >> x >> y >> type1 >> j1 >> k1 >> type2 >> j2 >> k2;
    Index1D ind1(j1, k1, (XType)type1);
    Index1D ind2(j2, k2, (XType)type2);
    Index2D ind(ind1, ind2);
    indexset.insert(ind); 
  }
}

} // namespace lawa
