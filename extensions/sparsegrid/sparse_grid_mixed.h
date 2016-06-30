# include <string>

namespace webbur
{
  void sparse_grid_mixed_index ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[], 
    FLENS_DEFAULT_INDEXTYPE point_num, FLENS_DEFAULT_INDEXTYPE point_total_num, FLENS_DEFAULT_INDEXTYPE sparse_unique_index[], 
    FLENS_DEFAULT_INDEXTYPE sparse_order[], FLENS_DEFAULT_INDEXTYPE sparse_index[] );

  void sparse_grid_mixed_point ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[], 
    double alpha[], double beta[], FLENS_DEFAULT_INDEXTYPE point_num, FLENS_DEFAULT_INDEXTYPE sparse_order[], 
    FLENS_DEFAULT_INDEXTYPE sparse_index[], double sparse_point[] );

  FLENS_DEFAULT_INDEXTYPE sparse_grid_mixed_size ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[], double alpha[],
    double beta[], double tol );

  FLENS_DEFAULT_INDEXTYPE sparse_grid_mixed_size_total ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[] );

  void sparse_grid_mixed_unique_index ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[], 
    double alpha[], double beta[], double tol, FLENS_DEFAULT_INDEXTYPE point_num, FLENS_DEFAULT_INDEXTYPE point_total_num,
   FLENS_DEFAULT_INDEXTYPE sparse_unique_index[] );

  void sparse_grid_mixed_weight ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level_max, FLENS_DEFAULT_INDEXTYPE rule[], 
    double alpha[], double beta[], FLENS_DEFAULT_INDEXTYPE point_num, FLENS_DEFAULT_INDEXTYPE point_total_num, 
    FLENS_DEFAULT_INDEXTYPE sparse_unique_index[], double sparse_weight[] );

  void sparse_grid_mixed_write ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE rule[], double alpha[], 
    double beta[], FLENS_DEFAULT_INDEXTYPE point_num, double sparse_weight[], double sparse_point[], 
    std::string file_name );
}
  

