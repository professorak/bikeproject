//#include <Rcpp.h>
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>     /* atoi */
#include <assert.h>
#include <ctime>

using namespace Rcpp;
using namespace std;
using namespace arma;

typedef unsigned int uint;


void print_vec(colvec);
void print_vec(uvec);
void print_vec(vector<string>);
void print_vec(vector<int>);

/*  functions for cartesian product
 * 
 * 
 */
// Cartesion product of vector of vectors

#include <vector>
#include <iostream>
#include <iterator>

// Types to hold vector-of-ints (Vi) and vector-of-vector-of-ints (Vvi)
typedef std::vector<int> Vi;
typedef std::vector<Vi> Vvi;
typedef std::vector<double> Vd;
typedef std::vector<Vd> Vvd;
// Just for the sample -- populate the intput data set
Vvi build_input() {
  Vvi vvi;
  
  for(int i = 0; i < 3; i++) {
    Vi vi;
    for(int j = 0; j < 3; j++) {
      vi.push_back(i*10+j);
    }
    vvi.push_back(vi);
  }
  return vvi;
}

// just for the sample -- print the data sets
std::ostream&
  operator<<(std::ostream& os, const Vi& vi)
{
    os << "(";
    std::copy(vi.begin(), vi.end(), std::ostream_iterator<int>(os, ", "));
    os << ")";
    return os;
  }
std::ostream&
  operator<<(std::ostream& os, const Vvi& vvi)
{
    os << "(\n";
    for(Vvi::const_iterator it = vvi.begin();
        it != vvi.end();
        it++) {
      os << "  " << *it << "\n";
    }
    os << ")";
    return os;
  }
  std::ostream&
  operator<<(std::ostream& os, const Vd& vd)
{
    os << "(";
    std::copy(vd.begin(), vd.end(), std::ostream_iterator<double>(os, ", "));
    os << ")";
    return os;
  }
std::ostream&
  operator<<(std::ostream& os, const Vvd& vvd)
{
    os << "(\n";
    for(Vvd::const_iterator it = vvd.begin();
        it != vvd.end();
        it++) {
      os << "  " << *it << "\n";
    }
    os << ")";
    return os;
  }


// recursive algorithm to to produce cart. prod.
// At any given moment, "me" points to some Vi in the middle of the
// input data set. 
//   for int i in *me:
  //      add i to current result
//      recurse on next "me"
// 
  void cart_product(
    Vvi& rvvi,  // final result
    Vi&  rvi,   // current result 
    Vvi::const_iterator me, // current input
    Vvi::const_iterator end) // final input
{
  if(me == end) {
    // terminal condition of the recursion. We no longer have
    // any input vectors to manipulate. Add the current result (rvi)
    // to the total set of results (rvvvi).
    rvvi.push_back(rvi);
    return;
  }
  
  // need an easy name for my vector-of-ints
  const Vi& mevi = *me;
  for(Vi::const_iterator it = mevi.begin();
      it != mevi.end();
      it++) {
    // final rvi will look like "a, b, c, ME, d, e, f"
    // At the moment, rvi already has "a, b, c"
    rvi.push_back(*it);  // add ME
    cart_product(rvvi, rvi, me+1, end); //add "d, e, f"
    rvi.pop_back(); // clean ME off for next round
  }
}
  void cart_product(
    Vvd& rvvi,  // final result
    Vd&  rvi,   // current result 
    Vvd::const_iterator me, // current input
    Vvd::const_iterator end) // final input
{
  if(me == end) {
    // terminal condition of the recursion. We no longer have
    // any input vectors to manipulate. Add the current result (rvi)
    // to the total set of results (rvvvi).
    rvvi.push_back(rvi);
    return;
  }
  
  // need an easy name for my vector-of-ints
  const Vd& mevi = *me;
  for(Vd::const_iterator it = mevi.begin();
      it != mevi.end();
      it++) {
    // final rvi will look like "a, b, c, ME, d, e, f"
    // At the moment, rvi already has "a, b, c"
    rvi.push_back(*it);  // add ME
    cart_product(rvvi, rvi, me+1, end); //add "d, e, f"
    rvi.pop_back(); // clean ME off for next round
  }
}
//// sample only, to drive the cart_product routine.
//int main_test() {
//  Vvi input(build_input());
//  std::cout << input << "\n";
//  
//  Vvi output;
//  Vi outputTemp;
//  cart_product(output, outputTemp, input.begin(), input.end());
//  std::cout << output << "\n";
//}


/*
 * */



colvec latlondistance(colvec,colvec,double,double);
vector<int>  &split(const std::string &s, char delim, vector<int> &elems);
vector<int>  split(const std::string &s, char delim);
vector<uint>  which_r(vector<int> j_loc_st, vector<int> st_point_list);
uvec  which_r_str(vector<string> j_loc_st, vector<string> st_point_list);
vector<string> covert_row_str(umat mat_st_state);
vector<string> unique_str(vector<string> myvector);
uvec unique_idx(uvec myvector);
vector<mat>  compute_prob(uint i, mat station_data, NumericMatrix xpoints, uint wdclat1_col, uint wdclon1_col,
    uint pointslat1_col, uint pointslon1_col, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, NumericVector xv0_vec,
    uint points_density_col);
vector<mat>  compute_prob_2(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, double pointsden_i);
mat  compute_prob_theta(uint i, mat station_data, NumericMatrix xpoints, uint wdclat1_col, uint wdclon1_col,
    uint pointslat1_col, uint pointslon1_col, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, NumericVector xv0_vec,
    uint points_density_col);
mat  compute_prob_theta_2(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, double pointsden_i);
rowvec  compute_prob_unweighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec);
mat compute_hessian_delta_sq(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index);
double compute_hessian_beta1_sq(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index);
mat compute_hessian_theta1_sq_weighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index, uint xtheta1_size,
    double point_density_i, rowvec points_den_covariates);
mat compute_hessian_theta1_delta_weighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index, uint xtheta1_size,
    double point_density_i, rowvec points_den_covariates);
void construct_mat_st_state_str_unq(uint i, colvec &lat1, colvec &lon1,double pointslat1_i,double pointslon1_i,
  std::vector<string> &wdc_sto_state_local, std::vector<string> &wdc_local_stations,
  std::vector<string> &points_local_stations,
  vector<int> &st_point_list_org, vector<int> &st_point_list, uvec &list_obs, umat &mat_st_state, 
  imat &obs_st_state, uvec &st_point_list_uvec, vector<string> &mat_st_state_str, 
  vector<string> &mat_st_state_str_unq, uvec &station_id_index_r );
void construct_a_prob_delta_vec(vector< vector<int> > &a_vec, vector< vector<double> > &prob_a,
  vector< vector<double> > &delta_a, imat &obs_st_state, umat &mat_st_state, vector<string> &mat_st_state_str, 
  vector<string> &str_vec, vector<int> &st_point_list, colvec &xdeltain,
  urowvec &col_na, rowvec &delta_avg, urowvec &mat_st_state_row, uvec &station_id_index_r, 
  colvec &wdcobswt);

// [[Rcpp::export]]
SEXP eval_lambda_delta_list_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints(points);
    uint points_density_col = 2;
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
    
    colvec lambda_t(xwdcMergedday.nrow(),fill::zeros); 
    arma::mat grad_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    
    for(uint i=0;i<xpoints.nrow();i++) {
    //  cout << "point no" << i << endl;
    //for(uint i=2;i<3;i++) {
        
        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              vector<mat> ret = compute_prob(i, station_data, xpoints, wdclat1_col, wdclon1_col, pointslat1_col, 
              pointslon1_col, beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, points_density_col); 
                            
              rowvec lambda_st_t = ret[0];              
              lambda_t(a_vec[k][l]) +=  lambda_st_t(k);

              mat util_grad = ret[1];              

              rowvec grad_temp = util_grad.row(k);

              vector< vector<double> > prob_a_temp = prob_a;
              vector<double> temp_vec1(1); temp_vec1[0]=1;
              prob_a_temp[k] = temp_vec1;
              vector< vector<int> > a_vec_temp = a_vec;
              vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
              a_vec_temp[k] = temp_vec2;
              vector< vector<double> > grad_temp_list = prob_a_temp;              
              vector<int> a_vec_unlisted;
              vector<double> grad_temp_unlisted;
              for(uint m=0; m <prob_a_temp.size(); ++m) {                
                std::transform(grad_temp_list[m].begin(), grad_temp_list[m].end(), 
                  grad_temp_list[m].begin(), std::bind1st(std::multiplies<double>(),grad_temp[m]));                
                a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                grad_temp_unlisted.insert(grad_temp_unlisted.end(),grad_temp_list[m].begin(),grad_temp_list[m].end());
              }
           
              rowvec grad_temp_unlisted_rowvec = conv_to<rowvec>::from(grad_temp_unlisted);
              uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
              uvec rowno(1); rowno(0) = a_vec_temp[k][0];
              grad_t(rowno,a_vec_unlisted_uvec)  += mat(grad_temp_unlisted_rowvec);
            }
            }
          }

        }
    }//end of points loop      
    colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    mat wdcobswt_mat = repmat(wdcobswt_colvec,1,xwdcMergedday.nrow());
    lambda_t = lambda_t%wdcobswt_colvec;    
    grad_t = grad_t%wdcobswt_mat;    
    
    mat obj_ret = join_rows(lambda_t,grad_t); 
    return(wrap(obj_ret));  
      
}
      
// [[Rcpp::export]]
SEXP eval_lambda_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints(points);
    uint points_density_col = 2;
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
    
    colvec lambda_t(xwdcMergedday.nrow(),fill::zeros); 
    //arma::mat grad_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])

    for(uint i=0;i<xpoints.nrow();i++) {

        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
               
              vector<mat> ret = compute_prob(i, station_data, xpoints, wdclat1_col, wdclon1_col, pointslat1_col, 
              pointslon1_col, beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, points_density_col); 
                            
              rowvec lambda_st_t = ret[0];              

              lambda_t(a_vec[k][l]) +=  lambda_st_t(k);
//              mat util_grad = ret[1];              
//
//              rowvec grad_temp = util_grad.row(k);
//
//              vector< vector<double> > prob_a_temp = prob_a;
//              vector<double> temp_vec1(1); temp_vec1[0]=1;
//              prob_a_temp[k] = temp_vec1;
//              vector< vector<int> > a_vec_temp = a_vec;
//              vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
//              a_vec_temp[k] = temp_vec2;
//              vector< vector<double> > grad_temp_list = prob_a_temp;              
//              vector<int> a_vec_unlisted;
//              vector<double> grad_temp_unlisted;
//              for(uint m=0; m <prob_a_temp.size(); ++m) {                
//                std::transform(grad_temp_list[m].begin(), grad_temp_list[m].end(), 
//                  grad_temp_list[m].begin(), std::bind1st(std::multiplies<double>(),grad_temp[m]));                
//                a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
//                grad_temp_unlisted.insert(grad_temp_unlisted.end(),grad_temp_list[m].begin(),grad_temp_list[m].end());
//              }
//           
//              rowvec grad_temp_unlisted_rowvec = conv_to<rowvec>::from(grad_temp_unlisted);
//              uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
//              uvec rowno(1); rowno(0) = a_vec_temp[k][0];
//              grad_t(rowno,a_vec_unlisted_uvec)  += mat(grad_temp_unlisted_rowvec);
            }
            }
          }

        }
    }//end of points loop      

    // colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    // mat wdcobswt_mat = repmat(wdcobswt_colvec,1,xwdcMergedday.nrow());
    // lambda_t = lambda_t%wdcobswt_colvec;    
    //grad_t = grad_t%wdcobswt_mat;    
    
    //mat obj_ret = join_rows(lambda_t,grad_t); 
    return(wrap(lambda_t));  
      
}

// [[Rcpp::export]]
SEXP eval_lambda_multiple_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
    
    arma::mat lambda_t(xwdcMergedday.nrow(), xpoints.n_cols-2, fill::zeros);
    //arma::mat grad_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {

        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.

          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              rowvec lambda_st_t = compute_prob_unweighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec); 
                            
              lambda_t.row(a_vec[k][l]) +=  lambda_st_t(k)*xpoints( i, span(min_points_col,max_points_col));

            }
            }
          }

        }
    }//end of points loop  

    colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    mat wdcobswt_mat = repmat(wdcobswt_colvec,1,xpoints.n_cols-2);
    lambda_t = lambda_t % wdcobswt_mat;    
    //grad_t = grad_t%wdcobswt_mat;    
    
    //mat obj_ret = join_rows(lambda_t,grad_t); 
    return(wrap(lambda_t));  
      
}

// [[Rcpp::export]]
SEXP eval_grad_lambda_theta_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints(points);
    uint points_density_col = 2;
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
    
    arma::mat grad_t(xwdcMergedday.nrow(),2,fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    
    for(uint i=0;i<xpoints.nrow();i++) {

        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              mat util_grad = compute_prob_theta(i, station_data, xpoints, wdclat1_col, wdclon1_col, pointslat1_col, 
              pointslon1_col, beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, points_density_col); 

              rowvec grad_temp = util_grad.row(k);

              uvec rowno(1); rowno(0) = a_vec[k][l];
              grad_t.rows(rowno)  += util_grad.row(k);
            }
            }
          }

        }
    }//end of points loop      
    colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    mat wdcobswt_mat = repmat(wdcobswt_colvec,1,grad_t.n_cols);
    grad_t = grad_t%wdcobswt_mat;    
    
    return(wrap(grad_t));  
      
}

void print_vec(colvec z) {
  /* Print st_id_index vector to console */
    copy(z.begin(), z.end(), ostream_iterator<double>(cout, " "));
  cout << "\n";
} 
void print_vec(uvec z) {
  /* Print st_id_index vector to console */
    copy(z.begin(), z.end(), ostream_iterator<uword>(cout, " "));    
  cout << "\n";
} 
void print_vec(vector<string> z) {
  /* Print st_id_index vector to console */
    copy(z.begin(), z.end(), ostream_iterator<string>(cout, " "));    
  cout << "\n";
} 
void print_vec(vector<int> z) {
  /* Print st_id_index vector to console */
    copy(z.begin(), z.end(), ostream_iterator<int>(cout, " "));    
  cout << "\n";
} 

colvec latlondistance(colvec lat1, colvec lon1, double lat2, double lon2) {
  //check Rcpp sugar why this cant be done on Rccp vectors
  lat1 *= 3.14/180;
  lon1 *= 3.14/180;
  lat2 *= 3.14/180;
  lon2 *= 3.14/180;
  double R = 6371;
  
  colvec disx = (lon2-lon1) % cos((lat1+lat2)/2);
  colvec disy = (lat2-lat1);
  colvec dis_v  = R*sqrt(disx%disx + disy%disy);
  return(dis_v);
  
}

//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  //    std::stringstream ss(s);
  //    std::string item;
  //    while (std::getline(ss, item, delim)) {
    //        elems.push_back(item);
    //    }
  //    return elems;
  //}
//
  //
  //std::vector<std::string> split(const std::string &s, char delim) {
    //    std::vector<std::string> elems;
    //    split(s, delim, elems);
    //    return elems;
    //}

vector<int> &split(const std::string &s, char delim, vector<int> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(atoi(item.c_str()));
  }
  return elems;
}


vector<int> split(const std::string &s, char delim) {
  vector<int> elems;
  split(s, delim, elems);
  return elems;
}

//returns which element from second list do each element of first list matches to
vector<uint> which_r(vector<int> j_loc_st, vector<int> st_point_list) {
  vector<uint> elems;
  for(uint i=0;i<j_loc_st.size() ;i++) {
    if(find(st_point_list.begin(), st_point_list.end(), j_loc_st[i])!=st_point_list.end()) {
      elems.push_back(i);
      //cout << i << endl;
    }
  }
  return elems;
}
uvec which_r_str(vector<string> j_loc_st, vector<string> st_point_list) {
  uvec elems(j_loc_st.size());
  uint count = 0;
  for(uint i=0;i<j_loc_st.size() ;i++) {
    if(find(st_point_list.begin(), st_point_list.end(), j_loc_st[i])!=st_point_list.end()) {
      elems(count)= i;
      count++;
    }
  }
  
  elems.resize(count);
  return elems;
}



vector<string> covert_row_str(umat mat_st_state) {
  vector<string> mat_st_state_str;
  for(uint k=0; k < mat_st_state.n_rows; ++k) {
    std::ostringstream oss;  
    urowvec v = mat_st_state.row(k);
    copy(v.begin(), v.end(), ostream_iterator<uword>(oss, ""));    
    //cout << oss.str() << endl;   
    mat_st_state_str.push_back(oss.str());
  }
  return mat_st_state_str;
}

vector<string> unique_str(vector<string> myvector) {
  std::sort (myvector.begin(), myvector.end()); 
  std::vector<string>::iterator it;
  it = std::unique (myvector.begin(), myvector.end());   // 10 20 30 20 10 ?  ?  ?  ?    
  myvector.resize( std::distance(myvector.begin(),it) ); // 10 20 30 20 10      
  return myvector;
}

uvec unique_idx(uvec myvector) {
  //return the indexes of unique elements of myvector
  vector<uword> uvector = conv_to< vector<uword> >::from (myvector);      
  std::sort (uvector.begin(), uvector.end()); 
  
  vector<uword> u_idx;
  u_idx.push_back(0);
  
  uint result = uvector[0];
  for(uint j=1; j<uvector.size(); ++j) {
    if(result!=uvector[j]) {
      u_idx.push_back(j);
      result=uvector[j];
    }
  }
  
  return conv_to< uvec >::from(u_idx);
}

//http://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors      
//#include <algorithm>
//#include <functional>
//
//template <typename T>
//std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
//{
//    assert(a.size() == b.size());
//
//    std::vector<T> result;
//    result.reserve(a.size());
//
//    std::transform(a.begin(), a.end(), b.begin(), 
//                   std::back_inserter(result), std::plus<T>());
//    return result;
//}    
//      

vector<mat>  compute_prob(uint i, mat station_data, NumericMatrix xpoints, uint wdclat1_col, uint wdclon1_col,
    uint pointslat1_col, uint pointslon1_col, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint points_density_col) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), xpoints(i,pointslat1_col), xpoints(i,pointslon1_col)));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          rowvec lambda_st_t(no_t_st,fill::zeros);
          mat util_grad(no_t_st,no_t_st,fill::zeros);                    
          mat A(no_t_st,no_t_st,fill::zeros);
          
          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
          //for(int m=0; m<1; ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              rowvec util_prob_t =  util/denutil_t;
              
              lambda_st_t += util_prob_t;
              for(uint l=0; l<no_t_st; ++l) {
                A( l, l ) = util_prob_t(l);          
              }
              mat B = repmat(util_prob_t,no_t_st,1);
              mat C = repmat(vectorise( util_prob_t, 0 ),1,no_t_st);              
              util_grad += A - (B%C);
          }
          lambda_st_t *= (xpoints(i,points_density_col)/xv0_vec.size());
          util_grad *= (xpoints(i,points_density_col)/xv0_vec.size());
          
        
//          lambda_st_t = lambda_st_t % prob_row;
//          mat prob_row_mat = repmat(vectorise( prob_row, 0 ),1,no_t_st);          
//          util_grad = util_grad % prob_row_mat;

          vector<mat> obj_ret(2);
          obj_ret[0] = lambda_st_t;
          obj_ret[1] = util_grad;
          return((obj_ret));  
}

vector<mat>  compute_prob_2(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, double pointsden_i) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          rowvec lambda_st_t(no_t_st,fill::zeros);
          mat util_grad(no_t_st,no_t_st,fill::zeros);                    
          mat A(no_t_st,no_t_st,fill::zeros);
          
          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
          //for(int m=0; m<1; ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              rowvec util_prob_t =  util/denutil_t;
              
              lambda_st_t += util_prob_t;
              for(uint l=0; l<no_t_st; ++l) {
                A( l, l ) = util_prob_t(l);          
              }
              mat B = repmat(util_prob_t,no_t_st,1);
              mat C = repmat(vectorise( util_prob_t, 0 ),1,no_t_st);              
              util_grad += A - (B%C);
          }
          lambda_st_t *= (pointsden_i/xv0_vec.size());
          util_grad *= (pointsden_i/xv0_vec.size());
          
        
//          lambda_st_t = lambda_st_t % prob_row;
//          mat prob_row_mat = repmat(vectorise( prob_row, 0 ),1,no_t_st);          
//          util_grad = util_grad % prob_row_mat;

          vector<mat> obj_ret(2);
          obj_ret[0] = lambda_st_t;
          obj_ret[1] = util_grad;
          return((obj_ret));  
}

rowvec  compute_prob_unweighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          rowvec lambda_st_t(no_t_st,fill::zeros);
          // mat util_grad(no_t_st,no_t_st,fill::zeros);                    
          // mat A(no_t_st,no_t_st,fill::zeros);
          
          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
          //for(int m=0; m<1; ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              rowvec util_prob_t =  util/denutil_t;
              
              lambda_st_t += util_prob_t;
              // for(uint l=0; l<no_t_st; ++l) {
              //   A( l, l ) = util_prob_t(l);          
              // }
              // mat B = repmat(util_prob_t,no_t_st,1);
              // mat C = repmat(vectorise( util_prob_t, 0 ),1,no_t_st);              
              // util_grad += A - (B%C);
          }
          lambda_st_t *= (1/xv0_vec.size());
          // util_grad *= (1/xv0_vec.size());
          
        
//          lambda_st_t = lambda_st_t % prob_row;
//          mat prob_row_mat = repmat(vectorise( prob_row, 0 ),1,no_t_st);          
//          util_grad = util_grad % prob_row_mat;

//          vector<mat> obj_ret(2);
//          obj_ret[0] = lambda_st_t;
//          obj_ret[1] = util_grad;
          return((lambda_st_t));  
}

mat  compute_prob_theta(uint i, mat station_data, NumericMatrix xpoints, uint wdclat1_col, uint wdclon1_col,
    uint pointslat1_col, uint pointslon1_col, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint points_density_col) {
  
    rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
      station_data.col(wdclon1_col), xpoints(i,pointslat1_col), xpoints(i,pointslon1_col)));                    
    
    rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
    double den_util = sum(util);
    uint no_t_st = util.size();          
    rowvec grad_beta1(no_t_st,fill::zeros);
    rowvec grad_sigma0(no_t_st,fill::zeros);
    
    for(int m=0; m<xv0_vec.size(); ++m) {
    //for(int m=0; m<1; ++m) {
        double out = exp(-xv0_vec(m)*sigma0);
        double denutil_t = den_util+out;        
        rowvec util_prob_t =  util/denutil_t;
          
        rowvec disP = util_prob_t%station_data_dis_vIdx;
        double disP_sum = sum( disP);
        rowvec disP_sum_vec(no_t_st);
        disP_sum_vec.fill(disP_sum);            
        //#gradient wrt beta1 is given as dis_ij*Pij-Pij(\sum_k dis_ik*Pik) for a t
        grad_beta1 += disP - util_prob_t%disP_sum_vec;
        //#gradient wrt sigma0 is given as vio*pij*pio
        double prob0 = 1-sum( util_prob_t);              
        grad_sigma0 += xv0_vec(m)*prob0*util_prob_t;  

    }
    grad_beta1 *= (xpoints(i,points_density_col)/xv0_vec.size());
    grad_sigma0 *= (xpoints(i,points_density_col)/xv0_vec.size());

//    grad_beta1 = grad_beta1 % prob_row;
//    grad_sigma0 = grad_sigma0 % prob_row;
    
    mat obj_ret = join_rows(trans(grad_beta1),trans(grad_sigma0)); 
    return(obj_ret);  

}

mat  compute_prob_theta_2(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, double pointsden_i) {
  
    rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
      station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
    
    rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
    double den_util = sum(util);
    uint no_t_st = util.size();          
    rowvec grad_beta1(no_t_st,fill::zeros);
    rowvec grad_sigma0(no_t_st,fill::zeros);
    
    for(int m=0; m<xv0_vec.size(); ++m) {
    //for(int m=0; m<1; ++m) {
        double out = exp(-xv0_vec(m)*sigma0);
        double denutil_t = den_util+out;        
        rowvec util_prob_t =  util/denutil_t;
          
        rowvec disP = util_prob_t%station_data_dis_vIdx;
        double disP_sum = sum( disP);
        rowvec disP_sum_vec(no_t_st);
        disP_sum_vec.fill(disP_sum);            
        //#gradient wrt beta1 is given as dis_ij*Pij-Pij(\sum_k dis_ik*Pik) for a t
        grad_beta1 += disP - util_prob_t%disP_sum_vec;
        //#gradient wrt sigma0 is given as vio*pij*pio
        double prob0 = 1-sum( util_prob_t);              
        grad_sigma0 += xv0_vec(m)*prob0*util_prob_t;  

    }
    grad_beta1 *= (pointsden_i/xv0_vec.size());
    grad_sigma0 *= (pointsden_i/xv0_vec.size());

//    grad_beta1 = grad_beta1 % prob_row;
//    grad_sigma0 = grad_sigma0 % prob_row;
    
    mat obj_ret = join_rows(trans(grad_beta1),trans(grad_sigma0)); 
    return(obj_ret);  

}

/********************/
//hessian implementation
// [[Rcpp::export]]
SEXP eval_hessian_lambda_delta_sq_cpp(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in,
         SEXP lambda_multiplers_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    NumericVector lambda_multiplers(lambda_multiplers_in);
    assert(lambda_multiplers.size()==xwdcMergedday.nrow());
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
        
    //arma::mat lambda_t(xwdcMergedday.nrow(), xpoints.n_cols-2, fill::zeros);
    arma::mat hessian_delta_sq_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {
        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              // rowvec lambda_st_t = compute_prob_unweighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
              //   xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              // xv0_vec); 
                            
              // lambda_t.row(a_vec[k][l]) +=  lambda_st_t(k)*xpoints( i, span(min_points_col,max_points_col));
              if(lambda_multiplers(a_vec[k][l])==0) continue;

              mat hessian_delta_sq_kl = compute_hessian_delta_sq(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, k);
              //multiply with point weight and lambda_multiplers_in(a_vec[k][l])
              
              hessian_delta_sq_kl *=  wdcobswt(a_vec[k][l])*lambda_multiplers(a_vec[k][l]) * xpoints( i, min_points_col); 

              //need to expand hessian_delta_sq_kl to reflect gradients wrt
              //a_vec columns to reflect the deltaaveraged gradients.
              //test in a seperate Rcpp file how to repeat rows and columns and then 
              //multiply rows and columns with prob_a values.
              //creating version of a_vec and prob_a which have the focal station 
              //with only one entry and rest of the stations with actual list.
              vector< vector<double> > prob_a_temp = prob_a;
              vector<double> temp_vec1(1); temp_vec1[0]=1;
              prob_a_temp[k] = temp_vec1;
              vector< vector<int> > a_vec_temp = a_vec;
              vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
              a_vec_temp[k] = temp_vec2;
              //cout << "simplify above lines, there should be way of direclty assigning\
              //instead of creating temp vecs" << endl; 
              //unlisting above lists
              vector<int> a_vec_unlisted;
              vector<double> prob_a_unlisted;
              //create list of hessian_delta_sq_kl indexes to expand
              uvec hessian_expand_index;

              for(uint m=0; m <prob_a_temp.size(); ++m) {  
                uvec hessian_expand_index_temp(prob_a_temp[m].size());
                hessian_expand_index_temp.fill(m);
                hessian_expand_index.insert_rows( hessian_expand_index.size(), hessian_expand_index_temp ); 
                a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                prob_a_unlisted.insert(prob_a_unlisted.end(),prob_a_temp[m].begin(),prob_a_temp[m].end());                
              }              
              mat weights_mat(prob_a_unlisted.size(),prob_a_unlisted.size(),fill::zeros);
              weights_mat.diag()  = conv_to<vec>::from(prob_a_unlisted);
              mat hessian_delta_sq_kl_expanded = hessian_delta_sq_kl.rows(hessian_expand_index);
              hessian_delta_sq_kl_expanded = hessian_delta_sq_kl_expanded.cols(hessian_expand_index);
              hessian_delta_sq_kl_expanded = weights_mat * hessian_delta_sq_kl_expanded * weights_mat;
              uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
              hessian_delta_sq_t(a_vec_unlisted_uvec,a_vec_unlisted_uvec) += hessian_delta_sq_kl_expanded;
              
            }
            }
          }

        }
    }//end of points loop  

    // colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    // mat wdcobswt_mat = repmat(wdcobswt_colvec,1,xpoints.n_cols-2);
    // lambda_t = lambda_t % wdcobswt_mat;    
    //grad_t = grad_t%wdcobswt_mat;    
    
    //mat obj_ret = join_rows(lambda_t,grad_t); 
    return(wrap(hessian_delta_sq_t));  
      
}


mat compute_hessian_delta_sq(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          //rowvec lambda_st_t(no_t_st,fill::zeros);
          mat hessian_delta_sq_t(no_t_st,no_t_st,fill::zeros);

          uvec no_focal_indexes(no_t_st,fill::zeros);
          
          //fill  no_focal_indexes with index sequence
          //find more efficient way to do this
          for(uint m=0; m<no_focal_indexes.size(); ++m) {
            no_focal_indexes(m)=m;
          }
          no_focal_indexes.shed_row(focal_station_index);
          

          mat A(no_t_st,no_t_st,fill::zeros);
          urowvec e_f(no_t_st, fill::zeros);
          e_f(focal_station_index) = 1;          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
          //for(int m=0; m<1; ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              
              rowvec util_prob_t =  util/denutil_t;
              
              rowvec util_prob_t_nofocal = util_prob_t;
              
              util_prob_t_nofocal.shed_col(focal_station_index);
               
              mat B_1 = repmat(vectorise( util_prob_t_nofocal, 0),1,no_t_st-1);
              mat B_2 = repmat(util_prob_t_nofocal,no_t_st-1,1);
              mat I = eye<mat>(no_t_st-1,no_t_st-1);
              mat A_2 = -B_1 % (I - 2*B_2);
              A(no_focal_indexes,no_focal_indexes) = A_2;

              rowvec A_1 = (1-2*util_prob_t(focal_station_index)) * (e_f-util_prob_t);
              A.row(focal_station_index) = A_1;
              A.col(focal_station_index) = A_1.t();
              hessian_delta_sq_t += util_prob_t(focal_station_index)*A;
          }
          hessian_delta_sq_t *= (1/xv0_vec.size());
          return((hessian_delta_sq_t));  
}

// [[Rcpp::export]]
SEXP eval_hessian_lambda_beta1_sq_cpp(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in,
         SEXP lambda_multiplers_in) {
 
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    NumericVector lambda_multiplers(lambda_multiplers_in);
    assert(lambda_multiplers.size()==xwdcMergedday.nrow());
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
        
    //arma::mat lambda_t(xwdcMergedday.nrow(), xpoints.n_cols-2, fill::zeros);
    double hessian_beta1_sq=0;
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {
 
         vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              if(lambda_multiplers(a_vec[k][l])==0) continue;

              double hessian_beta1_sq_kl = compute_hessian_beta1_sq(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, k);
              //multiply with point weight and lambda_multiplers_in(a_vec[k][l])
              
              hessian_beta1_sq_kl *=  wdcobswt(a_vec[k][l])*lambda_multiplers(a_vec[k][l]) * xpoints( i, min_points_col); 

              hessian_beta1_sq += hessian_beta1_sq_kl;
              
            }
            }
          }

        }
    }//end of points loop  

    return(wrap(hessian_beta1_sq));  
      
}



double compute_hessian_beta1_sq(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          //rowvec lambda_st_t(no_t_st,fill::zeros);
          double hessian_beta1_sq_t=0;

          uvec no_focal_indexes(no_t_st,fill::zeros);
          
          //fill  no_focal_indexes with index sequence
          //find more efficient way to do this
          for(uint m=0; m<no_focal_indexes.size(); ++m) {
            no_focal_indexes(m)=m;
          }
          no_focal_indexes.shed_row(focal_station_index);
          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              
              rowvec util_prob_t =  util/denutil_t;

              rowvec disP = util_prob_t%station_data_dis_vIdx;
              double disP_sum = sum( disP);

              double val1 = util_prob_t(focal_station_index) *\
                pow(station_data_dis_vIdx(focal_station_index)-disP_sum,2.0); 

              rowvec disP_sum_vec(no_t_st);
              disP_sum_vec.fill(disP_sum);

              double val2 = util_prob_t(focal_station_index)*sum((station_data_dis_vIdx-disP_sum_vec) % disP);
              hessian_beta1_sq_t += val1 - val2;
          }
          hessian_beta1_sq_t *= (1/xv0_vec.size());
          return((hessian_beta1_sq_t));  
}

// [[Rcpp::export]]
SEXP eval_hessian_lambda_theta1_sq_cpp(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in,
         SEXP lambda_multiplers_in) {
 
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    assert(max_points_col-min_points_col==xtheta1.size()-2); //to assert other points covariates supplied correspond to density vector.
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    NumericVector lambda_multiplers(lambda_multiplers_in);
    assert(lambda_multiplers.size()==xwdcMergedday.nrow());
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
        
    //arma::mat lambda_t(xwdcMergedday.nrow(), xpoints.n_cols-2, fill::zeros);
    mat hessian_theta1_sq(xtheta1.size(),xtheta1.size(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {
         vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              if(lambda_multiplers(a_vec[k][l])==0) continue;

              mat hessian_theta1_sq_kl = compute_hessian_theta1_sq_weighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, k, xtheta1.size(), xpoints(i, min_points_col), xpoints(i,span(min_points_col+1,max_points_col)));
              //multiply with observation wt & lambda_multiplers_in(a_vec[k][l])
              
              hessian_theta1_sq_kl *=  wdcobswt(a_vec[k][l])*lambda_multiplers(a_vec[k][l]);
              

              hessian_theta1_sq += hessian_theta1_sq_kl;
              
            }
            }
          }

        }
    }//end of points loop  

    return(wrap(hessian_theta1_sq));  
      
}



mat compute_hessian_theta1_sq_weighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index, uint xtheta1_size,
    double point_density_i, rowvec points_den_covariates) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          //rowvec lambda_st_t(no_t_st,fill::zeros);
          mat hessian_theta1(xtheta1_size,xtheta1_size,fill::zeros);
          double hessian_beta1_sq_t=0;
          double hessian_sigma_sq_t=0;
          double hessian_beta1_thetaden_t=0;  

          uvec no_focal_indexes(no_t_st,fill::zeros);          
          
          for(int m=0; m<xv0_vec.size(); ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              
              rowvec util_prob_t =  util/denutil_t;

              rowvec disP = util_prob_t%station_data_dis_vIdx;
              double disP_sum = sum( disP);

              double val1 = util_prob_t(focal_station_index) *\
                pow(station_data_dis_vIdx(focal_station_index)-disP_sum,2.0); 

              rowvec disP_sum_vec(no_t_st);
              disP_sum_vec.fill(disP_sum);
              
              double val2 = util_prob_t(focal_station_index)*sum((station_data_dis_vIdx-disP_sum_vec) % disP);
              hessian_beta1_sq_t += val1 - val2;

              hessian_beta1_thetaden_t += util_prob_t(focal_station_index)*\
                (station_data_dis_vIdx(focal_station_index)-disP_sum);
          }
          hessian_beta1_sq_t *= (1/xv0_vec.size());

          hessian_beta1_thetaden_t *= (1/xv0_vec.size());

          vec hessian_beta1_thetaden_vec_t(xtheta1_size-2);
          hessian_beta1_thetaden_vec_t.fill(hessian_beta1_thetaden_t);
          //multilpy hessian_beta1_thetaden_vec_t with points_den_covariates
          hessian_beta1_thetaden_vec_t = hessian_beta1_thetaden_vec_t % points_den_covariates.t();

          hessian_theta1(0,0) = hessian_beta1_sq_t * point_density_i;
          hessian_theta1(0,span(2,xtheta1_size-1)) = hessian_beta1_thetaden_vec_t.t();
          hessian_theta1(span(2,xtheta1_size-1),0) = hessian_beta1_thetaden_vec_t;

          return((hessian_theta1));  
}







// [[Rcpp::export]]
SEXP eval_hessian_lambda_theta1_delta_cpp(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in,
         SEXP lambda_multiplers_in) {
 
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    assert(max_points_col-min_points_col==xtheta1.size()-2); //to assert other points covariates supplied correspond to density vector.
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    NumericVector lambda_multiplers(lambda_multiplers_in);
    assert(lambda_multiplers.size()==xwdcMergedday.nrow());
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
        
    mat hessian_theta1_delta(xtheta1.size(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {
        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
           
              if(lambda_multiplers(a_vec[k][l])==0) continue;

              mat hessian_theta1_delta_kl = compute_hessian_theta1_delta_weighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, k, xtheta1.size(), xpoints(i, min_points_col), xpoints(i,span(min_points_col+1,max_points_col)));
              //multiply with observation wt & lambda_multiplers_in(a_vec[k][l])
              hessian_theta1_delta_kl *=  wdcobswt(a_vec[k][l])*lambda_multiplers(a_vec[k][l]);
              
              //need to expand hessian_delta_sq_kl to reflect gradients wrt
              //a_vec columns to reflect the deltaaveraged gradients.
              //test in a seperate Rcpp file how to repeat rows and columns and then 
              //multiply rows and columns with prob_a values.
              //creating version of a_vec and prob_a which have the focal station 
              //with only one entry and rest of the stations with actual list.
              vector< vector<double> > prob_a_temp = prob_a;
              vector<double> temp_vec1(1); temp_vec1[0]=1;
              prob_a_temp[k] = temp_vec1;
              vector< vector<int> > a_vec_temp = a_vec;
              vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
              a_vec_temp[k] = temp_vec2;
              //cout << "simplify above lines, there should be way of direclty assigning\
              //instead of creating temp vecs" << endl; 
              //unlisting above lists
              vector<int> a_vec_unlisted;
              vector<double> prob_a_unlisted;
              //create list of hessian_delta_sq_kl indexes to expand
              uvec hessian_expand_index;

              for(uint m=0; m <prob_a_temp.size(); ++m) {  
                uvec hessian_expand_index_temp(prob_a_temp[m].size());
                hessian_expand_index_temp.fill(m);
                hessian_expand_index.insert_rows( hessian_expand_index.size(), hessian_expand_index_temp ); 
                a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                prob_a_unlisted.insert(prob_a_unlisted.end(),prob_a_temp[m].begin(),prob_a_temp[m].end());                
              }              
              mat weights_mat(prob_a_unlisted.size(),prob_a_unlisted.size(),fill::zeros);
              weights_mat.diag()  = conv_to<vec>::from(prob_a_unlisted);
              mat hessian_theta1_delta_kl_expanded = hessian_theta1_delta_kl.cols(hessian_expand_index);
              //hessian_delta_sq_kl_expanded = hessian_delta_sq_kl_expanded.cols(hessian_expand_index);
              hessian_theta1_delta_kl_expanded = hessian_theta1_delta_kl_expanded * weights_mat;
              uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
              hessian_theta1_delta.cols(a_vec_unlisted_uvec) += hessian_theta1_delta_kl_expanded;
            }
            }
          }

        }
    }//end of points loop  

    return(wrap(hessian_theta1_delta));  
      
}



mat compute_hessian_theta1_delta_weighted(uint i, mat station_data, uint wdclat1_col, uint wdclon1_col,
    double pointslat1_i, double pointslon1_i, double beta1, double sigma0, colvec xdeltain, 
    uvec st_point_list_uvec, rowvec deltain_row, urowvec mat_st_state_row, 
    NumericVector xv0_vec, uint focal_station_index, uint xtheta1_size,
    double point_density_i, rowvec points_den_covariates) {
  
          rowvec station_data_dis_vIdx = conv_to< rowvec >::from(latlondistance(station_data.col(wdclat1_col), 
            station_data.col(wdclon1_col), pointslat1_i, pointslon1_i));                    
          
          rowvec util = exp(beta1*station_data_dis_vIdx + deltain_row)% (mat_st_state_row==0);
          double den_util = sum(util);
          uint no_t_st = util.size();          
          //rowvec lambda_st_t(no_t_st,fill::zeros);
          rowvec hessian_beta1_delta_t(no_t_st,fill::zeros);
          mat hessian_theta1_delta_t(xtheta1_size,no_t_st,fill::zeros);
          rowvec grad_delta(no_t_st,fill::zeros);

          uvec no_focal_indexes(no_t_st,fill::zeros);
          
          //fill  no_focal_indexes with index sequence
          //find more efficient way to do this
          for(uint m=0; m<no_focal_indexes.size(); ++m) {
            no_focal_indexes(m)=m;
          }
          no_focal_indexes.shed_row(focal_station_index);

          for(int m=0; m<xv0_vec.size(); ++m) {
              double out = exp(-xv0_vec(m)*sigma0);
              double denutil_t = den_util+out;        
              
              rowvec util_prob_t =  util/denutil_t;

              rowvec disP = util_prob_t%station_data_dis_vIdx;
              double disP_sum = sum( disP);
              // rowvec disP_sum_vec(no_t_st);
              // disP_sum_vec.fill(disP_sum);

              rowvec val1(no_t_st,fill::zeros);
              val1 = station_data_dis_vIdx; 

              val1 += station_data_dis_vIdx(focal_station_index) - 2*disP_sum;
              
              val1 = val1 % util_prob_t;
              val1 *= -util_prob_t(focal_station_index);
              //remove focal_station_index from val1 as it is incorrect.
              val1.shed_col(focal_station_index);
              hessian_beta1_delta_t(no_focal_indexes) += val1;
              hessian_beta1_delta_t(focal_station_index) += util_prob_t(focal_station_index) * (1-2*util_prob_t(focal_station_index))*\
                (station_data_dis_vIdx(focal_station_index)-disP_sum);

              grad_delta -= util_prob_t(focal_station_index)*util_prob_t;
              grad_delta(focal_station_index) += util_prob_t(focal_station_index);
          }
          grad_delta *= (1/xv0_vec.size());
          hessian_beta1_delta_t *= (1/xv0_vec.size())* point_density_i;
          
          mat hessian_thetaden_delta_t = points_den_covariates.t() * grad_delta;
          assert(hessian_thetaden_delta_t.n_rows==points_den_covariates.size());
          assert(hessian_thetaden_delta_t.n_cols==grad_delta.size());

          hessian_theta1_delta_t.row(0)=hessian_beta1_delta_t;
          hessian_theta1_delta_t.rows(span(2,xtheta1_size-1))=hessian_thetaden_delta_t;

          return((hessian_theta1_delta_t));  
}


void construct_mat_st_state_str_unq(uint i, colvec &lat1, colvec &lon1,double pointslat1_i,double pointslon1_i,
  std::vector<string> &wdc_sto_state_local, std::vector<string> &wdc_local_stations,
  std::vector<string> &points_local_stations,
  vector<int> &st_point_list_org, vector<int> &st_point_list, uvec &list_obs, umat &mat_st_state, 
  imat &obs_st_state, uvec &st_point_list_uvec, vector<string> &mat_st_state_str, 
  vector<string> &mat_st_state_str_unq, uvec &station_id_index_r ) {

  colvec dis_vIdx = latlondistance(lat1, lon1, 
                                   pointslat1_i, pointslon1_i);
  // ref for find - http://arma.sourceforge.net/docs.html#find
  //        uvec list_obs = find(dis_vIdx <=xmax_walking_dis);
  //        if(list_obs.size()==0) continue;
  //        print_vec(list_obs);
  
  
  //for the observations that are within range of points, 
  //find the list of states of stations in neighbourhood of points
  //compute share of points for each of those states and add to correponding lambda
  //list of stations near point
  //uvec st_point_list = conv_to<uvec>::from(split(points_local_stations[i],'_'));
  st_point_list_org = (split(points_local_stations[i],'_'));
  st_point_list = st_point_list_org;
  for(uint j=0; j<st_point_list.size();++j) {
    st_point_list[j] -=1; //subtracting 1 to create 0 based indexes in c++
  
  }
  list_obs = which_r(conv_to< Vi >::from(station_id_index_r),st_point_list);
  //print_vec(list_obs);
  //st_point_list = st_point_list - vector<int>(st_point_list.size(),fill::ones); //subtracting 1 to create 0 based indexes in c++
  //alternatively could select list_obs from st_point_list instead of calculating distances above
  
  mat_st_state = umat(list_obs.size(),st_point_list.size());
  obs_st_state = imat(list_obs.size(),st_point_list.size());
  st_point_list_uvec = conv_to<uvec>::from(st_point_list);        

  //for each station compute the states of local stations of points 
  for(uint j=0; j<list_obs.size();++j) {
    //for(uint j=0; j<1;++j) {
    vector<int> j_loc_st = (split(wdc_local_stations[list_obs(j)],'_'));
    uvec j_st_state = conv_to<uvec>::from(split(wdc_sto_state_local[list_obs(j)],'_'));
    uvec Idx = conv_to<uvec>::from(which_r(j_loc_st, st_point_list_org));
    mat_st_state.row(j) = conv_to<urowvec>::from(j_st_state.rows(Idx));            
    irowvec obs_st_state_row(st_point_list.size());
    obs_st_state_row.fill(-1); // -1 serves a equivalent of NA in R
    uint st_id_temp = station_id_index_r(list_obs(j));
    
    vector<int> st_vec(1);
    st_vec[0] = st_id_temp;
    assert(which_r(st_point_list,st_vec).size()==1);
    obs_st_state_row(which_r(st_point_list,st_vec)[0]) = list_obs(j);
    obs_st_state.row(j) = obs_st_state_row;
    //print_vec(conv_to<vec>::from(obs_st_state.row(j)));
    //join_cols(mat_st_state1,j_st_state.rows(Idx));
    //mat_st_state <- rbind(mat_st_state,j_st_state[which(j_loc_st %in% st_point_list)])
  }
  mat_st_state_str = covert_row_str(mat_st_state);
  mat_st_state_str_unq = unique_str(mat_st_state_str);        

}

  
void construct_a_prob_delta_vec(vector< vector<int> > &a_vec, vector< vector<double> > &prob_a,
  vector< vector<double> > &delta_a, imat &obs_st_state, umat &mat_st_state, vector<string> &mat_st_state_str, 
  vector<string> &str_vec, vector<int> &st_point_list, colvec &xdeltain,
  urowvec &col_na, rowvec &delta_avg, urowvec &mat_st_state_row, uvec &station_id_index_r, 
  colvec &wdcobswt) {

  imat a = obs_st_state.rows(which_r_str(mat_st_state_str,str_vec));
  mat_st_state_row = mat_st_state.row(which_r_str(mat_st_state_str,str_vec)[0]);
  //cout << a << endl;
  //now for each column of a remove the -1's and get the expanded grid
  //convert to vector of vectors from columns of a, remove -1 
  //and then apply the recursive strategy to create all combinations
  
  
  for(uint k=0; k <a.n_cols; ++k) {
    //cout << a.col(k)(find(a.col(k)>=0)) << endl;
    ivec a_sub = a.col(k);
    vector<int> a_vec_elem = conv_to< vector<int> >::from(a_sub.elem(find(a_sub>=0)));
    //if a_vec_elem is empty add the first element of that station_id to this 
    //also make correponding entry in col_na =1, which will restrict us from updating those rows, but only
    //use those delta;s for computation of other station probabilities

    if(a_vec_elem.size()==0) {
      uint st_k = which_r(conv_to< vector<int> >::from(station_id_index_r),
      vector<int>(1,st_point_list[k]))[0];
      a_vec_elem.push_back(st_k);
      col_na(k) = 1;
    }
    a_vec.push_back(a_vec_elem);            
  }

//          Vvi obs_no_table;
//          Vi outputTemp;
//          cart_product(obs_no_table, outputTemp, a_vec.begin(), a_vec.end());
//          //to get delta_table, will have to take all elements of obs_no_table and take corresponding elements from 
//          //deltain
//          //first convert obs_no_table to mat
//          umat obs_no_table_mat(obs_no_table.size(),obs_no_table[0].size());
//          for(uint k=0; k<obs_no_table.size(); ++k) {
//            obs_no_table_mat.row(k) = conv_to<urowvec>::from(obs_no_table[k]);
//          }
  
  //cout << obs_no_table_mat << endl;
//          mat delta_table(obs_no_table_mat.n_rows,obs_no_table_mat.n_cols);          
//          for(uint k=0;k<delta_table.n_cols;++k) {
//            uvec obs_no_table_col = obs_no_table_mat.col(k);            
//            delta_table.col(k) = xdeltain.elem(obs_no_table_col);
//          }
  //create prob_a to store probabilities of observatiosn corresponding to 
  //obs_no in a_vec and           
  for(uint k=0; k <a_vec.size(); ++k) {
    uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
    //get obs_wt corresponding to this vector
    colvec obs_wt_col = wdcobswt.elem(a_vec_col);
    colvec prob_col = obs_wt_col/sum(obs_wt_col);
    prob_a.push_back(conv_to< vector<double> >::from(prob_col));
    colvec delta_col = xdeltain.elem(a_vec_col);
    delta_a.push_back(conv_to< vector<double> >::from(delta_col));
  }

  
  for(uint k=0; k <a_vec.size(); ++k) {
    uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
    colvec delta_col = xdeltain.elem(a_vec_col);
    colvec prob_a_col = conv_to<colvec>::from(prob_a[k]);
    delta_avg(k) = sum(delta_col % prob_a_col);  //this is wrong as it woudl a do a term by term prod without summing
  }

}



// [[Rcpp::export]]
SEXP eval_grad_lambda_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in) {
  
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint points_density_col = 2;
    uint min_density_points_col = 3;
    uint max_density_points_col = xpoints.n_cols-1;
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
    
    arma::mat grad_nondensity_cols(xwdcMergedday.nrow(),2,fill::zeros);
    arma::mat grad_density_cols(xwdcMergedday.nrow(), xpoints.n_cols-3, fill::zeros);
    arma::mat grad_delta(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    //arma::mat grad_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {

        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.

          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
              {
                mat util_grad = compute_prob_theta_2(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec, xpoints(i,points_density_col)); 

                rowvec grad_temp = util_grad.row(k);

                uvec rowno(1); rowno(0) = a_vec[k][l];
                grad_nondensity_cols.rows(rowno)  += util_grad.row(k);
              }  
              {
                rowvec grad_density_cols_int = compute_prob_unweighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec);                             
                grad_density_cols.row(a_vec[k][l]) +=  grad_density_cols_int(k)*xpoints( i, span(min_density_points_col,max_density_points_col));                
              }
              {
                vector<mat> ret = compute_prob_2(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec, xpoints(i,points_density_col)); 
                              
                mat util_grad = ret[1];              

                rowvec grad_temp = util_grad.row(k);
                //CLEAN THIS. YOU HAVE ALTERNATE WAY IN HESSIAN IMPLEMENTATION  
                vector< vector<double> > prob_a_temp = prob_a;
                vector<double> temp_vec1(1); temp_vec1[0]=1;
                prob_a_temp[k] = temp_vec1;
                vector< vector<int> > a_vec_temp = a_vec;
                vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
                a_vec_temp[k] = temp_vec2;
                vector< vector<double> > grad_temp_list = prob_a_temp;              
                vector<int> a_vec_unlisted;
                vector<double> grad_temp_unlisted;
                for(uint m=0; m <prob_a_temp.size(); ++m) {                
                  std::transform(grad_temp_list[m].begin(), grad_temp_list[m].end(), 
                    grad_temp_list[m].begin(), std::bind1st(std::multiplies<double>(),grad_temp[m]));                
                  a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                  grad_temp_unlisted.insert(grad_temp_unlisted.end(),grad_temp_list[m].begin(),grad_temp_list[m].end());
                }
             
                rowvec grad_temp_unlisted_rowvec = conv_to<rowvec>::from(grad_temp_unlisted);
                uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
                uvec rowno(1); rowno(0) = a_vec_temp[k][0];
                grad_delta(rowno,a_vec_unlisted_uvec)  += mat(grad_temp_unlisted_rowvec);
              }

            }
            }
          }

        }
    }//end of points loop  

    // colvec wdcobswt_colvec = xwdcmat.col(wdcobswt_col);    
    // grad_nondensity_cols = grad_nondensity_cols% (repmat(wdcobswt_colvec,1,grad_nondensity_cols.n_cols));    
    // grad_density_cols = grad_density_cols % (repmat(wdcobswt_colvec,1,grad_density_cols.n_cols));    
    // grad_delta = grad_delta % (repmat(wdcobswt_colvec,1,grad_delta.n_cols));  
    //grad_t = grad_t%wdcobswt_mat;    
    
    //mat obj_ret = join_rows(grad_nondensity_cols,grad_density_cols,grad_delta); 
    // mat obj_ret_1 = join_rows(grad_nondensity_cols,grad_density_cols);
    // mat obj_ret = join_rows(obj_ret_1,grad_delta);  
    mat obj_ret = join_rows(grad_nondensity_cols,join_rows(grad_density_cols,grad_delta)); 
    return(wrap(obj_ret));  
      
}



// [[Rcpp::export]]
SEXP eval_hessian_lambda_cpp(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in,
         SEXP lambda_multiplers_in) {
 
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1(theta1);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints_temp(points);
    arma::mat xpoints(xpoints_temp.begin(), xpoints_temp.nrow(), xpoints_temp.ncol(), 
                      true);  

    uint min_points_col = 2;
    uint max_points_col = xpoints.n_cols-1;
    assert(max_points_col-min_points_col==xtheta1.size()-2); //to assert other points covariates supplied correspond to density vector.
    
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    NumericVector xv0_vec(v0_vec);
    NumericVector lambda_multiplers(lambda_multiplers_in);
    assert(lambda_multiplers.size()==xwdcMergedday.nrow());
    double beta1 = xtheta1(0);
    double sigma0 = xtheta1(1);
    
    uint wdclat1_col = 3;
    uint wdclon1_col = 4;
    uint wdcobswt_col = 5;
    uint wdcstation_id_index_col = 2;
    
    NumericVector lat1_r = xwdcMergedday(_,wdclat1_col);
    NumericVector lon1_r = xwdcMergedday(_,wdclon1_col);
    NumericVector wdcobswt_r = xwdcMergedday(_,wdcobswt_col);      
    colvec lat1(lat1_r.begin(),lat1_r.size(),true);
    colvec lon1(lon1_r.begin(),lon1_r.size(),true);
    colvec wdcobswt(wdcobswt_r.begin(),wdcobswt_r.size(),true);
        
    mat hessian_theta1_delta(xtheta1.size(),xwdcMergedday.nrow(),fill::zeros);
    mat hessian_delta_sq_t(xwdcMergedday.nrow(),xwdcMergedday.nrow(),fill::zeros);
    mat hessian_theta1_sq(xtheta1.size(),xtheta1.size(),fill::zeros);
    
    uint pointslat1_col = 0;
    uint pointslon1_col = 1;
    
    arma::mat xwdcmat(xwdcMergedday.begin(), xwdcMergedday.nrow(), xwdcMergedday.ncol(), 
                      true);  
    
    uvec station_id_index_r = conv_to< uvec >::from(xwdcmat.col(wdcstation_id_index_col));
    station_id_index_r -=1; //subtracting 1 to create 0 based indexes in c++
    mat station_data_all = xwdcmat.rows(unique_idx(station_id_index_r));
    
    //station_data_all <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
    
    for(uint i=0;i<xpoints.n_rows;i++) {
        vector<int> st_point_list_org;
        vector<int> st_point_list;
        uvec list_obs;
        umat mat_st_state;
        imat obs_st_state;
        uvec st_point_list_uvec;
        vector<string> mat_st_state_str;
        vector<string> mat_st_state_str_unq;  
        construct_mat_st_state_str_unq(i, lat1, lon1,  xpoints(i,pointslat1_col), xpoints(i,pointslon1_col),
        wdc_sto_state_local, wdc_local_stations,
        points_local_stations, st_point_list_org, st_point_list, list_obs, mat_st_state, 
        obs_st_state,  st_point_list_uvec,  mat_st_state_str, mat_st_state_str_unq, station_id_index_r);
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          vector< vector<int> > a_vec;
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          urowvec col_na(st_point_list.size(),fill::zeros);
          rowvec delta_avg(st_point_list.size(),fill::zeros);
          urowvec mat_st_state_row;

          construct_a_prob_delta_vec(a_vec, prob_a, delta_a, obs_st_state, mat_st_state, mat_st_state_str, str_vec, st_point_list,
            xdeltain, col_na, delta_avg, mat_st_state_row, station_id_index_r, wdcobswt);
          //and prob_table by expanding similar to obs_no_table          
//          Vvd prob_table;
//          Vd outputTemp2;
//          cart_product(prob_table, outputTemp2, prob_a.begin(), prob_a.end());          
          //now for each row of delta_table compute the probabilities for each station and gradients,
          //multiply them with weights from
          //prob_table and add to corresponding row.
          mat station_data = station_data_all.rows(st_point_list_uvec);
                              
          for(uint k=0; k <a_vec.size(); ++k) {
            if(col_na(k)==0) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            for(uint l=0; l <a_vec_col.size(); ++l) {
              rowvec deltain_row = delta_avg;

              deltain_row(k) = delta_a[k][l];
              if(lambda_multiplers(a_vec[k][l])==0) continue;
              {                
                mat hessian_theta1_delta_kl = compute_hessian_theta1_delta_weighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec, k, xtheta1.size(), xpoints(i, min_points_col), xpoints(i,span(min_points_col+1,max_points_col)));
                //multiply with observation wt & lambda_multiplers_in(a_vec[k][l])
                hessian_theta1_delta_kl *=  lambda_multiplers(a_vec[k][l]);
                
                //need to expand hessian_delta_sq_kl to reflect gradients wrt
                //a_vec columns to reflect the deltaaveraged gradients.
                //test in a seperate Rcpp file how to repeat rows and columns and then 
                //multiply rows and columns with prob_a values.
                //creating version of a_vec and prob_a which have the focal station 
                //with only one entry and rest of the stations with actual list.
                vector< vector<double> > prob_a_temp = prob_a;
                vector<double> temp_vec1(1); temp_vec1[0]=1;
                prob_a_temp[k] = temp_vec1;
                vector< vector<int> > a_vec_temp = a_vec;
                vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
                a_vec_temp[k] = temp_vec2;
                //cout << "simplify above lines, there should be way of direclty assigning\
                //instead of creating temp vecs" << endl; 
                //unlisting above lists
                vector<int> a_vec_unlisted;
                vector<double> prob_a_unlisted;
                //create list of hessian_delta_sq_kl indexes to expand
                uvec hessian_expand_index;

                for(uint m=0; m <prob_a_temp.size(); ++m) {  
                  uvec hessian_expand_index_temp(prob_a_temp[m].size());
                  hessian_expand_index_temp.fill(m);
                  hessian_expand_index.insert_rows( hessian_expand_index.size(), hessian_expand_index_temp ); 
                  a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                  prob_a_unlisted.insert(prob_a_unlisted.end(),prob_a_temp[m].begin(),prob_a_temp[m].end());                
                }              
                mat weights_mat(prob_a_unlisted.size(),prob_a_unlisted.size(),fill::zeros);
                weights_mat.diag()  = conv_to<vec>::from(prob_a_unlisted);
                mat hessian_theta1_delta_kl_expanded = hessian_theta1_delta_kl.cols(hessian_expand_index);
                //hessian_delta_sq_kl_expanded = hessian_delta_sq_kl_expanded.cols(hessian_expand_index);
                hessian_theta1_delta_kl_expanded = hessian_theta1_delta_kl_expanded * weights_mat;
                uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
                hessian_theta1_delta.cols(a_vec_unlisted_uvec) += hessian_theta1_delta_kl_expanded;
              }
              {
                mat hessian_delta_sq_kl = compute_hessian_delta_sq(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec, k);
                //multiply with point weight and lambda_multiplers_in(a_vec[k][l])
                
                hessian_delta_sq_kl *=  lambda_multiplers(a_vec[k][l]) * xpoints( i, min_points_col); 

                //need to expand hessian_delta_sq_kl to reflect gradients wrt
                //a_vec columns to reflect the deltaaveraged gradients.
                //test in a seperate Rcpp file how to repeat rows and columns and then 
                //multiply rows and columns with prob_a values.
                //creating version of a_vec and prob_a which have the focal station 
                //with only one entry and rest of the stations with actual list.
                vector< vector<double> > prob_a_temp = prob_a;
                vector<double> temp_vec1(1); temp_vec1[0]=1;
                prob_a_temp[k] = temp_vec1;
                vector< vector<int> > a_vec_temp = a_vec;
                vector<int> temp_vec2(1); temp_vec2[0]=a_vec[k][l];
                a_vec_temp[k] = temp_vec2;
                //cout << "simplify above lines, there should be way of direclty assigning\
                //instead of creating temp vecs" << endl; 
                //unlisting above lists
                vector<int> a_vec_unlisted;
                vector<double> prob_a_unlisted;
                //create list of hessian_delta_sq_kl indexes to expand
                uvec hessian_expand_index;

                for(uint m=0; m <prob_a_temp.size(); ++m) {  
                  uvec hessian_expand_index_temp(prob_a_temp[m].size());
                  hessian_expand_index_temp.fill(m);
                  hessian_expand_index.insert_rows( hessian_expand_index.size(), hessian_expand_index_temp ); 
                  a_vec_unlisted.insert(a_vec_unlisted.end(),a_vec_temp[m].begin(),a_vec_temp[m].end());
                  prob_a_unlisted.insert(prob_a_unlisted.end(),prob_a_temp[m].begin(),prob_a_temp[m].end());                
                }              
                mat weights_mat(prob_a_unlisted.size(),prob_a_unlisted.size(),fill::zeros);
                weights_mat.diag()  = conv_to<vec>::from(prob_a_unlisted);
                mat hessian_delta_sq_kl_expanded = hessian_delta_sq_kl.rows(hessian_expand_index);
                hessian_delta_sq_kl_expanded = hessian_delta_sq_kl_expanded.cols(hessian_expand_index);
                hessian_delta_sq_kl_expanded = weights_mat * hessian_delta_sq_kl_expanded * weights_mat;
                uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
                hessian_delta_sq_t(a_vec_unlisted_uvec,a_vec_unlisted_uvec) += hessian_delta_sq_kl_expanded;
              }
              {
                mat hessian_theta1_sq_kl = compute_hessian_theta1_sq_weighted(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                  xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
                xv0_vec, k, xtheta1.size(), xpoints(i, min_points_col), xpoints(i,span(min_points_col+1,max_points_col)));
                //multiply with observation wt & lambda_multiplers_in(a_vec[k][l])
                
                hessian_theta1_sq_kl *=  lambda_multiplers(a_vec[k][l]);
                

                hessian_theta1_sq += hessian_theta1_sq_kl;
              }
            }
            }
          }

        }
    }//end of points loop  

    mat hessian_all = join_cols(join_rows(hessian_theta1_sq,hessian_theta1_delta),
                                join_rows(hessian_theta1_delta.t(),hessian_delta_sq_t));

    return(wrap(hessian_all));  
      
}

