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



// [[Rcpp::export]]
void test_main() {
  cout << "hi" << endl;
  mat A = randu<mat>(2,2);
  A.print("A:");
  cout << endl; 
  urowvec row_indices(3);
  row_indices(0)=0;
  row_indices(1)=0;
  row_indices(2)=1;
  mat A1 = A.rows(row_indices);
  A1 = A1.cols(row_indices);
  A1.print("A1:");
  cout << endl; 
  
  vec weights(2);
  weights(0)=1;
  weights(1)=2;
  mat weights_mat(2,2,fill::zeros);
  
  weights_mat.diag()  = weights;

  //multiply cols with weights
  mat B = A * weights_mat;
  B.print("B:");
  cout << endl; 
  
  //multiply rows with weights
  mat C = weights_mat * A;
  C.print("C:");
  cout << endl; 
  
  uvec hessian_expand_index;
  for(uint m=0; m <2; ++m) {  
    uvec hessian_expand_index_temp(m+1);
    hessian_expand_index_temp.fill(m+1);
    hessian_expand_index_temp.print("hessian_expand_index_temp: ");
    hessian_expand_index.insert_rows( hessian_expand_index.size(), hessian_expand_index_temp ); 
    //hessian_expand_index = join_cols(hessian_expand_index,hessian_expand_index_temp);
  }
  hessian_expand_index.print("hessian_expand_index: ");           
}



