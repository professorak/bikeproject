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
        
    arma::mat lambda_t(xwdcMergedday.nrow(), xpoints.n_cols-2, fill::zeros);
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
    //  cout << "point no" << i << endl;
    //for(uint i=2;i<3;i++) {
        colvec dis_vIdx = latlondistance(lat1, lon1, 
                                         xpoints(i,pointslat1_col), xpoints(i,pointslon1_col));
        // ref for find - http://arma.sourceforge.net/docs.html#find
//        uvec list_obs = find(dis_vIdx <=xmax_walking_dis);
//        if(list_obs.size()==0) continue;
//        print_vec(list_obs);
        
        
        //for the observations that are within range of points, 
        //find the list of states of stations in neighbourhood of points
        //compute share of points for each of those states and add to correponding lambda
        //list of stations near point
        //uvec st_point_list = conv_to<uvec>::from(split(points_local_stations[i],'_'));
        vector<int> st_point_list_org = (split(points_local_stations[i],'_'));
        vector<int> st_point_list = st_point_list_org;
        for(uint j=0; j<st_point_list.size();++j) {
          st_point_list[j] -=1; //subtracting 1 to create 0 based indexes in c++
        }
        
        uvec list_obs = which_r(conv_to< Vi >::from(station_id_index_r),st_point_list);
        //print_vec(list_obs);
        //st_point_list = st_point_list - vector<int>(st_point_list.size(),fill::ones); //subtracting 1 to create 0 based indexes in c++
        //alternatively could select list_obs from st_point_list instead of calculating distances above
        
        umat mat_st_state(list_obs.size(),st_point_list.size());
        imat obs_st_state(list_obs.size(),st_point_list.size());
        uvec st_point_list_uvec = conv_to<uvec>::from(st_point_list);        

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
        vector<string> mat_st_state_str = covert_row_str(mat_st_state);
        vector<string> mat_st_state_str_unq = unique_str(mat_st_state_str);        
        
        for(uint j=0; j< mat_st_state_str_unq.size(); ++j) {
        //for(uint j=0; j< 1; ++j) {
          vector<string> str_vec(1);
          str_vec[0] = mat_st_state_str_unq[j];          
          //cout << which_r_str(mat_st_state_str,str_vec) << endl;
          imat a = obs_st_state.rows(which_r_str(mat_st_state_str,str_vec));
          urowvec mat_st_state_row = mat_st_state.row(which_r_str(mat_st_state_str,str_vec)[0]);
          //cout << a << endl;
          //now for each column of a remove the -1's and get the expanded grid
          //convert to vector of vectors from columns of a, remove -1 
          //and then apply the recursive strategy to create all combinations
          vector< vector<int> > a_vec;
          urowvec col_na(st_point_list.size(),fill::zeros);
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


          //create prob_a to store probabilities of observatiosn corresponding to 
          //obs_no in a_vec and           
          vector< vector<double> > prob_a;
          vector< vector<double> > delta_a;
          for(uint k=0; k <a_vec.size(); ++k) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            //get obs_wt corresponding to this vector
            colvec obs_wt_col = wdcobswt.elem(a_vec_col);
            colvec prob_col = obs_wt_col/sum(obs_wt_col);
            prob_a.push_back(conv_to< vector<double> >::from(prob_col));
            colvec delta_col = xdeltain.elem(a_vec_col);
            delta_a.push_back(conv_to< vector<double> >::from(delta_col));
          }

          rowvec delta_avg(st_point_list.size(),fill::zeros);
          for(uint k=0; k <a_vec.size(); ++k) {
            uvec a_vec_col = conv_to<uvec>::from(a_vec[k]);
            colvec delta_col = xdeltain.elem(a_vec_col);
            colvec prob_a_col = conv_to<colvec>::from(prob_a[k]);
            delta_avg(k) = sum(delta_col % prob_a_col);  //this is wrong as it woudl a do a term by term prod without summing
          }

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
              
              mat hessian_delta_sq_kl = compute_hessian_delta_sq(i, station_data, wdclat1_col, wdclon1_col, xpoints(i,pointslat1_col), 
                xpoints(i,pointslon1_col), beta1, sigma0, xdeltain, st_point_list_uvec, deltain_row, mat_st_state_row,  
              xv0_vec, k);
              //multiply with point weight and lambda_multiplers_in(a_vec[k][l])
              hessian_delta_sq_kl *=  wdcobswt(a_vec[k][l])*lambda_multiplers_in(a_vec[k][l]) * xpoints( i, span(min_points_col)); 
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
              cout << "simplify above lines, there should be way of direclty assigning
              instead of creating temp vecs" << endl; 
              //unlisting above lists
              vector<int> a_vec_unlisted;
              vector<int> prob_a_unlisted;
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
              weights_mat.diag()  = prob_a_unlisted;
              mat hessian_delta_sq_kl_expanded = hessian_delta_sq_kl.rows(hessian_expand_index);
              hessian_delta_sq_kl_expanded = hessian_delta_sq_kl_expanded.cols(hessian_expand_index);
              hessian_delta_sq_kl_expanded = weights_mat * hessian_delta_sq_kl_expanded * weights_mat;

              uvec a_vec_unlisted_uvec = conv_to<uvec>::from(a_vec_unlisted);
              hessian_delta_sq_t(a_vec_unlisted_uvec,a_vec_unlisted_uvec) = hessian_delta_sq_kl_expanded;
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
          urowvec no_focal_indexes(no_t_st,fill::zeros);
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
              util_prob_t_nofocal = util_prob_t(no_focal_indexes);
              mat B_1 = repmat(vectorise( util_prob_t_nofocal, 0),1,no_t_st-1);
              mat B_2 = repmat(util_prob_t_nofocal,no_t_st-1,1);
              mat I = eye<mat>(no_t_st-1,no_t_st-1);
              mat A_2 = -B_1 % (I - 2*B_2);
              A(no_focal_indexes,no_focal_indexes) = A_2;

              rowvec A_1 = (1-2*util_prob_t(focal_station_index)) * (e_f-util_prob_t);
              A.row(focal_station_index) = A_1;
              A.col(focal_station_index) = A_1;
              hessian_delta_sq_t += util_prob_t(focal_station_index)*A;
          }
          hessian_delta_sq_t *= (1/xv0_vec.size());
          return((hessian_delta_sq_t));  
}











