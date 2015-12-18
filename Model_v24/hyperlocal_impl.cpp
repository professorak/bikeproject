// [[Rcpp::export]]
SEXP eval_lambda_delta_list_cpp_new(SEXP deltain , SEXP  theta1 ,SEXP wdcMergedday , SEXP points,
         SEXP no_st, SEXP max_walking_dis,
         SEXP v0_vec, SEXP v0_vec_weights, SEXP wdc_sto_state_local_in, 
         SEXP wdc_local_stations_in, SEXP points_local_stations_in, SEXP nonden_ceoflength_in) {
  // cout << "hi" << endl;
  //try {
    std::vector<string> wdc_sto_state_local = Rcpp::as< std::vector<string> >(wdc_sto_state_local_in); 
    std::vector<string> wdc_local_stations = Rcpp::as< std::vector<string> >(wdc_local_stations_in); 
    std::vector<string> points_local_stations = Rcpp::as< std::vector<string> >(points_local_stations_in);    
    NumericVector deltain_r(deltain);
    colvec xdeltain(deltain_r.begin(),deltain_r.size(),true);
    NumericVector xtheta1_in(theta1);  rowvec xtheta1(xtheta1_in.begin(),xtheta1_in.size(),true);
    NumericMatrix xwdcMergedday(wdcMergedday);
    NumericMatrix xpoints(points);
    uint points_density_col = 2;
    int xno_st = as<int>(no_st); 
    double xmax_walking_dis = as<double>(max_walking_dis);
    double nonden_ceoflength = as<double>(nonden_ceoflength_in);

    NumericVector xv0_vec(v0_vec);
    NumericVector xv0_vec_weights(v0_vec_weights);
    
    rowvec beta1_vec = xtheta1(span(0,nonden_ceoflength-1));
    beta1_vec.shed_col(1);
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
    
    //construct obs_no_vec_all, prob_vec_all and delta_avg_all
    vector< vector<uint> > obs_no_vec_all;
    vector< vector<double> > prob_vec_all;
    vector<double> delta_avg_all;
    
    construct_obs_no_prob_delta_avg_all_vec(obs_no_vec_all, prob_vec_all, delta_avg_all, 
        xdeltain, wdcobswt, xno_st, station_id_index_r);
    
    for(uint i=0;i<xpoints.nrow();i++) {
    //  cout << "point no" << i << endl;
    //for(uint i=2;i<3;i++) {
        
        vector<int> st_point_list_org = (split(points_local_stations[i],'_'));
        vector<int> st_point_list = st_point_list_org;
        for(uint j=0; j<st_point_list.size();++j) {
          st_point_list[j] -=1; //subtracting 1 to create 0 based indexes in c++        
        } //how to write above efficiently
        uvec st_point_list_uvec = conv_to<uvec>::from(st_point_list);

        mat station_data = station_data_all.rows(st_point_list_uvec);

        //subset obs_no_vec_all, prob_vec_all and delta_avg_all
        vector< vector<uint> > obs_no_vec;
        vector< vector<double> > prob_vec;
        vector<double> delta_avg;

        for(uint j=0; j<st_point_list.size();++j) {
          obs_no_vec.push_back(obs_no_vec_all[st_point_list[j]]); 
          prob_vec.push_back(prob_vec_all[st_point_list[j]]); 
          delta_avg.push_back(delta_avg_all[st_point_list[j]]); 
        }
        
        //now loop over obs_no_vec and for each observation construct stockout state vector
        //and compute lamba, gradients, hessian etc.
        for(uint k=0; k< obs_no_vec.size(); ++k) {
          vector<uint> obs_no_vec_col = obs_no_vec[k];
          for(uint l=0; l <obs_no_vec_col.size(); ++l) {
            rowvec deltain_row = delta_avg;
            deltain_row(k) = xdeltain(obs_no_vec_col[l]);
            vector<int> station_local_stations = (split(wdc_local_stations[obs_no_vec_col[l]],'_'));
            uvec station_local_state = conv_to<uvec>::from(split(wdc_sto_state_local[obs_no_vec_col[l]],'_'));
            uvec station_point_intersection_index = conv_to<uvec>::from(which_r(station_local_stations, st_point_list_org));
            urowvec station_point_stkt_state = conv_to<urowvec>::from(station_local_state.rows(station_point_intersection_index));

            vector<mat> ret = compute_prob(i, station_data, xpoints, wdclat1_col, wdclon1_col, pointslat1_col, 
              pointslon1_col, beta1_vec, sigma0, xdeltain, st_point_list_uvec, deltain_row, station_point_stkt_state,  
              xv0_vec, xv0_vec_weights, points_density_col); 

            rowvec lambda_st_t = ret[0];              
            lambda_t(obs_no_vec_col[l]) +=  lambda_st_t(k);

            mat util_grad = ret[1];              

            rowvec grad_temp = util_grad.row(k);

            vector< vector<double> > prob_vec_temp = prob_vec;
            vector<double> temp_vec1(1); temp_vec1[0]=1;
            prob_vec_temp[k] = temp_vec1;
            vector< vector<uint> > obs_no_vec_temp = obs_no_vec;
            vector<uint> temp_vec2(1); temp_vec2[0]=obs_no_vec[k][l];
            obs_no_vec_temp[k] = temp_vec2;
            vector< vector<double> > grad_temp_list = prob_vec_temp;              
            vector<uint> obs_no_vec_unlisted;
            vector<double> grad_temp_unlisted;
            for(uint m=0; m <prob_vec_temp.size(); ++m) {                
              std::transform(grad_temp_list[m].begin(), grad_temp_list[m].end(), 
                grad_temp_list[m].begin(), std::bind1st(std::multiplies<double>(),grad_temp[m]));

              obs_no_vec_unlisted.insert(obs_no_vec_unlisted.end(),obs_no_vec_temp[m].begin(),obs_no_vec_temp[m].end());
              grad_temp_unlisted.insert(grad_temp_unlisted.end(),grad_temp_list[m].begin(),grad_temp_list[m].end());
            }

            rowvec grad_temp_unlisted_rowvec = conv_to<rowvec>::from(grad_temp_unlisted);
            uvec obs_no_vec_unlisted_uvec = conv_to<uvec>::from(obs_no_vec_unlisted);
            uvec rowno(1); rowno(0) = obs_no_vec_temp[k][0];
            
            grad_t(rowno,obs_no_vec_unlisted_uvec)  += mat(grad_temp_unlisted_rowvec);

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
