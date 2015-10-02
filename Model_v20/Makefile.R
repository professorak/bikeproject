PKG_CXXFLAGS=`Rscript -e "RcppArmadillo:::CxxFlags()"`
PKG_LIBS=`Rscript -e "RcppArmadillo:::LdFlags()"`
#PKG_CXXFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`
#PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`
R CMD SHLIB eval_func_2_new_deltaaveraged.cpp


