#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double marginal_column_bin(NumericVector y, arma::field<arma::vec> log_Vnt, int n_rep, 
                            double a0, double b0, double gamma_mfm){
  // should let 
  // for univariate
  int n=y.length();
  int i, i_rep, i_clst, n_clst, idx;
  // for integrals and cluster-wise integrals:
  double log_int, log_int_clst, log_int_clst_max;
  // for Monte Carlo replications
  double log_lik_rep, log_lik_rep_max;
  double weight_temp, weight_sum;
  arma::vec log_Vnt_temp;
  IntegerVector ind;
  NumericVector an, bn, n_size;
  double an_temp, bn_temp;
  NumericVector log_int_list, range;
  log_lik_rep=R_NegInf;
  for(i_rep=0; i_rep < n_rep; i_rep++){
    // initialize
    an=rep(a0,n);
    bn=rep(b0,n);
    
    an_temp=a0+y[0];
    bn_temp=b0+1;
    
    // for the first observation, it forms a cluster by itself
    n_size=rep(0,n);
    n_size[0]=1;
    n_clst=1;
    
    log_int=lgamma(a0 + y[0]) + lgamma(1 - y[0] + b0) - lgamma(1 + a0 + b0) -
      (lgamma(a0) + lgamma(b0) - lgamma(a0 + b0));
    an[0]=an_temp;
    bn[0]=bn_temp;
    
    for(i=1; i < n; i++){
      log_Vnt_temp=log_Vnt(i);
      log_int_list=rep(0,n_clst+1);
      weight_sum=0;
      // initially 0, NegInf after log
      log_int_clst=R_NegInf;
      log_int_clst_max=R_NegInf;
      // for each cluster, calculate the integral
      for(i_clst=0; i_clst < n_clst; i_clst++){
        an_temp=an[i_clst] + y[i];
        bn_temp=bn[i_clst] + 1 - y[i];
        
        log_int_list[i_clst] = lgamma(an_temp) + lgamma(bn_temp) - lgamma(an_temp + bn_temp) -
          (lgamma(an[i_clst]) + lgamma(bn[i_clst]) - lgamma(an[i_clst] + bn[i_clst]));
        // get the weight for the cluster i_clst
        weight_temp=n_size[i_clst]+gamma_mfm;
        // calculate the summation for the sum of the log integrals
        weight_sum=weight_sum+weight_temp;
        // get the maximal integral
        log_int_clst_max=max(NumericVector::create(log_int_list[i_clst],log_int_clst_max));
        // log of the sum of the weighted integrals
        log_int_clst=log(exp(log_int_clst-log_int_clst_max)+
          weight_temp*exp(log_int_list[i_clst]-log_int_clst_max))+log_int_clst_max;
        // reassign the log_int_list to incorporate weight:
        log_int_list[i_clst]=log_int_list[i_clst]+log(weight_temp);
      }
      an_temp = a0 + y[i];
      bn_temp = b0 + 1 - y[i];
      
      log_int_list[n_clst] = lgamma(an_temp) + lgamma(bn_temp) - lgamma(an_temp + bn_temp) -
        (lgamma(a0) + lgamma(b0) - lgamma(a0 + b0));
      
      weight_temp=exp(log_Vnt_temp[n_clst]-log_Vnt_temp[n_clst-1]+log(gamma_mfm));
      weight_sum=weight_sum+weight_temp;
      log_int_clst_max=max(NumericVector::create(log_int_list[n_clst],log_int_clst_max));
      log_int_clst=log(exp(log_int_clst-log_int_clst_max)+
        weight_temp*exp(log_int_list[n_clst]-log_int_clst_max))+log_int_clst_max;
      // subtract the sum of the weights:
      log_int_clst=log_int_clst-log(weight_sum);
      log_int_list[n_clst]=log_int_list[n_clst]+log(weight_temp);
      // accumulate u_i
      log_int=log_int+log_int_clst;
      // assign the label of the observation i:
      log_int_clst_max=max(log_int_list);
      for (i_clst=0; i_clst < (n_clst+1); i_clst++) {
        log_int_list[i_clst]=exp(log_int_list[i_clst]-log_int_clst_max);
      }
      range=seq(0,n_clst);
      ind=sample(range,1,true,log_int_list);
      idx=ind[0];
      
      if(idx < n_clst){
        n_size[idx]=n_size[idx]+1;
        
        an_temp = an[idx] + y[i];
        bn_temp = bn[idx] + 1 - y[i];
        
        an[idx]=an_temp;
        bn[idx]=bn_temp;
      } else {
        n_size[n_clst]=1;
        
        an_temp = a0 + y[i];
        bn_temp = b0 + 1 - y[i];
        
        an[idx]=an_temp;
        bn[idx]=bn_temp;
        n_clst=n_clst+1;
      }
    }
    log_lik_rep_max=max(NumericVector::create(log_lik_rep,log_int));
    log_lik_rep=log(exp(log_lik_rep-log_lik_rep_max)+exp(log_int-log_lik_rep_max))+log_lik_rep_max;
  }
  log_lik_rep=log_lik_rep-log(n_rep);
  return(log_lik_rep);
}

// [[Rcpp::export]]
void columnAssign(int i, NumericMatrix X,arma::mat Y){
  // Y is a column matrix
  // assign Y to the i-th column of X
  int k,n=X.nrow();
  for(k=0;k<n;k++){
    X(k,i)=Y(k);
  }
}

// [[Rcpp::export]]
void i_col_remove(int i_col, int n_col, int n_row,
                  NumericVector col_assign, NumericVector col_size, NumericVector col_n_clst,
                  NumericMatrix row_assign, NumericMatrix row_size, NumericVector row_n_clst){
  int cur_col_clst, i;
  
  cur_col_clst = col_assign[i_col];// current cluster assignment
  col_assign[i_col] = (-1);
  //print(col_assign);
  --col_size[cur_col_clst - 1]; //decrease the number of grids by 1 in the corresponding cluster.
  
  //If a grid is clustered by itself
  if (col_size[cur_col_clst - 1] == 0) {
    for (i = 0; i < n_col; i++) {
      if (col_assign[i] > cur_col_clst ) {
        --col_assign[i];
      }
    }
    for (i = cur_col_clst; i < col_n_clst[0]; i++) {
      col_size[i-1] = col_size[i];
      row_n_clst[i-1] = row_n_clst[i];
      row_size(_,i-1) = row_size(_,i);
      row_assign(_,i-1) = row_assign(_,i);
    }
    col_size[col_n_clst[0] - 1] = 0;
    row_n_clst[col_n_clst[0] - 1] = 0;
    row_size(_,col_n_clst[0]-1) = rep(0,n_row);
    row_assign(_,col_n_clst[0]-1) = rep(0,n_row);
    --col_n_clst[0];
  }
}

// [[Rcpp::export]]
void i_row_remove(int i_row, int n_row, NumericVector col_n_clst,
                  NumericMatrix row_assign, NumericMatrix row_size, NumericVector row_n_clst){
  int cur_row_clst, i, j;
  for(i=0; i < col_n_clst[0]; i++){
    cur_row_clst = row_assign(i_row,i);// current cluster assignment
    row_assign(i_row,i) = (-1);
    --row_size(cur_row_clst-1,i);//decrease the number of grids by 1 in the corresponding cluster.
    //If a grid is clustered by itself
    if (row_size(cur_row_clst-1,i) == 0) {
      for (j = 0; j < n_row; j++) {
        if (row_assign(j,i) > cur_row_clst ) {
          --row_assign(j,i);
        }
      }
      for (j = cur_row_clst; j < row_n_clst[i]; j++) {
        row_size(j-1,i) = row_size(j,i);
      }
      row_size(row_n_clst[i] - 1,i) = 0;
      --row_n_clst[i];
    }
  }
}

// [[Rcpp::export]]
void partition_samp(int n, double gamma_mfm, arma::field<arma::vec> log_Vnt,
                    NumericVector partition, NumericVector size_list, NumericVector n_clst){
  NumericVector prob_list(n+0L,0), tprob;
  arma::vec log_Vnt_temp;
  IntegerVector range, ind;
  int i, k;
  double prob_max;
  
  partition[0]=1;
  size_list[0]=1;
  n_clst[0]=1;
  for(i=1; i<n; i++){
    std::fill(prob_list.begin(), prob_list.end(), 0);
    log_Vnt_temp=log_Vnt(i);
    for(k=0; k<n_clst[0]; k++){
      prob_list[k]=log(size_list[k]+gamma_mfm);
    }
    prob_list[n_clst[0]]=log_Vnt_temp[n_clst[0]]-log_Vnt_temp[n_clst[0]-1]+log(gamma_mfm);
    prob_max=max(prob_list);
    for(k=0; k<(n_clst[0]+1); k++){
      prob_list[k]=exp(prob_list[k]-prob_max);
    }
    range=seq(0,n_clst[0]);
    tprob=prob_list[range];
    ind=sample(range,1,true,tprob);
    partition[i]=ind[0]+1;
    if(ind[0] < n_clst[0]){
      ++size_list[ind[0]];
    } else {
      size_list[ind[0]]=1;
      ++n_clst[0];
    }
  }
}

// [[Rcpp::export]]
double log_ml_bin(arma::mat y_mat, arma::uvec col_assign_temp_idx, arma::uvec col_assign_union_idx,
                  double row_n_clst, arma::vec row_assign_temp,
                  double a0, double b0){
  int k;
  double n_row_temp = 1, n_col_temp = 1, log_lik = 0;
  double a0_temp, b0_temp, y_sum;
  arma::uvec row_assign_temp_idx;
  arma::mat y_sub;
  
  for(k = 0; k < row_n_clst; k++){
    // find the sub-matrix
    row_assign_temp_idx = arma::find(row_assign_temp==(k+1));
    
    // adding back the i_col_k+1 - th column, calculate the marginal likelihood
    y_sub = y_mat.submat(row_assign_temp_idx,col_assign_temp_idx);
    
    // excluding the i_col_k+1 - th column, calculate the marginal likelihood
    n_col_temp = y_sub.n_cols;
    n_row_temp = y_sub.n_rows;
    y_sum = arma::accu(y_sub);
    
    // marginal prior on sub matrix
    a0_temp = a0 + y_sum;
    b0_temp = b0 + n_col_temp * n_row_temp - y_sum;
    
    log_lik = log_lik - // minus here
      (lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp)); 
    
    // adding back the i_col_k+1 - th column, calculate the marginal likelihood
    y_sub = y_mat.submat(row_assign_temp_idx,col_assign_union_idx);
    
    // adding back the i_col_k+1 - th column, calculate the marginal likelihood
    n_col_temp = y_sub.n_cols;
    n_row_temp = y_sub.n_rows;
    y_sum = arma::accu(y_sub);
    
    // marginal prior on sub matrix
    a0_temp = a0 + y_sum;
    b0_temp = b0 + n_col_temp * n_row_temp - y_sum;
    
    log_lik = log_lik + // plus here
      (lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp)); 
  }
  
  return(log_lik);
}

// [[Rcpp::export]]
double log_ml_bin_col(arma::mat y_mat,
                      double row_n_clst, arma::vec row_assign_temp,
                      double a0, double b0){
  int k;
  double n_row_temp = 1, n_col_temp = 1, log_lik = 0;
  double a0_temp, b0_temp, y_sum;
  arma::uvec row_assign_temp_idx;
  arma::mat y_sub;
  
  n_col_temp = y_mat.n_cols;
  for(k = 0; k < row_n_clst; k++){
    // find the sub-matrix
    row_assign_temp_idx = arma::find(row_assign_temp==(k+1));
    
    // adding back the i_col_k+1 - th column, calculate the marginal likelihood
    y_sub = y_mat.rows(row_assign_temp_idx);
    
    // excluding the i_col_k+1 - th column, calculate the marginal likelihood
    n_row_temp = y_sub.n_rows;
    y_sum = arma::accu(y_sub);
    
    // marginal prior on sub matrix
    a0_temp = a0 + y_sum;
    b0_temp = b0 + n_col_temp * n_row_temp - y_sum;
    
    log_lik = lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp) -
      (lgamma(a0) + lgamma(b0) - lgamma(a0 + b0)); 
  }
  
  return(log_lik);
}

// [[Rcpp::export]]
double log_ml_bin_row(arma::mat y_sub,
                      arma::uvec row_assign_temp_idx, 
                      arma::uvec row_assign_union_idx,
                      double a0, double b0){
  
  double n_row_temp = 1, n_col_temp = 1, log_lik = 0;
  double a0_temp, b0_temp, y_sum;
  arma::mat y_sub_rows,y_sub_sum;
  
  // before adding the i_col_k+1 - th column, calculate the marginal likelihood
  y_sub_rows = y_sub.rows(row_assign_temp_idx);
  
  n_row_temp = y_sub_rows.n_rows;
  n_col_temp = y_sub_rows.n_cols;
  y_sum = arma::accu(y_sub_rows);
  
  a0_temp = a0 + y_sum;
  b0_temp = b0 + n_row_temp * n_col_temp - y_sum;
  
  log_lik = log_lik - //minus here
    (lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp));
  
  // adding back the i_col_k+1 - th column, calculate the marginal likelihood
  y_sub_rows = y_sub.rows(row_assign_union_idx);
  
  n_row_temp = y_sub_rows.n_rows;
  n_col_temp = y_sub_rows.n_cols;
  y_sum = arma::accu(y_sub_rows);
  
  a0_temp = a0 + y_sum;
  b0_temp = b0 + n_row_temp * n_col_temp - y_sum;
  
  log_lik = log_lik + //minus here
    (lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp)); 
  
  return(log_lik);
}

// [[Rcpp::export]]
double log_ml_bin_row_i(arma::mat y_i, double a0, double b0){
  
  double n_row_temp = 1, n_col_temp = 1, log_lik = 0;
  double a0_temp, b0_temp, y_sum;
  
  // before adding the i_col_k+1 - th column, calculate the marginal likelihood
  n_row_temp = y_i.n_rows;
  n_col_temp = y_i.n_cols;
  y_sum = arma::accu(y_i);
  
  a0_temp = a0 + y_sum;
  b0_temp = b0 + n_row_temp * n_col_temp - y_sum;
  
  log_lik = log_lik +
    (lgamma(a0_temp) + lgamma(b0_temp) - lgamma(a0_temp + b0_temp)) -
    (lgamma(a0) + lgamma(b0) - lgamma(a0 + b0));
  
  return(log_lik);
}


// [[Rcpp::export]]
void consBinClstK(arma::mat y_mat, NumericVector col_ml, 
                 arma::field<arma::vec> log_Vnt, arma::field<arma::vec> log_Vnt_Cons_K,
                 NumericVector col_assign, NumericMatrix row_assign,
                 NumericVector col_size, NumericMatrix row_size,
                 NumericVector col_n_clst, NumericVector row_n_clst,
                 double a0, double b0,
                 int n_iter, int n_rep, double gamma_mfm){
  // col_assign is the column assignment, of which the length should be n_col
  // row_assign is a n_row * n_col matrix, with the i-th column recording the row assignment of i-th column cluster
  // col_size is a vector,of which the length is n_col, with the first col_n_clst[0] entries being non-zero
  // row_size is a n_row * n_col matrix, with the i-th column recording the i-th column cluster size
  // col_n_clst is a length 1 vector, with the first entry recording the number of column clusters
  // row_n_clst is a length n_col vector, with the first col_n_clst[0] entries recording the number of row clusters within each column cluster
  int n_row=y_mat.n_rows,n_col=y_mat.n_cols;
  int i_iter, i_row, i_col, i_rep, i_col_temp;
  int i_row_k, i_col_k, k;
  int logVntIdx;
  
  IntegerVector range;
  IntegerVector ind;
  // for column reassignment
  arma::vec col_assign_temp, col_assign_union, row_assign_temp, row_assign_union;// inorder to convert NumericVector to arma::vec
  arma::uvec col_assign_temp_idx, col_assign_union_idx, row_assign_temp_idx, row_assign_union_idx;// record the index of a specific column and row cluster
  
  arma::mat y_sub,y_sub_sum;// submatrix of y given column index and row index
  
  NumericVector col_log_lik(n_col+0L,0),col_log_lik_temp(2+0L,0);
  NumericVector log_lik_temp;
  double log_lik_max, log_lik_comp;
  arma::vec log_Vnt_last;
  NumericVector partition_row(n_row+0L,0), size_row(n_row+0L,0), clst_row(1+0L,0);
  // for row reassignment
  arma::mat y_i;
  NumericVector row_log_lik(n_row+0L,0);
  
  arma::mat y_mat_sub, y_mat_subPart;
  NumericMatrix row_assign_single(n_row,1), row_size_single(n_row,1);
  NumericVector row_n_clst_single(1+0L,0), col_n_clst_single(1+0L,0);
  
  for(i_iter=0; i_iter < n_iter; i_iter++){
    for(i_col=0; i_col < n_col; i_col++){
      // get the label of i_col
      i_col_temp = col_assign[i_col];
      // check if the constraint is violated when one of the columns is removed.
      if(row_n_clst[i_col_temp] > col_size[i_col_temp]/2){
        // initialization. Assume that all the observations are labelled the same
        std::fill(row_assign_single.begin(),row_assign_single.end(),1);
        std::fill(row_size_single.begin(),row_size_single.end(),0);
        row_size_single(0,0) = n_row;
        row_n_clst_single[0] = 1;
        col_n_clst_single[0] = 1;
        
        // duplicate the current column assign and remove column i_col
        col_assign_temp = col_assign;
        col_assign_temp[i_col] = -1;
        
        // get the remainder index
        col_assign_temp_idx=arma::find(col_assign_temp==(i_col_temp+1));
        
        y_mat_sub = y_mat.cols(col_assign_temp_idx);
        
        // for each column cluster, sample a new row cluster:
        logVntIdx = (col_assign_temp_idx.n_elem+1)/2;
        
        // update row assign, use log_Vnt(n_row-2):
        log_Vnt_last = log_Vnt_Cons_K(logVntIdx-1);
        
        for(i_rep=0; i_rep < n_rep; i_rep++){
          for(i_row=0; i_row < n_row; i_row++){
            
            i_row_remove(i_row, n_row, col_n_clst_single, row_assign_single, 
                         row_size_single, row_n_clst_single);
            
            // get the i-th row for the marginal likelihood
            y_i = y_mat_sub.row(i_row);
            
            std::fill(row_log_lik.begin(), row_log_lik.end(), 0);
            
            row_assign_temp = row_assign_single(_,0);
            row_assign_union = row_assign_single(_,0);
            
            for(i_row_k=0; i_row_k < row_n_clst_single[0]; i_row_k++){
              
              row_assign_temp_idx = arma::find(row_assign_temp == (i_row_k + 1));
              
              row_assign_union[i_row] = i_row_k + 1;
              row_assign_union_idx = arma::find(row_assign_union == (i_row_k + 1));
              
              row_log_lik[i_row_k] = row_log_lik[i_row_k] + 
                log_ml_bin_row(y_mat_sub,
                               row_assign_temp_idx, 
                               row_assign_union_idx,
                               a0, b0) + 
                                 log(row_size_single(i_row_k,0)+gamma_mfm);
              
            }
            
            if(row_n_clst_single[0] < logVntIdx){
              row_log_lik[row_n_clst_single[0]] = log_ml_bin_row_i(y_i, a0, b0)+
                log_Vnt_last[row_n_clst_single[0]] - log_Vnt_last[row_n_clst_single[0]-1] + 
                log(gamma_mfm);
              range = seq(0,row_n_clst_single[0]);
              
              log_lik_temp = row_log_lik[range];
              log_lik_max = max(log_lik_temp);
              
              for(i_row_k=0; i_row_k < (row_n_clst_single[0]+1); i_row_k++){
                log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
              }
              
            }else{
              range = seq(0,row_n_clst_single[0]-1);
              
              log_lik_temp = row_log_lik[range];
              log_lik_max = max(log_lik_temp);
              
              for(i_row_k=0; i_row_k < (row_n_clst_single[0]); i_row_k++){
                log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
              }
              
            }
            
            ind = sample(range,1,true,log_lik_temp);
            k = ind[0];
            row_assign_single(i_row,0) = k+1;
            
            if(k < row_n_clst_single[0]){
              ++row_size_single(k,0);
            }else{
              row_size_single(row_n_clst_single[0],0) = 1;
              ++row_n_clst_single[0];
            }
          }      
        }
        // given the last row_assign, calculate the difference in the row-wise marginal likelihood
        // now convert numericMatrix into arma::vec
        row_assign_temp = row_assign_single(_,0);
        row_assign_union = row_assign(_,i_col_temp);
          
        log_lik_comp = log_ml_bin_col(y_mat_sub, row_n_clst_single[0], 
                                      row_assign_temp, a0, b0) - 
                      log_ml_bin_col(y_mat_sub, row_n_clst[i_col_temp], 
                                     row_assign_union, a0, b0);
        // also need to consider the difference in the prior mass
        for(k = 0; k < row_n_clst_single[0]; k++){
          log_lik_comp = log_lik_comp + 
            (lgamma(gamma_mfm + row_size_single(k,0)) - lgamma(gamma_mfm));
        }
        log_lik_comp = log_lik_comp + log_Vnt_last[row_n_clst_single[0] - 1];
        
        for(k = 0; k < row_n_clst[i_col_temp]; k++){
          log_lik_comp = log_lik_comp - 
            (lgamma(gamma_mfm + row_size(k,i_col_temp)) - lgamma(gamma_mfm));
        }
        log_lik_comp = log_lik_comp - log_Vnt_last[row_n_clst[i_col_temp] - 1];
        
        // normalizing
        col_log_lik_temp[0] = 0;
        col_log_lik_temp[1] = log_lik_comp;
        log_lik_max = max(col_log_lik_temp);
        col_log_lik_temp[0] = exp(col_log_lik_temp[0] - log_lik_max);
        col_log_lik_temp[1] = exp(col_log_lik_temp[1] - log_lik_max);
        
        range = seq(0,1);
        
        ind = sample(range,1,true,col_log_lik_temp);
        k = ind[0];
        
        if(k == 1){
          // assign the sampled row-wise label back into the labelling matrix
          row_assign(_,i_col_temp) = row_assign_single(_,0);
          row_size(_,i_col_temp) = row_size_single(_,0);
          row_n_clst[i_col_temp] = row_n_clst_single[0];
          continue;
        }
      }else{
          // update column assign, use log_Vnt(n_col-2):
          log_Vnt_last=log_Vnt(n_col-1);
          
          i_col_remove(i_col, n_col, n_row,
                       col_assign, col_size, col_n_clst,
                       row_assign, row_size, row_n_clst);
          
          std::fill(col_log_lik.begin(), col_log_lik.end(), 0);
          
          // evaluate the marginal likelihood given the row partition
          for(i_col_k=0; i_col_k < col_n_clst[0]; i_col_k++){
            
            row_assign_temp=row_assign(_,i_col_k);
            
            // record the current column assign
            col_assign_temp = col_assign;
            col_assign_union = col_assign;
            
            // put this column to the column cluster
            col_assign_union[i_col] = i_col_k + 1;
            
            // get the index of i_col_k+1 - th column cluster
            col_assign_temp_idx=arma::find(col_assign_temp==(i_col_k+1));
            col_assign_union_idx=arma::find(col_assign_union==(i_col_k+1));
            
            // sample row clusters given col_assign_union
            y_mat_sub = y_mat.cols(col_assign_union_idx);
            
            // when the number of row clusters in the old column cluster reaches the maximum limit
            // need to resample the old row clustering structure.
            
            col_log_lik[i_col_k] = log_ml_bin(y_mat, col_assign_temp_idx, col_assign_union_idx,
                                              row_n_clst[i_col_k], row_assign_temp,
                                              a0, b0) + 
                                                log(col_size[i_col_k]+gamma_mfm);
          }        
          
          // new cluster term
          col_log_lik[col_n_clst[0]] = col_ml[i_col] + 
            log_Vnt_last[col_n_clst[0]]-log_Vnt_last[col_n_clst[0]-1]+log(gamma_mfm);
          
          range=seq(0,col_n_clst[0]);
          
          log_lik_temp=col_log_lik[range];
          // sample the cluster
          log_lik_max=max(log_lik_temp);
          for (k = 0; k < (col_n_clst[0]+1); k++) {
            log_lik_temp[k]=exp(log_lik_temp[k]-log_lik_max);
          }
          ind=sample(range,1,true,log_lik_temp);
          k=ind[0];
          
          col_assign[i_col]=k+1;
          if(k < col_n_clst[0]){
            ++col_size[k];
          } else {
            std::fill(partition_row.begin(), partition_row.end(), 1);
            std::fill(size_row.begin(), size_row.end(), 0);
            std::fill(clst_row.begin(), clst_row.end(), 0);
            col_size[k]=1;
            ++col_n_clst[0];
            
            size_row[0] = n_row;
            clst_row[0] = 1;
            row_assign(_,k)=partition_row;
            row_size(_,k)=size_row;
            row_n_clst[k]=clst_row[0];
          }      
        }
      }
    
    // conditioning on column assignment, do replication:
    col_assign_temp = col_assign;// record the current column assign
    for(i_rep=0; i_rep < n_rep; i_rep++){
      for(i_row=0; i_row < n_row; i_row++){
        i_row_remove(i_row, n_row, col_n_clst, row_assign, row_size, row_n_clst);
        
        // for each column cluster, sample a new row cluster:
        for(i_col=0; i_col < col_n_clst[0]; i_col++){
          
          col_assign_temp_idx = arma::find(col_assign_temp==(i_col+1));
          logVntIdx = (col_assign_temp_idx.n_elem+1)/2;
          
          // update row assign, use log_Vnt(n_row-2):
          log_Vnt_last=log_Vnt_Cons_K(logVntIdx-1);
          
          y_sub = y_mat.cols(col_assign_temp_idx);
          // get the i-th row for the marginal likelihood
          y_i = y_sub.row(i_row);
          
          std::fill(row_log_lik.begin(), row_log_lik.end(), 0);
          
          row_assign_temp = row_assign(_,i_col);
          row_assign_union = row_assign(_,i_col);
          
          
          for(i_row_k=0; i_row_k < row_n_clst[i_col]; i_row_k++){
            
            row_assign_temp_idx = arma::find(row_assign_temp == (i_row_k + 1));
            
            row_assign_union[i_row] = i_row_k + 1;
            row_assign_union_idx = arma::find(row_assign_union == (i_row_k + 1));
            
            row_log_lik[i_row_k] = row_log_lik[i_row_k] + 
              log_ml_bin_row(y_sub,
                             row_assign_temp_idx, 
                             row_assign_union_idx,
                             a0, b0) + 
                               log(row_size(i_row_k,i_col)+gamma_mfm);
            
          }
          
          if(row_n_clst[i_col] < logVntIdx){
            row_log_lik[row_n_clst[i_col]] = log_ml_bin_row_i(y_i, a0, b0)+
              log_Vnt_last[row_n_clst[i_col]] - log_Vnt_last[row_n_clst[i_col]-1] + 
              log(gamma_mfm);
            range = seq(0,row_n_clst[i_col]);
            
            log_lik_temp = row_log_lik[range];
            log_lik_max = max(log_lik_temp);
            
            for(i_row_k=0; i_row_k < (row_n_clst[i_col]+1); i_row_k++){
              log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
            }
            
          }else{
            range = seq(0,row_n_clst[i_col]-1);
            
            log_lik_temp = row_log_lik[range];
            log_lik_max = max(log_lik_temp);
            
            for(i_row_k=0; i_row_k < (row_n_clst[i_col]); i_row_k++){
              log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
            }
            
          }
          
          ind=sample(range,1,true,log_lik_temp);
          k=ind[0];
          
          row_assign(i_row,i_col)=k+1;
          if(k < row_n_clst[i_col]){
            ++row_size(k,i_col);
          } else {
            row_size(row_n_clst[i_col],i_col)=1;
            ++row_n_clst[i_col];
          }        
        }
      }      
    }
  }
}

// [[Rcpp::export]]
void consBinClstRow(arma::mat y_mat, arma::field<arma::vec> log_Vnt_Cons_K,
                    NumericVector col_assign, NumericMatrix row_assign,
                    NumericMatrix row_size,
                    NumericVector col_n_clst, NumericVector row_n_clst,
                    double a0, double b0,
                    int n_iter, int n_rep, double gamma_mfm){
  // col_assign is the column assignment, of which the length should be n_col
  // row_assign is a n_row * n_col matrix, with the i-th column recording the row assignment of i-th column cluster
  // col_size is a vector,of which the length is n_col, with the first col_n_clst[0] entries being non-zero
  // row_size is a n_row * n_col matrix, with the i-th column recording the i-th column cluster size
  // col_n_clst is a length 1 vector, with the first entry recording the number of column clusters
  // row_n_clst is a length n_col vector, with the first col_n_clst[0] entries recording the number of row clusters within each column cluster
  int n_row=y_mat.n_rows,n_col=y_mat.n_cols;
  int i_iter, i_row, i_col, i_rep;
  int i_row_k, k;
  int logVntIdx;
  
  IntegerVector range;
  IntegerVector ind;
  // for column reassignment
  arma::vec col_assign_temp, col_assign_union, row_assign_temp, row_assign_union;// inorder to convert NumericVector to arma::vec
  arma::uvec col_assign_temp_idx, col_assign_union_idx, row_assign_temp_idx, row_assign_union_idx;// record the index of a specific column and row cluster
  
  arma::mat y_sub, y_sub_sum;// submatrix of y given column index and row index
  
  NumericVector col_log_lik(n_col+0L,0);
  NumericVector log_lik_temp;
  double log_lik_max;
  arma::vec log_Vnt_last;
  NumericVector partition_row(n_row+0L,0), size_row(n_row+0L,0), clst_row(1+0L,0);
  // for row reassignment
  arma::mat y_i;
  NumericVector row_log_lik(n_row+0L,0);
  
  arma::mat y_mat_sub, y_mat_subPart;
  
  for(i_iter=0; i_iter < n_iter; i_iter++){
    
    // conditioning  on column assignment, do replication:
    col_assign_temp=col_assign;// record the current column assign
    for(i_rep=0; i_rep < n_rep; i_rep++){
      for(i_row=0; i_row < n_row; i_row++){
        i_row_remove(i_row, n_row, col_n_clst, row_assign, row_size, row_n_clst);
        
        // for each column cluster, sample a new row cluster:
        for(i_col=0; i_col < col_n_clst[0]; i_col++){
          
          col_assign_temp_idx = arma::find(col_assign_temp==(i_col+1));
          logVntIdx = (col_assign_temp_idx.n_elem+1)/2;
          
          // update row assign, use log_Vnt(n_row-2):
          log_Vnt_last=log_Vnt_Cons_K(logVntIdx-1);
          
          y_sub = y_mat.cols(col_assign_temp_idx);
          // get the i-th row for the marginal likelihood
          y_i = y_sub.row(i_row);
          
          std::fill(row_log_lik.begin(), row_log_lik.end(), 0);
          
          row_assign_temp = row_assign(_,i_col);
          row_assign_union = row_assign(_,i_col);
          
          
          for(i_row_k=0; i_row_k < row_n_clst[i_col]; i_row_k++){
            
            row_assign_temp_idx = arma::find(row_assign_temp == (i_row_k + 1));
            
            row_assign_union[i_row] = i_row_k + 1;
            row_assign_union_idx = arma::find(row_assign_union == (i_row_k + 1));
            
            row_log_lik[i_row_k] = row_log_lik[i_row_k] + 
              log_ml_bin_row(y_sub,
                             row_assign_temp_idx, 
                             row_assign_union_idx,
                             a0, b0) + 
                               log(row_size(i_row_k,i_col)+gamma_mfm);
            
          }
          
          if(row_n_clst[i_col] < logVntIdx){
            row_log_lik[row_n_clst[i_col]] = log_ml_bin_row_i(y_i, a0, b0)+
              log_Vnt_last[row_n_clst[i_col]] - log_Vnt_last[row_n_clst[i_col]-1] + 
              log(gamma_mfm);
            range = seq(0,row_n_clst[i_col]);
            
            log_lik_temp = row_log_lik[range];
            log_lik_max = max(log_lik_temp);
            
            for(i_row_k=0; i_row_k < (row_n_clst[i_col]+1); i_row_k++){
              log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
            }
            
          }else{
            range = seq(0,row_n_clst[i_col]-1);
            
            log_lik_temp = row_log_lik[range];
            log_lik_max = max(log_lik_temp);
            
            for(i_row_k=0; i_row_k < (row_n_clst[i_col]); i_row_k++){
              log_lik_temp[i_row_k] = exp(log_lik_temp[i_row_k]-log_lik_max);
            }
            
          }
          
          ind=sample(range,1,true,log_lik_temp);
          k=ind[0];
          
          row_assign(i_row,i_col)=k+1;
          if(k < row_n_clst[i_col]){
            ++row_size(k,i_col);
          } else {
            row_size(row_n_clst[i_col],i_col)=1;
            ++row_n_clst[i_col];
          }        
        }
      }      
    }
  }
}