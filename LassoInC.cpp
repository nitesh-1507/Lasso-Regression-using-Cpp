#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  if (a > lambda){
    
    return(a-lambda);
    
  }else if (a < -lambda){
    
    return(a + lambda);
    
  }else{
    
    return(0);
    
  }
}


// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  
  // Number of data points
  int n = Xtilde.n_rows;
  
  // Initialize parameters
  arma::mat X = Xtilde;
  arma::colvec Y = Ytilde;
  arma::colvec param = beta;
  
  // Calculate and return objective value
  double obj = accu(square(Y - X * param)/(2 * n)) + lambda * accu(abs(param));
  return obj;
  
}

// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.0001){
  // Your function code goes here
  // Create beta_last and beta_new vectors
  int row = Xtilde.n_rows;
  int col = Xtilde.n_cols;
  arma::colvec beta_last = beta_start;
  arma::colvec beta_new = beta_start;
  
  // Loss difference variable
  double loss_diff = 100;
  
  // Full initial residual
  arma::colvec r = Ytilde - Xtilde * beta_start;
  
  // Update beta
  while (loss_diff >= eps) {
    
    beta_last = beta_new;
    
    // Objective value old
    double loss_old = lasso_c(Xtilde, Ytilde, beta_last, lambda);
    
    for(int i = 0; i < col; i++)
    {
      // Update beta
      beta_new(i) = soft_c(arma::as_scalar(beta_last(i) + (Xtilde.col(i).t() * r / row)), lambda);
      
      // Partial residual
      r = r + Xtilde.col(i) * (beta_last(i) - beta_new(i));
    }
    
    // Objective value new
    double loss_new = lasso_c(Xtilde, Ytilde, beta_new, lambda);
    
    // Objective value difference
    loss_diff = loss_old - loss_new;
  }
  
  // Return updated beta
  return beta_new;
  
}  

// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.0001){
  // Your function code goes here
  // Creates beta matrix 
  int row = Xtilde.n_cols;
  int col = lambda_seq.size();
  arma::mat beta_mat(row, col);
  
  // Initialize beta start
  arma::colvec beta_start(row);
  beta_start.fill(0.0);
  
  // Calculate beta for each lambda
  for(int i = 0; i < col; i++)
  {
    
    arma::colvec beta_new = fitLASSOstandardized_c(Xtilde,  Ytilde, lambda_seq(i), beta_start, eps);
    beta_mat.col(i) = beta_new;
    beta_start = beta_new;
    
  }
  
  // Return beta matrix
  return beta_mat;
  
}