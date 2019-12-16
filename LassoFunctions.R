# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean = mean(Y)
  Ytilde = Y - Ymean
  
  # [ToDo] Center and scale X
  Xmeans = colMeans(X)
  X_centred = t(t(X) - Xmeans)
  
  # Scale X
  weights = sqrt(colSums(X_centred * X_centred)/nrow(X_centred))
  Xtilde = X_centred %*% diag(1/weights)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  
  # calculate new beta
  new_beta = sign(a) * max(abs(a) - lambda, 0)
  
  return(new_beta)

}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  
  return((crossprod(Ytilde - Xtilde %*% beta)) /(2 * length(Ytilde)) + lambda * sum(abs(beta)))
 
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001){
  
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    
    stop('The length of X and Y cannot be different')
    
  }
  
  #[ToDo]  Check that lambda is non-negative
  if(lambda < 0){
    
    stop('Regularization parameter (lambda) cannot be negative')
    
  }
  
  #[ToDo]  Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  if(length(beta_start) == 0){
    
    beta_start = rep(0, ncol(Xtilde))
    
  }else{
    
    if(length(beta_start) != ncol(Xtilde)){
      
      stop('The number of features in beta_init and X cannot be different')
    }
  }
  
  #[ToDo]  Coordinate-descent implementation. Stop when the difference between objective functions is less than eps for the first time.
  beta_last = beta_start
  beta_new = beta_start
  loss_diff = 100
  r = Ytilde - Xtilde %*% beta_last
  # write a while loop to optimize beta until the loss function does not change larger than eps
  while (loss_diff >= eps) {
    beta_last = beta_new
    for (i in 1:ncol(Xtilde)) {
      beta_new[i] = soft((beta_last[i] + crossprod(r, Xtilde[, i])/nrow(Xtilde)), lambda )
      r = r + Xtilde[, i] * (beta_last[i] - beta_new[i])
    }
    
    loss_old = lasso(Xtilde, Ytilde, beta_last, lambda)
    loss_new = lasso(Xtilde, Ytilde, beta_new, lambda)
    loss_diff = loss_old - loss_new
    
  }
  fmin = lasso(Xtilde, Ytilde, beta_new, lambda)
  
  
  # For example, if you had 3 iterations with objectives 3, 1, 0.99999, your should return fmin = 0.99999, and not have another iteration
  
  # Return 
  # beta - the solution
  # fmin - optimal function value (value of objective at beta)
  return(list(beta = beta_new, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 50, eps = 0.0001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    
    stop('The length of X and Y cannot be different')
    
  }
  
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >=0, and make sure the values are sorted from largest to smallest. If none of the supplied values satisfy the requirement, print the warning message and proceed as if the values were not supplied.
  if(length(lambda_seq) != 0){
    
    lambda_seq = sort(lambda_seq[lambda_seq >= 0], decreasing = TRUE)
  }
  
  # If lambda_seq is not supplied, calculate lambda_max (the minimal value of lambda that gives zero solution), and create a sequence of length n_lambda as
  if(length(lambda_seq) == 0){
    
    lambda_max = max(abs((t(Xtilde) %*% Ytilde)/nrow(Xtilde)))
    
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
    
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda (make sure supplied eps is carried over). Use warm starts strategy discussed in class for setting the starting values.
  # set a container for beta_mat, the final output
  beta_mat = matrix(NA, ncol(Xtilde), length(lambda_seq))
  # set a container for fmin_vec for holding the loss value for each lambda
  fmin_vec = rep(NA, length(lambda_seq))
  
  # firstly initiate beta_start as 0's to use warm start strategy
  beta_start = rep(0, ncol(Xtilde))
  
  for (i in 1:length(lambda_seq)) {
    result = fitLASSOstandardized(Xtilde, Ytilde, lambda = lambda_seq[i], beta_start = beta_start, eps = eps)
    beta_mat[, i] = result$beta
    fmin_vec[i] = result$fmin
    # update beta hat (lambda) as a new starting point
    beta_start = result$beta
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 50, eps = 0.0001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardized = standardizeXY(X, Y)
  Xtilde = standardized$Xtilde
  Ytilde = standardized$Ytilde
  weights = standardized$weights
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
  lasso_fit = fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)
  
  # Getting lambdas
  lambda_seq = lasso_fit$lambda_seq
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector for each lambda
  beta_mat = diag(1/weights) %*% lasso_fit$beta_mat
  beta0_vec = rep(standardized$Ymean, length(lambda_seq)) - (standardized$Xmeans %*% beta_mat)
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifiying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.0001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 50, k = 5, fold_ids = NULL, eps = 0.0001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  Lasso = fitLASSO(X ,Y, lambda_seq, n_lambda, eps)
  lambda_seq = Lasso$lambda_seq
  beta_mat = Lasso$beta_mat
  beta0_vec = Lasso$beta0_vec
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds. If fold_ids is not NULL, split the data according to supplied fold_ids.
  if(length(fold_ids) == 0){
    
    fold_ids = cut(1:nrow(X), breaks = k, labels = FALSE)
    
  }
  
  # errors container of size k by len(lambda_sec)
  errors = matrix(NA, length(unique(fold_ids)), length(lambda_seq))
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO, and perform any additional calculations needed for CV(lambda) and SE_CV(lambda
  for (j in unique(fold_ids)) {
    # Xtrain and Ytrain for each fold
    Xtrain = X[fold_ids != j, ]
    Ytrain = Y[fold_ids != j]
    
    # Xtest and Ytest for each fold
    Xtest = X[fold_ids == j, ]
    Ytest = Y[fold_ids == j]
    
    result = fitLASSO(Xtrain, Ytrain, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
    result_beta = result$beta_mat
    result_beta0 = result$beta0_vec
    result_lambda = result$lambda_seq
    
    # Vectorize the calculation of losses per fold. store the results in a corresponding row in the errors_all matrix
    y_hat_matrix =  matrix(rep(Ytest, each = length(result_lambda)), length(Ytest), byrow = TRUE)
    
    # subtract the intercept from the matrix
    loss_intercept = sweep(y_hat_matrix, 2, result_beta0)
    
    # calculate X'beta in matrix form
    loss_xbeta = Xtest %*% result_beta
    errors[j, ] = colMeans( (loss_intercept - loss_xbeta) ^ 2)
  }
  
  # calculate cvlambda
  cvm = colMeans(errors)
  
  # [ToDo] Find lambda_min
  lambda_min = lambda_seq[which.min(cvm)]
  
  # calculate cvse
  cvse = apply(errors, 2, function(x) sd(x) / sqrt(length(unique(fold_ids))))
  
  # [ToDo] Find lambda_1SE
  upper_limit = cvm[which.min(cvm)] + cvse[which.min(cvm)]
  lambda_1se = lambda_seq[which.max(cvm <= upper_limit)]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

