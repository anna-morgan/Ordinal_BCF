#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec rmvt_cpp(arma::mat Sigma,
                   arma::vec delta, 
                   int df){ // Generate random MVT vector
  arma::rowvec delta_r = arma::conv_to<arma::rowvec>::from(delta);
  // sample MVN
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  Y = Y * arma::chol(Sigma);
  // Sample chi-square df
  double u = R::rchisq(df);
  // Compute X (MVT)
  return(arma::conv_to<arma::vec>::from(delta_r + sqrt(df / u) * Y));
}

// [[Rcpp::export]]
double dmvt_cpp(arma::vec x,
                    arma::mat Sigma, 
                    arma::vec delta, 
                    int df,
                    bool log_bool){ // Evaluate MVT density at x (log or not)
  int p = Sigma.n_cols;
  double d1 = tgamma((df + p)/2.0) / tgamma(df/2.0) / pow(df, p/2.0) / 
                pow(arma::datum::pi, p / 2.0) / sqrt(det(Sigma));
  arma::mat m1 = arma::conv_to<arma::mat>::from(x - delta);
  arma::mat m2 = trans(m1) * inv(Sigma) * m1;
  double d2 = pow(1.0 + arma::conv_to<double>::from(m2) / df, -(df + p) / 2.0);
  if(log_bool == true){
   return(log(d1 * d2));
  } else {
   return(d1 * d2);
  }  
}

// [[Rcpp::export]]
arma::mat Hessian_cpp(arma::vec y, 
                      arma::vec mu_vect, 
                      arma::vec cuts) { // compute hessian of joint likelihood
  int J = cuts.size();
  int N = y.size();
  arma::mat hessian(J - 2, J - 2);
  for(int k = 2; k <= (J - 1); k++){
    double diag_entry = 0.0;
    for(int i = 0; i < N; i++){
      if(y[i] == k){
        double v1 = cuts(k - 1) - mu_vect(i);
        double v2 = cuts(k - 2) - mu_vect(i);
        double v3 = R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false);
        double comp = (- v3 * R::dnorm(v1, 0.0, 1.0, false) * v1 - 
                       pow(R::dnorm(v1, 0.0, 1.0, false), 2.0)) / pow(v3, 2.0);
        diag_entry = diag_entry + comp;
      }
      if(y[i] == k + 1){
        double v1 = cuts(k) - mu_vect(i);
        double v2 = cuts(k - 1) - mu_vect(i);
        double v3 = R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false);
        double comp = (v3 * R::dnorm(v2, 0.0, 1.0, false) * v2 - 
                       pow(R::dnorm(v2, 0.0, 1.0, false), 2.0)) / pow(v3, 2.0);
        diag_entry = diag_entry + comp;
      }
    }
    hessian(k-2, k-2) = diag_entry;
  }q
  for(int k = 3; k <= (J - 1); k++){
    double diag_entry = 0.0;
    for(int i = 0; i < N; i++){
      if(y[i] == k){
        double v1 = cuts(k - 1) - mu_vect(i);
        double v2 = cuts(k - 2) - mu_vect(i);
        double comp = R::dnorm(v1, 0.0, 1.0, false) * R::dnorm(v2, 0.0, 1.0, false) /
          pow(R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false), 2);
        diag_entry = diag_entry + comp;
      }
    }
    hessian(k-3, k-2) = diag_entry;
  }
  for(int k = 2; k <= (J - 2); k++){
    double diag_entry = 0.0;
    for(int i = 0; i < N; i++){
      if(y[i] == k + 1){
        double v1 = cuts(k) - mu_vect(i);
        double v2 = cuts(k - 1) - mu_vect(i);
        double comp = R::dnorm(v1, 0.0, 1.0, false) * R::dnorm(v2, 0.0, 1.0, false) / 
          pow(R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false), 2);
        diag_entry = diag_entry + comp;
      }
    }
    hessian(k-1, k-2) = diag_entry;
  }
  
  return hessian;
}

// [[Rcpp::export]]
arma::rowvec Gradient_cpp(arma::vec y, 
                          arma::vec mu_vect, 
                          arma::vec cuts,
                          bool d_prior) { // compute gradient of joint likelihood
  int J = cuts.size();
  int N = y.size();
  arma::rowvec gradient(J - 2);
  for(int k = 2; k <= (J - 1); k++){
    double grad_entry = 0;
    for(int i = 0; i < N; i++){
      if(y[i] == k){
        double v1 = cuts(k - 1) - mu_vect(i);
        double v2 = cuts(k - 2) - mu_vect(i);
        grad_entry = grad_entry + R::dnorm(v1, 0.0, 1.0, false) /
          (R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false));
      }
      
      if(y[i] == k + 1){
        double v1 = cuts(k) - mu_vect(i);
        double v2 = cuts(k - 1) - mu_vect(i);
        grad_entry = grad_entry - R::dnorm(v2, 0.0, 1.0, false) /
          (R::pnorm(v1, 0.0, 1.0, true, false) - R::pnorm(v2, 0.0, 1.0, true, false));
      }
    }
    
    if(d_prior == true){
      arma::vec alpha = 2.0 * arma::vec(J, arma::fill::ones);
      grad_entry = grad_entry + (alpha(k - 1) - 1) * R::dnorm(cuts(k - 1), 0.0, 1.0, false) / 
        (R::pnorm(cuts(k - 1), 0.0, 1.0, true, false) - R::pnorm(cuts(k - 2), 0.0, 1.0, true, false));
    }
    
    gradient(k - 2) = grad_entry;
  }
  return gradient;
}

// [[Rcpp::export]]
double Jacob_dir (arma::vec cuts_free){
  int J = cuts_free.size() + 2;
  arma::vec cuts_ext = arma::shift(arma::resize(cuts_free, J-1, 1), 1);
  arma::mat Jacob = arma::mat(J, J, arma::fill::zeros);
  for(int j = 0; j <= J - 2; j++){
    Jacob(j, 0) = 1;
    Jacob(j, j + 1) = R::pnorm(cuts_ext(j), 0.0, 1.0, true, false);
    Jacob(j + 1, j + 1) = - R::pnorm(cuts_ext(j), 0.0, 1.0, true, false);
  }
  Jacob(J - 1, 0) = 1;
  return(log(abs(arma::det(Jacob))));
}

// [[Rcpp::export]]
double num_deriv_Jacob(arma::vec cuts_free,
                       int ind){
  double h = sqrt(DBL_EPSILON);
  arma::vec cuts_free_p = cuts_free;
  cuts_free_p(ind) = cuts_free(ind) + h;
  cuts_free(ind) = cuts_free(ind) - h;
  return((Jacob_dir(cuts_free_p) - Jacob_dir(cuts_free)) / (2 * h));
}

// [[Rcpp::export]]
double second_num_deriv_Jacob(arma::vec cuts_free,
                              int ind1,
                              int ind2){
  double h = sqrt(DBL_EPSILON);
  arma::vec cuts_free_p = cuts_free;
  cuts_free_p(ind2) = cuts_free(ind2) + h;
  cuts_free(ind2) = cuts_free(ind2) - h;
  return((num_deriv_Jacob(cuts_free_p, ind1) - num_deriv_Jacob(cuts_free, ind1)) / (2 * h));
}

// [[Rcpp::export]]
double log_post_cpp (arma::vec cuts_free,
                        arma::vec y, 
                        arma::vec mu_vect,
                        bool d_prior) { // compute log posterior of joint likelihood
  int J = cuts_free.size() + 2;
  int N = y.size();
  double lp = 0;
  for(int i = 0; i < N; i++){
    if(y(i) == 2){
      double a = log(R::pnorm(cuts_free(0) - mu_vect(i), 0.0, 1.0, true, false) -
                     R::pnorm(0 - mu_vect(i), 0.0, 1.0, true, false));
      lp = lp + a;
    }
    
    for(int j = 3; j <= (J - 1); j++){
      if(y(i) == j){
        double a = log(R::pnorm(cuts_free(j - 2) - mu_vect(i), 0.0, 1.0, true, false) -
                       R::pnorm(cuts_free(j - 3) - mu_vect(i), 0.0, 1.0, true, false));
        lp = lp + a;
      }
    }
    
    if(y(i) == J){
      double a = log(1 - R::pnorm(cuts_free(J - 3) - mu_vect(i), 0.0, 1.0, true, false));
      lp = lp + a;
    }
  }
  if(d_prior == true){
    arma::vec cuts_ext = arma::shift(arma::resize(cuts_free, J-1, 1), 1);
    arma::vec alpha = 2.0 * arma::vec(J, arma::fill::ones);
    
    lp = lp + Jacob_dir(cuts_free);
     
     double comp = 0;
     comp = comp + (alpha(0) - 1) * log(R::pnorm(cuts_ext(0), 0.0, 1.0, true, false));
     for(int k = 1; k < J - 1; k++){
       comp = comp + alpha(k) * log(R::pnorm(cuts_ext(k), 0.0, 1.0, true, false) - R::pnorm(cuts_ext(k - 1), 0.0, 1.0, true, false));
     }
     comp = comp + alpha(J - 1) * log(1 - R::pnorm(cuts_ext(J - 2), 0.0, 1.0, true, false));
     lp = lp + comp;
     return(lp);
  } else {return(lp);
  }
}

// [[Rcpp::export]]
arma::vec z_update(arma::vec cuts,
                   arma::vec y,
                   arma::vec mu_vect){ // update Z's
  int N = y.size();
  arma::vec z(N);
  for(int i = 0; i < N; i++){
    arma::vec range = {cuts(y(i) - 1), cuts(y(i))};
    double mu = mu_vect(i);
    double F_a = R::pnorm(min(range), mu, 1, true, false);
    double F_b = R::pnorm(max(range), mu, 1, true, false);
    double u = R::runif(F_a, F_b);
    z(i) = R::qnorm(u, mu, 1, true, false);
  } 
  return(z);
}

//[[Rcpp::export]]
arma::mat beta_sample(arma::mat X,
                      arma::mat B0,
                      arma::mat beta0,
                      arma::mat z){ // Sample beta's
  arma::mat B_hat = inv(inv(B0) + trans(X)*X); // covariance matrix
  arma::mat beta_hat = B_hat * ((inv(B0)*beta0 + trans(X)*z)); // mean
  // MVN sampler:
  int ncols = B0.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return arma::repmat(beta_hat, 1, 1).t() + Y * arma::chol(B_hat);
  //return (beta_hat);
}

//// [[Rcpp::export]]
// arma::vec Newt_Raph_cpp(arma::vec y, 
//                         arma::vec mu_vect, 
//                         arma::vec cuts) {  // Newton Raphson optimization to find delta hat for cutpoints
//   int J = cuts.size();
//   arma::vec cutpoints_temp = cuts;
//   int B = 10;
//   for(int i = 0; i <= B; i++){
//     arma::mat inv_hess = arma::inv(Hessian_cpp(y, mu_vect, cutpoints_temp));
//     arma::rowvec grad = Gradient_cpp(y, mu_vect, cutpoints_temp);
//     arma::mat diff = inv_hess * trans(arma::conv_to<arma::mat>::from(grad));
//     for(int j = 0; j <= (J - 3); j ++){
//     cutpoints_temp(j + 1) = cutpoints_temp(j + 1) - diff(j);
//     }
//   }
//   return(cutpoints_temp);
// }

//// [[Rcpp::export]]
//
//
//arma::vec cutpoints_Update(arma::vec y, 
//                           arma::vec mu_vect, 
//                           arma::vec cuts) { // Update cutpoints
//  int J = cuts.size(); 
//  
//  // current cutpoints (0, ..., infty):
//  arma::vec current = arma::resize(shift(cuts,(J - 1)), J-2 ,1);
//  
//  // generate proposal:
//  arma::vec g = arma::resize(shift(Newt_Raph_cpp(y, mu_vect, cuts),(J - 1)), J-2 ,1);
//  arma::vec g_full = arma::resize(g, J, 1);
//  g_full(J - 2) = pow(10.0, 9);
//  g_full = shift(g_full, 1);
//  arma::mat D = arma::inv(- Hessian_cpp(y, mu_vect, g_full));
//  double df = 100.0;
//  arma::vec proposal = rmvt_cpp(D, g, df);
//  
//  // compute acceptance prob.
//
//  arma::vec test = current; /////////////////////
//  double a = log_post_cpp(proposal, y, mu_vect) - log_post_cpp(current, y, mu_vect);
//  double b = dmvt_cpp(current, D, g, df, true);
//  double c = dmvt_cpp(proposal, D, g, df, true);
//    
//  double alpha = 0;
//  if(a + b - c <= 0){
//    alpha = a + b - c;
//  }
//  alpha = exp(alpha);
//  
//  // accept / reject proposal:
//  double p = 0;
//  arma::vec prob = {1 - alpha, alpha};///
//  arma::vec x = {0.0, 1.0};///
//  const int size = 1;///
//  if(isnan(exp(a + b - c)) == false) {
//    arma::vec prob = {1 - alpha, alpha};
//    arma::vec x = {0.0, 1.0};
//    const int size = 1;
//    p = Rcpp::as<double>(wrap(Rcpp::RcppArmadillo::sample(x, size, false, prob)));
//    if(p == 1) {
//      current = proposal;
//    }
//  }
//  
//  // resize vector to include fixed cutpoints:
//  current = arma::resize(current, J, 1);
//  current(J - 2) = pow(10.0, 9);
//  current = shift(current, 1);
//  NumericVector test2 = {alpha, p, Rcpp::as<double>(wrap(Rcpp::RcppArmadillo::sample(x, size, false, prob)))};
//  
//  return(current);
//}
//
//// [[Rcpp::export]]
//List Algorithm_2_cpp(arma::vec y, 
//                          arma::vec mu_vect, 
//                          arma::vec cuts,
//                          arma::mat X,
//                          arma::mat B0,
//                          arma::vec beta0) { // one iteration of the algorithm
//  int J = cuts.size();
// // Step 1(a):
//  arma::vec cutpoints = cutpoints_Update(y, mu_vect, cuts);
//  
// // Step 1(b):
//  arma::vec cutpoints_z = shift(arma::resize(cutpoints, J + 1, 1), 1);
//  cutpoints_z(0) = - pow(10.0, 9);
//  arma::vec z = z_update(cutpoints_z, y, mu_vect);
//    
// // Step 2:
//  arma::mat beta = beta_sample(X, B0, arma::conv_to<arma::mat>::from(beta0), 
//                               arma::conv_to<arma::mat>::from(z));
//  
// // Compute log posterior
//  double lp = log_post_cpp(resize(shift(cutpoints, J - 1), J - 2, 1), y, mu_vect); // need to update mu_vect
//
//  List L = List::create(cutpoints, z, beta, lp, 1);
//  
//  return(L);
//}
//
//// [[Rcpp::export]]
//List ordinal_Sampler_cpp(List trace_init,
//                          int M,
//                          arma::vec y,
//                          arma::mat X,
//                          arma::mat B0,
//                          arma::vec beta0){ // Run algorithm M times and save parameters
//  
//   arma::vec cuts = trace_init[0];
//   arma::vec z = trace_init[1];
//   arma::vec beta = trace_init[2];
//   double lp = trace_init[3];
//   int mh_accept = trace_init[4];
//   
//   arma::vec mus = arma::conv_to<arma::vec>::from(X * beta);
//   
//   Rcpp::List T(5 * M);
//   for(int m = 0; m < M; m++){
//     T[5 * m] = cuts;
//     T[5 * m + 1] = z;
//     T[5 * m + 2] = beta;
//     T[5 * m + 3] = lp;
//     T[5 * m + 4] = mh_accept;
//     
//     mus = arma::conv_to<arma::vec>::from(X * beta);
//     
//     List L = Algorithm_2_cpp(y, mus, cuts, X, B0, beta0);
//     cuts = as<arma::vec>(L[0]);
//     z = as<arma::vec>(L[1]);
//     beta = as<arma::vec>(L[2]);
//     lp = L[3];
//     mh_accept = L[4];
//   }
//   return(T);
// }

/*** R
 N <- 200
 J <- 4
 set.seed(2020)
 x0 <- rep(1, N)
 x1 <- rnorm(N)
 x2 <- rnorm(N, sd = 0.4)
 X <- matrix(c(x0, x1, x2), ncol = 3)
 beta_true <- c(1.1, 0.14, 0.05)
 z_true <- as.vector(X %*% beta_true) + rnorm(N, 0, 1)
 #hist(z_true)
 cutpoints_true <- c(-1000000, 0, 1.25, 2.25, 1000000)
 y <- rep(0, N)
 for(j in 1:length(cutpoints_true)){
   for(i in 1:N){
     if(z_true[i] >= cutpoints_true[j]){ 
       y[i] <- j
     }
   }
 }
 #hist(y)
 
 # Initiate parameters for sampler
 cutpoints0 <- c(0, 1, 2, 100)
 cutpoints0[4] <- 10 ^ 9
 beta0 <- c(1, .1, .06)
 B0 <- matrix(c(c(0.4,0,0), 
                c(0,0.4,0),
                c(0,0,0.4)), nrow = 3)
 Bhat <- solve(B0 + t(X) %*% X)
 
 # function arguments - temporary
 beta <- beta0
 cutpoints <- cutpoints0
 cp_vect <- cutpoints[2:(J-1)]
 beta <- as.vector(beta)
 mus <- as.vector(X %*% beta)
# 
# Hessian_cpp(y, mus, cutpoints)
# Gradient_cpp(y, mus, cutpoints)
log_post_cpp(cp_vect, y, mus, TRUE)
# cutpoints_z <- c(-1 * 10^9, cutpoints)
# z_update(cutpoints_z, y, mus)
# Newt_Raph_cpp(y, mus, cutpoints)
# cutpoints_Update(y, mus, cutpoints)
# Algorithm_2_cpp(y, mus, cutpoints, X, B0, beta0)
# 
# start_time <- Sys.time()
# T_list <- ordinal_Sampler_cpp(trace, M, y, X, B0, beta0)
# end_time <- Sys.time()
# end_time - start_time

*/

