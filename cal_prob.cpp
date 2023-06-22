// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace Rcpp;


// inputs
// 	arma::mat zlimits_lower = t(zlimits[ , (CnsIndx+1):D, 1]) -- (D-CnsIndx) by N
// 	arma::mat zlimits_upper = t(zlimits[ , (CnsIndx+1):D, 2]) -- (D-CnsIndx) by N
//  arma::mat mu_noncns = t(mu.noncns) -- (D-CnsIndx) by N
// 	arma::mat Sigma_noncns_sqrt = sqrt(Sigma.noncns) -- (D-CnsIndx) by (D-CnsIndx)

// 	int Nnorms

// 	arma::colvec ?
// 	double ?

// outputs
// arma::colvec tau.noncns

// [[Rcpp::export]]
List cal_prob (
	arma::mat zlimits_lower, arma::mat zlimits_upper,
	arma::mat mu_noncns, arma::mat Sigma_noncns_sqrt,
	int Nnorms, int dim, int sampleSize
) {

	int i = 0, j = 0, D = dim, N = sampleSize ;
	int ind1 = 0, ind2 = 0 ;
	double prob = 0 ;
	arma::mat samples(D, Nnorms) ;
	arma::mat transformed_samples(D, Nnorms) ;
	arma::colvec z_l(D) ;
	arma::colvec z_u(D) ;
	arma::colvec oneSample(D) ;
	arma::colvec inds(Nnorms) ;
	// output
	arma::vec resprob(N) ;

	// generate sample
	samples = arma::randn(D, Nnorms) ;
	transformed_samples = Sigma_noncns_sqrt * samples ;	
	for (i = 0 ; i < N ; i++) {
		// lower limit and upper limit
		z_l = zlimits_lower.col(i) - mu_noncns.col(i) ;
		z_u = zlimits_upper.col(i) - mu_noncns.col(i) ;

		// truncate		
		for (j = 0 ; j < Nnorms ; j++) {
			oneSample = transformed_samples.col(j) ;
			ind1 = all(z_l < oneSample) ;
			ind2 = all(oneSample < z_u) ;
			inds(j) = (ind1 * ind2) ;
		}
		prob = mean(inds) ;
		resprob(i) = prob ;
	}
	
	List res ;
	res["prob"] = resprob ;
	return(res) ;
}