// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


//' Graph Laplacian calculation
//' 
//' Calculate graph Laplacian of a symmetrix matrix
//' 
//' @param A symmetric matrix
//' @export
// [[Rcpp::export]]
arma::mat norm_laplacian(arma::mat A) {
    A = exp(-A/A.max());
    arma::rowvec D_row = pow(sum(A), -0.5);
    A.each_row() %= D_row; // here is Schur product: element-wise multiplication of two objects ! not equal to that in python
    arma::colvec D_col = arma::conv_to< arma::colvec >::from(D_row);
    A.each_col() %= D_col;
    arma::mat res = arma::eye(A.n_cols, A.n_cols) - A;
    return(res);
}

//' Matrix left-multiplied by its transpose
//' 
//' Given matrix A, the procedure returns A'A.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
    return(x.t()*x);
}

//' Compute Euclidean distance matrix by columns
//' 
//' Used in sc3-funcs.R distance matrix calculation
//' and within the consensus clustering.
//' 
//' @param x A numeric matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix & x) {
	unsigned int outcols = x.ncol(), i = 0, j = 0;
	double d;
	Rcpp::NumericMatrix out(outcols, outcols);

	for (j = 0; j < outcols - 1; j++) {
	    Rcpp::NumericVector v1 = x.column(j);
		for (i = j + 1; i < outcols; i++) {
			d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
			out(i, j) = d;
			out(j, i) = d;
		}
	}

	return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
