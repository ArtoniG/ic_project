#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix genes_matrix_pm(NumericMatrix data, NumericMatrix genes_pm, NumericMatrix by_gene_id){
	int nrows = genes_pm.nrow(), ncols = genes_pm.ncol(), aux;
	
	for (int i = 0; i < nrows; i++){
		for (int j = 0; j < ncols; j++){
			aux = by_gene_id(i,0);
			genes_pm(i,j) = data(aux, j);
		}
	}
	return genes_pm;
}
