#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector rowWhichMax(NumericMatrix m) {
	int nr = m.nrow();
	int nc = m.ncol();
	IntegerVector index(nr);

	for(int i = 0; i < nr; i ++) {
		double max = -99999999;
		int which_max = -1;
		for(int j = 0; j < nc; j ++) {
			if(max < m(i, j)) {
				max = m(i, j);
				which_max = j;
			}
		}
		index[i] = which_max + 1;
	}
	return index;
}
