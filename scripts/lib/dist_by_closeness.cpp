#include <Rcpp.h>
using namespace Rcpp;


float closeness(NumericVector x1, NumericVector x2) {
	float y = 0;
	int y_n = 0;
	int n = x1.size();

	for(int i = 0; i < n; i ++) {
		for(int j = 0; j < n; j ++) {
			if(std::abs(x1[i] - 0) < 1e-6 | std::abs(x2[j] - 0) < 1e-6) {
				continue;
			}
			y += std::abs(i - j);
			y_n += 1;
		}
	}

	return y/y_n;
}


int min_v(IntegerVector x) {
	int min = 1e8;

	for(int i = 0; i < x.size(); i ++) {
		if(x[i] != -1) {
			if(x[i] < min) {
				min = x[i];
			}
		}
	}
	return min;
}

float closeness_with_min_dist(NumericVector x1, NumericVector x2) {
	int n1 = x1.size();
	int n2 = x2.size();

	IntegerMatrix dist(n1, n2);

	for(int i = 0; i < n1; i ++) {
		for(int j = 0; j < n2; j ++) {
			if(std::abs(x1[i] - 0) < 1e-6 | std::abs(x2[j] - 0) < 1e-6) {
				dist(i, j) = -1;
			} else {
				dist(i, j) = std::abs(i - j);
			} 
		}
	}

	IntegerVector d1(n1);
	IntegerVector d2(n2);
	for(int i = 0; i < n1; i ++) {
		d1[i] = min_v(dist(i, _));
	}
	for(int j = 0; j < n2; j ++) {
		d2[j] = min_v(dist(_, j));
	}

	float y = 0;
	int y_n = 0;
	for(int i = 0; i < n1; i ++) {
		if(d1[i] != -1) {
			y = y + d1[i];
			y_n ++;
		}
	}
	for(int j = 0; j < n2; j ++) {
		if(d2[j] != -1) {
			y = y + d2[j];
			y_n ++;
		}
	}

	return y/y_n;
}


// [[Rcpp::export]]
NumericMatrix dist_by_closeness(NumericMatrix mat) {
	int n = mat.nrow();

	NumericMatrix dist(n, n);

	for(int i = 0; i < n - 1; i ++) {
		for(int j = i+1; j < n; j ++) {
			dist(i, j) = closeness(mat(i, _), mat(j, _));
			dist(j, i) = dist(i, j);
		}
	}
	for(int i = 0; i < n; i ++) {
		dist(i, i) = 0;
	}
	return dist;
}

