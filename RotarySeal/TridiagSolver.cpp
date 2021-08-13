#include "TridiagSolver.h"


void TDMASolver(VECTORTYPE& x, VECTORTYPE& a, VECTORTYPE& b, VECTORTYPE& c, VECTORTYPE& d)
{
	int n= a.size()-1;
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);
	x(n) = d(n);

	for (int i = n; i-- > 0;) {
		x[i] = d[i]- c[i] * x[i + 1];
	}
}

void CTDMASolver(VECTORTYPE& x, VECTORTYPE& a, VECTORTYPE& b, VECTORTYPE& c, VECTORTYPE& d) {
	DATATYPE temp, p, q, r;
	int n = a.size() - 1;
	b[1] /= a[1];
	c[1] /= a[1];
	d[1] /= a[1];

	for (int i = 2; i <= n-2; i++) {
		temp = a[i] - c[i] * b[i - 1];
		b[i] /= temp;
		d[i] = (d[i] + c[i] * d[i - 1]) / temp;
		c[i] *= c[i - 1] / temp;
	}

	p = a[n-1];
	q = b[n-1];
	r = d[n-1];

	for (int i = 2; i <= n - 2; i++) {
		p = p - q * c[i - 1];
		r = r + q * d[i - 1];
		q = q * b[i - 1];
	}

	x[n - 1] = ((q + c[n - 1]) * d[n - 2] + r) / (p - (q + c[n - 1]) * (b[n - 2] + c[n - 2]));

	for (int i = n - 2; i >= 1; i--) {

		x[i] = b[i] * x[i + 1] + c[i] * x[n - 1] + d[i];
	}
	x[0] = x[n - 1];
	x[n] = x[1];

}