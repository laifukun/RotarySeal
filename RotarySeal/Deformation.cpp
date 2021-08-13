#include "Deformation.h"
#include <cmath>

int influenceMatrix(MATRIXTYPE& Ik, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD)
{
	DATATYPE Iii, LD2, d2, c2, sqrtdc;
	int M = y.size();
	int N = x.size();
	LD2 = LD * LD;
	d2 = d * d;
	c2 = c * c;
	sqrtdc = sqrt(LD2 * d2 + c2);
	Iii = (LD * (-d) * log(((-c) + sqrtdc) / ((+c) + sqrtdc)) + LD * (d) * log(((+c) + sqrtdc) / ((-c) + sqrtdc)) + (+c) * log((LD * (+d) + sqrtdc) / (LD * (-d) + sqrtdc)) + (-c) * log((LD * (-d) + sqrtdc) / (LD * (+d) + sqrtdc)));

	for (int i = 0; i < N; i++) 
		for (int j = 0; j < M; j++) {
			int o = j + M * i;
			Ik(o, o) = Iii;
			for (int k=0; k<N; k++)
				for (int l = 0; l < M; l++) {
					int p = l + M * k;
					if (p > o) {
						DATATYPE xikpc = x(i) - x(k) + c;
						DATATYPE xikmc = x(i) - x(k) - c;
						DATATYPE yjlpd = y(j) - y(l) + d;
						DATATYPE yjlmd = y(j) - y(l) - d;
						DATATYPE xikpc2 = xikpc* xikpc;
						DATATYPE xikmc2 = xikmc* xikmc;
						DATATYPE yjlpd2 = yjlpd* yjlpd;
						DATATYPE yjlmd2 = yjlmd* yjlmd;
						DATATYPE x1 = (xikmc + sqrt(LD2 * yjlmd2 + xikmc2)) / ((xikpc)+sqrt(LD2 * yjlmd2 + xikpc2));
						DATATYPE x2 = (xikpc + sqrt(LD2 * yjlpd2 + xikpc2)) / (xikmc + sqrt(LD2 * yjlpd2 + xikmc2));
						DATATYPE x3 = (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikpc2)) / (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikpc2));
						DATATYPE x4 = (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikmc2)) / (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikmc2));
						Ik(o, p) = LD * (yjlmd)* log(x1) + LD * (yjlpd)* log(x2) + (xikpc)* log(x3) + (xikmc)* log(x4);
						Ik(p, o) = Ik(o, p);
					}
				}
		}
	return 1;
}

int influenceMatrix(MATRIXTYPE& Ik, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& xi, const DATATYPE& yi)
{
	DATATYPE  LD2=LD*LD;
	int M = y.size();
	int N = x.size();
	DATATYPE xikpc, xikmc, yjlpd, yjlmd, x1, x2, x3, x4, xikpc2, xikmc2, yjlpd2, yjlmd2;
			for (int k = 0; k < N; k++)
				for (int l = 0; l < M; l++) {

					//if (xi - x(k) <= 20 * c && yi - y(l) <= 20 * d) {

						xikpc = xi - x(k) + c;
						xikmc = xi - x(k) - c;
						yjlpd = yi - y(l) + d;
						yjlmd = yi - y(l) - d;
						xikpc2 = xikpc * xikpc;
						xikmc2 = xikmc * xikmc;
						yjlpd2 = yjlpd * yjlpd;
						yjlmd2 = yjlmd * yjlmd;
						x1 = (xikmc + sqrt(LD2 * yjlmd2 + xikmc2)) / (xikpc + sqrt(LD2 * yjlmd2 + xikpc2));
						x2 = (xikpc + sqrt(LD2 * yjlpd2 + xikpc2)) / (xikmc + sqrt(LD2 * yjlpd2 + xikmc2));
						x3 = (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikpc2)) / (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikpc2));
						x4 = (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikmc2)) / (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikmc2));
						Ik(l, k) = LD * (yjlmd)* log(x1) + LD * (yjlpd)* log(x2) + (xikpc)* log(x3) + (xikmc)* log(x4);
					//}
					//else
					//	Ik(l, k) = 0.0;
				}

	return 1;
}

DATATYPE influenceMatrix(const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& xi, const DATATYPE& yi, const MATRIXTYPE& P)
{
	DATATYPE  LD2 = LD * LD;
	int M = y.size();
	int N = x.size();
	DATATYPE xikpc, xikmc, yjlpd, yjlmd, x1, x2, x3, x4, xikpc2, xikmc2, yjlpd2, yjlmd2;
	DATATYPE Hi = 0.0;
	for (int k = 0; k < N; k++)
		for (int l = 0; l < M; l++) {

			//if (xi - x(k) <= 20 * c && yi - y(l) <= 20 * d) {

			xikpc = xi - x(k) + c;
			xikmc = xi - x(k) - c;
			yjlpd = yi - y(l) + d;
			yjlmd = yi - y(l) - d;
			xikpc2 = xikpc * xikpc;
			xikmc2 = xikmc * xikmc;
			yjlpd2 = yjlpd * yjlpd;
			yjlmd2 = yjlmd * yjlmd;
			x1 = (xikmc + sqrt(LD2 * yjlmd2 + xikmc2)) / (xikpc + sqrt(LD2 * yjlmd2 + xikpc2));
			x2 = (xikpc + sqrt(LD2 * yjlpd2 + xikpc2)) / (xikmc + sqrt(LD2 * yjlpd2 + xikmc2));
			x3 = (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikpc2)) / (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikpc2));
			x4 = (LD * yjlmd + sqrt(LD2 * yjlmd2 + xikmc2)) / (LD * yjlpd + sqrt(LD2 * yjlpd2 + xikmc2));
			//Ik(l, k) = LD * (yjlmd)* log(x1) + LD * (yjlpd)* log(x2) + (xikpc)* log(x3) + (xikmc)* log(x4);
			Hi+= (LD * (yjlmd)* log(x1) + LD * (yjlpd)* log(x2) + (xikpc)* log(x3) + (xikmc)* log(x4))*P(l,k);
			//}
			//else
			//	Ik(l, k) = 0.0;
		}

	return Hi;
}

int Deformation(MATRIXTYPE& H, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& deformPara, const MATRIXTYPE& P)
{
	DATATYPE  LD2 = LD * LD;
	int M = y.size();
	int N = x.size();
	DATATYPE xikpc, xikmc, yjlpd, yjlmd, x1, x2, x3, x4, xikpc2, xikmc2, yjlpd2, yjlmd2,xj,yi,sqrtyxmm, sqrtyxmp, sqrtyxpm, sqrtyxpp;
	DATATYPE Iop;
	DATATYPE sqrtdc = sqrt(LD2 * d*d + c*c);
	DATATYPE Iii = (LD * (-d) * log(((-c) + sqrtdc) / ((+c) + sqrtdc)) + LD * (d)* log(((+c) + sqrtdc) / ((-c) + sqrtdc)) + (+c) * log((LD * (+d) + sqrtdc) / (LD * (-d) + sqrtdc)) + (-c) * log((LD * (-d) + sqrtdc) / (LD * (+d) + sqrtdc)));
	int o, p,i,j,k,l;
	//for (i = 0; i < N; i++)
		//for (j = 0; j < M; j++) {
			H = Iii * P;
		//}
	//for (int i = 0; i < N; i++)
		//for (int j = 0; j < M; j++) 
	//		int im = (N - 1) / 2;
	for (o = 0; o < M*N; o++) {
			
		i = o % N;
		j = o / N;
			xj = x[i];
			yi = y[j];
			//o = j + M * i;
			//for (int k = 0; k < N; k++)
				//for (int l = 0; l < M; l++) 
			for (p = o+1; p < M*N; p++) {

				//	p = l + M * k;
				k = p % N;
				l = p / N;
				//	if (p > o) {
						xikpc = xj - x(k) + c;
						xikmc = xj - x(k) - c;
						yjlpd = yi - y(l) + d;
						yjlmd = yi - y(l) - d;
						//xikpc2 = xikpc * xikpc;
						//xikmc2 = xikmc * xikmc;
						//yjlpd2 = yjlpd * yjlpd;
						//yjlmd2 = yjlmd * yjlmd;
						sqrtyxmm=sqrt(LD2 * yjlmd * yjlmd + xikmc * xikmc);
						sqrtyxmp=sqrt(LD2 * yjlmd * yjlmd + xikpc * xikpc);
						sqrtyxpm = sqrt(LD2 * yjlpd * yjlpd + xikmc * xikmc);
						sqrtyxpp = sqrt(LD2 * yjlpd * yjlpd + xikpc * xikpc);
						x1 = (xikmc + sqrtyxmm) / (xikpc + sqrtyxmp);
						x2 = (xikpc + sqrtyxpp) / (xikmc + sqrtyxpm);
						x3 = (LD * yjlpd + sqrtyxpp) / (LD * yjlmd + sqrtyxmp);
						x4 = (LD * yjlmd + sqrtyxmm) / (LD * yjlpd + sqrtyxpm);

						Iop = LD * (yjlmd* log(x1) + yjlpd* log(x2)) + xikpc* log(x3) + xikmc* log(x4);
						H(j, i) += Iop * P(l, k);
						H(l, k) += Iop * P(j, i);
				//	}

				}
		}
	
	H *= deformPara;
//	H += 
//	for (i = 0; i < N; i++)
//		for (j = 0; j < M; j++) {
//			H(j, i) = H(j, i)/Ebar + A1;
//		}
	return 1;
	
}

int DeformationPolar(MATRIXTYPE& H, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& deformPara, const MATRIXTYPE& P)
{
	DATATYPE  LD2 = LD * LD;
	int M = y.size();
	int N = x.size();
	DATATYPE xikpc, xikmc, yjlpd, yjlmd, x1, x2, x3, x4, xikpc2, xikmc2, yjlpd2, yjlmd2, xj, yi, sqrtyxmm, sqrtyxmp, sqrtyxpm, sqrtyxpp;
	DATATYPE Iop;
	DATATYPE sqrtdc = sqrt(LD2 * d * d + c * c);
	DATATYPE Iii = 0.0;
		//(LD * (-d) * log(((-c) + sqrtdc) / ((+c) + sqrtdc)) + LD * (d)* log(((+c) + sqrtdc) / ((-c) + sqrtdc)) + (+c) * log((LD * (+d) + sqrtdc) / (LD * (-d) + sqrtdc)) + (-c) * log((LD * (-d) + sqrtdc) / (LD * (+d) + sqrtdc)));
	int o, p, i, j, k, l;
	//for (i = 0; i < N; i++)
		//for (j = 0; j < M; j++) {
	H = Iii * P;
	//}
	int im = (N - 1) / 2;
	//for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
	
	//for (o = 0; o < M * N; o++) {

		//i = o % N;
		//j = o / N;
		xj = x[im];
		yi = y[j];
		//o = j + M * i;
		for (int k = 0; k < N; k++)
			for (int l = 0; l < M; l++) {
		//for (p = 0; p < M * N; p++) {

			//	p = l + M * k;
			//k = p % N;
			//l = p / N;
			//	if (p > o) {
			xikpc = xj - x(k) + c;
			xikmc = xj - x(k) - c;
			yjlpd = yi - y(l) + d;
			yjlmd = yi - y(l) - d;
			//xikpc2 = xikpc * xikpc;
			//xikmc2 = xikmc * xikmc;
			//yjlpd2 = yjlpd * yjlpd;
			//yjlmd2 = yjlmd * yjlmd;
			sqrtyxmm = sqrt(LD2 * yjlmd * yjlmd + xikmc * xikmc);
			sqrtyxmp = sqrt(LD2 * yjlmd * yjlmd + xikpc * xikpc);
			sqrtyxpm = sqrt(LD2 * yjlpd * yjlpd + xikmc * xikmc);
			sqrtyxpp = sqrt(LD2 * yjlpd * yjlpd + xikpc * xikpc);
			x1 = (xikmc + sqrtyxmm) / (xikpc + sqrtyxmp);
			x2 = (xikpc + sqrtyxpp) / (xikmc + sqrtyxpm);
			x3 = (LD * yjlpd + sqrtyxpp) / (LD * yjlmd + sqrtyxmp);
			x4 = (LD * yjlmd + sqrtyxmm) / (LD * yjlpd + sqrtyxpm);

			Iop = LD * (yjlmd * log(x1) + yjlpd * log(x2)) + xikpc * log(x3) + xikmc * log(x4);
			H(j, im) += Iop * P(l, k)* deformPara;
			//H(l, k) += Iop * P(j, i);
	//	}

		}
	}
		for (int i = 0; i < N; i++)
			//for (int j = 0; j < M; j++)
				H.col(i) = H.col(im);
	//H *= deformPara;
	//	H += 
	//	for (i = 0; i < N; i++)
	//		for (j = 0; j < M; j++) {
	//			H(j, i) = H(j, i)/Ebar + A1;
	//		}
	return 1;

}

int asperityAngle(MATRIXTYPE& theta, const DATATYPE& ybbar, const DATATYPE& eps, const DATATYPE& alpha, const VECTORTYPE& y) {
	int M = theta.rows();
	int N = theta.cols();

	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			if (y[j] <= ybbar)
				theta(j, i) = atan(1.0/(2.0*alpha*eps*(ybbar-y(j))/ybbar/ybbar));
			else
				theta(j,i)= atan(1.0 / (2.0 * alpha * eps * (ybbar - y(j)) /(1.0-ybbar) / (1.0-ybbar)));
		}

	return 1;

}
