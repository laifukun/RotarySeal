#include "flowFactor.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

DATATYPE shearFlowFactorEquation(const int& Equ, const DATATYPE& m, const DATATYPE& c2, const DATATYPE& Fai, const DATATYPE& H) {
	//DATATYPE fsc;
	if (H > 1) {
		switch (Equ) {
		case 1: return -m * Fai / H / H + c2; 
		case 2: return 1.0 / (c2 + m * sqrt(Fai) / H);
		case 3: return 1.0 / (c2 + m * Fai / H / H);
		default: return c2; 
		}
	}
	else {
		switch (Equ) {
		case 1: return c2 - m*Fai;
		case 2: return 1.0 / (c2 + m * sqrt(Fai));
		case 3: return 1.0 / (c2 + m * Fai);
		default: return c2;
		}
	}
	
	

}
DATATYPE densityFlowFactor(const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai) {
	int iFC0, iFC90;
	DATATYPE a2, b2, fc0, fc90, FaiH;
	DATATYPE fcc;
	std::vector<DATATYPE>::const_iterator it = find(FCgam.begin(), FCgam.end(), gam);
	if (it != FCgam.end()) {
		iFC0 = std::distance(FCgam.begin(),it);
		iFC90 = 2 + (iFC0 - 2) * (-1);
		if (Fai > 0) 
			FaiH = pow(Fai, -1.5) * H*H*H*H*H*H;
		else
			FaiH = 0;
		a2 = FCa2.at(iFC0);
		b2 = FCb2.at(iFC0);
		fc0 = 1.0 / (a2 + b2 * FaiH);

		a2 = FCa2.at(iFC90);
		b2 = FCb2.at(iFC90);
		fc90 = 1.0 / (a2 + b2 * FaiH);
		DATATYPE cost = cos(theta), sint = sin(theta);
		fcc = cost*cost * fc0 + sint*sint * fc90;
	}
	else
		return 1.0;

	return fcc;

}

int densityFlowFactor(MATRIXTYPE& fcc, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& H, const MATRIXTYPE& Fai) {
	int iFC0, iFC90;
	DATATYPE a2, b2, fc0, fc90, FaiH, cost, sint;
	//DATATYPE fcc;
	int M = fcc.rows(), N = fcc.cols();
	std::vector<DATATYPE>::const_iterator it = find(FCgam.begin(), FCgam.end(), gam);
	if (it != FCgam.end()) {
		iFC0 = std::distance(FCgam.begin(), it);
		iFC90 = 2 + (iFC0 - 2) * (-1);
		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
			
				FaiH = Fai(j, i) > 0 ? pow(Fai(j, i), -1.5) * H(j, i) * H(j, i) * H(j, i) * H(j, i) * H(j, i) * H(j, i) : 0.0;

				a2 = FCa2.at(iFC0);
				b2 = FCb2.at(iFC0);
				fc0 = 1.0 / (a2 + b2 * FaiH);

				a2 = FCa2.at(iFC90);
				b2 = FCb2.at(iFC90);
				fc90 = 1.0 / (a2 + b2 * FaiH);
				cost = costheta(j,i), sint = sintheta(j,i);
				fcc(j,i) = cost * cost * fc0 + sint * sint * fc90;
			}
			
	}
	else
		fcc =MATRIXTYPE::Constant(M,N,1.0);
	

	return 1;

}

int pressureFlowFactor(DATATYPE& fxx, DATATYPE& fyy, DATATYPE& fxy , const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai) {
	int iFXC, iFYC;
	DATATYPE c1, r;
	DATATYPE fx, fy;
	
	std::vector<DATATYPE>::const_iterator it = find(FXCgam.begin(), FXCgam.end(), gam);
	if (it != FXCgam.end()) {
		iFXC = std::distance(FXCgam.begin(), it);
		iFYC = 2 + (iFXC - 2) * (-1);
		c1 = FXCc1.at(iFXC);
		r = FXCr.at(iFXC);
		if (gam <= 1)
			fx = 1.0 - c1 * exp(-r * H);
		else
			fx = 1.0 + c1 * pow(H, -r);

		c1 = FXCc1.at(iFYC);
		r = FXCr.at(iFYC);
		if (1.0/gam <= 1)
			fy = 1.0 - c1 * exp(-r * H);
		else
			fy = 1.0 + c1 * pow(H, -r);
	} else
	{
		fx = 1;
		fy = 1;
	}
	DATATYPE cost = cos(theta), sint = sin(theta);
	fxx = cost * cost * fx + sint * sint * fy;
	fyy = cost * cost * fy + sint * sint * fx;
	fxy = sint * cost * (fx - fy);
	return 1;

}

int pressureFlowFactor(MATRIXTYPE& fxx, MATRIXTYPE& fyy, MATRIXTYPE& fxy, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& H, const MATRIXTYPE& Fai) {
	int iFXC, iFYC;
	DATATYPE c1, r,c1y,ry;
	DATATYPE fx, fy, cost,sint;
	std::vector<DATATYPE>::const_iterator it = find(FXCgam.begin(), FXCgam.end(), gam);

	int M = H.rows(), N = H.cols();
	if (it != FXCgam.end()) {
		iFXC = std::distance(FXCgam.begin(), it);
		iFYC = 2 + (iFXC - 2) * (-1);
		c1 = FXCc1.at(iFXC);
		r = FXCr.at(iFXC);
		c1y = FXCc1.at(iFYC);
		ry = FXCr.at(iFYC);
		if (gam <= 1) {
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++) {

					//	gam <= 1 ? ( fx = 1.0 - c1 * exp(-r * H(j, i)), fy = 1.0 + c1y * pow(H(j, i), -ry) ) : (fx = 1.0 + c1 * pow(H(j, i), -r), fy = 1.0 - c1y * exp(-ry * H(j, i)));
					fx = 1.0 - c1 * exp(-r * H(j, i)), fy = 1.0 + c1y * pow(H(j, i), -ry);
					cost = costheta(j, i), sint = sintheta(j, i);
					fxx(j, i) = cost * cost * fx + sint * sint * fy;
					fyy(j, i) = cost * cost * fy + sint * sint * fx;
					fxy(j, i) = sint * cost * (fx - fy);
				}
		}
		else
		{
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++) {

					//	gam <= 1 ? ( fx = 1.0 - c1 * exp(-r * H(j, i)), fy = 1.0 + c1y * pow(H(j, i), -ry) ) : (fx = 1.0 + c1 * pow(H(j, i), -r), fy = 1.0 - c1y * exp(-ry * H(j, i)));
					fx = 1.0 + c1 * pow(H(j, i), -r), fy = 1.0 - c1y * exp(-ry * H(j, i));
					cost = costheta(j, i), sint = sintheta(j, i);
					fxx(j, i) = cost * cost * fx + sint * sint * fy;
					fyy(j, i) = cost * cost * fy + sint * sint * fx;
					fxy(j, i) = sint * cost * (fx - fy);
				}
		}
		

		
	}
	else
	{
		fxx = MATRIXTYPE::Constant(M, N, 1.0);
		fyy = MATRIXTYPE::Constant(M, N, 1.0);
		fxy = MATRIXTYPE::Constant(M, N, 1.0);
	}
	
	return 1;

}
int shearFlowFactor(DATATYPE& fscx, DATATYPE& fscy, const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai) {

	DATATYPE HL, HU;
	DATATYPE FaiHHL, FaiHHU;

	int iFSC0L, iFSC0U, iFSC90L, iFSC90U, iFSC0, iFSC90;
	DATATYPE fsc0L, fsc0U, fsc0, fsc90L, fsc90U, fsc90;

	HL = floor(H);
	HU = ceil(H);

	if (HL > 1)
		FaiHHL = Fai / HL / HL;
	else
		FaiHHL = Fai;

	if (HU > 1)
		FaiHHU = Fai / HU / HU;
	else
		FaiHHU = Fai;

	if (HL != HU) {

		iFSC0L = findShearFactorPosition(gam, FaiHHL, HL);
		iFSC0U = findShearFactorPosition(gam, FaiHHU, HU);
		iFSC90L = findShearFactorPosition(1.0/gam, FaiHHL, HL);
		iFSC90U = findShearFactorPosition(1.0/gam, FaiHHU, HU);

		fsc0L = shearFlowFactorEquation(FSCEqu.at(iFSC0L),FSCm.at(iFSC0L),FSCc2.at(iFSC0L),Fai,HL);
		fsc0U = shearFlowFactorEquation(FSCEqu.at(iFSC0U), FSCm.at(iFSC0U), FSCc2.at(iFSC0U), Fai, HU);
		fsc0 = fsc0L + (H - HL) * (fsc0U - fsc0L) / (HU - HL);
		fsc90L = shearFlowFactorEquation(FSCEqu.at(iFSC90L), FSCm.at(iFSC90L), FSCc2.at(iFSC90L), Fai, HL);
		fsc90U = shearFlowFactorEquation(FSCEqu.at(iFSC90U), FSCm.at(iFSC90U), FSCc2.at(iFSC90U), Fai, HU);
		fsc90 = fsc90L + (H - HL) * (fsc90U - fsc90L) / (HU - HL);
	} 
	else
	{
		iFSC0 = findShearFactorPosition(gam, FaiHHL, H);
		iFSC90 = findShearFactorPosition(1.0/gam, FaiHHL, H);
		fsc0 = shearFlowFactorEquation(FSCEqu.at(iFSC0), FSCm.at(iFSC0), FSCc2.at(iFSC0), Fai, H);
		fsc90 = shearFlowFactorEquation(FSCEqu.at(iFSC90), FSCm.at(iFSC90), FSCc2.at(iFSC90), Fai, H);
	}
	DATATYPE cost = cos(theta), sint = sin(theta);
	fscx = cost* cost* fsc0 + sint * sint * fsc90;
	fscy = sint * cost * (fsc90 - fsc0);
	return 1;
}

int shearFlowFactor(MATRIXTYPE& fscx, MATRIXTYPE& fscy, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& Hin, const MATRIXTYPE& Faiin) {

	DATATYPE HL, HU;
	DATATYPE FaiHHL, FaiHHU;

	int iFSC0L, iFSC0U, iFSC90L, iFSC90U, iFSC0, iFSC90;
	DATATYPE fsc0L, fsc0U, fsc0, fsc90L, fsc90U, fsc90, H, Fai, cost, sint;
	int M = Hin.rows(), N = Hin.cols();

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			H = Hin(j, i);
			Fai = Faiin(j, i);
			HL = int(H);
			HU = H>HL ? HL+1:H;

			FaiHHL=HL > 1 ?  Fai / HL / HL : Fai;
			FaiHHU=HU > 1 ?  Fai / HU / HU : Fai;

			if (HU > HL) {

				iFSC0L = findShearFactorPosition(gam, FaiHHL, HL);
				iFSC0U = findShearFactorPosition(gam, FaiHHU, HU);
				iFSC90L = findShearFactorPosition(1.0 / gam, FaiHHL, HL);
				iFSC90U = findShearFactorPosition(1.0 / gam, FaiHHU, HU);
				//findShearFactorPosition(iFSC0L,iFSC90L,gam, FaiHHL, HL);
				//findShearFactorPosition(iFSC0U, iFSC90U, gam, FaiHHU, HU);
				fsc0L = shearFlowFactorEquation(FSCEqu.at(iFSC0L), FSCm.at(iFSC0L), FSCc2.at(iFSC0L), Fai, HL);
				fsc0U = shearFlowFactorEquation(FSCEqu.at(iFSC0U), FSCm.at(iFSC0U), FSCc2.at(iFSC0U), Fai, HU);
				fsc0 = fsc0L + (H - HL) * (fsc0U - fsc0L) / (HU - HL);
				fsc90L = shearFlowFactorEquation(FSCEqu.at(iFSC90L), FSCm.at(iFSC90L), FSCc2.at(iFSC90L), Fai, HL);
				fsc90U = shearFlowFactorEquation(FSCEqu.at(iFSC90U), FSCm.at(iFSC90U), FSCc2.at(iFSC90U), Fai, HU);
				fsc90 = fsc90L + (H - HL) * (fsc90U - fsc90L) / (HU - HL);
			}
			else
			{
				iFSC0 = findShearFactorPosition(gam, FaiHHL, H);
				iFSC90 = findShearFactorPosition(1.0 / gam, FaiHHL, H);
				//findShearFactorPosition(iFSC0,iFSC90,gam, FaiHHL, H);
				fsc0 = shearFlowFactorEquation(FSCEqu.at(iFSC0), FSCm.at(iFSC0), FSCc2.at(iFSC0), Fai, H);
				fsc90 = shearFlowFactorEquation(FSCEqu.at(iFSC90), FSCm.at(iFSC90), FSCc2.at(iFSC90), Fai, H);
			}
			cost = costheta(j,i), sint = sintheta(j,i);
			fscx(j,i) = cost * cost * fsc0 + sint * sint * fsc90;
			fscy(j,i) = sint * cost * (fsc90 - fsc0);
		}
	
	return 1;
}

int findShearFactorPosition(const DATATYPE& gam, const DATATYPE& FaiHH, const DATATYPE & H) {

	int pos0, pos1;
	
	for (std::vector<DATATYPE>::const_iterator it = FSCgam.begin(); it != FSCgam.end(); ++it) {
		pos1 = std::distance( FSCgam.begin(),it);
		if (*it == gam && FSCH.at(pos1) == H && FaiHH < FSCFaiHH.at(pos1))  return pos1;

		if (*it == gam && FSCH.at(pos1) == 100) pos0 = pos1;
	}
//	if (it != FSCgam.end())
//		return pos1;
//	else
		return pos0;
}

void findShearFactorPosition(int& iFSC0, int& iFSC90, const DATATYPE& gam, const DATATYPE& FaiHH, const DATATYPE& H) {

	int pos0, pos1;

	for (std::vector<DATATYPE>::const_iterator it = FSCgam.begin(); it != FSCgam.end(); ++it) {
		pos1 = std::distance(FSCgam.begin(), it);
		if (*it == gam && FSCH.at(pos1) == H && FaiHH < FSCFaiHH.at(pos1))  iFSC0= pos1;
		if (*it == 1.0/gam && FSCH.at(pos1) == H && FaiHH < FSCFaiHH.at(pos1))  iFSC90 = pos1;
		
		if (*it == gam && FSCH.at(pos1) == 0) iFSC0 = pos1;
		if (*it == 1.0/gam && FSCH.at(pos1) == 0) iFSC90 = pos1;
	}
	//	if (it != FSCgam.end())
	//		return pos1;
	//	else

}