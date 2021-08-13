// RotarySeal.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "RotarySeal.h"
#include "Deformation.h"
#include "flowFactor.h"
#include "TridiagSolver.h"

#define RUNTIMETEST 1

constexpr
unsigned int str2int(const char* str, int h = 0)
{
	return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

int main()
{
	ofstream logfile, datafile;
	ofstream pressuefile, deformationfile, initialdeformfile;
	ifstream inputFile("input.txt");



	DATATYPE DP = 6.875 * 25.4 / 1000;
	DATATYPE h = 2e-6;
	DATATYPE Lx = 600E-3;
	DATATYPE Ly = 50E-3;

	DATATYPE sigma = 1e-6;
	DATATYPE alpha = 0.001;
	DATATYPE lamday = 5e-6;
	DATATYPE lamdax = 0.556E-6;
	DATATYPE yb = Ly / 2;
	DATATYPE ps = 1E5;
	DATATYPE pa = 1E5;
	DATATYPE pc = 0;
	DATATYPE Pcmax = 2E5;
	DATATYPE mu = 0.01;
	DATATYPE U = 160 * PI * DP / 60.;

	DATATYPE A1 = 1.4;
	DATATYPE gam = 9.0;



	DATATYPE ESeal = 6.2E6;
	DATATYPE EShaft = 6.2E6;

	DATATYPE Ralpha = 1.75;
	DATATYPE Rbeta = 1.0;
	DATATYPE RR = Ralpha * sigma;
	DATATYPE eita = Rbeta / RR / RR;
	DATATYPE sigmabar = sigma * pow(RR, 1 / 3) * pow(eita, 2 / 3);
	DATATYPE sealpoissonratio = 0.45;
	DATATYPE shaftpoissonratio = 0.33;
	DATATYPE eps;


	DATATYPE pin ;
	DATATYPE pout ;

	DATATYPE ptol = 1E-5;
	DATATYPE dtol = 1E-3;
	//epsilon=mu*U*Lx/((ps-pc)*sigma^2);

	//theta = %pi/4;

	string UnitSystem;
	string LogFileName, DataFileName, PressureFileName,DeformationFileName;
	DATATYPE URP = 0.75;
	DATATYPE URD = 0.75;

	int M = 81;
	int N = 31;
	int PrintInt;
	int MaxPressIter = 10000;
	int MaxDeformIter = 100;
	bool DeformCal;
	
	if (inputFile.is_open())
	{
		std::string line;
		while (getline(inputFile, line)) {
			line.erase(std::remove_if(line.begin(), line.end(), isspace),
				line.end());
			if (line[0] == '#' || line.empty()) continue;
			auto delimiterPos = line.find("=");
			auto name = line.substr(0, delimiterPos);
			auto value = line.substr(delimiterPos + 1);

			switch (str2int(name.data()))
			{
			case str2int("UNIT"): UnitSystem = value; break;
			case str2int("LogFileName"): LogFileName = value; break;
			case str2int("DataFileName"): DataFileName = value; break;
			case str2int("PressureFileName"): PressureFileName = value; break;
			case str2int("DeformationFileName"): DeformationFileName = value; break;
			case str2int("ShaftDiameter"): DP = stod(value); break;
			case str2int("FilmThickness"): h = stod(value); break;

			case str2int("SealingLength"): Lx = stod(value); break;
			case str2int("SealingWidth"): Ly = stod(value); break;
			case str2int("RMSRoughnessofSealSurface"): sigma = stod(value); break;

			case str2int("ElasticParameterofLip"): alpha = stod(value); break;
			case str2int("AutocorrelationLengthinY"): lamday = stod(value); break;
			case str2int("AutocorrelationLengthinX"): lamdax = stod(value); break;

			case str2int("MaximumLipSurfaceDisplacement"): yb = stod(value)*Ly; break;
			case str2int("SealedPressure"): ps = stod(value); break;
			case str2int("AmbientPressure"): pa = stod(value); break;

			case str2int("CavitationPressure"): pc = stod(value); break;
			case str2int("MaximumContactPressure"): Pcmax = stod(value); break;
			case str2int("Viscosity"): mu = stod(value); break;

			case str2int("SlidingSpeed"): U = stod(value); break;
			case str2int("DimensionlessStaticFilmThickness"): A1 = stod(value); break;
			case str2int("AspectRatioofAsperity"): 
				gam = stod(value); 
				gam = gam >= 1.0 ? gam : round(gam * 10) / 9;
				break;

			case str2int("SealYoungsModulus"): ESeal = stod(value); break;
			case str2int("SealPoissonRatio"): sealpoissonratio = stod(value); break;
			case str2int("ShaftYoungsModulus"): EShaft = stod(value); break;
			case str2int("shaftPoissonRatio"): shaftpoissonratio = stod(value); break;
			case str2int("SealInsidePressure"): pin = stod(value); break;

			case str2int("SealOutsidePressure"): pout = stod(value); break;
			case str2int("PressureTolerance"): ptol = stod(value); break;
			case str2int("DisplacementTolerance"): dtol = stod(value); break;

			case str2int("PressureUnderRelaxationFactor"): URP = stod(value); break;
			case str2int("DisplacementUnderRelaxationFactor"): URP = stod(value); break;
			case str2int("NumberofPointinY"): M = stoi(value); break;
			case str2int("NumberofPointinX"): N = stoi(value); break;
			case str2int("PrintInterval"): PrintInt = stoi(value); break;
			case str2int("DeformationCalculation"): 
				DeformCal = (value=="Yes"||value=="YES")?true:false; 
				break;
			case str2int("MaximumPressureIteration"): MaxPressIter = stoi(value); break;
			case str2int("MaximumDeformationIteration"): MaxDeformIter = stoi(value); break;



			}

			std::cout << name << ": " << value << '\n';
		}

	}
	else {
		std::cerr << "Couldn't open config file for reading.\n";
	}

	if (UnitSystem == "MKS") {
		U *= PI * DP / 60;
	//	U = U;
	}
	cout << "Sliding Speed: " << U << "\n";
	//cout << "pin: " << pin;
	logfile.open(LogFileName, ios::out);
	datafile.open(DataFileName, ios::out);

	DATATYPE Lxy = Lx / Ly;
	DATATYPE Lxy2 = Lxy * Lxy;
	DATATYPE H0 = h / sigma;
	DATATYPE dx = 1.0 / (N - 1);
	DATATYPE dy = 1.0 / (M - 1);
	DATATYPE dxy = dx / dy;
	DATATYPE dyx = dy / dx;
	DATATYPE LD = 1.0/Lxy;
	DATATYPE ybbar = yb / Ly;
	DATATYPE psbar = (pa - pc);

	DATATYPE Ebar = 1.0/((1-sealpoissonratio* sealpoissonratio)/ESeal+ (1 - shaftpoissonratio * shaftpoissonratio) / EShaft )  / psbar;
	//ESeal / (1 - sealpoissonratio * sealpoissonratio) / psbar;//
	DATATYPE pex = psbar;

	DATATYPE inL = 1000 / 25.4;
	DATATYPE milL = 1E6 / 25.4;

	DATATYPE c = 0.5 / N;
	DATATYPE d = 0.5 / M;
	DATATYPE deformPara = Lx / sigma / Ebar/PI;
		//6.0*mu*U*Ly*Ly/Ey/sigma/sigma/sigma;

	VECTORTYPE A(M), B(M), C(M), D(M), X(M);
	VECTORTYPE AA(N), BB(N), CC(N), DD(N),XX(N);
	DATATYPE AE1, AW1, AE2, AW2, AN1, AS1, AN2, AS2, SP, AE, AW, ANE, ANW, ASW, ASE, AN, AS;

	
	cout << "******************************Start The Program " <<  " **********************"<<"\n";
	VECTORTYPE x = VECTORTYPE::LinSpaced(N, 0, 1);
	VECTORTYPE y= VECTORTYPE::LinSpaced(M, 0, 1);

	MATRIXTYPE Ik(M , N);
	MATRIXTYPE H=MATRIXTYPE::Zero(M,N);
	MATRIXTYPE P = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE Psc = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE pact = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE Fai = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE Theta = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE cosTheta = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE sinTheta = MATRIXTYPE::Zero(M, N);

	MATRIXTYPE fxx = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE fyy = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE fxy = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE fcc = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE fscx = MATRIXTYPE::Zero(M, N);
	MATRIXTYPE fscy = MATRIXTYPE::Zero(M, N);

	//Initialization
	P= MATRIXTYPE::Constant(M,N,(pex-pc) / psbar);
	P.row(0) = MATRIXTYPE::Constant(1,N,(pin-pc) / psbar);
	P.row(M-1) = MATRIXTYPE::Constant(1,N,(pout-pc) / psbar);

	MATRIXTYPE Pold = P;
	DATATYPE tempang;
	for (int i=0; i<N; i++)
		for (int j = 0; j < M; j++) {
			
			tempang=(y[j] <= ybbar ? PI / 2.0 / ybbar * y[j] : (PI /2.0* (1-y[j])/(1-ybbar)));
			//cout << tempang;
			Psc(j, i) = (Pcmax * sin(tempang) - pc+pa) / psbar;
			//Psc(j, i) = Pcmax/psbar;
		}

	MATRIXTYPE Px = P - Psc;
	MATRIXTYPE P0 = Psc;
	//cout << Px;
	/*
	auto t0 = std::chrono::high_resolution_clock::now();
	MATRIXTYPE Ik2(M*N, N*N);
	influenceMatrix(Ik2, x, y, c, d, LD);
	Map<VECTORTYPE> Pxv(Px.data(), Px.size());
//	Px.reshaped();
	MATRIXTYPE temp1 = Ik2 * Pxv / Ebar;
	Map<MATRIXTYPE> temp2(temp1.data(), M, N);
	auto t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<DATATYPE> elapsed = t1 - t0;

	std::cout << "Old InfluenceMatrix time: " << elapsed.count() << " s\n";
	*/
	MATRIXTYPE temp(M, N);

	cout << "****************************** Initial Deformation " << " **********************" << "\n";
	auto t0 = std::chrono::high_resolution_clock::now();
/*	for (int i =0; i<N; i++)
		for (int j = 0; j < M; j++) {
			influenceMatrix(Ik, x, y, c, d, LD, x[j], y[i]);
			//cout << Ik;
//			cout << Px;
			temp=Ik.cwiseProduct(Px);
			H(j,i) = temp.sum()/Ebar + A1;
		}
*/	
	auto t1=std::chrono::high_resolution_clock::now();
	std::chrono::duration<DATATYPE> elapsed = t1 - t0;

	std::cout << "Initial deformation time: " << elapsed.count() << " s\n";

	t0 = std::chrono::high_resolution_clock::now();
/*	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			;
			//temp = Ik.cwiseProduct(Px);
			H(j, i) = influenceMatrix( x, y, c, d, LD, x[j], y[i],Px) / Ebar + A1;
		}
*/
	t1 = std::chrono::high_resolution_clock::now();
	elapsed = t1 - t0;

	std::cout << "Initial deformation time new: " << elapsed.count() << " s\n";

	t0 = std::chrono::high_resolution_clock::now();
	//cout << Psc<<"\n";
	//DeformationPolar(H, x, y, c, d, LD, deformPara, Px);
	//H = Px / Ebar;
	//H *= deformPara;
	cout << deformPara << " Ebar:  " <<Ebar<<"\n";
	//H += MATRIXTYPE::Constant(M,N,A1);
	
	t1 = std::chrono::high_resolution_clock::now();
	elapsed = t1 - t0;

	std::cout << "Initial deformation time best: " << elapsed.count() << " s\n";

	initialdeformfile.open("initialdeform.txt", ios::out);

	initialdeformfile << "Initial Deformation" << "\n";
	initialdeformfile << H << "\n" << "\n";
	initialdeformfile.close();
	MATRIXTYPE Hold(M,N);	
	MATRIXTYPE H3(M,N);
	MATRIXTYPE dP(M, N), PHold(M, N);
	H = MATRIXTYPE::Constant(M, N, H0);

	for (int t = 0; t <= 0; t++) {

		cout << "****************************Start Case " << t << " **********************"<<"\n";
		eps = mu * U * Lx / psbar / sigma / sigma;
		cout << "eps: " << eps << "\n";

		asperityAngle(Theta, ybbar, eps, alpha, y);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++) {
			
				cosTheta(j, i) = cos(Theta(j, i));
				sinTheta(j, i) = sin(Theta(j, i));
			}

		DATATYPE wd = 1;
		DATATYPE wp = 1;
		std::vector<DATATYPE> wds, wps;

		//cout << Theta;
		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++) 
				pact(j,i) = P(j,i) * psbar + pc;

		for (int k = 0; k < MaxDeformIter; k++) {
			
			Hold = H;
			PHold = P;
			wp = 1;
			int piter = 0;
			

			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					H3(j, i) = H(j,i)*H(j,i)*H(j,i);

			DATATYPE FaiPara = mu * U * lamdax / sigma / sigma;

			cout << "*************Fai Parameter: " << FaiPara << "**************\n";

			for (piter = 0; piter < MaxPressIter; piter++) {
				Pold = P;

				//t0 = std::chrono::high_resolution_clock::now();


				for (int i = 0; i < N; i++)
					for (int j = 0; j < M; j++) {
						//					pact(j, i) = P(j, i) * (ps - pc) + pc;
						Fai(j, i) = FaiPara / (P(j, i) * psbar+pc);
						//	Fai(j+1, i+1) = FaiPara / (P(j+1, i+1) * (ps - pc));
							//fcc(j, i) = densityFlowFactor(gam, Theta(j, i), H(j, i), Fai(j, i));
							//shearFlowFactor(fscx(j,i), fscy(j,i), gam, Theta(j, i), H(j,i), Fai(j,i));
							//pressureFlowFactor(fxx(j, i), fyy(j, i), fxy(j, i), gam, Theta(j, i), H(j, i), Fai(j, i));						
					}
				shearFlowFactor(fscx, fscy, gam, cosTheta, sinTheta, H, Fai);
				densityFlowFactor(fcc, gam, cosTheta, sinTheta, H, Fai);
				pressureFlowFactor(fxx, fyy, fxy, gam, cosTheta, sinTheta, H, Fai);
				/*
				initialdeformfile.open("initialdeform.txt", ios::out);

				initialdeformfile << "Initial Deformation" << "\n";
				initialdeformfile << fscx << "\n" << "\n";
				initialdeformfile.close();
				*/
				//fscx=MATRIXTYPE::Constant(M, N, 1.0);
				//fscy = MATRIXTYPE::Constant(M, N, 1.0);
				//fcc = MATRIXTYPE::Constant(M,N,1.0);
				//fxx = MATRIXTYPE::Constant(M, N, 1.0);
				//fyy = MATRIXTYPE::Constant(M, N, 1.0);
				//fxy = MATRIXTYPE::Constant(M, N, 1.0);
				//t1 = std::chrono::high_resolution_clock::now();
				//elapsed = t1 - t0;

				//std::cout << "Flow factor: " << elapsed.count() << " s\n";
				//t0 = std::chrono::high_resolution_clock::now();

				for (int j = 1; j < N - 1; j++) {
					for (int i = 1; i < M - 1; i++) {

						AE1 = 2. * H3(i, j + 1) * fxx(i, j + 1) * H3(i, j) * fxx(i, j) * dyx / (H3(i, j + 1) * fxx(i, j + 1) + H3(i, j) * fxx(i, j));
						AW1 = 2. * H3(i, j - 1) * fxx(i, j - 1) * H3(i, j) * fxx(i, j) * dyx / (H3(i, j - 1) * fxx(i, j - 1) + H3(i, j) * fxx(i, j));
						AE2 = H3(i, j + 1) * fxy(i, j + 1) * H3(i, j) * fxy(i, j) * Lxy / 2.0 / (H3(i, j + 1) * fxy(i, j + 1) + H3(i, j) * fxy(i, j));
						AW2 = H3(i, j - 1) * fxy(i, j - 1) * H3(i, j) * fxy(i, j) * Lxy / 2.0 / (H3(i, j - 1) * fxy(i, j - 1) + H3(i, j) * fxy(i, j));

						AN1 = Lxy2 * 2.0 * H3(i + 1, j) * fyy((i)+1, j) * H3(i, j) * fyy(i, j) / (H3(i + 1, j) * fyy((i)+1, j) + H3(i, j) * fyy(i, j)) * dxy;
						AS1 = Lxy2 * 2.0 * H3(i - 1, j) * fyy((i)-1, j) * H3(i, j) * fyy(i, j) / (H3(i - 1, j) * fyy((i)-1, j) + H3(i, j) * fyy(i, j)) * dxy;
						AN2 = Lxy / 2.0 * H3(i + 1, j) * fxy((i)+1, j) * H3(i, j) * fxy(i, j) / (H3(i + 1, j) * fxy((i)+1, j) + H3(i, j) * fxy(i, j));
						AS2 = Lxy / 2.0 * H3(i - 1, j) * fxy((i)-1, j) * H3(i, j) * fxy(i, j) / (H3(i - 1, j) * fxy((i)-1, j) + H3(i, j) * fxy(i, j));
						SP = 6.0 * eps * (dy * (2.0 * H(i, j) * (H(i, j + 1) / (H(i, j) + H(i, j + 1)) - H(i, j - 1) / (H(i, j) + H(i, j - 1))) + (fcc(i, j - 1) - fcc(i, j + 1)) / 2 + (fscx(i, j + 1) - fscx(i, j - 1)) / 2) + Lxy * dx * (fscy(i + 1, j) - fscy(i - 1, j)) / 2);
						B[i] = (AE1 + AW1 + AN1 + AS1);
						C[i] = -(AN1 + AE2 - AW2);
						A[i] = -(AS1 + AW2 - AE2);
						AE = AE1 + AN2 - AS2;
						AW = AW1 + AS2 - AN2;
						ANE = AE2 + AN2;
						ASE = -AE2 - AS2;
						ANW = -AW2 - AN2;
						ASW = AW2 + AS2;
						D[i] = AE * P(i, j + 1) + AW * P(i, j - 1) + ANE * P(i + 1, j + 1) + ASE * P(i - 1, j + 1) + ANW * P(i + 1, j - 1) + ASW * P(i - 1, j - 1) - SP;

					}
					D[0] = P(0, j);
					D[M - 1] = P(M - 1, j);
					B[0] = 1;
					B[M - 1] = 1;
					A[0] = 0.0;
					A[M - 1] = 0.0;
					C[0] = 0.0;
					C[M - 1] = 0.0;
					//t1 = std::chrono::high_resolution_clock::now();
					//elapsed = t1 - t0;

					//std::cout << "TDMA assign: " << elapsed.count() << " s\n";
					//t0 = std::chrono::high_resolution_clock::now();
					//cout << A<<"\n"<<B<<"\n"<<C<<"\n";
					TDMASolver(X, A, B, C, D);

					//t1 = std::chrono::high_resolution_clock::now();
					// elapsed = t1 - t0;

					//std::cout << "TDMA solver: " << elapsed.count() << " s\n";
					//cout << X;

					for (int i = 1; i < M - 1; i++)
					{
						P(i, j) = (1 - URP) * P(i, j) + URP * X[i];
						//		P(i+1, j) = (1 - URP) * P(i+1, j) + URP * X[i+1];
					//	P(i, j) = P(i, j) < P0(i, j) ? P0(i, j) : P(i, j);
					}

				}
			/*
				dP = P - Pold;
				wp = abs(dP.sum() / P.sum());
				//logfile << "       ****Pressure iteration number: " << piter << "; Iteration tolerance: " << wp << "\n";
				if (piter % PrintInt == 0)
					cout << "\n" << "       ****Pressure iteration number: " << piter << "; Iteration tolerance: " << wp << "\n";
				//wps.push_back(wp);
				if (wp < ptol)	break;
			}

			for (piter = 0; piter < 10000; piter++) {
				Pold = P;
				
			*/
				/*
				for (int i = 0; i < N; i ++)
					for (int j = 0; j < M ; j ++) {
						//					pact(j, i) = P(j, i) * (ps - pc) + pc;
						Fai(j, i) = FaiPara / (P(j, i) * psbar+pc);
	
						//fcc(j, i) = densityFlowFactor(gam, Theta(j, i), H(j, i), Fai(j, i));
					//	shearFlowFactor(fscx(j, i), fscy(j, i), gam, Theta(j, i), H(j, i), Fai(j, i));
					//	pressureFlowFactor(fxx(j, i), fyy(j, i), fxy(j, i), gam, Theta(j, i), H(j, i), Fai(j, i));						
					}
				shearFlowFactor(fscx, fscy, gam, cosTheta, sinTheta, H, Fai);
				densityFlowFactor(fcc, gam, cosTheta, sinTheta, H, Fai);
				pressureFlowFactor(fxx, fyy, fxy, gam, cosTheta, sinTheta, H, Fai);
			*/	
				//fscx = MATRIXTYPE::Constant(M, N, 1.0);
				//fscy = MATRIXTYPE::Constant(M, N, 1.0);
				//fcc = MATRIXTYPE::Constant(M, N, 1.0);
				//fxx = MATRIXTYPE::Constant(M, N, 1.0);
				//fyy = MATRIXTYPE::Constant(M, N, 1.0);
				//fxy = MATRIXTYPE::Constant(M, N, 1.0);
				for (int i = 1; i < M-1; i++) {
					for (int j = 1; j < N - 1; j++) {

						AE1 = 2. * H3(i,j+1) * fxx(i, j + 1)*H3(i,j) * fxx(i, j)*dyx / (H3(i,j) * fxx(i, j+1) + H3(i,j) * fxx(i, j));
						AW1 = 2. * H3(i,j - 1) * fxx(i, j-1)*H3(i,j) * fxx(i, j)*dyx / (H3(i,j - 1) * fxx(i, j-1) + H3(i,j) * fxx(i, j));
						AE2 =  H3(i,j+1) * fxy(i, j+1)*H3(i,j) * fxy(i, j)*Lxy / 2.0 / (H3(i,j+1) * fxy(i, j+1) + H3(i,j) * fxy(i, j));
						AW2 =  H3(i,j - 1) * fxy(i, j-1)*H3(i,j) * fxy(i, j)*Lxy / 2.0 / (H3(i,j - 1) * fxy(i, j-1) + H3(i,j) * fxy(i, j));

						AN1 = Lxy2*2.0 * H3(i + 1, j) * fyy(i + 1, j)*H3(i,j) * fyy(i, j) / (H3(i + 1, j) * fyy(i + 1, j) + H3(i,j) * fyy(i, j))*dxy;
						AS1 = Lxy2*2.0 * H3(i - 1, j) * fyy(i - 1, j)*H3(i,j) * fyy(i, j) / (H3(i - 1, j) * fyy(i - 1, j) + H3(i,j) * fyy(i, j))*dxy;
						AN2 = Lxy / 2.  * H3(i + 1, j) * fxy(i + 1, j)*H3(i,j) * fxy(i, j) / (H3(i + 1, j) * fxy(i + 1, j) + H3(i,j) * fxy(i, j));
						AS2 = Lxy/ 2.  * H3(i - 1, j)* fxy(i - 1, j)*H3(i,j) * fxy(i, j) / (H3(i - 1, j) * fxy(i - 1, j) + H3(i,j) * fxy(i, j));
						SP = 6. * eps*(dy*(2. * H(i, j)*(H(i, j+1) / (H(i, j) + H(i, j+1)) - H(i, j-1) / (H(i, j) + H(i, j-1))) + (fcc(i, j-1) - fcc(i, j+1)) / 2 + (fscx(i, j+1) - fscx(i, j-1)) / 2) + Lxy*dx*(fscy(i + 1, j) - fscy(i - 1, j)) / 2);
						AA[j] = AE1 + AW1 + AN1 + AS1;
						BB[j] = AE1 + AN2 - AS2;
						CC[j] = AW1 + AS2 - AN2;

						AN = AN1 + AE2 - AW2;
						AS = AS1 + AW2 - AE2;
						ANE = AE2 + AN2;
						ASE = -AE2 - AS2;
						ANW = -AW2 - AN2;
						ASW = AW2 + AS2;
						DD[j] = AN*P(i + 1, j) + AS*P(i - 1, j) + ANE*P(i + 1, j+1) + ASE*P(i - 1, j+1) + ANW*P(i + 1, j-1) + ASW*P(i - 1, j-1) - SP;

						
					}
					//cout << AA << "\n" << BB << "\n" << CC << "\n"<<DD<<"\n";
					//t0 = std::chrono::high_resolution_clock::now();

					CTDMASolver(XX,AA, BB, CC, DD);

					//t1 = std::chrono::high_resolution_clock::now();
					//elapsed = t1 - t0;

					//std::cout << "CTDMA solver: " << elapsed.count() << " s\n";
					//cout << XX<<"\n";
					for (int j = 0; j < N ; j++)
					{
						P(i, j) = (1 - URP) * P(i, j) + URP * XX[j];
					//	P(i, j+1) = (1 - URP) * P(i, j+1) + URP * XX[j+1];
						P(i, j) = P(i, j) < P0(i, j) ? P0(i, j) : P(i, j);
					}
						
				}
				dP = P - Pold;
				wp = abs(dP.sum() / P.sum());
				//logfile << "       ****Pressure iteration number: " << piter << "; Iteration tolerance: " << wp << "\n";
				if (piter%PrintInt==0)
					cout << "\n" << "       ****Pressure iteration number: " << piter << "; Iteration tolerance: " << wp << "\n";
				//wps.push_back(wp);
				if (wp < ptol)	break;

			}

			if (DeformCal) {
				Px = P-Psc;
				//cout << Psc;
				//Map<VECTORTYPE> Pxv(Px.data(), Px.size());

				//temp = Ik * Pxv / Ebar;
				//Map<MATRIXTYPE> temp2(temp.data(), M, N);
				/*for (int i = 0; i < N; i++)
					for (int j = 0; j < M; j++) {
						influenceMatrix(Ik, x, y, c, d, LD, x[j], y[i]);
						temp=Ik.cwiseProduct(Px);
						H(j, i) =  A1+(1.0 - URD) * (H(j, i) - A1) + URD* temp.sum()/Ebar;
					}
					*/
				/*for (int j = 0; j < N; j++)
					for (int i = 0; i < M; i++)
					{
						//temp2(i, j) = 1.9 * Px(i, j) * dx * (1.0 - pissonratio * pissonratio) / Ey * (ps - pc);
						//H(i, j) = A1 + (1.0 - URD) * (H(i, j) - A1) + URD * temp2(i, j);
						Px(j,i) = Px(j, i) < 0 ? 0.0 : Px(j, i);
					}
					*/
				DeformationPolar(H, x, y, c, d, LD, deformPara,Px);
				//H += Hold;
				//H += MATRIXTYPE::Constant(M, N, URD * A1);
				//H *= URD;
				//H += Hold;
				
					
				//H = Px / Ebar*deformPara*dy;
				H *= URD;
				//H += (1.0 - URD) * (Hold);
				//H += MATRIXTYPE::Constant(M, N, A1);
				H += (1.0 - URD) * (Hold)+MATRIXTYPE::Constant(M, N, URD * A1);
				//H= MATRIXTYPE::Constant(M, N, 1.74);
				/*for (int j=0; j<N; j++)
					for (int i = 0; i < M; i++)
					{
						//temp2(i, j) = 1.9 * Px(i, j) * dx * (1.0 - pissonratio * pissonratio) / Ey * (ps - pc);
						//H(i, j) = A1 + (1.0 - URD) * (H(i, j) - A1) + URD * temp2(i, j);
						H(i, j) = H(i, j) < 1 ? 1.0 : H(i, j);
					}
				*/
				MATRIXTYPE dH = H - Hold;
				wd = abs(dH.sum() / H.sum());
				wds.push_back(wd);

			}
			else
				wd = dtol / 10;
			
			pressuefile.open(PressureFileName, ios::out);
			deformationfile.open(DeformationFileName, ios::out);
			logfile <<"\n"<<"**********Deformation iteration number: " << k << "; Iteration tolerance: " << wd << "\n" << "\n";
			cout << "\n" << "**********Deformation iteration number: " << k << "; Iteration tolerance: " << wd << "\n" << "\n";
			pressuefile << "Pressure Results" << "\n";
			pressuefile << P << "\n" << "\n";
			deformationfile << "Defromation Results" << "\n";
			deformationfile << H << "\n" << "\n";
			pressuefile.close();
			deformationfile.close();
			if (wd < dtol) break;
		}

		for (int i = 0; i < wds.size(); i++)
			logfile << wds[i]<<"\n";

		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++) {
				pact(j, i) = P(j, i) * psbar + pc;
				Fai(j, i) = mu * U * lamdax / (pact(j, i) - pc) / sigma / sigma;
				fcc(j, i) = densityFlowFactor(gam, Theta(j, i), H(j, i), Fai(j, i));
				shearFlowFactor(fscx(j, i), fscy(j, i), gam, Theta(j, i), H(j, i), Fai(j, i));
				pressureFlowFactor(fxx(j, i), fyy(j, i), fxy(j, i), gam, Theta(j, i), H(j, i), Fai(j, i));

		}



		DATATYPE workLoad = dx * dy * (P.sum()*psbar+pc)*Lx*Ly;
		VECTORTYPE DPDY(N), DPDX(N);
		int j = 0;
		for (int i = 0; i < N; i++) {

			//DPDY[i] = (P(M - 1, i) - P(M - 2, i)) / dy;

			DPDY[i] =  -(-25 * P(M - 1, i) + 48 * P(M - 2, i) - 36 * P(M - 3, i) + 16 * P(M - 4, i) - 3 * P(M - 5, i)) / 12.0 / dy;
			DPDY[i] = (-25 * P(0, i) + 48 * P(1, i) - 36 * P(2, i) + 16 * P(3, i) - 3 * P(4, i)) / 12.0 / dy;

			if (i >= 2 && i <= N - 3)
			{
				//DPDX[i] = (P(M - 1, i) - P(M - 1, i - 1)) / dx;
				DPDX[i] = (-P(M - 1, i + 2) + 8 * P(M - 1, i + 1) - 8 * P(M - 1, i - 1) + P(M - 1, i - 2)) / 12 / dx;

			}

			//DPDX(0) = (P(M - 1, 0) - P(M - 1, N - 1)) / dx;
			DPDX(0) =  (-P(M - 1, 2) + 8 * P(M - 1, 1) - 8 * P(M - 1, N - 2) + P(M - 1, N - 3)) / 12 / dx;
			DPDX(N - 1) = DPDX(0);
			DPDX(1) = (-P(M - 1, 3) + 8 * P(M - 1, 2) - 8 * P(M - 1, 0) + P(M - 1, N - 1)) / 12 / dx;
			DPDX(N - 2) =  (-P(M - 1, 0) + 8 * P(M - 1, N - 1) - 8 * P(M - 1, N - 3) + P(M - 1, N - 4)) / 12 / dx;
		}
		VECTORTYPE qy(N);
		for (int i=0; i<N; i++)
			qy[i]= -pow(H(j,i),3)*(fxy(j,i)*DPDX[i]+Lxy*fyy(j,i)*DPDY[i]) + 6.0 * eps * fscy(j,i);

		DATATYPE Qflow = qy(0) ;//dx*qy(0)/2.0+dx*qy(qy.size()-1)/2;

		for (int i = 1; i < N; i++)
		{
			
			Qflow += (i%2==0?2*qy(i):4*qy(i));
		}

		Qflow *= dx / 3;
		Qflow *= psbar * sigma * sigma * sigma / 12.0 / mu;

		

		datafile << "**************************Case " << t <<"************************"<< "\n" << "\n";
		datafile << "******* Dimensionless value*********** "<<"\n" << "\n";
		datafile << "Speed" << "\t" <<"Work Load"<<"\t" <<"Flow Rate" <<"\n";

		datafile << eps << "\t" << workLoad/psbar/Lx/Ly << "\t" << Qflow*12*mu/psbar/sigma/sigma/sigma << "\n"<<"\n";

		datafile << "******* Actual value*********** " << "\n" << "\n";
		datafile << "Speed (m/s)" << "\t" << "Work Load (N)" << "\t" << "Flow Rate (m^3/s)" << "\n";

		datafile << U << "\t" << workLoad  << "\t" << Qflow  << "\n" << "\n";

		cout << "**************************Case " << t << "************************" << "\n" << "\n";
		cout << "******* Dimensionless value*********** " << "\n" << "\n";
		cout << "Speed" << "\t" << "Work Load" << "\t" << "Flow Rate" << "\n";

		cout << eps << "\t" << workLoad / psbar / Lx / Ly << "\t" << Qflow * 12 * mu / psbar / sigma / sigma / sigma << "\n" << "\n";

		cout << "******* Actual value*********** " << "\n" << "\n";
		cout << "Speed (m/s)" << "\t" << "Work Load (N)" << "\t" << "Flow Rate (m^3/s)" << "\n";

		cout << U << "\t" << workLoad << "\t" << Qflow << "\n" << "\n";


	}

	logfile.close();
	datafile.close();
	system("pause");
	return 0;
	
}


