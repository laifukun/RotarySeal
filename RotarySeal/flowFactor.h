#pragma once
#include "RotarySeal.h"

const std::vector<DATATYPE> FXCgam = { 1.0/9, 1.0/3,  1,     3,     9 };
const std::vector<DATATYPE> FXCc1 = {  1.714, 1.284,  0.611, 0.186, 1.078 };
const std::vector<DATATYPE> FXCr = {   0.430, 0.440,  0.540, 1.330, 1.640};

const std::vector<DATATYPE> FCgam = { 1.0/9, 1.0/3, 1,     3,     9 };
const std::vector<DATATYPE> FCa2 = {  0.926, 0.686, 0.657, 0.730, 0.927 };
const std::vector<DATATYPE> FCb2 = {  0.005, 0.019, 0.081, 0.503, 3.704 };

const std::vector<DATATYPE> FSCgam =   { 1.0/9,   1.0/9,   1.0/9, 1.0/9, 1.0/9,  1.0/9,  1.0/9,  1.0/9,  1.0/9,  1.0/9, 1.0/3,   1.0/3,   1.0/3, 1.0/3,  1.0/3,  1.0/3,  1.0/3,  1.0/3, 1.0/3,   1.0/3, 1.0,      1,        1,     1,      1,      1,      1,      1,     1,       1,     3,     3,      3,      3,      3,      3,      3,      3,      3,     9,       9,       9,      9,      9,      9,      9,      9,      9,      9 };
const std::vector<DATATYPE> FSCH =     { 0.0,     1,       1,     1,     2,      3,      4,      5,      6,      100,   0,       1,       1,     1,      2,      3,      4,      5,      6,      100,   0,        1,        1,     1,      2,      3,      4,      5,     6,       100,   0,     1,      1,      2,      3,      4,      5,      6,      100,   0,       1,       1,      1,      2,      3,      4,      5,      6,      100 };
const std::vector<DATATYPE> FSCm =     { 0,       0,       0.899, 2.121, 0,      0.4715, 0.4164, 0.2053, 0.1227, 2.121, 0,       0,       0.968, 1.706,  0,      0.5051, 0.2291, 0.1217, 0.0782,  1.706, 0,        0,        1.007, 1.4,    0,      0.1295, 0.0593, 0.0363, 0.0239, 1.400, 0,     0,      1.554,  0,      0.0341, 0.0110, 0.0084, 0.0046, 1.554, 0,       0,       3.171,  1.443,  0,      0.0059, 0.0029, 0.0015, 0.0010, 1.443 };
const std::vector<DATATYPE> FSCc2 =    { 1.247,   1.247,   0.802, 0.896, 1.10,   0.95,   0.74,   0.55,   0.45,   0.896, 1.172,   1.172,   0.853, 0.99,   0.98,   0.80,   0.60,   0.47,  0.39,    0.990, 0.867,    0.867,    1.153, 1.334,  0.68,   0.51,   0.38,   0.31,  0.26,    1.334, 0.478, 0.478,  1.936,  0.37,   0.26,   0.18,   0.15,   0.12,   1.936, 0.331,   0.331,   3.021,  4.762,  0.21,   0.11,   0.08,   0.06,   0.05,   4.762 };
const std::vector<DATATYPE> FSCFaiHH = { 7.03E-9, 7.03E-9, 0.03,  100,   0.01,   0.1,    0.4,    0.8,    1.0,    100,   7.61E-8, 7.61E-8, 0.07,  100,    0.02,   0.4,    1.0,    1.3,   1.7,     100,   1.26E-11, 1.26E-11, 0.12,  100,    0.1,    0.8,    1.5,    2.1,   2.7,     100,   0.1,   0.1,    100,    0.5,    2.3,    3.2,    4.5,    5.6,    100,   2.18E-9, 2.18E-9, 1.15,   100,    0.1,    6,      12,     21,     27,     100 };
const std::vector<int> FSCEqu =        { 1,       1,       2,     3,     1,      1,      1,      1,      1,      3,     1,       1,       2,     3,      1,      1,      1,      1,     1,       3,     1,        1,        2,     3,      1,      1,      1,      1,     1,       3,     1,     1,      3,      1,      1,      1,      1,      1,      3,     1,       1,       2,      3,      1,      1,      1,      1,      1,      3 };


DATATYPE shearFlowFactorEquation(const int& Equ, const DATATYPE& m, const DATATYPE& c2, const DATATYPE& Fai, const DATATYPE& H);
DATATYPE densityFlowFactor(const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai);
int densityFlowFactor(MATRIXTYPE& fcc, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& H, const MATRIXTYPE& Fai);
int pressureFlowFactor(DATATYPE& fxx, DATATYPE& fyy, DATATYPE& fxy, const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai);
int pressureFlowFactor(MATRIXTYPE& fxx, MATRIXTYPE& fyy, MATRIXTYPE& fxy, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& H, const MATRIXTYPE& Fai);
int shearFlowFactor(DATATYPE& fscx, DATATYPE& fscy, const DATATYPE& gam, const DATATYPE& theta, const DATATYPE& H, const DATATYPE& Fai);
int shearFlowFactor(MATRIXTYPE& fscx, MATRIXTYPE& fscy, const DATATYPE& gam, const MATRIXTYPE& costheta, const MATRIXTYPE& sintheta, const MATRIXTYPE& Hin, const MATRIXTYPE& Faiin);
//int shearFlowFactor(MATRIXTYPE& fscx, MATRIXTYPE& fscy, const DATATYPE& gam, const MATRIXTYPE& theta, const MATRIXTYPE& H, const MATRIXTYPE& Fai);
int findShearFactorPosition(const DATATYPE& gam, const DATATYPE& FaiHH, const DATATYPE& H);
void findShearFactorPosition(int& iFSC0, int& iFSC90, const DATATYPE& gam, const DATATYPE& FaiHH, const DATATYPE& H);

