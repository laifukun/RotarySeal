#pragma once
#include "RotarySeal.h"

int influenceMatrix(MATRIXTYPE& Ik, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD);
int influenceMatrix(MATRIXTYPE& Ik, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& xi, const DATATYPE& yi);
DATATYPE influenceMatrix(const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& xi, const DATATYPE& yi, const MATRIXTYPE& P);
int asperityAngle(MATRIXTYPE& theta, const DATATYPE& ybbar, const DATATYPE& eps, const DATATYPE& alpha, const VECTORTYPE& y);
int Deformation(MATRIXTYPE& H, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& deformPara, const MATRIXTYPE& P);
int DeformationPolar(MATRIXTYPE& H, const VECTORTYPE& x, const VECTORTYPE& y, const DATATYPE& c, const DATATYPE& d, const DATATYPE& LD, const DATATYPE& deformPara, const MATRIXTYPE& P);