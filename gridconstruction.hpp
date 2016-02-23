#ifndef GCONST
#define GCONST

#include "vectorXYZ.hpp"
#include "vectorIJK.hpp"
 
int gridfsle2d(vector<vectorXYZ> *itracer, vectorIJK *dimension, vector<int> *neighbor, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax);

vector<int> neighbors(vectorIJK dimension);

#endif
