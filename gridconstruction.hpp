#ifndef GRIDCONSTRUCTION
#define GRIDCONSTRUCTION

#include "vectorXYZ.hpp"
#include "vectorIJK.hpp"
 
int MakeRegularGrid(vector <vectorXYZ> *point, vectorIJK *dimension, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax);

vector<int> neighbors(vectorIJK dimension);

#endif
