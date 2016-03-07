#ifndef VELOCITY
#define VELOCITY

#include <ctime>
#include "vectorXYZ.hpp"
#include "vectorIJK.hpp"



// Functions
int LoadVGrid(struct tm rdate, vectorXYZ domainmin, vectorXYZ domainmax, vectorXYZ meanvel, double duration);
int LoadVFlow(struct tm seeddate, int ntime);
void FreeVFlow(unsigned int ntime);

int GetVFlow(double t,vectorXYZ point, vectorXYZ *vint);
int GetVPartialDeriv(double t,vectorXYZ point, vectorIJK dir, vectorXYZ *dvdr);

#endif
