#ifndef VELOCITY
#define VELOCITY

#include <ctime>
#include "vectorXYZ.hpp"

int loadvgrid(struct tm rdate,int vfield);
int loadvflow(struct tm seeddate, int ntime, int vfield);

int getvflow(double t,vectorXYZ point, vectorXYZ *vint, int vfield);

#endif
