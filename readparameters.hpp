#ifndef DATE
#define DATE

#include <ctime>
#include "vectorXYZ.hpp"

extern int verbose;
extern int vfield;
extern vectorXYZ domainmin;
extern vectorXYZ intergrid;
extern vectorXYZ domainmax;
extern struct tm seeddate;
extern double  intstep;
extern int tau;
extern double deltamax;

int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams);
int GetfileParameters(string nfileparameters);

#endif
