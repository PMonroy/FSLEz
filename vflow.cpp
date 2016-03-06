#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>

using namespace std;

#include "vectorXYZ.hpp"
#include "vectorIJK.hpp"
#include "constants.hpp"

static const int NC_ERR = 2;

//EXTERN VARIABLES
extern int verbose;
extern int tau;

//GLOBAL VARIABLES
int nvlon, nvlat, nvdepth;

vector<double > vlon;
vector<double > vlat;
vector<double > vdepth;
vector<vectorXYZ > vflow;

// FLOW MODEL PARAMETERS
enum vmodels {  MYOCEAN,fi
	      NUMMODELS
};
struct velocitymodel {
string dir;

string slatdim;
string slondim;
string sdepthdim;

string slatvar;
string slonvar;
string sdepthvar;

double latstep;//must be in degrees
double lonstep;//must be in degrees
double depthpositive;//it is used to define the direction positive

string svvar;
string suvar;
double vscalefactor;// factor to get velocity in meters per day
double uscalefactor;// factor to get velocity in meters per day
double vfillvalue;// velocity meters per day
double ufillvalue;// velocity meters per day
};
velocitymodel vmodel[NUMMODELS];

//FUNCTIONS PROTOTYPES
int GetIndices(vectorXYZ point, vectorIJK *index, int vfield);
vectorIJK GetIndices(vectorXYZ point, int vfield);
void locate(double xx[], unsigned int n, double x, unsigned int *j);
void loadparamsmodel(void);

//FUNCTIONS 
void loadparamsmodel(void ){
/**************************************
 * MYOCEAN PARAMETERS
 *************************************/
vmodel[MYOCEAN].dir="/net/argos/data/peps/dflod/data/FSLEz/MODEL_FORECAST/netcdf/";
vmodel[MYOCEAN].slatdim="latitude";
vmodel[MYOCEAN].slondim="longitude";
vmodel[MYOCEAN].sdepthdim="depth";

vmodel[MYOCEAN].slatvar="latitude";
vmodel[MYOCEAN].slonvar="longitude";
vmodel[MYOCEAN].sdepthvar="depth";

vmodel[MYOCEAN].latstep=1.0/12.0;//must be in degrees
vmodel[MYOCEAN].lonstep=1.0/12.0;//must be in degrees

vmodel[MYOCEAN].depthpositive=-1.0;

vmodel[MYOCEAN].svvar="v";
vmodel[MYOCEAN].suvar="u";

vmodel[MYOCEAN].vscalefactor=0.000610370188951492;
vmodel[MYOCEAN].vscalefactor*=secondsday;// factor to get velocity in meters per day

vmodel[MYOCEAN].uscalefactor=0.000610370188951492;
vmodel[MYOCEAN].uscalefactor*=secondsday;// factor to get velocity in meters per day

vmodel[MYOCEAN].vfillvalue= -32767.0;
vmodel[MYOCEAN].vfillvalue*=vmodel[MYOCEAN].vscalefactor;// velocity meters per day

vmodel[MYOCEAN].ufillvalue= -32767.0;
vmodel[MYOCEAN].ufillvalue*=vmodel[MYOCEAN].uscalefactor;// velocity meters per day
}
int loadvgrid(struct tm rdate, int vfield){
  
char ncfile[256];
NcError err(NcError::verbose_nonfatal);

//Loading parameters of the model
loadparamsmodel();

// Open the first Netcdf file
sprintf(ncfile, "%s%04d-%02d-%02d.nc",vmodel[vfield].dir.c_str(), rdate.tm_year,rdate.tm_mon+1,rdate.tm_mday);
NcFile dataFile(ncfile, NcFile::ReadOnly);  

// Check to see if the file was opened.
if(!dataFile.is_valid())
  return NC_ERR;

NcDim *nvlonDim;
if (!(nvlonDim = dataFile.get_dim(vmodel[vfield].slondim.c_str())))
  return NC_ERR;
nvlon = nvlonDim->size();

NcDim *nvlatDim;
if (!(nvlatDim = dataFile.get_dim(vmodel[vfield].slatdim.c_str())))
  return NC_ERR;
nvlat = nvlatDim->size();

NcDim *nvdepthDim;
if (!(nvdepthDim = dataFile.get_dim(vmodel[vfield].sdepthdim.c_str())))
  return NC_ERR;
nvdepth = nvdepthDim->size();

NcVar *vlonVar;
if (!(vlonVar = dataFile.get_var(vmodel[vfield].slonvar.c_str())))
    return NC_ERR;

NcVar *vlatVar;
if (!(vlatVar = dataFile.get_var(vmodel[vfield].slatvar.c_str())))
  return NC_ERR;

NcVar *vdepthVar;
if (!(vdepthVar = dataFile.get_var(vmodel[vfield].sdepthvar.c_str())))
  return NC_ERR;

vlon.resize(nvlon);
vlat.resize(nvlat);
vdepth.resize(nvdepth);

// Get the lat/lon data from the nc file.
if (!vlonVar->get(&vlon[0], nvlon))
  return NC_ERR;

if (!vlatVar->get(&vlat[0], nvlat))
  return NC_ERR;

if (!vdepthVar->get(&vdepth[0], nvdepth))
  return NC_ERR;

//for(int k=0; k<nvdepth; k++){
//vdepth[k]*=vmodel[vfield].depthpositive;
//}
  

return 0;
}
int loadvflow(struct tm seeddate, int ntime, int vfield, vectorXYZ domainmin, vectorXYZ domainmax){

  struct tm *date = {0};
  time_t seedtime, time; // Date in seconds UTC
  char ncfile[256];
  NcVar *uVar, *vVar;
  
  int t;
  
  vector <long int> ubuffer;
  vector <long int> vbuffer;

  vectorIJK indexmin,indexmax;
  vectorXYZ vdomainmin,vdomainmax;
  double avehorvel = 0.2; // meter per seconds

  vdomainmin.x=domainmin.x-(avehorvel*secondsday*abs(tau)*degrees)/(rearth*cos(rads*domainmin.y));
  vdomainmin.y=domainmin.y-(avehorvel*secondsday*abs(tau)*degrees)/rearth ;
  vdomainmin.z=domainmin.z;

  vdomainmax.x=domainmax.x+(avehorvel*secondsday*abs(tau)*degrees)/(rearth*cos(rads*domainmax.y));
  vdomainmax.y=domainmax.y+(avehorvel*secondsday*abs(tau)*degrees)/rearth;
  vdomainmax.z=domainmax.z;

  indexmin=GetIndices(vdomainmin,vfield);
  indexmax=GetIndices(vdomainmax,vfield);
  indexmax+=vectorIJK(1,1,1);

  vdepth.erase(vdepth.begin()+indexmax.k+1,vdepth.end());
  vdepth.erase(vdepth.begin(),vdepth.begin()+indexmin.k);

  vlon.erase(vlon.begin()+indexmax.i+1,vlon.end());
  vlon.erase(vlon.begin(),vlon.begin()+indexmin.i);

  vlat.erase(vlat.begin()+indexmax.j+1,vlat.end());
  vlat.erase(vlat.begin(),vlat.begin()+indexmin.j);

  int nv = vdepth.size()*vlat.size()*vlon.size();

  nvlon = vlon.size();
  nvlat = vlat.size();
  nvdepth = vdepth.size();

  ubuffer.resize(nv);
  vbuffer.resize(nv);
  vflow.reserve(nv*ntime);

  // Get the Velocity data from the nc file
  seedtime = mktime(&seeddate); // convert date to time in seconds    
  for(t=0; t<ntime; t++){
    time = seedtime + t*secondsday;                                  
    date = gmtime(&time);

    sprintf(ncfile, "%s%04d-%02d-%02d.nc",vmodel[vfield].dir.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday);
    NcFile dataFile(ncfile, NcFile::ReadOnly);

    if(verbose==1){ //Reading nc file 
      cout << "Reading nc file: " << ncfile << " ";
    }
    if(!dataFile.is_valid()){ //Check to see if the file was opened.
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }
    if (!(uVar = dataFile.get_var(vmodel[vfield].suvar.c_str()))){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }
    if (!(vVar = dataFile.get_var(vmodel[vfield].svvar.c_str()))){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }

    if (!uVar->set_cur(0, indexmin.k, indexmin.j, indexmin.i)){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }
    if (!vVar->set_cur(0, indexmin.k, indexmin.j, indexmin.i)){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }
     
    if (!uVar->get(&ubuffer[0], 1, nvdepth, nvlat, nvlon)){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }
    if (!vVar->get(&vbuffer[0], 1, nvdepth, nvlat, nvlon)){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }

    for(int q=0; q<nv; q++){       
      vflow.push_back(vectorXYZ(((double) ubuffer[q])*vmodel[vfield].uscalefactor,
				((double) vbuffer[q])*vmodel[vfield].vscalefactor,
				0.0));
    }
    if(verbose==1){//Success file read
      cout << "[Complete]" << endl;
    }
  }

  return 0;
}

int GetIndices(vectorXYZ point, vectorIJK *index, int vfield){

  /* Locate index longitude*/
  index->i = floor((point.x-vlon[0])/vmodel[vfield].lonstep);
  if(index->i < 0 || index->i >= (int)(vlon.size()-1))
      return 1;

  /* Locate index latitude*/
  index->j = floor((point.y-vlat[0])/vmodel[vfield].latstep);
  if(index->j < 0 || index->j >= (int)(vlat.size()-1))
      return 1;

  /* Locate index latitude*/
  double *vvdepth = &vdepth[0]-1;
  unsigned int index_depth;
  locate(vvdepth, vdepth.size(), point.z, &index_depth);
  if(index_depth == 0 || index_depth == vdepth.size())
    return 1;
  else
    index->k = index_depth - 1;

  return 0;
}
vectorIJK GetIndices(vectorXYZ point, int vfield){

  vectorIJK index;

  /* Locate index longitude*/
  index.i = floor((point.x-vlon[0])/vmodel[vfield].lonstep);
  if(index.i < 0)
    index.i=0;
  else if(index.i>=(nvlon-1))
    index.i=(nvlon-1);

  /* Locate index latitude*/
  index.j = floor((point.y-vlat[0])/vmodel[vfield].latstep);
  if(index.j < 0)
    index.j=0;
  else if(index.j>=(nvlat-1))
    index.j=(nvlat-1);

  /* Locate index latitude*/
  double *vvdepth = &vdepth[0]-1;
  unsigned int index_depth;
  locate(vvdepth, nvdepth, point.z, &index_depth);
  if(index_depth == 0)
    index.k=0;
  else
    index.k=index_depth-1;

  return index;
}

int getvflow(double t,vectorXYZ point, vectorXYZ *vint, int vfield){

  vectorXYZ vgrid[16];
  vectorXYZ vcomp[16];
  vectorIJK index;
  unsigned int time;

  double alpha, beta;
  int i,j,k,l;
  int deltai,deltaj,deltak,deltatime;
  unsigned int q;
  unsigned int p, pl,pk,pj;

  time = (unsigned long) t;

  if(GetIndices(point, &index, vfield)==1)
      return 1;

  /* Vectors and points with only one int index*/
  
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++){
    l = time+deltatime;
    pl = l*nvdepth;

    for(deltak=0; deltak<2; deltak++){
      k = index.k + deltak;
      pk = (pl + k) * nvlat; 

      for(deltaj=0; deltaj<2; deltaj++){
	j = index.j + deltaj;
	pj = (pk + j) * nvlon; 

	for(deltai=0; deltai<2; deltai++){
	  i = index.i + deltai;
	  p = pj + i;
	
	  vgrid[q].x = vlon[i];
	  vgrid[q].y = vlat[j];
	  vgrid[q].z = vdepth[k];
	  vcomp[q] = vflow[p];
	  
	  /* COAST CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	  if(vcomp[q].x == vmodel[vfield].ufillvalue || vcomp[q].y == vmodel[vfield].vfillvalue)
	    return 1;
	
	  q++; 
	}
      }
    }
  }

  /* Longitude Interpolation: */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1 - alpha;
  for(q=0; q<16; q+=2){
    vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
    vgrid[q>>1] = vgrid[q]; 
  }
  /* Latitude Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 8; q+=2){
    vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
    vgrid[q>>1] = vgrid[q];
  }
  /* Longitude Interpolation */
  alpha = (vgrid[1].z - point.z)/(vgrid[1].z - vgrid[0].z);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q+=2){
    vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
    vgrid[q>>1] = vgrid[q];
  }
  /* Time Interpolation: */ 
  alpha = ((double) (time + 1)) - t;  
  beta = 1.0 - alpha;   
  vcomp[0] = alpha * vcomp[0] + beta * vcomp[1];
  /* Interpolated V*/
  *vint = vcomp[0];

  return 0;
}
void locate(double xx[], unsigned int n, double x, unsigned int *j){

  /* Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j]
   * and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned
   * to indicate that x is out of range.
   */

  unsigned long ju,jm,jl;
  int ascnd;
  
  jl=0;
  ju=n+1;
  ascnd=(xx[n] > xx[1]);
  while ((ju-jl) > 1) 
    {
      jm=(ju+jl) >> 1;
      if ((x > xx[jm]) == ascnd)
	jl=jm;
      else
	ju=jm;
    }
  if (x == xx[1]) 
    *j=1;
  else if(x == xx[n]) 
    *j=n-1;
  else 
    *j=jl;
}
