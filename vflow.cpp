#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

#include "vflow.hpp"
#include "constants.hpp"
#include "rparameters.hpp"
static const int NC_ERR = 2;

//GLOBAL VARIABLES
unsigned int nvlon, nvlat, nvdepth;
unsigned int nv;

vector <double> vlon;
vector <double> vlat;
vector <double> vdepth;

vectorIJK indexmin,indexmax;

vectorXYZ**** vflow;

struct vflowParameters {

  const string dir;

  const string slatdim;
  const string slondim;
  const string sdepthdim;
  
  const string slatvar;
  const string slonvar;
  const string sdepthvar;
  
  const double latstep;//must be in degrees
  const double lonstep;//must be in degrees
  const double depthpositive;//it is used to define the direction positive
  
  const string suvar;
  const string svvar;
  const double uscalefactor;// factor to get velocity in meters per day
  const double vscalefactor;// factor to get velocity in meters per day
  const double ufillvalue;// velocity meters per day
  const double vfillvalue;// velocity meters per day


  // The constructor:
  vflowParameters(const string & vflowParamsFileName)
    :dir(getStringParam(vflowParamsFileName, "dir"))
    ,slatdim(getStringParam(vflowParamsFileName, "slatdim"))
    ,slondim(getStringParam(vflowParamsFileName, "slondim"))
    ,sdepthdim(getStringParam(vflowParamsFileName, "sdepthdim"))
    ,slatvar(getStringParam(vflowParamsFileName, "slatvar"))
    ,slonvar(getStringParam(vflowParamsFileName, "slonvar"))
    ,sdepthvar(getStringParam(vflowParamsFileName, "sdepthvar"))
    ,latstep(getDoubleParam(vflowParamsFileName, "latstep"))
    ,lonstep(getDoubleParam(vflowParamsFileName, "lonstep"))
    ,depthpositive(getDoubleParam(vflowParamsFileName, "depthpositive"))
    ,suvar(getStringParam(vflowParamsFileName, "suvar"))
    ,svvar(getStringParam(vflowParamsFileName, "svvar"))
    ,uscalefactor(getDoubleParam(vflowParamsFileName, "uscalefactor"))
    ,vscalefactor(getDoubleParam(vflowParamsFileName, "vscalefactor"))
    ,ufillvalue(getDoubleParam(vflowParamsFileName, "ufillvalue"))
    ,vfillvalue(getDoubleParam(vflowParamsFileName, "vfillvalue"))
{}

};
const vflowParameters vflowparams("vflow.conf");

//FUNCTIONS PROTOTYPES
vectorIJK GetIndices(vectorXYZ point);

//INLINE FUNCTIONS 
inline int LocateIndex(double x, double start, double step){
  return floor((x-start)/step);
}
inline int LocateIndex(double x, const vector <double> &xx){
  /* Given a vector xx[0...n-1], and given a value x, returns a value i such that x is between xx[i]
   * and xx[i+1]. xx must be monotonic, either increasing or decreasing. -1 or n is returned
   * to indicate that x is out of range.
   */

  if (x < xx.front()) 
    return -1;
  else if(x >= xx.back()) 
    return xx.size();
  else { 		
    unsigned long jl,jm,ju;
    jl=0;
    ju=xx.size()+1;
    int ascnd = (xx.back()>xx.front());
    while((ju-jl)>1) {
      jm =(ju+jl)>>1;
      if((x>=xx.at(jm-1))==ascnd) jl=jm;
      else ju=jm;
    }
    return int(jl-1);
  }
}

//FUNCTIONS DECLARATIONS
int LoadVGrid(struct tm rdate, vectorXYZ domainmin, vectorXYZ domainmax, vectorXYZ meanvel, double duration){
 
  NcError err(NcError::verbose_nonfatal);
  
  // Open the first Netcdf file
  char ncfile[256];
  sprintf(ncfile, "%s%04d-%02d-%02d.nc",vflowparams.dir.c_str(), rdate.tm_year,rdate.tm_mon+1,rdate.tm_mday);
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);  
  if(!dataFile.is_valid()){
    cout <<ncfile <<endl;
    cout << "Error opening file:"<< ncfile <<endl;
    return NC_ERR;
  }

 // Get grid variable dimensions:
  NcDim *nvlonDim;
  NcDim *nvlatDim;
  NcDim *nvdepthDim;

  if(!(nvlonDim = dataFile.get_dim(vflowparams.slondim.c_str()))
    || !(nvlatDim = dataFile.get_dim(vflowparams.slatdim.c_str()))
    || !(nvdepthDim = dataFile.get_dim(vflowparams.sdepthdim.c_str()))){
    cout << "Error getting grid variable dimensions"<<endl;
    return NC_ERR;
  }

  NcVar *vlonVar;
  NcVar *vlatVar;
  NcVar *vdepthVar;
  
  nvlon = nvlonDim->size();
  nvlat = nvlatDim->size();
  nvdepth = nvdepthDim->size();
  vlon.resize(nvlon);
  vlat.resize(nvlat);
  vdepth.resize(nvdepth);

  if (!(vlonVar = dataFile.get_var(vflowparams.slonvar.c_str()))
      || !vlonVar->get(&vlon[0], nvlon)){
    cout << "Error reading longitude variable"<<endl;
    return NC_ERR;
  }

  if (!(vlatVar = dataFile.get_var(vflowparams.slatvar.c_str()))
      || !vlatVar->get(&vlat[0], nvlat)){
    cout << "Error reading latitude variable"<<endl;
    return NC_ERR;
  }

  if (!(vdepthVar = dataFile.get_var(vflowparams.sdepthvar.c_str()))
      || !vdepthVar->get(&vdepth[0], nvdepth)){
    cout << "Error reading depth variable"<<endl;
    return NC_ERR;
  }
  
#ifdef DEBUG
    cout << "Ncfile grid:" << endl;
    cout << " Longitude = "<< vlon.front()   <<","<< vlon.back()  <<" ("<< 0  <<","<< vlon.size()   <<")"<<endl;
    cout << " Latitude = " << vlat.front()   <<","<< vlat.back()  <<" ("<< 0  <<","<< vlat.size()   <<")"<<endl;
    cout << " Depth = "    << vdepth.front() <<","<< vdepth.back()<<" ("<< 0  <<","<< vdepth.size() <<")"<<endl;
#endif

  vectorXYZ vdomainmin,vdomainmax;
  //vectorXYZ avevelocity=vectorXYZ(0.1,0.1,0.0); // meter per seconds
  double scalefactor=(secondsday*duration*degrees)/rearth;

  vdomainmin=domainmin-meanvel*vectorXYZ(scalefactor/cos(rads*domainmin.y),scalefactor,1);
  vdomainmax=domainmax+meanvel*vectorXYZ(scalefactor/cos(rads*domainmax.y),scalefactor,1);

  indexmin=GetIndices(vdomainmin);
  indexmax=GetIndices(vdomainmax);
  indexmax+=vectorIJK(2,2,2);

  vdepth.erase(vdepth.begin()+indexmax.k,vdepth.end());
  vdepth.erase(vdepth.begin(),vdepth.begin()+indexmin.k);

  vlon.erase(vlon.begin()+indexmax.i,vlon.end());
  vlon.erase(vlon.begin(),vlon.begin()+indexmin.i);

  vlat.erase(vlat.begin()+indexmax.j,vlat.end());
  vlat.erase(vlat.begin(),vlat.begin()+indexmin.j);

  nvlon = vlon.size();
  nvlat = vlat.size();
  nvdepth = vdepth.size();
  nv = nvdepth*nvlat*nvlon;
  
#ifdef DEBUG
    cout << "Velocity grid:" << endl;
    cout << " Longitude = "<< vlon.front()   <<","<< vlon.back()  <<" ("<< 0  <<","<< vlon.size()   <<")"<<endl;
    cout << " Latitude = " << vlat.front()   <<","<< vlat.back()  <<" ("<< 0  <<","<< vlat.size()   <<")"<<endl;
    cout << " Depth = "    << vdepth.front() <<","<< vdepth.back()<<" ("<< 0  <<","<< vdepth.size() <<")"<<endl;
#endif

  return 0;
}
int LoadVFlow(struct tm seeddate, int ntime){

  struct tm *date = {0};
  time_t seedtime, time; // Date in seconds UTC
  char ncfile[256];
  NcVar *uVar, *vVar;
  
  int t;
  
  vector <long int> vflowx_buffer(nv);
  vector <long int> vflowy_buffer(nv);
  //vector <double*> vflowz_buffer(nv);

  vflow = new vectorXYZ ***[ntime];

  // Get the Velocity data from the nc file
  seedtime = mktime(&seeddate); // convert date to time in seconds    
  for(t=0; t<ntime; t++){
    time = seedtime + t*secondsday;                                  
    date = gmtime(&time);

    sprintf(ncfile, "%s%04d-%02d-%02d.nc",vflowparams.dir.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday);
    NcFile dataFile(ncfile, NcFile::ReadOnly);

#ifdef DEBUG
      //Reading nc file 
      cout << "Reading nc file: " << ncfile << " ";
#endif

    if(!dataFile.is_valid() 
       || !(uVar = dataFile.get_var(vflowparams.suvar.c_str()))
       || !(vVar = dataFile.get_var(vflowparams.svvar.c_str()))
       || !uVar->set_cur(0, indexmin.k, indexmin.j, indexmin.i)
       || !vVar->set_cur(0, indexmin.k, indexmin.j, indexmin.i)
       || !uVar->get(&vflowx_buffer[0], 1, nvdepth, nvlat, nvlon)
       || !vVar->get(&vflowy_buffer[0], 1, nvdepth, nvlat, nvlon)){
      cout << "[Fail]" <<endl;
      return NC_ERR;
    }

    unsigned int fk,fj,q;

    vflow[t] = new vectorXYZ **[nvdepth];
    for(unsigned int k=0; k<nvdepth; k++){
      fk=k*nvlat;
      vflow[t][k] = new vectorXYZ *[nvlat];      
      for(unsigned int j=0; j<nvlat; j++){
	fj=(fk+j)*nvlon;
	vflow[t][k][j] = new vectorXYZ [nvlon];
	for(unsigned int i=0; i<nvlon; i++){
	  q=fj+i;
	  vflow[t][k][j][i].x= double(vflowx_buffer[q])*vflowparams.uscalefactor;
	  vflow[t][k][j][i].y= double(vflowy_buffer[q])*vflowparams.vscalefactor;
	  vflow[t][k][j][i].z= 0.0;
	}
      }
    }

#ifdef DEBUG
    //Success file read
    cout << "[Complete]" << endl;
#endif
  }

  return 0;
}

void FreeVFlow(unsigned int ntime){

  unsigned int t, k, j;

  for(t=0; t<ntime; t++) {
    for(k=0; k<nvdepth; k++) {
      for(j=0; j<nvlat;j++) {
	delete[] vflow[t][k][j];
      }   
      delete[] vflow[t][k];
    }
    delete[] vflow[t];
  }
}

vectorIJK GetIndices(vectorXYZ point){

  vectorIJK index;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);
  index.k = LocateIndex(point.z, vdepth);

  if(index.i<0)
    index.i=0;
  else if((unsigned int) index.i>vlon.size())
    index.i=vlon.size()-1;

  if(index.j<0)
    index.j=0;
  else if((unsigned int) index.j>vlat.size())
    index.j=vlat.size()-1;

  if(index.k<0)
    index.k=0;
  else if((unsigned int) index.k>vdepth.size())
    index.k=vdepth.size()-1;

  return index;
}

int GetVFlow(double t,vectorXYZ point, vectorXYZ *vint){

  vectorXYZ vgrid[16];
  vectorXYZ vcomp[16];
  vectorIJK index;
  unsigned int time;

  double alpha, beta;
  unsigned int i,j,k,l;
  unsigned int deltai,deltaj,deltak,deltatime;
  unsigned int q;

  time = (unsigned long) t;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);
  index.k = LocateIndex(point.z, vdepth);

  /* Vectors and points with one flat index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++){
    l = time + deltatime;
    for(deltak=0; deltak<2; deltak++){
      k = index.k + deltak;
      for(deltaj=0; deltaj<2; deltaj++){
	j = index.j + deltaj;
	for(deltai=0; deltai<2; deltai++){
	  i = index.i + deltai;

	  try{
	    vgrid[q].x = vlon.at(i);
	    vgrid[q].y = vlat.at(j);
	    vgrid[q].z = vdepth.at(k);
	  }
	  catch(...){
	    return 1;
	    }

	  vcomp[q] = vflow[l][k][j][i];	
	  /* COAST CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	  if(vcomp[q].x <= vflowparams.ufillvalue || vcomp[q].y <= vflowparams.vfillvalue)
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
  /* Depth Interpolation */
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
int GetVPartialDeriv(double t,vectorXYZ point, vectorIJK dir, vectorXYZ *dvdr){

  unsigned int time;
  unsigned int p;
  vectorXYZ delta;
  vectorIJK index;

  vectorXYZ vgrid[16];
  vectorXYZ partialderiv[16];
  double alpha,beta;
  unsigned int i,j,k,l;
  unsigned int i0,j0,k0;
  unsigned int i1,j1,k1;
  unsigned int deltai,deltaj,deltak,deltatime;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);
  index.k = LocateIndex(point.z, vdepth);

  time = (unsigned int) t;

  p=0;
  for(deltatime=0; deltatime<2; deltatime++){
    l = time+deltatime;

    for(deltak=0; deltak<2; deltak++){
      k  = index.k + deltak;
      k0 = k - dir.k;
      k1 = k + dir.k;

      for(deltaj=0; deltaj<2; deltaj++){
	j  =  index.j + deltaj;
	j0 = j - dir.j;
	j1 = j + dir.j;
	
	delta.x=rearth*cos(rads*vlat[j])*vflowparams.lonstep; 
	delta.y=rearth*vflowparams.latstep;
	delta.z=1.0;

	for(deltai=0; deltai<2; deltai++){

	  i  =  index.i + deltai;
	  i0 = i - dir.i;
	  i1 = i + dir.i;

	  try{
	    vgrid[p].x = vlon.at(i);
	    vgrid[p].y = vlat.at(j);
	    vgrid[p].z = vdepth.at(k);

	    if(vgrid[p].x==vlon.front() 
	       || vgrid[p].y==vlon.back()
	       || vgrid[p].x==vlat.front() 
	       || vgrid[p].y==vlat.back()
	       || vgrid[p].x==vdepth.front() 
	       || vgrid[p].y==vdepth.back()){
	      return 1;
	    }
	  }
	  catch(...){
	    return 1;
	  }

	  //COASTAL CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    if(vflow[l][k][j][i].x == vflowparams.ufillvalue 
	     || vflow[l][k][j][i].y == vflowparams.vfillvalue
	     || vflow[l][k0][j0][i0].x == vflowparams.ufillvalue 
	     || vflow[l][k0][j0][i0].y == vflowparams.vfillvalue
	     || vflow[l][k1][j1][i1].x == vflowparams.ufillvalue 
	     || vflow[l][k1][j1][i1].y == vflowparams.vfillvalue){
	    return 1;
	    }

	  partialderiv[p] = (vflow[l][k1][j1][i1] - vflow[l][k0][j0][i0])/(2.0*delta);
	  p++;
	}
      }
    }
  }
  /* Longitude Interpolation: */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1 - alpha;
  for(p=0; p<16; p+=2){
    partialderiv[p>>1] = alpha * partialderiv[p] + beta * partialderiv[p+1];
    vgrid[p>>1] = vgrid[p]; 
  }

  /* Latitude Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(p=0; p<8; p+=2){
    partialderiv[p>>1] = alpha * partialderiv[p] + beta * partialderiv[p+1];
    vgrid[p>>1] = vgrid[p];
  }
 
  /* Depth Interpolation */
  alpha = (vgrid[1].z - point.z)/(vgrid[1].z - vgrid[0].z);
  beta = 1.0 - alpha;
  for(p=0; p<4; p+=2){
    partialderiv[p>>1] = alpha * partialderiv[p] + beta * partialderiv[p+1];
    vgrid[p>>1] = vgrid[p];
  }
 
  /* Time Interpolation: */

  alpha = ((double) (time + 1)) - t;
  beta = 1 - alpha;

  partialderiv[0] = alpha*partialderiv[0]+beta*partialderiv[1];

  *dvdr = partialderiv[0];
  return 0;
}
