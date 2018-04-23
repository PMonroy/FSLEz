#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw
#include <string>

using namespace std;

#include "rparameters.hpp"
#include "gridconstruction.hpp" 
#include "constants.hpp"
#include "vflow.hpp" 
#include "integration.hpp"
#include "vectorXYZ.hpp"

string numprintf(int ndigits, int ndecimals, double number);
struct rtimeParameters {

  const vectorXYZ domainmin;
  const vectorXYZ domainmax;
  const vectorXYZ intergrid;
  const struct tm seeddate;
  const double intstep;
  const int tau;
  const double fthreshold;

  // Define a constructor that will load stuff from a configuration file.
  rtimeParameters(const string & rtimeParamsFileName)
  :domainmin(getVectorXYZParam(rtimeParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(rtimeParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(rtimeParamsFileName, "intergrid")) 
  ,seeddate(getDateParam(rtimeParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(rtimeParamsFileName, "intstep")) 
  ,tau(getIntParam(rtimeParamsFileName, "tau"))
  ,fthreshold(getDoubleParam(rtimeParamsFileName, "fthreshold"))
{}
};

int main(int argc, char **argv){
  
  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/
  
  string fnameparams;// File name that stores the parameters
  int namefileflag;

  if(GetcmdlineParameters(argc, argv, &fnameparams, &namefileflag)){//Get cmd line parameters
    cout << "Error getting parameters from file" <<endl;
    return 1;
  }
  const rtimeParameters rtimeParams(fnameparams);//Set struct parameters

#ifdef DEBUG
  cout<<"Retention time parameters from file ";
  cout <<"\""<<fnameparams<<": "<<endl; 
  cout<<" domainmin = "<<rtimeParams.domainmin<<endl;
  cout<<" intergrid = "<<rtimeParams.intergrid<<endl;
  cout<<" domainmax = "<<rtimeParams.domainmax<<endl;
  cout<<" seeddate = "<<rtimeParams.seeddate.tm_mday<<"-";
  cout<<rtimeParams.seeddate.tm_mon+1<<"-";
  cout<<rtimeParams.seeddate.tm_year<<endl;
  cout<<" intstep = "<<rtimeParams.intstep<<endl;
  cout<<" tau = "<<rtimeParams.tau<<endl;
  cout<<" fthreshold = "<<rtimeParams.fthreshold<<endl;
  cout<<" [Complete]"<<endl;
  cout<< endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector <vectorXYZ> grid;
  vectorIJK dimgrid;
  if(MakeRegularGrid(&grid, &dimgrid, 
		     rtimeParams.domainmin, 
		     rtimeParams.intergrid, 
		     rtimeParams.domainmax)){//Grid construction 
    cout << "[Fail]" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size();
  vector<int> qflag(numgridpoints, 1);

#ifdef DEBUG
  cout << "Grid construction:"<<endl; 
  cout << " num. nodes = "<< grid.size() <<endl;
  cout << " dim(x) = "<< dimgrid.i <<endl;
  cout << " dim(y) = "<< dimgrid.j <<endl;
  cout << " dim(z) = "<< dimgrid.k <<endl;
  cout << "[Complete]" << endl;
  cout << endl;
#endif

  /**********************************************
   * SETUP TIME PARAMETERS
   **********************************************/

  struct tm *inidate = {0};

  double tstart;
  double tend;
  double h;

  int ntime = abs(rtimeParams.tau);
  int ascnd = rtimeParams.tau > 0;
  struct tm datebuff = rtimeParams.seeddate; // Buffer var to avoid error "cast const tm* to tm*"
  time_t seedtime = mktime(&datebuff); // Convert date to time in seconds (UTC) 

  if(ascnd){
    time_t initime = seedtime;
    inidate = gmtime(&initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=rtimeParams.intstep;
  }
  else{
    time_t initime = seedtime - ntime*secondsday;
    inidate = gmtime(&initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*rtimeParams.intstep;
  }

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

#ifdef DEBUG
  cout<<"Loading velocity grid:"<<endl; //Verbose: loading vel. grid
#endif

  vectorXYZ meanvel(0.1,0.1,0.0);
  if(LoadVGrid(rtimeParams.seeddate,
               rtimeParams.domainmin,
               rtimeParams.domainmax, 
               meanvel, 
               abs(rtimeParams.tau))!=0){//Load velocity grid
      cout<<"[Fail]"<<endl;
      return 1;
  }

#ifdef DEBUG
  cout<<"[Complete]" << endl;//Success loading vel. grid
  cout<<"Loading velocities from model:" << endl;//Verbose: loading velocities from model 
#endif

  if(LoadVFlow(*inidate,ntime+2)!=0){// Load velocity field
    cout << "[Fail]"<< endl;
    return 1;
  }

#ifdef DEBUG//Succes loading velocities
  cout << "[Complete]"<<endl;
#endif

  
  /**********************************************
   * TIME LOOP
   **********************************************/
  
  vector <vectorXYZ > tracer = grid;
  vector <vectorXYZ > dvdx;
  vector <vectorXYZ > dvdy;

  vector<double> ow(numgridpoints,0.0);
  vector<double> owt0(numgridpoints,0.0);
  vector<double> std_owt0(dimgrid.k,0.0);
  vector<double> owthreshold(dimgrid.k,0.0);
  vector <double> rtime(numgridpoints,0.0);
  
  vectorIJK dirx(1,0,0),diry(0,1,0);

  dvdx.resize(numgridpoints);
  dvdy.resize(numgridpoints);


  unsigned int step;
  double t;
  unsigned int q;


#ifdef DEBUG//Verbose: Response and exit time calc. 
  cout << "Calculation of initial Okubo Weiss:" << endl;
#endif

  /* Compute the intial OKUBO-WEISS */
  
  for(q=0; q<numgridpoints; q++){
    
    if(GetVPartialDeriv(tstart,tracer[q], dirx, &dvdx[q])==1){
	qflag[q]=-1;
	continue;
    }
    if(GetVPartialDeriv(tstart,tracer[q], diry, &dvdy[q])==1){
      qflag[q]=-1;
      continue;
    }
    
    owt0[q]=(dvdx[q].x-dvdy[q].y)*(dvdx[q].x-dvdy[q].y)+4.0*(dvdy[q].x*dvdx[q].y);
  }
  
  /* COMPUTE:Threshold, compute standard deviation OW at initial date */
  unsigned int nlayer=dimgrid.i*dimgrid.j;

  double sum_owt0=0.0, sum2_owt0=0.0;
  double nsize_owt0=0.0;

  unsigned int k;
  unsigned int dimk=dimgrid.k;

  for(k=0; k<dimk; k++){
    sum_owt0=0.0;
    sum2_owt0=0.0;
    nsize_owt0=0.0;
    for(q=k*nlayer; q<(k+1)*nlayer; q++){
      if(qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point
	continue;
      }
      nsize_owt0+=1.0;
      sum_owt0+=owt0[q];
      sum2_owt0+=(owt0[q]*owt0[q]);
    }
    std_owt0[k]=sqrt((sum2_owt0/nsize_owt0)-((sum_owt0/nsize_owt0)*(sum_owt0/nsize_owt0)));
    owthreshold[k]=rtimeParams.fthreshold*std_owt0[k];

  }


  for(q=0; q<numgridpoints; q++){

    if(qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
      continue;
    }
    k=(unsigned int)(q/nlayer);
    if(owt0[q]>=owthreshold[k]){
      rtime[q]=0.0;
      qflag[q]=0;
    }
  }

#ifdef DEBUG//Verbose: Response and exit time calc. 
  cout << "Calculation of eddy retention time:" << endl;
#endif

  for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++){

#ifdef DEBUG//Verbose: show the current step 
    cout << "step=" << step << "(" <<(tend-tstart)/h <<")"<<endl;
#endif

    /* Compute the position of tracer in time t=t+h*/
    for (q=0; q<numgridpoints; q++){
      if(qflag[q]==0 || qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
      if(RK4(t, h, &tracer[q], GetVFlow)==1){
	qflag[q]=-1; // q is a non-integrable grid point
      }
    }

    /* Compute the OKUBO-WEISS */
    for (q=0; q<numgridpoints; q++){

      if(qflag[q]==0 || qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
      if(GetVPartialDeriv(t,tracer[q], dirx, &dvdx[q])==1){
	qflag[q]=-1;
	continue;
      }
      if(GetVPartialDeriv(t,tracer[q], diry, &dvdy[q])==1){
	qflag[q]=-1;
	continue;
      }
      

      ow[q]=(dvdx[q].x-dvdy[q].y)*(dvdx[q].x-dvdy[q].y)+4.0*(dvdy[q].x*dvdx[q].y);//Okubo-Weiss Parameters
      k=(unsigned int)(q/nlayer);
      if(ow[q]>=owthreshold[k]){
	rtime[q]=step*h;
	qflag[q]=0;
      }
    }
  }

  //free vflow
  FreeVFlow((unsigned int)(ntime+2));

#ifdef DEBUG//Verbose: Success response and exit time calculation
    cout << "[Complete]" << endl;
#endif


  /**********************************************************
   * WRITE RESULTS
   **********************************************************/

  // Save fslez grid in a file
  string nfilegridrtime = 
    "rtime_lon"    + numprintf(4,0,rtimeParams.domainmin.x) 
    + numprintf(4,0,rtimeParams.domainmax.x)
    + numprintf(4,3,rtimeParams.intergrid.x)
    + "_lat"       + numprintf(4,0,rtimeParams.domainmin.y)
    + numprintf(4,0,rtimeParams.domainmax.y)
    + numprintf(4,3,rtimeParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,rtimeParams.domainmin.z)
    + numprintf(4,0,rtimeParams.domainmax.z)
    + numprintf(4,0,rtimeParams.intergrid.z)
    + ".grid";

  #ifdef DEBUG
  cout<<"Saving grid in file: " << nfilegridrtime <<endl;
  #endif

  if(!ifstream(nfilegridrtime.c_str())){ // Check if file exists
    ofstream ofilegridrtime(nfilegridrtime.c_str());
    for(q=0; q<numgridpoints; q++){
      ofilegridrtime<<grid[q]<<endl;
    }
    ofilegridrtime.close();
  }

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  // save rtime values in a file
  
  string rawname;
  
  if(namefileflag==1){
    size_t lastdot=fnameparams.find_last_of(".");
    if(lastdot==string::npos){
      rawname="rtime3d_"+fnameparams;
    }
    else{
      rawname="rtime3d_"+fnameparams.substr(0,lastdot);
    }
  }
  else{
    rawname= "rtime3d_lon" + numprintf(4,0,rtimeParams.domainmin.x) 
      + numprintf(4,0,rtimeParams.domainmax.x)
      + numprintf(4,3,rtimeParams.intergrid.x)
      + "_lat"       + numprintf(4,0,rtimeParams.domainmin.y)
      + numprintf(4,0,rtimeParams.domainmax.y)
      + numprintf(4,3,rtimeParams.intergrid.y)
      + "_dpt"       + numprintf(4,0,rtimeParams.domainmin.z)
      + numprintf(4,0,rtimeParams.domainmax.z)
      + numprintf(4,0,rtimeParams.intergrid.z)
      + "_h"        + numprintf(4,3,rtimeParams.intstep) 
      + "_t"         + numprintf(4,0,rtimeParams.tau)
      + "_f"         + numprintf(3,1,rtimeParams.fthreshold) 
      + "_d"         + numprintf(2,0,rtimeParams.seeddate.tm_mday)
      + numprintf(2,0,rtimeParams.seeddate.tm_mon+1)
      + numprintf(2,0,rtimeParams.seeddate.tm_year);
  }

  string nfilertime = rawname + ".data";

#ifdef DEBUG// Verbose: Save fsle values in ascii file
  cout << "Save rtime/owt0/qflag values in file: " << nfilertime <<endl;
#endif

  ofstream ofilertime(nfilertime.c_str());
  for(q=0; q<numgridpoints; q++){
    k=(unsigned int)(q/nlayer);
    ofilertime<<rtime[q]<<" ";
    ofilertime<<owt0[q]<<" ";
    ofilertime<<qflag[q]<<" ";
    ofilertime<<std_owt0[k]<<" ";
    ofilertime<<owthreshold[k]<<endl;
  }
  ofilertime.close();

#ifdef DEBUG// Verbose: success file saved  
  cout << "[Complete]" << endl;
#endif


 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/

  string vtkfilertime = rawname + ".vtk";

#ifdef DEBUG
  cout << "Save rtime field in vtk file: " << vtkfilertime <<endl;
#endif

  ofstream vtkfile(vtkfilertime.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Retention time 3D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_GRID"<< endl;
  vtkfile<<"DIMENSIONS "<< dimgrid.i <<" "<< dimgrid.j <<" "<< dimgrid.k<<endl;
  vtkfile<<"POINTS "<< numgridpoints <<" double"<<endl;
  vectorXYZ scalegrid(1.0,1.0,-0.01);
  for(q=0; q<numgridpoints; q++){
    vtkfile<<scalegrid*grid[q]<<endl;
  }  
  vtkfile<<"POINT_DATA "<< numgridpoints <<endl;
  vtkfile<<"SCALARS rtime double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++){
    vtkfile<<rtime[q]<<endl;
  }
  vtkfile<<"SCALARS owt0 double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++){
    vtkfile<<owt0[q]<<endl;
  }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++){
    vtkfile<<qflag[q]<<endl;
  }
  vtkfile.close();

  #ifdef DEBUG
    cout << "[Complete]" <<endl;
  #endif

  return 0;
}

string numprintf(int ndigits, int ndecimals, double number){
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
