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

  // Define a constructor that will load stuff from a configuration file.
  rtimeParameters(const string & rtimeParamsFileName)
  :domainmin(getVectorXYZParam(rtimeParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(rtimeParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(rtimeParamsFileName, "intergrid")) 
  ,seeddate(getDateParam(rtimeParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(rtimeParamsFileName, "intstep")) 
  ,tau(getIntParam(rtimeParamsFileName, "tau"))
{}
};

int main(int argc, char **argv){
  
  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/
  
  string fnameparams;// File name that stores the parameters
  if(GetcmdlineParameters(argc, argv, &fnameparams)){//Get cmd line parameters
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

  vector <double> ow;
  vector <double> rtime;
  
  vectorIJK dirx(1,0,0),diry(0,1,0);

  dvdx.resize(numgridpoints);
  dvdy.resize(numgridpoints);
  ow.resize(numgridpoints);
  rtime.resize(numgridpoints);

  unsigned int step;
  double t;
  unsigned int q;

#ifdef DEBUG//Verbose: Response and exit time calc. 
  cout << "Calculation of retention time:" << endl;
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
      

      ow[q]=(dvdx[q].x-dvdy[q].y)*(dvdx[q].x -dvdy[q].y)+4.0*(dvdy[q].x*dvdx[q].y);//Okubo-Weiss Parameters
      if(ow[q]>=0.0){
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
   * METRICS : must be moved to a separated file!!!!!!!!!!!!           
   **********************************************************/
  string nmeandpt = 
    "rtime_date"  + numprintf(2,0,rtimeParams.seeddate.tm_year) + 
    "-"           + numprintf(2,0,rtimeParams.seeddate.tm_mon+1)+
    "-"           + numprintf(2,0,rtimeParams.seeddate.tm_mday) +
    "_dmin"       + numprintf(4,0,rtimeParams.domainmin.x) +
    "_"           + numprintf(4,0,rtimeParams.domainmin.y) +
    "_"           + numprintf(4,0,rtimeParams.domainmin.z) +
    "_dmax"       + numprintf(4,0,rtimeParams.domainmax.x) +
    "_"           + numprintf(4,0,rtimeParams.domainmax.y) +
    "_"           + numprintf(4,0,rtimeParams.domainmax.z) +
    "_res"        + numprintf(4,3,rtimeParams.intergrid.x) +
    "_"           + numprintf(4,3,rtimeParams.intergrid.y) +
    "_"           + numprintf(3,0,rtimeParams.intergrid.z) +
    "_tau"        + numprintf(4,0,rtimeParams.tau) +
    "_h"          + numprintf(4,3,rtimeParams.intstep) +
    ".meandpt";

  #ifdef DEBUG// Verbose: Save fsle values in ascii file
    cout << "Save mean values in file: " << nmeandpt <<endl;
  #endif

  ofstream ofilemeandpt(nmeandpt.c_str());

  unsigned int idxdepth;
  
  double sumdepth,sum2depth;
  double meandepth,mean2depth;
  double variancedepth;

  for(idxdepth=0; idxdepth<(unsigned int) dimgrid.k; idxdepth++)
    {
      sumdepth=0.0;
      sum2depth=0.0;
      for(q=idxdepth*dimgrid.i*dimgrid.j; q<(idxdepth+1)*dimgrid.i*dimgrid.j; q++){

	sumdepth+=rtime[q];
	sum2depth+=(rtime[q]*rtime[q]);
      }
      meandepth=sumdepth/double(dimgrid.i*dimgrid.j);
      mean2depth=sum2depth/double(dimgrid.i*dimgrid.j);

      variancedepth=mean2depth-(meandepth*meandepth); 
      ofilemeandpt << grid[idxdepth*dimgrid.i*dimgrid.j].z << " ";
      ofilemeandpt << meandepth << " ";
      ofilemeandpt << variancedepth << endl;
      //cout << grid[idxdepth*dimgrid.i*dimgrid.j].z << " " << meandepth << endl;
    }
  ofilemeandpt.close();

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
   ofstream ofilegridrtime(nfilegridrtime.c_str());

  #ifdef DEBUG
  cout<<"Saving grid in file: " << nfilegridrtime <<endl;
  #endif

  for(q=0; q<numgridpoints; q++){
      ofilegridrtime<<grid[q]<<endl;
  }
  ofilegridrtime.close();

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  // save fslez values in a file
  string nfilertime = 
     "rtime_lon"    + numprintf(4,0,rtimeParams.domainmin.x) 
    + numprintf(4,0,rtimeParams.domainmax.x)
    + numprintf(4,3,rtimeParams.intergrid.x)
    + "_lat"       + numprintf(4,0,rtimeParams.domainmin.y)
    + numprintf(4,0,rtimeParams.domainmax.y)
    + numprintf(4,3,rtimeParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,rtimeParams.domainmin.z)
    + numprintf(4,0,rtimeParams.domainmax.z)
    + numprintf(4,0,rtimeParams.intergrid.z)
    + "_ts"        + numprintf(4,0,rtimeParams.intstep) 
    + "_t"         + numprintf(4,0,rtimeParams.tau) 
    + "_d"         + numprintf(2,0,rtimeParams.seeddate.tm_mday)
    + numprintf(2,0,rtimeParams.seeddate.tm_mon+1)
    + numprintf(2,0,rtimeParams.seeddate.tm_year)
    + ".data";

#ifdef DEBUG// Verbose: Save fsle values in ascii file
  cout << "Save rtime values in file: " << nfilertime <<endl;
#endif

  ofstream ofilertime(nfilertime.c_str());
  for(q=0; q<numgridpoints; q++){
    ofilertime<<rtime[q]<<endl;
  }
  ofilertime.close();

#ifdef DEBUG// Verbose: success file saved  
  cout << "[Complete]" << endl;
#endif


 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/

  string vtkfilertime = 
    "rtime_lon"    + numprintf(4,0,rtimeParams.domainmin.x) 
    + numprintf(4,0,rtimeParams.domainmax.x)
    + numprintf(4,3,rtimeParams.intergrid.x)
    + "_lat"       + numprintf(4,0,rtimeParams.domainmin.y)
    + numprintf(4,0,rtimeParams.domainmax.y)
    + numprintf(4,3,rtimeParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,rtimeParams.domainmin.z)
    + numprintf(4,0,rtimeParams.domainmax.z)
    + numprintf(4,0,rtimeParams.intergrid.z)
    + "_ts"        + numprintf(4,0,rtimeParams.intstep) 
    + "_t"         + numprintf(4,0,rtimeParams.tau) 
    + "_d"         + numprintf(2,0,rtimeParams.seeddate.tm_mday)
    + numprintf(2,0,rtimeParams.seeddate.tm_mon+1)
    + numprintf(2,0,rtimeParams.seeddate.tm_year)
    + ".data";

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
  vectorXYZ scalegrid(1.0,1.0,-1.0);
  for(q=0; q<numgridpoints; q++){
    vtkfile<<scalegrid*grid[q]<<endl;
  }  
  vtkfile<<"POINT_DATA "<< numgridpoints <<endl;
  vtkfile<<"SCALARS rtime double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++){
    vtkfile<<rtime[q]<<endl;
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
