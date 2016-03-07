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

string numprintf(int ndigits, int ndecimals, double number);
struct fslezParameters {

  const vectorXYZ domainmin;
  const vectorXYZ domainmax;
  const vectorXYZ intergrid;
  const struct tm seeddate;
  const double intstep;
  const int tau;
  const double deltamax;

  // Define a constructor that will load stuff from a configuration file.
  fslezParameters(const string & rtimeParamsFileName)
  :domainmin(getVectorXYZParam(rtimeParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(rtimeParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(rtimeParamsFileName, "intergrid")) 
  ,seeddate(getDateParam(rtimeParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(rtimeParamsFileName, "intstep")) 
  ,tau(getIntParam(rtimeParamsFileName, "tau"))
  ,deltamax(getDoubleParam(rtimeParamsFileName, "deltamax"))
{}
};

int main(int argc, char **argv){

  /*********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string fnameparams;// File name that stores the parameters
  if(GetcmdlineParameters(argc, argv, &fnameparams)){//Get cmd line parameters
    cout << "Error getting parameters from file" <<endl;
    return 1;
  }
  const fslezParameters fslezParams(fnameparams);


#ifdef DEBUG
  cout << "FSLEz Parameters from file ";
  cout <<"\""<<fnameparams<<": "<<endl; 
  cout << " domainmin = "<< fslezParams.domainmin<<endl;
  cout << " intergrid = " <<fslezParams.intergrid<<endl;
  cout << " domainmax = "<< fslezParams.domainmax<<endl;
  cout << " seeddate = "<< fslezParams.seeddate.tm_mday<<"-";
  cout<<fslezParams.seeddate.tm_mon+1<<"-";
  cout<<fslezParams.seeddate.tm_year<<endl;
  cout << " intstep = "<<fslezParams.intstep<<endl ;
  cout << " tau = "<<fslezParams.tau<<endl ;
  cout << " deltamax = "<<fslezParams.deltamax<<endl ;
  cout << " [Complete]" <<endl;
  cout << endl;
#endif
      
  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  vectorIJK dimgrid;
  vector<int> neighbor;
  if(MakeRegularGrid(&grid, &dimgrid, 
		     fslezParams.domainmin, 
		     fslezParams.intergrid, 
		     fslezParams.domainmax)){//Grid construction 
    cout << "[Fail]" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size();

  neighbor=neighbors(dimgrid);
  unsigned int numneighbors=neighbor.size();

#ifdef DEBUG
  cout << "Grid construction: "<<endl; 
  cout << " num. nodes = "<< grid.size() <<endl;
  cout << " dim(x) = "<< dimgrid.i <<endl;
  cout << " dim(y) = "<< dimgrid.j <<endl;
  cout << " dim(z) = "<< dimgrid.k <<endl;
  cout << "[Complete]" << endl;
  cout << endl;
#endif

  /**********************************************
   * COMPUTE INITIAL SEPARATION
   **********************************************/

  vectorXYZ delta,scalefactor;
  vector<double> ilength;
  unsigned int q;
  unsigned int p;
  unsigned int dir;

  
  ilength.reserve(numneighbors);

  for(q=0; q<numgridpoints; q++){
    p=6*q;
    for(dir=0; dir<6; dir++){
      if(neighbor[p+dir]!=-1){
	delta=grid[neighbor[p+dir]]-grid[q];
	
	delta.x=rads*delta.x;
	delta.y=rads*delta.y;
	
	scalefactor.x=rearth*cos(rads*grid[q].y); 
	scalefactor.y=rearth; 
	scalefactor.z=1.0; 
	
	delta=delta*scalefactor;
	delta*=delta;
	
	ilength.push_back(sqrt(delta.x+delta.y+delta.z));
      }
      else
	ilength.push_back(-1.0);
    }
  }

  /**********************************************
   * INITIALIZE VARIABLES
   **********************************************/

  vector<double> exit_time;
  vector<double> response;
  vector<int> qcore;
  vector<int> qflag;

  //Reserve memory to improve efficiency
  exit_time.reserve(numgridpoints);
  response.reserve(numgridpoints);
  qcore.reserve(numgridpoints);
  qflag.reserve(numgridpoints);

  for(q=0; q<numgridpoints; q++){
    p=6*q;
    if(neighbor[p]!=-1 && neighbor[p+1]!=-1 && neighbor[p+3]!=-1 && neighbor[p+4]!=-1){
      qflag.push_back(1);// Initially we can compute fsle in this points (num. hor. neighbours = 4)
      qcore.push_back(q);// qcore contains all the index of the later points
    }
    else{
      qflag.push_back(0); // Border points: we can not compute fsle on it (num. hor. neighbours <4)
    }
    //Response and exit_time are initialized to 1.0 (so initially fsle is 0.0)
    response.push_back(1.0);
    exit_time.push_back(1.0);
  }

  /**********************************************
   * SETUP TIME PARAMETERS
   **********************************************/

  struct tm *inidate = {0};

  double tstart;
  double tend;
  double h;

  int ntime = abs(fslezParams.tau);
  int ascnd = fslezParams.tau > 0;
  struct tm datebuff = fslezParams.seeddate; // Buffer var to avoid error "cast const tm* to tm*"
  time_t seedtime = mktime(&datebuff); // Convert date to time in seconds (UTC) 

  if(ascnd){
    time_t initime = seedtime;
    inidate = gmtime(&initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=fslezParams.intstep;
  }
  else{
    time_t initime = seedtime - ntime*secondsday;
    inidate = gmtime(&initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*fslezParams.intstep;
  }

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

#ifdef DEBUG
  cout<<"Loading velocity grid:"<<endl; //Verbose: loading vel. grid
#endif

  vectorXYZ meanvel(0.2,0.2,0.0);
  if(LoadVGrid(fslezParams.seeddate,
               fslezParams.domainmin,
               fslezParams.domainmax, 
               meanvel, 
               abs(fslezParams.tau))!=0){//Load velocity grid
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
  
  vector<vectorXYZ> tracer = grid;  
  vector<double> length = ilength;

  double lengthmax;
  int dirmax;
  
  int q0, q1, q3, q4; 
  int qdir;
  int pdir;

  unsigned int step;
  double t;

#ifdef DEBUG//Verbose: Response and exit time calc. 
  cout << "Calculation of response and exit_time:" << endl;
#endif

  for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++) {

#ifdef DEBUG//Verbose: show the current step 
    cout << "step=" << step << "(" <<(tend-tstart)/h <<")"<<endl;
#endif

    /* Compute the position of tracer in time t=t+h*/
    for (q=0; q<numgridpoints; q++)
      {
	//index of four q-neighbors
	p=6*q;
	q0=neighbor[p];
	q1=neighbor[p+1];
	q3=neighbor[p+3];
	q4=neighbor[p+4];
	
	if(qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
	if((qflag[q]==1) ||
	   (q0!=-1 && qflag[q0]==1) ||
	   (q1!=-1 && qflag[q1]==1) ||
	   (q3!=-1 && qflag[q3]==1) ||
	   (q4!=-1 && qflag[q4]==1)){
	  if(RK4(t, h, &tracer[q], GetVFlow)==1){
	    qflag[q]=-1; // q is a non-integrable grid point
	    if(q0!=-1 && qflag[q0]!=-1) qflag[q0]=0;// We keep non-integrable neighbor grid points
	    if(q1!=-1 && qflag[q1]!=-1) qflag[q1]=0;
	    if(q3!=-1 && qflag[q3]!=-1) qflag[q3]=0;
	    if(q4!=-1 && qflag[q4]!=-1) qflag[q4]=0;
	  }
	}
      }

    /* Compute the relative distances */
    for(q=0; q<numgridpoints; q++)
      {
	p=6*q;
	for(dir=0; dir<2; dir++){
	  qdir = neighbor[p+dir];
	  if(qdir==-1) // if grid point hasn't neighbour in direction dir,
	    continue;  // jump to the next step in the loop
	  
	  if(qflag[q]==1 || qflag[qdir]==1){
	    delta=tracer[qdir]-tracer[q];
	    
	    delta.x=rads*delta.x;
	    delta.y=rads*delta.y;
	    
	    scalefactor.x=rearth*cos(rads*grid[q].y); 
	    scalefactor.y=rearth;
	    scalefactor.z=1.0;
	    
	    delta=delta*scalefactor;
	    delta*=delta;

	    pdir=p+dir; 
	    length[pdir]=sqrt(delta.x+delta.y+delta.z);
	    length[6*qdir+3]=length[pdir];
	  }
	}
      }

    /* Compute the max length*/
    for(q=0; q<numgridpoints; q++)
      {
	if(qflag[q]==1){
	  p=6*q;
	  lengthmax=length[p];// May we initialize lengthmax before?
	  dirmax=0;
	  for(dir=0; dir<2; dir++){
	    pdir = p+dir;
	    if(length[pdir]>lengthmax){
	      lengthmax=length[pdir];
	      dirmax=dir;
	    }
	    if(length[pdir+3]>lengthmax){
	      lengthmax=length[pdir+3];
	      dirmax=dir+3;
	    }
	  }
	  
	  if(lengthmax>fslezParams.deltamax){
	    qflag[q]=0;
	    exit_time[q]=(step+1.0)*h;
	    response[q]=lengthmax/ilength[p+dirmax];
	  }
	}
      }
  }

  //free VFlow
#ifdef DEBUG//Verbose: Success response and exit time calculation
  cout << "[Complete]" << endl;
#endif

  /**********************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   **********************************************************/

  vector<double> fsle;
  unsigned int numcorepoints=qcore.size();

  fsle.reserve(numgridpoints);
  for(q=0; q<numcorepoints; q++){
    fsle.push_back(log(response[qcore[q]])/exit_time[qcore[q]]);
  }

  // Save fslez grid in a file
  string nfilegridfsle2d = 
    "fslez_lon"    + numprintf(4,0,fslezParams.domainmin.x) 
    + numprintf(4,0,fslezParams.domainmax.x)
    + numprintf(4,3,fslezParams.intergrid.x)
    + "_lat"       + numprintf(4,0,fslezParams.domainmin.y)
    + numprintf(4,0,fslezParams.domainmax.y)
    + numprintf(4,3,fslezParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,fslezParams.domainmin.z)
    + numprintf(4,0,fslezParams.domainmax.z)
    + numprintf(4,0,fslezParams.intergrid.z)
    + ".grid";
  ofstream ofilegridfsle2d(nfilegridfsle2d.c_str());

#ifdef DEBUG 
  cout << "Save grid in file: " << nfilegridfsle2d <<endl;
#endif

  for(q=0; q<numcorepoints; q++){
      ofilegridfsle2d<<grid[q]<<endl;
  }
  ofilegridfsle2d.close();

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  // save fslez values in a file
  string nfilefsle2d = 
    "fslez_lon"    + numprintf(4,0,fslezParams.domainmin.x) 
    + numprintf(4,0,fslezParams.domainmax.x)
    + numprintf(4,3,fslezParams.intergrid.x)
    + "_lat"       + numprintf(4,0,fslezParams.domainmin.y)
    + numprintf(4,0,fslezParams.domainmax.y)
    + numprintf(4,3,fslezParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,fslezParams.domainmin.z)
    + numprintf(4,0,fslezParams.domainmax.z)
    + numprintf(4,0,fslezParams.intergrid.z)
    + "_dl"        + numprintf(3,0,fslezParams.deltamax/1000.0)
    + "_ts"        + numprintf(4,0,fslezParams.intstep) 
    + "_t"         + numprintf(4,0,fslezParams.tau) 
    + "_d"         + numprintf(2,0,fslezParams.seeddate.tm_mday)
    + numprintf(2,0,fslezParams.seeddate.tm_mon+1)
    + numprintf(2,0,fslezParams.seeddate.tm_year)
    + ".data";

#ifdef DEBUG
  // Verbose: Save fsle values in ascii file
  cout << "Save fsle values in file: " << nfilefsle2d <<endl;
#endif

  ofstream ofilefsle2d(nfilefsle2d.c_str());
  for(q=0; q<numcorepoints; q++){
    ofilefsle2d<<fsle[q]<<endl;
  }
  ofilefsle2d.close();

#ifdef DEBUG
  // Verbose: success file saved  
  cout << "[Complete]" << endl;
#endif

 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/

  string vtkfilefsle2d = 
    "fslez_lon"       + numprintf(4,0,fslezParams.domainmin.x) 
    + numprintf(4,0,fslezParams.domainmax.x)
    + numprintf(4,3,fslezParams.intergrid.x)
    + "_lat"       + numprintf(4,0,fslezParams.domainmin.y)
    + numprintf(4,0,fslezParams.domainmax.y)
    + numprintf(4,3,fslezParams.intergrid.y)
    + "_dpt"       + numprintf(4,0,fslezParams.domainmin.z)
    + numprintf(4,0,fslezParams.domainmax.z)
    + numprintf(4,0,fslezParams.intergrid.z)
    + "_dl"        + numprintf(3,0,fslezParams.deltamax/1000.0)
    + "_ts"        + numprintf(4,0,fslezParams.intstep) 
    + "_t"         + numprintf(4,0,fslezParams.tau) 
    + "_d"         + numprintf(2,0,fslezParams.seeddate.tm_mday)
    + numprintf(2,0,fslezParams.seeddate.tm_mon+1)
    + numprintf(2,0,fslezParams.seeddate.tm_year)
    + ".vtk";

#ifdef DEBUG
  cout << "Save ftle field in vtk file: " << vtkfilefsle2d <<endl;
#endif

  ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 3D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< dimgrid.i-1 <<" "<< dimgrid.j-1 <<" "<< dimgrid.k+1 <<endl;
  vtkfile<<"ORIGIN ";
  vtkfile<<grid[qcore[0]].x-(fslezParams.intergrid.x/2.0)<<" ";
  vtkfile<<grid[qcore[0]].y-(fslezParams.intergrid.y/2.0)<<" ";
  vtkfile<<grid[qcore[0]].z-(fslezParams.intergrid.z/2.0)<<" "; 
  vtkfile<<"SPACING ";
  vtkfile<<fslezParams.intergrid.x<<" ";
  vtkfile<<fslezParams.intergrid.y<<" ";
  vtkfile<<-fslezParams.intergrid.z <<endl;
  vtkfile<<"CELL_DATA "<<(dimgrid.i-2)*(dimgrid.j-2)*(dimgrid.k)<<endl;
  vtkfile<<"SCALARS fslez double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numcorepoints; q++){
    vtkfile<<fsle[q]<<endl;
  }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numcorepoints; q++){
    vtkfile<<qflag[qcore[q]]<<endl;
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
