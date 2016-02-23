/* fslez.cpp*/
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw

using namespace std;

#include "readparameters.hpp"
#include "gridconstruction.hpp" 
#include "constants.hpp"
#include "vflow.hpp" 
#include "integration.hpp"

// Maybe one day I use this:
//int (*velocity)(double ,vectorXYZ , vectorXYZ* ); 

string numprintf(int ndigits, int ndecimals, double number);

int main(int argc, char **argv){

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/

  string fnameparams;// File name that stores the parameters
  if(verbose==1){//Verbose: reading command line parameters
      cout << "Reading command line parameters:" <<endl;
      cout << endl;
    }
  if(GetcmdlineParameters(argc, argv, &fnameparams)){//Get cmd parameters
    cout << "[Fail]" <<endl;
    return 1;
  }
  if(verbose==1){//Verbose: success reading cmd parameters
      cout << " parameters file: " << fnameparams <<endl;
      cout << " [Complete]" <<endl;
      cout << endl;
    }

  /**********************************************
   * READ PARAMETERS FROM FILE
   **********************************************/

  if(verbose==1){//Verbose: reading parameters from file
    cout << "Reading parameters from file:"<< fnameparams <<endl;
    }
  if(GetfileParameters(fnameparams)){//Get parameters from file
    cout << "[Fail]" <<endl;    
    return 1;
  }
  if(verbose==1){//Verbose: success reading file parameters
      cout << "Parameters: "<<endl; 
      cout << " vfield = "<<vfield<<endl ;
      cout << " domainmin = "<< domainmin<<endl;
      cout << " intergrid = " <<intergrid<<endl;
      cout << " domainmax = "<< domainmax<<endl;
      cout << " seeddate = "<< seeddate.tm_mday<<"-"<<seeddate.tm_mon+1<<"-"<<seeddate.tm_year<<endl ;
      cout << " intstep = "<<intstep<<endl ;
      cout << " tau = "<<tau<<endl ;
      cout << " deltamax = "<<deltamax<<endl ;
      cout << " [Complete]" <<endl;
      cout << endl;
    }

  /**********************************************
   * GRID CONSTRUCTION
   **********************************************/

  vector<vectorXYZ> grid;
  vectorIJK dimgrid;
  vector<int> neighbor;

  if(verbose==1){//Verbose: grid construction 
    cout << "Grid construction:"<<endl;
  }  
  if(gridfsle2d(&grid, &dimgrid, &neighbor, domainmin, intergrid, domainmax)){//Grid construction 
    cout << "[Fail]" << endl;
    return 1;
  }
  if(verbose == 1){//Verbose: success grid construction
    cout << "Grid parameters: "<<endl; 
    cout << " num. nodes = "<< grid.size() <<endl;
    cout << " dim(x) = "<< dimgrid.i <<endl;
    cout << " dim(y) = "<< dimgrid.j <<endl;
    cout << " dim(z) = "<< dimgrid.k <<endl;
    cout << "[Complete]" << endl;
    cout << endl;
  }

  /**********************************************
   * COMPUTE INITIAL SEPARATION
   **********************************************/

  vectorXYZ delta,scalefactor;
  vector<double> ilength;
  unsigned int p;

  ilength.reserve(6*grid.size());

  for(unsigned int q=0; q<grid.size(); q++){
    p=6*q;
    for(int dir=0; dir<6; dir++){
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
  exit_time.reserve(grid.size());
  response.reserve(grid.size());
  qcore.reserve(grid.size());
  qflag.reserve(grid.size());

  for(unsigned int q=0; q<grid.size(); q++){
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

  int ntime = abs(tau);
  int ascnd = tau > 0;
  time_t seedtime = mktime(&seeddate); // Convert date to time in seconds (UTC) 

  if(ascnd){
    time_t initime = seedtime;
    inidate = gmtime(&initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=intstep;
  }
  else{
    time_t initime = seedtime - ntime*secondsday;
    inidate = gmtime(&initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*intstep;
  }

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

  if(verbose==1){//Verbose: loading vel. grid
    cout << "Loading velocity grid:" << endl;
  }
  if(loadvgrid(seeddate,vfield)!=0){//Load velocity grid
      cout << "[Fail]" << endl;
      return 1;
  }
  if(verbose == 1){//Success loading vel. grid
    cout << "[Complete]" << endl;
  }

  if(verbose == 1){//Verbose: loading velocities from model 
    cout << "Loading velocities from model:" << endl;
  }
  if(loadvflow(*inidate, ntime+2, vfield)!=0){// Load velocity field
    cout << "[Fail]"<< endl;
    return 1;
  }
  if(verbose == 1){//Succes loading velocities
    cout << "[Complete]"<<endl;
  }

  /**********************************************
   * TIME LOOP
   **********************************************/
  
  vector<vectorXYZ> tracer = grid;  
  vector<double> length = ilength;

  double lengthmax;
  int dirmax;
  
  int q0, q1, q3, q4; 
  int qdir;

  unsigned int step;
  double t;

  if(verbose == 1){//Verbose: Response and exit time calc. 
    cout << "Calculation of response and exit_time:" << endl;
  }
  for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++) {

    if(verbose == 1){//Verbose: show the current step 
      cout << "step=" << step << "(" <<(tend-tstart)/h <<")"<<endl;
    }

    /* Compute the position of tracer in time t=t+h*/
    for (unsigned int q = 0; q<tracer.size() ; q++)
      {
	//index of four q-neighbors
	q0 = neighbor[6*q];
	q1 = neighbor[6*q+1];
	q3 = neighbor[6*q+3];
	q4 = neighbor[6*q+4];
	
	if(qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
	if((qflag[q]==1) ||
	   (q0!=-1 && qflag[q0]==1) ||
	   (q1!=-1 && qflag[q1]==1) ||
	   (q3!=-1 && qflag[q3]==1) ||
	   (q4!=-1 && qflag[q4]==1)){
	  if(RK4(t, h, &tracer[q], getvflow, vfield)==1){
	    qflag[q]=-1; // q is a non-integrable grid point
	    if(q0!=-1 && qflag[q0]!=-1) qflag[q0]=0;// We keep non-integrable neighbor grid points
	    if(q1!=-1 && qflag[q1]!=-1) qflag[q1]=0;
	    if(q3!=-1 && qflag[q3]!=-1) qflag[q3]=0;
	    if(q4!=-1 && qflag[q4]!=-1) qflag[q4]=0;
	  }
	}
      }

    /* Compute the relative distances */
    for(unsigned int q=0; q<tracer.size(); q++)
      {
	p=6*q;
	for(int dir=0; dir<2; dir++){
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
	  
	    length[p+dir]=sqrt(delta.x+delta.y+delta.z);
	    length[6*qdir+3]=length[p+dir];
	  }
	}
      }

    /* Compute the max length*/
    for(unsigned int q=0; q<tracer.size(); q++)
      {
	if(qflag[q]==1){
	  p=6*q;
	  lengthmax=length[p];// May we initialize lengthmax before?
	  for(int dir=0; dir<2; dir++){
	    if(length[p+dir]>lengthmax){
	      lengthmax=length[p+dir];
	      dirmax=dir;
	    }
	    if(length[p+dir+3]>lengthmax){
	      lengthmax=length[p+dir+3];
	      dirmax=dir+3;
	    }
	  }
	  
	  if(lengthmax>deltamax){
	    qflag[q]=0;
	    exit_time[q]=(step+1.0)*h;
	    response[q]=lengthmax/ilength[p+dirmax];
	  }
	}
      }
  }
  if(verbose==1){//Verbose: Success response and exit time calculation
    cout << "[Complete]" << endl;
  }

  /**********************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   **********************************************************/

  vector<double> fsle;

  fsle.reserve(grid.size());
  for(unsigned int q=0; q<qcore.size(); q++){
    fsle.push_back(log(response[qcore[q]])/exit_time[qcore[q]]);
  }

  // Save fslez grid in a file
  string nfilegridfsle2d = 
    "fslez_dmin" + numprintf(4,0,domainmin.x) +
    "_"          + numprintf(4,0,domainmin.y) +
    "_"          + numprintf(4,0,domainmin.z) +
    "_dmax"      + numprintf(4,0,domainmax.x) +
    "_"          + numprintf(4,0,domainmax.y) +
    "_"          + numprintf(4,0,domainmax.z) +
    "_res"       + numprintf(4,3,intergrid.x) +
    "_"          + numprintf(4,3,intergrid.y) +
    "_"          + numprintf(4,3,intergrid.z) +
    ".grid";
  ofstream ofilegridfsle2d(nfilegridfsle2d.c_str());

  if(verbose==1) cout << "Save grid in file: " << nfilegridfsle2d <<endl;
  for(unsigned int q=0; q<qcore.size(); q++){
      ofilegridfsle2d<<grid[q]<<endl;
  }
  ofilegridfsle2d.close();

  if(verbose==1) cout << "[Complete]" << endl;

  // save fslez values in a file
  string nfilefsle2d = 
    "fslez_vm"   + numprintf(1,0,vfield) +
    "_date"       + numprintf(2,0,seeddate.tm_year) + 
    "-"           + numprintf(2,0,seeddate.tm_mon+1)+
    "-"           + numprintf(2,0,seeddate.tm_mday) +
    "-dmin"       + numprintf(4,0,domainmin.x) +
    "_"           + numprintf(4,0,domainmin.y) +
    "_"           + numprintf(4,0,domainmin.z) +
    "_dmax"       + numprintf(4,0,domainmax.x) +
    "_"           + numprintf(4,0,domainmax.y) +
    "_"           + numprintf(4,0,domainmax.z) +
    "_res"        + numprintf(4,3,intergrid.x) +
    "_"           + numprintf(4,3,intergrid.y) +
    "_"           + numprintf(3,0,intergrid.z) +
    "_tau"        + numprintf(4,0,tau) +
    "_h"          + numprintf(4,3,intstep) +
    "_dmax"       + numprintf(3,0, deltamax/1000.0) +
    ".data";

  if(verbose==1){// Verbose: Save fsle values in ascii file
    cout << "Save fsle values in file: " << nfilefsle2d <<endl;
  }
  ofstream ofilefsle2d(nfilefsle2d.c_str());
  for(unsigned int q=0; q<fsle.size(); q++){
    ofilefsle2d<<fsle[q]<<endl;
  }
  ofilefsle2d.close();
  if(verbose==1){// Verbose: success file saved  
    cout << "[Complete]" << endl;
  }

 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/

  string vtkfilefsle2d = 
    "fslez_vm"   + numprintf(1,0,vfield) +
    "_date"       + numprintf(2,0,seeddate.tm_year) + 
    "-"           + numprintf(2,0,seeddate.tm_mon+1)+
    "-"           + numprintf(2,0,seeddate.tm_mday) +
    "-dmin"       + numprintf(4,0,domainmin.x) +
    "_"           + numprintf(4,0,domainmin.y) +
    "_"           + numprintf(4,0,domainmin.z) +
    "_dmax"       + numprintf(4,0,domainmax.x) +
    "_"           + numprintf(4,0,domainmax.y) +
    "_"           + numprintf(4,0,domainmax.z) +
    "_res"        + numprintf(4,3,intergrid.x) +
    "_"           + numprintf(4,3,intergrid.y) +
    "_"           + numprintf(3,0,intergrid.z) +
    "_tau"        + numprintf(4,0,tau) +
    "_h"          + numprintf(4,3,intstep) +
    "_dmax"       + numprintf(3,0, deltamax/1000.0) +
    ".vtk";

  if(verbose==1)  cout << "Save ftle field in vtk file: " << vtkfilefsle2d <<endl;

  /*ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< ni-2 <<" "<< nj-2 <<" "<< 1 <<endl;
  vtkfile<<"ORIGIN "<<grid[qcore[0]].x<<" "<<grid[qcore[0]].y<<" "<< 0.0 <<endl;
  vtkfile<<"SPACING "<<intergrid.x<<" "<<intergrid.y<<" "<< 0.0 <<endl;
  vtkfile<<"POINT_DATA "<<(ni-2)*(nj-2)<<endl;
  vtkfile<<"SCALARS fsle2d double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      vtkfile<<fsle[q]<<endl;
    }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      vtkfile<<qflag[qcore[q]]<<endl;
    }
    vtkfile.close();*/

  ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< dimgrid.i-1 <<" "<< dimgrid.j-1 <<" "<< dimgrid.k+1 <<endl;
  vtkfile<<"ORIGIN ";
  vtkfile<<grid[qcore[0]].x-(intergrid.x/2.0)<<" ";
  vtkfile<<grid[qcore[0]].y-(intergrid.y/2.0)<<" ";
  vtkfile<<grid[qcore[0]].z-(intergrid.z/2.0) <<endl;
  vtkfile<<"SPACING "<<intergrid.x<<" "<<intergrid.y<<" "<< intergrid.z <<endl;
  vtkfile<<"CELL_DATA "<<(dimgrid.i-2)*(dimgrid.j-2)*(dimgrid.k)<<endl;
  vtkfile<<"SCALARS fslez double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      vtkfile<<fsle[q]<<endl;
    }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      vtkfile<<qflag[qcore[q]]<<endl;
    }
    vtkfile.close();

  if(verbose==1)  cout << "[Complete]" <<endl;

  return 0;
}

string numprintf(int ndigits, int ndecimals, double number)
{
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
