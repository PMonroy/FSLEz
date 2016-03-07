#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#include "vectorXYZ.hpp"
#include "vectorIJK.hpp"

int MakeRegularGrid(vector<vectorXYZ> *point, vectorIJK *dimension, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax){

  double x,y,z;

  int i=0;
  int j=0;
  int k=0;

  int numpoints;

  // GRIDS POINTS
  numpoints  = (int)((domainmax.z-domainmin.z)/intergrid.z)+1;
  numpoints *= (int)((domainmax.y-domainmin.y)/intergrid.y)+1;
  numpoints *= (int)((domainmax.x-domainmin.x)/intergrid.x)+1;
 
  if(numpoints<1)
    return 1;
  
  (*point).reserve(numpoints);
  
  for(z=domainmin.z,k=0; z<=domainmax.z; z+=intergrid.z,k++){
    for(y=domainmin.y,j=0; y<=domainmax.y; y+=intergrid.y,j++){
      for(x=domainmin.x,i=0; x<=domainmax.x; x+=intergrid.x,i++){
	(*point).push_back(vectorXYZ(x,y,z));
      }
    }
  }

  // DIMENSIONS
  *dimension=vectorIJK(i,j,k);

  return 0;
}

vector<int> neighbors(vectorIJK dimension){

  int i,j,k;
  
  int nlayer=dimension.i*dimension.j;
  int ntracers=nlayer*dimension.k;
  
  int q;
  int qmin[6],qincr[6],qmax[6];
  int qneighbor;
  int dir;

  vector<int> neighbors;
  neighbors.reserve(6*ntracers);

  qincr[0]=1;
  qincr[1]=dimension.i;
  qincr[2]=nlayer;
  qincr[3]=-1;
  qincr[4]=-dimension.i;
  qincr[5]=-nlayer;

  qmin[2]=qmin[5]=0;
  qmax[2]=qmax[5]=ntracers;

  for(k=0; k<dimension.k; k++){
    qmin[1]=qmin[4]=k*nlayer;
    qmax[1]=qmax[4]=(k+1)*nlayer;

    for(j=0; j<dimension.j; j++){
      qmin[0]=qmin[3]=j*dimension.i+qmin[1];
      qmax[0]=qmax[3]=(j+1)*dimension.i+qmin[1];

      for(i=0; i<dimension.i; i++){
	q=i+j*dimension.i+k*nlayer;

	for(dir=0; dir<6; dir++){
	  qneighbor=q+qincr[dir];

	  if(qneighbor>=qmin[dir] && qneighbor<qmax[dir])
	    neighbors.push_back(qneighbor);
	  else
	    neighbors.push_back(-1);
	}
      }
    }
  }

  return neighbors;
}
