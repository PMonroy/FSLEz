#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <string>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <ctime>

using namespace std;
#include "vectorXYZ.hpp"

//GLOBALS PARAMETERS
int verbose=0; //Verbosity is disabled by default
int vfield;
vectorXYZ domainmin;
vectorXYZ intergrid;
vectorXYZ domainmax;
struct tm seeddate= {0};                                                
/* The structure tm contains nine members of type int, which are:
   * tm_sec       int     seconds after the minute 0-61
   * tm_min       int     minutes after the hour   0-59
   * tm_hour      int     hours since midnight     0-23
   * tm_mday      int     day of the month         1-31
   * tm_mon       int     months since January     0-11
   * tm_year      int     years since 1900
   * tm_wday      int     days since Sunday        0-6
   * tm_yday      int     days since January 1     0-365
   * tm_isdst     int     Daylight Saving Time flag     
 */    
double  intstep;
int tau;
double deltamax;

//FUNCTIONS
int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams){ 
  string me="GetcmdlineParameters";
  int opt;

  int verboseflag=0; // Verbose flag option
  int fnameparamsflag=0; // File parameters flag option

  while ((opt = getopt(narg, cmdarg, "f:vh")) != -1) 
    {
      switch(opt)
	{

	case 'f':
	  fnameparamsflag++;
	  if (optarg)
	    *fnameparams = optarg;
	  else
	    return 1;
	  break;

	case 'v':
	  verboseflag++;
	  verbose=1;
	  break;

	case 'h':
	  cout<<"Usage: "<< cmdarg[0] <<" [OPTIONS]"<<endl;
	  cout<<" -f [file]      Where [file] is the input file parmeters (mandatory)" <<endl;
	  cout<<" -v             Verbose mode enable"<<endl;
	  cout<<" -h             Print this help and exit (optional)"<<endl;
         return 1;

	case '?':
	  cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
	  return 1;
	}
    }

  /* Check mandatory or repeated options*/

  if(fnameparamsflag==0){//Check if file option is set  
    cout << "Error("<<me <<"): Option -f <file> is mandatory" <<endl;
    cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
    return 1;
  }
  else if (fnameparamsflag>1){//Check if file option is repeated
    cout << "Warning("<<me <<"): Option -f is repeated. Parameters file now is the last one ("<<*fnameparams<<")"<<endl;
  }

  if(verboseflag>1){//Check if verbose option is repeated
      cout << "Warning("<<me <<"): Option -v is repeated. Verbose mode continues to be enable. "<<endl;
    }

  return 0;
}
int GetfileParameters(string nfileparameters){

  string me="readparams()";
  ifstream fparameters(nfileparameters.c_str());

  enum  enum_parameters { VFIELD,
			  DOMAINMIN,
			  INTERGRID,
			  DOMAINMAX,
			  SEEDDATE,
			  INTSTEP,
			  TAU,
			  DELTAMAX,
			  NPARAMETERS
  };
  string pname[NPARAMETERS];

  pname[VFIELD]="vfield";
  pname[DOMAINMIN]="domainmin";
  pname[INTERGRID]="intergrid";
  pname[DOMAINMAX]="domainmax";
  pname[SEEDDATE]="seeddate";
  pname[INTSTEP]="intstep";
  pname[TAU]="tau";
  pname[DELTAMAX]="deltamax";

  int delimiter,end;
  string name, value;
  stringstream svalue;
  string line;

  int *pflag;
  int parameter;
  
  if (!fparameters.is_open())
    {
      cout << me <<": Skipping unreadable file \"" << nfileparameters.c_str() << "\" "<<endl; 
      return 1;
    }

  pflag = (int*) calloc (NPARAMETERS,sizeof(int));

  while(!fparameters.eof()){
    getline(fparameters,line);
    
    if (line[0] == '#') 
      continue;  /* ignore comment line which starts with #*/
    
    if (isspace(line[0])) 
      continue; /* ignore blank line */
    
    delimiter = line.find("=");
    end = line.length();
    
    if (end==0) 
      continue; /* ignore blank line */
    
    value = line.substr(delimiter+1, end);
    name = line.substr(0, delimiter);
    
    for(parameter=0; parameter<NPARAMETERS; parameter++){
      if (name.compare(pname[parameter]) == 0){
	pflag[parameter]++;
	if(pflag[parameter]>1){
	  cout  << me << ": Parameter "<< name << " repeated"<<endl;
	  return 1;
	}
	if ((delimiter+1) == end){
	  cout << me << ": Parameter "<< name << " has no value"<<endl;
	  return 1;
	}
	switch (parameter) {
	case VFIELD:
	  vfield = atoi(value.c_str());
	  break;
	case DOMAINMIN:
	  svalue<<value;
	  if((svalue>>domainmin)==0){
	    cout  << me << ": Format of domainmin is incorrect" << endl;
	    return 1;
	  }
	  svalue.clear();//clear any bits set
	  svalue.str(string());
	  break;
	case INTERGRID:
	  svalue << value;		  
	  if((svalue>>intergrid)==0){
	    cout  << me << ": Format of domainmin is incorrect" << endl;
	    return 1;
	  }
	  svalue.clear();//clear any bits set
	  svalue.str(string());
	  break;
	case DOMAINMAX:
	  svalue << value;		  
	  if((svalue>>domainmax)==0){
	    cout  << me << ": Format of domainmin is incorrect" << endl;
	    return 1;
	  }
	  svalue.clear();//clear any bits set
	  svalue.str(string());
	  break;
	case SEEDDATE:                                           
	  if(sscanf(value.c_str(),"%d-%d-%d",&seeddate.tm_mday,&seeddate.tm_mon,&seeddate.tm_year)!=3){
	    cout  << me << ": Date format in start_date is incorrect" << endl;
	    return 1;                                                
	  }                    
	  seeddate.tm_mon = seeddate.tm_mon - 1;
	  seeddate.tm_hour= 11;
	  seeddate.tm_min = 0;
	  seeddate.tm_sec = 0;
	  break; 
	case INTSTEP:
	  intstep = atof(value.c_str());
	  break;
	case TAU:
	  tau = atoi(value.c_str());
	  break;  
	case DELTAMAX:
	  deltamax = atof(value.c_str());
	  break;  
	default:
	  cout  << me << ": Unknown parameter "<< name <<endl;
	  return 1;
	}
	break;
      }
    }
    if(parameter==NPARAMETERS){
      cout  << me << ": Unknown parameter "<< name <<endl ;
      return 1;
    }
  }

  for(parameter=0; parameter<NPARAMETERS; parameter++){
    if(pflag[parameter]==0){
      cout  << me << ": Parameter "<< pname[parameter] << " is not defined"<<endl;
      return 1;
    }
  }  
  fparameters.close();
  free(pflag);

  return 0;
}
