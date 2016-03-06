#CC_MONO=g++ -ggdb -Wall -mcmodel=medium #TO DEBUG
CC_MONO=g++ -O2 -Wall -mcmodel=medium #TO OPTIMIZE

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: fslez

readparameters.o: readparameters.cpp readparameters.hpp
	$(CC) -c readparameters.cpp 

vectorXYZ.o: vectorXYZ.cpp vectorXYZ.hpp
	$(CC) -c vectorXYZ.cpp 

vectorIJK.o: vectorIJK.cpp vectorIJK.hpp
	$(CC) -c vectorIJK.cpp 

constants.o: constants.cpp constants.hpp
	$(CC) -c constants.cpp 

gridconstruction.o: gridconstruction.cpp gridconstruction.hpp vectorXYZ.o vectorIJK.o
	$(CC) -c gridconstruction.cpp

integration.o: integration.cpp integration.hpp vectorXYZ.o constants.o
	$(CC) -c integration.cpp

vflow.o: vflow.cpp vflow.hpp vectorXYZ.o constants.o
	$(CC) -c vflow.cpp -lnetcdf_c++

fslez.o: fslez.cpp
	$(CC) -c fslez.cpp 

fslez: fslez.o readparameters.o gridconstruction.o vflow.o constants.o vectorXYZ.o vectorIJK.o integration.o
	$(CC)  fslez.o readparameters.o gridconstruction.o vflow.o integration.o constants.o vectorXYZ.o vectorIJK.o -o fslez -lnetcdf_c++
clean:
	$(RM) *.o 
