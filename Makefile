CPP=g++
CPPFLAGS=-O2 -Wall

LDFLAGS= -lnetcdf_c++ 
RM=rm -rf

DESTDIR=$(HOME)/bin

common_src=  constants.cpp vflow.cpp vectorXYZ.cpp vectorIJK.cpp rparameters.cpp gridconstruction.cpp constants.cpp integration.cpp
common_obj=$(common_src:.cpp=.o) 
common_dep=$(common_obj:.o=.d)  # one dependency file for each source

#RTIME 
rtime_src= rtime.cpp
rtime_obj=$(rtime_src:.cpp=.o) 
rtime_dep=$(rtime_obj:.o=.d)  # one dependency file for each source

#FSLEZ
fslez_src= fslez.cpp
fslez_obj=$(fslez_src:.cpp=.o) 
fslez_dep=$(fslez_obj:.o=.d)  # one dependency file for each source

.PHONY: all rtime fslez

all: rtime fslez

rtime: $(common_obj) $(rtime_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

fslez: $(common_obj) $(fslez_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

-include $(common_dep) # include all dep files in makefile
-include $(rtime_dep) 
-include $(fselz_dep) 


%.d: %.cpp 	# rule to generate a dep file by using the g++ prepocesor
	$(CPP) $(CPPFLAGS) -MM -MT $(@:.d=.o) $< -MF $@

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

.PHONY: debug rtime fslez
debug: CPPFLAGS+= -DDEBUG -ggdb # debug with gdb
debug: rtime fslez

.PHONY: clean
clean:
	$(RM) $(common_obj) $(rtime_obj) $(fslez_obj) *.d *~ *# rtime

.PHONY: install
install: rtime fslez
	install $^ $(DESTDIR)

# Trick to run comfortably in emacs (emulating compile command)
#ARGS="--default"
#run:	# For include other args use -> make run ARGS="-arg1 -arg2=foo" 
#	rtime $(ARGS)
