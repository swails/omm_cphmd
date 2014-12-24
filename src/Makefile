include ../config.h

.NOTPARALLEL: clean install all

OBJS = amberparm.o readparm.o

all: lib

.SUFFIXES: .F90 .cpp .o

depends::
	../makedepends

.F90.o:
	$(F90) -I../include $(F90FLAGS) -c $*.F90

.cpp.o:
	$(CXX) -I../include $(CXXFLAGS) -c $*.cpp

clean:
	/bin/rm -f $(OBJS)

lib: $(OBJS)
	$(CXX) -shared $(OBJS) -o libommcphmd.so
