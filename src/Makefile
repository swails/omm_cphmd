include ../config.h

.NOTPARALLEL: clean install all

OBJS = amberparm.o readparm.o ambercrd.o string_manip.o NetCDFFile.o gbmodels.o \
	   unitcell.o

install: all
	/bin/mv libamber$(SHARED_EXT) libamber.a $(PREFIX)/lib

clean:
	/bin/rm -f $(OBJS) libamber$(SHARED_EXT) libamber.a

all: lib

.SUFFIXES: .F90 .cpp .o

depends::
	../makedepends

.F90.o:
	$(F90) -I../include $(F90FLAGS) -c $*.F90

.cpp.o:
	$(CXX) -I../include $(CXXFLAGS) -c $*.cpp

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) $(MAKESHARED) $(OBJS) -o libamber$(SHARED_EXT) $(LDFLAGS)
	ar -rcs libamber.a $(OBJS)

include depends
