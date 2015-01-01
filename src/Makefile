include ../config.h

.NOTPARALLEL: clean install all

OBJS = amberparm.o readparm.o ambercrd.o string_manip.o NetCDFFile.o

install: all
	/bin/mv libcphmd$(SHARED_EXT) libcphmd.a $(PREFIX)/lib

clean:
	/bin/rm -f $(OBJS) libcphmd$(SHARED_EXT) libcphmd.a

all: lib

.SUFFIXES: .F90 .cpp .o

depends::
	../makedepends

.F90.o:
	$(F90) -I../include $(F90FLAGS) -c $*.F90

.cpp.o:
	$(CXX) -I../include $(CXXFLAGS) -c $*.cpp

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) $(MAKESHARED) $(OBJS) -o libcphmd$(SHARED_EXT) $(LDFLAGS)
	ar -rcs libcphmd.a $(OBJS)

include depends
