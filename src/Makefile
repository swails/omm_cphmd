include ../config.h

.NOTPARALLEL: clean install all

OBJS = amberparm.o readparm.o

install: all
	/bin/mv libcphmd$(SHARED_EXT) libcphmd.a $(PREFIX)/lib

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
	$(CXX) $(MAKESHARED) $(OBJS) -o libcphmd$(SHARED_EXT)
	ar -rcs libcphmd.a $(OBJS)
