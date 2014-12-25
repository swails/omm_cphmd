include ../config.h

.NOTPARALLEL: clean install all

OBJS = amberparm.o readparm.o

install: all
	/bin/mv libcphmd.so libcphmd.a ../lib

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
	$(CXX) -shared $(OBJS) -o libcphmd.so
	ar -rcs libcphmd.a $(OBJS)
