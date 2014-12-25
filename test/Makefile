include ../config.h

test:: TopologyTest
	./TopologyTest && /bin/rm ./TopologyTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libcphmd.a

depends::
	../makedepends

include depends
