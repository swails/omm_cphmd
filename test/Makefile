include ../config.h

test:: TopologyTest AmberParmTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libcphmd.a

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libcphmd.a

depends::
	../makedepends

include depends
