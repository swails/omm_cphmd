include ../config.h

test:: clean TopologyTest AmberParmTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libcphmd.a

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libcphmd.a

clean:
	/bin/rm -f TopologyTest AmberParmTest

depends::
	../makedepends

include depends
