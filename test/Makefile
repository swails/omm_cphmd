include ../config.h

test:: clean TopologyTest AmberParmTest OpenMMTest CoordinateFileTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest
	./OpenMMTest && /bin/rm ./OpenMMTest serialized_system.xml
	./CoordinateFileTest && /bin/rm ./CoordinateFileTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libcphmd.a $(LDFLAGS)

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libcphmd.a $(LDFLAGS)

OpenMMTest: OpenMMTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o OpenMMTest OpenMMTest.cpp ../lib/libcphmd.a $(LDFLAGS)

CoordinateFileTest: CoordinateFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o CoordinateFileTest CoordinateFileTest.cpp ../lib/libcphmd.a $(LDFLAGS)

clean:
	/bin/rm -f TopologyTest AmberParmTest OpenMMTest

depends::
	../makedepends

include depends
