include ../config.h

test:: clean TopologyTest AmberParmTest OpenMMTest CoordinateFileTest \
       NetCDFCoordinateFileTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest
	./OpenMMTest && /bin/rm ./OpenMMTest serialized_system.xml
	./CoordinateFileTest && /bin/rm ./CoordinateFileTest test.rst7
	./NetCDFCoordinateFileTest && /bin/rm ./NetCDFCoordinateFileTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libcphmd.a $(LDFLAGS)

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libcphmd.a $(LDFLAGS)

OpenMMTest: OpenMMTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o OpenMMTest OpenMMTest.cpp ../lib/libcphmd.a $(LDFLAGS)

CoordinateFileTest: CoordinateFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o CoordinateFileTest CoordinateFileTest.cpp ../lib/libcphmd.a $(LDFLAGS)

NetCDFCoordinateFileTest: NetCDFCoordinateFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o NetCDFCoordinateFileTest NetCDFCoordinateFileTest.cpp ../lib/libcphmd.a $(LDFLAGS)

clean:
	/bin/rm -f TopologyTest AmberParmTest OpenMMTest NetCDFCoordinateFileTest

depends::
	../makedepends

include depends
