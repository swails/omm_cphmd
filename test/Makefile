include ../config.h

test:: clean TopologyTest AmberParmTest OpenMMTest CoordinateFileTest \
       NetCDFCoordinateFileTest NetCDFFileTest UnitCellTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest
	./OpenMMTest && /bin/rm ./OpenMMTest serialized_system.xml
	./CoordinateFileTest && /bin/rm ./CoordinateFileTest test.rst7
	./NetCDFCoordinateFileTest && /bin/rm ./NetCDFCoordinateFileTest
	./NetCDFFileTest && /bin/rm ./NetCDFFileTest
	./UnitCellTest && /bin/rm ./UnitCellTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libamber.a $(LDFLAGS)

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libamber.a $(LDFLAGS)

OpenMMTest: OpenMMTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o OpenMMTest OpenMMTest.cpp ../lib/libamber.a $(LDFLAGS)

CoordinateFileTest: CoordinateFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o CoordinateFileTest CoordinateFileTest.cpp ../lib/libamber.a $(LDFLAGS)

NetCDFCoordinateFileTest: NetCDFCoordinateFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o NetCDFCoordinateFileTest NetCDFCoordinateFileTest.cpp ../lib/libamber.a $(LDFLAGS)

NetCDFFileTest: NetCDFFileTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o NetCDFFileTest NetCDFFileTest.cpp ../lib/libamber.a $(LDFLAGS)

UnitCellTest: UnitCellTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o UnitCellTest UnitCellTest.cpp ../lib/libamber.a $(LDFLAGS)

clean:
	/bin/rm -f TopologyTest AmberParmTest OpenMMTest NetCDFCoordinateFileTest
	/bin/rm -f CoordinateFileTest NetCDFFileTest UnitCellTest

depends::
	../makedepends

include depends
