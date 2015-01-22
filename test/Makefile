include ../config.h

test:: clean TopologyTest AmberParmTest OpenMMTest CoordinateFileTest \
       NetCDFCoordinateFileTest NetCDFFileTest UnitCellTest
	./TopologyTest && /bin/rm ./TopologyTest
	./AmberParmTest && /bin/rm ./AmberParmTest
	./OpenMMTest && /bin/rm ./OpenMMTest
	$(MAKE) XmlDiff
	./CoordinateFileTest && /bin/rm ./CoordinateFileTest test.rst7
	./NetCDFCoordinateFileTest && /bin/rm ./NetCDFCoordinateFileTest
	./NetCDFFileTest && /bin/rm -f ./NetCDFFileTest files/tmp12345.nc
	./UnitCellTest && /bin/rm ./UnitCellTest

TopologyTest: TopologyTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o TopologyTest TopologyTest.cpp ../lib/libamber.a $(LDFLAGS)

AmberParmTest: AmberParmTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o AmberParmTest AmberParmTest.cpp ../lib/libamber.a $(LDFLAGS)

OpenMMTest: OpenMMTest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o OpenMMTest OpenMMTest.cpp ../lib/libamber.a $(LDFLAGS)

# Dependencies for this are wrong -- make sure it is run after OpenMMTest
XmlDiff:
	@if [ -z `/usr/bin/diff -w serialized_system.xml files/serialized_system.xml` ]; then \
		echo "XML File comparison... OK" && /bin/rm serialized_system.xml; \
	else \
		echo "XML File comparison... FAILED" && false; \
	fi

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
