include ../config.h

test:: unittest
	./unittest && /bin/rm ./unittest

unittest: unittest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o unittest unittest.cpp

depends::
	../makedepends

include depends
