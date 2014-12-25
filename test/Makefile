include ../config.h

test:: unittest
	./unittest && /bin/rm ./unittest

unittest: unittest.cpp
	$(CXX) $(CXXFLAGS) -I../include -o unittest unittest.cpp ../lib/libcphmd.a

depends::
	../makedepends

include depends
