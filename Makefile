include config.h

test:: unittest
	./unittest && /bin/rm -f ./unittest

unittest: unittest.cpp
	$(CXX) -o unittest unittest.cpp

include depends
