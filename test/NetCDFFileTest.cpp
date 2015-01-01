// NetCDFFileTest.cpp -- tests the NetCDFFile class
#include <cassert>
#include <cmath>
#include <iostream>

#include "Amber.h"

using namespace std;

void test_ncrst_read(void) {

    Amber::AmberNetCDFFile testFile(Amber::AmberNetCDFFile::RESTART);

    testFile.readFile("files/amber.ncrst");

    assert(testFile.getApplication() == "AMBER");
    assert(testFile.getProgram() == "sander");
    assert(testFile.getProgramVersion() == "11.0");
    assert(testFile.getTitle() == "ACE");

    assert(testFile.hasCoordinates());
    assert(testFile.hasVelocities());
    assert(!testFile.hasForces());
    assert(testFile.hasBox());
    assert(!testFile.hasREMD());
    assert(testFile.getREMDDimension() == 0);
    assert(testFile.getNumFrames() == 1);
    assert(testFile.getNatom() == 2101);

    vector<OpenMM::Vec3> *coords = testFile.getCoordinates();

    assert(coords->size() == 2101);
    assert(abs((*coords)[0][0] - 6.82122493) < 1e-5);
    assert(abs((*coords)[0][1] - 6.62762507) < 1e-5);
    assert(abs((*coords)[0][2] - -8.51669) < 1e-5);
    assert(abs((*coords)[2100][0] - 3.64159598) < 1e-5);
    assert(abs((*coords)[2100][1] - 8.62796977) < 1e-5);
    assert(abs((*coords)[2100][2] - -8.56491885) < 1e-5);

    vector<OpenMM::Vec3> *vels = testFile.getVelocities();

    assert(vels->size() == 2101);
    assert(abs((*vels)[0][0] - -0.14058448*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vels)[0][1] - -0.14768776*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vels)[0][2] - 0.18278307*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vels)[2100][0] - 0.3497022*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vels)[2100][1] - 0.39152533*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vels)[2100][2] - 0.41794168*Amber::AMBER_TIME_PER_PS) < 1e-5);

    delete coords;
    delete vels;
}

int main() {
    cout << "Testing NetCDF reading of NetCDF restart file...";
    test_ncrst_read();
    cout << " OK." << endl;

    return 0;
}
