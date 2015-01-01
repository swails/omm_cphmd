/** NetCDFCoordinateFileTest -- This file contains tests for the NetCDF
  * coordinate file parsing and writing routines
  */
#include <cassert>
#include <cmath>
#include <iostream>

#include "Amber.h"

using namespace std;

void check_ncrst_read(void) {

    Amber::AmberCoordinateFrame frame;

    frame.readRst7("files/amber.ncrst");

    assert(frame.getNatom() == 2101);
    assert(frame.getPositions().size() == 2101);
    assert(abs(frame.getPositions()[0][0] - 6.82122493) < 1e-5);
    assert(abs(frame.getPositions()[0][1] - 6.62762507) < 1e-5);
    assert(abs(frame.getPositions()[0][2] - -8.51669) < 1e-5);
    assert(abs(frame.getPositions()[2100][0] - 3.64159598) < 1e-5);
    assert(abs(frame.getPositions()[2100][1] - 8.62796977) < 1e-5);
    assert(abs(frame.getPositions()[2100][2] - -8.56491885) < 1e-5);

    assert(frame.getVelocities().size() == 2101);
    assert(abs(frame.getVelocities()[0][0] - -0.14058448*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[0][1] - -0.14768776*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[0][2] - 0.18278307*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[2100][0] - 0.3497022*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[2100][1] - 0.39152533*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[2100][2] - 0.41794168*Amber::AMBER_TIME_PER_PS) < 1e-5);

    assert(abs(frame.getBoxA() - 30.2642725) < 1e-6);
    assert(abs(frame.getBoxB() - 30.2642725) < 1e-6);
    assert(abs(frame.getBoxC() - 30.2642725) < 1e-6);
    assert(abs(frame.getBoxAlpha() - 109.471219) < 1e-6);
    assert(abs(frame.getBoxBeta() - 109.471219) < 1e-6);
    assert(abs(frame.getBoxGamma() - 109.471219) < 1e-6);
}

int main() {

    cout << "Checking NetCDF restart file reading...";
    check_ncrst_read();
    cout << " OK." << endl;
}
