/// Tests the coordinate file reading/writing

#include <cassert>
#include <cmath>
#include <iostream>

#include "amber_constants.h"
#include "ambercrd.h"

using namespace std;

void check_inpcrd_crd_parsing(void) {

    Amber::AmberCoordinateFrame frame;

    frame.readRst7("crdsonly.rst7");

    assert(frame.getNatom() == 28);
    assert(frame.getPositions().size() == 28);
    assert(frame.getVelocities().size() ==  0);

    assert(abs(frame.getPositions()[0][0] - 1.8662825) < 1e-5);
    assert(abs(frame.getPositions()[0][1] - 1.2455949) < 1e-5);
    assert(abs(frame.getPositions()[0][2] - 0.6672110) < 1e-5);
}

void check_inpcrd_crdvel_parsing(void) {

    Amber::AmberCoordinateFrame frame;

    frame.readRst7("trx.inpcrd");

    assert(frame.getNatom() == 1654);
    assert(frame.getPositions().size() == 1654);
    assert(frame.getVelocities().size() == 1654);

    assert(abs(frame.getPositions()[0][0] - -12.6649955) < 1e-7);
    assert(abs(frame.getPositions()[0][1] - 7.5329535) < 1e-7);
    assert(abs(frame.getPositions()[0][2] - -3.3805141) < 1e-7);

    assert(abs(frame.getVelocities()[1653][0] - -0.1679970 * Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[1653][1] - -0.0931307 * Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[1653][2] - 0.1423713 * Amber::AMBER_TIME_PER_PS) < 1e-5);
}

void check_inpcrd_crdbox_parsing(void) {

    Amber::AmberCoordinateFrame frame;

    frame.readRst7("crds_box.rst7");

    assert(frame.getNatom() == 18660);
    assert(frame.getPositions().size() == 18660);
    assert(frame.getVelocities().size() == 0);

    assert(abs(frame.getPositions()[0][0] - 42.4969453) < 1e-5);
    assert(abs(frame.getPositions()[0][1] - 26.5038927) < 1e-5);
    assert(abs(frame.getPositions()[0][2] - 24.4848153) < 1e-5);

    assert(abs(frame.getBoxA() - 65.3721047) < 1e-5);
    assert(abs(frame.getBoxB() - 65.3721047) < 1e-5);
    assert(abs(frame.getBoxC() - 65.3721047) < 1e-5);
    assert(abs(frame.getBoxAlpha() - 109.4712190) < 1e-5);
    assert(abs(frame.getBoxBeta() - 109.4712190) < 1e-5);
    assert(abs(frame.getBoxGamma() - 109.4712190) < 1e-5);
}

void check_inpcrd_crdvelbox_parsing(void) {

    Amber::AmberCoordinateFrame frame;

    frame.readRst7("crds_vels_box.rst7");

    assert(frame.getNatom() == 2101);
    assert(frame.getPositions().size() == 2101);
    assert(frame.getVelocities().size() == 2101);

    assert(abs(frame.getPositions()[0][0] - 6.6528795) < 1e-5);
    assert(abs(frame.getPositions()[0][1] - 6.6711416) < 1e-5);
    assert(abs(frame.getPositions()[0][2] - -8.5963255) < 1e-5);

    assert(abs(frame.getVelocities()[2100][0] - 0.3091332*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[2100][1] - 0.7355925*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs(frame.getVelocities()[2100][2] - (-0.1206961)*Amber::AMBER_TIME_PER_PS) < 1e-5);

    assert(abs(frame.getBoxA() - 30.2642725) < 1e-5);
    assert(abs(frame.getBoxB() - 30.2642725) < 1e-5);
    assert(abs(frame.getBoxC() - 30.2642725) < 1e-5);
    assert(abs(frame.getBoxAlpha() - 109.4712190) < 1e-5);
    assert(abs(frame.getBoxBeta() - 109.4712190) < 1e-5);
    assert(abs(frame.getBoxGamma() - 109.4712190) < 1e-5);
}

int main() {

    cout << "Testing amber inpcrd file reading with just coordinates...";
    check_inpcrd_crd_parsing();
    cout << " OK." << endl;

    cout << "Testing amber inpcrd file reading with coordinates and velocities...";
    check_inpcrd_crdvel_parsing();
    cout << " OK." << endl;

    cout << "Testing amber inpcrd file reading with coordinates and box...";
    check_inpcrd_crdbox_parsing();
    cout << " OK." << endl;

    cout << "Testing amber inpcrd file reading with coordinates, velocities, and box...";
    check_inpcrd_crdvelbox_parsing();
    cout << " OK." << endl;

}
