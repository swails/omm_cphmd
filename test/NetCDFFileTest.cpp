// NetCDFFileTest.cpp -- tests the NetCDFFile class
#include <cassert>
#include <cmath>
#include <iostream>

#include "Amber.h"
#include "OpenMM.h"

#define ASSERT_RAISES(statement, exctype) \
    try { statement; assert(false);} \
    catch (exctype &e) {assert(true);}

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

    OpenMM::Vec3 abc = testFile.getCellLengths();

    assert(abs(abc[0] - 30.2642725) < 1e-6);
    assert(abs(abc[1] - 30.2642725) < 1e-6);
    assert(abs(abc[2] - 30.2642725) < 1e-6);

    OpenMM::Vec3 ang = testFile.getCellAngles();

    assert(abs(ang[0] - 109.471219) < 1e-6);
    assert(abs(ang[1] - 109.471219) < 1e-6);
    assert(abs(ang[2] - 109.471219) < 1e-6);

    delete coords;
    delete vels;
}

void test_bad_ncrst_read() {
    {
        Amber::AmberNetCDFFile file(Amber::AmberNetCDFFile::TRAJECTORY);
        ASSERT_RAISES(file.readFile("files/amber.ncrst"), Amber::AmberCrdError)
    }
    {
        Amber::AmberNetCDFFile file;
        file.readFile("files/amber.ncrst");
        ASSERT_RAISES(file.getCoordinates(1), Amber::AmberCrdError)
        ASSERT_RAISES(file.getCoordinates(-1), Amber::AmberCrdError)
        ASSERT_RAISES(file.getVelocities(1), Amber::AmberCrdError)
        ASSERT_RAISES(file.getVelocities(-1), Amber::AmberCrdError)
        ASSERT_RAISES(file.getForces(), Amber::AmberCrdError)
    }
    {
        Amber::AmberNetCDFFile file(Amber::AmberNetCDFFile::RESTART);
        ASSERT_RAISES(file.readFile("files/trx.prmtop"), Amber::NotNetcdf)
    }
    {
        Amber::AmberNetCDFFile file(Amber::AmberNetCDFFile::TRAJECTORY);
        ASSERT_RAISES(file.readFile("files/amber.ncrst"), Amber::AmberCrdError)
    }
}

void test_nctraj_read(void) {

    Amber::AmberNetCDFFile traj(Amber::AmberNetCDFFile::TRAJECTORY);

    traj.readFile("files/crdvelfrc.nc");

    assert(traj.getNatom() == 12288);
    assert(traj.getNumFrames() == 5);
    assert(traj.hasCoordinates());
    assert(traj.hasVelocities());
    assert(traj.hasForces());
    assert(traj.hasBox());
    assert(!traj.hasREMD());

    vector<OpenMM::Vec3> *crd;
    vector<OpenMM::Vec3> *vel;
    vector<OpenMM::Vec3> *frc;

    crd = traj.getCoordinates(0);
    vel = traj.getVelocities(0);
    frc = traj.getForces(0);

    assert(abs((*crd)[0][0] - 0.70450461) < 1e-5);
    assert(abs((*crd)[0][1] - 4.95110083) < 1e-5);
    assert(abs((*crd)[0][2] - 3.75570297) < 1e-5);
    assert(abs((*crd)[12287][0] - 48.6208725) < 1e-5);
    assert(abs((*crd)[12287][1] - 45.05132675) < 1e-5);
    assert(abs((*crd)[12287][2] - 2.33579612) < 1e-5);

    assert(abs((*vel)[0][0] - 0.12262645*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[0][1] - -0.09410848*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[0][2] - 0.18947887*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][0] - 0.63712376*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][1] - -0.17904748*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][2] - -0.01668696*Amber::AMBER_TIME_PER_PS) < 1e-5);

    assert(abs((*frc)[0][0] - 27.66975403) < 1e-5);
    assert(abs((*frc)[0][1] - 19.00367165) < 1e-5);
    assert(abs((*frc)[0][2] - 2.24443197) < 1e-5);
    assert(abs((*frc)[12287][0] - -4.94580793) < 1e-5);
    assert(abs((*frc)[12287][1] - -9.5912199) < 1e-5);
    assert(abs((*frc)[12287][2] - -5.21790218) < 1e-5);

    delete crd;
    delete vel;
    delete frc;

    crd = traj.getCoordinates(4);
    vel = traj.getVelocities(4);
    frc = traj.getForces(4);

    assert(abs((*crd)[0][0] - 0.71575785) < 1e-5);
    assert(abs((*crd)[0][1] - 4.94342136) < 1e-5);
    assert(abs((*crd)[0][2] - 3.77085066) < 1e-5);
    assert(abs((*crd)[12287][0] - 48.67422104) < 1e-5);
    assert(abs((*crd)[12287][1] - 45.03504562) < 1e-5);
    assert(abs((*crd)[12287][2] - 2.34541345) < 1e-5);

    assert(abs((*vel)[0][0] - 0.1462021*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[0][1] - -0.09275416*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[0][2] - 0.18184757*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][0] - 0.65629011*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][1] - -0.21025696*Amber::AMBER_TIME_PER_PS) < 1e-5);
    assert(abs((*vel)[12287][2] - 0.19651674*Amber::AMBER_TIME_PER_PS) < 1e-5);

    assert(abs((*frc)[0][0] - 27.13381577) < 1e-5);
    assert(abs((*frc)[0][1] - 20.83939171) < 1e-5);
    assert(abs((*frc)[0][2] - 1.19384193) < 1e-5);
    assert(abs((*frc)[12287][0] - -4.78174591) < 1e-5);
    assert(abs((*frc)[12287][1] - -10.30693913) < 1e-5);
    assert(abs((*frc)[12287][2] - -5.04377699) < 1e-5);

    delete crd;
    delete vel;
    delete frc;

    // Check some errors
    ASSERT_RAISES(traj.getCoordinates(5), Amber::AmberCrdError)
    ASSERT_RAISES(traj.getCoordinates(-1), Amber::AmberCrdError)
    ASSERT_RAISES(traj.getVelocities(5), Amber::AmberCrdError)
    ASSERT_RAISES(traj.getVelocities(-1), Amber::AmberCrdError)
    ASSERT_RAISES(traj.getForces(5), Amber::AmberCrdError)
    ASSERT_RAISES(traj.getForces(-1), Amber::AmberCrdError)
}

void test_ncrst_write(void) {
    Amber::AmberNetCDFFile testFile(Amber::AmberNetCDFFile::RESTART);

    testFile.writeFile("files/tmp12345.nc", 10, true, true, false, true,
                       false, 0, "test dummy file", "NetCDFFileTest");

    vector<OpenMM::Vec3> positions, velocities;
    for (int i = 0; i < 10; i++) {
        double f = (double) i;
        positions.push_back(OpenMM::Vec3(f, f+1, f+2));
        velocities.push_back(OpenMM::Vec3(f*10, f*10+10, f*10+20));
    }
    testFile.setCoordinates(positions);
    testFile.setVelocities(velocities);
    testFile.setUnitCell(OpenMM::Vec3(1.0, 0.0, 0.0),
                         OpenMM::Vec3(0.0, 1.0, 0.0),
                         OpenMM::Vec3(0.0, 0.0, 1.0));
    testFile.setTime(0.0);
}

int main() {
    cout << "Testing NetCDF reading of NetCDF restart file...";
    test_ncrst_read();
    cout << " OK." << endl;

    cout << "Testing error checking of NetCDF restart file...";
    test_bad_ncrst_read();
    cout << " OK." << endl;

    cout << "Testing NetCDF reading of trajectory with coordinates, forces, and velocities...";
    test_nctraj_read();
    cout << " OK." << endl;

    cout << "Testing NetCDF restart file writing...";
    test_ncrst_write();
    cout << " OK." << endl;
    return 0;
}
