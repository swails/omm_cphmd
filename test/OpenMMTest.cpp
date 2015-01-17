// This tests the OpenMM functionality here
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "Amber.h"
#include "OpenMM.h"

using namespace std;

void check_omm_system(void) {
    Amber::AmberParm parm;

    parm.rdparm("files/trx.prmtop");

    OpenMM::System *system = parm.createSystem();

    ofstream output;
    output.open("serialized_system.xml");
    OpenMM::XmlSerializer::serialize<OpenMM::System>(system, string("System"),
                                                     output);
    output.close();
    delete system;
}

void check_gas_energy(void) {
    Amber::AmberParm parm;
    Amber::AmberCoordinateFrame frame;

    parm.rdparm("files/trx.prmtop");
    frame.readRst7("files/trx.inpcrd");

    OpenMM::System *system = 0;

    // Get the coordinates in nanometers
    vector<OpenMM::Vec3> positions = frame.getPositions();
    for (size_t i = 0; i < positions.size(); i++)
        positions[i] *= Amber::NANOMETER_PER_ANGSTROM;

    // Create the system, integrator, and context with the Reference platform
    system = parm.createSystem();
    OpenMM::VerletIntegrator integrator(0.002);
    OpenMM::Context *context = new OpenMM::Context(*system, integrator, 
                OpenMM::Platform::getPlatformByName(string("Reference")));

    // Set the starting coordinates
    context->setPositions(positions);

    // Now compare the energies to Amber energies

    // Now get the bond energy
    OpenMM::State s = context->getState(OpenMM::State::Energy, false,
                                        1<<Amber::AmberParm::BOND_FORCE_GROUP);
    double e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 631.8992707) < 1e-5);
    // Now get the angle energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::ANGLE_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 898.2542659) < 1e-5);
    // Now get the dihedral energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::DIHEDRAL_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 566.4454191) < 1e-5);
    // Now get the nonbonded energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::NONBONDED_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(1 - e/-2313.5905426) < 1e-6); // relative error...
}

void check_omm_pme(void) {
    Amber::AmberParm parm;
    Amber::AmberCoordinateFrame frame;

    parm.rdparm("files/4096wat.parm7");
    frame.readRst7("files/4096wat.rst7");

    OpenMM::System *system = 0;

    // Get the coordinates in nanometers
    vector<OpenMM::Vec3> positions = frame.getPositions();
    for (size_t i = 0; i < positions.size(); i++)
        positions[i] *= Amber::NANOMETER_PER_ANGSTROM;

    // Create the system, integrator, and context with the Reference platform
    system = parm.createSystem(OpenMM::NonbondedForce::PME, 8.0);
    OpenMM::VerletIntegrator integrator(0.002);
    OpenMM::Context *context = new OpenMM::Context(*system, integrator, 
                OpenMM::Platform::getPlatformByName(string("Reference")));

    // Set up a unit cell from our coordinate file
    Amber::UnitCell cell(frame.getBoxA(), frame.getBoxB(), frame.getBoxC(),
                         frame.getBoxAlpha(), frame.getBoxBeta(),
                         frame.getBoxGamma());

    // Set the starting coordinates and box
    context->setPositions(positions);
    context->setPeriodicBoxVectors(cell.getVectorA()/10,
                                   cell.getVectorB()/10,
                                   cell.getVectorC()/10);

    // Now compare the energies to Amber energies

    // Now get the bond energy
    OpenMM::State s = context->getState(OpenMM::State::Energy, false,
                                        1<<Amber::AmberParm::BOND_FORCE_GROUP);
    double e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 0) < 2e-8);
    // Now get the angle energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::ANGLE_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(e == 0);
    // Now get the dihedral energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::DIHEDRAL_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(e == 0);
    // Now get the nonbonded energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::NONBONDED_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(1 - e/-39344.3259106) < 2e-4); // relative error...
}

void check_omm_gb1(void) {
    Amber::AmberParm parm;
    Amber::AmberCoordinateFrame frame;

    parm.rdparm("files/trx.prmtop");
    frame.readRst7("files/trx.inpcrd");

    OpenMM::System *system = 0;

    // Get the coordinates in nanometers
    vector<OpenMM::Vec3> positions = frame.getPositions();
    for (size_t i = 0; i < positions.size(); i++)
        positions[i] *= Amber::NANOMETER_PER_ANGSTROM;

    // Create the system, integrator, and context with the Reference platform
    system = parm.createSystem(
            OpenMM::NonbondedForce::NoCutoff, 0.0, string("None"), false,
            string("HCT"));
    OpenMM::VerletIntegrator integrator(0.002);
//  cout << "I have " << OpenMM::Platform::getNumPlatforms() << " plats" << endl;
    OpenMM::Context *context = new OpenMM::Context(*system, integrator, 
                OpenMM::Platform::getPlatformByName(string("CPU")));

    // Set the starting coordinates
    context->setPositions(positions);

    // Now compare the energies to Amber energies

    // Now get the bond energy
    OpenMM::State s = context->getState(OpenMM::State::Energy, false,
                                        1<<Amber::AmberParm::BOND_FORCE_GROUP);
    double e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 631.8992707) < 1e-5);
    // Now get the angle energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::ANGLE_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 898.2542659) < 1e-5);
    // Now get the dihedral energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::DIHEDRAL_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(e - 566.4454191) < 1e-5);
    // Now get the nonbonded energy
    s = context->getState(OpenMM::State::Energy, false,
                          1<<Amber::AmberParm::NONBONDED_FORCE_GROUP);
    e = s.getPotentialEnergy() * Amber::CALORIE_PER_JOULE;
    assert(abs(1 - e/-4383.2214985) < 1e-6); // relative error...
}

int main() {

    cout << "Testing OpenMM System creation and serialization...";
    check_omm_system();
    cout << " OK." << endl;

    cout << "Testing OpenMM gas phase energy...";
    check_gas_energy();
    cout << " OK." << endl;

    cout << "Testing OpenMM PME potential energies...";
    check_omm_pme();
    cout << " OK." << endl;

    cout << "Testing OpenMM GB HCT energy...";
    check_omm_gb1();
    cout << " OK." << endl;
}
