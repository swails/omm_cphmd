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
    delete system;
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
    delete system;
}

void check_omm_gb(string const& model, double cutoff,
                  double saltcon, double nonbe) {
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
    OpenMM::NonbondedForce::NonbondedMethod meth;
    if (cutoff <= 0)
        meth = OpenMM::NonbondedForce::NoCutoff;
    else
        meth = OpenMM::NonbondedForce::CutoffNonPeriodic;
    system = parm.createSystem(
            meth, cutoff, string("None"), false,
            model, 0.0, saltcon);
    OpenMM::VerletIntegrator integrator(0.002);
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
    assert(abs(1 - e/nonbe) < 5e-7); // relative error...
    delete system;
}

// So we can pass a const char*
void check_omm_gb(const char* model, double cutoff,
                  double saltcon, double nonbe) {
    check_omm_gb(string(model), cutoff, saltcon, nonbe);
}

int main() {

    // Load the main plugins
    OpenMM::Platform::loadPluginsFromDirectory(
            OpenMM::Platform::getDefaultPluginsDirectory());

    // Run the tests
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
    check_omm_gb("HCT", 0.0, 0.0, -4380.6377735);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB HCT w/ 0.1M salt...";
    check_omm_gb("HCT", 0.0, 0.1, -4383.2215249);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC1 energy...";
    check_omm_gb("OBC1", 0.0, 0.0, -4430.6048991);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC1 w/ 0.1M salt...";
    check_omm_gb("OBC1", 0.0, 0.1, -4433.1897402);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC2 energy...";
    check_omm_gb("OBC2", 0.0, 0.0, -4317.4276516);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC2 w/ 0.1M salt...";
    check_omm_gb("OBC2", 0.0, 0.1, -4319.9948287);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn energy...";
    check_omm_gb("GBn", 0.0, 0.0, -4252.4065109);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn w/ 0.1M salt...";
    check_omm_gb("GBn", 0.0, 0.1, -4254.9660314);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn2 energy...";
    check_omm_gb("GBn2", 0.0, 0.0, -4324.7676537);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn2 w/ 0.1M salt...";
    check_omm_gb("GBn2", 0.0, 0.1, -4327.3449966);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB HCT energy (15A cutoff)...";
    check_omm_gb("HCT", 15.0, 0.0, -4541.2393356);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB HCT w/ 0.1M salt (15A cutoff)...";
    check_omm_gb("HCT", 15.0, 0.1, -4548.0488375);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC1 energy (15A cutoff)...";
    check_omm_gb("OBC1", 15.0, 0.0, -4593.9490639);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC1 w/ 0.1M salt (15A cutoff)...";
    check_omm_gb("OBC1", 15.0, 0.1, -4600.9552586);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC2 energy (15A cutoff)...";
    check_omm_gb("OBC2", 15.0, 0.0, -4481.5135247);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB OBC2 w/ 0.1M salt (15A cutoff)...";
    check_omm_gb("OBC2", 15.0, 0.1, -4488.8354190);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn energy (15A cutoff)...";
    check_omm_gb("GBn", 15.0, 0.0, -4410.0631502);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn w/ 0.1M salt (15A cutoff)...";
    check_omm_gb("GBn", 15.0, 0.1, -4417.7292305);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn2 energy (15A cutoff)...";
    check_omm_gb("GBn2", 15.0, 0.0, -4480.8521321);
    cout << " OK." << endl;

    cout << "Testing OpenMM GB GBn2 w/ 0.1M salt (15A cutoff)...";
    check_omm_gb("GBn2", 15.0, 0.1, -4488.0561661);
    cout << " OK." << endl;
}
