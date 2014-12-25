
// Testing program driving the unit tests for the AmberParm class

#include <cassert>
#include <iostream>

#include "amberparm.h"
#include "exceptions.h"

using namespace std;

void check_add_atoms(void) {
    CpHMD::AmberParm parm;

    CpHMD::Atom atom(0, "N", "N3", 7, 14.01, -1.0, 0.95,
                     0.1, 1.3, 0.85);
    parm.addAtom(atom);
    parm.addAtom("H1", "H", 1, 1.008, 1.0, 0.5, 0.1, 0.8, 0.85);

    assert(parm.Atoms().size() == 2);
}

void check_bad_add_atoms(void) {
    CpHMD::AmberParm parm;

    bool caught = false;
    try {
        CpHMD::Atom atom(1, "N", "N3", 7, 14.01, -1.0, 0.95,
                         0.1, 1.3, 0.85);
        parm.addAtom(atom);
    } catch (CpHMD::AmberParmError &e) {
        caught = true;
    }
    assert(caught);
}

void check_add_bonds(void) {
    CpHMD::AmberParm parm;

    parm.addAtom("N", "N3", 7, 14.01, -1.0, 0.95, 0.1, 1.3, 0.85);
    parm.addAtom("H1", "H", 1, 1.008, 0.5, 0.5, 0.05, 0.8, 0.85);
    parm.addAtom("H2", "H", 1, 1.008, 0.5, 0.5, 0.05, 0.8, 0.85);

    CpHMD::Bond bond(0, 1, 500.0, 0.8);
    parm.addBond(bond);
    parm.addBond(1, 2, 500.0, 0.8);

    assert(parm.Atoms().size() == 3);
    assert(parm.Bonds().size() == 2);
    assert(parm.Bonds()[0].getForceConstant() == 500);
    assert(parm.Bonds()[1].getForceConstant() == 500);
    assert(parm.Bonds()[0].getEquilibriumDistance() == 0.8);
    assert(parm.Bonds()[1].getEquilibriumDistance() == 0.8);

    assert(parm.Bonds()[0].getAtomI() == 0);
    assert(parm.Bonds()[0].getAtomJ() == 1);
    assert(parm.Bonds()[1].getAtomI() == 1);
    assert(parm.Bonds()[1].getAtomJ() == 2);
}

void check_bad_add_bonds(void) {
    CpHMD::AmberParm parm;

    bool caught = false;
    try {
        parm.addBond(0, 1, 500.0, 1.0);
    } catch (CpHMD::AmberParmError &e) {
        caught = true;
    }

    assert(caught);
}

void check_add_angles(void) {
    CpHMD::AmberParm parm;

    CpHMD::Angle angle(0, 1, 2, 50.0, 109.47);

    parm.addAtom("N", "N", 7, 14.01, -1.0, 0.8, 0.1, 1.2, 0.85);
    parm.addAtom("CA", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CB", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CG", "CX", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);

    parm.addAngle(angle);
    parm.addAngle(1, 2, 3, 60.0, 120.0);

    assert(parm.Atoms().size() == 4);
    assert(parm.Angles().size() == 2);
    assert(parm.Bonds().size() == 0);
    assert(parm.Angles().size() == 2);

    assert(parm.Angles()[0].getAtomI() == 0);
    assert(parm.Angles()[0].getAtomJ() == 1);
    assert(parm.Angles()[0].getAtomK() == 2);
    assert(parm.Angles()[1].getAtomI() == 1);
    assert(parm.Angles()[1].getAtomJ() == 2);
    assert(parm.Angles()[1].getAtomK() == 3);

    assert(parm.Angles()[0].getForceConstant() == 50.0);
    assert(parm.Angles()[1].getForceConstant() == 60.0);
    assert(parm.Angles()[0].getEquilibriumAngle() == 109.47);
    assert(parm.Angles()[1].getEquilibriumAngle() == 120.0);
}

int main() {

    cout << "Checking adding atoms to AmberParm...";
    check_add_atoms();
    cout << " OK" << endl;

    cout << "Checking error catching in adding atoms to AmberParm...";
    check_bad_add_atoms();
    cout << " OK" << endl;

    cout << "Checking adding bonds to AmberParm...";
    check_add_bonds();
    cout << " OK" << endl;

    cout << "Checking error catching in adding bonds to AmberParm...";
    check_bad_add_bonds();
    cout << " OK" << endl;

    cout << "Checking adding angles to AmberParm...";
    check_add_angles();
    cout << " OK" << endl;

    return 0;
}
