
// Testing program driving the unit tests for the functionality contained here

#include <cassert>
#include <iostream>
#include "topology.h"

using namespace std;

void check_atoms(void) {
    CpHMD::Atom my_atom(0, string("CA"), string("CX"), 6, 12.01,
                        0.3, 1.0, 0.5, 1.2, 0.85);
    CpHMD::Atom my_atom2(0, "CB", "CT", 6, 12.01, 0.4, 1.05, 0.5, 1.18, 0.80);

    // Check the first atom
    assert(my_atom.getName() == "CA");
    assert(my_atom.getType() == "CX");
    assert(my_atom.getElement() == 6);
    assert(my_atom.getMass() == 12.01);
    assert(my_atom.getCharge() == 0.3);
    assert(my_atom.getLJRadius() == 1);
    assert(my_atom.getLJEpsilon() == 0.5);
    assert(my_atom.getGBRadius() == 1.2);
    assert(my_atom.getGBScreen() == 0.85);

    // Check the second atom
    assert(my_atom2.getName() == "CB");
    assert(my_atom2.getType() == "CT");
    assert(my_atom2.getElement() == 6);
    assert(my_atom2.getMass() == 12.01);
    assert(my_atom2.getCharge() == 0.4);
    assert(my_atom2.getLJRadius() == 1.05);
    assert(my_atom2.getLJEpsilon() == 0.5);
    assert(my_atom2.getGBRadius() == 1.18);
    assert(my_atom2.getGBScreen() == 0.80);

    // Now check the atom list
    CpHMD::AtomList atoms;
    for (int i = 0; i < 10; i++) {
        atoms.push_back(CpHMD::Atom(i, "CA", "CX", 6, 12.01, 0.3,
                                    1.0, 0.5, 1.2, 0.85));
    }

    int i = 0;
    for (CpHMD::AtomList::const_iterator it = atoms.begin(); 
         it != atoms.end(); it++) {
        assert(it->getIndex() == i++);
        assert(it->getName() == "CA");
        assert(it->getType() == "CX");
        assert(it->getMass() == 12.01);
        assert(it->getCharge() == 0.3);
        assert(it->getLJRadius() == 1);
        assert(it->getLJEpsilon() == 0.5);
        assert(it->getGBRadius() == 1.2);
        assert(it->getGBScreen() == 0.85);
    }
}

void check_bonds(void) {
    CpHMD::Bond my_bond(0, 1, 100.0, 1.2);
    CpHMD::Bond my_bond2(1, 2, 150.0, 1.0);

    // Check the first bond
    assert(my_bond.getForceConstant() == 100.0);
    assert(my_bond.getEquilibriumDistance() == 1.2);
    assert(my_bond.getAtomI() == 0);
    assert(my_bond.getAtomJ() == 1);

    // Check the second bond
    assert(my_bond.getForceConstant() == 100.0);
    assert(my_bond.getEquilibriumDistance() == 1.2);
    assert(my_bond.getAtomI() == 0);
    assert(my_bond.getAtomJ() == 1);

    // Now check the bond list
    CpHMD::BondList bonds;
    for (int i = 0; i < 10; i++) {
        bonds.push_back(CpHMD::Bond(i, i+1, 100.0, 1.2));
    }

    int i = 0;
    for (CpHMD::BondList::const_iterator it = bonds.begin(); it != bonds.end();
         it++) {
        assert(it->getAtomI() == i++);
        assert(it->getAtomJ() == i);
        assert(it->getForceConstant() == 100.0);
        assert(it->getEquilibriumDistance() == 1.2);
    }
}

void check_angles(void) {
    CpHMD::Angle my_angle(0, 1, 2, 20.0, 109.5);

    // Check the angle
    assert(my_angle.getForceConstant() == 20.0);
    assert(my_angle.getEquilibriumAngle() == 109.5);
    assert(my_angle.getAtomI() == 0);
    assert(my_angle.getAtomJ() == 1);
    assert(my_angle.getAtomK() == 2);

    // Now check the angle list
    CpHMD::AngleList angles;
    for (int i = 0; i < 10; i++) {
        angles.push_back(CpHMD::Angle(i, i+1, i+2, 50.0, 109.47));
    }

    int i = 0;
    for (CpHMD::AngleList::const_iterator it = angles.begin();
         it != angles.end(); it++) {
        assert(it->getAtomI() == i++);
        assert(it->getAtomJ() == i);
        assert(it->getAtomK() == i+1);
        assert(it->getForceConstant() == 50.0);
        assert(it->getEquilibriumAngle() == 109.47);
    }
}

void check_dihedrals(void) {
    CpHMD::Dihedral my_dihedral(0, 1, 2, 3, 10.0, 180.0, 2, false);
    CpHMD::Dihedral my_dihedral2(1, 2, 3, 4, 20.0, 0.0, 4, true);

    // Check the dihedral
    assert(my_dihedral.getForceConstant() == 10);
    assert(my_dihedral.getPhase() == 180);
    assert(my_dihedral.getPeriodicity() == 2);
    assert(!my_dihedral.ignoreEndGroups());
    assert(my_dihedral.getAtomI() == 0);
    assert(my_dihedral.getAtomJ() == 1);
    assert(my_dihedral.getAtomK() == 2);
    assert(my_dihedral.getAtomL() == 3);

    // Check the second dihedral
    assert(my_dihedral2.getForceConstant() == 20);
    assert(my_dihedral2.getPhase() == 0);
    assert(my_dihedral2.getPeriodicity() == 4);
    assert(my_dihedral2.ignoreEndGroups());
    assert(my_dihedral2.getAtomI() == 1);
    assert(my_dihedral2.getAtomJ() == 2);
    assert(my_dihedral2.getAtomK() == 3);
    assert(my_dihedral2.getAtomL() == 4);

    // Now check the dihedral list
    CpHMD::DihedralList dihedrals;
    for (int i = 0; i < 10; i++) {
        dihedrals.push_back(CpHMD::Dihedral(i, i+1, i+2, i+3,
                                            10.0, 180.0, 2, false));
    }

    int i = 0;
    for (CpHMD::DihedralList::const_iterator it = dihedrals.begin();
         it != dihedrals.end(); it++) {
        assert(it->getAtomI() == i++);
        assert(it->getAtomJ() == i);
        assert(it->getAtomK() == i+1);
        assert(it->getAtomL() == i+2);
        assert(it->getForceConstant() == 10.0);
        assert(it->getPhase() == 180.0);
        assert(!it->ignoreEndGroups());
        assert(it->getPeriodicity() == 2);
    }
}

int main(int argc, char** argv) {

    cout << "Checking Atom class with list...";
    check_atoms();
    cout << " OK" << endl;

    cout << "Checking Bond class with list...";
    check_bonds();
    cout << " OK" << endl;

    cout << "Checking Angle class with list...";
    check_angles();
    cout << " OK" << endl;

    cout << "Checking Dihedral class with list...";
    check_dihedrals();
    cout << " OK" << endl;

    return 0;
}
