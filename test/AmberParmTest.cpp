
// Testing program driving the unit tests for the AmberParm class

#include <cassert>
#include <cmath>
#include <iostream>

#include "amberparm.h"
#include "exceptions.h"

using namespace std;

void check_add_atoms(void) {
    Amber::AmberParm parm;

    Amber::Atom atom(0, "N", "N3", 7, 14.01, -1.0, 0.95,
                     0.1, 1.3, 0.85);
    parm.addAtom(atom);
    parm.addAtom("H1", "H", 1, 1.008, 1.0, 0.5, 0.1, 0.8, 0.85);

    assert(parm.Atoms().size() == 2);
}

void check_bad_add_atoms(void) {
    Amber::AmberParm parm;

    bool caught = false;
    try {
        Amber::Atom atom(1, "N", "N3", 7, 14.01, -1.0, 0.95,
                         0.1, 1.3, 0.85);
        parm.addAtom(atom);
    } catch (Amber::AmberParmError &e) {
        caught = true;
    }
    assert(caught);
}

void check_add_bonds(void) {
    Amber::AmberParm parm;

    parm.addAtom("N", "N3", 7, 14.01, -1.0, 0.95, 0.1, 1.3, 0.85);
    parm.addAtom("H1", "H", 1, 1.008, 0.5, 0.5, 0.05, 0.8, 0.85);
    parm.addAtom("H2", "H", 1, 1.008, 0.5, 0.5, 0.05, 0.8, 0.85);

    Amber::Bond bond(0, 1, 500.0, 0.8);
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
    Amber::AmberParm parm;

    bool caught = false;
    try {
        parm.addBond(0, 1, 500.0, 1.0);
    } catch (Amber::AmberParmError &e) {
        caught = true;
    }

    assert(caught);
}

void check_add_angles(void) {
    Amber::AmberParm parm;

    Amber::Angle angle(0, 1, 2, 50.0, 109.47);

    parm.addAtom("N", "N", 7, 14.01, -1.0, 0.8, 0.1, 1.2, 0.85);
    parm.addAtom("CA", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CB", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CG", "CX", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);

    parm.addAngle(angle);
    parm.addAngle(1, 2, 3, 60.0, 120.0);

    assert(parm.Atoms().size() == 4);
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

void check_bad_add_angles(void) {
    Amber::AmberParm parm;

    bool caught = false;
    try {
        parm.addAngle(0, 1, 2, 500.0, 1.0);
    } catch (Amber::AmberParmError &e) {
        caught = true;
    }

    assert(caught);
}

void check_add_dihedrals(void) {
    Amber::AmberParm parm;

    Amber::Dihedral dihed(0, 1, 2, 3, 50.0, 180.0, 2, 1.2, 2.0, true);
    Amber::Dihedral dihed2(0, 1, 2, 3, 20.0, 0.0, 3, 1.2, 2.0, false);

    parm.addAtom("N", "N", 7, 14.01, -1.0, 0.8, 0.1, 1.2, 0.85);
    parm.addAtom("CA", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CB", "CT", 6, 12.01, 0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("CG", "CX", 6, 12.01, -0.5, 0.9, 0.15, 1.3, 0.85);
    parm.addAtom("HG", "CX", 1, 1.01, 0.5, 0.6, 0.05, 0.8, 0.80);

    parm.addDihedral(dihed);
    parm.addDihedral(dihed2);
    parm.addDihedral(1, 2, 3, 4, 40.0, 0.0, 1, 1.2, 2.0, false);

    assert(parm.Atoms().size() == 5);
    assert(parm.Bonds().size() == 0);
    assert(parm.Angles().size() == 0);
    assert(parm.Dihedrals().size() == 3);

    assert(parm.Dihedrals()[0].getAtomI() == 0);
    assert(parm.Dihedrals()[0].getAtomJ() == 1);
    assert(parm.Dihedrals()[0].getAtomK() == 2);
    assert(parm.Dihedrals()[0].getAtomL() == 3);
    assert(parm.Dihedrals()[1].getAtomI() == 0);
    assert(parm.Dihedrals()[1].getAtomJ() == 1);
    assert(parm.Dihedrals()[1].getAtomK() == 2);
    assert(parm.Dihedrals()[1].getAtomL() == 3);
    assert(parm.Dihedrals()[2].getAtomI() == 1);
    assert(parm.Dihedrals()[2].getAtomJ() == 2);
    assert(parm.Dihedrals()[2].getAtomK() == 3);
    assert(parm.Dihedrals()[2].getAtomL() == 4);

    assert(parm.Dihedrals()[0].getForceConstant() == 50.0);
    assert(parm.Dihedrals()[1].getForceConstant() == 20.0);
    assert(parm.Dihedrals()[2].getForceConstant() == 40.0);
    assert(parm.Dihedrals()[0].getPhase() == 180);
    assert(parm.Dihedrals()[1].getPhase() == 0);
    assert(parm.Dihedrals()[2].getPhase() == 0);
    assert(parm.Dihedrals()[0].getPeriodicity() == 2);
    assert(parm.Dihedrals()[1].getPeriodicity() == 3);
    assert(parm.Dihedrals()[2].getPeriodicity() == 1);
}

void check_bad_add_dihedrals(void) {
    Amber::AmberParm parm;

    bool caught = false;
    try {
        parm.addDihedral(0, 1, 2, 3, 50.0, 180.0, 2, 1.2, 2.0, false);
    } catch (Amber::AmberParmError &e) {
        caught = true;
    }

    assert(caught);
}

void check_rdparm(void) {
    Amber::AmberParm parm;

    parm.rdparm("trx.prmtop");

    // Check atoms and atom properties

    assert(parm.Atoms().size() == 1654);
    assert(parm.Atoms()[0].getName() == "N");
    assert(parm.Atoms()[0].getType() == "N3");
    assert(parm.Atoms()[0].getMass() == 14.01);
    assert(abs(parm.Atoms()[0].getCharge() - 0.1849) < 1e-4);
    assert(parm.Atoms()[0].getGBRadius() == 1.55);
    assert(parm.Atoms()[0].getGBScreen() == 0.79);
    assert(abs(parm.Atoms()[0].getLJRadius() - 1.824) < 1e-4);
    assert(abs(parm.Atoms()[0].getLJEpsilon() - 0.17) < 1e-4);

    assert(parm.Atoms()[1653].getName() == "OXT");
    assert(parm.Atoms()[1653].getType() == "O2");
    assert(parm.Atoms()[1653].getMass() == 16.00);
    assert(abs(parm.Atoms()[1653].getCharge() - -0.8055) < 1e-4);
    assert(parm.Atoms()[1653].getGBRadius() == 1.5);
    assert(parm.Atoms()[1653].getGBScreen() == 0.85);
    assert(abs(parm.Atoms()[1653].getLJRadius() - 1.6612) < 1e-4);
    assert(abs(parm.Atoms()[1653].getLJEpsilon() - 0.21) < 1e-4);

    assert(parm.Bonds().size() == 1670);
    assert(parm.Bonds()[0].getAtomI() == 9);
    assert(parm.Bonds()[0].getAtomJ() == 10);
    assert(parm.Bonds()[0].getForceConstant() == 553);
    assert(parm.Bonds()[0].getEquilibriumDistance() == 0.96);

    assert(parm.Angles().size() == 3049);
    assert(parm.Angles()[0].getAtomI() == 11);
    assert(parm.Angles()[0].getAtomJ() == 13);
    assert(parm.Angles()[0].getAtomK() == 14);
    assert(parm.Angles()[0].getForceConstant() == 30);
    assert(abs(parm.Angles()[0].getEquilibriumAngle() - 120) < 5e-4);

    assert(parm.Dihedrals().size() == 5402);
    assert(parm.Dihedrals()[0].getAtomI() == 12);
    assert(parm.Dihedrals()[0].getAtomJ() == 11);
    assert(parm.Dihedrals()[0].getAtomK() == 13);
    assert(parm.Dihedrals()[0].getAtomL() == 14);
    assert(!parm.Dihedrals()[0].ignoreEndGroups());
    assert(parm.Dihedrals()[0].getForceConstant() == 2.0);
    assert(parm.Dihedrals()[0].getPhase() == 0);
    assert(parm.Dihedrals()[0].getPeriodicity() == 1);
    assert(parm.Dihedrals()[0].getScee() == 1.2);
    assert(parm.Dihedrals()[0].getScnb() == 2.0);

    assert(parm.Dihedrals()[1].getAtomI() == 12);
    assert(parm.Dihedrals()[1].getAtomJ() == 11);
    assert(parm.Dihedrals()[1].getAtomK() == 13);
    assert(parm.Dihedrals()[1].getAtomL() == 14);
    assert(parm.Dihedrals()[1].ignoreEndGroups());
    assert(parm.Dihedrals()[1].getForceConstant() == 2.5);
    assert(abs(parm.Dihedrals()[1].getPhase() - 180) < 1e-4);
    assert(parm.Dihedrals()[1].getPeriodicity() == 2);
    assert(parm.Dihedrals()[1].getScee() == 1.2);
    assert(parm.Dihedrals()[1].getScnb() == 2.0);

    // Check the residue properties
    assert(parm.ResidueLabels()[0] == "SER");
    assert(parm.ResidueLabels()[107] == "ALA");
    assert(parm.ResiduePointers()[1] - parm.ResiduePointers()[0] == 13);
    assert(parm.ResiduePointers()[0] == 0);
    assert(parm.ResiduePointers()[1] == 13);
    assert(parm.ResiduePointers()[108] - parm.ResiduePointers()[107] == 11);

    // Check the exclusions
    for (Amber::AmberParm::bond_iterator it = parm.BondBegin();
            it != parm.BondEnd(); it++) {
        assert(parm.isExcluded(it->getAtomI(), it->getAtomJ()));
        assert(parm.isExcluded(it->getAtomJ(), it->getAtomI()));
    }
    for (Amber::AmberParm::angle_iterator it = parm.AngleBegin();
            it != parm.AngleEnd(); it++) {
        assert(parm.isExcluded(it->getAtomI(), it->getAtomJ()));
        assert(parm.isExcluded(it->getAtomJ(), it->getAtomI()));
        assert(parm.isExcluded(it->getAtomI(), it->getAtomK()));
        assert(parm.isExcluded(it->getAtomK(), it->getAtomI()));
        assert(parm.isExcluded(it->getAtomJ(), it->getAtomK()));
        assert(parm.isExcluded(it->getAtomK(), it->getAtomJ()));
    }
    /* Make sure the exceptions (i.e., 1-4s that are NOT excluded) are NOT part
     * of the exceptions.
     */
    for (Amber::AmberParm::dihedral_iterator it = parm.DihedralBegin();
            it != parm.DihedralEnd(); it++) {
        bool ignore_end = it->ignoreEndGroups();
        if (!ignore_end) {
            assert(!parm.isExcluded(it->getAtomI(), it->getAtomL()));
            assert(!parm.isExcluded(it->getAtomL(), it->getAtomI()));
        }
    }

    assert(!parm.isPeriodic());
    assert(parm.IfBox() == 0);
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

    cout << "Checking error catching in adding angles to AmberParm...";
    check_bad_add_angles();
    cout << " OK" << endl;

    cout << "Checking adding dihedrals to AmberParm...";
    check_add_dihedrals();
    cout << " OK" << endl;

    cout << "Checking error catching in adding dihedrals to AmberParm...";
    check_bad_add_dihedrals();
    cout << " OK" << endl;

    cout << "Checking Amber topology file reading...";
    check_rdparm();
    cout << " OK" << endl;

    return 0;
}
