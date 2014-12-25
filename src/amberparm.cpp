// amberparm.cpp -- contains the functionality

#include "amberparm.h"
#include "exceptions.h"

using namespace std;
using namespace CpHMD;

AmberParm::AmberParm(string const& filename) {
    rdparm(filename);
}

AmberParm::AmberParm(const char* filename) {
    rdparm(filename);
}

/// Implement the add**** methods
void AmberParm::addAtom(Atom& new_atom) {
    if (new_atom.getIndex() != (int) atoms_.size())
        throw AmberParmError("Atoms must be added sequentially!");
    atoms_.push_back(new_atom);
}

void AmberParm::addAtom(std::string const& name, std::string const& type,
                        int element, double mass, double charge, double lj_rad,
                        double lj_eps, double gb_rad, double gb_screen) {
    int i = (int) atoms_.size();
    atoms_.push_back(Atom(i, name, type, element, mass, charge, lj_rad,
                          lj_eps, gb_rad, gb_screen));
}

void AmberParm::addAtom(const char* name, const char* type,
                        int element, double mass, double charge, double lj_rad,
                        double lj_eps, double gb_rad, double gb_screen) {
    int i = (int) atoms_.size();
    atoms_.push_back(Atom(i, name, type, element, mass, charge, lj_rad,
                          lj_eps, gb_rad, gb_screen));
}

// addBonds
void AmberParm::addBond(Bond& new_bond) {
    int i = new_bond.getAtomI();
    int j = new_bond.getAtomJ();
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom)
        throw AmberParmError("Bond atom index out of range");
    bonds_.push_back(new_bond);
}

void AmberParm::addBond(int i, int j, double kf, double req) {
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom)
        throw AmberParmError("Bond atom index out of range");
    bonds_.push_back(Bond(i, j, kf, req));
}

// addAngles
void AmberParm::addAngle(Angle& new_angle) {
    int i = new_angle.getAtomI();
    int j = new_angle.getAtomJ();
    int k = new_angle.getAtomK();
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom || k < 0 || k >= natom)
        throw AmberParmError("Angle atom index out of range");
    angles_.push_back(new_angle);
}

void AmberParm::addAngle(int i, int j, int k, double kf, double theteq) {
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom || k < 0 || k >= natom)
        throw AmberParmError("Angle atom index out of range");
    angles_.push_back(Angle(i, j, k, kf, theteq));
}

// addDihedrals
void AmberParm::addDihedral(Dihedral& new_dihedral) {
    int i = new_dihedral.getAtomI();
    int j = new_dihedral.getAtomJ();
    int k = new_dihedral.getAtomK();
    int l = new_dihedral.getAtomL();
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom || k < 0 || k >= natom ||
        l < 0 || l >= natom)
        throw AmberParmError("Dihedral atom index out of range");
    dihedrals_.push_back(new_dihedral);
}

void AmberParm::addDihedral(int i, int j, int k, int l, double kf, double phase,
                            int periodicity, bool ignore_end) {
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom || k < 0 || k >= natom ||
        l < 0 || l >= natom)
        throw AmberParmError("Dihedral atom index out of range");
    dihedrals_.push_back(Dihedral(i, j, k, l, kf, phase,
                                  periodicity, ignore_end));
}

void AmberParm::rdparm(string const& filename) {
    throw AmberParmError("This feature is not yet implemented");
}

void AmberParm::rdparm(const char* filename) {
    rdparm(string(filename));
}
