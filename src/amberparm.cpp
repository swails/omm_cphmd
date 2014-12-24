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
    atoms_.push_back(new_atom);
}

void AmberParm::addAtom(int i, std::string const& name, std::string const& type,
                        int element, double mass, double charge, double lj_rad,
                        double lj_eps, double gb_rad, double gb_screen) {
    atoms_.push_back(Atom(i, name, type, element, mass, charge, lj_rad,
                          lj_eps, gb_rad, gb_screen));
}

void AmberParm::addAtom(int i, const char* name, const char* type,
                        int element, double mass, double charge, double lj_rad,
                        double lj_eps, double gb_rad, double gb_screen) {
    atoms_.push_back(Atom(i, name, type, element, mass, charge, lj_rad,
                          lj_eps, gb_rad, gb_screen));
}

// addBonds
void AmberParm::addBond(Bond& new_bond) {
    bonds_.push_back(new_bond);
}

void AmberParm::addBond(int i, int j, double kf, double req) {
    bonds_.push_back(Bond(i, j, kf, req));
}

// addAngles
void AmberParm::addAngle(Angle& new_angle) {
    angles_.push_back(new_angle);
}

void AmberParm::addAngle(int i, int j, int k, double kf, double theteq) {
    angles_.push_back(Angle(i, j, k, kf, theteq));
}

// addDihedrals
void AmberParm::addDihedral(Dihedral& new_dihedral) {
    dihedrals_.push_back(new_dihedral);
}

void AmberParm::addDihedral(int i, int j, int k, int l, double kf, double phase,
                            int periodicity, bool ignore_end) {
    dihedrals_.push_back(Dihedral(i, j, k, l, kf, phase,
                                  periodicity, ignore_end));
}

void AmberParm::rdparm(string const& filename) {
    throw AmberParmError("This feature is not yet implemented");
}

void AmberParm::rdparm(const char* filename) {
    rdparm(string(filename));
}
