// amberparm.cpp -- contains the functionality

#include "amberparm.h"
#include "exceptions.h"

#include <cmath>
#include <iostream>

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

    // Data structures for parsing data from the prmtop
    ParmDataMap parmData;
    ParmStringMap parmComments, unkParmData;
    ParmFormatMap parmFormats;
    vector<string> flagList;
    string version;

    ExitStatus retval = readparm(filename, flagList, parmData, parmComments,
                                 unkParmData, parmFormats, version);

    if (retval == NOOPEN) {
        string msg = "Could not open " + filename + " for reading";
        throw AmberParmError(msg.c_str());
    }

    if (retval == NOVERSION) {
        string msg = "No %VERSION in " + filename + ". Bad parm format";
        throw AmberParmError(msg.c_str());
    }

    if (retval == EMPTY) {
        string msg = filename + " was empty";
        throw AmberParmError(msg.c_str());
    }

    if (retval == ERR) {
        string msg = "Prmtop parsing error parsing " + filename;
        throw AmberParmError(msg.c_str());
    }

    // If we got here, we must have been successful in our parsing

    // First add our atoms
    string atomic_number = "ATOMIC_NUMBER";
    string atom_name = "ATOM_NAME";
    string atom_type = "AMBER_ATOM_TYPE";
    string radii = "RADII";
    string screen = "SCREEN";
    string acoef = "LENNARD_JONES_ACOEF";
    string bcoef = "LENNARD_JONES_BCOEF";
    string ljtype = "ATOM_TYPE_INDEX";
    string nbidx = "NONBONDED_PARM_INDEX";
    string pointers = "POINTERS";
    string mass = "MASS";
    string charge = "CHARGE";

    if (parmData.count(pointers) < 1)
        throw AmberParmError("Missing POINTERS in prmtop");
    int N = parmData[pointers][NATOM].i;

    if (parmData.count(atomic_number) < 1 || parmData[atomic_number].size() != N)
        throw AmberParmError(
                "Missing ATOMIC_NUMBER in prmtop or wrong # of elements");
    if (parmData.count(atom_name) < 1 || parmData[atom_name].size() != N)
        throw AmberParmError("Missing  ATOM_NAME in prmtop");
    if (parmData.count(atom_type) < 1 || parmData[atom_type].size() != N)
        throw AmberParmError(
                "Missing AMBER_ATOM_TYPE in prmtop or wrong # of elements");
    if (parmData.count(radii) < 1 || parmData[radii].size() != N)
        throw AmberParmError("Missing RADII in prmtop or wrong # of elements");
    if (parmData.count(screen) < 1 || parmData[screen].size() != N)
        throw AmberParmError("Missing SCREEN in prmtop or wrong # of elements");
    if (parmData.count(acoef) < 1)
        throw AmberParmError("Missing LENNARD_JONES_ACOEF in prmtop");
    if (parmData.count(bcoef) < 1)
        throw AmberParmError("Missing LENNARD_JONES_BCOEF in prmtop");
    if (parmData.count(ljtype) < 1 || parmData[ljtype].size() != N)
        throw AmberParmError(
                "Missing ATOM_TYPE_INDEX in prmtop or wrong # of elements");
    if (parmData.count(nbidx) < 1)
        throw AmberParmError("Missing NONBONDED_PARM_INDEX in prmtop");
    if (parmData.count(mass) < 1 || parmData[mass].size() != N)
        throw AmberParmError("Missing MASS in prmtop or wrong # of elements");
    if (parmData.count(charge) < 1 || parmData[charge].size() != N)
        throw AmberParmError("Missing CHARGE in prmtop or wrong # of elements");

    // Now extract the per-particle LJ rmin/epsilon parameters
    vector<double> lj_rmin, lj_eps;
    int ntypes = parmData[pointers][NTYPES].i;
    const double ONE_SIXTH = 1.0 / 6.0;
    for (int i = 0; i < ntypes; i++) {
        int lj_index = parmData[nbidx][ntypes*i+i].i - 1;
        if (parmData[acoef][lj_index].f < 1.0e-10) {
            lj_rmin.push_back(0.5);
            lj_eps.push_back(0);
        } else {
            double factor = 2 * parmData[acoef][lj_index].f /
                                parmData[bcoef][lj_index].f;
            lj_rmin.push_back(pow(factor, ONE_SIXTH) * 0.5);
            lj_eps.push_back(parmData[bcoef][lj_index].f * 0.5 / factor);
        }
    }

    // Now add the atoms
    for (int i = 0; i < parmData[pointers][NATOM].i; i++) {
        int typ = parmData[ljtype][i].i;
        double chg = parmData[charge][i].f / 18.2223;
        addAtom(parmData[atom_name][i].c, parmData[atom_type][i].c,
                parmData[atomic_number][i].i, parmData[mass][i].f,
                chg, lj_rmin[typ-1], lj_eps[typ-1],
                parmData[radii][i].f, parmData[screen][i].f);
    }

    // Now add the bonds
    string bonds = "BONDS_WITHOUT_HYDROGEN";
    string bondsh = "BONDS_INC_HYDROGEN";
    string bondk = "BOND_FORCE_CONSTANT";
    string bondeq = "BOND_EQUIL_VALUE";
    int nbonh = parmData[pointers][NBONH].i;
    int mbona = parmData[pointers][MBONA].i;
    int numbnd = parmData[pointers][NUMBND].i;

    if (parmData.count(bonds) < 1 || parmData[bonds].size() != mbona*3)
        throw AmberParmError("Bad (or missing) BONDS_WITHOUT_HYDROGEN section");
    if (parmData.count(bondsh) < 1 || parmData[bondsh].size() != nbonh*3)
        throw AmberParmError("Bad (or missing) BONDS_INC_HYDROGEN section");
    if (parmData.count(bondk) < 1 || parmData[bondk].size() != numbnd)
        throw AmberParmError("Bad (or missing) BOND_FORCE_CONSTANT section");
    if (parmData.count(bondeq) < 1 || parmData[bondeq].size() != numbnd)
        throw AmberParmError("Bad (or missing) BOND_EQUIL_VALUE section");
    for (int i = 0; i < nbonh; i++) {
        int i3 = i * 3;
        int ii = parmData[bondsh][i3  ].i / 3;
        int jj = parmData[bondsh][i3+1].i / 3;
        int bi = parmData[bondsh][i3+2].i - 1;
        addBond(ii, jj, parmData[bondk][bi].f, parmData[bondeq][bi].f);
    }
    for (int i = 0; i < mbona; i++) {
        int i3 = i * 3;
        int ii = parmData[bonds][i3  ].i / 3;
        int jj = parmData[bonds][i3+1].i / 3;
        int bi = parmData[bonds][i3+2].i - 1;
        addBond(ii, jj, parmData[bondk][bi].f, parmData[bondeq][bi].f);
    }

    // Now add the angles
    string angles = "ANGLES_WITHOUT_HYDROGEN";
    string anglesh = "ANGLES_INC_HYDROGEN";
    string anglek = "ANGLE_FORCE_CONSTANT";
    string angleeq = "ANGLE_EQUIL_VALUE";
    int ntheth = parmData[pointers][NTHETH].i;
    int mtheta = parmData[pointers][MTHETA].i;
    int numang = parmData[pointers][NUMANG].i;

    if (parmData.count(angles) < 1 || parmData[angles].size() != mtheta*4)
        throw AmberParmError("Bad (or missing) ANGLES_WITHOUT_HYDROGEN section");
    if (parmData.count(anglesh) < 1 || parmData[anglesh].size() != ntheth*4)
        throw AmberParmError("Bad (or missing) ANGLES_INC_HYDROGEN section");
    if (parmData.count(anglek) < 1 || parmData[anglek].size() != numang)
        throw AmberParmError("Bad (or missing) ANGLE_FORCE_CONSTANT section");
    if (parmData.count(angleeq) < 1 || parmData[angleeq].size() != numang)
        throw AmberParmError("Bad (or missing) ANGLE_EQUIL_VALUE section");

    for (int i = 0; i < ntheth; i++) {
        int i4 = i * 4;
        int ii = parmData[anglesh][i4  ].i / 3;
        int jj = parmData[anglesh][i4+1].i / 3;
        int kk = parmData[anglesh][i4+2].i / 3;
        int ai = parmData[anglesh][i4+3].i - 1;
        double ang = parmData[angleeq][ai].f * 180.0 / M_PI;
        addAngle(ii, jj, kk, parmData[anglek][ai].f, ang);
    }
    for (int i = 0; i < mtheta; i++) {
        int i4 = i * 4;
        int ii = parmData[angles][i4  ].i / 3;
        int jj = parmData[angles][i4+1].i / 3;
        int kk = parmData[angles][i4+2].i / 3;
        int ai = parmData[angles][i4+3].i - 1;
        double ang = parmData[angleeq][ai].f * 180.0 / M_PI;
        addAngle(ii, jj, kk, parmData[anglek][ai].f, ang);
    }
}

void AmberParm::rdparm(const char* filename) {
    rdparm(string(filename));
}
