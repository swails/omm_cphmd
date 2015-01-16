/* amberparm.cpp -- contains the functionality for reading an Amber topology
 * file and using that to instantiate an OpenMM System
 */

#include "amber/amber_constants.h"
#include "amber/amberparm.h"
#include "amber/exceptions.h"
#include "amber/gbmodels.h"
#include "amber/unitcell.h"

#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;
using namespace Amber;

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
                            int periodicity, double scee, double scnb,
                            bool ignore_end) {
    int natom = (int)atoms_.size();
    if (i < 0 || i >= natom || j < 0 || j >= natom || k < 0 || k >= natom ||
        l < 0 || l >= natom)
        throw AmberParmError("Dihedral atom index out of range");
    dihedrals_.push_back(Dihedral(i, j, k, l, kf, phase, periodicity,
                                  scee, scnb, ignore_end));
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
    string resptr = "RESIDUE_POINTER";
    string reslab = "RESIDUE_LABEL";

    if (parmData.count(pointers) < 1)
        throw AmberParmError("Missing POINTERS in prmtop");
    int N = parmData[pointers][NATOM].i;
    int NR = parmData[pointers][NRES].i;
    ifbox_ = parmData[pointers][IFBOX].i;

    // Allocate space for the exclusion and exception lists
    exclusion_list_.reserve(N);
    vector<set<int> > exception_list(N);
    for (int i = 0; i < N; i++) {
        set<int> s1;
        set<int> s2;
        exclusion_list_.push_back(s1);
        exception_list.push_back(s2);
    }

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
    if (parmData.count(reslab) < 1 || parmData[reslab].size() != NR)
        throw AmberParmError(
                "Missing RESIDUE_LABEL in prmtop or wrong # of elements");
    if (parmData.count(resptr) < 1 || parmData[resptr].size() != NR)
        throw AmberParmError(
                "Missing RESIDUE_POINTER in prmtop or wrong # of elements");

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

    // Now add the residues
    for (int i = 0; i < parmData[pointers][NRES].i; i++) {
        residue_pointers_.push_back(parmData[resptr][i].i-1);
        residue_labels_.push_back(string(parmData[reslab][i].c));
    }
    // Make it so the trick of subtracting pointers[n+1]-pointers[n] works even
    // for the last residue
    residue_pointers_.push_back(N);

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

    // Now add the dihedrals
    string dihedrals = "DIHEDRALS_WITHOUT_HYDROGEN";
    string dihedralsh = "DIHEDRALS_INC_HYDROGEN";
    string dihedralk = "DIHEDRAL_FORCE_CONSTANT";
    string dihedralphase = "DIHEDRAL_PHASE";
    string dihedralperiodicity = "DIHEDRAL_PERIODICITY";
    string scee = "SCEE_SCALE_FACTOR";
    string scnb = "SCNB_SCALE_FACTOR";
    int nphih = parmData[pointers][NPHIH].i;
    int mphia = parmData[pointers][MPHIA].i;
    int nptra = parmData[pointers][NPTRA].i;

    vector<double> sceefac(nptra, 1.2);
    vector<double> scnbfac(nptra, 2.0);

    if (parmData.count(dihedrals) < 1 || parmData[dihedrals].size() != mphia*5)
        throw AmberParmError("Bad (or missing) DIHEDRALS_WITHOUT_HYDROGEN section");
    if (parmData.count(dihedralsh) < 1 || parmData[dihedralsh].size() != nphih*5)
        throw AmberParmError("Bad (or missing) DIHEDRALS_INC_HYDROGEN section");
    if (parmData.count(dihedralk) < 1 || parmData[dihedralk].size() != nptra)
        throw AmberParmError("Bad (or missing) DIHEDRAL_FORCE_CONSTANT section");
    if (parmData.count(dihedralphase) < 1 || 
                parmData[dihedralphase].size() != nptra)
        throw AmberParmError("Bad (or missing) DIHEDRAL_PHASE section");
    if (parmData.count(dihedralperiodicity) < 1 || 
                parmData[dihedralperiodicity].size() != nptra)
        throw AmberParmError("Bad (or missing) DIHEDRAL_PERIODICITY section");
    if (parmData.count(scee) > 0) {
        for (int i = 0; i < nptra; i++)
            sceefac[i] = parmData[scee][i].f;
    }
    if (parmData.count(scnb) > 0) {
        for (int i = 0; i < nptra; i++)
            scnbfac[i] = parmData[scnb][i].f;
    }

    for (int i = 0; i < nphih; i++) {
        int i5 = i * 5;
        int ii = parmData[dihedralsh][i5  ].i / 3;
        int jj = parmData[dihedralsh][i5+1].i / 3;
        int kk = parmData[dihedralsh][i5+2].i / 3;
        int ll = parmData[dihedralsh][i5+3].i / 3;
        int ai = parmData[dihedralsh][i5+4].i - 1;
        double phase = parmData[dihedralphase][ai].f * 180.0 / M_PI;
        int per = (int) parmData[dihedralperiodicity][ai].f;
        bool ignore_end = kk < 0 || ll < 0;
        addDihedral(ii, jj, abs(kk), abs(ll), parmData[dihedralk][ai].f, phase,
                    per, sceefac[ai], scnbfac[ai], ignore_end);
        // Add this to the exception list (NOT the exclusion list)
        if (!ignore_end) {
            if (ii < ll) {
                exception_list[ii].insert(ll);
            } else {
                exception_list[ll].insert(ii);
            }
        }
    }
    for (int i = 0; i < mphia; i++) {
        int i5 = i * 5;
        int ii = parmData[dihedrals][i5  ].i / 3;
        int jj = parmData[dihedrals][i5+1].i / 3;
        int kk = parmData[dihedrals][i5+2].i / 3;
        int ll = parmData[dihedrals][i5+3].i / 3;
        int ai = parmData[dihedrals][i5+4].i - 1;
        double phase = parmData[dihedralphase][ai].f * 180.0 / M_PI;
        int per = (int) parmData[dihedralperiodicity][ai].f;
        bool ignore_end = kk < 0 || ll < 0;
        addDihedral(ii, jj, abs(kk), abs(ll), parmData[dihedralk][ai].f, phase,
                    per, sceefac[ai], scnbfac[ai], ignore_end);
        // Add this to the exception list (NOT the exclusion list)
        if (!ignore_end) {
            if (ii < ll) {
                exception_list[ii].insert(ll);
            } else {
                exception_list[ll].insert(ii);
            }
        }
    }

    // Now go through and build the exclusion list
    string num_exclusions = "NUMBER_EXCLUDED_ATOMS";
    string exclusions = "EXCLUDED_ATOMS_LIST";

    if (parmData.count(num_exclusions) < 1 ||
            parmData[num_exclusions].size() != N)
        throw AmberParmError("Bad (or missing) NUMBER_EXCLUDED_ATOMS section");
    int exclptr = 0;
    for (int i = 0; i < N; i++) {
        int nexcl = parmData[num_exclusions][i].i;
        for (int j = exclptr; j < exclptr + nexcl; j++) {
            int e = parmData[exclusions][j].i - 1;
            if (e < 0) continue;
            if (exception_list[i].count(e) > 0) continue; // it is an exception
            exclusion_list_[i].insert(e);
        }
        exclptr += nexcl;
    }
}

void AmberParm::rdparm(const char* filename) {
    rdparm(string(filename));
}

void AmberParm::printExclusions(int i) {
    cout << "The atoms excluded from atom " << i << " are:" << endl << "\t";
    if (exclusion_list_[i].size() == 0) {
        cout << "None." << endl;
    } else {
        for (set<int>::const_iterator it = exclusion_list_[i].begin();
                it != exclusion_list_[i].end(); it++)
            cout << *it << " ";
        cout << endl;
    }
}

OpenMM::System* AmberParm::createSystem(
                OpenMM::NonbondedForce::NonbondedMethod nonbondedMethod,
                double nonbondedCutoff,
                std::string constraints,
                bool rigidWater,
                std::string implicitSolvent,
                double implicitSolventKappa,
                double implicitSolventSaltConc,
                double temperature,
                double soluteDielectric,
                double solventDielectric,
                bool removeCMMotion,
                double ewaldErrorTolerance,
                bool flexibleConstraints,
                bool useSASA) {

    OpenMM::System* system = new OpenMM::System();

    // Make sure we have a legal choice for constraints and implicitSolvent
    if (constraints != "None" && constraints != "HBonds" &&
            constraints != "AllBonds") {
        string msg = "constraints must be None, HBonds, or AllBonds; not " +
                     constraints;
        throw AmberParmError(msg.c_str());
    }

    if (implicitSolvent != "None" && implicitSolvent != "HCT" &&
            implicitSolvent != "OBC1" && implicitSolvent != "OBC2" &&
            implicitSolvent != "GBn" && implicitSolvent != "GBn2") {
        string msg = "implicitSolvent must be None, HCT, OBC1, OBC2, GBn, or "
                     "GBn2; not " + implicitSolvent;
        throw AmberParmError(msg.c_str());
    }

    if (rigidWater && (constraints != "HBonds" && constraints != "AllBonds")) {
        cerr << "rigidWater is incompatible with constraints=None; setting to "
             << "false" << endl;
    }

    // Catch illegal nonbonded choice for system

    if (isPeriodic()) {
        if (nonbondedMethod == OpenMM::NonbondedForce::CutoffNonPeriodic ||
                nonbondedMethod == OpenMM::NonbondedForce::NoCutoff)
            throw AmberParmError("Illegal nonbondedForce choice for periodic system");
    } else {
        if (nonbondedMethod == OpenMM::NonbondedForce::CutoffPeriodic ||
                nonbondedMethod == OpenMM::NonbondedForce::PME ||
                nonbondedMethod == OpenMM::NonbondedForce::Ewald)
            throw AmberParmError("Illegal nonbondedForce choice for non-periodic system");
    }

    // Add all particles
    for (atom_iterator it = AtomBegin(); it != AtomEnd(); it++) {
        system->addParticle(it->getMass());
    }

    // Add constraints
    bool hcons = constraints == "HBonds" || constraints == "AllBonds";
    bool allcons = constraints == "AllBonds";
    for (bond_iterator it = BondBegin(); it != BondEnd(); it++) {
        Atom a1 = atoms_[it->getAtomI()];
        Atom a2 = atoms_[it->getAtomJ()];
        if (hcons && (a1.getElement() == 1 || a2.getElement() == 1)) {
            system->addConstraint(it->getAtomI(), it->getAtomJ(),
                                  it->getEquilibriumDistance()*NANOMETER_PER_ANGSTROM);
        } else if (allcons) {
            system->addConstraint(it->getAtomI(), it->getAtomJ(),
                                  it->getEquilibriumDistance()*NANOMETER_PER_ANGSTROM);
        }
    }
    
    // Add all bonds
    if (!allcons || flexibleConstraints) {
        OpenMM::HarmonicBondForce *bond_force = new OpenMM::HarmonicBondForce();
        bond_force->setForceGroup(BOND_FORCE_GROUP);
        double conv = ANGSTROM_PER_NANOMETER*ANGSTROM_PER_NANOMETER*JOULE_PER_CALORIE;
        for (bond_iterator it=BondBegin(); it != BondEnd(); it++) {
            // See if this bond needs to be skipped due to constraints
            Atom a1 = atoms_[it->getAtomI()];
            Atom a2 = atoms_[it->getAtomJ()];
            if (hcons && (a1.getElement() == 1 || a2.getElement() == 1) &&
                        !flexibleConstraints) continue;

            bond_force->addBond(it->getAtomI(), it->getAtomJ(),
                                it->getEquilibriumDistance()*NANOMETER_PER_ANGSTROM,
                                2*it->getForceConstant()*conv);
        }
        system->addForce(bond_force);
    }

    // Add all angles
    OpenMM::HarmonicAngleForce *angle_force = new OpenMM::HarmonicAngleForce();
    angle_force->setForceGroup(ANGLE_FORCE_GROUP);
    for (angle_iterator it = AngleBegin(); it != AngleEnd(); it++) {
        angle_force->addAngle(it->getAtomI(), it->getAtomJ(), it->getAtomK(),
                              it->getEquilibriumAngle()*RADIAN_PER_DEGREE,
                              2*it->getForceConstant()*JOULE_PER_CALORIE);
    }
    system->addForce(angle_force);

    // Add all torsions
    OpenMM::PeriodicTorsionForce *dihedral_force = new OpenMM::PeriodicTorsionForce();
    dihedral_force->setForceGroup(DIHEDRAL_FORCE_GROUP);
    for (dihedral_iterator it = DihedralBegin(); it != DihedralEnd(); it++) {
        dihedral_force->addTorsion(it->getAtomI(), it->getAtomJ(),
                                   it->getAtomK(), it->getAtomL(),
                                   it->getPeriodicity(),
                                   it->getPhase()*RADIAN_PER_DEGREE,
                                   it->getForceConstant()*JOULE_PER_CALORIE);
    }
    system->addForce(dihedral_force);

    // Add nonbonded force
    OpenMM::NonbondedForce *nonb_frc = new OpenMM::NonbondedForce();
    nonb_frc->setForceGroup(NONBONDED_FORCE_GROUP);
    nonb_frc->setNonbondedMethod(nonbondedMethod);
    if (nonbondedMethod != OpenMM::NonbondedForce::NoCutoff)
        nonb_frc->setCutoffDistance(nonbondedCutoff*NANOMETER_PER_ANGSTROM);
    const double ONE_SIXTH = 1.0 / 6.0;
    const double SIGMA_SCALE = pow(2, -ONE_SIXTH) * 2 * NANOMETER_PER_ANGSTROM;
    for (atom_iterator it = AtomBegin(); it != AtomEnd(); it++) {
        nonb_frc->addParticle(it->getCharge(),
                              it->getLJRadius()*SIGMA_SCALE,
                              it->getLJEpsilon()*JOULE_PER_CALORIE);
    }
    // Now do exceptions
    const double SIGMA_SCALE2 = pow(2, -ONE_SIXTH) * NANOMETER_PER_ANGSTROM;
    for (dihedral_iterator it = DihedralBegin(); it != DihedralEnd(); it++) {
        if (it->ignoreEndGroups()) continue;
        Atom a1 = atoms_[it->getAtomI()];
        Atom a2 = atoms_[it->getAtomL()];
        double eps = sqrt(a1.getLJEpsilon() * a2.getLJEpsilon()) *
                     JOULE_PER_CALORIE / it->getScnb();
        double sig = (a1.getLJRadius() + a2.getLJRadius()) * SIGMA_SCALE2;
        nonb_frc->addException(it->getAtomI(), it->getAtomL(),
                               a1.getCharge()*a2.getCharge()/it->getScee(),
                               sig, eps);
    }
    // Now do exclusions
    for (int i = 0; i < atoms_.size(); i++) {
        if (exclusion_list_[i].size() == 0) continue; // No exclusions here
        for (set<int>::const_iterator it = exclusion_list_[i].begin();
                it != exclusion_list_[i].end(); it++) {
            nonb_frc->addException(i, *it, 0.0, 1.0, 0.0);
        }
    }
    // Set the ewald error tolerance
    if (nonbondedMethod == OpenMM::NonbondedForce::PME ||
            nonbondedMethod == OpenMM::NonbondedForce::Ewald)
        nonb_frc->setEwaldErrorTolerance(ewaldErrorTolerance);
    system->addForce(nonb_frc);

    // See about removing the center of mass motion
    if (removeCMMotion)
        system->addForce(new OpenMM::CMMotionRemover());

    // If no implicit solvent, we can return system now
    if (implicitSolvent == "None")
        return system;

    // Otherwise, we need to add the GB force
    if (implicitSolventKappa <= 0 && implicitSolventSaltConc > 0)
        implicitSolventKappa = 50.33355 * 0.73 *
                sqrt(implicitSolventSaltConc/(solventDielectric*temperature));

    OpenMM::CustomGBForce *gb_frc = 0;
    if (implicitSolvent == "HCT") {
        gb_frc = GB_HCT(*this, solventDielectric, soluteDielectric, useSASA,
                        nonbondedCutoff, implicitSolventKappa);
    } else if (implicitSolvent == "OBC1") {
        gb_frc = GB_OBC1(*this, solventDielectric, soluteDielectric, useSASA,
                         nonbondedCutoff, implicitSolventKappa);
    } else if (implicitSolvent == "OBC2") {
        gb_frc = GB_OBC2(*this, solventDielectric, soluteDielectric, useSASA,
                         nonbondedCutoff, implicitSolventKappa);
    } else if (implicitSolvent == "GBn") {
        gb_frc = GB_GBn(*this, solventDielectric, soluteDielectric, useSASA,
                        nonbondedCutoff, implicitSolventKappa);
    } else if (implicitSolvent == "GBn2") {
        gb_frc = GB_GBn2(*this, solventDielectric, soluteDielectric, useSASA,
                         nonbondedCutoff, implicitSolventKappa);
    } else {
        stringstream iss;
        iss << "Should not be here; bad GB model " << implicitSolvent;
        throw InternalError(iss.str());
    }

    if (nonbondedMethod == OpenMM::NonbondedForce::NoCutoff)
        gb_frc->setNonbondedMethod(OpenMM::CustomGBForce::NoCutoff);
    else if (nonbondedMethod == OpenMM::NonbondedForce::CutoffNonPeriodic)
        gb_frc->setNonbondedMethod(OpenMM::CustomGBForce::CutoffNonPeriodic);
    else // all remaining options are periodic cutoff...
        gb_frc->setNonbondedMethod(OpenMM::CustomGBForce::CutoffPeriodic);

    system->addForce(gb_frc);

    return system;
}

