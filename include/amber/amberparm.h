/** amberparm.h
  * 
  * This file contains the functionality for reading Amber topology files and
  * instantiating molecular structure classes based on them.
  *
  */

#ifndef AMBERPARM_H
#define AMBERPARM_H

#include <cmath>
#include <string>
#include <set>

#include "topology.h"
#include "readparm.h"

#include "OpenMM.h"

namespace Amber {

class AmberParm {
    public:
        /**
         * Optionally parses an Amber topology file
         *
         * \param filename Name of the Amber prmtop file to parse, if provided
         */
        AmberParm(void) : ifbox_(0) {}
        AmberParm(std::string const& filename);
        AmberParm(const char* filename);

        // Iterators
        typedef AtomList::const_iterator atom_iterator;
        typedef BondList::const_iterator bond_iterator;
        typedef AngleList::const_iterator angle_iterator;
        typedef DihedralList::const_iterator dihedral_iterator;

        /// Iterator through atoms of this system
        atom_iterator AtomBegin(void) const {return atoms_.begin();}
        atom_iterator AtomEnd(void) const {return atoms_.end();}
        /// Iterator through bonds of this system
        bond_iterator BondBegin(void) const {return bonds_.begin();}
        bond_iterator BondEnd(void) const {return bonds_.end();}
        /// Iterator through angles of this system
        angle_iterator AngleBegin(void) const {return angles_.begin();}
        angle_iterator AngleEnd(void) const {return angles_.end();}
        /// Iterator through dihedrals of this system
        dihedral_iterator DihedralBegin(void) const {return dihedrals_.begin();}
        dihedral_iterator DihedralEnd(void) const {return dihedrals_.end();}

        /* To add parameters. They can either be added using the parameter type
         * directly or by providing the atom indexes and the parameters. If the
         * indexes of the provided atoms are not sequential, an AmberParmError
         * will be thrown. If atom properties are provided, the indexes will be
         * assigned sequentially.
         *
         * If any bonds, angles, or dihedrals are assigned with atom indexes
         * out of range of the main atom list, an AmberParmError will be thrown
         * (meaning that atoms MUST be added first).
         */
        /**
         * \brief Adds an atom to the system
         *
         * \param new_atom An Amber::Atom to add to this system
         *
         * Atoms must be added sequentially, meaning that if the index of
         * new_atom is not one more than the index of the last atom added, a
         * Amber::AmberParmError will be thrown
         */
        void addAtom(Atom& new_atom);
        /**
         * \brief Adds an atom to the system with the given attributes
         *
         * \param name The name of the atom to add
         * \param type The name of the atom type to add
         * \param element The atomic number of the new atom
         * \param mass The molecular mass of the atom in atomic mass units
         * \param charge The partial atomic charge of the atom in units of
         *               electrons
         * \param lj_rad The Lennard-Jones radius (Rmin/2) of this atom in
         *               angstroms
         * \param lj_eps The Lennard-Jones epsilon parameter of this atom in
         *               kcal/mol
         * \param gb_rad The GB intrinsic solvent radius of this atom
         * \param gb_screen The GB screening factor for this atom
         */
        void addAtom(std::string const& name, std::string const& type,
                     int element, double mass, double charge, double lj_rad,
                     double lj_eps, double gb_rad, double gb_screen);
        void addAtom(const char *name, const char *type, int element,
                     double mass, double charge, double lj_rad, double lj_eps,
                     double gb_rad, double gb_screen);

        /**
         * \brief Adds a bond to the system
         *
         * \param new_bond A reference to a Amber::Bond instance to be added to this
         *                 system
         *
         * If the bond references any atom indices that are out of range of the
         * number of atoms already added to this system, an AmberParmError is
         * thrown
         */
        void addBond(Bond& new_bond);
        /**
         * \brief Adds a new bond to the system based on atom indices and
         *        parameters
         *
         * \param i Index of the first atom in the bond
         * \param j Index of the second atom in the system
         * \param kf Force constant of this bond, in kcal/mole/angstrom^2
         * \param req Equilibrium distance of this bond in angstroms
         *
         * If any of the indexes are out of range for the atoms already added,
         * an Amber::AmberParmError is thrown
         */
        void addBond(int i, int j, double kf, double req);
        /**
         * \brief Adds an angle to the system
         *
         * \param new_angle Reference to an Amber::Angle instance to be added to
         *                  the system
         *
         * If the atom indexes in the angle are out of range of the number of
         * atoms currently added to the system, an Amber::AmberParmError is
         * thrown
         */
        void addAngle(Angle& new_angle);
        /**
         * \brief Adds an angle to the system
         *
         * \param i Index of the first atom in the angle to add
         * \param j Index of the second atom in the angle to add
         * \param k Index of the third atom in the angle to add
         * \param kf Force constant of the angle in kcal/mol/radians^2
         * \param theteq Equilibrum angle in degrees
         */
        void addAngle(int i, int j, int k, double kf, double theteq);
        /**
         * \brief Adds a dihedral to the system
         *
         * \param new_dihedral Reference to an Amber::Dihedral to add to the
         *                     system
         *
         * If any of the indexes in new_dihedral are out of range for the total
         * number of atoms added to the system so far, an Amber::AmberParmError
         * is thrown
         */
        void addDihedral(Dihedral& new_dihedral);
        /**
         * \brief Adds a dihedral to the system
         *
         * \param i Index of the first atom in the dihedral to add
         * \param j Index of the second atom in the dihedral to add
         * \param k Index of the third atom in the dihedral to add
         * \param l Index of the fourth atom in the dihedral to add
         * \param kf Force constant of the angle in kcal/mol
         * \param phase Phase shift of the dihedral in degrees
         * \param periodicity Periodicity of the dihedral
         * \param scee 1-4 electrostatic scaling factor for this pair
         * \param scnb 1-4 van der Waals scaling factor for this pair
         * \param ignore_end If true, the -14 nonbonded parameters are not
         *                   computed for this pair (e.g., in ring systems,
         *                   multiterm dihedrals, and impropers)
         *
         * The scee/scnb terms are effectively ignored when ignore_end is true.
         * If any atom indices are out of range of the currently added atoms, an
         * Amber::AmberParmError is thrown.
         */
        void addDihedral(int i, int j, int k, int l, double kf, double phase,
                         int periodicity, double scee, double scnb,
                         bool ignore_end);

        /// Returns a reference to the list of atoms in the system
        AtomList Atoms(void) const {return atoms_;}
        /// Returns a reference to the list of bonds in the system
        BondList Bonds(void) const {return bonds_;}
        /// Returns a reference to the list of angles in the system
        AngleList Angles(void) const {return angles_;}
        /// Returns a reference to the list of dihedrals in the system
        DihedralList Dihedrals(void) const {return dihedrals_;}
        /// Returns a reference to the beginning of each residue in the system
        std::vector<int> ResiduePointers(void) const {return residue_pointers_;}
        /// Returns the list of residue names of each residue in the sysetm
        std::vector<std::string> ResidueLabels(void) const {
            return residue_labels_;
        }
        /**
         * \brief Determines if 2 atoms are excluded
         *
         * \param i Index of the first atom
         * \param j Index of the second atom
         *
         * \return true if atoms i and j are excluded, false if they are not (or
         *         if they are 1-4 exceptions)
         */
        bool isExcluded(int i, int j) const {
            if (i == j) return true;
            if (i < j) {
                return exclusion_list_[i].count(j) > 0;
            }
            return exclusion_list_[j].count(i) > 0;
        }

        /**
         * \brief Indicates whether this system is periodic or not
         *
         * \return true if the system has periodic boundaries, false otherwise
         */
        bool isPeriodic(void) const {return ifbox_ > 0;}
        /**
         * \brief Tells what kind, if any, periodic boundary conditions exist
         *
         * \return 0 for no box, 1 for orthorhombic box, 2 for generalized
         *         triclinic cells
         */
        int IfBox(void) const {return ifbox_;}
        /**
         * \brief Prints all atoms that are excluded from a single atom to
         *        stdout
         *
         * \param i Index of the atom for which to print exclusions
         */
        void printExclusions(int i);

        /**
         * Read a prmtop file and instantiate a structure from it
         *
         * \param filename Name of the prmtop file to read
         */
        void rdparm(std::string const& filename);
        void rdparm(const char* filename);

        /**
         * \brief Creates and returns a pointer to an OpenMM::System
         *
         * \param nonbondedMethod The nonbonded method to use for the system.
         *                        Can be OpenMM::NonbondedForce::NoCutoff,
         *                        CutoffNonPeriodic, CutoffPeriodic, or PME
         * \param nonbondedCutoff Cutoff of nonbonded interactions in angstroms
         * \param constraints Can be "None", "HBonds" (for typical SHAKE), or
         *                    "AllBonds" to constrain all bonds
         * \param rigidWater Not used right now, must be false
         * \param implicitSolvent Name of the GB model to use, if any. Can be
         *                        "None", "HCT", "OBC1", "OBC2", "GBn" or "GBn2"
         * \param implicitSolventKappa kappa in the Debye salt concentration
         *      equation (in angstrom^-1). This takes precedence over
         *      implicitSolventSaltConc (below)
         * \param implicitSolventSaltConc Alternative way to specify salt
         *      concentration (converted internally to kappa). Given in Molar
         * \param temperature Temperature used for converting salt concentration
         *                    to kappa
         * \param soluteDielectric Dielectric constant to use for solute in GB
         * \param solventDielectric Dielectric constant to use for solvent in GB
         * \param removeCMMotion If true, remove COM motion. If false, don't
         * \param ewaldErrorTolerance Ewald error tolerance for PME and Ewald
         * \param flexibleConstraints If true, compute energy of constrained
         *                            bonds. If false, don't.
         * \param useSASA If true, use the ACE SASA-based non-polar solvation
         *                free energy model for the SA part of GBSA calculations
         */
        OpenMM::System* createSystem(
            OpenMM::NonbondedForce::NonbondedMethod nonbondedMethod=OpenMM::NonbondedForce::NoCutoff,
            double nonbondedCutoff=10.0,
            std::string constraints=std::string("None"),
            bool rigidWater=false,
            std::string implicitSolvent=std::string("None"),
            double implicitSolventKappa=0.0,
            double implicitSolventSaltConc=0.0,
            double temperature=298.15,
            double soluteDielectric=1.0,
            double solventDielectric=78.5,
            bool removeCMMotion=true,
            double ewaldErrorTolerance=0.0005,
            bool flexibleConstraints=true,
            bool useSASA=false);

    private:
        int ifbox_;
        AtomList atoms_;
        BondList bonds_;
        AngleList angles_;
        DihedralList dihedrals_;
        std::vector<int> residue_pointers_;
        std::vector<std::string> residue_labels_;
        std::vector<std::set<int> > exclusion_list_;
};

enum ForceGroup {BOND_FORCE_GROUP=0, ANGLE_FORCE_GROUP, DIHEDRAL_FORCE_GROUP,
                 NONBONDED_FORCE_GROUP};

}; // namespace Amber

#endif /* AMBERPARM_H */
