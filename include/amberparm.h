/** amberparm.h
  * 
  * This file contains the functionality for reading Amber topology files and
  * instantiating molecular structure classes based on them.
  *
  */

#ifndef AMBERPARM_H
#define AMBERPARM_H

#include <string>

#include "topology.h"
#include "readparm.h"

namespace CpHMD {

class AmberParm {
    public:
        /** Reads an Amber topology file and sets up data structures defining
          * the topology and parameters. If a filename is passed to the
          * constructor, it is read.
          */
        AmberParm(void) {}
        AmberParm(std::string const& filename);
        AmberParm(const char* filename);

        // Iterators
        typedef AtomList::const_iterator atom_iterator;
        typedef BondList::const_iterator bond_iterator;
        typedef AngleList::const_iterator angle_iterator;
        typedef DihedralList::const_iterator dihedral_iterator;

        atom_iterator AtomBegin(void) const {return atoms_.begin();}
        atom_iterator AtomEnd(void) const {return atoms_.end();}
        bond_iterator BondBegin(void) const {return bonds_.begin();}
        bond_iterator BondEnd(void) const {return bonds_.end();}
        angle_iterator AngleBegin(void) const {return angles_.begin();}
        angle_iterator AngleEnd(void) const {return angles_.end();}
        dihedral_iterator DihedralBegin(void) const {return dihedrals_.begin();}
        dihedral_iterator DihedralEnd(void) const {return dihedrals_.end();}

        /* To add parameters. They can either be added using the parameter type
         * directly or by providing the atom indexes and the parameters
         */
        void addAtom(Atom& new_atom);
        void addAtom(int i, std::string const& name, std::string const& type,
                     int element, double mass, double charge, double lj_rad,
                     double lj_eps, double gb_rad, double gb_screen);
        void addAtom(int i, const char *name, const char *type, int element,
                     double mass, double charge, double lj_rad, double lj_eps,
                     double gb_rad, double gb_screen);

        void addBond(Bond& new_bond);
        void addBond(int i, int j, double kf, double req);

        void addAngle(Angle& new_angle);
        void addAngle(int i, int j, int k, double kf, double theteq);

        void addDihedral(Dihedral& new_dihedral);
        void addDihedral(int i, int j, int k, int l, double kf, double phase,
                         int periodicity, bool ignore_end);

        // Getters for the parameter types -- cannot modify them
        AtomList Atoms(void) const {return atoms_;}
        BondList Bonds(void) const {return bonds_;}
        AngleList Angles(void) const {return angles_;}
        DihedralList Dihedrals(void) const {return dihedrals_;}

        /// Read a prmtop file and instantiate a structure from it.
        void rdparm(std::string const& filename);
        void rdparm(const char* filename);

    private:
        AtomList atoms_;
        BondList bonds_;
        AngleList angles_;
        DihedralList dihedrals_;
};

};

#endif /* AMBERPARM_H */
