/** topology.h
 *
 * This header contains classes used for storing parameter types like bonds,
 * angles, and dihedrals, read from an Amber topology file
 */
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <string>
#include <vector>

namespace CpHMD {

class Bond {
    public:
        /** A bond between two atoms
          * \param i : Index (from 0) of the first atom involved in the bond
          * \param j : Index (from 0) of the second atom involved in the bond
          * \param kf : Force constant (kcal/mol/A^2) of this bond
          * \param req : Equilibrium distance (A) of this bond
          */
        Bond(int i, int j, double kf, double req) :
            i_(i), j_(j), k_(kf), req_(req) {}

        int getAtomI(void) const {return i_;}
        int getAtomJ(void) const {return j_;}

        double getForceConstant(void) const {return k_;}
        double getEquilibriumDistance(void) const {return req_;}

    private:
        int i_, j_;      //< Atom indices
        double k_, req_; //< Bond parameters
};

class Angle {
    public:
        /** An angle between 3 atoms
          * \param i : Index of the first atom involved in the angle
          * \param j : Index of the second atom involved in the angle
          * \param k : Index of the third atom involved in the angle
          * \param kf : Force constant of the angle (kcal/mol/radians^2)
          * \param theteq : Equilibrium angle (in degrees)
          */
        Angle(int i, int j, int k, double kf, double theteq) :
            i_(i), j_(j), k_(k), kf_(kf), theteq_(theteq) {}

        int getAtomI(void) const {return i_;}
        int getAtomJ(void) const {return j_;}
        int getAtomK(void) const {return k_;}

        double getForceConstant(void) const {return kf_;}
        double getEquilibriumAngle(void) const {return theteq_;}

    private:
        int i_, j_, k_; ///< Atom indices
        double kf_, theteq_; ///< Angle parameters
};

class Dihedral {
    public:
        /** A torsion between 4 atoms around a single bond
          * \param i : Index of the first atom involved in the torsion
          * \param j : Index of the second atom involved in the torsion
          * \param k : Index of the third atom involved in the torsion
          * \param l : Index of the fourth atom involved in the torsion
          * \param kf : Force constant (kcal/mol)
          * \param phase : The phase shift of the torsion angle (degrees)
          * \param periodicity : The periodicity of the torsion parameter
          * \param ignore_end : Whether or not the end-group interactions for
          *                     this torsion are ignored or not
          */
        Dihedral(int i, int j, int k, int l, double kf, double phase,
                 int periodicity, bool ignore_end) :
            i_(i), j_(j), k_(k), l_(l), kf_(kf), phase_(phase),
            periodicity_(periodicity), ignore_end_(ignore_end) {}

        int getAtomI(void) const {return i_;}
        int getAtomJ(void) const {return j_;}
        int getAtomK(void) const {return k_;}
        int getAtomL(void) const {return l_;}

        double getForceConstant(void) const {return kf_;}
        double getPhase(void) const {return phase_;}
        int getPeriodicity(void) const {return periodicity_;}
        bool ignoreEndGroups(void) const {return ignore_end_;}

    private:
        int i_, j_, k_, l_; ///< Atom indices
        double kf_, phase_; ///< Dihedral parameters
        int periodicity_;   ///< Phase shift angle
        bool ignore_end_;   ///< Whether end-group 1-4 ixns are ignored or not
};

class Atom {
    public:
        /** An atom that contains the atomic properties of a biomolecular system
          * \param i : Index of this particular atom in the system
          * \param name : The name of the atom
          * \param type : The name of the atom type
          * \param mass : The atomic mass of this atom (daltons)
          * \param charge : The partial atomic charge (electron fractions)
          * \param lj_rad : Rmin/2 Lennard-Jones parameter (Angstroms)
          * \param lj_eps : epsilon Lennard-Jones parameter (kcal/mol)
          * \param gb_rad : GB intrinsic radius (Angstroms)
          * \param gb_screen : The GB screening factor
          */
        Atom(int i, std::string const& name, std::string const& type,
             double mass, double charge, double lj_rad, double lj_eps,
             double gb_rad, double gb_screen) :
            index_(i), name_(name), type_(type), mass_(mass), charge_(charge),
            lj_rad_(lj_rad), lj_eps_(lj_eps), gb_rad_(gb_rad),
            gb_screen_(gb_screen) {}

        Atom(int i, const char* name, const char* type, double mass,
             double charge, double lj_rad, double lj_eps, double gb_rad,
             double gb_screen) : 
            index_(i), name_(std::string(name)), type_(std::string(type)),
            mass_(mass), charge_(charge), lj_rad_(lj_rad), lj_eps_(lj_eps),
            gb_rad_(gb_rad), gb_screen_(gb_screen) {}

        std::string getName(void) const {return name_;}
        std::string getType(void) const {return type_;}

        int getIndex(void) const {return index_;}

        double getMass(void) const {return mass_;}
        double getCharge(void) const {return charge_;}
        double getLJRadius(void) const {return lj_rad_;}
        double getLJEpsilon(void) const {return lj_eps_;}
        double getGBRadius(void) const {return gb_rad_;}
        double getGBScreen(void) const {return gb_screen_;}

    private:
        int index_;               ///< Index of atom in the system
        std::string name_, type_; ///< Name and type of the atom
        double mass_, charge_, lj_rad_,      ///< Various atomic properties
               lj_eps_, gb_rad_, gb_screen_; ///< and parameters
};

// Containers for various parameter types

typedef std::vector<Bond> BondList;
typedef std::vector<Angle> AngleList;
typedef std::vector<Dihedral> DihedralList;
typedef std::vector<Atom> AtomList;

}; // namespace CpHMD
#endif /* TOPOLOGY_H */
