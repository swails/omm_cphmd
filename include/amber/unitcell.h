/** unitcell.h
 *
 * Contains functionality for manipulating and dealing with unit cell vectors
 */
#ifndef UNITCELL_H
#define UNITCELL_H

#include "OpenMM.h"

namespace Amber {

class UnitCell {
    public:
        /**
         * \brief Constructs a unit cell where all unit cell vectors are null
         *        vectors (i.e., (0, 0, 0))
         */
        UnitCell(void) :
            a_(OpenMM::Vec3()),
            b_(OpenMM::Vec3()),
            c_(OpenMM::Vec3()) {}
        /**
         * \brief Constructs a periodic box from the lengths of the unit cell
         *        and angles between the vectors
         *
         * \param a Length of the first unit cell vector in Angstroms
         * \param b Length of the second unit cell vector in Angstroms
         * \param c Length of the third unit cell vector in Angstroms
         * \param alpha Angle between vectors b and c in degrees
         * \param beta Angle between vectors a and c in degrees
         * \param gama Angle between vectors b and c in degrees
         */
        UnitCell(double a, double b, double c,
                 double alpha, double beta, double gama);
        /**
         * \brief Constructs a periodic box from the three unit cell vectors
         *
         * \param a First unit cell vector in Angstroms
         * \param b Second unit cell vector in Angstroms
         * \param c Third unit cell vector in Angstroms
         */
        UnitCell(const OpenMM::Vec3 &a, const OpenMM::Vec3 &b,
                 const OpenMM::Vec3 &c);
        /**
         * \brief Gets the length of the first unit cell vector
         *
         * \return Length of unit cell vector a in Angstroms
         */
        double getLengthA(void) const;
        /**
         * \brief Gets the length of the second unit cell vector
         *
         * \return Length of unit cell vector b in Angstroms
         */
        double getLengthB(void) const;
        /**
         * \brief Gets the length of the third unit cell vector
         *
         * \return Length of unit cell vector c in Angstroms
         */
        double getLengthC(void) const;
        /**
         * \brief Gets the angle between unit cell vectors B and C
         *
         * \return Angle between unit cell vectors B and C in degrees
         */
        double getAlpha(void) const;
        /**
         * \brief Gets the angle between unit cell vectors A and C
         *
         * \return Angle between unit cell vectors A and C in degrees
         */
        double getBeta(void) const;
        /**
         * \brief Gets the angle between unit cell vectors A and B
         *
         * \return Angle between unit cell vectors A and B in degrees
         */
        double getGamma(void) const;
        /**
         * \brief Gets the first unit cell vector
         *
         * \return Gets the first unit cell vector in Angstroms
         */
        OpenMM::Vec3 getVectorA(void) const {return a_;}
        /**
         * \brief Gets the second unit cell vector
         *
         * \return Gets the second unit cell vector in Angstroms
         */
        OpenMM::Vec3 getVectorB(void) const {return b_;}
        /**
         * \brief Gets the third unit cell vector
         *
         * \return Gets the third unit cell vector in Angstroms
         */
        OpenMM::Vec3 getVectorC(void) const {return c_;}
        /**
         * \brief Sets the first unit cell vector
         *
         * \param a The first unit cell vector in Angstroms
         */
        void setVectorA(OpenMM::Vec3 const& a) {a_ = a;}
        /**
         * \brief Sets the second unit cell vector
         *
         * \param b The second unit cell vector in Angstroms
         */
        void setVectorB(OpenMM::Vec3 const& b) {b_ = b;}
        /**
         * \brief Sets the third unit cell vector
         *
         * \param c The third unit cell vector in Angstroms
         */
        void setVectorC(OpenMM::Vec3 const& c) {c_ = c;}
        /**
         * \brief Sets the unit cell from the cell lengths and cell angles
         *
         * \param a Length of the first unit cell vector in Angstroms
         * \param b Length of the second unit cell vector in Angstroms
         * \param c Length of the third unit cell vector in Angstroms
         * \param alpha Angle between vectors b and c in degrees
         * \param beta Angle between vectors a and c in degrees
         * \param gama Angle between vectors b and c in degrees
         */
        void setUnitCell(double a, double b, double c,
                         double alpha, double beta, double gama);
    private:
        OpenMM::Vec3 a_, b_, c_;
};

};

#endif /* UNITCELL_H */
