/** ambercrd.h
  *
  * This file contains functions for reading and writing Amber coordinate files
  * (both restart/inpcrd and trajectory files
  */

#ifndef AMBERCRD_H
#define AMBERCRD_H

#include <vector>

#include "exceptions.h"

#include "OpenMM.h"

namespace Amber {

class AmberCoordinateFrame {
    public:

        AmberCoordinateFrame(void) :
                natom_(-1), temp0_(0), has_remd_(false), a_(0), b_(0), c_(0),
                alpha_(0), beta_(0), gama_(0), coordinates_(NULL),
                velocities_(NULL) {
            coordinates_ = new std::vector<OpenMM::Vec3>;
            velocities_ = new std::vector<OpenMM::Vec3>;
        }

        ~AmberCoordinateFrame(void);

        /**
         * \brief Sets the input coordinates (in angstroms)
         *
         * If the number of atoms is already known, Amber::AmberCrdError will be
         * thrown if the size of the input array is not equal to the information
         * that is already set.
         *
         * \param crd A pointer to the coordinate array assigned to the
         *            instance. Ownership of this pointer is claimed, and it
         *            will be freed upon destruction of this instance
         */
        void setPositions(std::vector<OpenMM::Vec3> *crd) {
            if (natom_ == -1) {
                natom_ = (int) crd->size();
                coordinates_ = crd;
            } else if (natom_ != crd->size()) {
                throw AmberCrdError("coordinate size mismatch");
            } else {
                coordinates_ = crd;
            }
        }

        /**
         * \brief Returns a reference to the coordinate array
         *
         * \return const Reference to coordinate array
         */
        std::vector<OpenMM::Vec3>& getPositions(void) const {return *coordinates_;}

        /**
         * \brief Sets the input velocities (in angstroms/picosecond).
         *
         * If the number of atoms is already known, Amber::AmberCrdError will be
         * thrown if the size of the input array is not equal to the information
         * that is already set.
         *
         * \param vel A pointer to the velocity array assigned to the instance.
         *            Ownership of this pointer is claimed, and it will be freed
         *            upon destruction of this instance
         */
        void setVelocities(std::vector<OpenMM::Vec3> *vel) {
            if (natom_ == -1) {
                natom_ = (int) vel->size();
                velocities_ = vel;
            } else if (natom_ != vel->size()) {
                throw AmberCrdError("velocity size mismatch");
            } else {
                velocities_ = vel;
            }
        }

        /**
         * \brief Returns a reference to the velocity array
         *
         * \return const Reference to velocity array
         */
        std::vector<OpenMM::Vec3>& getVelocities(void) const {return *velocities_;}

        /**
         * Sets the periodic boundary conditions
         *
         * \param a Length of the first unit cell vector in angstroms
         * \param b Length of the second unit cell vector in angstroms
         * \param c Length of the third unit cell vector in angstroms
         * \param alpha Angle between unit cell vectors b and c in degrees
         * \param beta Angle between unit cell vectors a and c in degrees
         * \param gama Angle between unit cell vectors a and b in degrees
         */
        void setBox(double a, double b, double c,
                    double alpha=90, double beta=90, double gama=90) {
            a_ = a; b_ = b; c_ = c;
            alpha_ = alpha; beta_ = beta; gama_ = gama;
        }

        /// Returns the length of the first unit cell vector in angstroms
        double getBoxA(void) const {return a_;}
        /// Returns the length of the second unit cell vector in angstroms
        double getBoxB(void) const {return b_;}
        /// Returns the length of the third unit cell vector in angstroms
        double getBoxC(void) const {return c_;}
        /// Returns angle between unit cell vectors b and c in degrees
        double getBoxAlpha(void) const {return alpha_;}
        /// Returns angle between unit cell vectors a and c in degrees
        double getBoxBeta(void) const {return beta_;}
        /// Returns angle between unit cell vectors a and b in degrees
        double getBoxGamma(void) const {return gama_;}
        /// Returns the number of atoms defined in this coordinate frame
        int getNatom(void) const {return natom_;}

        /**
         * \brief Reads an Amber coordinate file
         *
         * \param filename Name of the file to read and instantiate data from
         */
        void readRst7(std::string const& filename);
        void readRst7(const char* filename);

        /**
         * \brief Writes an Amber coordinate file
         *
         * \param filename Name of the restart file to write
         * \param netcdf If true, write the restart as NetCDF. If false, use
         *               the ASCII format
         */
        void writeRst7(std::string const& filename, bool netcdf = false);
        void writeRst7(const char* filename, bool netcdf = false);

    private:
        int natom_;
        double temp0_; // Temperature or pH
        bool has_remd_;
        double a_, b_, c_, alpha_, beta_, gama_;
        std::vector<OpenMM::Vec3> *coordinates_;
        std::vector<OpenMM::Vec3> *velocities_;

        void readASCII_(std::string const& filename);
        int readNetCDF_(std::string const& filename);

        void writeASCII_(const char* filename);
        void writeNetCDF_(const char* filename);
};

}; // namespace Amber
#endif /* AMBERCRD_H */
