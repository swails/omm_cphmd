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
                alpha_(0), beta_(0), gama_(0)
            { }

        /// Sets the input coordinates (in angstroms)
        void setPositions(std::vector<OpenMM::Vec3> const &crd) {
            if (natom_ == -1) {
                natom_ = (int) crd.size();
                coordinates_ = crd;
            } else if (natom_ != crd.size()) {
                throw AmberCrdError("coordinate size mismatch");
            } else {
                coordinates_ = crd;
            }
        }

        std::vector<OpenMM::Vec3> getPositions(void) const {return coordinates_;}

        /// Sets the input velocities (in angstroms/picosecond)
        void setVelocities(std::vector<OpenMM::Vec3> const &vel) {
            if (natom_ == -1) {
                natom_ = (int) vel.size();
                velocities_ = vel;
            } else if (natom_ != vel.size()) {
                throw AmberCrdError("velocity size mismatch");
            } else {
                velocities_ = vel;
            }
        }

        std::vector<OpenMM::Vec3> getVelocities(void) const {return velocities_;}

        /// Gets and sets periodic boundary parameters
        void setBox(double a, double b, double c,
                    double alpha=90, double beta=90, double gama=90) {
            a_ = a; b_ = b; c_ = c;
            alpha_ = alpha; beta_ = beta; gama_ = gama;
        }

        double getBoxA(void) const {return a_;}
        double getBoxB(void) const {return b_;}
        double getBoxC(void) const {return c_;}
        double getBoxAlpha(void) const {return alpha_;}
        double getBoxBeta(void) const {return beta_;}
        double getBoxGamma(void) const {return gama_;}

        int getNatom(void) const {return natom_;}

        /// Reads an Amber coordinate file
        void readRst7(std::string const& filename);
        void readRst7(const char* filename);

        /// Writes an Amber coordinate file
        void writeRst7(std::string const& filename, bool netcdf = false);
        void writeRst7(const char* filename, bool netcdf = false);

    private:
        int natom_;
        double temp0_; // Temperature or pH
        bool has_remd_;
        double a_, b_, c_, alpha_, beta_, gama_;
        std::vector<OpenMM::Vec3> coordinates_;
        std::vector<OpenMM::Vec3> velocities_;
};

}; // namespace Amber
#endif /* AMBERCRD_H */
