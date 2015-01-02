/**
 * This file contains the class definition for reading and manipulating NetCDF
 * restart and trajectory files in the Amber format
 */

#ifndef NETCDFFILE_H
#define NETCDFFILE_H

#include <string>
#include <vector>

#include "OpenMM.h"

namespace Amber {

class AmberNetCDFFile {

    public:
        enum FileType {RESTART, TRAJECTORY, AUTOMATIC};

        AmberNetCDFFile(FileType type=AUTOMATIC);
        ~AmberNetCDFFile();

        /**
         * Open a NetCDF file for reading. If the file does not exist, or is not
         * a NetCDF file, a Amber::NotNetcdf exception is thrown
         *
         * \param filename Name of the file to open for reading
         */
        void readFile(std::string const& filename);
        void readFile(const char* filename);

        /** Returns the number of frames in this coordinate file
         */
        int getNumFrames(void) const {return (int)num_frames_;}
        /** Returns the number of atoms in this coordinate file
         */
        int getNatom(void) const {return (int)natom_;}
        /** Returns the number of REMD dimensions in this file
         */
        int getREMDDimension(void) const {return (int)remd_dimension_;}
        /** Returns whether or not a box is present
         */
        bool hasBox(void) const {
            return cell_lengthsVID_ > -1 && cell_anglesVID_ > -1;
        }
        /** Returns whether or not velocities are present
         */
        bool hasVelocities(void) const {return velocitiesVID_ > -1;}
        /** Returns whether or not coordinates are present
         */
        bool hasCoordinates(void) const {return coordinatesVID_ > -1;}
        /** Returns whether or not the file has force information
         */
        bool hasForces(void) const {return forcesVID_ > -1;}
        /** Returns whether or not the file has REMD information
         */
        bool hasREMD(void) const {
            if (temp0VID_ != -1) return true;
            return remdDID_ != -1 && remd_indicesVID_ != -1 &&
                   remd_dimtypeVID_ != -1;
        }
        /**
         * Returns the coordinates in a given frame (must be 0 for restart), if
         * coordinates are not present, AmberCrdError is thrown
         */
        std::vector<OpenMM::Vec3> *getCoordinates(int frame=0) const;
        /**
         * Returns the velocities in a given frame (must be 0 for restart), if
         * velocities are not present, AmberCrdError is thrown
         */
        std::vector<OpenMM::Vec3> *getVelocities(int frame=0) const;
        /**
         * Returns the forces in a given frame. If forces are not present, an
         * AmberCrdError is thrown.
         */
        std::vector<OpenMM::Vec3> *getForces(int frame=0) const;
        /**
         * Returns the 3 cell lengths in Angstroms in a given frame. If cell
         * lengths are not present, an AmberCrdError is thrown.
         */
        OpenMM::Vec3 getCellLengths(int frame=0) const;
        /**
         * Returns the 3 cell angles in degrees in a given frame. If cell angles
         * are not present, an AmberCrdError is thrown.
         */
        OpenMM::Vec3 getCellAngles(int frame=0) const;
        /**
         * Returns the program that created the NetCDF file
         */
        std::string getProgram(void) const {return program_;}
        /**
         * Returns the program version that created the NetCDF file
         */
        std::string getProgramVersion(void) const {return programVersion_;}
        /**
         * Returns the appplication that created the NetCDF file
         */
        std::string getApplication(void) const {return application_;}
        /**
         * Returns the title of the NetCDF file
         */
        std::string getTitle(void) const {return title_;}
    private:
        // Utility functions
        /**
         * Gets a textual attribute and returns it as a string. If not present
         * but required, an AmberCrdError is thrown. If not present but not
         * required, an empty string is returned
         *
         * \param varID The ID of the variable to get the attribute from
         *              (NC_GLOBAL for a global attribute)
         * \param attr The name of the attribute whose text should be gotten
         * \param required (optional) Default is false -- if required and not
         *                            present, an exception is thrown
         * \return string containing the attribute text
         */
        std::string GetAttributeText_(int varID, const char* attr,
                                      bool required=false) const;
        /**
         * Gets a floating point number attribute and returns it as a string. If
         * not present, the default value is returned
         *
         * \param varID The ID of the variable to get the attribute from
         *              (NC_GLOBAL for a global attribute)
         * \param attr The name of the attribute whose value should be gotten
         * \param default_ The default value to return if not found
         *
         * \return value of the attribute (or default_ if it does not exist)
         */
        double GetAttributeFloat_(int varID, const char* attr,
                                  double default_=1.0) const;
        /**
         * Gets the dimension ID for a particular dimension. If the dimension
         * does not exist, a -1 is returned and value is assigned to 0
         *
         * \param name Name of the dimension to retrieve
         * \param value Size of the dimension returned to caller (by reference)
         *
         * \return dimension ID number (or -1 if dimension does not exist)
         */
        int GetDimensionID_(const char* name, size_t &value);
        /**
         * Gets the variable ID for a particular variable. If the variable does
         * not exist, a -1 is returned.
         *
         * \param name Name of the variable to retrieve
         *
         * \return variable ID number
         */
        int GetVariableID_(const char* name);

        // File descriptor and dimensions
        int ncid_, atomDID_, frameDID_, spatialDID_, cell_spatialDID_,
            cell_angularDID_, labelDID_, remdDID_;
        // Variables
        int spatialVID_, cell_spatialVID_, cell_angularVID_, timeVID_,
            coordinatesVID_, cell_lengthsVID_, cell_anglesVID_, velocitiesVID_,
            forcesVID_, temp0VID_, remd_dimtypeVID_, remd_indicesVID_;
        bool is_open_, is_old_;

        FileType type_;
        // Internal variables
        size_t num_frames_, natom_, remd_dimension_, label_;

        std::string program_, programVersion_, application_, title_;
};

};
#endif /* NETCDFFILE_H */
