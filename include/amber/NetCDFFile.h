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

        /**
         * Sets the file type of the object. It can either be
         * AmberNetCDFFile::RESTART, AmberNetCDFFile::TRAJECTORY, or
         * AmberNetCDFFile::AUTOMATIC.
         *
         * \param type If AmberNetCDFFile::AUTOMATIC, the file type will be
         *             determined when AmberNetCDFFile::readFile is called
         *             (attempting to write will result in a thrown exception)
         */
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
         *
         * \param frame The frame of the NetCDF file to extract coordinates for
         *
         * \return A newly allocated pointer to a std::vector of OpenMM::Vec3.
         *         The caller is responsible for deallocating the result
         */
        std::vector<OpenMM::Vec3> *getCoordinates(int frame=0) const;
        /**
         * Returns the velocities in a given frame (must be 0 for restart), if
         * velocities are not present, AmberCrdError is thrown.
         *
         * \param frame The farme of the NetCDF file to extract velocities for
         *
         * \return A newly allocated pointer to a std::vector of OpenMM::Vec3.
         *         The caller is responsible for deallocating the result
         */
        std::vector<OpenMM::Vec3> *getVelocities(int frame=0) const;
        /**
         * Returns the forces in a given frame. If forces are not present, an
         * AmberCrdError is thrown.
         *
         * \param frame The farme of the NetCDF file to extract forces for
         *
         * \return A newly allocated pointer to a std::vector of OpenMM::Vec3.
         *         The caller is responsible for deallocating the result
         */
        std::vector<OpenMM::Vec3> *getForces(int frame=0) const;
        /**
         * Returns the 3 cell lengths in Angstroms in a given frame. If cell
         * lengths are not present, an AmberCrdError is thrown.
         *
         * \param frame The frame of the NetCDF file to extract the cell lengths
         *              for
         *
         * \return An OpenMM::Vec3 containing the 3 cell lengths in Angstroms
         */
        OpenMM::Vec3 getCellLengths(int frame=0) const;
        /**
         * Returns the 3 cell angles in degrees in a given frame. If cell angles
         * are not present, an AmberCrdError is thrown.
         *
         * \param frame The frame of the NetCDF file to extract the cell angles
         *              for
         *
         * \return An OpenMM::Vec3 containing the 3 cell angles in degrees
         */
        OpenMM::Vec3 getCellAngles(int frame=0) const;
        /**
         * Returns the program that created the NetCDF file
         *
         * \return string containing the text in the "program" attribute
         */
        std::string getProgram(void) const {return program_;}
        /**
         * Returns the program version that created the NetCDF file
         *
         * \return the string containing the text in the "programVersion"
         *         attribute
         */
        std::string getProgramVersion(void) const {return programVersion_;}
        /**
         * Returns the appplication that created the NetCDF file
         *
         * \return the string containing the text in the "application" attribute
         */
        std::string getApplication(void) const {return application_;}
        /**
         * Returns the title of the NetCDF file
         *
         * \return the string containing the text in the "title" attribute
         */
        std::string getTitle(void) const {return title_;}

        // Functions for writing new files
        /**
         * \brief Opens a new file to be written with data
         *
         * \param filename Name of the file to write
         * \param natom Number of atoms present in this system
         * \param hasCrd If true, coordinates will be provided. Not otherwise
         * \param hasVel If true, velocities will be provided. Not otherwise
         * \param hasFrc If true, forces will be provided. Not otherwise
         * \param hasBox If true, unit cells will be provided. Not otherwise
         * \param hasRemd If true, REMD information needs to be written
         * \param remdDimension If hasRemd, this is the number of REMD
         *      dimensions that will be present. Dimension of 0 means multi-D
         *      REMD info will not be written (Do not set this for 1-D T-REMD or
         *      1-D pH-REMD, as that will use the standard temp0 machinery)
         * \param title The title of the file written to the global attributes
         * \param application The name of the program calling these functions
         *      (will be written as an attribute to the NetCDF file)
         */
        void writeFile(std::string const &filename, int natom, bool hasCrd,
                bool hasVel, bool hasFrc, bool hasRemd, bool hasBox,
                int remdDimension, std::string const& title,
                std::string const& application);
        /**
         * See the other definition for writeFile (but with const char* instead
         * of std::string for textual input)
         */
        void writeFile(const char* filename, int natom, bool hasCrd,
                bool hasVel, bool hasFrc, bool hasBox, bool hasRemd,
                int remdDimension, const char* title, const char* application);
        /**
         * \brief Writes a set of coordinates to the current file
         *
         * \param coordinates The coordinates to write to the NetCDF file in
         *                    angstroms
         */
        void setCoordinates(std::vector<OpenMM::Vec3> const &coordinates);
        /**
         * \brief Writes a set of coordinates to the current file
         *
         * \param coordinates The coordinates to write to the NetCDF file in
         *                    nanometers
         *
         * This function is provided for convenience to work with programs that
         * use the kJ-nm-s unit system (e.g., OpenMM)
         */
        void setCoordinatesNm(std::vector<OpenMM::Vec3> const &coordinates);
        /**
         * \brief Writes a set of velocities to the current file
         *
         * \param velocities The velocities to write to the NetCDF file in
         *                   Angstroms/picosecond/20.455
         */
        void setVelocities(std::vector<OpenMM::Vec3> const& velocities);
        /**
         * \brief Writes a set of velocities to the current file
         *
         * \param velocities The velocities to write to the NetCDF file in
         *                   nanometers/picosecond
         *
         * This function is provided for convenience to work with programs that
         * use the kJ-nm-s unit system (e.g., OpenMM)
         */
        void setVelocitiesNmPerPs(std::vector<OpenMM::Vec3> const& velocities);
        /**
         * \brief Writes a set of forces to the current file
         *
         * \param forces The forces to write to the NetCDF file in kcal/Ang
         */
        void setForces(std::vector<OpenMM::Vec3> const& forces);
        /**
         * \brief Writes a set of forces to the current file
         *
         * \param forces The forces to write to the NetCDF file in kJ/nm
         *
         * This function is provided for convenience to work with programs that
         * use the kJ-nm-s unit system (e.g., OpenMM)
         */
        void setForcesKJPerNm(std::vector<OpenMM::Vec3> const& forces);
        /**
         * \brief Writes cell lengths and angles from a set of unit cell vectors
         *
         * \param a The first unit cell vector (in Angstroms)
         * \param b The second unit cell vector (in Angstroms)
         * \param c The third unit cell vector (in Angstroms)
         */
        void setUnitCell(OpenMM::Vec3 const &a, OpenMM::Vec3 const &b,
                         OpenMM::Vec3 const &c);
        /**
         * \brief Writes cell lengths and angles from a set of unit cell vectors
         *
         * \param a The first unit cell vector (in nm)
         * \param b The second unit cell vector (in nm)
         * \param c The third unit cell vector (in nm)
         *
         * This function is provided for convenience to work with programs that
         * use the kJ-nm-s unit system (e.g., OpenMM)
         */
        void setUnitCellNm(OpenMM::Vec3 const &a, OpenMM::Vec3 const &b,
                           OpenMM::Vec3 const &c);
        /**
         * \brief Writes a set of unit cell lengths to the current file
         *
         * \param a Length of the first unit cell vector in Angstroms
         * \param b Length of the second unit cell vector in Angstroms
         * \param c Length of the third unit cell vector in Angstroms
         */
        void setCellLengths(double a, double b, double c);
        /**
         * \brief Write cell lengths to the current file
         *
         * \param a Length of the first unit cell vector in nm
         * \param b Length of the second unit cell vector in nm
         * \param c Length of the third unit cell vector in nm
         *
         * This function is provided for convenience to work with programs that
         * use the kJ-nm-s unit system (e.g., OpenMM)
         */
        void setCellLengthsNm(double a, double b, double c);
        /**
         * \brief Write the cell angles to the current file
         *
         * \param alpha Angle between the 2nd and 3rd box vectors (in degrees)
         * \param beta Angle between the 1st and 3rd box vectors (in degrees)
         * \param gama Angle between the 1st and 2nd box vectors (in degrees)
         */
        void setCellAngles(double a, double b, double c);
        /**
         * \brief Write the cell angles to the current file
         *
         * \param alpha Angle between the 2nd and 3rd box vectors (in radians)
         * \param beta Angle between the 1st and 3rd box vectors (in radians)
         * \param gama Angle between the 1st and 2nd box vectors (in radians)
         *
         * This function is provided for convenience to work with programs that
         * specify angles in radians
         */
        void setCellAnglesRad(double a, double b, double c);
        /**
         * \brief Write the time to the current file
         *
         * \param time The time to write, in ps
         */
        void setTime(double time);
        /**
         * \brief Write the temperature to the current file (REMD)
         *
         * \param temp The temperature to write (in Kelvin) for T-REMD
         *             simulations (or the pH for pH-REMD simulations)
         */
        void setTemp(double temp);
        /**
         * \brief Sets information about the replica exchange dimension types.
         *
         * \param remdTypes The type of REMD move in each dimension. 1 for
         *      temperature REMD, 3 for Hamiltonian REMD. Should have the same
         *      length as the REMD dimensionality
         *
         * This should only be called once for each file.
         */
        void setRemdTypes(std::vector<int> remdTypes);
        /**
         * \brief Sets the multi-D REMD indices to the current file
         *
         * \param remdIndices The multi-D REMD indices. Length must be as large
         *                    as the REMD-dimension variable
         */
        void setRemdIndices(std::vector<int> remdIndices);
        /**
         * \brief Closes the file
         */
        void close(void);
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
        size_t num_frames_, natom_, natom3_, remd_dimension_, label_;

        // Counters
        size_t coordinate_frame_, velocity_frame_, cell_length_frame_,
               cell_angle_frame_, force_frame_, time_frame_, temp0_frame_,
               remd_indices_frame_;
        bool remd_types_set_;
        std::string program_, programVersion_, application_, title_;
};

};
#endif /* NETCDFFILE_H */
