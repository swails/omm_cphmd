/** ambercrd.cpp -- Contains the functionality to read and write Amber
  * coordinate files
  */
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>

#include "amber/NetCDFFile.h"
#include "amber/amber_constants.h"
#include "amber/ambercrd.h"
#include "amber/readparm.h"
#include "amber/string_manip.h"
#include "OpenMM.h"

using namespace Amber;
using namespace std;
using namespace OpenMM;

void AmberCoordinateFrame::readRst7(string const& filename) {
    if (readNetCDF_(filename) != 0)
        readASCII_(filename);
}

AmberCoordinateFrame::~AmberCoordinateFrame() {
    if (coordinates_ != NULL) delete coordinates_;
    if (velocities_ != NULL) delete velocities_;
}

void AmberCoordinateFrame::readASCII_(string const& filename) {

    ifstream input;
    input.open(filename.c_str());
    string line;

    bool has_velocities = false;
    bool has_box = false;

    // Parse all lines and save in memory
    vector<string> lines;

    while (getline(input, line)) {
        lines.push_back(line);
    }

    if (lines.size() < 2) {
        throw AmberCrdError("Too few lines in inpcrd/restart");
    }

    vector<string> words = split(lines[1]);

    if (words.size() > 2) {
        has_remd_ = true;
        natom_ = StringToDouble(words[0]);
//      double time = StringToDouble(words[1]);
        temp0_ = StringToDouble(words[2]);
    } else if (words.size() == 2) {
        natom_ = StringToDouble(words[0]);
//      double time = StringToDouble(words[1]);
    } else if (words.size() == 1) {
        natom_ = StringToDouble(words[0]);
    } else {
        throw AmberCrdError("Unknown format of Amber inpcrd/restrt");
    }

    if (lines.size() == (natom_+1)/2 + 2) {
    } else if (lines.size() == (natom_+1)/2 + 3) {
        has_box = true;
    } else if (lines.size() == ((natom_+1)/2) * 2 + 2) {
        has_velocities = true;
    } else if (lines.size() == ((natom_+1)/2) * 2 + 3) {
        has_velocities = true;
        has_box = true;
    } else {
        throw AmberCrdError("Unknown number of lines in inpcrd/restrt");
    }

    // Now parse the coordinates
    coordinates_->clear();
    velocities_->clear();
    coordinates_->reserve(natom_);
    if (has_velocities)
        velocities_->reserve(natom_);

    int current_line = 2;
    int atomno = 0;
    for (int i = 0; i < (natom_+1)/2; i++) {
        int cl = current_line + i;
        coordinates_->push_back(OpenMM::Vec3(StringToDouble(lines[cl].substr(0, 12)),
                                             StringToDouble(lines[cl].substr(12, 12)),
                                             StringToDouble(lines[cl].substr(24, 12))));
        if (++atomno < natom_) {
            coordinates_->push_back(OpenMM::Vec3(StringToDouble(lines[cl].substr(36, 12)),
                                                 StringToDouble(lines[cl].substr(48, 12)),
                                                 StringToDouble(lines[cl].substr(60, 12))));
            ++atomno;
        }
    }

    current_line += (natom_+1) / 2;

    if (has_velocities) {
        atomno = 0;
        for (int i = 0; i < (natom_+1)/2; i++) {
            int cl = current_line + i;
            velocities_->push_back(
                    OpenMM::Vec3(StringToDouble(lines[cl].substr(0, 12)) * AMBER_TIME_PER_PS,
                                 StringToDouble(lines[cl].substr(12, 12)) * AMBER_TIME_PER_PS,
                                 StringToDouble(lines[cl].substr(24, 12)) * AMBER_TIME_PER_PS));
            if (++atomno < natom_) {
                velocities_->push_back(
                        OpenMM::Vec3(StringToDouble(lines[cl].substr(36, 12)) * AMBER_TIME_PER_PS,
                                     StringToDouble(lines[cl].substr(48, 12)) * AMBER_TIME_PER_PS,
                                     StringToDouble(lines[cl].substr(60, 12)) * AMBER_TIME_PER_PS));
                ++atomno;
            }
        }
        current_line += (natom_+1) / 2;
    }

    if (has_box) {
        a_ = StringToDouble(lines[current_line].substr(0, 12));
        b_ = StringToDouble(lines[current_line].substr(12, 12));
        c_ = StringToDouble(lines[current_line].substr(24, 12));
        alpha_ = StringToDouble(lines[current_line].substr(36, 12));
        beta_ = StringToDouble(lines[current_line].substr(48, 12));
        gama_ = StringToDouble(lines[current_line].substr(60, 12));
    }
}

int AmberCoordinateFrame::readNetCDF_(string const& filename) {
    AmberNetCDFFile ncfile(AmberNetCDFFile::RESTART);
    try {
        ncfile.readFile(filename);
    } catch (NotNetcdf &e) {
        return 1;
    }
    if (!ncfile.hasCoordinates())
        throw AmberCrdError("NetCDF file does not have coordinates");
    vector<OpenMM::Vec3> *crd;
    vector<OpenMM::Vec3> *vel;

    crd = ncfile.getCoordinates();
    delete coordinates_;
    coordinates_ = crd;
    natom_ = ncfile.getNatom();

    if (ncfile.hasVelocities()) {
        vel = ncfile.getVelocities();
        delete velocities_;
        velocities_ = vel;
    }

    if (ncfile.hasBox()) {
        OpenMM::Vec3 cell_lengths = ncfile.getCellLengths();
        OpenMM::Vec3 cell_angles = ncfile.getCellAngles();
        a_ = cell_lengths[0];
        b_ = cell_lengths[1];
        c_ = cell_lengths[2];
        alpha_ = cell_angles[0];
        beta_ = cell_angles[1];
        gama_ = cell_angles[2];
    }

    return 0;
}

void AmberCoordinateFrame::readRst7(const char* filename) {
    readRst7(string(filename));
}

void AmberCoordinateFrame::writeRst7(string const& filename, bool netcdf) {
    writeRst7(filename.c_str(), netcdf);
}

void AmberCoordinateFrame::writeRst7(const char* filename, bool netcdf) {

    if (netcdf) {
        writeNetCDF_(filename);
    } else {
        writeASCII_(filename);
    }

}

void AmberCoordinateFrame::writeASCII_(const char* filename) {

    FILE *file = NULL;

    file = fopen(filename, "w");

    if (!file) {
        string msg = "Could not open ";
        msg = msg + filename + " for writing";
        throw AmberCrdError(msg.c_str());
    }

    // Write the header
    fprintf(file, "\n%5d\n", natom_);

    for (int i = 0; i < natom_; i++) {
        fprintf(file, "%12.7f%12.7f%12.7f",
                (*coordinates_)[i][0], (*coordinates_)[i][1], (*coordinates_)[i][2]);
        if (i % 2 == 1) fprintf(file, "\n");
    }
    if (natom_ % 2 == 1) fprintf(file, "\n");

    if (velocities_->size() > 0) {
        for (int i = 0; i < natom_; i++) {
            fprintf(file, "%12.7f%12.7f%12.7f",
                    (*velocities_)[i][0]*PS_PER_AMBER_TIME,
                    (*velocities_)[i][1]*PS_PER_AMBER_TIME,
                    (*velocities_)[i][2]*PS_PER_AMBER_TIME);
            if (i % 2 == 1) fprintf(file, "\n");
        }
        if (natom_ % 2 == 1) fprintf(file, "\n");
    }

    if (a_ > 0 && b_ > 0 && c_ > 0 && alpha_ > 0 && beta_ > 0 && gama_ > 0) {
        fprintf(file, "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n",
                a_, b_, c_, alpha_, beta_, gama_);
    }
    fclose(file);
}

void AmberCoordinateFrame::writeNetCDF_(const char* filename) {
#if 0
#ifdef HAS_NETCDF
    int ncid;
    int spatialDID, cell_spatialDID, cell_angularDID, atomDID, labelDID;
    int timeVID, coordinateVID, velocityVID, cell_anglesVID, cell_lengthVID,
        spatialVID, cell_angularVID, cell_spatialVID;

    // Create the file
    if (nc_create(filename, NC_64BIT_OFFSET, &ncid) != NC_NOERR) {
        stringstream iss;
        iss << "Could not open " << filename << " for writing";
        throw AmberCrdError(iss.str().c_str());
    }
    // Put into define mode
    if (nc_redef(ncid) != NC_NOERR)
        throw AmberCrdError("Could not put NetCDF file into define mode");
    // Define the dimensions
    if (nc_def_dim(ncid, "atom", (size_t)natom_, &atomDID) != NC_NOERR) {
        nc_close(ncid);
        throw AmberCrdError("Could not define atom dimension");
    }
    if (nc_def_dim(ncid, "spatial", 3, &spatialDID) != NC_NOERR) {
        nc_close(ncid);
        throw AmberCrdError("Could not define spatial dimension");
    }
    if (a_ > 0 && b_ > 0 && c_ > 0 && alpha_ > 0 && beta_ > 0 && gama_ > 0) {
        if (nc_def_dim(ncid, "cell_spatial", 3, &cell_spatialDID) != NC_NOERR) {
            nc_close(ncid);
            throw AmberCrdError("Could not define cell_spatial dimension");
        }
        if (nc_def_dim(ncid, "cell_angular", 3, &cell_angularDID) != NC_NOERR) {
            nc_close(ncid);
            throw AmberCrdError("Could not define cell_angular dimension");
        }
        if (nc_def_dim(ncid, "label", 5, &labelDID) != NC_NOERR) {
            nc_close(ncid);
            throw AmberCrdError("Could not define label dimension");
        }
    }
    // Define the variables
    int dimensionID[NC_MAX_VAR_DIMS];
    if (nc_def_var(ncid, "time", NC_DOUBLE, 0, dimensionID, &timeVID) != NC_NOERR) {
        nc_close(ncid);
        throw AmberCrdError("Could not create the time variable");
    }
    dimensionID[0] = atomDID;
    dimensionID[1] = spatialDID;
    if (nc_def_var(ncid, "coordinates", NC_DOUBLE, 2, dimensionID,
                   &coordinateVID) != NC_NOERR) {
        nc_close(ncid);
        throw AmberCrdError("Could not create the coordinate variable");
    }
    dimensionID[0] = spatialDID;
    if (nc_def_var(ncid, "spatial", NC_CHAR, 1, dimensionID, &spatialVID) != NC_NOERR) {
        nc_close(ncid);
        throw AmberCrdError("Could not create the spatial variable");
    }
    if (a_ > 0 && b_ > 0 && c_ > 0 && alpha_ > 0 && beta_ > 0 && gama_ > 0) {
        dimensionID[0] = cell_angularDID;
        dimensionID[1] = labelDID;
        if (nc_def_var(ncid, "cell_angular", NC_CHAR, 2, dimensionID,
                       &cell_angularVID) != NC_NOERR) {
            nc_close(ncid);
            throw AmberCrdError("Could not create the cell_angular variable");
        }
        dimensionID[0] = cell_spatialDID;
        if (nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, dimensionID,
                       &cell_angularVID) != NC_NOERR) {
            nc_close(ncid);
            throw AmberCrdError("Could not create the cell_spatial variable");
        }
        dimensionID[0] = cell_spatialDID;
    }
#else
    throw AmberCrdError("libcphmd not built with NetCDF support");
#endif /* HAS_NETCDF */
#endif
}
