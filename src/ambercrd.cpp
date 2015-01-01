/** ambercrd.cpp -- Contains the functionality to read and write Amber
  * coordinate files
  */
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>

#include "amber/amber_constants.h"
#include "amber/ambercrd.h"
#include "amber/readparm.h"
#include "amber/string_manip.h"
#include "OpenMM.h"

// TODO: Get rid of this
#ifdef HAS_NETCDF
#include "netcdf.h"
#endif /* HAS_NETCDF */

using namespace Amber;
using namespace std;
using namespace OpenMM;

#ifdef HAS_NETCDF
// Some private NetCDF functions, dump it in the Amber namespace
namespace Amber {
int GetNCDimension(int ncid, const char* name, size_t *len) {
    int dimID;
    if (nc_inq_dimid(ncid, name, &dimID) != NC_NOERR) {
        stringstream iss;
        iss << "Could not find dimension ID for " << name;
        throw AmberCrdError(iss.str().c_str());
    }

    if (nc_inq_dimlen(ncid, dimID, len) != NC_NOERR) {
        stringstream iss;
        iss << "Could not find length of dimension " << name;
        throw AmberCrdError(iss.str().c_str());
    }

    return dimID;
}

};
#endif /* HAS_NETCDF */

void AmberCoordinateFrame::readRst7(string const& filename) {
    if (readNetCDF_(filename) != 0)
        readASCII_(filename);
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
    coordinates_.clear();
    velocities_.clear();
    coordinates_.reserve(natom_);
    if (has_velocities)
        velocities_.reserve(natom_);

    int current_line = 2;
    int atomno = 0;
    for (int i = 0; i < (natom_+1)/2; i++) {
        int cl = current_line + i;
        coordinates_.push_back(OpenMM::Vec3(StringToDouble(lines[cl].substr(0, 12)),
                                            StringToDouble(lines[cl].substr(12, 12)),
                                            StringToDouble(lines[cl].substr(24, 12))));
        if (++atomno < natom_) {
            coordinates_.push_back(OpenMM::Vec3(StringToDouble(lines[cl].substr(36, 12)),
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
            velocities_.push_back(
                    OpenMM::Vec3(StringToDouble(lines[cl].substr(0, 12)) * AMBER_TIME_PER_PS,
                                 StringToDouble(lines[cl].substr(12, 12)) * AMBER_TIME_PER_PS,
                                 StringToDouble(lines[cl].substr(24, 12)) * AMBER_TIME_PER_PS));
            if (++atomno < natom_) {
                velocities_.push_back(
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
#ifdef HAS_NETCDF
    int ncid;
    if (nc_open(filename.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
        return 1; // It is not a NetCDF file

    /* If we got here, it is a NetCDF file. Check that it is a restart file by
     * checking the conventions
     */
    size_t attr_len;
    if (nc_inq_attlen(ncid, NC_GLOBAL, "Conventions", &attr_len) != NC_NOERR)
        throw AmberCrdError("Problem getting length of Conventions NetCDF attribute");

    char *conventions = new char[attr_len+1];
    if (nc_get_att_text(ncid, NC_GLOBAL, "Conventions", conventions) != NC_NOERR)
        throw AmberCrdError("Problem getting Conventions NetCDF attribute");
    conventions[attr_len] = '\0'; // null-terminate

    if (string("AMBERRESTART") != conventions) {
        stringstream iss;
        iss << "NetCDF Conventions (" << conventions << ") is not AMBERRESTART";
        throw AmberCrdError(iss.str().c_str());
    }

    if (nc_inq_attlen(ncid, NC_GLOBAL, "ConventionVersion", &attr_len) != NC_NOERR)
        throw AmberCrdError("Problem getting length of ConventionVersion "
                            "NetCDF attribute");
    delete conventions;
    conventions = new char[attr_len+1];
    if (nc_get_att_text(ncid, NC_GLOBAL, "ConventionVersion", conventions) != NC_NOERR)
        throw AmberCrdError("Problem getting ConventionVersion NetCDF attribute");
    conventions[attr_len] = '\0';

    if (string("1.0") != conventions) {
        stringstream iss;
        iss << "NetCDF ConventionVersion (" << conventions << ") is not 1.0";
        throw AmberCrdError(iss.str().c_str());
    }
    delete conventions;

    // OK, I'm satisfied that it's an Amber restart file. Now check for natom
    int atomDID, spatialDID, angularDID;
    size_t atom, spatial, angular;
    bool has_box = true;
    atomDID = GetNCDimension(ncid, "atom", &atom);
    try {
        spatialDID = GetNCDimension(ncid, "cell_spatial", &spatial);
        angularDID = GetNCDimension(ncid, "cell_angular", &angular);
        if (spatial != 3 || angular != 3) {
            cerr << "WARNING: cell_spatial (" << spatial
                 << ") and cell_angular (" << angular
                 << ") should both be 3. Ignoring box info." << endl;
            has_box = false;
        }
    } catch (AmberCrdError &e) {
        has_box = false;
    }
    natom_ = (int) atom;

    // Now get the coordinates (and make sure units are correct)
    int coordVID;
    if (nc_inq_varid(ncid, "coordinates", &coordVID) != NC_NOERR)
        throw AmberCrdError("NetCDF restart does not have coordinates");
    double *coords = new double[3*natom_];
    size_t start[] = {0, 0};
    size_t count[] = {natom_, 3};
    if (nc_get_vara_double(ncid, coordVID, start, count, coords) != NC_NOERR) {
        delete[] coords;
        throw AmberCrdError("Trouble getting coordinates from NetCDF restart");
    }
    // Make sure the units are set and set to "angstrom"
    if (nc_inq_attlen(ncid, coordVID, "units", &attr_len) != NC_NOERR) {
        cerr << "WARNING: coordinates variable does not have units attribute"
             << endl;
    } else {
        char *attr = new char[attr_len+1];
        if (nc_get_att_text(ncid, coordVID, "units", attr) != NC_NOERR)
            cerr << "WARNING: trouble getting coordinates units attribute"
                 << endl;
        else {
            attr[attr_len] = '\0';
            if (string("angstrom") != attr)
                cerr << "WARNING: coordinates units are " << attr
                     << ", not angstroms" << endl;
        }
        delete attr;
    }
    coordinates_.clear();
    coordinates_.reserve(natom_);
    for (int i = 0; i < natom_*3; i+=3) {
        coordinates_.push_back(
                OpenMM::Vec3(coords[i], coords[i+1], coords[i+2]));
    }

    // Velocities
    int velVID;
    if (nc_inq_varid(ncid, "velocities", &velVID) == NC_NOERR) {
        // Has velocities
        if (nc_get_vara_double(ncid, velVID, start, count, coords) != NC_NOERR) {
            delete[] coords;
            throw AmberCrdError("Trouble getting velocities from NetCDF restart");
        }
        if (nc_inq_attlen(ncid, velVID, "units", &attr_len) != NC_NOERR)
            cerr << "WARNING: velocities variable does not have units attribute"
                 << endl;
        else {
            char *attr = new char[attr_len+1];
            if (nc_get_att_text(ncid, velVID, "units", attr) != NC_NOERR)
                cerr << "WARNING: trouble getting velocities units attribute"
                     << endl;
            else {
                attr[attr_len] = '\0';
                if (string("angstrom/picosecond") != attr)
                    cerr << "WARNING: velocities units are " << attr
                         << ", not angstroms/picoseconds" << endl;
            }
            delete attr;
        }
        // See if we have a scale factor
        nc_type type;
        double scale_factor = 1;
        if (nc_inq_atttype(ncid, velVID, "scale_factor", &type) == NC_NOERR) {
            if (type != NC_DOUBLE && type != NC_FLOAT) {
                cerr << "WARNING: velocities scale_factor is not a double"
                     << endl;
            } else if (type == NC_DOUBLE) {
                if (nc_get_att_double(ncid, velVID, "scale_factor",
                                      &scale_factor) != NC_NOERR) {
                    cerr << "WARNING: could not get velocity scale_factor"
                         << endl;
                }
            } else if (type == NC_FLOAT) {
                float sf;
                if (nc_get_att_float(ncid, velVID, "scale_factor", &sf) != NC_NOERR)
                    cerr << "WARNING: could not get velocity scale_factor"
                         << endl;
                else
                    scale_factor = (double) sf;
            }
        }
        velocities_.clear();
        velocities_.reserve(natom_);
        for (int i = 0; i < natom_*3; i+=3) {
            velocities_.push_back(
                    OpenMM::Vec3(coords[i]*scale_factor,
                                 coords[i+1]*scale_factor,
                                 coords[i+2]*scale_factor)
            );
        }
    }
    // END VELOCITIES
    delete[] coords;

    // Box
    if (has_box) {
        int abcVID, angleVID;
        if (nc_inq_varid(ncid, "cell_lengths", &abcVID) != NC_NOERR) {
            cerr << "WARNING: Could not find expected cell_lengths variable"
                 << endl;
            return 0;
        }

        if (nc_inq_varid(ncid, "cell_angles", &angleVID) != NC_NOERR) {
            cerr << "WARNING: Could not find expected cell_angles variable"
                 << endl;
            return 0;
        }

        size_t start[] = {0};
        size_t count[] = {3};
        double box[3];
        if (nc_get_vara_double(ncid, abcVID, start, count, box) != NC_NOERR) {
            cerr << "WARNING: Could not get cell_lengths from NetCDF file"
                 << endl;
            return 0;
        }
        a_ = box[0]; b_ = box[1]; c_ = box[2];

        if (nc_get_vara_double(ncid, angleVID, start, count, box) != NC_NOERR) {
            cerr << "WARNING: Could not get cell_angles from NetCDF file"
                 << endl;
            return 0;
        }
        alpha_ = box[0]; beta_ = box[1]; gama_ = box[2];

    }
    return 0;
#else /* !HAS_NETCDF */
    return 1; // no NetCDF
#endif /* HAS_NETCDF */
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
                coordinates_[i][0], coordinates_[i][1], coordinates_[i][2]);
        if (i % 2 == 1) fprintf(file, "\n");
    }
    if (natom_ % 2 == 1) fprintf(file, "\n");

    if (velocities_.size() > 0) {
        for (int i = 0; i < natom_; i++) {
            fprintf(file, "%12.7f%12.7f%12.7f",
                    velocities_[i][0]*PS_PER_AMBER_TIME,
                    velocities_[i][1]*PS_PER_AMBER_TIME,
                    velocities_[i][2]*PS_PER_AMBER_TIME);
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
}
