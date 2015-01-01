/* NetCDFFile.cpp -- contains the NetCDF file parsing routines for restarts and
 * trajectories
 */

#include <iostream>
#include <sstream>

#include "NetCDFFile.h"
#include "exceptions.h"

using namespace std;
using namespace Amber;

#ifdef HAS_NETCDF
#include "netcdf.h"

AmberNetCDFFile::AmberNetCDFFile(AmberNetCDFFile::FileType type) :
        ncid_(-1), atomDID_(-1), frameDID_(-1), spatialDID_(-1),
        cell_spatialDID_(-1), cell_angularDID_(-1), labelDID_(-1),
        remdDID_(-1), spatialVID_(-1), cell_spatialVID_(-1),
        cell_angularVID_(-1), timeVID_(-1), coordinatesVID_(-1),
        cell_lengthsVID_(-1), cell_anglesVID_(-1), velocitiesVID_(-1),
        forcesVID_(-1), temp0VID_(-1), remd_dimtypeVID_(-1),
        remd_indicesVID_(-1), is_open_(false), is_old_(false), type_(type),
        num_frames_(0), natom_(0), remd_dimension_(0), label_(0) {}

AmberNetCDFFile::~AmberNetCDFFile(void) {
    if (is_open_) nc_close(ncid_);
}

void AmberNetCDFFile::readFile(string const& filename) {
    // Get all of the attributes
    if (nc_open(filename.c_str(), NC_NOWRITE, &ncid_) != NC_NOERR) {
        stringstream iss;
        iss << "Could not open " << filename << " for reading.";
        throw AmberCrdError(iss.str().c_str());
    }
    is_open_ = true;
    is_old_ = true;
    string conventions = GetAttributeText_(NC_GLOBAL, "Conventions", true);
    if (type_ == AUTOMATIC) {
        if (conventions == "AMBER")
            type_ = TRAJECTORY;
        else if (conventions == "AMBERRESTART")
            type_ = RESTART;
        else {
            stringstream iss;
            iss << "Unknown Conventions attribute (" << conventions << ")";
            throw AmberCrdError(iss.str().c_str());
        }
    } else if (type_ == RESTART) {
        if (conventions != "AMBERRESTART") {
            stringstream iss;
            iss << "Conventions " << conventions << " does not match expected "
                << "AMBERRESTART";
            throw AmberCrdError(iss.str().c_str());
        }
    } else if (type_ == TRAJECTORY) {
        if (conventions != "TRAJECTORY") {
            stringstream iss;
            iss << "Conventions " << conventions << " does not match expected "
                << "AMBER";
            throw AmberCrdError(iss.str().c_str());
        }
    } else {
        throw AmberCrdError("Should not be here");
    }
    string cv = GetAttributeText_(NC_GLOBAL, "ConventionVersion", true);
    if (cv != "1.0")
        cerr << "WARNING: ConventionVersion (" << cv << ") does not match "
             << "expected \"1.0\"" << endl;
    program_ = GetAttributeText_(NC_GLOBAL, "program");
    programVersion_ = GetAttributeText_(NC_GLOBAL, "programVersion");
    application_ = GetAttributeText_(NC_GLOBAL, "application");
    title_ = GetAttributeText_(NC_GLOBAL, "title");
    // Now get all of the dimensions
    atomDID_ = GetDimensionID_("atom", natom_);
    if (atomDID_ == -1) throw AmberCrdError("Could not get atom dimension");
    if (type_ == TRAJECTORY) {
        frameDID_ = GetDimensionID_("frame", num_frames_);
        if (frameDID_ == -1) throw AmberCrdError("Could not get frame dimension");
    } else {
        num_frames_ = 1;
    }
    remdDID_ = GetDimensionID_("remd_dimension", remd_dimension_);
    size_t value;
    cell_spatialDID_ = GetDimensionID_("cell_spatial", value);
    if (cell_spatialDID_ != -1 && value != 3) {
        cerr << "WARNING: cell_spatial was expected to be 3" << endl;
    }
    cell_angularDID_ = GetDimensionID_("cell_angular", value);
    if (cell_angularDID_ != -1 && value != 3) {
        cerr << "WARNING: cell_angular was expected to be 3" << endl;
    }
    labelDID_ = GetDimensionID_("label", label_);
    spatialDID_ = GetDimensionID_("spatial", value);
    if (spatialDID_ == -1) {
        cerr << "WARNING: spatial dimension not found" << endl;
    } else if (value != 3) {
        stringstream iss;
        iss << "spatial dimension (" << value << ") is not 3";
        throw AmberCrdError(iss.str().c_str());
    }
    // Now get all of the variables
    spatialVID_ = GetVariableID_("spatial");
    cell_spatialVID_ = GetVariableID_("cell_spatial");
    cell_angularVID_ = GetVariableID_("cell_angular");
    timeVID_ = GetVariableID_("time");
    coordinatesVID_ = GetVariableID_("coordinates");
    cell_lengthsVID_ = GetVariableID_("cell_lengths");
    cell_anglesVID_ = GetVariableID_("cell_angles");
    velocitiesVID_ = GetVariableID_("velocities");
    forcesVID_ = GetVariableID_("forces");
    temp0VID_ = GetVariableID_("temp0");
    remd_dimtypeVID_ = GetVariableID_("remd_dimtypeVID_");
}

void AmberNetCDFFile::readFile(const char* filename) {
    readFile(string(filename));
}

vector<OpenMM::Vec3> AmberNetCDFFile::getCoordinates(int frame) const {
    if (!is_old_)
        throw AmberCrdError("Cannot get coordinates from a new NetCDF file");
    if (type_ == RESTART) {
        if (frame != 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range for NetCDF restart";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {0, 0};
        size_t count[] = {(size_t)natom_, 3};
        double *coords = new double[3*natom_];
        if (nc_get_vara_double(ncid_, coordinatesVID_, start, count, coords) != NC_NOERR) {
            delete[] coords;
            throw AmberCrdError("Could not get coordinates from NetCDF restart file");
        }
        string units = GetAttributeText_(coordinatesVID_, "units");
        if (units.empty()) {
            cerr << "WARNING: No units attached to coordinates" << endl;
        } else if (units != "angstrom") {
            cerr << "WARNING: Coordinate units (" << units << ") not angstroms"
                 << endl;
        }
        vector<OpenMM::Vec3> ret;
        ret.reserve(natom_);
        for (int i = 0; i < natom_*3; i+=3) {
            ret.push_back(OpenMM::Vec3(coords[i], coords[i+1], coords[i+2]));
        }
        delete[] coords;
        return ret;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_) {
            stringstream iss;
            iss << "Frame " << frame << " is out of range of the total number"
                << "of frames (" << num_frames_ << ")";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0, 0};
        size_t count[] = {1, (size_t)natom_, 3};
        float *coords = new float[3*natom_];
        if (nc_get_vara_float(ncid_, coordinatesVID_, start, count, coords) != NC_NOERR) {
            delete[] coords;
            throw AmberCrdError("Could not get coordinates from NetCDF trajectory file");
        }
        vector<OpenMM::Vec3> ret;
        ret.reserve(natom_);
        for (int i = 0; i < natom_*3; i+=3) {
            ret.push_back(OpenMM::Vec3(coords[i], coords[i+1], coords[i+2]));
        }
        delete[] coords;
        return ret;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
}

// Private routines
string
AmberNetCDFFile::GetAttributeText_(int varID, const char* attr, bool required) const {
    size_t attr_len;
    string ret;
    if (nc_inq_attlen(ncid_, varID, attr, &attr_len) == NC_NOERR) {
        char *attrtext = new char[attr_len+1];
        if (nc_get_att_text(ncid_, varID, attr, attrtext) != NC_NOERR) {
            stringstream iss;
            iss << "Trouble getting text attribute " << attr;
            throw AmberCrdError(iss.str().c_str());
        }
        attrtext[attr_len] = '\0';
        ret = attrtext;
        delete attrtext;
    } else if (required) {
        stringstream iss;
        iss << "Required attribute " << attr << " was not found";
        throw AmberCrdError(iss.str().c_str());
    }
    return ret;
}

double
AmberNetCDFFile::GetAttributeFloat_(int varID, const char* attr, double default_) {
    double ret = default_;
    nc_type type;
    if (nc_inq_atttype(ncid_, varID, attr, &type) == NC_NOERR) {
        if (type == NC_DOUBLE) {
            if (nc_get_att_double(ncid_, varID, attr, &ret) != NC_NOERR)
                cerr << "WARNING: Could not get double " << attr << endl;
        } else if (type == NC_FLOAT) {
            float tmp;
            if (nc_get_att_float(ncid_, varID, attr, &tmp) != NC_NOERR)
                cerr << "WARNING: Could not get float " << attr << endl;
            else
                ret = (double) tmp;
        } else {
            cerr << "WARNING: Attribute " << attr << " is not a float" << endl;
        }
    }

    return ret;
}

int AmberNetCDFFile::GetDimensionID_(const char* name, size_t &value) {
    int dimID = -1;
    value = 0;
    if (nc_inq_dimid(ncid_, name, &dimID) == NC_NOERR) {
        if (nc_inq_dimlen(ncid_, dimID, &value) != NC_NOERR) {
            stringstream iss;
            iss << "Could not find length of dimension " << name;
            throw AmberCrdError(iss.str().c_str());
        }
    }
    return dimID;
}

int AmberNetCDFFile::GetVariableID_(const char* name) {
    int varID;
    if (nc_inq_varid(ncid_, name, &varID) == NC_NOERR) {
        return varID;
    }
    return -1;
}

#else
AmberNetCDFFile::AmberNetCDFFile(FileType type) :
        ncid_(-1), atomDID_(-1), frameDID_(-1), spatialDID_(-1),
        cell_spatialDID_(-1), cell_angularDID_(-1), labelDID_(-1),
        remdDID_(-1), spatialVID_(-1), cell_spatialVID_(-1),
        cell_angularVID_(-1), timeVID_(-1), coordinatesVID_(-1),
        cell_lengthsVID_(-1), cell_anglesVID_(-1), velocitiesVID_(-1),
        forcesVID_(-1), temp0VID_(-1), remd_dimtypeVID_(-1),
        remd_indicesVID_(-1), is_open_(false), is_old_(false), type_(type),
        num_frames_(0), natom_(0), remd_dimension_(0), label_(0) {}
    throw AmberCrdError("Compiled without NetCDF support. Cannot use NetCDF functionality");
}
#endif /* HAS_NETCDF */
