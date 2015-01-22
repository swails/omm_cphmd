/* NetCDFFile.cpp -- contains the NetCDF file parsing routines for restarts and
 * trajectories
 */

#include <iostream>
#include <sstream>

#include "amber/amber_constants.h"
#include "amber/exceptions.h"
#include "amber/NetCDFFile.h"
#include "amber/unitcell.h"
#include "amber/version.h"

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
        num_frames_(0), natom_(0), natom3_(0), remd_dimension_(0), label_(0),
        coordinate_frame_(0), velocity_frame_(0), cell_length_frame_(0),
        cell_angle_frame_(0), force_frame_(0), time_frame_(0), temp0_frame_(0),
        remd_indices_frame_(0), remd_types_set_(false) {}

AmberNetCDFFile::~AmberNetCDFFile(void) {
    if (is_open_) nc_close(ncid_);
}

void AmberNetCDFFile::readFile(string const& filename) {
    // This can only be called once (and not with writeFile)
    if (ncid_ != -1)
        throw AmberCrdError("AmberNetCDFFile already in use!");
    // Get all of the attributes
    if (nc_open(filename.c_str(), NC_NOWRITE, &ncid_) != NC_NOERR) {
        stringstream iss;
        iss << "Could not open " << filename << " for reading.";
        throw NotNetcdf(iss.str().c_str());
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
        if (conventions != "AMBER") {
            stringstream iss;
            iss << "Conventions " << conventions << " does not match expected "
                << "AMBER";
            throw AmberCrdError(iss.str().c_str());
        }
    } else {
        throw InternalError("Should not be here");
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
    natom3_ = natom_ * 3;
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

vector<OpenMM::Vec3> *AmberNetCDFFile::getCoordinates(int frame) const {
    // Basic error checking
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get coordinates from a new NetCDF file");
    if (!hasCoordinates())
        throw AmberCrdError("NetCDF file does not contain coordinates");
    // Handle restart and trajectory files differently
    if (type_ == RESTART) {
        if (frame != 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range for NetCDF restart";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {0, 0};
        size_t count[] = {natom_, 3};
        double *coords = new double[natom3_];
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
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3(coords[i], coords[i+1], coords[i+2]));
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
        size_t count[] = {1, natom_, 3};
        float *coords = new float[natom3_];
        if (nc_get_vara_float(ncid_, coordinatesVID_, start, count, coords) != NC_NOERR) {
            delete[] coords;
            throw AmberCrdError("Could not get coordinates from NetCDF trajectory file");
        }
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3((double) coords[i  ],
                                        (double) coords[i+1],
                                        (double) coords[i+2]));
        }
        delete[] coords;
        return ret;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
    throw AmberCrdError("Should not be here");
}

vector<OpenMM::Vec3> *AmberNetCDFFile::getVelocities(int frame) const {
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get velocities from a new NetCDF file");
    if (!hasVelocities())
        throw AmberCrdError("NetCDF file does not contain velocities");
    if (type_ == RESTART) {
        if (frame != 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range for NetCDF restart";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {0, 0};
        size_t count[] = {natom_, 3};
        double *vels = new double[natom3_];
        if (nc_get_vara_double(ncid_, velocitiesVID_, start, count, vels) != NC_NOERR) {
            delete[] vels;
            throw AmberCrdError("Could not get velocities from NetCDF restart file");
        }
        string units = GetAttributeText_(velocitiesVID_, "units");
        if (units.empty()) {
            cerr << "WARNING: No units attached to velocities" << endl;
        } else if (units != "angstrom/picosecond") {
            cerr << "WARNING: Coordinate units (" << units
                 << ") not angstroms/picosecond" << endl;
        }
        double scale = GetAttributeFloat_(velocitiesVID_, "scale_factor", 1.0);
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3(vels[i  ] * scale, vels[i+1] * scale,
                                        vels[i+2] * scale));
        }
        delete[] vels;
        return ret;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_) {
            stringstream iss;
            iss << "Frame " << frame << " is out of range of the total number"
                << "of frames (" << num_frames_ << ")";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0, 0};
        size_t count[] = {1, natom_, 3};
        float *vels = new float[natom3_];
        if (nc_get_vara_float(ncid_, velocitiesVID_, start, count, vels) != NC_NOERR) {
            delete[] vels;
            throw AmberCrdError("Could not get coordinates from NetCDF trajectory file");
        }
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        double scale = GetAttributeFloat_(velocitiesVID_, "scale_factor", 1.0);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3((double) vels[i  ] * scale,
                                        (double) vels[i+1] * scale,
                                        (double) vels[i+2] * scale));
        }
        delete[] vels;
        return ret;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
    throw AmberCrdError("Should not be here");
}

vector<OpenMM::Vec3> *AmberNetCDFFile::getForces(int frame) const {
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get forces from a new NetCDF file");
    if (!hasForces())
        throw AmberCrdError("NetCDF file does not contain forces");
    if (type_ == RESTART) {
        if (frame != 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range for NetCDF restart";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {0, 0};
        size_t count[] = {natom_, 3};
        double *frcs = new double[natom3_];
        if (nc_get_vara_double(ncid_, forcesVID_, start, count, frcs) != NC_NOERR) {
            delete[] frcs;
            throw AmberCrdError("Could not get forces from NetCDF restart file");
        }
        string units = GetAttributeText_(forcesVID_, "units");
        if (units.empty()) {
            cerr << "WARNING: No units attached to coordinates" << endl;
        } else if (units != "kilocalorie/mole/angstrom") {
            cerr << "WARNING: Force units (" << units
                 << ") not kilocalorie/mole/angstrom" << endl;
        }
        double scale = GetAttributeFloat_(forcesVID_, "scale_factor", 1.0);
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3(frcs[i  ] * scale, frcs[i+1] * scale,
                                        frcs[i+2] * scale));
        }
        delete[] frcs;
        return ret;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_) {
            stringstream iss;
            iss << "Frame " << frame << " is out of range of the total number"
                << "of frames (" << num_frames_ << ")";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0, 0};
        size_t count[] = {1, natom_, 3};
        float *frcs = new float[natom3_];
        if (nc_get_vara_float(ncid_, forcesVID_, start, count, frcs) != NC_NOERR) {
            delete[] frcs;
            throw AmberCrdError("Could not get coordinates from NetCDF trajectory file");
        }
        vector<OpenMM::Vec3> *ret = new vector<OpenMM::Vec3>;
        ret->reserve(natom_);
        double scale = GetAttributeFloat_(forcesVID_, "scale_factor", 1.0);
        for (int i = 0; i < natom3_; i+=3) {
            ret->push_back(OpenMM::Vec3((double) frcs[i  ] * scale,
                                        (double) frcs[i+1] * scale,
                                        (double) frcs[i+2] * scale));
        }
        delete[] frcs;
        return ret;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
    throw InternalError("Should not be here");
}

OpenMM::Vec3 AmberNetCDFFile::getCellLengths(int frame) const {
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get cell lengths from a new NetCDF file");
    if (cell_lengthsVID_ == -1)
        throw AmberCrdError("NetCDF file does not contain cell lengths");
    if (type_ == RESTART) {
        if (frame != 0)
            throw AmberCrdError("Restart files only have a single frame");
        size_t start[] = {0};
        size_t count[] = {3};
        double box[3];
        if (nc_get_vara_double(ncid_, cell_lengthsVID_, start, count, box) != NC_NOERR)
            throw AmberCrdError("Could not get cell_lengths from NetCDF file");
        return OpenMM::Vec3(box[0], box[1], box[2]);
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_ || frame < 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range; only " << num_frames_
                << "are present.";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0};
        size_t count[] = {1, 3};
        float box[3];
        if (nc_get_vara_float(ncid_, cell_lengthsVID_, start, count, box) != NC_NOERR)
            throw AmberCrdError("Could not get cell_lengths form NetCDF file");
        return OpenMM::Vec3((double)box[0], (double)box[1], (double)box[2]);
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
}

OpenMM::Vec3 AmberNetCDFFile::getCellAngles(int frame) const {
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get cell lengths from a new/nonexistent NetCDF file");
    if (cell_anglesVID_ == -1)
        throw AmberCrdError("NetCDF file does not contain cell lengths");
    if (type_ == RESTART) {
        if (frame != 0)
            throw AmberCrdError("Restart files only have a single frame");
        size_t start[] = {0};
        size_t count[] = {3};
        double box[3];
        if (nc_get_vara_double(ncid_, cell_anglesVID_, start, count, box) != NC_NOERR)
            throw AmberCrdError("Could not get cell_lengths from NetCDF file");
        return OpenMM::Vec3(box[0], box[1], box[2]);
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_ || frame < 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range; only " << num_frames_
                << "are present.";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0};
        size_t count[] = {1, 3};
        float box[3];
        if (nc_get_vara_float(ncid_, cell_anglesVID_, start, count, box) != NC_NOERR)
            throw AmberCrdError("Could not get cell_lengths form NetCDF file");
        return OpenMM::Vec3((double)box[0], (double)box[1], (double)box[2]);
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
}

double AmberNetCDFFile::getTime(int frame) const {
    // Basic error checking
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get time from a new NetCDF file");
    if (timeVID_ == -1)
        throw AmberCrdError("No time present in this NetCDF file");
    if (type_ == RESTART) {
        if (frame != 0)
            throw AmberCrdError("Restart files only have a single frame");
        double time;
        if (nc_get_var_double(ncid_, timeVID_, &time) != NC_NOERR)
            throw AmberCrdError("Could not get time from NetCDF file");
        return time;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_ || frame < 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range; only " << num_frames_
                << "are present.";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0};
        size_t count[] = {1, 1};
        float time;
        if (nc_get_vara_float(ncid_, cell_anglesVID_, start, count, &time) != NC_NOERR)
            throw AmberCrdError("Could not get time form NetCDF file");
        return (double)time;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
}

double AmberNetCDFFile::getTemp(int frame) const {
    // Basic error checking
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get time from a new NetCDF file");
    if (temp0VID_ == -1)
        throw AmberCrdError("No time present in this NetCDF file");
    if (type_ == RESTART) {
        if (frame != 0)
            throw AmberCrdError("Restart files only have a single frame");
        double temp;
        if (nc_get_var_double(ncid_, temp0VID_, &temp) != NC_NOERR)
            throw AmberCrdError("Could not get time from NetCDF file");
        return temp;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_ || frame < 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range; only " << num_frames_
                << "are present.";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0};
        size_t count[] = {1, 1};
        double temp;
        if (nc_get_vara_double(ncid_, temp0VID_, start, count, &temp) != NC_NOERR)
            throw AmberCrdError("Could not get time form NetCDF file");
        return temp;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
}

vector<int> AmberNetCDFFile::getRemdIndices(int frame) const {
    // Basic error checking
    if (!is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot get REMD indices from a new NetCDF file");
    if (remd_indicesVID_ == -1)
        throw AmberCrdError("NetCDF file does not contain REMD indices");
    // Handle restart and trajectory files differently
    if (type_ == RESTART) {
        if (frame != 0) {
            stringstream iss;
            iss << "Frame " << frame << " out of range for NetCDF restart";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {0};
        size_t count[] = {remd_dimension_};
        int *ind = new int[remd_dimension_];
        if (nc_get_vara_int(ncid_, remd_indicesVID_, start, count, ind) != NC_NOERR) {
            delete[] ind;
            throw AmberCrdError("Could not get REMD indices from NetCDF restart file");
        }
        vector<int> ret;
        ret.reserve(remd_dimension_);
        for (int i = 0; i < natom3_; i+=3) {
            ret.push_back(ind[i]);
        }
        delete[] ind;
        return ret;
    } else if (type_ == TRAJECTORY) {
        if (frame >= num_frames_) {
            stringstream iss;
            iss << "Frame " << frame << " is out of range of the total number"
                << "of frames (" << num_frames_ << ")";
            throw AmberCrdError(iss.str().c_str());
        }
        size_t start[] = {frame, 0};
        size_t count[] = {1, remd_dimension_};
        int *ind = new int[remd_dimension_];
        if (nc_get_vara_int(ncid_, remd_indicesVID_, start, count, ind) != NC_NOERR) {
            delete[] ind;
            throw AmberCrdError("Could not get REMD indices from NetCDF file");
        }
        vector<int> ret;
        ret.reserve(remd_dimension_);
        for (int i = 0; i < natom3_; i+=3) {
            ret.push_back(ind[i]);
        }
        delete[] ind;
        return ret;
    } else {
        throw AmberCrdError("Unrecognized NetCDF file type");
    }
    throw AmberCrdError("Should not be here");
}

void AmberNetCDFFile::writeFile(const char* filename, int natom, bool hasCrd,
                bool hasVel, bool hasFrc, bool hasBox, bool hasRemd,
                int remdDimension, const char* title, const char* application) {
    writeFile(string(filename), natom, hasCrd, hasVel, hasFrc, hasBox,
              hasRemd, remdDimension, string(title), string(application));
}

void AmberNetCDFFile::writeFile(string const& filename, int natom, bool hasCrd,
                bool hasVel, bool hasFrc, bool hasBox, bool hasRemd,
                int remdDimension, string const& title, string const& application) {

    // Check and clean the string input
    if (filename.empty())
        throw AmberCrdError("Cannot write to blank file");
    string myTitle, myApp;
    if (title.empty()) {
        stringstream iss;
        iss << "NetCDF file created by " << AMBER_LIBRARY_NAME
            << " v." << AMBER_LIBRARY_VERSION;
        myTitle = iss.str();
    } else {
        myTitle = title;
    }
    if (application.empty())
        myApp = AMBER_LIBRARY_NAME;
    else
        myApp = application;

    // Error checking
    if (type_ == RESTART && !hasCrd)
        throw AmberCrdError("RESTART files MUST have coordinates");
    if (remdDimension > 0 && !hasRemd)
        throw AmberCrdError("remdDimension > 0 requires hasRemd");

    int NDIM;
    nc_type data_type;
    switch (type_) {
        case TRAJECTORY:
            NDIM = 3;
            data_type = NC_FLOAT;
            break;
        case RESTART:
            NDIM = 2;
            data_type = NC_DOUBLE;
            break;
        case AUTOMATIC:
            throw AmberCrdError("Cannot write to an AUTOMATIC NetCDF file");
    }

    if (type_ == AUTOMATIC)
        throw AmberCrdError("Cannot write to an AUTOMATIC file type");

    // Do not write to a file that is already open for reading
    if (ncid_ != -1)
        throw AmberCrdError("AmberNetCDFFile instance already initialized");

    // Create and open the file
    if (nc_create(filename.c_str(), NC_64BIT_OFFSET, &ncid_) != NC_NOERR) {
        stringstream iss;
        iss << "Could not open " << filename << " for writing";
        throw AmberCrdError(iss.str());
    }
    is_open_ = true;

    // Define the global attributes

    string program = string(AMBER_LIBRARY_NAME);
    string programVersion = string(AMBER_LIBRARY_VERSION);
    string conventions = "AMBER";
    if (type_ == RESTART) conventions = "AMBERRESTART";
    if (nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", conventions.size(),
                        conventions.c_str()) != NC_NOERR)
        throw AmberCrdError("Error writing Conventions attribute");
    if (nc_put_att_text(ncid_, NC_GLOBAL, "ConventionVersion", 3, "1.0") != NC_NOERR)
        throw AmberCrdError("Error writing ConventionVersion attribute");
    if (nc_put_att_text(ncid_, NC_GLOBAL, "title",
                        myTitle.size(), myTitle.c_str()) != NC_NOERR)
        throw AmberCrdError("Error writing title attribute");
    if (nc_put_att_text(ncid_, NC_GLOBAL, "program",
                        program.size(), program.c_str()) != NC_NOERR)
        throw AmberCrdError("Error writing program attribute");
    if (nc_put_att_text(ncid_, NC_GLOBAL, "programVersion", programVersion.size(),
                        programVersion.c_str()) != NC_NOERR)
        throw AmberCrdError("Error writing programVersion attribute");
    if (nc_put_att_text(ncid_, NC_GLOBAL, "application",
                        myApp.size(), myApp.c_str()) != NC_NOERR)
        throw AmberCrdError("Error writing application attribute");

    // Define dimensions

    natom_ = (size_t) natom;
    natom3_ = natom_ * 3;
    num_frames_ = 0;

    if (type_ == TRAJECTORY) {
        if (nc_def_dim(ncid_, "frame", NC_UNLIMITED, &frameDID_) != NC_NOERR)
            throw AmberCrdError("Error defining frame dimension");
    }
    if (nc_def_dim(ncid_, "atom", natom_, &atomDID_) != NC_NOERR)
        throw AmberCrdError("Error defining atom dimension");
    if (nc_def_dim(ncid_, "spatial", 3, &spatialDID_) != NC_NOERR)
        throw AmberCrdError("Error defining spatial dimension");
    if (remdDimension > 0) {
        if (nc_def_dim(ncid_, "remd_dimension", (size_t)remdDimension, &remdDID_) != NC_NOERR)
            throw AmberCrdError("Error defining REMD dimension");
    }
    if (hasBox) {
        if (nc_def_dim(ncid_, "cell_spatial", 3, &cell_spatialDID_) != NC_NOERR)
            throw AmberCrdError("Error defining cell_spatial dimension");
        if (nc_def_dim(ncid_, "cell_angular", 3, &cell_angularDID_) != NC_NOERR)
            throw AmberCrdError("Error defining cell_angular dimension");
        if (nc_def_dim(ncid_, "label", 5, &labelDID_) != NC_NOERR)
            throw AmberCrdError("Error defining label dimension");
    }

    // Define variables

    int dimensionID[NC_MAX_VAR_DIMS];
    if (type_ == TRAJECTORY) dimensionID[0] = frameDID_;

    // time variable
    if (nc_def_var(ncid_, "time", data_type, NDIM-2, dimensionID, &timeVID_) != NC_NOERR)
        throw AmberCrdError("Error defining the time variable");
    if (nc_put_att_text(ncid_, timeVID_, "units", 10, "picosecond") != NC_NOERR)
        throw AmberCrdError("Error defining the time unit attribute");

    // Coordinate, velocity, and force variables
    if (type_ == RESTART) {
        dimensionID[0] = atomDID_;
        dimensionID[1] = spatialDID_;
    } else if (type_ == TRAJECTORY) {
        dimensionID[0] = frameDID_;
        dimensionID[1] = atomDID_;
        dimensionID[2] = spatialDID_;
    } else {
        throw InternalError("Should not be here");
    }
    if (hasCrd) {
        if (nc_def_var(ncid_, "coordinates", data_type, NDIM,
                       dimensionID, &coordinatesVID_) != NC_NOERR)
            throw AmberCrdError("Error defining coordinates variable");
        if (nc_put_att_text(ncid_, coordinatesVID_, "units", 8, "angstrom") != NC_NOERR)
            throw AmberCrdError("Error defining coordinates units attribute");
    }
    if (hasVel) {
        if (nc_def_var(ncid_, "velocities", data_type, NDIM,
                       dimensionID, &velocitiesVID_) != NC_NOERR)
            throw AmberCrdError("Error defining velocities variable");
        if (nc_put_att_text(ncid_, velocitiesVID_, "units", 19,
                            "angstrom/picosecond") != NC_NOERR)
            throw AmberCrdError("Error defining velocities units attribute");
        double scale_factor = 1;
        if (nc_put_att_double(ncid_, velocitiesVID_, "scale_factor", NC_DOUBLE,
                              1, &scale_factor) != NC_NOERR)
            throw AmberCrdError("Error defining velocities scale_factor");
    }
    if (hasFrc) {
        if (nc_def_var(ncid_, "forces", data_type, NDIM,
                       dimensionID, &forcesVID_) != NC_NOERR)
            throw AmberCrdError("Error defining forces variable");
        if (nc_put_att_text(ncid_, forcesVID_, "units", 25,
                            "kilocalorie/mole/angstrom") != NC_NOERR)
            throw AmberCrdError("Error defining forces units attribute");
    }

    // spatial variable
    if (nc_def_var(ncid_, "spatial", NC_CHAR, 1, dimensionID, &spatialVID_) != NC_NOERR)
        throw AmberCrdError("Error defining spatial variable");

    // Unit cell-related variable
    if (hasBox) {
        // Labels
        dimensionID[0] = cell_spatialDID_;
        if (nc_def_var(ncid_, "cell_spatial", NC_CHAR, 1,
                       dimensionID, &cell_spatialVID_) != NC_NOERR)
            throw AmberCrdError("Error defining cell_spatial variable");
        dimensionID[0] = cell_angularDID_;
        dimensionID[1] = labelDID_;
        if (nc_def_var(ncid_, "cell_angular", NC_CHAR, 2, dimensionID,
                       &cell_angularVID_) != NC_NOERR)
            throw AmberCrdError("Error defining the cell_angular variable");

        // Cell dimensions
        size_t boxDim;
        if (type_ == RESTART) {
            boxDim = 0;
        } else if (type_ == TRAJECTORY) {
            dimensionID[0] = frameDID_;
            boxDim = 1;
        } else {
            throw InternalError("Should not be here");
        }
        dimensionID[boxDim] = cell_spatialDID_;
        if (nc_def_var(ncid_, "cell_lengths", NC_DOUBLE, NDIM-1,
                       dimensionID, &cell_lengthsVID_) != NC_NOERR)
            throw AmberCrdError("Error defining cell_lengths variable");
        if (nc_put_att_text(ncid_, cell_lengthsVID_,
                            "units", 8, "angstroms") != NC_NOERR)
            throw AmberCrdError("Error defining cell_lengths units");
        dimensionID[boxDim] = cell_angularDID_;
        if (nc_def_var(ncid_, "cell_angles", NC_DOUBLE, NDIM-1,
                       dimensionID, &cell_anglesVID_) != NC_NOERR)
            throw AmberCrdError("Error defining cell_angles variable");
        if (nc_put_att_text(ncid_, cell_anglesVID_, "units", 6, "degree") != NC_NOERR)
            throw AmberCrdError("Error defining degree variable");
    }

    if (hasRemd) {
        if (remdDimension > 0) {
            dimensionID[0] = remdDID_;
            if (nc_def_var(ncid_, "remd_dimtype", NC_INT, 1, dimensionID,
                           &remd_dimtypeVID_) != NC_NOERR)
                throw AmberCrdError("Error creating remd_dimtype variable");
            if (type_ == TRAJECTORY) {
                dimensionID[0] = frameDID_;
                dimensionID[1] = remdDID_;
            } else if (type_ == RESTART) {
                dimensionID[0] = remdDID_;
            } else {
                throw InternalError("Should not be here");
            }
            if (nc_def_var(ncid_, "remd_indices", NC_INT, NDIM-1,
                           dimensionID, &remd_indicesVID_) != NC_NOERR)
                throw AmberCrdError("Error creating remd_indices variable");
        } else {
            // Define temperature
            dimensionID[0] = frameDID_;
            if (nc_def_var(ncid_, "temp0", NC_DOUBLE, NDIM-2,
                           dimensionID, &temp0VID_) != NC_NOERR)
                throw AmberCrdError("Error defining temp0 variable");
            if (nc_put_att_text(ncid_, temp0VID_, "units", 6, "kelvin") != NC_NOERR)
                throw AmberCrdError("Error defining temp0 units");
        }
    }

    // We are done with all definitions

    if (nc_enddef(ncid_) != NC_NOERR)
        throw AmberCrdError("Error taking NetCDF file out of define mode");

    // Fill the variables we already know the values of

    size_t start[3], count[3];

    // spatial variable
    start[0] = 0; count[0] = 3;
    char str[] = {'x', 'y', 'z'};
    dimensionID[0] = spatialDID_;
    if (nc_put_vara_text(ncid_, spatialVID_, start, count, str) != NC_NOERR)
        throw AmberCrdError("Error filling spatial variable");

    if (hasBox) {
        str[0] = 'a'; str[1] = 'b'; str[2] = 'c';
        if (nc_put_vara_text(ncid_, cell_spatialVID_, start, count, str) != NC_NOERR)
            throw AmberCrdError("Error filling cell_spatial variable");
        char labels[] = {'a', 'l', 'p', 'h', 'a',
                         'b', 'e', 't', 'a', ' ',
                         'g', 'a', 'm', 'm', 'a'};
        start[0] = 0; start[1] = 0;
        count[0] = 3; count[1] = 5;
        if (nc_put_vara_text(ncid_, cell_angularVID_, start, count, labels) != NC_NOERR)
            throw AmberCrdError("Error filling cell_angular variable");
    }
}

void AmberNetCDFFile::setCoordinates(vector<OpenMM::Vec3> const &coordinates) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set coordinates on an old file");
    if (coordinates.size() != natom_)
        throw AmberCrdError("Wrong number of coordinates");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {coordinate_frame_, 0, 0};
                size_t count[] = {1, natom_, 3};
                float coords[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    coords[i3  ] = (float) coordinates[i][0];
                    coords[i3+1] = (float) coordinates[i][1];
                    coords[i3+2] = (float) coordinates[i][2];
                }
                if (nc_put_vara_float(ncid_, coordinatesVID_,
                                      start, count, coords) != NC_NOERR)
                    throw AmberCrdError("Error writing coordinates to NetCDF file");
            }
            break;
        case RESTART:
            if (coordinate_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0, 0};
                size_t count[] = {natom_, 3};
                double coords[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    coords[i3  ] = coordinates[i][0];
                    coords[i3+1] = coordinates[i][1];
                    coords[i3+2] = coordinates[i][2];
                }
                if (nc_put_vara_double(ncid_, coordinatesVID_,
                                       start, count, coords) != NC_NOERR)
                    throw AmberCrdError("Error writing coordinates to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    coordinate_frame_++;
}

void AmberNetCDFFile::setCoordinatesNm(vector<OpenMM::Vec3> const &coordinates) {
    vector<OpenMM::Vec3> coords = coordinates;
    for (size_t i = 0; i < coords.size(); i++)
        coords[i] *= ANGSTROM_PER_NANOMETER;
    setCoordinates(coords);
}

void AmberNetCDFFile::setVelocities(vector<OpenMM::Vec3> const &velocities) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set velocities on an old file");
    if (velocities.size() != natom_)
        throw AmberCrdError("Wrong number of velocities");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {velocity_frame_, 0, 0};
                size_t count[] = {1, natom_, 3};
                float vels[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    vels[i3  ] = (float) velocities[i][0];
                    vels[i3+1] = (float) velocities[i][1];
                    vels[i3+2] = (float) velocities[i][2];
                }
                if (nc_put_vara_float(ncid_, velocitiesVID_,
                                      start, count, vels) != NC_NOERR)
                    throw AmberCrdError("Error writing velocities to NetCDF file");
            }
            break;
        case RESTART:
            if (velocity_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0, 0};
                size_t count[] = {natom_, 3};
                double vels[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    vels[i3  ] = velocities[i][0];
                    vels[i3+1] = velocities[i][1];
                    vels[i3+2] = velocities[i][2];
                }
                if (nc_put_vara_double(ncid_, velocitiesVID_,
                                       start, count, vels) != NC_NOERR)
                    throw AmberCrdError("Error writing velocities to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    velocity_frame_++;
}

void AmberNetCDFFile::setVelocitiesNmPerPs(vector<OpenMM::Vec3> const &velocities) {
    vector<OpenMM::Vec3> vels = velocities;
    for (size_t i = 0; i < vels.size(); i++)
        vels[i] *= ANGSTROM_PER_NANOMETER;
    setVelocities(vels);
}

void AmberNetCDFFile::setForces(vector<OpenMM::Vec3> const &forces) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set forces on an old file");
    if (forces.size() != natom_)
        throw AmberCrdError("Wrong number of forces");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {force_frame_, 0, 0};
                size_t count[] = {1, natom_, 3};
                float frcs[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    frcs[i3  ] = (float) forces[i][0];
                    frcs[i3+1] = (float) forces[i][1];
                    frcs[i3+2] = (float) forces[i][2];
                }
                if (nc_put_vara_float(ncid_, forcesVID_,
                                      start, count, frcs) != NC_NOERR)
                    throw AmberCrdError("Error writing forces to NetCDF file");
            }
            break;
        case RESTART:
            if (force_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0, 0};
                size_t count[] = {natom_, 3};
                double frcs[natom3_];
                for (size_t i = 0; i < natom_; i++) {
                    size_t i3 = i*3;
                    frcs[i3  ] = forces[i][0];
                    frcs[i3+1] = forces[i][1];
                    frcs[i3+2] = forces[i][2];
                }
                if (nc_put_vara_double(ncid_, forcesVID_,
                                       start, count, frcs) != NC_NOERR)
                    throw AmberCrdError("Error writing forces to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    force_frame_++;
}

void AmberNetCDFFile::setForcesKJPerNm(vector<OpenMM::Vec3> const &forces) {
    vector<OpenMM::Vec3> frcs = forces;
    for (size_t i = 0; i < frcs.size(); i++)
        frcs[i] *= CALORIE_PER_JOULE * NANOMETER_PER_ANGSTROM;
    setForces(frcs);
}

void AmberNetCDFFile::setUnitCell(OpenMM::Vec3 const &a, OpenMM::Vec3 const &b,
                                  OpenMM::Vec3 const &c) {
    UnitCell cell(a, b, c);
    setCellLengths(cell.getLengthA(), cell.getLengthB(), cell.getLengthC());
    setCellAngles(cell.getAlpha(), cell.getBeta(), cell.getGamma());
}

void AmberNetCDFFile::setUnitCellNm(OpenMM::Vec3 const &a, OpenMM::Vec3 const &b,
                                    OpenMM::Vec3 const &c) {
    setUnitCell(a*ANGSTROM_PER_NANOMETER,
                b*ANGSTROM_PER_NANOMETER,
                c*ANGSTROM_PER_NANOMETER);
}

void AmberNetCDFFile::setCellLengths(double a, double b, double c) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set cell lengths on an old file");
    if (cell_lengthsVID_ == -1)
        throw AmberCrdError("NetCDF file defined without box; "
                            "cannot write cell lengths");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {cell_length_frame_, 0};
                size_t count[] = {1, 3};
                double lengths[] = {a, b, c};
                if (nc_put_vara_double(ncid_, cell_lengthsVID_,
                                       start, count, lengths) != NC_NOERR)
                    throw AmberCrdError("Error writing cell lengths to NetCDF file");
            }
            break;
        case RESTART:
            if (cell_length_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0};
                size_t count[] = {3};
                double lengths[] = {a, b, c};
                if (nc_put_vara_double(ncid_, cell_lengthsVID_,
                                       start, count, lengths) != NC_NOERR)
                    throw AmberCrdError("Error writing cell lengths to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    cell_length_frame_++;
}

void AmberNetCDFFile::setCellLengthsNm(double a, double b, double c) {
    setCellLengths(a*ANGSTROM_PER_NANOMETER,
                   b*ANGSTROM_PER_NANOMETER,
                   c*ANGSTROM_PER_NANOMETER);
}

void AmberNetCDFFile::setCellAngles(double alpha, double beta, double gama) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set cell angles on an old file");
    if (cell_anglesVID_ == -1)
        throw AmberCrdError("NetCDF file defined without box; "
                            "cannot write cell angles");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {cell_angle_frame_, 0};
                size_t count[] = {1, 3};
                double angles[] = {alpha, beta, gama};
                if (nc_put_vara_double(ncid_, cell_anglesVID_,
                                       start, count, angles) != NC_NOERR)
                    throw AmberCrdError("Error writing cell angles to NetCDF file");
            }
            break;
        case RESTART:
            if (cell_angle_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0};
                size_t count[] = {3};
                double angles[] = {alpha, beta, gama};
                if (nc_put_vara_double(ncid_, cell_anglesVID_,
                                       start, count, angles) != NC_NOERR)
                    throw AmberCrdError("Error writing cell angles to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    cell_angle_frame_++;
}

void AmberNetCDFFile::setCellAnglesRad(double alpha, double beta, double gama) {
    setCellAngles(alpha*DEGREE_PER_RADIAN,
                  beta*DEGREE_PER_RADIAN,
                  gama*DEGREE_PER_RADIAN);
}

void AmberNetCDFFile::setTime(double time) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set time on an old file");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {time_frame_};
                size_t count[] = {1};
                float t = (float) time;
                if (nc_put_vara_float(ncid_, timeVID_, start, count, &t) != NC_NOERR)
                    throw AmberCrdError("Error writing time to NetCDF file");
            }
            break;
        case RESTART:
            if (time_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                if (nc_put_var_double(ncid_, timeVID_, &time) != NC_NOERR)
                    throw AmberCrdError("Error writing time to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    time_frame_++;
}

void AmberNetCDFFile::setTemp(double temp) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set temperature on an old file");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {temp0_frame_};
                size_t count[] = {1};
                if (nc_put_vara_double(ncid_, temp0VID_, start, count, &temp) != NC_NOERR)
                    throw AmberCrdError("Error writing temperature to NetCDF file");
            }
            break;
        case RESTART:
            if (temp0_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                if (nc_put_var_double(ncid_, temp0VID_, &temp) != NC_NOERR)
                    throw AmberCrdError("Error writing temperature to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    temp0_frame_++;
}

void AmberNetCDFFile::setRemdTypes(vector<int> remdTypes) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set REMD types on an old file");
    if (remdTypes.size() != remd_dimension_)
        throw AmberCrdError("REMD types is not the correct size");
    if (remd_types_set_)
        throw AmberCrdError("REMD dimtypes should only be set once");

    size_t start[] = {0};
    size_t count[] = {remd_dimension_};

    int array[remd_dimension_];

    for (int i = 0; i < remd_dimension_; i++)
        array[i] = remdTypes[i];

    if (nc_put_vara_int(ncid_, remd_dimtypeVID_, start, count, array) != NC_NOERR)
        throw AmberCrdError("Error writing REMD dimtypes to NetCDF file.");

    remd_types_set_ = true;
}

void AmberNetCDFFile::setRemdIndices(vector<int> remdIndices) {
    if (is_old_ || ncid_ == -1)
        throw AmberCrdError("Cannot set remd indices on an old file");
    if (remdIndices.size() != remd_dimension_)
        throw AmberCrdError("remdIndices is the wrong size");
    switch (type_) {
        case TRAJECTORY:
            {
                size_t start[] = {remd_indices_frame_, 0};
                size_t count[] = {1, remd_dimension_};
                int indices[remd_dimension_];
                for (int i = 0; i < remd_dimension_; i++)
                    indices[i] = remdIndices[i];
                if (nc_put_vara_int(ncid_, remd_indicesVID_,
                                    start, count, indices) != NC_NOERR)
                    throw AmberCrdError("Error writing REMD indices to NetCDF file");
            }
            break;
        case RESTART:
            if (cell_angle_frame_ > 0)
                throw AmberCrdError("Restart files can only have 1 frame!");
            {
                size_t start[] = {0};
                size_t count[] = {remd_dimension_};
                int indices[remd_dimension_];
                for (int i = 0; i < remd_dimension_; i++)
                    indices[i] = remdIndices[i];
                if (nc_put_vara_int(ncid_, remd_indicesVID_,
                                    start, count, indices) != NC_NOERR)
                    throw AmberCrdError("Error writing REMD indices to NetCDF file");
            }
            break;
        default:
            throw InternalError("Should not be here");
            break;
    }
    remd_indices_frame_++;
}

void AmberNetCDFFile::close(void) {
    if (!is_open_)
        throw AmberCrdError("Cannot close file that is not open");
    if (nc_close(ncid_) != NC_NOERR)
        throw AmberCrdError("Error closing file");
    is_open_ = false;
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
AmberNetCDFFile::GetAttributeFloat_(int varID, const char* attr, double default_) const {
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
        num_frames_(0), natom_(0), natom3_(0), remd_dimension_(0), label_(0),
        coordinate_frame_(0), velocity_frame_(0), cell_length_frame_(0),
        cell_angle_frame_(0), force_frame_(0), time_frame_(0), temp0_frame_(0),
        remd_indices_frame_(0), remd_types_set_(false) {
    throw NotNetcdf("Compiled without NetCDF support. Cannot use NetCDF functionality");
}
#endif /* HAS_NETCDF */
