/** ambercrd.cpp -- Contains the functionality to read and write Amber
  * coordinate files
  */
#include <cstdio>
#include <fstream>

#include "amber_constants.h"
#include "ambercrd.h"
//#include "netcdf.h"
#include "OpenMM.h"
#include "readparm.h"
#include "string_manip.h"

using namespace Amber;
using namespace std;
using namespace OpenMM;

void AmberCoordinateFrame::readRst7(string const& filename) {

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

void AmberCoordinateFrame::readRst7(const char* filename) {
    readRst7(string(filename));
}

void AmberCoordinateFrame::writeRst7(std::string const& filename, bool netcdf) {
    writeRst7(filename.c_str(), netcdf);
}

void AmberCoordinateFrame::writeRst7(const char* filename, bool netcdf) {

    if (netcdf) {
        throw AmberCrdError("NetCDF file writing not yet supported.");
    }

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
