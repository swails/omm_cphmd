/// Contains some constants for unit conversions b/w Amber and OpenMM

#ifndef AMBER_CONSTANTS_H
#define AMBER_CONSTANTS_H

#include <cmath>

namespace Amber {
// Some conversion constants
static const double ANGSTROM_PER_NANOMETER = 10.0;
static const double JOULE_PER_CALORIE = 4.184;
static const double DEGREE_PER_RADIAN = 180.0 / M_PI;

static const double NANOMETER_PER_ANGSTROM = 1 / ANGSTROM_PER_NANOMETER;
static const double CALORIE_PER_JOULE = 1 / JOULE_PER_CALORIE;
static const double RADIAN_PER_DEGREE = 1 / DEGREE_PER_RADIAN;

static const double AMBER_TIME_PER_PS = 20.455;
static const double PS_PER_AMBER_TIME = 1 / AMBER_TIME_PER_PS;

}; // namespace Amber

#endif /* AMBER_CONSTANTS_H */
