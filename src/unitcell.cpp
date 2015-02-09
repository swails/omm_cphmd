/// unitcell.cpp

#include <cmath>

#include "amber/exceptions.h"
#include "amber/unitcell.h"

using namespace std;
using namespace Amber;

#define RAD_PER_DEG  0.0174532925199432954743716805978692718782
#define DEG_PER_RAD 57.2957795130823228646477218717336654663086

#define TINY 1e-6

UnitCell::UnitCell(double a, double b, double c,
                   double alpha, double beta, double gama) {
    setUnitCell(a, b, c, alpha, beta, gama);
}

UnitCell::UnitCell(const OpenMM::Vec3 &a, const OpenMM::Vec3 &b,
                   const OpenMM::Vec3 &c) : a_(a), b_(b), c_(c) {}

double UnitCell::getLengthA(void) const {
    return sqrt(a_.dot(a_));
}

double UnitCell::getLengthB(void) const {
    return sqrt(b_.dot(b_));
}

double UnitCell::getLengthC(void) const {
    return sqrt(c_.dot(c_));
}

double UnitCell::getAlpha(void) const {
    double lb = getLengthB();
    double lc = getLengthC();
    if (lb == 0 || lc == 0)
        throw UnitCellError("Cell vectors of 0 detected");
    return acos(b_.dot(c_) / (lb*lc)) * DEG_PER_RAD;
}

double UnitCell::getBeta(void) const {
    double la = getLengthA();
    double lc = getLengthC();
    if (la == 0 || lc == 0)
        throw UnitCellError("Cell vectors of 0 detected");
    return acos(a_.dot(c_) / (la*lc)) * DEG_PER_RAD;
}

double UnitCell::getGamma(void) const {
    double la = getLengthA();
    double lb = getLengthB();
    if (la == 0 || lb == 0)
        throw UnitCellError("Cell vectors of 0 detected");
    return acos(a_.dot(b_) / (la*lb)) * DEG_PER_RAD;
}
void UnitCell::setUnitCell(double a, double b, double c,
                           double alpha, double beta, double gama) {
    double al = alpha * RAD_PER_DEG;
    double be = beta * RAD_PER_DEG;
    double ga = gama * RAD_PER_DEG;

    a_[0] = a;
    a_[1] = 0;
    a_[2] = 0;

    b_[0] = b * cos(ga);
    b_[1] = b * sin(ga);
    b_[2] = 0;

    c_[0] = c * cos(be);
    c_[1] = c * (cos(al) - cos(be)*cos(ga)) / sin(ga);
    c_[2] = sqrt(abs(c*c - c_[0]*c_[0] - c_[1]*c_[1]));

    // Now make sure if anything is CLOSE to 0, we set it to 0 exactly
    b_[0] = abs(b_[0]) < TINY ? 0 : b_[0];
    b_[1] = abs(b_[1]) < TINY ? 0 : b_[1];
    c_[0] = abs(c_[0]) < TINY ? 0 : c_[0];
    c_[1] = abs(c_[1]) < TINY ? 0 : c_[1];
    c_[2] = abs(c_[2]) < TINY ? 0 : c_[2];
}

