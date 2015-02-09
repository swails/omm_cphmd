// This tests the UnitCell class here
#include <cassert>
#include <cmath>
#include <iostream>

#include "Amber.h"
#include "OpenMM.h"

using namespace std;
using namespace Amber;

void diff_vecs(OpenMM::Vec3 const &a, OpenMM::Vec3 const&b, int places=6) {
    OpenMM::Vec3 diff = a - b;
    double delta = pow(10, -(double)places);
    for (int i = 0; i < 3; i++)
        assert(diff[i] < delta);
}

void check_default_constructor(void) {

    UnitCell cell;
    OpenMM::Vec3 nullvec;

    assert(cell.getVectorA() == nullvec);
    assert(cell.getVectorB() == nullvec);
    assert(cell.getVectorC() == nullvec);

    assert(cell.getLengthA() == 0);
    assert(cell.getLengthB() == 0);
    assert(cell.getLengthC() == 0);

    try {
        cell.getAlpha();
        assert(0);
    } catch (UnitCellError &e) {
        assert(1);
    }
    try {
        cell.getBeta();
        assert(0);
    } catch (UnitCellError &e) {
        assert(1);
    }
    try {
        cell.getGamma();
        assert(0);
    } catch (UnitCellError &e) {
        assert(1);
    }
}

void check_vec_constructor(void) {

    UnitCell cell(OpenMM::Vec3(1, 0, 0), OpenMM::Vec3(0, 1, 0),
                  OpenMM::Vec3(0, 0, 1));

    assert(cell.getVectorA() == OpenMM::Vec3(1, 0, 0));
    assert(cell.getVectorB() == OpenMM::Vec3(0, 1, 0));
    assert(cell.getVectorC() == OpenMM::Vec3(0, 0, 1));

    assert(cell.getLengthA() == 1);
    assert(cell.getLengthB() == 1);
    assert(cell.getLengthC() == 1);

    assert(abs(cell.getAlpha() - 90) < 1e-5);
    assert(abs(cell.getBeta() - 90) < 1e-5);
    assert(abs(cell.getGamma() - 90) < 1e-5);
}

void check_dim_constructor(void) {

    UnitCell cell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0);

    assert(cell.getVectorA() == OpenMM::Vec3(1, 0, 0));
    assert(cell.getVectorB() == OpenMM::Vec3(0, 1, 0));
    assert(cell.getVectorC() == OpenMM::Vec3(0, 0, 1));

    assert(cell.getLengthA() == 1);
    assert(cell.getLengthB() == 1);
    assert(cell.getLengthC() == 1);

    assert(abs(cell.getAlpha() - 90) < 1e-5);
    assert(abs(cell.getBeta() - 90) < 1e-5);
    assert(abs(cell.getGamma() - 90) < 1e-5);
}

void check_vec_set(void) {

    UnitCell cell;

    cell.setVectorA(OpenMM::Vec3(1, 0, 0));
    cell.setVectorB(OpenMM::Vec3(0, 1, 0));
    cell.setVectorC(OpenMM::Vec3(0, 0, 1));

    assert(cell.getVectorA() == OpenMM::Vec3(1, 0, 0));
    assert(cell.getVectorB() == OpenMM::Vec3(0, 1, 0));
    assert(cell.getVectorC() == OpenMM::Vec3(0, 0, 1));

    assert(cell.getLengthA() == 1);
    assert(cell.getLengthB() == 1);
    assert(cell.getLengthC() == 1);

    assert(abs(cell.getAlpha() - 90) < 1e-5);
    assert(abs(cell.getBeta() - 90) < 1e-5);
    assert(abs(cell.getGamma() - 90) < 1e-5);
}

void check_dim_set(void) {

    UnitCell cell;

    cell.setUnitCell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0);

    assert(cell.getVectorA() == OpenMM::Vec3(1, 0, 0));
    assert(cell.getVectorB() == OpenMM::Vec3(0, 1, 0));
    assert(cell.getVectorC() == OpenMM::Vec3(0, 0, 1));

    assert(cell.getLengthA() == 1);
    assert(cell.getLengthB() == 1);
    assert(cell.getLengthC() == 1);

    assert(abs(cell.getAlpha() - 90) < 1e-5);
    assert(abs(cell.getBeta() - 90) < 1e-5);
    assert(abs(cell.getGamma() - 90) < 1e-5);
}

void check_truncoct(void) {

    UnitCell cell(44.8903851, 44.8903851, 44.8903851,
                  109.471219, 109.471219, 109.471219);

    diff_vecs(cell.getVectorA(), OpenMM::Vec3(44.8903851, 0, 0));
    diff_vecs(cell.getVectorB(),
              OpenMM::Vec3(-14.963460492639706, 42.323061379247044, 0.0));
    diff_vecs(cell.getVectorC(),
              OpenMM::Vec3(-14.963460492639706, -21.161528128425644, 36.65284779906417));

    assert(abs(cell.getLengthA() - 44.8903851) < 1e-6);
    assert(abs(cell.getLengthB() - 44.8903851) < 1e-6);
    assert(abs(cell.getLengthC() - 44.8903851) < 1e-6);

    assert(abs(cell.getAlpha() - 109.471219) < 1e-5);
    assert(abs(cell.getBeta() - 109.471219) < 1e-5);
    assert(abs(cell.getGamma() - 109.471219) < 1e-5);

    UnitCell cell2;
    cell2.setVectorA(OpenMM::Vec3(44.8903851, 0, 0));
    cell2.setVectorB(OpenMM::Vec3(-14.963460492639706, 42.323061379247044, 0.0));
    cell2.setVectorC(OpenMM::Vec3(-14.963460492639706, -21.161528128425644,
                                  36.65284779906417));

    assert(abs(cell2.getLengthA() - 44.8903851) < 1e-6);
    assert(abs(cell2.getLengthB() - 44.8903851) < 1e-6);
    assert(abs(cell2.getLengthC() - 44.8903851) < 1e-6);
    assert(abs(cell2.getAlpha() - 109.471219) < 1e-5);
    assert(abs(cell2.getBeta() - 109.471219) < 1e-5);
    assert(abs(cell2.getGamma() - 109.471219) < 1e-5);
}

int main() {
    cout << "Checking the default UnitCell constructor... ";
    check_default_constructor();
    cout << " OK." << endl;

    cout << "Checking the Vec3 UnitCell constructor... ";
    check_vec_constructor();
    cout << " OK." << endl;

    cout << "Checking the dimension UnitCell constructor... ";
    check_dim_constructor();
    cout << " OK." << endl;

    cout << "Checking vector setting... ";
    check_vec_set();
    cout << " OK." << endl;

    cout << "Check dimension setting... ";
    check_dim_set();
    cout << " OK." << endl;

    cout << "Check truncated octahedron... ";
    check_truncoct();
    cout << " OK." << endl;
}
