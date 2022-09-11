#pragma once

#include <cmath>
#include <shared/matrix.h>
#include <shared/math/VFour.h>
#include <shared/math/VThree.h>


class rotation : public matrix<double> {

public:
    rotation() : matrix<double>(3, 3) { ; }

    rotation(double, double, double);

    ~rotation() { ; }

    rotation set(double, double, double);

    rotation set(const math::VThree &);

    friend math::VThree operator*=(math::VThree &, const rotation &);

};


class lorentzTransform : public matrix<double> {
    friend class math::VFour;

private:
    double _gamma;

public:
    lorentzTransform() : matrix<double>(4, 4) {
        _gamma = 1;
    }

    lorentzTransform(double, double, double);

    lorentzTransform(const math::VThree &beta);

    lorentzTransform(const math::VFour &);

    lorentzTransform(const rotation &);

    ~lorentzTransform() { ; }

    lorentzTransform set(double, double, double); // rotation
    lorentzTransform set(const math::VThree &beta); // boost
    lorentzTransform set(const math::VFour &); //boost to rest frame
    lorentzTransform set(const rotation &);

    friend math::VFour operator*=(math::VFour &, const lorentzTransform &);

};
