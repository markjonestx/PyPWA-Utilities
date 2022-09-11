//
// Created by Mark Jones on 9/10/22.
//

#pragma once
#include <math/VThree.h>

#include <compare>

namespace math {

    class VFour {
    private:

        double t;
        VThree vector;

    public:
        VFour(): t(0.0), vector() {};
        VFour(double _t, VThree v): t(_t), vector(v) {};
        VFour(const VFour &vector) = default;
        ~VFour() = default;

        double getT() const { return t; }
        double getX() const { return vector.getX(); }
        double getY() const { return vector.getY(); }
        double getZ() const { return vector.getZ(); }
        VThree getVector() const { return vector; }
        VFour &setT(double t) { this->t = t; return *this; }
        VFour &setX(double x) { vector.setX(x); return *this; }
        VFour &setY(double y) { vector.setY(y); return *this; }
        VFour &setZ(double z) { vector.setZ(z); return *this; }
        VFour &setVector(VThree vector) { this->vector = vector; return *this; }

        double getR() const;
        double getTheta() const;
        double getPhi() const;
        double getLenSq() const;
        double getLen() const;
        double getMass() const;

        double dot(const VFour&) const;
        VThree crossMultiply(const VFour&) const;

        VFour &setPolar(double r, double theta, double phi);

        double &at(int i);

        bool operator==(const VFour &) const;
        std::partial_ordering operator<=>(const VFour &) const;

        VFour &operator+=(const VFour &);
        VFour &operator-=(const VFour &);
        VFour &operator*=(double);

        VFour operator+(const VFour &) const;
        VFour operator-(const VFour &) const;
        VFour operator-() const;

        friend VFour operator*(double, const VFour &);

    };

} // math
