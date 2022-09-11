//
// Created by Mark Jones on 9/9/22.
//
#pragma once

#include <compare>


namespace math {

    class VThree {
    private:

        double x, y, z;

    public:

        VThree(): x(0.0), y(0.0), z(0.0) {};
        VThree(double _x, double _y, double _z):
            x(_x), y(_y), z(_z) {}

        ~VThree() = default;

        double getX() const {return x;}
        double getY() const {return y;}
        double getZ() const {return z;}
        VThree &setX(double x) { this->x = x; return *this; }
        VThree &setY(double y) { this->y = y; return *this; }
        VThree &setZ(double z) { this->z = z; return *this; }

        double getR() const;
        double getTheta() const;
        double getCosTheta() const;
        double getPhi() const;
        double getLenSq() const;
        double getLen() const;

        double dot(const VThree&) const;
        VThree crossMultiply(const VThree&) const;

        VThree &setPolar(double r, double theta, double phi);

        double &at(int i);

        bool operator==(const VThree &) const;
        std::partial_ordering operator<=>(const VThree &) const;

        VThree& operator+=(const VThree&);
        VThree& operator-=(const VThree&);
        VThree& operator*=(double);

        VThree operator+(const VThree&) const;
        VThree operator-(const VThree&) const;
        VThree operator-() const;
        friend VThree operator*(double, const VThree&);

    };
}
