//
// Created by Mark Jones on 9/10/22.
//

#include "VFour.h"

#include <compare>
#include <cmath>
#include <stdexcept>
#include <shared/math/VThree.h>


namespace math {

    double VFour::getR() const {
        return sqrt(vector.getLenSq());
    }

    double VFour::getTheta() const {
        return acos(vector.getCosTheta());
    }

    double VFour::getPhi() const {
        return vector.getPhi();
    }

    double VFour::getLen() const {
        return sqrt(getLenSq());
    }

    double VFour::getLenSq() const {
        return (t * t - vector.getLenSq());
    }

    double VFour::getMass() const {
        double rsq = getLenSq();
        double r = sqrt(abs(rsq));
        return ((rsq < 0) ? -r : r);
    }

    double VFour::dot(const VFour &other) const {
        return t * other.t - vector.dot(other.vector);
    }

    VThree VFour::crossMultiply(const VFour &other) const {
        return vector.crossMultiply(other.vector);
    }

    VFour &VFour::setPolar(double r, double theta, double phi) {
        this->vector.setPolar(r, theta, phi);
        return *this;
    }

    double &VFour::at(int index) {
        switch (index) {
            case 0:
                return t;
            case 1:
            case 2:
            case 3:
                return vector.at(index - 1);
            default:
                throw std::out_of_range("Four vector only has 4 elements");
        }
    }

    bool VFour::operator==(const VFour &other) const {
        return this->t == other.t && (this->vector == other.vector);
    }

    std::partial_ordering VFour::operator<=>(const VFour &other) const {
        return getLenSq() <=> other.getLenSq();
    }

    VFour &VFour::operator+=(const VFour &other) {
        this->t += other.t;
        this->vector += other.vector;
        return *this;
    }

    VFour &VFour::operator-=(const VFour &other) {
        this->t -= other.t;
        this->vector -= other.vector;
        return *this;
    }

    VFour &VFour::operator*=(double a) {
        this->t *= a;
        this->vector *= a;
        return *this;
    }

    VFour VFour::operator+(const VFour &other) const {
        return {this->t + other.t, this->vector + other.vector};
    }

    VFour VFour::operator-(const VFour &other) const {
        return {this->t - other.t, this->vector - other.vector};
    }

    VFour VFour::operator-() const {
        return {-this->t, -this->vector};
    }

    VFour operator*(double a, const VFour &other) {
        return {a * other.t, a * other.vector};
    }

} // math