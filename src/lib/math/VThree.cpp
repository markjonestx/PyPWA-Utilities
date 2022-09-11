//
// Created by Mark Jones on 9/9/22.
//

#include "VThree.h"

#include <cmath>
#include <compare>
#include <stdexcept>


namespace math {

    double VThree::getR() const {
        return pow(this->getLenSq(), 0.5);
    }

    double VThree::getTheta() const {
        return acos(this->getCosTheta());
    }

    double VThree::getCosTheta() const {
        return (this->z / this->getR());
    }

    double VThree::getPhi() const {
        return atan2(this->y, this->x);
    }

    double VThree::getLenSq() const {
        return (
                this->x * this->x + this->y * this->y + this->z * this->z
        );
    }

    double VThree::getLen() const {
        return pow(this->getLenSq(), 0.5);
    }

    double VThree::dot(const VThree &V) const {
        return (this->x * V.x + this->y * V.y + this->z * V.z);
    }

    VThree VThree::crossMultiply(const VThree &V) const {
        return {
                this->y * V.z - this->z * V.y,
                this->z * V.x - this->x * V.z,
                this->x * V.y - this->y * V.x
        };
    }

    VThree &VThree::setPolar(double r, double theta, double phi) {
        this->x = r * sin(theta) * cos(phi);
        this->y = r * sin(theta) * sin(phi);
        this->z = r * cos(theta);
        return *this;
    }

    double &VThree::at(int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                throw std::out_of_range("Three Vector on has 3 elements");
        }
    }

    bool VThree::operator==(const VThree &V) const {
        return (this->x == V.x
                && this->y == V.y
                && this->z == V.z);
    }

    std::partial_ordering VThree::operator<=>(const VThree &V) const {
        return this->getLenSq() <=> V.getLenSq();
    }

    VThree &VThree::operator+=(const VThree &V) {
        this->x += V.x;
        this->y += V.y;
        this->z += V.z;
        return *this;
    }

    VThree &VThree::operator-=(const VThree &) {
        this->x -= x;
        this->y -= y;
        this->z -= z;
        return *this;
    }

    VThree &VThree::operator*=(double a) {
        this->x *= a;
        this->y *= a;
        this->z *= a;
        return *this;
    }

    VThree VThree::operator+(const VThree &other) const {
        return {
            this->x + other.x, this->y + other.y, this->z + other.z
        };
    }

    VThree VThree::operator-(const VThree &other) const {
        return {
            this->x - other.x, this->y - other.y, this->z - other.z
        };
    }

    VThree VThree::operator-() const {
        return { -this->x, -this->y, -this->z };
    }

    VThree operator*(double a, const VThree &V) {
        return {a * V.x, a * V.y, a * V.z};
    }
}