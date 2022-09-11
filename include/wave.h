#pragma once

#include <complex>
#include <string>

#include <math/VFour.h>
#include <lorentz.h>
#include <particle.h>
#include <event.h>
#include <pputil.h>


#define getsign(x) (((x)>0)?"+":"-")

using namespace std;

class wave : public particle {
private:
    int _m;
    int _epsilon;
    math::VFour _beam;
    math::VFour _target;
    string _channel;
    double _b;
    double _t;

public:
    wave() : particle() {
        this->_b = 0.0;
        this->_t = 0.0;
    }

    wave(const wave &wv);

    ~wave() { ; }

    wave &operator=(const wave &);

    wave setM(int m) {
        this->_m = m;
        return (*this);
    }

    wave setSlope(double b) {
        this->_b = b;
        return (*this);
    }

    wave setT(double t) {
        this->_t = t;
        return (*this);
    }

    math::VFour getBeam() const { return this->_beam; }

    math::VFour getTarget() const { return this->_target; }

    wave channel(char *c) {
        this->_channel = c;
        return (*this);
    }

    string channel() { return (this->_channel); }

    wave fill(const event &, int debug = 0);

    wave &setupFrames(int debug = 0);

    complex<double> decayAmp(int debug = 0);

    wave &operator*=(const lorentzTransform &);

    void print() const;
};
