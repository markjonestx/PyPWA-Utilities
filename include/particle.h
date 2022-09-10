#pragma once

#include <list>
#include <string>
#include <Vec.h>
#include <matrix.h>
#include <lorentz.h>
#include <particleData.h>
#include <event.h>
#include <pputil.h>
#include <massDep.h>


class Decay;

class massDep;

class particle : public particleData {

private:
    static int _particle_debug;
    int _lambda;
    int _charge;
    fourVec _p;
    Decay *_decay;
    int _index;
    std::list<int> _helicities;
    int _inRestFrame;
    massDep *_massDep;

public:

    particle();

    particle(const particleData &, int);

    particle(const particle &);

    ~particle();

    particle &operator=(const particle &);

    particle &set4P(const fourVec &);

    particle &Index(int);

    particle &setDecay(const Decay &);

    particle &setMassDep(massDep *md) {
        _massDep = md;
        return *this;
    }

    friend particle operator*(const lorentzTransform &, const particle &);

    particle &operator*=(const lorentzTransform &);

    Decay *get_decay() const;

    int Stable() const;

    int operator==(const particle &);

    int operator!=(const particle &);

    int operator<(const particle &);

    fourVec get4P() const;

    fourVec *get4P(particle *p, int debug = 0) const;

    threeVec get3P() const;

    int Index() const;

    int Charge() const;

    std::list<int> &helicities();

    particle &addHelicity(int lam);

    int is(std::string) const;

    fourVec setupFrames(int debug = 0);

    std::string sprint(std::string space = " ");

    double q() const;

    double q0() const;

    std::complex<double> decayAmp(int, int debug = 0);

    void print() const;

    void printFrames() const;

    void debug(int d = 1) {
        _particle_debug = d;
    }

};


class event;

class Decay {

private:
    // list<particle> _children;
    static int _decay_debug;
    std::list<particle> _childrenInFrames;
    int _l;
    int _s;
    double _mass;

    void _init(const std::list<particle> &, int, int, double);


public:
    std::list<particle> _children;

    Decay();

    Decay(const Decay &);

    ~Decay();

    Decay &addChild(const particle &);

    Decay &setL(int);

    Decay &setS(int);

    int L() const;

    int S() const;

    Decay &calculateS();

    fourVec fill(const event &, int debug = 0);

    fourVec *get4P(particle *part, int debug = 0);

    Decay &setupFrames(lorentzTransform T, int debug = 0);

    std::complex<double> amp(int, int, int debug = 0) const;

    std::complex<double> expt_amp(double b, double t, int debug = 0) const;

    Decay &operator*=(const lorentzTransform &L);

    Decay &operator=(const Decay &);

    void print() const;

    void printFrames() const;

    void debug(int d = 1) {
        _decay_debug = d;
    }

};
