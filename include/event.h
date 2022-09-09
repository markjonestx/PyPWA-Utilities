#pragma once

#include <iostream>
#include <string>
#include <list>
#include <pputil.h>
#include <Vec.h>
#include <lorentz.h>
#include <particle.h>

extern particleDataTable PDGtable;

class particle;

class event {
protected:
    std::list<particle> _initial;
    std::list<particle> _final;
    particle *_beam;
    particle *_target;
    int _ioversion;
public:


    event();

    event(const event &);

    ~event();

    event &operator=(const event &e);


    event &addfinal(const particle &);

    event &addinitial(const particle &p);

    event &erase();

    event &beam(const particle &);

    event &target(const particle &);

    particle beam() const;

    particle target() const;

    fourVec
    getPartPFinal(std::string name, int charge, int index, int debug = 0) const;


    std::list<particle> f_mesons() const;

    std::list<particle> f_baryons() const;

    std::list<particle> i_mesons() const;

    std::list<particle> i_baryons() const;

    void print() const;

    friend std::istream &operator>>(std::istream &is, event &e);

    friend std::ostream &operator<<(std::ostream &os, event &e);

    std::istream &read1(std::istream &is);

    std::ostream &write1(std::ostream &os);

    std::istream &read2(std::istream &is);

    std::ostream &write2(std::ostream &os);

    friend event operator*(const lorentzTransform &, const event &);

};
