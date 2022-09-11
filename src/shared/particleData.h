#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <shared/pputil.h>

#define getsign(x) (((x)>0)?"+":"-")


class particleData {
private:
    static int _particle_data_count;
    static int _particle_data_debug;
    std::string _name;
    double _mass;
    double _width;
    int _isospin;
    int _gparity;
    int _spin;
    int _parity;
    int _cparity;

    void _init(std::string, double, double, int, int, int, int, int);

public:

    particleData();

    particleData(const particleData &);

    particleData(std::string n, double m, double w, int i, int g, int j, int p,
                 int c);

    ~particleData();

    particleData &setJ(int j) {
        this->_spin = j;
        return (*this);
    }

    particleData &setP(int p) {
        this->_parity = p;
        return (*this);
    }

    std::string Name() const;

    double Mass() const {
        return (this->_mass);
    }

    double Width() const {
        return (this->_width);
    }

    int I() const {
        return (this->_isospin);
    }

    int G() const {
        return (this->_gparity);
    }

    int J() const {
        return (this->_spin);
    }

    int P() const {
        return (this->_parity);
    }

    int C() const {
        return (this->_cparity);
    }

    void print() const;

    void dump() const;

    particleData &operator=(const particleData &);

    void debug(int d = 1) {
        _particle_data_debug = d;
    }

};


class tableEntry {

private:
    particleData particle;
    tableEntry *nextparticle;

    void _init(const particleData &p, tableEntry *n);

public:
    tableEntry(particleData p, tableEntry *n);

    tableEntry(const tableEntry &te);

    ~tableEntry();

    tableEntry &operator=(const tableEntry &te);

    tableEntry *next() const;

    particleData Particle() const;

    void print() const;

    void dump() const;
};


class particleDataTable {

private:
    tableEntry *head;

public:
    particleDataTable(tableEntry *p = NULL);

    void initialize();

    void initialize(char *PDTfile);

    void insert(particleData p);

    particleData get(std::string _name) const;

    double mass(std::string _name) const;

    void print() const;

    void dump() const;
};
