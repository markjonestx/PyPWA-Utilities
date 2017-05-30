#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <complex>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <cstdlib>
#include <matrix.h>

using namespace std;

class integral {
    private:
        matrix<complex<double> > _sum;
        map<string, int> _index;
        int _nwaves;
        int _nevents;
        int _maxEvents;
    public:
        integral();
        integral(char**);
        integral(const integral&);
        ~integral();
        integral& operator=(const integral&);
        integral& integrate();
        integral& files(char**);
        integral& files(list<string>);
        integral& renormalize(int n);
        integral& max(int m);
        integral& events(int n);
        complex<double>& el(string, string);
        int nevents() const;
        list<string> files() const;
        char** files_c_str() const;
        complex<double> val(string, string);
        integral get(char** flist);
        integral get(list<string> flist);
        int index(string s);
        int index(char* s);
        matrix<complex<double> > mat();
        const integral& print(ostream& os = cout) const;
        const integral& print_events(ostream& os = cout) const;
        integral& scan(istream& is = cin);

};

#endif

