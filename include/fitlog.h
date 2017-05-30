#ifndef FITLOG_H
#define FITLOG_H

#include <iostream>
#include <fstream>
#include <ext/hash_map>
#include <cstdlib>
#include <unistd.h>
#include <fit.h>
#include <fitparse.h>

namespace std {
  using namespace __gnu_cxx;
}

using namespace std;

int yyparse(void);
extern FILE* yyin;

class fitlog {
	private:
		int _nevents;
		double _minus_log_like;
		double _lowEdge;
		double _highEdge;
		hash_map<const char*, int> _parSet;
		hash_map<const char*, int> _parSet1;
		hash_map<const char*, int> _parSet2;
		hash_map<const char*, parameter> _pars;
		hash_map<const char*, integral, hash<const char*>, eqstr> _integrals;
		matrix<double> _emat;
		FILE* _progFile;

	public:
		fitlog();
		fitlog(const fitlog&);
		fitlog(int nev, double le, double he, double mll,
			hash_map<const char*, parameter> pars,
			hash_map<const char*, integral, hash<const char*>, eqstr> ints,
			matrix<double> emat);
		~fitlog();

		fitlog& operator=(const fitlog&);

		list<string> waveNames() const;
		double binCenter() const { return _lowEdge+binWidth()/2; }
		double binWidth() const {return _highEdge-_lowEdge;}
		hash_map<const char*, parameter> parameters() const;
		hash_map<const char*, integral, hash<const char*>, eqstr> integrals() const;
		fitlog& progFile(const char*);
		fitlog& progParse();
		fitlog& setNevents();
		fitlog& setParams();
		fitlog& setParams(list<string>);
		fitlog& setParams(list<string>,int listno);
		fitlog& setIntegrals();
		double nevents() const;
		double phase();
		double delta(double (*)(fitlog&));
		matrix<double> ddp(double (*)(fitlog&));

		void print(ofstream&) const;
		fitlog& scan(ifstream&);

		friend ifstream& operator>>(ifstream& ifs, fitlog& fl);
		friend ofstream& operator<<(ofstream& ofs, fitlog& fl);

		//Added by PL COltharp
		void get_emat(matrix<double>& emat ){ emat=_emat;}
};


double call_fitlog_nevents(fitlog& fl);
double call_fitlog_phase(fitlog& fl);


#endif

