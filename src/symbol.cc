#include <cstdlib>
#include <fstream>
#include <complex>
#include <ext/hash_map>
#include <fit.h>
#include <minuit.h>
#include <integral.h>
#include <fitparse.h>

namespace std {
  using namespace __gnu_cxx;
}

using namespace std;

int debug = 0;

static hash_map<const char*, Symbol*, hash<const char*>, eqstr> symtable;

int n_amps = 0;

void print_stable() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		cout << sym->first << " <-> " << sym->second << "->name=" << sym->second->name << " ->val=" << sym->second->val << endl;
		sym++;
	}
}

Symbol* lookup(char* s) {
	Symbol* ret;
	if (debug) {
		cout << "begin lookup of " << s << endl;
		cout << "sym table before lookup:" << endl;
		print_stable();
	}
	// ret = symtable[s];
	if( symtable.find(s) != symtable.end() )
		ret = symtable.find(s)->second;
	else
		ret = NULL;
	if (debug) {
		if (ret) {
			cout << "symbol found is " << ret->name << endl;
			cout << "symbol is nevents: val is " << ret->val <<endl;
		} else {
			cout << "symbol " << s << " not found, returning " << ret << endl;
		}
		cout << "sym table after lookup:" << endl;
		print_stable();
	}
	if (debug) {
		cout << "end lookup of " << s << " returning " << ret << endl;
		cout << flush;
	}
	return ret;
}

Symbol *install(char *s, int t, complex<double> d) {
	Symbol* sym_p;
	char *emalloc(unsigned n);
	
	
	if (debug) 
		cout << "begin install of " << s << endl;
	sym_p = (Symbol*) emalloc(sizeof(Symbol));
	sym_p->name = emalloc(strlen(s)+1);
	strcpy(sym_p->name, s);
	sym_p->type = t;
	sym_p->val = d;
	
	if (debug) {
		cout << "installing " << sym_p->name << " as " << s << endl;
		cout << "sym table before install:" << endl;
		print_stable();
	}
	symtable[sym_p->name] = sym_p;
	if (debug) {
		cout << "sym table after install:" << endl;
		print_stable();
		cout << "end install of " << s << endl;
	}
	
	return symtable[s];

}

void count_amps() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();

	n_amps = 0;
	while (sym != end) {
		if( (sym->second)->type == ATYPE ) n_amps++;
		sym++;
	}
}

static hash_map<const char*, vector< complex<double> >, hash<const char*>, eqstr > amp;

static int amp_ptr = 0;
static int namps = 0;

void read_amps() {
	int nevent = 0;
	complex<double> a;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == ATYPE) {
			nevent = 0;
			ifstream ampfile(sym->second->name);
			while (ampfile.read((char*)&a,sizeof(complex<double>)) ) {
				amp[sym->second->name].push_back(a);
				nevent++;
			}
			if (!namps) namps = nevent;
		}
		sym++;
	}
	lookup("nevents")->val = nevent;
}

int update_amps() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == ATYPE) {
			sym->second->val = amp[sym->second->name][amp_ptr];
		}
		sym++;
	}
	amp_ptr++;
	return (amp_ptr <= namps);
}

static hash_map<const char*, ifstream*> ampfiles;

void open_ampfiles() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == ATYPE) {
			ampfiles[sym->second->name] = new ifstream(sym->second->name);
		}
		sym++;
	}
}

int update_amps_from_files() {
	int read_ok = 0;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == ATYPE) {
		  read_ok = !(ampfiles[sym->second->name]->read((char*)&(sym->second->val),sizeof(sym->second->val))).eof();
		}
		sym++;
	}
	return (read_ok);
}

void rewind_amps() {
	amp_ptr = 0;
}



/* rmat stuff */

static hash_map<const char*, vector< matrix<double> >, hash<const char*>, eqstr > rmat;

static int rmat_ptr = 0;
static int nrmat = 0;

void read_rmat() {
	int nevent = 0;
	matrix<double> a;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == RMTYPE) {
			nevent = 0;
			ifstream rmatfile(sym->second->name);
			while (!(rmatfile >> a).eof() ) {
			  //		  cout << a;
				rmat[sym->second->name].push_back(a);
				nevent++;
			}
			if (!nrmat) nrmat = nevent;
		}
		sym++;
		
	}

}

int update_rmat() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == RMTYPE) {
			sym->second->rmatval = rmat[sym->second->name][rmat_ptr];
		}
		sym++;
	}
	rmat_ptr++;
	return (rmat_ptr <= nrmat);
}

int update_symbols() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == ATYPE) {
			sym->second->val = amp[sym->second->name][amp_ptr];
		}
		else if (sym->second->type == RMTYPE) {
			sym->second->rmatval = rmat[sym->second->name][rmat_ptr];
		}
		sym++;
	}
	rmat_ptr++;
	amp_ptr++;
	return (amp_ptr <= namps);
}

static hash_map<const char*, ifstream*, hash<const char*>, eqstr> rmatfiles;

void open_rmatfiles() {
	matrix<double> a;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == RMTYPE) {
			rmatfiles[sym->second->name] = new ifstream(sym->second->name);
		}
		sym++;
	}
}

int update_rmat_from_files() {
	int read_ok = 0;
	int found_rmat = 0;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == RMTYPE) {
		  found_rmat = 1;
		  read_ok = !(*rmatfiles[sym->second->name] >> sym->second->rmatval).eof();
		  //		  cout << sym->second->rmatval;
		}
		sym++;
	}
	return (found_rmat?read_ok:1);
}


void rewind_rmat() {
	rmat_ptr = 0;
}




void set_param(string pname, parameter p) {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;

	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == PTYPE) {
			if(string(sym->second->name) == pname) {
				sym->second->par = p;
			}
		} else if (sym->second->type == RPTYPE) {
			if(string(sym->second->name) == pname) {
				if (p.val.imag() == 0.0 && p.imag_index == -1) {
					sym->second->par = p;
				} else {
					cerr << "error: setting real parameters imaginary part: setting val to " << p.val << " with MINUIT indices (" << p.real_index << "," << p.imag_index << ")" << endl;
					throw("PARAM_RANGE");
				}
			}
		}
		sym++;
	}
}

void set_params() {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	int parno = 1;
	string pname = "parameter";
	double start = 1.0;
	double step = 1.0;
	double lbound = 0.0, ubound = 0.0;

	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == PTYPE) {
			pname = string(sym->second->name) + ".re";
			mnparm(parno,pname,start,step,lbound,ubound);
			sym->second->par.real_index = parno++;
			pname = string(sym->second->name) + ".im";
			mnparm(parno,pname,start,step,lbound,ubound);
			sym->second->par.imag_index = parno++;
		} 
		else if (sym->second->type == RPTYPE) {
			pname = string(sym->second->name) + ".re";
			mnparm(parno,pname,start,step,lbound,ubound);
			sym->second->par.real_index = parno++;
			sym->second->par.imag_index = -1;
		}
		sym++;
	}

}

void update_params(double *minval) {
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
	hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
	sym = symtable.begin();
	end = symtable.end();
	while (sym != end) {
		if (sym->second->type == PTYPE) {
		  sym->second->par.val = complex<double>(minval[sym->second->par.real_index-1],minval[sym->second->par.imag_index-1]);	
		}
		else if (sym->second->type == RPTYPE) {
		  sym->second->par.val = complex<double>(minval[sym->second->par.real_index-1],0.0);	
		}
		sym++;
	}
}

hash_map<const char*, parameter>
get_params() {

    hash_map<const char*, parameter> ret;
    
    hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
    hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
    sym = symtable.begin();
    end = symtable.end();
    while (sym != end) {
        if (sym->second->type == PTYPE || sym->second->type == RPTYPE) {
            ret[sym->first] = sym->second->par;
        }
        sym++;
    }
    return ret;
}

void
print_params(ofstream& os) {

    hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator sym;
    hash_map<const char*, Symbol*, hash<const char*>, eqstr>::iterator end;
    sym = symtable.begin();
    end = symtable.end();
    while (sym != end) {
        if (sym->second->type == PTYPE || sym->second->type == RPTYPE) {
            os << sym->first << " = " << sym->second->par.val << endl;
        }
        sym++;
    }
}

char *emalloc(unsigned n) {
	char *p;

	p = (char*) malloc(n);
	if (p == 0)
		execerror("out of memory", (char*) 0);
	return p;
}


static hash_map<const char*, integral, hash<const char*>, eqstr> integrals;

hash_map<const char*, integral, hash<const char*>, eqstr> get_integrals() {
	return (integrals);
}

void set_integral(string s, integral ni) {
	char* ss = (char*) malloc(s.length()*sizeof(char));
	strcpy(ss,s.c_str());

	integrals[ss] = ni;
}

void read_integral(string s) {
	integral ni;
	ifstream ni_file(s.c_str());
	char* ss = (char*) malloc(s.length()*sizeof(char));

	strcpy(ss,s.c_str());

	integrals[ss] = ni.scan(ni_file);
}

void read_integral(string s, string file) {
	integral ni;
	ifstream ni_file(file.c_str());
	char* ss = (char*) malloc(s.length()*sizeof(char));

	strcpy(ss,s.c_str());

	integrals[ss] = ni.scan(ni_file);
}

integral get_integral(string s) {
	integral ni;
	ni = integrals[s.c_str()];
	return ni;
}

complex<double> get_integral_val(string s, string wv1, string wv2) {
	complex<double> ret;
	integral ni;
	ni = integrals[s.c_str()];
	ret = ni.val(wv1,wv2);
	return ret;
}

void print_integral(string s) {
	integrals[s.c_str()].print();
}

void print_integrals(ofstream& os) {
	hash_map<const char*, integral, hash<const char*>, eqstr>::iterator norm_int;
	norm_int = integrals.begin();
	while (norm_int != integrals.end()) {
		os << "norm_int(" << norm_int->first << "):" << endl;
		norm_int->second.print(os);
		norm_int++;
	}
}
