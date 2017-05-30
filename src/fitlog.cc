#include <parseutils.h>
#include <fitlog.h>

double call_fitlog_nevents(fitlog& fl) {
  return fl.nevents();
}

double call_fitlog_phase(fitlog& fl) {
  return fl.phase();
}

fitlog::fitlog() {
}

fitlog::fitlog(const fitlog& fl) {
	_minus_log_like = fl._minus_log_like;
	_nevents = fl._nevents;
	_lowEdge = fl._lowEdge;
	_highEdge = fl._highEdge;
	_pars = fl._pars;
	_integrals = fl._integrals;
	_emat = fl._emat;
}

fitlog::fitlog(int nev, double le, double he, double mll,
	hash_map<const char*, parameter> pars,
	hash_map<const char*, integral, hash<const char*>, eqstr> ints,
	matrix<double> emat) {
		_minus_log_like = mll;
		_nevents = nev;
		_lowEdge = le;
		_highEdge = he;
		_pars = pars;
		_integrals = ints;
		_emat = emat;

}

fitlog::~fitlog() {
}

fitlog& fitlog::operator=(const fitlog& fl) {
	_minus_log_like = fl._minus_log_like;
	_nevents = fl._nevents;
	_lowEdge = fl._lowEdge;
	_highEdge = fl._highEdge;
	_pars = fl._pars;
	_integrals = fl._integrals;
	_emat = fl._emat;
	
	return *this;
}

list<string> fitlog::waveNames() const {
	list<string> p_names;
	hash_map<const char*, parameter>::const_iterator p = _pars.begin();

	while(p != _pars.end()) {
		p_names.push_back(p->first);
		p++;
	}

	return p_names;
}

hash_map<const char*, parameter> fitlog::parameters() const {
	hash_map<const char*, parameter> ret;
	ret = _pars;
	return ret;
}

hash_map<const char*, integral, hash<const char*>, eqstr> fitlog::integrals() const {
	hash_map<const char*, integral, hash<const char*>, eqstr> ret;
	ret = _integrals;
	return ret;
}

fitlog& fitlog::progFile(const char* pname) {
	_progFile = fopen(pname,"r");
	yyin = _progFile;
	return *this;
}

fitlog& fitlog::progParse() {

	install("events",CTYPE,0.0);

	// --- set the values of the integrals to the ones from logfile
	//     this must be done before parsing since intgeral values are
	//     looked up at parse time.

	// changed how integrals are evaluated in fitparse.yy
	// hash_map<const char*, integral, hash<const char*>, eqstr>::iterator ni = _integrals.begin();
	// while(ni != _integrals.end()) {
//  		set_integral(ni->first,ni->second);
 // 		ni++;
 // 	}

	init();
	init_norm_code();
	yyparse();
	return *this;
}

fitlog& fitlog::setNevents() {
	// --- set variable "nevents" to the number of observed events in
	//     this bin.  this is available for use in the nevents program
	Symbol *sp;
	sp = lookup("nevents");
	//cerr << "looked up nevents: found " << sp->name << " " << sp->val <<endl;
	sp->val = this->_nevents;
	return *this;
}

fitlog& fitlog::setParams() {
	// --- set the values of the parameters to the ones from logfile
	hash_map<const char*, parameter>::iterator p = _pars.begin();
	while(p != _pars.end()) {
		set_param(p->first,p->second);
		_parSet[p->first] = 1;
		p++;
	}
	return *this;
}

fitlog& fitlog::setIntegrals() {
        hash_map<const char*, integral, hash<const char*>, eqstr>::iterator ni =
 _integrals.begin();
        while(ni != _integrals.end()) {
                set_integral(ni->first,ni->second);
                ni++;
        }
	return *this;
}
fitlog& fitlog::setParams(list<string> wvlist) {

	hash_map<const char*, parameter>::const_iterator p = _pars.begin();
	list<string>::const_iterator wv = wvlist.begin();
	parameter par;
	char s[200];

	par.val = complex<double>(0.0,0.0);
	par.real_index = -1;
	par.imag_index = -1;

	// loops over _pars and set all params in symbol table to zero
	while(p != _pars.end()) {
		//cerr << "setting " << p->first << " to zero" << endl;
		set_param(p->first,par);
		_parSet[p->first] = 0;
		p++;
	}
	// now set the requested pars to their value
	while(wv != wvlist.end()) {
		strcpy(s,(*wv).c_str());
		p = _pars.begin();
		while ( p != _pars.end() ) {
			if(string(p->first)==*wv) {
				set_param(*wv,p->second);
				_parSet[p->first] = 1;
			}
			p++;
		}
		wv++;
	}
	return *this;
}

fitlog& fitlog::setParams(list<string> wvlist,int listno) {

	hash_map<const char*, parameter>::const_iterator p = _pars.begin();
	list<string>::const_iterator wv = wvlist.begin();
	parameter par;
	char s[200];

	par.val = complex<double>(0.0,0.0);
	par.real_index = -1;
	par.imag_index = -1;

	// loops over _pars and set all params in symbol table to zero
	while(p != _pars.end()) {
		//cerr << "setting " << p->first << " to zero" << endl;
		set_param(p->first,par);
		if (listno==1)
			_parSet1[p->first] = 0;
		else if (listno==2)
			_parSet2[p->first] = 0;
		p++;
	}
	// now set the requested pars to their value
	while(wv != wvlist.end()) {
		strcpy(s,(*wv).c_str());
		p = _pars.begin();
		while ( p != _pars.end() ) {
			if(string(p->first)==*wv) {
				set_param(*wv,p->second);
				if (listno==1)
					_parSet1[p->first] = 1;
				else if (listno==2)
					_parSet2[p->first] = 1;
			}
			p++;
		}
		wv++;
	}
	return *this;
}

double fitlog::phase() {
	complex<double> z1, z2;
	// loop over pars and add pars in set1 and set2
	hash_map<const char*, parameter>::const_iterator p = _pars.begin();
	while (p != _pars.end()) {
		if(_parSet1[p->first]) z1+=p->second.val;
		if(_parSet2[p->first]) z2+=p->second.val;
		p++;
	}
	return arg(z1*conj(z2));
}

double fitlog::nevents() const{
	double nev = 0;

	//cout << "before nevents() calc." << endl;
	//print_stable();
	execute(norm_prog);
	nev = real(lookup("events")->val);
	return nev;
}

double fitlog::delta(double (*func)(fitlog&)) {

	double deltaf = 0.0;
	matrix<double> dfdp;

	dfdp = ddp(func);
	deltaf = pow( (dfdp*this->_emat*dfdp.transpose()).el(0,0), 0.5 );

	return deltaf;

}

matrix<double> fitlog::ddp(double (*func)(fitlog&)) {

	double fx, fxdx;
	matrix<double> deriv(1,_emat.nrows());
	parameter par, ptmp;
	double eps = 0.0000001;

	fx = func(*this);
	
	hash_map<const char*,int>::const_iterator wp = _parSet.begin();
	while( wp != _parSet.end() ) {
		// par = _pars.find(wp->first)->second;
		par = _pars[wp->first];
		if (wp->second == 1) {
			// wave is on -- set to derivative
			ptmp = par;
			// ``bump'' real part of param
			ptmp.val = complex<double>( par.val.real() + eps, par.val.imag() );
			set_param(wp->first, ptmp);
			_pars[wp->first] = ptmp;
			fxdx = func(*this);
			deriv.el(0, par.real_index-1) = (fxdx - fx)/eps;
			// ``bump'' imag part of param
			if ( par.imag_index != -1 ) {
				ptmp.val = complex<double>( par.val.real(), par.val.imag() + eps );
				set_param(wp->first, ptmp);
				_pars[wp->first] = ptmp;
				fxdx = func(*this);
				deriv.el(0, par.imag_index-1) = (fxdx - fx)/eps;
			}
			// restore parameters original value
			_pars[wp->first] = par;
			set_param(wp->first, par);
		}
		else {
			// wave is off -- set to 0
			deriv.el(0, par.real_index-1) = 0.0;
			deriv.el(0, par.imag_index-1) = 0.0;
		}
		wp++;
	}
	
	return deriv;
	
}

void fitlog::print(ofstream& os) const {
	os << "@begin bin" << endl;
	os << "@begin bin_stats" << endl;
	os << "@val integer n " << _nevents << endl;
	os << "@val float low_edge " << _lowEdge << endl;
	os << "@val float high_edge " << _highEdge << endl;
	os << "@val float minus_log_like " << _minus_log_like << endl;
	os << "@end bin_stats" << endl;

	hash_map<const char*, parameter>::const_iterator p = _pars.begin();
	os << "@begin prod_amps" << endl;
	    while(p != _pars.end()) {
		os << "@val amp " << p->first << " " << p->second.val << " " << p->second.real_index << " " << p->second.imag_index << endl;
		p++;
	    }
	os << "@end prod_amps" << endl;

	hash_map<const char*, integral, hash<const char*>, eqstr>::const_iterator i = _integrals.begin();
	os << "@begin integrals" << endl;
	    while(i != _integrals.end()) {
		os << "@val integral " << i->first << endl;
		i->second.print(os);
		i++;
	    }
	os << "@end integrals" << endl;

	os << "@begin error_matrix" << endl;
		os << "@val matrix<double> emat" << endl;
		os << _emat;
	os << "@end error_matrix" << endl;
	os << "@end bin" << endl;
}

fitlog & fitlog::scan(ifstream & ifs)
{
    string word, kind;
    string type, name;
    int ival;
    double fval;
    complex < double >cpar;
    int r_index;
    int i_index;
    integral integral_val;
    matrix < double >mat_val;

    while (!(ifs >> word).eof() && word != "@begin") ifs >> word;
    if(ifs.eof()) return *this;
    ifs >> kind;
    if (kind != "bin") {
	cerr << "cannot find start of bin: missing @begin bin" << endl;
	exit(1);
    }
    while (!(ifs >> word).eof()) {
	if (word == "@end") {
	    ifs >> kind;
	    if (kind == "bin") return *this;
	}
	if (word == "@begin")
	    ifs >> kind;
	if (kind == "bin_stats") {
	    do {
		ifs >> word;
		if (word == "@val") {
		    ifs >> type;
		    ifs >> name;
		    if (type == "integer" && name == "n") {
			ifs >> ival;
			_nevents = ival;
		    } else 
		    if (type == "float" && name == "low_edge") {
			ifs >> fval;
			_lowEdge = fval;
		    } else 
		    if (type == "float" && name == "high_edge") {
			ifs >> fval;
			_highEdge = fval;
		    } else 
		    if (type == "float" && name == "minus_log_like") {
			ifs >> fval;
			_minus_log_like = fval;
		    } else {
			cerr << "unknown @val type name: " <<
			    type << " " << name << endl;
			abort();
		    }
		}
	    } while (word != "@end");
	    ifs >> word;
	    if (word != "bin_stats") {
		cerr << "@end " << word <<
		    " found, @end bin_stats expected" << endl;
		abort();
	    }
	} else if (kind == "prod_amps") {
	    _pars.erase(_pars.begin(),_pars.end());
	    do {
		ifs >> word;
		if (word == "@val") {
		    ifs >> type;
		    ifs >> name;
		    ifs >> cpar;
		    ifs >> r_index;
		    ifs >> i_index;

		    parameter p;
		    p.val = cpar;
		    p.real_index = r_index;
		    p.imag_index = i_index;
		    char *cname =
			(char *) malloc((name.size() + 1) * sizeof(char));
		    strcpy(cname, name.c_str());
		    _pars[cname] = p;
		}
	    } while (word != "@end");
	    ifs >> word;
	    if (word != "prod_amps") {
		cerr << "@end " << word <<
		    " found, @end prod_amps expected" << endl;
		abort();
	    }
	} else if (kind == "integrals") {
	    _integrals.erase(_integrals.begin(),_integrals.end());
	    do {
		ifs >> word;
		if (word == "@val") {
		    ifs >> type;
		    ifs >> name;
		    integral_val.scan(ifs);
		    char *cname =
			(char *) malloc((name.size() + 1) * sizeof(char));
		    strcpy(cname, name.c_str());
		    _integrals[cname] = integral_val;
		}
	    } while (word != "@end");
	    ifs >> word;
	    if (word != "integrals") {
		cerr << "@end " << word <<
		    " found, @end integrals expected" << endl;
		abort();
	    }
	} else if (kind == "error_matrix") {
	    do {
		ifs >> word;
		if (word == "@val") {
		    ifs >> type;
		    ifs >> name;
		    ifs >> mat_val;
		    _emat = mat_val;
		}
	    } while (word != "@end");
	    ifs >> word;
	    if (word != "error_matrix") {
		cerr << "@end " << word <<
		    " found, @end error_matrix expected" << endl;
		abort();
	    }
	} else {
	    //ERROR
	    cerr << "unrecognized kind " << kind << endl;
	    abort();
	}
    }
    return *this;
}

ifstream& operator>>(ifstream& ifs, fitlog& fl) {
	fl.scan(ifs);
	return ifs;
}

ofstream& operator<<(ofstream& ofs, fitlog& fl) {
	fl.print(ofs);
	return ofs;
}

#ifdef UNDEFINED
char* progname;
int lineno = 1;

void warning(char* s, char* t) {
    fprintf(stderr, "%s: %s", progname, s);
    if (t)
        fprintf(stderr, " %s", t);
    fprintf(stderr, " near line %d\n", lineno);
}

void yyerror(char* s) {
    warning(s, (char*) 0);
    exit(1);
}

void execerror(char* s, char* t) {
    warning(s, t);
    exit(1);
}

void fpecatch(int err) {
    execerror("floating point exception", (char*) 0);
}
#endif

