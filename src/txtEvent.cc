#include <txtEvent.h>
	using std::list;
	using std::string;
	using std::cout;
	using std::endl;
	using std::ostream;
	using std::istream;
	txtEvent::txtEvent() {
		this->_beam = NULL;
		this->_target = NULL;	
		this->_run = this->_event = 0;
	}
	txtEvent::~txtEvent() {
		if (this->_beam) delete this->_beam;
		if (this->_target) delete this->_target;
	}
	txtEvent::txtEvent(const txtEvent& e) {
		this->_final = e._final;
		this->_initial = e._initial;
		this->_beam = new particle(*e._beam);
		this->_target = new particle(*e._target);
		this->_run = e.run();
		this->_event = e.event();
	}
	txtEvent& txtEvent::operator=(const txtEvent& e) {
	  this->_code = e.code();
	  this->_weight = e.weight();
	  this->_run = e.run();
	  this->_event = e.event();
		this->_final = e._final;
		this->_initial = e._initial;
		if (this->_beam) delete this->_beam;
		this->_beam = new particle(*e._beam);
		if (this->_target) delete this->_target;
		this->_target = new particle(*e._target);
		return *this;
	}
	txtEvent& txtEvent::addfinal(const particle& p) {
		this->_final.push_back(p);
		return *this;
	}
	txtEvent& txtEvent::addinitial(const particle& p) {
		this->_initial.push_back(p);
		return *this;
	}
	txtEvent& txtEvent::erase() {
		if( !_initial.empty() ) {
			_initial.erase(_initial.begin(),_initial.end());
		}
		if( !_final.empty() ) {
			_final.erase(_final.begin(),_final.end());
		}
		return *this;
	}
	txtEvent& txtEvent::beam(const particle& p) {
		this->_initial.push_front(p);
		if (this->_beam) {
			*(this->_beam) = p;
		}
		else {
			this->_beam = new particle(p);
		}
		return *this;
	}
	txtEvent& txtEvent::target(const particle& p) {
		this->_initial.push_back(p);
		if (this->_target) {
			*(this->_target) = p;
		}
		else {
			this->_target = new particle(p);
		}
		return *this;
	}
	int txtEvent::OK(double epsilon = 1e-6) const {
		list<particle>::const_iterator p;
		int q_initial = 0;
		int q_final = 0;
		fourVec p_initial, p_final;
		int q_conserved, p_conserved;
		p = this->_initial.begin();
		while ( p != this->_initial.end() ) {
			q_initial += p->Charge();
			p_initial += p->get4P();
			p++;
		}
		p = this->_final.begin();
		while ( p != this->_final.end() ) {
			q_final += p->Charge();
			p_final += p->get4P();
			p++;
		}
		q_conserved = q_initial == q_final;
		p_conserved = (p_initial - p_final).lenSq() < epsilon;
		return(q_conserved && p_conserved);
	}
	particle txtEvent::beam() const{
		return ( *(this->_beam) );
	}
	particle txtEvent::target() const{
		return ( *(this->_target) );
	}
	fourVec txtEvent::getPartPFinal(string name,int charge,int index,int debug) const {
		int i = 0;
		if (debug) {
			cout << "Looking for " << name << charge << "[" << index << "] in txtEvent" << endl;
		}
		list<particle>::const_iterator p = this->_final.begin();
		while (p != this->_final.end() ) {
			if (debug) {
				cout << "checking against " << p->Name() << p->Charge() << endl;
			}
			if ( p->Name() == name && p->Charge() == charge ) {
				i++;
				if (debug) {
					cout << "found one" << endl;
					cout << "checking against index " << i << endl;
				}
				if ( i == index ) {
					if (debug) {
						cout << "found the right one, getting 4p" << endl;
						cout << "4p:" << endl;
						p->get4P().print();
					}
					return p->get4P();
				}
			}
			p++;
		}
		throw("PartNotFound");
		return(fourVec(0,threeVec(0,0,0)));
	}
	fourVec txtEvent::getPartPInitial(string name,int charge,int index) const {
		int i = 1;
		list<particle>::const_iterator p = this->_initial.begin();
		while (p != this->_initial.end() ) {
			if ( p->Name() == name && i++ == index ) {
				return p->get4P();
			}
		}
		throw("PartNotFound");
		return(fourVec(0,threeVec(0,0,0)));
	}
	int txtEvent::f_charge() const{
		int q = 0;
		list<particle>::const_iterator p = _final.begin();
		while( p != _final.end() ) {
			q += p->Charge();
			p++;
		}
		return ( q );
	}
	list<particle> txtEvent::f_mesons() const{
		list<particle> l;
		list<particle>::const_iterator p = _final.begin();
		while( p != _final.end() ) {
			if ( p->J()%2 == 0 ) {
				l.push_back(*p);
			}
			p++;
		}
		return ( l );
	}
	list<particle> txtEvent::f_baryons() const{
		list<particle> l;
		list<particle>::const_iterator p = _final.begin();
		while( p != _final.end() ) {
			if ( p->J()%2 == 1 ) {
				l.push_back(*p);
			}
			p++;
		}
		return ( l );
	}
	list<particle> txtEvent::f_particles() const{
		return ( _final );
	}
	particle txtEvent::f_particle(const string& name, int charge, int index) const{
		int i = 0;
		list<particle>::const_iterator p = this->_final.begin();
		while (p != this->_final.end() ) {
			if ( p->Name() == name && p->Charge() == charge ) {
				i++;
				if ( i == index ) {
					return *p;
				}
			}
			p++;
		}
		throw("PartNotFound");
	}
	int txtEvent::i_charge() const{
		int q = 0;
		list<particle>::const_iterator p = _initial.begin();
		while( p != _initial.end() ) {
			q += p->Charge();
			p++;
		}
		return ( q );
	}
	list<particle> txtEvent::i_mesons() const{
		list<particle> l;
		list<particle>::const_iterator p = _initial.begin();
		while( p != _initial.end() ) {
			if ( p->J()%2 == 0 ) {
				l.push_back(*p);
			}
			p++;
		}
		return ( l );
	}
	list<particle> txtEvent::i_baryons() const{
		list<particle> l;
		list<particle>::const_iterator p = _initial.begin();
		while( p != _initial.end() ) {
			if ( p->J()%2 == 1 ) {
				l.push_back(*p);
			}
			p++;
		}
		return ( l );
	}
	list<particle> txtEvent::i_particles() const{
		return ( _initial );
	}
	particle txtEvent::i_particle(const string& name, int charge, int index) const{
		int i = 0;
		list<particle>::const_iterator p = this->_initial.begin();
		while (p != this->_initial.end() ) {
			if ( p->Name() == name && p->Charge() == charge ) {
				i++;
				if ( i == index ) {
					return *p;
				}
			}
			p++;
		}
		throw("PartNotFound");
	}
	threeVec txtEvent::mesonPlane() const {
		threeVec A,C,N;
		list<particle> i,f;
		list<particle>::const_iterator ip,fp;
		i = this->i_mesons();
		ip = i.begin();
		while (ip != i.end()) {
			A += ip->get3P();
			ip++;
		}
		f = this->f_mesons();
		fp = f.begin();
		while (fp != f.end()) {
			C += fp->get3P();
			fp++;
		}
		if ( (A < threeVec(1e-4,0,0)) || (C < threeVec(1e-4,0,0)) )
			return threeVec(0,0,0);
		N = A / C;
		N *= (1/N.len());
		return N;
	}
	threeVec txtEvent::baryonPlane() const {
		threeVec B,D,N;
		list<particle> i,f;
		list<particle>::const_iterator ip,fp;
		i = this->i_baryons();
		ip = i.begin();
		while (ip != i.end()) {
			B += ip->get3P();
			ip++;
		}
		f = this->f_baryons();
		fp = f.begin();
		while (fp != f.end()) {
			D += fp->get3P();
			fp++;
		}
		if ( (B < threeVec(1e-4,0,0)) || (D < threeVec(1e-4,0,0)) )
			return threeVec(0,0,0);
		N = B / D;
		N *= (1/N.len());
		return N;
	}
	void txtEvent::print() const {
		cout << "beam: ";
		this->_beam->get4P().print();
		cout << "target: ";
		this->_target->get4P().print();
		cout << "final particles: ";
		cout << endl;
		list<particle>::const_iterator p = this->_final.begin();
		while( p != this->_final.end() ) {
			p->print();
			p++;
		}
	}
	ostream& operator<<(ostream& os, txtEvent& e) {
	  os << e._final.size()+1 << " " << e.code() << " " << e.weight() << " " << e.run() << " " << e.event()<< endl;
		fourVec v = e.beam().get4P();
		os << name2id(e._beam->Name(),e._beam->Charge()) << " " 
			<< e._beam->Charge() << " " 
			<< v.x() << " " << v.y() << " " << v.z() << " " 
			<< v.t() << endl;
		list<particle>::iterator part = e._final.begin();
		while (part != e._final.end()) {
			v = part->get4P();
			os << name2id(part->Name(),part->Charge()) << " "
				<< part->Charge() << " "
				<< v.x() << " " << v.y() << " " << v.z() << " "
				<< v.t() << endl;
			part++;
		}
		return os;
	}
	istream& operator>>(istream& is, txtEvent& e) {
		int nparticles = 0;
		int ptype, q;
		double px, py, pz, t;
		int code;
		double w;
		int run,event;
		string name;
		e.erase();
		particle Target(PDGtable.get("p"),1);
		Target.set4P(fourVec(Target.Mass(),threeVec(0,0,0)));
		e.target(Target);
		is >> nparticles;
		is >> code >> w >> run >> event;
		e.weight(w);
		e.code(code);
		e.run(run);
		e.event(event);
		for (int i = 0; i < nparticles; i++ ) {
			is >> ptype >> q >> px >> py >> pz >> t;
			name = id2name( (Geant_ID) ptype );
			if ( i==0 ) {
				particle Beam(PDGtable.get(name),q);
				Beam.set4P(fourVec(t,threeVec(px,py,pz)));
				e.beam(Beam);
				
			}
			else {
				particle part(PDGtable.get(name),q);
				part.set4P(fourVec(t,threeVec(px,py,pz)));
				e.addfinal(part);
			}
		}
		return is;
	}	
	txtEvent operator*(const lorentzTransform& L,const txtEvent& e) {
		txtEvent r;
		list<particle>::const_iterator p;
		r.beam(L*e.beam());
		r.target(L*e.target());
		p = e._initial.begin();
		while (p != e._initial.end()) {
			r.addinitial(L*(*p));
			p++;
		}
		p = e._final.begin();
		while (p != e._final.end()) {
			r.addfinal(L*(*p));
			p++;
		}
		return r;
	}
