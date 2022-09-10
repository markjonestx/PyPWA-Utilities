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

txtEvent &txtEvent::operator=(const txtEvent &e) {
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

txtEvent &txtEvent::addfinal(const particle &p) {
    this->_final.push_back(p);
    return *this;
}

txtEvent &txtEvent::addinitial(const particle &p) {
    this->_initial.push_back(p);
    return *this;
}

txtEvent &txtEvent::erase() {
    if (!_initial.empty()) {
        _initial.erase(_initial.begin(), _initial.end());
    }
    if (!_final.empty()) {
        _final.erase(_final.begin(), _final.end());
    }
    return *this;
}

txtEvent &txtEvent::beam(const particle &p) {
    this->_initial.push_front(p);
    if (this->_beam) {
        *(this->_beam) = p;
    } else {
        this->_beam = new particle(p);
    }
    return *this;
}

txtEvent &txtEvent::target(const particle &p) {
    this->_initial.push_back(p);
    if (this->_target) {
        *(this->_target) = p;
    } else {
        this->_target = new particle(p);
    }
    return *this;
}

particle txtEvent::beam() const {
    return (*(this->_beam));
}

particle txtEvent::target() const {
    return (*(this->_target));
}

void txtEvent::print() const {
    cout << "beam: ";
    this->_beam->get4P().print();
    cout << "target: ";
    this->_target->get4P().print();
    cout << "final particles: ";
    cout << endl;
    list<particle>::const_iterator p = this->_final.begin();
    while (p != this->_final.end()) {
        p->print();
        p++;
    }
}

ostream &operator<<(ostream &os, txtEvent &e) {
    os << e._final.size() + 1 << " " << e.code() << " " << e.weight() << " "
       << e.run() << " " << e.event() << endl;
    fourVec v = e.beam().get4P();
    os << name2id(e._beam->Name(), e._beam->Charge()) << " "
       << e._beam->Charge() << " "
       << v.x() << " " << v.y() << " " << v.z() << " "
       << v.t() << endl;
    list<particle>::iterator part = e._final.begin();
    while (part != e._final.end()) {
        v = part->get4P();
        os << name2id(part->Name(), part->Charge()) << " "
           << part->Charge() << " "
           << v.x() << " " << v.y() << " " << v.z() << " "
           << v.t() << endl;
        part++;
    }
    return os;
}

istream &operator>>(istream &is, txtEvent &e) {
    int nparticles = 0;
    int ptype, q;
    double px, py, pz, t;
    int code;
    double w;
    int run, event;
    string name;
    e.erase();
    particle Target(PDGtable.get("p"), 1);
    Target.set4P(fourVec(Target.Mass(), threeVec(0, 0, 0)));
    e.target(Target);
    is >> nparticles;
    is >> code >> w >> run >> event;
    e.weight(w);
    e.code(code);
    e.run(run);
    e.event(event);
    for (int i = 0; i < nparticles; i++) {
        is >> ptype >> q >> px >> py >> pz >> t;
        name = id2name((Geant_ID) ptype);
        if (i == 0) {
            particle Beam(PDGtable.get(name), q);
            Beam.set4P(fourVec(t, threeVec(px, py, pz)));
            e.beam(Beam);

        } else {
            particle part(PDGtable.get(name), q);
            part.set4P(fourVec(t, threeVec(px, py, pz)));
            e.addfinal(part);
        }
    }
    return is;
}

txtEvent operator*(const lorentzTransform &L, const txtEvent &e) {
    txtEvent r;
    list<particle>::const_iterator p;
    r.beam(L * e.beam());
    r.target(L * e.target());
    p = e._initial.begin();
    while (p != e._initial.end()) {
        r.addinitial(L * (*p));
        p++;
    }
    p = e._final.begin();
    while (p != e._final.end()) {
        r.addfinal(L * (*p));
        p++;
    }
    return r;
}
