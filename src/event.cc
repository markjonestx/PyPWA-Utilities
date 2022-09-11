#include <event.h>
#include <math/VFour.h>
#include <math/VThree.h>
#include <math/io.h>

using std::list;
using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;

event::event() {
    this->_beam = NULL;
    this->_target = NULL;
    this->_ioversion = 1;
}

event::~event() {
    if (this->_beam) delete this->_beam;
    if (this->_target) delete this->_target;
}

event::event(const event &e) {
    this->_final = e._final;
    this->_initial = e._initial;
    this->_beam = new particle(*e._beam);
    this->_target = new particle(*e._target);
    this->_ioversion = e._ioversion;
}

event &event::operator=(const event &e) {
    this->_final = e._final;
    this->_initial = e._initial;
    if (this->_beam) delete this->_beam;
    this->_beam = new particle(*e._beam);
    if (this->_target) delete this->_target;
    this->_target = new particle(*e._target);
    this->_ioversion = e._ioversion;
    return *this;
}


event &event::addfinal(const particle &p) {
    this->_final.push_back(p);
    return *this;
}

event &event::addinitial(const particle &p) {
    this->_initial.push_back(p);
    return *this;
}

event &event::erase() {
    if (!_initial.empty()) {
        _initial.erase(_initial.begin(), _initial.end());
    }
    if (!_final.empty()) {
        _final.erase(_final.begin(), _final.end());
    }
    return *this;
}

event &event::beam(const particle &p) {
    this->_initial.push_front(p);
    if (this->_beam) {
        *(this->_beam) = p;
    } else {
        this->_beam = new particle(p);
    }
    return *this;
}

event &event::target(const particle &p) {
    this->_initial.push_back(p);
    if (this->_target) {
        *(this->_target) = p;
    } else {
        this->_target = new particle(p);
    }
    return *this;
}

particle event::beam() const {
    return (*(this->_beam));
}

particle event::target() const {
    return (*(this->_target));
}

math::VFour
event::getPartPFinal(const string &name, int charge, int index, int debug) const {
    int i = 0;
    if (debug) {
        cout << "Looking for " << name << charge << "[" << index << "] in event"
             << endl;
    }

    for (const auto &p: _final) {
        if (debug) {
            cout << "checking against " << p.Name() << p.Charge() << endl;
        }
        if (p.Name() == name && p.Charge() == charge) {
            i++;
            if (debug) {
                cout << "found one" << endl;
                cout << "checking against index " << i << endl;
            }
            if (i == index) {
                if (debug) {
                    cout << "found the right one, getting 4p" << endl;
                    cout << "4p:" << endl;
                    cout << p.get4P() << endl;
                }
                return p.get4P();
            }
        }
    }
    throw ("PartNotFound");
}


void event::print() const {
    cout << "beam: " << _beam->get4P() << endl;
    cout << "target: " << _target->get4P() << endl;
    cout << "final particles: "<< endl;

    for (const auto &p: _final) {
        p.print();
    }
}

ostream &operator<<(ostream &os, event &e) {
    switch (e._ioversion) {
        case 1:
            return e.write1(os);
        case 2:
            return e.write2(os);
        default:
            throw ("badIOVersion");
    }
}

istream &operator>>(istream &is, event &e) {
    switch (e._ioversion) {
        case 1:
            return e.read1(is);
        case 2:
            return e.read2(is);
        default:
            throw ("badIOVersion");
    }
}

ostream &event::write2(ostream &os) {
    math::VFour v = this->beam().get4P();
    os << "B " << name2id(this->_beam->Name(), this->_beam->Charge()) << " "
       << this->_beam->Charge() << " "
       << v.getT() << " "
       << v.getX() << " " << v.getY() << " " << v.getZ() << " "
       << endl;
    v = this->target().get4P();
    os << "T " << name2id(this->_target->Name(), this->_target->Charge()) << " "
       << this->_target->Charge() << " "
       << v.getT() << " "
       << v.getX() << " " << v.getY() << " " << v.getZ() << " "
       << endl;

    for (const auto &p: _final) {
        v = p.get4P();
        os << "F " << name2id(p.Name(), p.Charge()) << " "
           << p.Charge() << " "
           << v.getT() << " "
           << v.getX() << " " << v.getY() << " " << v.getZ() << " "
           << endl;
    }
    os << "E" << endl;
    return os;
}

istream &event::read2(istream &is) {
    int ptype, q;
    double px, py, pz, t;
    string name;

    char Tag = 0;
    this->erase();
    while (!(is >> Tag).eof()) {
        switch (Tag) {
            case 'I': {
                is >> ptype >> q >> t >> px >> py >> pz;
                name = id2name((Geant_ID) ptype);
                particle part(PDGtable.get(name), q);
                part.set4P(math::VFour(t, math::VThree(px, py, pz)));
                this->addinitial(part);
            }
                break;
            case 'F': {
                is >> ptype >> q >> t >> px >> py >> pz;
                name = id2name((Geant_ID) ptype);
                particle part(PDGtable.get(name), q);
                part.set4P(math::VFour(t, math::VThree(px, py, pz)));
                this->addfinal(part);
            }
                break;
            case 'B': {
                is >> ptype >> q >> t >> px >> py >> pz;
                name = id2name((Geant_ID) ptype);
                particle part(PDGtable.get(name), q);
                part.set4P(math::VFour(t, math::VThree(px, py, pz)));
                this->beam(part);
            }
                break;
            case 'T': {
                is >> ptype >> q >> t >> px >> py >> pz;
                name = id2name((Geant_ID) ptype);
                particle part(PDGtable.get(name), q);
                part.set4P(math::VFour(t, math::VThree(px, py, pz)));
                this->target(part);
            }
                break;
            case 'E':
                return is;
        }
    }
    return is;
}

ostream &event::write1(ostream &os) {
    os << this->_final.size() + 1 << endl;
    math::VFour v = this->beam().get4P();
    os << name2id(this->_beam->Name(), this->_beam->Charge()) << " "
       << this->_beam->Charge() << " "
       << v.getX() << " " << v.getY() << " " << v.getZ() << " "
       << v.getT() << endl;

    for (const auto &part: _final) {
        v = part.get4P();
        os << name2id(part.Name(), part.Charge()) << " "
           << part.Charge() << " "
           << v.getX() << " " << v.getY() << " " << v.getZ() << " "
           << v.getT() << endl;
    }
    return os;
}

istream &event::read1(istream &is) {
    int nparticles = 0;
    int ptype, q;
    double px, py, pz, t;
    string name;

    this->erase();

    particle Target(PDGtable.get("p"), 1);
    Target.set4P(math::VFour(Target.Mass(), math::VThree(0, 0, 0)));
    this->target(Target);

    is >> nparticles;
    for (int i = 0; i < nparticles; i++) {
        is >> ptype >> q >> px >> py >> pz >> t;
        name = id2name((Geant_ID) ptype);
        if (i == 0) {
            particle Beam(PDGtable.get(name), q);
            Beam.set4P(math::VFour(t, math::VThree(px, py, pz)));
            this->beam(Beam);
        } else {
            particle part(PDGtable.get(name), q);
            part.set4P(math::VFour(t, math::VThree(px, py, pz)));
            this->addfinal(part);
        }
    }
    return is;
}


event operator*(const lorentzTransform &L, const event &e) {
    event r;

    r.beam(L * e.beam());
    r.target(L * e.target());

    for (auto &p: e._initial)
        r.addinitial(L * p);

    for (auto &p: e._final)
        r.addfinal(L * p);


    return r;
}

