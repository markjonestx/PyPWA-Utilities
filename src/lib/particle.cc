#include "particle.h"

#include <lib/particleData.h>
#include <lib/math/io.h>


using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::complex;
using std::string;


particle::particle() : particleData() {
    if (_particle_debug) {
        cout << "in particle(" << this << ")::particle()" << endl;
    }
    this->_charge = 0;
    this->_index = 0;
    this->_lambda = 0;
    this->_inRestFrame = 0;
    this->_p = math::VFour(0, math::VThree(0, 0, 0));
    this->_decay = NULL;
    this->_massDep = NULL;
}

particle::particle(const particle &particle) : particleData(particle) {
    if (_particle_debug) {
        cout << "in particle(" << this
             << ")::particle(const particle& particle(" << &particle << ")="
             << particle.Name() << ")" << endl;
    }
    this->_charge = particle._charge;
    this->_p = particle._p;
    this->_index = particle._index;
    this->_lambda = particle._lambda;
    this->_inRestFrame = particle._inRestFrame;
    this->_helicities = particle._helicities;
    if (particle._decay) {
        this->_decay = new Decay(*(particle._decay));
    } else {
        this->_decay = NULL;
    }
    if (particle._massDep) {
        _massDep = particle._massDep->create();
    } else {
        this->_massDep = NULL;
    }
}

particle::particle(const particleData &data, int c) : particleData(data) {
    if (_particle_debug) {
        cout << "in particle(" << this
             << ")::particle(const particleData& data(" << &data << ")="
             << data.Name()
             << ",int c)" << endl;
    }
    this->_charge = c;
    this->_index = 0;
    this->_lambda = 0;
    this->_inRestFrame = 0;
    this->_p = math::VFour(0, math::VThree(0, 0, 0));
    this->_decay = NULL;
    this->_massDep = NULL;
}

particle::~particle() {
    if (_particle_debug) {
        cout << "in particle((" << this << ")=" << this->Name()
             << ")::~particle()" << endl;
    }
    // if (_decay) delete _decay;
    delete _decay;
    delete _massDep;
}

particle &particle::operator=(const particle &particle) {
    if (_particle_debug) {
        cout << "in particle(" << this
             << ")::operator=(const particle& particle(" << &particle << ")="
             << particle.Name() << ")" << endl;
    }
    if (this != &particle) {
        particleData::operator=(particle);
        this->_charge = particle._charge;
        this->_p = particle._p;
        this->_index = particle._index;
        this->_lambda = particle._lambda;
        this->_inRestFrame = particle._inRestFrame;
        this->_helicities = particle._helicities;
        delete this->_decay;
        if (particle._decay) {
            this->_decay = new Decay(*(particle._decay));
        } else {
            this->_decay = NULL;
        }
        if (particle._massDep) {
            massDep *md = particle._massDep->create();
            delete _massDep;
            _massDep = md;
        } else {
            _massDep = NULL;
        }
    }
    return (*this);
}

int particle::operator==(const particle &p) {
    return (this->Mass() == p.Mass());
}

int particle::operator!=(const particle &p) {
    return (this->Mass() != p.Mass());
}

int particle::operator<(const particle &p) {
    return (this->Mass() < p.Mass());
}

particle &particle::set4P(const math::VFour &p4) {
    this->_p = p4;
    return (*this);
}

particle &particle::Index(int i) {
    this->_index = i;
    return (*this);
}

particle &particle::setDecay(const Decay &d) {
    if (this->_decay) {
        delete this->_decay;
    }
    this->_decay = new Decay(d);
    return (*this);
}

particle operator*(const lorentzTransform &L, const particle &p) {
    particle part;
    part = p;
    part._p *= L;
    return part;
}

particle &particle::operator*=(const lorentzTransform &L) {
    this->_p *= L;
    if (this->_decay) {
        *(this->_decay) *= L;
    }
    return (*this);
}

Decay *particle::get_decay() const {
    return (this->_decay);
}


void Decay::_init(const list<particle> &cl, int l, int s, double m) {
    _children = cl;
    _l = l;
    _s = s;
    _mass = m;
}

Decay::Decay() {
    _l = 0;
    _s = 0;
}

Decay::Decay(const Decay &d) {
    if (_decay_debug) {
        cout << "in Decay(" << this << ")::Decay(const Decay& d)" << endl;
    }
    _init(d._children, d._l, d._s, d._mass);
}

Decay::~Decay() {
    if (_decay_debug) {
        cout << "in Decay(" << this << ")::Decay()" << endl;
    }
}

Decay &Decay::addChild(const particle &p) {
    _children.push_back(p);
    return (*this);
}

Decay &Decay::setL(int l) {
    _l = l;
    return (*this);
}

Decay &Decay::setS(int s) {
    _s = s;
    return (*this);
}

int Decay::L() const {
    return (this->_l);
}

int Decay::S() const {
    return (this->_s);
}

Decay &Decay::calculateS() {
    int spin = 0;
    int numNonZero = 0;
    list<particle>::iterator child;

    child = this->_children.begin();
    while (child != this->_children.end()) {
        if (child->J() != 0) {
            spin += child->J();
            numNonZero++;
        }
        child++;
    }
    if (numNonZero > 1) {
        cerr << "final state spin is undetermined in Decay: " << endl;
        this->print();
        exit(1);
    } else {
        this->_s = spin;
        return *this;
    }
}

Decay &Decay::operator=(const Decay &d) {
    if (_decay_debug) {
        cout << "in Decay(" << this << ")::operator=(const Decay& d)" << endl;
    }
    _init(d._children, d._l, d._s, d._mass);
    return (*this);
}

void Decay::print() const {
    ptab();
    cout << "children : {" << endl;

    for (const auto &c : _children)
        c.print();

    ptab();
    cout << "L: " << _l << " S: " << _s << endl;
    ptab();
    cout << "}" << endl;
}

void Decay::printFrames() const {
    ptab();
    cout << "i have " << this->_childrenInFrames.size() << " children." << endl;
    ptab();
    cout << "children in Decay frame : {" << endl;

    for (const auto &c: _children)
        c.printFrames();

    ptab();
    cout << "L: " << _l << " S: " << _s << endl;
    ptab();
    cout << "}" << endl;
}

math::VFour Decay::fill(const event &e, int debug) {
    math::VFour p;
    Decay *d;

    for (auto &c: _children) {
        math::VFour v;
        if (c.Stable()) {
            if (debug) {
                cout << "Found stable particle " << c.Name() << endl;
            }
            v = e.getPartPFinal(c.Name(), c.Charge(), c.Index(), debug);
            if (debug) {
                cout << "Setting p of " << c.Name() << " to:" << endl;
                cout << v << endl;
            }
            c.set4P(v);
            p += v;
        } else {
            if (debug) {
                cout << "Found unstable particle " << c.Name() << endl;
            }
            d = c.get_decay();
            if (debug) {
                cout << "Calling fill for " << c.Name() << endl;
            }
            v = d->fill(e, debug);
            if (debug) {
                cout << "Setting p of " << c.Name() << " to:" << endl;
                cout << v << endl;
            }
            c.set4P(v);
            p += v;
        }
    }

    return p;
}

math::VFour *Decay::get4P(particle *part, int debug) {

    math::VFour *p = NULL;

    list<particle>::iterator child;

    if (debug) {
        cout << "looking for " << part->Name() << endl;
    }
    child = this->_children.begin();
    while (child != this->_children.end()) {
        if (debug) {
            cout << "checking against " << child->Name() << endl;
        }
        p = child->get4P(part, debug);
        if (p != NULL) {
            if (debug) {
                cout << "found it! " << child->Name() << " == " << part->Name()
                     << endl;
            }
            break;
        }
        child++;
    }
    return p;
}

Decay &Decay::setupFrames(lorentzTransform T, int debug) {

    // assert(this->_children.size() == 2);
    math::VFour p(0, math::VThree(0, 0, 0));


    // boost children into correct frame
    for (auto &child: _children) {
        if (debug) {
            cout << "boosting child (" << child.Name()
                 << ") into correct frame" << endl;
            cout << "p before: " << endl;
            cout << child.get4P() << endl;
        }
        child *= T;
        if (debug) {
            cout << "p after: " << endl;
            cout << child.get4P() << endl;
        }
    }

    // setupDecay for children
    for (auto &child: _children) {
        p += child.setupFrames(debug);
    }
    this->_mass = p.getMass();

    if (debug) {
        cout << "Decay mass: " << this->_mass << endl;
        cout << "Decay analyzer: " << _children.front().Name()
             << _children.front().Charge() << "["
             << _children.front().Index() << "]" << endl;
        cout << "Decay angles: theta: " << _children.front().get3P().getTheta()
             << " " << "phi: " << _children.front().get3P().getPhi() << endl;
    }

    return *this;
}

complex<double> Decay::expt_amp(double b, double t, int debug) const {
    assert(b >= 0);

    complex<double> a, amp;
    list<particle>::const_iterator child;

    int s1, s2;
    particle child1, child2, analyzer;
    list<int> helicities1;
    list<int> helicities2;
    list<int>::iterator lambda1;
    list<int>::iterator lambda2;

    child = this->_children.begin();
    child1 = *child;
    child++;
    child2 = *child;

    s1 = child1.J();
    s2 = child2.J();
    helicities1 = child1.helicities();
    helicities2 = child2.helicities();
    analyzer = child1;

    addtab();

    if (debug) {
        ptab();
        cout << "My children are: " << child1.Name() << " and " << child2.Name()
             << endl;
        ptab();
        cout << "My childrens spins are: " << s1 << " and " << s2 << endl;
        ptab();
        cout << "child1's helicity list is " << helicities1.size()
             << " elements long" << endl;
        ptab();
        cout << "child1's helicities are: ";
        list<int>::const_iterator i = helicities1.begin();
        while (i != helicities1.end()) {
            cout << *i++ << " ";
        }
        cout << endl;
        ptab();
        cout << "child2's helicity list is " << helicities2.size()
             << " elements long" << endl;
        ptab();
        cout << "child2's helicities are: ";
        i = helicities2.begin();
        while (i != helicities2.end()) {
            cout << *i++ << " ";
        }
        cout << endl;
    }

    // calculate my amplitude

    amp = complex<double>(0, 0);
    a = complex<double>(0, 0);
    lambda1 = helicities1.begin();
    while (lambda1 != helicities1.end()) {
        lambda2 = helicities2.begin();
        while (lambda2 != helicities2.end()) {
            if (debug) {
                ptab();
                cout << "lambda1: " << *lambda1 << " lambda2: " << *lambda2
                     << endl;
                cout << "exp(-" << b << "|" << t << "|)" << endl;
            }
            a = complex<double>(exp(-(b / 2.0) * fabs(t)), 0.0);
            if (debug) {
                ptab();
                cout << "a = " << a << endl;
            }

            if (debug) {
                ptab();
                cout << "calculating amp for " << child1.Name() << endl;
            }
            a *= child1.decayAmp(*lambda1, debug);
            if (debug) {
                ptab();
                cout << "a *= child1.decayAmp = " << a << endl;
            }

            if (debug) {
                ptab();
                cout << "calculating amp for " << child2.Name() << endl;
            }
            a *= child2.decayAmp(*lambda2, debug);
            if (debug) {
                ptab();
                cout << "a *= child2.decayAmp = " << a << endl;
            }

            if (debug) {
                ptab();
                cout << "PLUS" << endl;
            }
            amp += a;
            if (debug) {
                ptab();
                cout << "amp += a = " << amp << endl;
            }
            lambda2++;
        }
        lambda1++;
    }
    subtab();
    return amp;
}

complex<double> Decay::amp(int j, int m, int debug) const {
    assert(j >= 0);
    assert(m <= j);

    complex<double> a, amp;
    int lambda;
    list<particle>::const_iterator child;

    int s1, s2;
    particle child1, child2, analyzer;
    list<int> helicities1;
    list<int> helicities2;
    list<int>::iterator lambda1;
    list<int>::iterator lambda2;

    child = this->_children.begin();
    child1 = *child;
    child++;
    child2 = *child;

    s1 = child1.J();
    s2 = child2.J();
    helicities1 = child1.helicities();
    helicities2 = child2.helicities();
    analyzer = child1;

    addtab();

    if (debug) {
        ptab();
        cout << "My children are: " << child1.Name() << " and " << child2.Name()
             << endl;
        ptab();
        cout << "My childrens spins are: " << s1 << " and " << s2 << endl;
        ptab();
        cout << "child1's helicity list is " << helicities1.size()
             << " elements long" << endl;
        ptab();
        cout << "child1's helicities are: ";
        list<int>::const_iterator i = helicities1.begin();
        while (i != helicities1.end()) {
            cout << *i++ << " ";
        }
        cout << endl;
        ptab();
        cout << "child2's helicity list is " << helicities2.size()
             << " elements long" << endl;
        ptab();
        cout << "child2's helicities are: ";
        i = helicities2.begin();
        while (i != helicities2.end()) {
            cout << *i++ << " ";
        }
        cout << endl;
    }

    // calculate my amplitude

    amp = complex<double>(0, 0);
    a = complex<double>(0, 0);
    lambda1 = helicities1.begin();
    while (lambda1 != helicities1.end()) {
        lambda2 = helicities2.begin();
        while (lambda2 != helicities2.end()) {
            lambda = *lambda1 - *lambda2;
            if (abs(lambda) <= j) {
                if (debug) {
                    ptab();
                    cout << "lambda1: " << *lambda1 << " lambda2: " << *lambda2
                         << endl;
                }
                complex<double> df;
                double CGcoefficient1, CGcoefficient2;
                double tildeFactor;
                double barrierFactor;
                double lambdaFactor;
                double phi = 0.0;
                double theta = 0.0;
                if (this->_children.size() == 2) {
                    phi = analyzer.get3P().getPhi();
                    theta = analyzer.get3P().getTheta();
                } else if (this->_children.size() == 3) {
                    // omega case, use normal to Decay plane
                    math::VThree normal = child1.get3P().crossMultiply(child2.get3P());
                    phi = normal.getPhi();
                    theta = normal.getTheta();
                }

                tildeFactor = tilde(this->_l);
                df = conj(D(phi, theta, 0, j, m, lambda));
                CGcoefficient1 = clebsch(_l, _s, j, 0, lambda, lambda);
                CGcoefficient2 = clebsch(s1, s2, _s, *lambda1, -*lambda2,
                                         lambda);
                if (this->_children.size() == 2) {
                    barrierFactor = F(this->_l, analyzer.get3P().getLen());
                    lambdaFactor = 1;
                } else if (this->_children.size() == 3) {
                    // omega case
                    // instead of barrier factor use lambda factor
                    particle piZero, piPlus, piMinus;
                    child = this->_children.begin();
                    while (child != this->_children.end()) {
                        if ((*child).is("pi0")) {
                            piZero = *child;
                        } else if ((*child).is("pi")) {
                            switch ((*child).Charge()) {
                                case -1:
                                    piMinus = *child;
                                    break;
                                case 1:
                                    piPlus = *child;
                                    break;
                                default:
                                    cerr << "bad child for omega: ";
                                    cerr << (*child).Name() << endl;
                            }
                        } else {
                            cerr << "bad child for omega: ";
                            cerr << (*child).Name() << endl;
                        }
                        child++;
                    }
                    if (debug) {
                        ptab();
                        cout << "calculate sqrt(lambda)" << endl;
                        ptab();
                        cout << "P_pi+: " << piPlus.get3P() << endl;
                        ptab();
                        cout << "P_pi-: " << piMinus.get3P() << endl;
                        ptab();
                        cout << "| P_pi+ X P_pi- |: " \
 << piPlus.get3P().crossMultiply(piMinus.get3P()).getLen() << endl;
                        ptab();
                        cout << "M(pi^+ pi^- pi^0): " << this->_mass << endl;
                        ptab();
                        cout << "M(pi^-): " << piMinus.get4P().getLenSq() << endl;
                        ptab();
                        cout << "numerator: " \
 << piPlus.get3P().crossMultiply(piMinus.get3P()).getLen() << endl;
                        ptab();
                        cout << "denominator: " \
 << (sqrt(3.0 / 4.0) * (pow(this->_mass / 3.0, 2.0) - piMinus.get4P().getLenSq()))
                             << endl;
                    }
                    lambdaFactor = piPlus.get3P().crossMultiply(piMinus.get3P()).getLen() \
 / (sqrt(3.0 / 4.0) * (pow(this->_mass / 3.0, 2.0) - piMinus.get4P().getLenSq()));
                    barrierFactor = 1;
                }

                if (debug) {
                    ptab();
                    cout << "tilde(" << this->_l << "){=" << tildeFactor << "}";
                    cout << "D[" << j << ", " << m << ", " << lambda << "]";
                    cout << "(" << phi << ", " << theta << ", 0){=" << df
                         << "}";
                    cout << "( " << _l << " 0 " << _s << " " << lambda << \
                        " | " << j << " " << lambda << " ){=" << CGcoefficient1
                         << "}";
                    cout << "( " << s1 << " " << *lambda1 << " " << s2 << " "
                         << -*lambda2 << \
                        " | " << _s << " " << lambda << " ){=" << CGcoefficient2
                         << "}";
                    cout << "F_" << this->_l << "(" << analyzer.get3P().getLen()
                         << "){=" << barrierFactor << "}";
                    cout << "sqrt(lambda)" << "{=" << lambdaFactor << "}";
                    cout << "Delta(" << this->_mass << ")";
                    cout << endl;
                }
                a = tildeFactor * df * CGcoefficient1 * CGcoefficient2 *
                    barrierFactor * lambdaFactor;
                if (debug) {
                    ptab();
                    cout << "a = " << a << endl;
                }

                if (debug) {
                    ptab();
                    cout << "calculating amp for " << child1.Name() << endl;
                }
                a *= child1.decayAmp(*lambda1, debug);
                if (debug) {
                    ptab();
                    cout << "a *= child1.decayAmp = " << a << endl;
                }

                if (debug) {
                    ptab();
                    cout << "calculating amp for " << child2.Name() << endl;
                }
                a *= child2.decayAmp(*lambda2, debug);
                if (debug) {
                    ptab();
                    cout << "a *= child2.decayAmp = " << a << endl;
                }

                if (debug) {
                    ptab();
                    cout << "PLUS" << endl;
                }
                amp += a;
                if (debug) {
                    ptab();
                    cout << "amp += a = " << amp << endl;
                }
            }
            lambda2++;
        }
        lambda1++;
    }
    subtab();
    return amp;
}

Decay &Decay::operator*=(const lorentzTransform &L) {
    for (auto &c: _children) {
        c *= L;
    }
    return *this;
}

int particle::_particle_debug = 0;
int Decay::_decay_debug = 0;

int particle::Stable() const {
    return ((this->_decay) ? 0 : 1);
}

math::VFour particle::get4P() const {
    return (this->_p);
}

math::VFour *particle::get4P(particle *part, int debug) const {
    math::VFour *ret = NULL;
    if (this->Name() == part->Name()
        && this->Charge() == part->Charge()
        && this->Index() == part->Index()) {
        ret = new math::VFour(this->_p);
        if (debug) {
            cout << "found particle " << part->Name() << part->Charge() << "["
                 << part->Index() << "]" << endl;
            cout << "returning math::VFour:" << ret << endl;
        }
    } else if (this->get_decay()) {
        if (debug) {
            cout << "I'm a " << this->Name() << this->Charge() << "["
                 << this->Index() << "]" << endl;
            cout << "not a " << part->Name() << part->Charge() << "["
                 << part->Index() << "]" << endl;
            cout << "checking my children..." << endl;
        }
        Decay *d = this->get_decay();
        ret = d->get4P(part, debug);
    }
    return ret;
}

math::VThree particle::get3P() const {
    return (this->_p.getVector());
}

int particle::Index() const {
    return (this->_index);
}

int particle::Charge() const {
    return (this->_charge);
}

list<int> &particle::helicities() {
    if (this->_helicities.size() == 0) {
        for (int lam = -(this->J()); lam <= this->J(); lam += 2) {
            this->_helicities.push_back(lam);
        }
    }
    return (this->_helicities);
}

particle &particle::addHelicity(int lam) {
    this->_helicities.push_back(lam);
    return (*this);
}

int particle::is(string nm) const {
    if (this->Name() == nm) {
        return (1);
    } else {
        return (0);
    }
}

math::VFour particle::setupFrames(int debug) {
    math::VFour ret;
    if (_decay) {
        if (debug) {
            cout << "put ";
            cout << this->Name();
            cout << this->Charge();
            cout << "[" << this->Index() << "]";
            cout << " into helicity frame:" << endl;
            cout << "momentum: " << _p << endl;
        }

        // i should be in my parents rest frame when this is called
        // so this->_p.get3P() should be a breakup momentum
        ret = this->_p;
        // this->_decay->Parent(*this);

        lorentzTransform L, T;
        matrix<double> X(4, 4);
        rotation R;
        math::VThree normal;
        math::VFour tempP;

        tempP = this->_p;

        // make y perpendicular to z_old and p
        normal = math::VThree(0, 0, 1).crossMultiply(tempP.getVector());
        T.set(R.set(normal.getPhi(), normal.getTheta() - M_PI / 2, -M_PI / 2));
        L = T;
        tempP *= T;

        // make z_new parallel to p
        T.set(R.set(0, signof (tempP.getX()) * tempP.getTheta(), 0));
        X = T * L;
        L = *((lorentzTransform *) &X);
        tempP *= T;

        // boost into p rest frame
        T.set(tempP);
        X = T * L;
        L = *((lorentzTransform *) &X);
        tempP *= T;

        this->_decay->setupFrames(L, debug);
        _inRestFrame = 1;
    } else {
        if (debug) {
            cout << "found stable particle ";
            cout << this->Name();
            cout << this->Charge();
            cout << "[" << this->Index() << "]" << endl;
        }
        ret = this->_p;
    }
    return ret;
}

double particle::q() const {
    list<particle>::const_iterator child = this->get_decay()->_children.begin();

    if (_inRestFrame) {
        return child->get3P().getLen();
    } else {
        assert(this->get_decay()->_children.size() == 2);
        particle child1 = *child;
        child++;
        particle child2 = *child;

        double lam;
        double Msq = this->get4P().getLenSq();
        double m1sq = child1.get4P().getLenSq();
        double m2sq = child2.get4P().getLenSq();

        lam = lambda(Msq, m1sq, m2sq);
        return (sqrt(fabs(lam / (4 * Msq))));
    }
}

double particle::q0() const {
    list<particle>::const_iterator child = this->get_decay()->_children.begin();

    assert(this->get_decay()->_children.size() == 2);
    particle child1 = *child;
    child++;
    particle child2 = *child;

    double lam;
    double Msq = pow(this->Mass(), 2);
    double m1sq = child1.get4P().getLenSq();
    double m2sq = child2.get4P().getLenSq();

    lam = lambda(Msq, m1sq, m2sq);
    return (sqrt(fabs(lam / (4 * Msq))));
}


string particle::sprint(string space) {
    string s;
    s = this->Name();
    s += chargetos(this->Charge());
    s += "[";
    s += itos(this->_index);
    s += "]";
    if (this->get_decay()) {
        s += space;
        s += "{";
        list<particle>::iterator c;
        for (c = this->get_decay()->_children.begin();
             c != this->get_decay()->_children.end(); c++) {
            s += space;
            s += c->sprint(space);
        }
        s += space;
        s += itos(this->get_decay()->L());
        s += space;
        s += itos(this->get_decay()->S());
        s += space;
        s += "}";

    }
    return s;
}

complex<double> particle::decayAmp(int lambda, int debug) {
    complex<double> a;
    complex<double> bw;
    if (Stable()) {
        a = complex<double>(1, 0);
    } else {
        bw = this->_massDep->val(*this);
        if (debug) {
            ptab();
            cout << "calculate Decay amplitude for " << this->Name()
                 << this->Charge() << "[" << this->Index()
                 << "] {bw=" << bw << "}" << endl;
        }
        Decay *d = this->get_decay();
        a = d->amp(this->J(), lambda, debug) * bw;
    }
    if (debug) {
        ptab();
        cout << this->Name() << " Decay amp = " << a << endl;
    }
    return a;
}


void particle::print() const {
    this->particleData::print();
    ptab();
    cout << "charge: " << _charge << "\tid: " << _index << endl;
    ptab();
    cout << "momentum: " << _p << endl;
    if (_decay) {
        addtab();
        cout << "mass dependance: ";
        _massDep->print();
        cout << endl;
        _decay->print();
        subtab();
    }
}

void particle::printFrames() const {
    this->particleData::print();
    ptab();
    cout << "charge: " << _charge << "\tid: " << _index << endl;
    ptab();
    cout << "momentum: " << _p << endl;
    if (_decay) {
        addtab();
        _decay->printFrames();
        subtab();
    }
}



