#include "waveH.h"

#include <lib/math/VFour.h>
#include <lib/math/io.h>


wave::wave(const wave &wv) : particle(wv) {
    this->_m = wv._m;
    this->_epsilon = wv._epsilon;
    this->_beam = wv._beam;
    this->_target = wv._target;
    this->_t = wv._t;
    this->_b = wv._b;
}

wave &wave::operator=(const wave &wv) {
    particle::operator=(wv);
    this->_m = wv._m;
    this->_epsilon = wv._epsilon;
    this->_beam = wv._beam;
    this->_target = wv._target;
    this->_t = wv._t;
    this->_b = wv._b;
    return (*this);
}

void wave::print() const {
    cout << "wave: " << endl;
    cout << "beam: " << _beam << endl;
    cout << "target: " << _target << endl;
    if (this->_b == 0.0) {
        cout << this->I() << getsign(this->G());
        cout << "(" << this->J() << getsign(this->P()) << getsign(this->C())
             << ")";
        cout << " m = " << this->_m << " eps = " << this->_epsilon << endl;
    } else {
        cout << "b: " << this->_b << " t: " << this->_t;
    }
    cout << "momentum: " << get4P() << endl;
    cout << "decays to:" << endl;
    if (!(this->Stable())) this->get_decay()->print();
}

wave wave::fill(const event &e, int debug) {
    Decay *d;
    math::VFour p;

    if (debug) {
        cout << "Filling wave" << endl;
        cout << "Setting beam p to:" << endl;
        cout << e.beam().get4P() << endl;
    }
    this->_beam = e.beam().get4P();
    if (debug) {
        cout << "Setting target p to:" << endl;
        cout << e.target().get4P() << endl;
    }
    this->_target = e.target().get4P();
    d = this->get_decay();
    if (debug) {
        cout << "Calling fill for wave: " << endl;
    }
    p = d->fill(e, debug);
    if (debug) {
        cout << "Setting p of wave to: " << endl;
        cout << p << endl;
    }
    this->set4P(p);

    return *this;
}

wave &wave::setupFrames(int debug) {
    if (debug) {
        cout << "put wave into HELICITY frame:" << endl;
    }
    if (Stable()) { ;
    } else {
        lorentzTransform L, T;
        matrix<double> X(4, 4);
        rotation R;
        math::VFour tempX, tempX_copy, tempBeam, tempTarget, tempChild1;
        list<particle>::iterator child1 = this->get_decay()->_children.begin();
        math::VThree N;

        tempX = this->get4P();
        tempBeam = this->_beam;
        tempTarget = this->_target;
        tempChild1 = (*child1).get4P();

        if (debug) {
            cout << "initially in lab:" << endl;
            cout << "tempX: " << tempX << endl;
            cout << "tempBeam: " << tempBeam << endl;
            cout << "tempTarget: " << tempTarget << endl;
            cout << "tempChild1: " << tempChild1 << endl;
        }

        if (this->channel() == "t") {
            // put normal to production plane along y
            N = tempBeam.getVector().crossMultiply(tempX.getVector());
        } else if (this->channel() == "s" || this->channel() == "expt") {
            // use lab y
            N = math::VThree(0.0, 1.0, 0.0);
        }
        T.set(R.set(N.getPhi(), N.getTheta() - M_PI / 2.0, -M_PI / 2.0));
        L = T;
        tempX *= T;
        tempBeam *= T;
        tempTarget *= T;
        tempChild1 *= T;
        if (debug) {
            cout << "put normal to PP along y:" << endl;
            cout << "tempX: " << tempX << endl;
            cout << "tempBeam: " << tempBeam << endl;
            cout << "tempTarget: " << tempTarget << endl;
            cout << "tempChild1: " << tempChild1 << endl;
        }

        // boost to CM frame
        T.set(tempBeam + tempTarget);
        X = T * L;
        L = *((lorentzTransform *) &X);
        tempX *= T;
        tempBeam *= T;
        tempTarget *= T;
        tempChild1 *= T;
        if (debug) {
            cout << "boost to CM frame:" << endl;
            cout << "tempX: " << tempX << endl;
            cout << "tempBeam: " << tempBeam << endl;
            cout << "tempTarget: " << tempTarget << endl;
            cout << "tempChild1: " << tempChild1 << endl;
        }

        //tempX_copy = new math::VFour();
        tempX_copy = math::VFour(tempX.getT(), tempX.getVector());

        // boost to X rest frame
        T.set(tempX);
        X = T * L;
        L = *((lorentzTransform *) &X);
        tempX *= T;
        tempBeam *= T;
        tempTarget *= T;
        tempChild1 *= T;
        if (debug) {
            cout << "boost to XRF:" << endl;
            cout << "tempX: " << tempX << endl;
            cout << "tempBeam: " << tempBeam << endl;
            cout << "tempTarget: " << tempTarget << endl;
            cout << "tempChild1: " << tempChild1 << endl;
        }

        // put z along X
        // T.set (R.set (tempBeam.V ()));
        T.set(R.set(0.0, signof(tempX_copy.getX()) * tempX_copy.getVector().getTheta(), 0.0));
        X = T * L;
        L = *((lorentzTransform *) &X);
        tempX *= T;
        tempBeam *= T;
        tempTarget *= T;
        tempChild1 *= T;
        if (debug) {
            cout << "put beam along z:" << endl;
            cout << "tempX: " << tempX << endl;
            cout << "tempBeam: " << tempBeam << endl;
            cout << "tempTarget: " << tempTarget << endl;
            cout << "tempChild1: " << tempChild1 << endl;
        }

        // boost the beam and the target
        this->_beam *= L;
        this->_target *= L;

        // setupFrames of children
        this->get_decay()->setupFrames(L, debug);
    }
    return *this;
}

complex<double> wave::decayAmp(int debug) {
    complex<double> a;
    if (Stable()) {
        a = complex<double>(1, 0);
    } else if (this->_b != 0.0) {
        if (debug) {
            cout << "calculate Decay amplitude for expt wave b=" << this->_b;
            cout << " t=" << this->_t << endl;
        }
        Decay *d = this->get_decay();
        a = d->expt_amp(this->_b, this->_t, debug);
    } else {
        if (debug) {
            cout << "calculate Decay amplitude for wave J=" << this->J();
            cout << " m=" << this->_m << endl;
        }
        Decay *d = this->get_decay();
        a = d->amp(this->J(), this->_m, debug);
    }
    return a;
}

wave &wave::operator*=(const lorentzTransform &L) {
    this->_beam *= L;
    this->_target *= L;
    particle::operator*=(L);
    return *this;
}
