#include <lorentz.h>

rotation rotation::set(double alpha, double beta, double gamma) {
    double ca, sa, cb, sb, cg, sg;

    ca = cos(alpha);
    sa = sin(alpha);
    cb = cos(beta);
    sb = sin(beta);
    cg = cos(gamma);
    sg = sin(gamma);

    this->el(0, 0) = ca * cb * cg - sa * sg;
    this->el(0, 1) = cb * cg * sa + ca * sg;
    this->el(0, 2) = -cg * sb;

    this->el(1, 0) = -sg * cb * ca - cg * sa;
    this->el(1, 1) = -sg * cb * sa + cg * ca;
    this->el(1, 2) = sb * sg;

    this->el(2, 0) = ca * sb;
    this->el(2, 1) = sa * sb;
    this->el(2, 2) = cb;

    return (*this);
}

rotation::rotation(double alpha, double beta, double gamma) : matrix<double>(3,
                                                                             3) {
    this->set(alpha, beta, gamma);
}

rotation rotation::set(const threeVec &V) {
    this->set(V.phi(), V.theta(), 0.0);
    return (*this);
}

threeVec operator*=(threeVec &V, const rotation &R) {
    V = R * V;
    return V;
}

lorentzTransform lorentzTransform::set(const rotation &r) {
    this->zero();
    for (int row = 1; row < 4; row++) {
        for (int col = 1; col < 4; col++) {
//				this->el(row,col) = r.el(row-1,col-1);
            this->el(row, col) = (const_cast<rotation *>(&r))->el(row - 1,
                                                                  col - 1);
        }
    }
    this->el(0, 0) = 1;
    return (*this);
}

lorentzTransform
lorentzTransform::set(double alpha, double beta, double gamma) {
    rotation r(alpha, beta, gamma);
    this->zero();
    this->set(r);
    return (*this);
}

lorentzTransform::lorentzTransform(double alpha, double beta, double gamma)
        : matrix<double>(4, 4) {
    this->set(alpha, beta, gamma);
}

lorentzTransform::lorentzTransform(const rotation &r) : matrix<double>(4, 4) {
    this->set(r);
}


lorentzTransform lorentzTransform::set(const threeVec &beta) {
    double gamma;
    double gFactor;

    this->_gamma = 1.0 / sqrt(1 - beta.lenSq());
    gamma = this->_gamma;
    gFactor = pow(gamma, 2) / (gamma + 1);

    this->el(0, 0) = gamma;
    this->el(0, 1) = gamma * beta.x();
    this->el(0, 2) = gamma * beta.y();
    this->el(0, 3) = gamma * beta.z();

    this->el(1, 1) = (pow(beta.x(), 2) * gFactor) + 1;
    this->el(1, 2) = beta.x() * beta.y() * gFactor;
    this->el(1, 3) = beta.x() * beta.z() * gFactor;

    this->el(2, 2) = (pow(beta.y(), 2) * gFactor) + 1;
    this->el(2, 3) = beta.y() * beta.z() * gFactor;

    this->el(3, 3) = (pow(beta.z(), 2) * gFactor) + 1;

    this->el(1, 0) = this->el(0, 1);
    this->el(2, 0) = this->el(0, 2);
    this->el(2, 1) = this->el(1, 2);
    this->el(3, 0) = this->el(0, 3);
    this->el(3, 1) = this->el(1, 3);
    this->el(3, 2) = this->el(2, 3);

    return (*this);
}

lorentzTransform::lorentzTransform(const threeVec &beta) : matrix<double>(4,
                                                                          4) {
    this->set(beta);
}


lorentzTransform lorentzTransform::set(const fourVec &p) {
    threeVec beta;

    beta.el(0) = -p.x() / p.t();
    beta.el(1) = -p.y() / p.t();
    beta.el(2) = -p.z() / p.t();

    this->set(beta);

    return (*this);
}

lorentzTransform::lorentzTransform(const fourVec &p) : matrix<double>(4, 4) {
    this->set(p);
}

fourVec operator*=(fourVec &v, const lorentzTransform &L) {  //new
    v = L * v;
    return v;
}

//	template class matrix<double>;


