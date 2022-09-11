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

rotation::rotation(
        double alpha, double beta, double gamma
    ) : matrix<double>(3,3)
{
    this->set(alpha, beta, gamma);
}

rotation rotation::set(const math::VThree &V) {
    this->set(V.getPhi(), V.getTheta(), 0.0);
    return (*this);
}

math::VThree operator*=(math::VThree &V, const rotation &R) {
    V = R * V;
    return V;
}

lorentzTransform lorentzTransform::set(const rotation &r) {
    this->zero();
    for (int row = 1; row < 4; row++) {
        for (int col = 1; col < 4; col++) {
            this->el(row, col) = (const_cast<rotation *>(&r))->el(
                    row - 1, col - 1
            );
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


lorentzTransform lorentzTransform::set(const math::VThree &beta) {
    double gamma;
    double gFactor;

    this->_gamma = 1.0 / sqrt(1 - beta.getLenSq());
    gamma = this->_gamma;
    gFactor = pow(gamma, 2) / (gamma + 1);

    this->el(0, 0) = gamma;
    this->el(0, 1) = gamma * beta.getX();
    this->el(0, 2) = gamma * beta.getY();
    this->el(0, 3) = gamma * beta.getZ();

    this->el(1, 1) = (pow(beta.getX(), 2) * gFactor) + 1;
    this->el(1, 2) = beta.getX() * beta.getY() * gFactor;
    this->el(1, 3) = beta.getX() * beta.getZ() * gFactor;

    this->el(2, 2) = (pow(beta.getY(), 2) * gFactor) + 1;
    this->el(2, 3) = beta.getY() * beta.getZ() * gFactor;

    this->el(3, 3) = (pow(beta.getZ(), 2) * gFactor) + 1;

    this->el(1, 0) = this->el(0, 1);
    this->el(2, 0) = this->el(0, 2);
    this->el(2, 1) = this->el(1, 2);
    this->el(3, 0) = this->el(0, 3);
    this->el(3, 1) = this->el(1, 3);
    this->el(3, 2) = this->el(2, 3);

    return (*this);
}

lorentzTransform::lorentzTransform(const math::VThree &beta) : matrix<double>(4,
                                                                          4) {
    this->set(beta);
}


lorentzTransform lorentzTransform::set(const math::VFour &p) {
    math::VThree beta;

    beta.setX(-p.getX() / p.getT());
    beta.setY(-p.getY() / p.getT());
    beta.setZ(-p.getZ() / p.getT());

    this->set(beta);

    return (*this);
}

lorentzTransform::lorentzTransform(const math::VFour &p) : matrix<double>(4, 4) {
    this->set(p);
}

math::VFour operator*=(math::VFour &v, const lorentzTransform &L) {  //new
    v = L * v;
    return v;
}

//	template class matrix<double>;


