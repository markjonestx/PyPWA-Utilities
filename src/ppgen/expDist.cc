#include "plib.h"

#include <cmath>


double expDist(double slope) {
    double r = randm(0, 1);
    double val = -log(r) / slope;
    return (val);
}
