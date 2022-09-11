//
// Created by Mark Jones on 9/10/22.
//

#include "io.h"
#include <shared/math/VFour.h>
#include <shared/math/VThree.h>


namespace math {

    std::ostream& operator<<(std::ostream &out, const VFour &v) {
        out << v.getT()
            << "\t" << v.getX() << "\t" << v.getY() << "\t" << v.getZ();
        return out;
    }

    std::ostream& operator<<(std::ostream &out, const VThree &v) {
        out << v.getX() << "\t" << v.getY() << "\t" << v.getZ();
        return out;
    }

} // math