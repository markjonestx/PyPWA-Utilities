//
// Created by Mark Jones on 9/10/22.
//

#pragma once

#include <iostream>
#include <lib/math/VFour.h>
#include <lib/math/VThree.h>

namespace math {

    std::ostream &operator<<(std::ostream &, const VFour &);
    std::ostream &operator<<(std::ostream &, const VThree &);


} // math
