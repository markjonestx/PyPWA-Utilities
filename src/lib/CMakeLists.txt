
add_subdirectory(math)

add_library(lib
        event.cc
        lorentz.cc
        particle.cc
        particleData.cc
        pputil.cc
        massDep.cc
)


target_link_libraries(lib
        PUBLIC
        libMath
)
