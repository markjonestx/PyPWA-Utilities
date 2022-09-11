
add_subdirectory(math)

add_library(shared
        event.cc
        lorentz.cc
        particle.cc
        particleData.cc
        pputil.cc
        massDep.cc
)


target_link_libraries(shared
        PUBLIC
        math
)
