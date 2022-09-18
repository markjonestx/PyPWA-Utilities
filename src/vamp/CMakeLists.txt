
add_executable(vamp vamp.cpp Arguments.cpp)
target_precompile_headers(vamp PRIVATE
        <complex>
        <iostream>
        <string>
)
target_link_libraries(vamp PRIVATE cxxopts)
install(TARGETS vamp DESTINATION bin)
