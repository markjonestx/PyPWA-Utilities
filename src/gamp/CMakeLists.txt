
set(GAMP_SOURCES
        gamp.cc
        keyfile.cc
        keyParse.cc
        keyScan.cc
)

add_executable(gamp ${GAMP_SOURCES} wave.cc)
add_executable (hamp ${GAMP_SOURCES} waveH.cc)

target_link_libraries(gamp lib)
target_link_libraries(hamp lib)

install(TARGETS gamp DESTINATION bin)
install(TARGETS hamp DESTINATION bin)
