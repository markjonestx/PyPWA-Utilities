cmake_minimum_required(VERSION 3.1)
project(PWA_UTILS VERSION 1.1
		DESCRIPTION "A Modern revision of PWA2000"
		LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)
include_directories(src)

set(SHARED_LIBPP
		src/event.cc
		src/lorentz.cc
		src/particle.cc
		src/particleData.cc
		src/pputil.cc
		src/math/VFour.cpp
		src/math/VThree.cpp
		src/math/io.cpp
)

set (
    GAMP_SOURCES
	src/gamp.cc
	src/keyfile.cc
	src/keyParse.cc
	src/keyScan.cc
	src/massDep.cc
)

set (
    PPGEN_SOURCES
	src/expDist.cc
	src/ppgen.cc
	src/txtEvent.cc
)

if (MSVC)
	# warning level 4 and all warnings as errors
	add_compile_options(/W4 /WX)
else()
	# lots of warnings and all warnings as errors
	add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

add_library (pp STATIC ${SHARED_LIBPP})

add_executable (gamp ${GAMP_SOURCES} src/wave.cc)
add_executable (hgamp ${GAMP_SOURCES} src/waveH.cc)
add_executable (ppgen ${PPGEN_SOURCES})
add_executable (vamp src/vamp.cc)

target_link_libraries (gamp pp)
target_link_libraries (hgamp pp)
target_link_libraries (ppgen pp)
target_link_libraries (vamp)

install (TARGETS gamp DESTINATION bin)
install (TARGETS hgamp DESTINATION bin)
install (TARGETS ppgen DESTINATION bin)
install (TARGETS vamp DESTINATION bin)

