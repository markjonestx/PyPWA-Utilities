cmake_minimum_required (VERSION 2.8.12)
project (PWA_UTILS)

# find_package (ROOT REQUIRED COMPONENTS MINUIT)
# include (${ROOT_USE_FILE})

include_directories (include)

set (
    SHARED_LIBPP
	src/event.cc
	src/lorentz.cc
	src/particle.cc
	src/particleData.cc
	src/pputil.cc
	src/Vec.cc
)

set (
    GAMP_SOURCES
	src/keyfile.cc
	src/keyParse.cc
	src/keyScan.cc
	src/massDep.cc
	src/gamp.cc
)

set (
    PPGEN_SOURCES
	src/expDist.cc
	src/matrix.cc
	src/ppgen.cc
	src/txtEvent.cc
)

set (
    SHARED_FIT
	src/code.cc
	src/fitlog.cc
	src/fitparse.cc
	src/init.cc
	src/integral.cc
	src/lex.cc
	src/math.cc
	src/parseutils.cc
	src/symbol.cc
)

set (
    FIT_SOURCES
	src/fit.cc
	src/fort_open.f
	src/minuit.cc
	src/matrix.cc
)


add_library (pp ${SHARED_LIBPP})
add_library (fits ${SHARED_FIT})

# Create executables
# add_executable (fit ${FIT_SOURCES})
add_executable (gamp ${GAMP_SOURCES} src/wave.cc)
add_executable (hgamp ${GAMP_SOURCES} src/waveH.cc)
add_executable (integral src/int.cc src/integral.cc)
add_executable (logview src/logview.cc)
add_executable (ppgen ${PPGEN_SOURCES})
add_executable (vamp src/vamp.cc)

# target_link_libraries (fit fits)
target_link_libraries (fit pp)
target_link_libraries (gamp pp)
target_link_libraries (hgamp pp)
target_link_libraries (logview fits)
target_link_libraries (ppgen pp)

