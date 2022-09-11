
set(PPGEN_SOURCES
        expDist.cc
        ppgen.cc
        txtEvent.cc
)

add_executable(ppgen ${PPGEN_SOURCES})
target_link_libraries(ppgen shared)

install(TARGETS ppgen DESTINATION bin)
