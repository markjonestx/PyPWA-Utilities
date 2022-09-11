
add_subdirectory(shared)
add_subdirectory(gamp)
add_subdirectory(ppgen)

add_executable(vamp vamp.cc)
install(TARGETS vamp DESTINATION bin)
