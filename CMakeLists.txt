cmake_minimum_required(VERSION 3.20)
project(PWA_UTILS VERSION 1.1
		DESCRIPTION "A Modern revision of PWA2000"
		LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

include_directories(src)

if (MSVC)
	# warning level 4 and all warnings as errors
	add_compile_options(/W4 /WX)
else()
	# lots of warnings and all warnings as errors
	add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

# Support IDE Folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


add_subdirectory(dependencies)
add_subdirectory(src)

