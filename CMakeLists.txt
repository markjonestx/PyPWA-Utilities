cmake_minimum_required(VERSION 3.1)
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

# These features should only be enabled when this is the project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

	# Support IDE Folders
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	# Here you can enable documentation and testing
endif()


add_subdirectory(src)

