include(FetchContent)

# Include the logger
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.10.0
)
FetchContent_MakeAvailable(spdlog)


# Include Argument Parser
FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG v3.0.0
)
FetchContent_MakeAvailable(cxxopts)

# Include Eigen Library
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)
