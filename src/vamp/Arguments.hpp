//
// Created by Mark Jones on 9/16/22.
//

#pragma once
#include <cxxopts.hpp>

namespace vamp {

    class Arguments {
    private:

        cxxopts::Options *options;
        cxxopts::ParseResult result;

    public:
        Arguments(int argc, char **argv);
        ~Arguments() { delete options; }

        int getMaxEvents();
        bool getBinary();
        bool askedForHelp();
        void printHelp();
    };
} // vamp
