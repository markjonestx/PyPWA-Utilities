//
// Created by Mark Jones on 9/16/22.
//

#include "Arguments.hpp"

namespace vamp {

    Arguments::Arguments(int argc, char **argv) {
        options = new cxxopts::Options("VAMP", "Converts amplitude files to binary and vice versa");
        options->add_options()
            ("m,max_events", "Max number of events to convert",
                 cxxopts::value<int>()->default_value("-1"))
            ("b,binary", "Changes the output file to binary",
                 cxxopts::value<bool>()->default_value("false"))
            ("h,help", "Prints this help message",
                 cxxopts::value<bool>()->default_value("false"));

        this->result = options->parse(argc, argv);
    }

    int Arguments::getMaxEvents() {
        return result["max_events"].as<int>();
    }

    bool Arguments::getBinary() {
        return result["binary"].as<bool>();
    }

    bool Arguments::askedForHelp() {
        return result["help"].as<bool>();
    }

    void Arguments::printHelp() {
        std::cout << options->help() << std::endl;
    }
}
