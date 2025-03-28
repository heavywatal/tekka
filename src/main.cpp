#include "program.hpp"

#include <iostream>
#include <exception>

//! Just instantiate and run pbf::Program
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        pbf::Program program(arguments);
        program.run();
        program.write();
    } catch (const pbf::exit_success&) {
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
