// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <blackthunnus/program.hpp>

//' Run C++ simulation
//' @return conf and population as strings
//' @rdname blackthunnus
// [[Rcpp::export]]
std::vector<std::string> cpp_blackthunnus(const std::vector<std::string>& args) {
    pbt::Program program(args);
    std::stringbuf bypass;
    std::streambuf* original_buf = std::cout.rdbuf(&bypass);
    program.run();
    std::cout.rdbuf(original_buf);
    return {bypass.str()};
}
