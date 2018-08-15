// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <blackthunnus/program.hpp>

//' Run C++ simulation
//' @rdname blackthunnus
// [[Rcpp::export]]
Rcpp::CharacterVector cpp_blackthunnus(const std::vector<std::string>& args) {
    try {
        pbt::Program program(args);
        program.run();
        return Rcpp::CharacterVector::create(
            Rcpp::Named("sample_family", program.sample_family()),
            Rcpp::Named("demography", program.demography())
        );
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
