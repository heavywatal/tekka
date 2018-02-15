/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"
#include "segment.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>
#include <wtl/exception.hpp>
#include <sfmt.hpp>

namespace pbt {

Population::Population(const size_t initial_size, std::random_device::result_type seed)
: engine_(std::make_unique<URBG>(seed)) {HERE;
    Individual::set_default_values();
    const size_t half = initial_size / 2UL;
    const size_t rest = initial_size - half;
    males_.reserve(half);
    females_.reserve(rest);
    auto founder = std::make_shared<Individual>();
    for (size_t i=0; i<half; ++i) {males_.emplace_back(new Individual(founder, founder, 0));}
    for (size_t i=0; i<rest; ++i) {females_.emplace_back(new Individual(founder, founder, 0));}
}

Population::~Population() {} // to allow forward declaration of Individual

void Population::run(const uint_fast32_t simulating_duration,
                     const double sample_rate,
                     const uint_fast32_t recording_duration) {HERE;
    auto recording_start = simulating_duration - recording_duration;
    for (year_ = 4u; year_ < simulating_duration; ++year_) {
        DCERR(year_ << ": " << sizes() << std::endl);
        reproduce();
        survive(0u);
        survive(1u);
        survive(2u);
        survive(3u);
        if (year_ >= recording_start) {
            sample(sample_rate);
        }
        migrate();
    }
    DCERR(year_ << ": " << sizes() << std::endl);
}

void Population::reproduce() {
    std::vector<std::shared_ptr<Individual>> boys;
    std::vector<std::shared_ptr<Individual>> girls;
    std::vector<std::vector<std::shared_ptr<Individual>>> males_located(Individual::num_locations());
    for (const auto& p: males_) {
        males_located[p->location()].emplace_back(p);
    }
    for (const auto& mother: females_) {
        if (!mother->is_in_breeding_place()) continue;
        const auto& potential_fathers = males_located[mother->location()];
        if (potential_fathers.empty()) continue;
        const uint_fast32_t num_juveniles = mother->recruitment(year_, *engine_);
        const std::shared_ptr<Individual> father = *wtl::choice(potential_fathers.begin(), potential_fathers.end(), *engine_);
        // TODO: multiple fathers
        for (uint_fast32_t i=0; i<num_juveniles; ++i) {
            if (wtl::generate_canonical(*engine_) < 0.5) {
                boys.emplace_back(new Individual(father, mother, year_));
            } else {
                girls.emplace_back(new Individual(father, mother, year_));
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

void Population::survive(const uint_fast32_t quarter) {
    auto impl = [quarter,this](const decltype(males_)& v) {
        decltype(males_) survivors;
        survivors.reserve(v.size());
        std::copy_if(v.begin(), v.end(), std::back_inserter(survivors),
                     [&](const auto& p) {return p->has_survived(year_, quarter, *engine_);});
        return survivors;
    };
    males_ = impl(males_);
    females_ = impl(females_);
}

void Population::migrate() {
    for (auto& p: males_) {p->migrate(year_, *engine_);}
    for (auto& p: females_) {p->migrate(year_, *engine_);}
}

void Population::sample(const double rate) {
    auto impl = [rate,this](decltype(males_)* individuals) {
        decltype(males_) survivors;
        survivors.reserve(individuals->size());
        std::vector<std::vector<size_t>> adults(Individual::num_breeding_places());
        std::vector<std::vector<size_t>> juveniles(Individual::num_breeding_places());
        size_t num_adults = 0u;
        for (size_t i=0; i<individuals->size(); ++i) {
            auto& p = individuals->at(i);
            if (p->is_in_breeding_place()) {
                if (p->birth_year() == year_) {
                    juveniles[p->location()].emplace_back(i);
                } else {
                    adults[p->location()].emplace_back(i);
                    ++num_adults;
                }
            } else {
                survivors.emplace_back(p);
            }
        }
        decltype(males_) samples;
        samples.reserve(static_cast<size_t>(std::round(3.0 * rate * num_adults)));
        auto sort = [&](const std::vector<size_t>& indices, const size_t k) {
            const auto chosen = wtl::sample(indices.size(), k, *engine_);
            for (size_t i=0; i<indices.size(); ++i) {
                if (chosen.find(i) == chosen.end()) {
                    survivors.emplace_back(individuals->at(indices[i]));
                } else {
                    samples.emplace_back(individuals->at(indices[i]));
                }
            }
        };
        for (uint_fast32_t loc=0; loc<Individual::num_breeding_places(); ++loc) {
            const size_t n_adults = adults[loc].size();
            const size_t n_adult_samples = static_cast<size_t>(std::round(rate * n_adults));
            const size_t n_juvenile_samples = 2u * n_adult_samples;
            sort(adults[loc], n_adult_samples);
            sort(juveniles[loc], n_juvenile_samples);
        }
        individuals->swap(survivors);
        return samples;
    };
    auto m_samples = impl(&males_);
    auto f_samples = impl(&females_);
    auto& samples = year_samples_[year_];
    samples.reserve(m_samples.size() + f_samples.size());
    std::copy(m_samples.begin(), m_samples.end(), std::back_inserter(samples));
    std::copy(f_samples.begin(), f_samples.end(), std::back_inserter(samples));
}

std::ostream& Population::write_sample_header(std::ostream& ost) const {
    if (year_samples_.empty()) return ost;
    return wtl::join(Individual::names(), ost, "\t") << "\tcapture_year\n";
}

std::ostream& Population::write_sample(std::ostream& ost) const {
    write_sample_header(ost);
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            ost << *p << "\t" << ys.first << "\n";
        }
    }
    return ost;
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    std::unordered_map<const Individual*, uint_fast32_t> id_year;
    std::set<const Individual*, Individual::less> nodes;
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            p->trace_back(&nodes);
            id_year.emplace(p.get(), ys.first);
        }
    }
    write_sample_header(ost);
    for (const auto& p: nodes) {
        ost << *p << "\t";
        auto it = id_year.find(p);
        if (it != id_year.end()) {
            ost << it->second;
        }
        ost << "\n";
    }
    return ost;
}

//! Generate [0,1) doubles; size ~ Poisson(lambda)
class RunifPoisson {
  public:
    RunifPoisson(double lambda): poisson_(lambda) {}
    std::vector<double> operator()(URBG& engine) {
        uint_fast32_t n = poisson_(engine);
        std::vector<double> v(n);
        for (uint_fast32_t i=0; i<n; ++i) {
            v[i] = wtl::generate_canonical(engine);
        }
        return v;
    }
  private:
    std::poisson_distribution<uint_fast32_t> poisson_;
};

std::ostream& Population::write_ms(double lambda, std::ostream& ost) const {
    std::set<Segment, Segment::less> nodes;
    std::vector<const Segment*> sampled_nodes;
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            for (bool b: {true, false}) {
                auto it = nodes.emplace(p.get(), b).first;
                sampled_nodes.push_back(&*it);
            }
        }
    }
    std::vector<const Segment*> current(sampled_nodes);
    RunifPoisson runif(lambda);
    uint_fast32_t num_coalesced = 0u;
    uint_fast32_t num_uncoalesced = 0u;
    while (!current.empty()) {
        std::vector<const Segment*> next;
        for (auto& x: current) {
            x->set_mutations(runif(*engine_));
            if (x->is_first_gen()) {
                ++num_uncoalesced;
                continue;
            }
            auto p = nodes.emplace(x->ancestor(), wtl::generate_canonical(*engine_) < 0.5);
            auto ancp = const_cast<Segment*>(&*p.first);
            x->set_ancestral_segment(ancp);
            if (p.second) {
                next.push_back(ancp);
            } else {
                ++num_coalesced;
            }
        }
        current.clear();
        next.swap(current);
    }
    DCERR("uncoalesced: " << num_uncoalesced << " / " << sampled_nodes.size() << "\n");
    WTL_ASSERT(num_coalesced + num_uncoalesced == sampled_nodes.size());

    std::set<double> positions_collector;
    for (const auto& s: nodes) {
        positions_collector.insert(s.begin(), s.end());
    }
    std::vector<double> positions(positions_collector.begin(), positions_collector.end());

    ost << "\n//\n"
        << "segsites: " << positions.size() << "\n"
        << "positions: ";
    wtl::join(positions, ost, " ") << "\n";
    for (const auto& p: sampled_nodes) {
        wtl::join(p->binary_genotype(positions), ost, "") << "\n";
    }
    return ost;
}

std::vector<size_t> Population::sizes() const {
    std::vector<size_t> counter(Individual::num_locations(), 0u);
    for (const auto& p: males_) {++counter[p->location()];}
    for (const auto& p: females_) {++counter[p->location()];}
    return counter;
}

std::ostream& Population::write(std::ostream& ost) const {
    for (const auto& p: males_) {ost << *p << "\n";}
    for (const auto& p: females_) {ost << *p << "\n";}
    return ost;
}

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbt
