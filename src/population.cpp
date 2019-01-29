/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"
#include "segment.hpp"
#include "random.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/exception.hpp>

namespace pbf {

Population::Population(const size_t initial_size, std::random_device::result_type seed)
: engine_(std::make_unique<URBG>(seed)) {HERE;
    const size_t half = initial_size / 2UL;
    const size_t rest = initial_size - half;
    males_.reserve(half);
    females_.reserve(rest);
    for (size_t i=0; i<half; ++i) {males_.emplace_back(std::make_shared<Individual>());}
    for (size_t i=0; i<rest; ++i) {females_.emplace_back(std::make_shared<Individual>());}
}

Population::~Population() = default;

void Population::run(const uint_fast32_t simulating_duration,
                     const double sample_rate,
                     const uint_fast32_t recording_duration) {HERE;
    auto recording_start = simulating_duration - recording_duration;
    for (year_ = 4u; year_ < simulating_duration; ++year_) {
        append_demography();
        DCERR(year_ << ": " << demography_.at(year_) << std::endl);
        reproduce();
        survive(0u, false);
        survive(1u, false);
        survive(2u, false);
        survive(3u, true);
        if (year_ >= recording_start) {
            sample(sample_rate);
        }
        migrate();
    }
    append_demography();
    DCERR(year_ << ": " << demography_.at(year_) << std::endl);
}

void Population::reproduce() {
    size_t num_potential_mothers = 0;
    for (const auto& mother: females_) {
        if (mother->is_in_breeding_place()) ++num_potential_mothers;
    }
    size_t num_potential_fathers = 0;
    std::vector<std::vector<std::shared_ptr<Individual>>> males_located(Individual::num_breeding_places());
    for (const auto& p: males_) {
        if (p->is_in_breeding_place()) {
            ++num_potential_fathers;
            males_located[p->location()].emplace_back(p);
        }
    }
    std::vector<std::discrete_distribution<unsigned>> mate_distrs;
    mate_distrs.reserve(males_located.size());
    for (size_t i=0; i<males_located.size(); ++i) {
        std::vector<double> fitnesses;
        fitnesses.reserve(males_located[i].size());
        for (const auto& male: males_located[i]) {
            fitnesses.push_back(male->weight(year_));
        }
        mate_distrs.emplace_back(fitnesses.begin(), fitnesses.end());
    }
    const size_t popsize = num_potential_mothers + num_potential_fathers;
    std::vector<std::shared_ptr<Individual>> boys;
    std::vector<std::shared_ptr<Individual>> girls;
    boys.reserve(num_potential_fathers * 60u);
    girls.reserve(num_potential_fathers * 60u);
    for (const auto& mother: females_) {
        if (!mother->is_in_breeding_place()) continue;
        const auto& potential_fathers = males_located[mother->location()];
        if (potential_fathers.empty()) continue;
        auto& mate_distr = mate_distrs[mother->location()];
        const uint_fast32_t num_juveniles = mother->recruitment(year_, popsize, *engine_);
        for (uint_fast32_t i=0; i<num_juveniles; ++i) {
            const auto father = potential_fathers[mate_distr(*engine_)];
            if (wtl::generate_canonical(*engine_) < 0.5) {
                boys.emplace_back(std::make_shared<Individual>(father, mother, year_));
            } else {
                girls.emplace_back(std::make_shared<Individual>(father, mother, year_));
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

void Population::survive(const uint_fast32_t quarter, bool shrink) {
    auto has_survived = [quarter,this](const auto& p) -> bool {
        return p && p->has_survived(year_, quarter, *engine_);
    };
    if (shrink) {
        auto impl = [&has_survived](decltype(males_)* v) {
            decltype(males_) survivors;
            survivors.reserve(v->size() / 4u);
            std::copy_if(v->begin(), v->end(), std::back_inserter(survivors), has_survived);
            survivors.swap(*v);
        };
        impl(&males_);
        impl(&females_);
    } else {
        for (auto& p: males_) {
            if (!has_survived(p)) p.reset();
        }
        for (auto& p: females_) {
            if (!has_survived(p)) p.reset();
        }
    }
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
        auto sort = [this,&individuals,&survivors,&samples]
                    (const std::vector<size_t>& indices, const size_t k) {
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

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    if (year_samples_.empty()) return ost;
    wtl::join(Individual::names(), ost, "\t") << "\tcapture_year\n";
    std::map<const Individual*, size_t> ids;
    ids.emplace(nullptr, 0u);
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            p->trace_back(ost, &ids, ys.first);
        }
    }
    return ost;
}

//! Generate [0,1) doubles; size ~ Poisson(lambda)
//! @cond
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

struct less_shptr_Segment {
    bool operator()(const std::shared_ptr<Segment>& lhs, const std::shared_ptr<Segment>& rhs) const noexcept {
        return *lhs < *rhs;
    }
};
//! @endcond

std::list<std::shared_ptr<Segment>> Population::coalesce() const {
    std::set<std::shared_ptr<Segment>, less_shptr_Segment> nodeset;
    std::list<std::shared_ptr<Segment>> nodes;
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            for (bool b: {true, false}) {
                auto it = nodeset.emplace(std::make_shared<Segment>(p.get(), b, true)).first;
                nodes.push_back(*it);
            }
        }
    }
    const size_t sample_size = nodes.size();
    uint_fast32_t num_coalesced = 0u;
    uint_fast32_t num_uncoalesced = 0u;
    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        auto& p = *it;
        if (p->individual_->is_first_gen()) {
            ++num_uncoalesced;
            continue;
        }
        auto res = nodeset.emplace(std::make_shared<Segment>(p->ancestor(), wtl::generate_canonical(*engine_) < 0.5, false));
        auto ancp = const_cast<Segment*>(res.first->get());
        p->set_ancestral_segment(ancp);
        if (res.second) {
            nodes.push_back(*res.first);
        } else {
            ++num_coalesced;
        }
    }
    DCERR("uncoalesced: " << num_uncoalesced << " / " << sample_size << "\n");
    WTL_ASSERT(num_coalesced + num_uncoalesced == sample_size);
    return nodes;
}

std::ostream& Population::write_ms(const std::list<std::shared_ptr<Segment>>& nodes, double lambda, std::ostream& ost) const {
    RunifPoisson runif(lambda);
    std::set<double> collector;
    for (const auto& x: nodes) {
        auto sites = runif(*engine_);
        collector.insert(sites.begin(), sites.end());
        x->set_mutations(std::move(sites));
    }
    std::vector<double> positions(collector.begin(), collector.end());

    ost << "\n//\n"
        << "segsites: " << positions.size() << "\n"
        << "positions: ";
    wtl::join(positions, ost, " ") << "\n";
    for (const auto& x: nodes) {
        if (x->is_sampled_) {
            wtl::join(x->binary_genotype(positions), ost, "") << "\n";
        }
    }
    return ost;
}

void Population::append_demography() {
    std::vector<std::map<uint_fast32_t, size_t>> counter(Individual::num_locations());
    for (const auto& p: males_) {++counter[p->location()][year_ - p->birth_year()];}
    for (const auto& p: females_) {++counter[p->location()][year_ - p->birth_year()];}
    demography_.emplace_hint(demography_.end(), year_, counter);
}

std::ostream& Population::write_demography(std::ostream& ost) const {
    ost << "year\tlocation\tage\tcount\n";
    for (const auto& year_structure: demography_) {
        const auto year = year_structure.first;
        for (size_t loc=0; loc<Individual::num_locations(); ++loc) {
            const auto structure = year_structure.second[loc];
            for (const auto& age_count: structure) {
                ost << year << "\t" << loc << "\t"
                    << age_count.first << "\t"
                    << age_count.second << "\n";
            }
        }
    }
    return ost;
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

} // namespace pbf
