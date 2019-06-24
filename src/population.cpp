/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/exception.hpp>
#include <sfmt.hpp>

namespace pbf {

Population::Population(const size_t initial_size, std::random_device::result_type seed)
: engine_(std::make_unique<URBG>(seed)) {
    const size_t half = initial_size / 2UL;
    const size_t rest = initial_size - half;
    males_.reserve(half);
    females_.reserve(rest);
    for (size_t i=0; i<half; ++i) {males_.emplace_back(std::make_shared<Individual>());}
    for (size_t i=0; i<rest; ++i) {females_.emplace_back(std::make_shared<Individual>());}
}

Population::~Population() = default;

void Population::run(const int_fast32_t simulating_duration,
                     const std::vector<size_t>& sample_size_adult,
                     const std::vector<size_t>& sample_size_juvenile,
                     const int_fast32_t recording_duration) {
    auto recording_start = simulating_duration - recording_duration;
    append_demography(0);
    for (year_ = 0; year_ < simulating_duration; ++year_) {
        reproduce();
        survive(0, false);
        survive(1, false);
        survive(2, false);
        survive(3, true);
        if (year_ >= recording_start) {
            sample(sample_size_adult, sample_size_juvenile);
        }
        migrate();
    }
    append_demography(0);
}

void Population::reproduce() {
    double female_biomass = 0.0;
    for (const auto& mother: females_) {
        if (mother->is_in_breeding_place()) {
            female_biomass += mother->weight(year_);
        }
    }
    double male_biomass = 0.0;
    std::vector<std::vector<std::shared_ptr<Individual>>> males_located(Individual::num_breeding_places());
    for (const auto& p: males_) {
        if (p->is_in_breeding_place()) {
            male_biomass += p->weight(year_);
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
    const double biomass = female_biomass + male_biomass;
    const double popsize = biomass / Individual::weight_for_age().back();
    const double density_effect = (1.0 - popsize / Individual::param().CARRYING_CAPACITY);
    const double exp_recruitment = density_effect * Individual::param().RECRUITMENT_COEF * female_biomass;
    std::vector<std::shared_ptr<Individual>> boys;
    std::vector<std::shared_ptr<Individual>> girls;
    boys.reserve(static_cast<size_t>(exp_recruitment * 0.6));
    girls.reserve(static_cast<size_t>(exp_recruitment * 0.6));
    for (const auto& mother: females_) {
        if (!mother->is_in_breeding_place()) continue;
        const auto& potential_fathers = males_located[mother->location()];
        if (potential_fathers.empty()) continue;
        auto& mate_distr = mate_distrs[mother->location()];
        const uint_fast32_t num_juveniles = mother->recruitment(year_, density_effect, *engine_);
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

void Population::survive(const int_fast32_t season, bool shrink) {
    auto has_survived = [season, this](const auto& p) -> bool {
        return p && p->has_survived(year_, season, *engine_);
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
    append_demography(season);
}

void Population::migrate() {
    for (auto& p: males_) {p->migrate(year_, *engine_);}
    for (auto& p: females_) {p->migrate(year_, *engine_);}
}

void Population::sample(const std::vector<size_t>& sample_size_adult,
                        const std::vector<size_t>& sample_size_juvenile) {
    std::vector<std::vector<size_t>> sample_size_female(2u, std::vector<size_t>(Individual::num_locations()));
    std::vector<std::vector<size_t>> sample_size_male = sample_size_female;
    size_t total_sampled_juvenile = 0u;
    size_t total_sampled_adult = 0u;
    size_t total_sampled_female = 0u;
    size_t total_sampled_male = 0u;
    for (size_t i = 0; i < sample_size_juvenile.size(); ++i) {
        const size_t n = sample_size_juvenile[i];
        total_sampled_juvenile += n;
        std::binomial_distribution<size_t> binom(n, 0.5);
        total_sampled_female += sample_size_female[0u][i] = binom(*engine_);
        total_sampled_male += sample_size_male[0u][i] = n - sample_size_female[0u][i];
    }
    for (size_t i = 0; i < sample_size_adult.size(); ++i) {
        const size_t n = sample_size_adult[i];
        total_sampled_adult += n;
        std::binomial_distribution<size_t> binom(n, 0.5);
        total_sampled_female += sample_size_female[1u][i] = binom(*engine_);
        total_sampled_male += sample_size_male[1u][i] = n - sample_size_female[1u][i];
    }
    WTL_ASSERT(total_sampled_adult + total_sampled_juvenile == total_sampled_female + total_sampled_male);
    year_samples_[year_].reserve(total_sampled_adult + total_sampled_juvenile);
    using VecPtrs = std::vector<std::shared_ptr<Individual>>;
    auto impl = [this](VecPtrs* individuals, size_t total,
                       std::vector<std::vector<size_t>>& sample_size) {
        VecPtrs& sampled = year_samples_[year_];
        VecPtrs unsampled;
        unsampled.reserve(std::min(individuals->size() - total, 0ul));
        std::shuffle(individuals->begin(), individuals->end(), *engine_);
        for (const auto p: *individuals) {
            const size_t idx = (p->birth_year() == year_) ? 0u : 1u;
            auto& to_be_sampled = sample_size[idx][p->location()];
            if (to_be_sampled) {
                sampled.emplace_back(p);
                --to_be_sampled;
            } else {
                unsampled.emplace_back(p);
            }
        }
        individuals->swap(unsampled);
    };
    year_samples_[year_].reserve(total_sampled_female + total_sampled_male);
    impl(&females_, total_sampled_female, sample_size_female);
    impl(&males_, total_sampled_male, sample_size_male);
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    if (year_samples_.empty()) return ost;
    wtl::join(Individual::names(), ost, "\t") << "\tcapture_year\n";
    std::map<const Individual*, size_t> ids;
    ids.emplace(nullptr, 0u);
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            ids.emplace(p.get(), ids.size());
        }
    }
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            p->trace_back(ost, &ids, ys.first);
        }
    }
    return ost;
}

std::vector<std::map<int_fast32_t, size_t>> Population::count() const {
    std::vector<std::map<int_fast32_t, size_t>> counter(Individual::num_locations());
    for (const auto& p: males_) {
        if (p) ++counter[p->location()][year_ - p->birth_year()];
    }
    for (const auto& p: females_) {
        if (p) ++counter[p->location()][year_ - p->birth_year()];
    }
    return counter;
}

void Population::append_demography(const int_fast32_t season) {
    demography_.emplace_hint(
      demography_.end(),
      std::pair<int_fast32_t, int_fast32_t>(year_, season),
      count()
    );
}

std::ostream& Population::write_demography(std::ostream& ost) const {
    ost << "year\tseason\tlocation\tage\tcount\n";
    for (const auto& time_structure: demography_) {
        const auto time = time_structure.first;
        for (size_t loc=0; loc<Individual::num_locations(); ++loc) {
            const auto structure = time_structure.second[loc];
            for (const auto& age_count: structure) {
                ost << time.first << "\t" << time.second << "\t"
                    << loc << "\t"
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
