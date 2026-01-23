#ifndef GCIS_STATISTICS_HPP
#define GCIS_STATISTICS_HPP
#include "util.hpp"
#include <cstdint>
#include <iostream>
#include <vector>

using std::vector;
class gcis_level_statistics {
  public:
    void serialize(std::ofstream &o) const {
        gcis::util::serialize(alphabet_size, o);
        gcis::util::serialize(string_size, o);
        gcis::util::serialize(dictionary_size, o);
        gcis::util::serialize(rules_n, o);
        gcis::util::serialize(discarded_rules_n, o);
        gcis::util::serialize(total_discarded_rules_len, o);
        gcis::util::serialize(total_rule_len, o);
        gcis::util::serialize(total_lcp_len, o);
        gcis::util::serialize(total_rule_suffix_len, o);
        gcis::util::serialize(total_rl_potential, o);
        gcis::util::serialize(rule_suffix_width, o);
        gcis::util::serialize(tail_len, o);
        gcis::util::serialize(tail_width, o);
        gcis::util::serialize(avg_discarded_rules_len, o);
        gcis::util::serialize(avg_rule_len, o);
        gcis::util::serialize(avg_lcp, o);
        gcis::util::serialize(avg_rule_suffix_len, o);
        gcis::util::serialize(size_in_bytes, o);
    }
    void load(std::ifstream &i) {
        gcis::util::load(alphabet_size, i);
        gcis::util::load(string_size, i);
        gcis::util::load(dictionary_size, i);
        gcis::util::load(rules_n, i);
        gcis::util::load(discarded_rules_n, i);
        gcis::util::load(total_discarded_rules_len, i);
        gcis::util::load(total_rule_len, i);
        gcis::util::load(total_lcp_len, i);
        gcis::util::load(total_rule_suffix_len, i);
        gcis::util::load(total_rl_potential, i);
        gcis::util::load(rule_suffix_width, i);
        gcis::util::load(tail_len, i);
        gcis::util::load(tail_width, i);
        gcis::util::load(avg_discarded_rules_len, i);
        gcis::util::load(avg_rule_len, i);
        gcis::util::load(avg_lcp, i);
        gcis::util::load(avg_rule_suffix_len, i);
        gcis::util::load(size_in_bytes, i);
    }

  public:
    uint64_t alphabet_size = 0;
    uint64_t string_size = 0;
    uint64_t dictionary_size = 0;
    uint64_t rules_n = 0;
    uint64_t discarded_rules_n = 0;
    uint64_t total_discarded_rules_len = 0;
    uint64_t total_rule_len = 0;
    uint64_t total_lcp_len = 0;
    uint64_t total_rule_suffix_len = 0;
    uint64_t total_rl_potential = 0;
    uint64_t rule_suffix_width = 0;
    uint64_t tail_len = 0;
    uint64_t tail_width = 0;
    double avg_discarded_rules_len = 0;
    double avg_rule_len = 0;
    double avg_lcp = 0;
    double avg_rule_suffix_len = 0;
    uint64_t size_in_bytes = 0;
};

export class gcis_statistics {
  public:
    void serialize(std::ofstream &o) const {
        gcis::util::serialize(level_n, o);
        gcis::util::serialize(premature_stop, o);
        size_t lvl_stats_sz = lvl_stats.size();
        gcis::util::serialize(lvl_stats_sz, o);
        for (const auto &stats : lvl_stats) {
            stats.serialize(o);
        }
    }
    void load(std::ifstream &i) {
        gcis::util::load(level_n, i);
        gcis::util::load(premature_stop, i);
        size_t lvl_stats_sz = 0;
        gcis::util::load(lvl_stats_sz, i);
        lvl_stats.resize(lvl_stats_sz);
        for (auto &stats : lvl_stats) {
            stats.load(i);
        }
    }

  public:
    int level_n = 0;
    bool premature_stop = false;
    vector<gcis_level_statistics> lvl_stats;
};

#endif