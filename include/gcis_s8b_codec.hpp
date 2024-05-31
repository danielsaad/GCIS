//
// Created by danielsaad on 10/1/17.
//

#ifndef GC_IS_GCIS_S8_CODEC_CPP_H
#define GC_IS_GCIS_S8_CODEC_CPP_H

#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "simple8b.hpp"
#include "util.hpp"
#include <cstdint>

class gcis_s8b_codec_level {
  public:
    gcis_s8b_codec_level() = default;
    gcis_s8b_codec_level(gcis_s8b_codec_level &) = default;
    gcis_s8b_codec_level(gcis_s8b_codec_level &&) = default;

    sdsl::int_vector<> rule;
    sdsl::bit_vector rule_delim;
    sdsl::bit_vector::select_1_type rule_delim_sel;

  public:
    void expand_rule(uint_t rule_num, sdsl::int_vector<> &r_string, uint_t &l);
    void expand_rule(uint_t rule_num, char *s, uint_t &l);
};

/**
 * @brief Same as gcis_s8b_codec_level, but stores rules
 * position in a plain array of pointers.
 */
class gcis_s8b_pointers_codec_level {
  public:
    // concatenaded rules
    sdsl::int_vector<> rule;
    // begin of each rule
    vector<uint_t> rule_pos;
    std::tuple<uint_t, uint_t> get_rule_info(uint_t rule_num);

  public:
    void expand_rule(uint_t rule_num, sdsl::int_vector<> &r_string, uint_t &l);
    void expand_rule(uint_t rule_num, char *s, uint_t &l);

    void expand_rule_bkt(uint_t rule_num, sdsl::int_vector<> &r_string,
                         uint_t &l, int_t *bkt);
    void expand_rule_bkt(uint_t rule_num, unsigned char *s, uint_t &l,
                         int_t *bkt);
};

class gcis_s8b_codec {
  public:
    uint_t string_size;
    uint_t alphabet_size;
    simple8b_codec lcp;
    simple8b_codec rule_suffix_length;
    sdsl::int_vector<> rule;
    sdsl::int_vector<> tail;

  public:
    uint64_t size_in_bytes();
    void serialize(std::ostream &o);
    void load(std::istream &i);
    // gcis_s8b_codec_level decompress();
    gcis_s8b_pointers_codec_level decompress();
};

#endif // GC_IS_GCIS_S8_CODEC_CPP_H
