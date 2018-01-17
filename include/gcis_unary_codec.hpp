//
// Created by danielsaad on 10/1/17.
//

#ifndef GC_IS_GCIS_UNARY_CODEC_HPP
#define GC_IS_GCIS_UNARY_CODEC_HPP

#include <cstdint>
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "util.hpp"


class gcis_unary_codec_level{
public:
    gcis_unary_codec_level() = default;
    gcis_unary_codec_level(gcis_unary_codec_level& ) = default;
    gcis_unary_codec_level(gcis_unary_codec_level&&) = default;

    sdsl::int_vector<> rule;
    sdsl::bit_vector rule_delim;
    sdsl::bit_vector::select_1_type rule_delim_sel;
public:
    void expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l);
    void expand_rule(uint64_t rule_num, char* s, uint64_t &l);
};

class gcis_unary_codec {
public:
    //Default constructor
    gcis_unary_codec() = default;
    //Default copy constructor
    gcis_unary_codec(gcis_unary_codec& ) = default;
    //Default move constructor
    gcis_unary_codec(gcis_unary_codec&& ) = default;

    uint_t  string_size;
    uint_t alphabet_size;
    sdsl::bit_vector lcp;
    sdsl::bit_vector rule_delim;
    sdsl::bit_vector::select_1_type rule_sel;
    sdsl::bit_vector::select_1_type lcp_sel;
    sdsl::int_vector<> rule;
    sdsl::int_vector<> tail;
public:
    //! \brief Give the number of bytes of this data structure
    //! \return The total number of bytes spent
    uint64_t size_in_bytes();

    void expand_rule(uint64_t rule,sdsl::int_vector<>& r_string);

    gcis_unary_codec_level decompress();

    void serialize(std::ostream& o);

    void load(std::istream& i);

    uint64_t get_lcp_length(uint64_t i);
    uint64_t get_s_length(uint64_t i);
    uint64_t get_s_start(uint64_t i);

    void extract_lms_substring_suffix(uint64_t i,
                                      uint64_t l,
                                      sdsl::int_vector<>& extracted_text,
                                      uint64_t& k);

private:
    void expand_lcp(uint64_t rule_num,
                    uint64_t lcp_length,
                    sdsl::int_vector<>& r_string,
                    uint64_t rule_start
    );
};



#endif //GC_IS_GCIS_UNARY_CODEC_HPP
