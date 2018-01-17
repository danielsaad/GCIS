#ifndef GC_IS_GCIS_ELIASFANO_CODEC_HPP
#define GC_IS_GCIS_ELIASFANO_CODEC_HPP


#include <cstdint>
#include "eliasfano.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "util.hpp"

class gcis_eliasfano_codec_level{
public:
    gcis_eliasfano_codec_level() = default;
    gcis_eliasfano_codec_level(gcis_eliasfano_codec_level& ) = default;
    gcis_eliasfano_codec_level(gcis_eliasfano_codec_level&&) = default;

    sdsl::int_vector<> rule;
    eliasfano_codec rule_delim;

public:
    void expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l);
    void expand_rule(uint64_t rule_num, char* s, uint64_t &l);
};

class gcis_eliasfano_codec{
public:
    gcis_eliasfano_codec() = default;
    gcis_eliasfano_codec(gcis_eliasfano_codec& rhs) = default;
    gcis_eliasfano_codec(gcis_eliasfano_codec&& rhs) = default;

    uint_t string_size;
    uint_t alphabet_size;
    eliasfano_codec lcp;
    sdsl::int_vector<> rule;
    eliasfano_codec rule_suffix_length;
    sdsl::int_vector<> tail;
    // sdsl::sd_vector<> lms_bv;
    // sdsl::sd_vector<>::rank_1_type lms_rnk1;
    // sdsl::sd_vector<>::select_1_type lms_sel1;
public:
    uint64_t get_lcp(uint64_t i);

    uint64_t get_rule_pos(uint64_t i);

    uint64_t get_rule_length(uint64_t i);

    uint64_t size_in_bytes();

    // void extract_rules(uint64_t l,
    //                    uint64_t r,
    //                    sdsl::int_vector<>& extracted_text,
    //                    sdsl::int_vector<>& tmp_text);
    //
    // void extract_lcp(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);
    // void extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);
    // void extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);


    void serialize(std::ostream& o);
    void load(std::istream& i);
    gcis_eliasfano_codec_level decompress();
};

#endif //GC_IS_GCIS_ELIASFANO_CODEC_HPP
