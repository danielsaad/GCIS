#ifndef GCIS_ELIASFANO_CODEC_NO_LCP_HPP
#define GCIS_ELIASFANO_CODEC_NO_LCP_HPP

#include <cstdint>
#include "eliasfano.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/dac_vector.hpp"
#include "util.hpp"

class gcis_eliasfano_codec_no_lcp_level{
public:
    gcis_eliasfano_codec_no_lcp_level() = default;
    gcis_eliasfano_codec_no_lcp_level(gcis_eliasfano_codec_no_lcp_level& ) = default;
    gcis_eliasfano_codec_no_lcp_level(gcis_eliasfano_codec_no_lcp_level&&) = default;

    sdsl::int_vector<> rule;
    eliasfano_codec rule_delim;

public:
    void expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l);
    void expand_rule(uint64_t rule_num, char* s, uint64_t &l);
};


class gcis_eliasfano_codec_no_lcp{
public:
    gcis_eliasfano_codec_no_lcp() = default;
    gcis_eliasfano_codec_no_lcp(gcis_eliasfano_codec_no_lcp& rhs) = default;
    gcis_eliasfano_codec_no_lcp(gcis_eliasfano_codec_no_lcp&& rhs) = default;

    // string size
    uint_t string_size;
    // alphabet size
    uint_t alphabet_size;
    // fixed-width of rule suffixes array
    sdsl::int_vector<> rule;
    // fixed-width of rule suffixes length
    eliasfano_codec rule_suffix_length;
    // vector of tails
    sdsl::int_vector<> tail;
    // A dac vector storing the fully decoded rule lengths
    sdsl::dac_vector_dp<> fully_decoded_rule_len;
    // A integer storing the fully decoded tail length
    uint64_t fully_decoded_tail_len;

public:
    uint64_t size_in_bytes();
    uint64_t get_rule_pos(uint64_t i);
    gcis_eliasfano_codec_no_lcp_level decompress();
    uint64_t get_rule_length(uint64_t i);

    void extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);
    void extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);

    void extract_rule_suffix(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k);
    void extract_rule(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k);



    void serialize(std::ostream& o);
    void load(std::istream& i);
};


#endif //GCIS_ELIASFANO_CODEC_NO_LCP_HPP
