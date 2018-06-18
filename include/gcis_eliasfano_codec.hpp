#ifndef GC_IS_GCIS_ELIASFANO_CODEC_HPP
#define GC_IS_GCIS_ELIASFANO_CODEC_HPP


#include <cstdint>
#include "eliasfano.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/dac_vector.hpp"
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

    //expand rules while counting symbol's frequency
    void expand_rule_bkt(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l, int_t *bkt);
    void expand_rule_bkt(uint64_t rule_num, char* s, uint64_t &l, int_t *bkt);
};


class gcis_eliasfano_codec{
public:
    gcis_eliasfano_codec() = default;
    gcis_eliasfano_codec(gcis_eliasfano_codec& rhs) = default;
    gcis_eliasfano_codec(gcis_eliasfano_codec&& rhs) = default;

    // string size
    uint_t string_size;
     // alphabet size
    uint_t alphabet_size;
    // elias-fano coded lcp information
    eliasfano_codec lcp;
    // fixed-width of rule suffixes array
    sdsl::int_vector<> rule;
    // fixed-width of rule suffixes length
    eliasfano_codec rule_suffix_length;
    // vector of tails
    sdsl::int_vector<> tail;
    // Fixed-width integer reduced string
    std::vector<uint64_t> reduced_string_ps;
    // A dac vector storing the fully decoded rule lengths
    sdsl::dac_vector_dp<> fully_decoded_rule_len;
    // A integer storing the fully decoded tail length
    uint64_t fully_decoded_tail_len;

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
     void extract_lcp(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);
     void extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);
     void extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k);

    void extract_lcp(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k);
    void extract_rule_suffix(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k);
    void extract_rule(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k);



    void serialize(std::ostream& o);
    void load(std::istream& i);
    gcis_eliasfano_codec_level decompress();
};

#endif //GC_IS_GCIS_ELIASFANO_CODEC_HPP
