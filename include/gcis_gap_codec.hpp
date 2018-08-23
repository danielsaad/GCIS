

#ifndef GCIS_GAP_CODEC_HPP
#define GCIS_GAP_CODEC_HPP

#include "sdsl/bit_vectors.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/sd_vector.hpp"
#include "sdsl/rrr_vector.hpp"
#include "util.hpp"
#include <cstdint>

class gcis_gap_codec_level {
  public:
    gcis_gap_codec_level() = default;
    gcis_gap_codec_level(gcis_gap_codec_level &) = default;
    gcis_gap_codec_level(gcis_gap_codec_level &&);

    sdsl::int_vector<> rule;
    sdsl::bit_vector rule_delim;
    sdsl::bit_vector::select_1_type rule_delim_sel;

  public:
    void expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string,
                     uint64_t &l);
    void expand_rule(uint64_t rule_num, char *s, uint64_t &l);
};

class gcis_gap_codec {
  public:
    // Default constructor
    gcis_gap_codec() = default;
    // Default copy constructor
    gcis_gap_codec(gcis_gap_codec &) = default;
    // Default move constructor
    gcis_gap_codec(gcis_gap_codec &&) = default;

    uint_t string_size;
    uint_t alphabet_size;
    sdsl::enc_vector<sdsl::coder::elias_delta> lcp;
    sdsl::enc_vector<sdsl::coder::elias_delta> rule_pos;
    sdsl::int_vector<> rule;
    sdsl::int_vector<> tail;
    // Fixed-width integer reduced string
    std::vector<uint64_t> reduced_string_ps;
    // A dac vector storing the fully decoded rule lengths
    sdsl::dac_vector_dp<> fully_decoded_rule_len;
    // A integer storing the fully decoded tail length
    uint64_t fully_decoded_tail_len;

  public:
    /**
     * @brief Tell the number of bytes consumed by the object.
     *
     * @return uint64_t The number of bytes consumed by the object.
     */
    uint64_t size_in_bytes();

    /**
     * @brief Get the first position of a rule suffix in the concatenated
     * data structure.
     *
     * @param i The rule id.
     * @return uint64_t The rule suffix position.
     */
    uint64_t get_rule_pos(uint64_t i);

    /**
     * @brief Decompress all the rules in order to do extraction more quickly.
     *
     * @return gcis_gap_codec_level A data structure containing all the rules
     * expanded in order to make decompression faster.
     */
    gcis_gap_codec_level decompress();

    /**
     * @brief Get the rule's suffix length
     * @param i The rule number
     * @return uint64_t The length of the rule part not encoded by the LCP.
     */
    uint64_t get_rule_length(uint64_t i);

    /**
     * @brief Get the lcp length between a rule and the previous one.
     *
     * @param i The rule number.
     * @return uint64_t The lcp length between rules i and i-1. lcp[0] = 0.
     */
    uint64_t get_lcp(uint64_t i);

    /**
     * @brief Let R be the expansion of a rule. This procedure extracts a rule
     * into a buffer regarding the substring R[l,r].
     *
     * @param rule_num Rule number
     * @param l The beggining of the rule's expansion substring.
     * @param r The end of the rules' expansion substring.
     * @param extracted_text  The buffer where the rule's expansion substring
     * shall be stored.
     * @param k A counter which marks the next available index of the buffer.
     */
    void extract_rule(uint64_t rule_num, int64_t l, int64_t r,
                      sdsl::int_vector<> &extracted_text, uint64_t &k);

    /**
     * @brief Extract a rule into a buffer.
     *
     * @param rule_num The rule number.
     * @param extracted_text  A buffer
     * @param k A counter which marks the next available index in the buffer.
     */
    void extract_rule(uint64_t rule_num, sdsl::int_vector<> &extracted_text,
                      uint64_t &k);

    /**
     * @brief Extracts the rule's suffix part into a buffer.
     *
     * @param rule_num Rule number.
     * @param extracted_text The buffer.
     * @param k A counter which marks the next available position
     * in the buffer.
     */
    void extract_rule_suffix(uint64_t rule_num,
                             sdsl::int_vector<> &extracted_text, uint64_t &k);

    /**
     * @brief Serializes the object into a file.
     *
     * @param o The ostream in which the information shall be stored.
     */
    void serialize(std::ostream &o);

    /**
     * @brief Loads the object from a file.
     *
     * @param i The istream from which the information shall be loaded.
     */
    void load(std::istream &i);

  private:
    /**
     * @brief Extract the LCP part of a rule into a buffer.
     *
     * @param rule_num Rule number.
     * @param extracted_text The buffer.
     * @param k A counter which marks the next available position
     * of the buffer.
     */
    void extract_lcp(uint64_t rule_num, sdsl::int_vector<> &extracted_text,
                     uint64_t &k);

    /**
     * @brief Let S be the expansion of a rule's suffix part.
     * This procedure S[l,r] into a buffer.
     *
     * @param rule_num Rule number.
     * @param l The beginning of the substring S.
     * @param r The end of substring S.
     * @param extracted_text The buffer.
     * @param k A counter which marks the next available position of
     * the buffer.
     */
    void extract_rule_suffix(uint64_t rule_num, int64_t l, int64_t r,
                             sdsl::int_vector<> &extracted_text, uint64_t &k);
};

#endif // GC_IS_GCIS_UNARY_CODEC_HPP
