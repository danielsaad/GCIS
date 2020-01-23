//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_ELIASFANO_INDEX_HPP
#define GC_IS_GCIS_ELIASFANO_INDEX_HPP

#include "grammar_builder.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/sd_vector.hpp"
#include "util.hpp"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <vector>

namespace gcis {

// Forward declaration
class elias_fano_grammar_builder;

/**
 * @brief Gives the grammar info for each level.
 */
class elias_fano_grammar_info {
  public:
    elias_fano_grammar_info() = default;
    // Number of rules in each level
    std::vector<uint_t> m_number_of_rules;
    // Text size in each level
    std::vector<uint_t> m_text_size;
    // alphabet size in each level
    std::vector<uint_t> m_alphabet_size;
    // Number of total levels in the grammar
    uint_t m_level_n = 0;
    // Grammar size (sum of right-hand)
    uint_t m_grammar_size = 0;
    // Length of the right-hand side concatenation from the rules in the first
    // level (the ones that derives terminals)
    uint_t m_first_level_expansion_len = 0;

    void serialize(ofstream &o) {
        util::serialize(m_number_of_rules, o);
        util::serialize(m_text_size, o);
        util::serialize(m_alphabet_size, o);
        util::serialize(m_level_n, o);
        util::serialize(m_grammar_size, o);
        util::serialize(m_first_level_expansion_len, o);
    }
    void load(ifstream &in) {
        util::load(m_number_of_rules, in);
        util::load(m_text_size, in);
        util::load(m_alphabet_size, in);
        util::load(m_level_n, in);
        util::load(m_grammar_size, in);
        util::load(m_first_level_expansion_len, in);
    }
    void print() {
        using std::cout;
        using std::endl;
        for (uint_t i = 0; i < m_level_n; i++) {
            cout << "Level = " << i << endl;
            cout << "Text size = " << m_text_size[i] << endl;
            cout << "Alphabet size = " << m_alphabet_size[i] << endl;
            cout << "Number of rules = " << m_number_of_rules[i] << endl;
            cout << "Grammar size = " << m_grammar_size << endl;
            cout << "First level expansion = " << m_first_level_expansion_len
                 << endl;
        }
    }
};

class elias_fano_grammar {
  public:
    typedef elias_fano_grammar_builder grammar_builder_t;

    elias_fano_grammar() = default;
    elias_fano_grammar(const elias_fano_grammar &rhs) = default;
    elias_fano_grammar(elias_fano_grammar &&rhs) = default;
    elias_fano_grammar &operator=(elias_fano_grammar &rhs) {
        using std::swap;
        elias_fano_grammar tmp(rhs);
        swap(*this, tmp);
        return *this;
    }
    elias_fano_grammar &operator=(elias_fano_grammar &&rhs) {
        if (this != &rhs) {
            using std::swap;
            swap(*this, rhs);
        }
        return *this;
    }

    void swap(elias_fano_grammar &lhs, elias_fano_grammar &rhs) {
        using std::swap;
        swap(lhs.rules_suffix, rhs.rules_suffix);
        swap(lhs.rules_lcp, rhs.rules_lcp);
        swap(lhs.rules_delim, rhs.rules_delim);
        swap(lhs.m_xs, rhs.m_xs);
        swap(lhs.rules_lcp_select, rhs.rules_lcp_select);
        swap(lhs.rules_delim_select, rhs.rules_delim_select);
        rules_lcp_select.set_vector(&rules_lcp);
        rules_delim_select.set_vector(&rules_delim);
    }

    char *decode() {
        m_info.print();
        char *str = new char[m_info.m_text_size[1]];
        // Decompress all the rules
        uint_t total_rules = std::accumulate(m_info.m_number_of_rules.begin(),
                                             m_info.m_number_of_rules.end(), 0);
        // fix for first level, if we have rules Xc->c this is not necessary
        auto width = max(sdsl::bits::hi(total_rules) + 1, sdsl::bits::hi(256));
        sdsl::int_vector<> rules_derivation(m_info.m_grammar_size, 0, width);
        width = sdsl::bits::hi(m_info.m_grammar_size) + 1;
        sdsl::int_vector<> rules_pos(total_rules + 1, 0, width);

        // Decompress sequentially each rule

        // Special case to avoid ifs
        uint_t idx = 0;
        rules_derivation[idx] = 0;
        rules_pos[idx++] = 0;
        uint_t rule_concat_idx = 1;

        uint_t prev_lcp_pos = 0;
        uint_t prev_rule_pos = 1;
        uint_t prev_rule_len = 0;

        for (uint_t i = 1; i < total_rules; i++) {
            rules_pos[i] = idx;
            auto cur_lcp_pos = rules_lcp_select(i + 1);
            uint_t lcp_len = cur_lcp_pos - prev_lcp_pos - 1;
            prev_lcp_pos = cur_lcp_pos;

            auto cur_rule_pos = rules_delim_select(i + 1);
            uint_t suffix_len = cur_rule_pos - prev_rule_pos - 1;
            prev_rule_pos = cur_rule_pos;

            auto cur_rule_len = lcp_len + suffix_len;
            copy_lcp(rules_derivation, lcp_len, prev_rule_len, idx);
            copy_suffix(rules_derivation, suffix_len, rule_concat_idx, idx);
            prev_rule_len = cur_rule_len;
        }
        rules_pos[total_rules] = m_info.m_grammar_size;

        idx = 0;
        // fix for first level, if we have rules Xc->c this is not necessary
        width = max(sdsl::bits::hi(total_rules) + 1, sdsl::bits::hi(256));
        int_t stack_idx = 0;
        sdsl::int_vector<> rule_stack(m_info.m_grammar_size, 0, width);
        rule_stack[stack_idx++] = m_xs;
        while (stack_idx > 0) {
            uint_t rule_idx = rule_stack[--stack_idx];
            uint_t pos = rules_pos[rule_idx];
            uint_t len = rules_pos[rule_idx + 1] - pos;
            // cout << "Decompressing rule " << rule_idx << endl;
            if (rule_idx < 256) {
                // Will generate terminal symbols
                // cout << "Generating terminal symbols from rule " << rule_idx
                //  << endl;
                str[idx++] = rule_idx;

                // cout << endl;
            } else {
                // Will expand into non-terminal symbols
                for (int_t i = len + pos - 1; i >= pos; i--) {
                    // cout << "Inserting rule " << rules_derivation[i]
                    //      << " with level " << rule_level + 1
                    //      << " into the stack" << endl;
                    rule_stack[stack_idx++] = rules_derivation[i];
                }
            }
        }

        return str;
    }

    std::string extract_rule(int_t rule_id) const {};
    int compare_suffix(int_t rule_id, const std::string &pattern,
                       size_t len) const {};
    int compare_preffix(int_t rule_id, const std::string &pattern,
                        size_t len) const {};

    void serialize(std::ofstream &o) {
        rules_suffix.serialize(o);
        rules_lcp.serialize(o);
        rules_delim.serialize(o);
        rules_lcp_select.serialize(o);
        rules_delim_select.serialize(o);
        util::serialize(m_xs, o);

        m_info.serialize(o);

        // cout << "Rules concat" << endl;
        // for (uint_t i = 0; i < rules_suffix.size(); i++) {
        //     cout << rules_suffix[i] << ' ';
        // }
        // cout << endl << "LCP " << endl;
        // for (uint_t i = 0; i < rules_lcp.size(); i++) {
        //     cout << rules_lcp[i] << ' ';
        // }
        // cout << endl << "Rules delimiters " << endl;
        // for (uint_t i = 0; i < rules_delim.size(); i++) {
        //     cout << rules_delim[i] << ' ';
        // }
        // cout << endl;
    }

    void load(std::ifstream &in) {
        rules_suffix.load(in);
        rules_lcp.load(in);
        rules_delim.load(in);
        rules_lcp_select.load(in);
        rules_delim_select.load(in);
        util::load(m_xs, in);

        rules_lcp_select.set_vector(&rules_lcp);
        rules_delim_select.set_vector(&rules_delim);

        m_info.load(in);

        // cout << "Rules concat" << endl;
        // for (uint_t i = 0; i < rules_suffix.size(); i++) {
        //     cout << rules_suffix[i] << ' ';
        // }
        // cout << endl << "LCP " << endl;
        // for (uint_t i = 0; i < rules_lcp.size(); i++) {
        //     cout << rules_lcp[i] << ' ';
        // }
        // cout << endl << "Rules delimiters " << endl;
        // for (uint_t i = 0; i < rules_delim.size(); i++) {
        //     cout << rules_delim[i] << ' ';
        // }
        // cout << endl;
    }

  public:
    // stores the rules suffix in a single concatenated array
    sdsl::int_vector<0> rules_suffix;
    // stores the LCP encoding between rules using Elias-Fano encoding
    sdsl::sd_vector<> rules_lcp;
    // stores the rules delimiters
    sdsl::sd_vector<> rules_delim;
    // stores the rule index which expansion generates the text
    uint_t m_xs;
    // Select support for lcp
    sdsl::sd_vector<>::select_1_type rules_lcp_select;
    // Select support for rules_delim
    sdsl::sd_vector<>::select_1_type rules_delim_select;
    // Grammar info
    elias_fano_grammar_info m_info;

    void copy_lcp(sdsl::int_vector<> &rules_derivation,
                  const uint_t cur_lcp_len, const uint_t prev_rule_len,
                  uint_t &idx) {
        for (uint_t i = 0; i < cur_lcp_len; i++) {
            rules_derivation[idx + i] =
                rules_derivation[idx + i - prev_rule_len];
            // cout << rules_derivation[idx + i] << " ";
            // cout.flush();
        }
        idx += cur_lcp_len;
    }

    void copy_suffix(sdsl::int_vector<> &rules_derivation,
                     const uint_t suffix_len, uint_t &rule_concat_idx,
                     uint_t &idx) {
        for (uint_t i = 0; i < suffix_len; i++) {
            rules_derivation[idx + i] = rules_suffix[rule_concat_idx + i];
            // cout << rules_derivation[idx + i] << " ";
        }
        // cout << endl;
        idx += suffix_len;
        rule_concat_idx += suffix_len;
    }
    friend class elias_fano_grammar_builder;
};

/**
 * @brief Holds information necessary to process
 * the elias fano grammar.
 */
class elias_fano_grammar_builder {

  public:
    elias_fano_grammar_builder() = default;

    void set_gcis_structures(int_t *s, uint_t *SA, int_t n, int_t K, int cs,
                             int level, unsigned char *t) {
        m_text = s;
        m_SA = SA;
        m_text_size = n;
        m_alphabet_size = K;
        m_cs = cs;
        m_level = level;
        m_type = t;
    }

    void update_grammar(std::vector<int_t> &lcp) {
        // account for the tail
        auto &info = m_grammar_info;
        update_lcp(lcp);
        update_rules_suffix(lcp);
        update_tail();
        info.m_number_of_rules.push_back(lcp.size() + 1);
        info.m_text_size.push_back(m_text_size);
        info.m_alphabet_size.push_back(m_alphabet_size);
        info.m_level_n++;
        m_last_tail_idx += lcp.size() + 1;
        // info.print();
    }

    elias_fano_grammar &get_grammar() { return m_g; }

    /**
     * @brief Creates a new rule
     * m_xs -> last_tail + reduced_string
     * Hence, we need to store just the ID m_xs to
     */
    void compute_reduced_string() {

        // Set LCP to 0
        lcp_unary.resize(lcp_unary.size() + 1);
        lcp_unary[lcp_unary.size() - 1] = 1;

        // We will insert a right hand side of size
        // |reduced_string| + 1 (we have to account for the tail rule)
        uint_t base = rules_delim_unary.size();
        rules_delim_unary.resize(base + m_text_size + 2);
        for (uint_t i = base; i < rules_delim_unary.size() - 1; i++) {
            rules_delim_unary[i] = 0;
        }
        rules_delim_unary[rules_delim_unary.size() - 1] = 1;

        // Copy the text to the right-hand side of m_xs
        base = rules_concat.size();
        rules_concat.resize(base + m_text_size + 1);
        rules_concat[base] = m_last_tail_idx;
        uint_t offset =
            std::accumulate(m_grammar_info.m_number_of_rules.begin(),
                            m_grammar_info.m_number_of_rules.end() - 1, 0);
        // cout << "Last tail idx = " << m_last_tail_idx << endl;
        for (uint_t j = 0; j < m_text_size; j++) {
            rules_concat[j + base + 1] = m_text[j] + offset;
            // cout << rules_concat[j + base + 1] << ' ';
        }
        // cout << endl;
        // Update grammar info
        m_grammar_info.m_alphabet_size.push_back(m_text_size);
        m_grammar_info.m_number_of_rules.push_back(1);
        m_grammar_info.m_text_size.push_back(m_text_size);
        m_grammar_info.m_level_n++;
        m_grammar_info.m_grammar_size += m_text_size + 1;
        // m_grammar_info.print();
    }

    /**
     * @brief Preprocessing adds rules of type Xc->c in the grammar.
     * Now, every non-terminal different from Xc has index >= 256. So we can
     * differentiate a non-terminal from a Xc rule by looking to its index.
     */
    void pre_process() {
        lcp_unary.resize(256);
        rules_delim_unary.resize(512);
        rules_concat.resize(256);
        for (uint_t i = 0; i < 256; i++) {
            lcp_unary[i] = 1;
            rules_concat[i] = i;
            rules_delim_unary[i << 1] = 0;
            rules_delim_unary[(i << 1) + 1] = 1;
        }
        m_grammar_info.m_number_of_rules.push_back(256);
        m_grammar_info.m_text_size.push_back(0);
        m_grammar_info.m_alphabet_size.push_back(256);
        m_grammar_info.m_grammar_size += 256;
        m_last_tail_idx = 255;
    }
    void post_process() {
        // Copy the data to the grammar and perform further compressions
        sdsl::util::bit_compress(rules_concat);
        m_g.rules_delim = sdsl::sd_vector<>(rules_delim_unary);
        m_g.rules_lcp = sdsl::sd_vector<>(lcp_unary);
        m_g.rules_suffix = std::move(rules_concat);
        m_g.m_info = m_grammar_info;
        sdsl::util::bit_compress(m_g.rules_suffix);
        sdsl::util::init_support(m_g.rules_lcp_select, &(m_g.rules_lcp));
        sdsl::util::init_support(m_g.rules_delim_select, &(m_g.rules_delim));
        // set xs label
        m_g.m_xs = std::accumulate(m_grammar_info.m_number_of_rules.begin(),
                                   m_grammar_info.m_number_of_rules.end(), 0) -
                   1;
    }

    ~elias_fano_grammar_builder() = default;

  private:
    elias_fano_grammar_info m_grammar_info;
    elias_fano_grammar m_g;
    sdsl::bit_vector lcp_unary;
    sdsl::bit_vector rules_delim_unary;
    sdsl::int_vector<> rules_concat;
    int_t m_last_tail_idx = -1;
    int_t *m_text;
    uint_t *m_SA;
    int_t m_text_size;
    int_t m_alphabet_size;
    int m_cs;
    int m_level;
    unsigned char *m_type;

    void update_lcp(const std::vector<int_t> lcp) {
        /**
         * Updates the LCP in unary form.
         * First it sums all the LCP to resize the bitvector and
         * then it represents each LCP value 'x' as a sequence of
         * x zeroes followed by a one **/
        uint_t lcp_total_size = std::accumulate(lcp.begin(), lcp.end(), 0);
        uint_t lcp_prev_sz = lcp_unary.size();
        lcp_unary.resize(lcp_prev_sz + lcp_total_size + lcp.size());
        for (size_t i = 0; i < lcp.size(); i++) {
            set_unary(lcp_unary, lcp_prev_sz, lcp[i]);
        }
    }

    void update_rules_suffix(const std::vector<int_t> &lcp) {
        /**
         * Stores the rules sufixes in a concatenated int_vector.
         * First it computer the total length of the concatenated suffixes.
         */
        int_t number_of_rules = lcp.size();
        int_t rules_suffix_total_size = 0;
        for (int_t i = 0; i < number_of_rules; i++) {
            int_t pos = m_SA[i];
            rules_suffix_total_size += get_lms_substring_len(pos) - lcp[i];
        }
        /**
         * Now it resizes the int_vector and copy the concatenated suffixes
         * considering the LCP between rules which cannot be copied.
         */
        uint_t rules_cat_idx = rules_concat.size();
        rules_concat.resize(rules_cat_idx + rules_suffix_total_size);
        uint_t rules_delim_idx = rules_delim_unary.size();
        rules_delim_unary.resize(rules_delim_unary.size() + number_of_rules +
                                 rules_suffix_total_size);
        uint_t offset =
            m_level == 0
                ? 0
                : std::accumulate(m_grammar_info.m_number_of_rules.begin(),
                                  m_grammar_info.m_number_of_rules.end() - 1,
                                  0);
        // cout << "Computing rules suffixes. Offset = " << offset << endl;
        for (int_t i = 0; i < number_of_rules; i++) {
            int_t pos = m_SA[i];
            int_t lcp_len = lcp[i];
            int_t rule_len = get_lms_substring_len(pos);
            // If it is the first level
            if (m_level == 0) {
                m_grammar_info.m_first_level_expansion_len += rule_len;
            }
            m_grammar_info.m_grammar_size += rule_len;
            copy_lms_substring(pos, rules_cat_idx, rule_len, lcp_len, offset);
            set_unary(rules_delim_unary, rules_delim_idx, rule_len - lcp_len);
        }
    }

    void update_tail() {
        // cout << "Update tail m_level = " << m_level << endl;
        uint_t tail_len = 0;
        while (!isLMS(m_type, tail_len)) {
            tail_len++;
        }
        // Resize LCP and set tail LCP to 0
        lcp_unary.resize(lcp_unary.size() + 1);
        lcp_unary[lcp_unary.size() - 1] = 1;

        if (m_level == 0) {
            // Resize rule delimiters bitvector and set the last entry to
            // tail_len
            uint_t base = rules_delim_unary.size();
            rules_delim_unary.resize(base + tail_len + 1);
            for (uint_t i = base; i < rules_delim_unary.size(); i++) {
                rules_delim_unary[i] = 0;
            }
            rules_delim_unary[rules_delim_unary.size() - 1] = 1;

            // Resize rule suffix concatenated array
            base = rules_concat.size();
            rules_concat.resize(base + tail_len);
            // cout << "Tail = ";
            for (uint_t i = 0; i < tail_len; i++) {
                rules_concat[base + i] = chr(m_text, m_cs, i);
                // cout << rules_concat[base + i];
            }
            // cout << endl;
            m_grammar_info.m_first_level_expansion_len += tail_len;
        } else {
            // Resize rule delimiters bitvector and set the last entry to
            // tail_len
            tail_len++;
            uint_t base = rules_delim_unary.size();
            rules_delim_unary.resize(base + tail_len + 1);
            for (uint_t i = base; i < rules_delim_unary.size(); i++) {
                rules_delim_unary[i] = 0;
            }
            rules_delim_unary[rules_delim_unary.size() - 1] = 1;

            // Resize rule suffix concatenated array
            base = rules_concat.size();
            rules_concat.resize(base + tail_len);
            rules_concat[base] = m_last_tail_idx;
            uint_t offset =
                std::accumulate(m_grammar_info.m_number_of_rules.begin(),
                                m_grammar_info.m_number_of_rules.end() - 1, 0);
            // cout << "offset = " << offset << endl;
            // cout << "Tail = " << rules_concat[base];

            for (uint_t i = 1; i < tail_len; i++) {
                rules_concat[base + i] = chr(m_text, m_cs, i - 1) + offset;
                // cout << rules_concat[base + i];
            }
            // cout << endl;
        }
        m_grammar_info.m_grammar_size += tail_len;
    }

    int_t get_lms_substring_len(int_t pos, int_t lcp_len = 0) {

        size_t len = 1;
        if (pos != m_text_size - 1)
            while (!isLMS(m_type, pos + len)) {
                len++;
            }
        return len - lcp_len;
    }

    void copy_lms_substring(int_t pos, uint_t &idx, int_t rule_len,
                            int_t lcp_len = 0, uint_t offset = 0) {

        // cout << "Offset = " << offset << endl;
        // cout << "Rule suffix: ";
        for (int_t j = 0;
             j < rule_len - lcp_len && j + pos + lcp_len < m_text_size; j++) {
            rules_concat[idx++] =
                (uint_t)chr(m_text, m_cs, j + pos + lcp_len) + offset;
            // cout << (uint_t)chr(m_text, m_cs, j + pos + lcp_len) + offset <<
            // " ";
        }
        // cout << endl;
    }

    void set_unary(sdsl::bit_vector &b, uint_t &idx, uint_t v) {
        for (uint_t i = idx; i < v + idx; i++) {
            b[i] = 0;
        }
        b[v + idx] = 1;
        idx += v + 1;
    }
};

} // end of namespace gcis

#endif // GC_IS_GCIS_ELIASFANO_INDEX_HPP
