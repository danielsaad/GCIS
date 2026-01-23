#include "gcis_s8b_no_lcp_codec.hpp"

uint64_t gcis_s8b_no_lcp_codec::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += 2 * sizeof(uint_t);
    total_bytes += rule_suffix_length.size_in_bytes();
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    return total_bytes;
}

void gcis_s8b_no_lcp_codec::serialize(std::ostream &o) {
    o.write((char *)&string_size, sizeof(string_size));
    o.write((char *)&alphabet_size, sizeof(alphabet_size));
    rule_suffix_length.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
}
void gcis_s8b_no_lcp_codec::load(std::istream &i) {
    i.read((char *)&string_size, sizeof(string_size));
    i.read((char *)&alphabet_size, sizeof(alphabet_size));
    rule_suffix_length.load(i);
    rule.load(i);
    tail.load(i);
}

gcis_s8b_no_lcp_pointers_codec_level gcis_s8b_no_lcp_codec::decompress() {
    gcis_s8b_no_lcp_pointers_codec_level gd;
    uint64_t number_of_rules = rule_suffix_length.size();
    uint64_t total_rule_suffix_length = 0;
    // Compute the total rule suffix length and the number of rules
    for (uint64_t i = 0; i < rule_suffix_length.size(); i++) {
        total_rule_suffix_length += rule_suffix_length.get_next();
    }

    // Resize data structures
    uint64_t total_length = total_rule_suffix_length;
    gd.rule.width(sdsl::bits::hi(alphabet_size - 1) + 1);
    gd.rule.resize(total_length);
    uint64_t rule_start = 0;
    uint64_t prev_rule_start = 0;
    uint64_t start = 0;
    rule_suffix_length.reset();

    for (uint64_t i = 0; i < number_of_rules; i++) {
        int64_t k;

        uint64_t rule_length = rule_suffix_length.get_next();

        total_length = rule_length;

        uint64_t j = 0;
        while (j < total_length) {
            gd.rule[rule_start + j] = rule[start++];
            j++;
        }

        gd.rule_pos.push_back(rule_start);
        prev_rule_start = rule_start;
        rule_start += total_length;
    }
    gd.rule_pos.push_back(rule_start);
    return gd;
}

void gcis_s8b_no_lcp_codec_level::expand_rule(uint_t rule_num,
                                       sdsl::int_vector<> &r_string,
                                       uint_t &l) {
    uint_t rule_start = rule_delim_sel(rule_num + 1);
    uint_t rule_length = rule_delim_sel(rule_num + 2) - rule_start;
    for (uint64_t i = 0; i < rule_length; i++) {
        r_string[l] = rule[rule_start + i];
        l++;
    }
}

void gcis_s8b_no_lcp_codec_level::expand_rule(uint_t rule_num, char *s, uint_t &l) {
    uint_t rule_start = rule_delim_sel(rule_num + 1);
    uint_t rule_length = rule_delim_sel(rule_num + 2) - rule_start;
    for (uint64_t i = 0; i < rule_length; i++) {
        s[l] = rule[rule_start + i];
        l++;
    }
}

void gcis_s8b_no_lcp_pointers_codec_level::expand_rule(uint_t rule_num,
                                                sdsl::int_vector<> &r_string,
                                                uint_t &l) {
    uint_t rule_start = rule_pos[rule_num];
    uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
    for (uint_t i = 0; i < rule_length; i++) {
        r_string[l] = rule[rule_start + i];
        l++;
    }
}

void gcis_s8b_no_lcp_pointers_codec_level::expand_rule(uint_t rule_num, char *s,
                                                uint_t &l) {
    uint_t rule_start = rule_pos[rule_num];
    uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
    for (uint_t i = 0; i < rule_length; i++) {
        s[l] = rule[rule_start + i];
        l++;
    }
}

void gcis_s8b_no_lcp_pointers_codec_level::expand_rule_bkt(
    uint_t rule_num, sdsl::int_vector<> &r_string, uint_t &l, int_t *bkt) {
    uint_t rule_start = rule_pos[rule_num];
    uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
    for (uint_t i = 0; i < rule_length; i++) {
        r_string[l++] = rule[rule_start + i];
        bkt[rule[rule_start + i]]++;
    }
}

void gcis_s8b_no_lcp_pointers_codec_level::expand_rule_bkt(uint_t rule_num,
                                                    unsigned char *s, uint_t &l,
                                                    int_t *bkt) {
    uint_t rule_start = rule_pos[rule_num];
    uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
    for (uint_t i = 0; i < rule_length; i++) {
        s[l++] = rule[rule_start + i];
        bkt[rule[rule_start + i]]++;
    }
}