#include "gcis_gap_codec.hpp"

gcis_gap_codec_level::gcis_gap_codec_level(gcis_gap_codec_level &&rhs) {
    rule = std::move(rhs.rule);
    rule_delim = std::move(rhs.rule_delim);
    rule_delim_sel = std::move(rhs.rule_delim_sel);
    rule_delim_sel.set_vector(&rule_delim);
}

uint64_t gcis_gap_codec::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += sizeof(string_size);
    total_bytes += sizeof(alphabet_size);
    total_bytes += sdsl::size_in_bytes(lcp);
    total_bytes += sdsl::size_in_bytes(rule_pos);
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    total_bytes += sdsl::size_in_bytes(fully_decoded_rule_len);
    total_bytes += sdsl::size_in_bytes(reduced_string_ps);
    total_bytes += sizeof(fully_decoded_tail_len);
    cout << "Size in bytes reduced_string_pc = "
         << sdsl::size_in_bytes(reduced_string_ps) << endl;
    cout << "Size in bytes lcp = " << sdsl::size_in_bytes(lcp) << endl;
    cout << "Size in bytes rule pos = " << sdsl::size_in_bytes(rule_pos)
         << endl;
    return total_bytes;
}

gcis_gap_codec_level gcis_gap_codec::decompress() {
    gcis_gap_codec_level gd;
    vector<uint64_t> lcp_values(lcp.get_sample_dens());
    vector<uint64_t> rule_pos_values(rule_pos.get_sample_dens());

    uint64_t lcp_idx = 0;
    uint64_t rule_pos_idx = 0;

    uint64_t lcp_bucket_idx = 0;
    uint64_t rule_pos_bucket_idx = 0;

    uint64_t last_lcp = 0;
    uint64_t last_rule_pos = 0;
    uint64_t total_lcp_length = 0;
    uint64_t total_rule_suffix_length = 0;

    uint64_t number_of_rules = lcp.size();

    total_lcp_length = lcp[lcp.size() - 1];
    total_rule_suffix_length = rule_pos[rule_pos.size() - 1];

    /* Once the LCP and rule suffix total lengths are computed, we resize
     * the data structures */

    uint64_t total_length = total_lcp_length + total_rule_suffix_length;
    sdsl::bit_vector rule_delim(total_length + 1);
    gd.rule.resize(total_length);
    gd.rule.width(sdsl::bits::hi(alphabet_size - 1) + 1);
    sdsl::util::set_to_value(rule_delim, 0);

    lcp_bucket_idx = lcp_idx = rule_pos_bucket_idx = rule_pos_idx = last_lcp =
        last_rule_pos = 0;

    uint64_t rule_start = 0;
    uint64_t prev_rule_start = 0;
    uint64_t start = 0;

    lcp.get_inter_sampled_values(0, lcp_values.data());
    uint64_t lcp_sample_value = lcp.sample(0);
    for (size_t i = 0; i < lcp_values.size(); i++) {
        lcp_values[i] += lcp_sample_value;
    }

    rule_pos.get_inter_sampled_values(0, rule_pos_values.data());
    uint64_t rule_pos_sample_value = rule_pos.sample(0);
    for (size_t i = 0; i < rule_pos_values.size(); i++) {
        rule_pos_values[i] += rule_pos_sample_value;
    }

    for (uint64_t i = 0; i < number_of_rules; i++) {
        /* Load more values from lcp and rule_pos */
        if (lcp_idx == lcp.get_sample_dens()) {
            lcp_bucket_idx++;
            lcp.get_inter_sampled_values(lcp_bucket_idx, lcp_values.data());
            lcp_sample_value = lcp.sample(lcp_bucket_idx);
            lcp_idx = 0;
            for (size_t i = 0; i < lcp_values.size(); i++) {
                lcp_values[i] += lcp_sample_value;
            }
        }
        if (rule_pos_idx == rule_pos.get_sample_dens()) {
            rule_pos_bucket_idx++;
            rule_pos.get_inter_sampled_values(rule_pos_bucket_idx,
                                              rule_pos_values.data());
            rule_pos_sample_value = rule_pos.sample(rule_pos_bucket_idx);
            rule_pos_idx = 0;
            for (size_t i = 0; i < rule_pos_values.size(); i++) {
                rule_pos_values[i] += rule_pos_sample_value;
            }
        }
        uint64_t lcp_length = lcp_values[lcp_idx] - last_lcp;
        uint64_t rule_suffix_length =
            rule_pos_values[rule_pos_idx] - last_rule_pos;

        uint64_t total_rule_length = lcp_length + rule_suffix_length;

        for (uint64_t j = 0; j < lcp_length; j++) {
            gd.rule[rule_start + j] = gd.rule[prev_rule_start + j];
        }
        for (uint64_t j = lcp_length; j < total_rule_length; j++) {
            gd.rule[rule_start + j] = rule[start++];
        }

        last_lcp = lcp_values[lcp_idx];
        last_rule_pos = rule_pos_values[rule_pos_idx];
        rule_delim[rule_start] = 1;
        prev_rule_start = rule_start;
        rule_start += total_rule_length;
        lcp_idx++;
        rule_pos_idx++;
    }
    rule_delim[rule_delim.size() - 1] = 1;
    gd.rule_delim = std::move(rule_delim);
    sdsl::util::init_support(gd.rule_delim_sel, &gd.rule_delim);
    return gd;
}

void gcis_gap_codec_level::expand_rule(uint64_t rule_num,
                                       sdsl::int_vector<> &r_string,
                                       uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num + 1);
    uint64_t rule_length;
    for (rule_length = 1; rule_delim[rule_start + rule_length] == 0;
         rule_length++)
        ;
    for (uint64_t i = 0; i < rule_length; i++) {
        r_string[l++] = rule[rule_start + i];
    }
}

void gcis_gap_codec_level::expand_rule(uint64_t rule_num, char *s,
                                       uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num + 1);
    uint64_t rule_length;
    for (rule_length = 1; rule_delim[rule_start + rule_length] == 0;
         rule_length++)
        ;
    for (uint64_t i = 0; i < rule_length; i++) {
        s[l++] = rule[rule_start + i];
    }
}

uint64_t gcis_gap_codec::get_rule_pos(uint64_t i) {
    return i > 0 ? rule_pos[i - 1] : 0;
}

uint64_t gcis_gap_codec::get_rule_length(uint64_t i) {
    return i > 0 ? rule_pos[i] - rule_pos[i - 1] : 1;
}

uint64_t gcis_gap_codec::get_lcp(uint64_t i) { return lcp[i]; }

/**
 * @brief Extracts the LCP from a rule into a buffer.
 *
 * @param rule_num Rule number.
 * @param extracted_text The buffer.
 * @param k Marks the buffer's next available position.
 */
void gcis_gap_codec::extract_lcp(uint64_t rule_num,
                                 sdsl::int_vector<> &extracted_text,
                                 uint64_t &k) {

    // Gap decoded lcp values
    vector<uint64_t> buffer(
        max(lcp.get_sample_dens(), rule_pos.get_sample_dens()));
    // Necessary LCP values
    vector<uint64_t> lcp_values;
    // Determine the bucket in which the rule_num falls
    uint64_t sample_bucket_lcp = (rule_num / lcp.get_sample_dens());
    // Determine the rule's LCP idx within the bucket
    uint64_t idx_lcp = rule_num % lcp.get_sample_dens();
    // Extract Gap encoded LCP Values
    lcp.get_inter_sampled_values(sample_bucket_lcp, buffer.data());
    // LCP sample from the bucket
    uint64_t lcp_sample_value = lcp.sample(sample_bucket_lcp);
    uint64_t lcp_len;

    // Extract the necessary LCP values
    do {
        if (idx_lcp > 0) {
            lcp_len = buffer[idx_lcp] - buffer[idx_lcp - 1] + lcp_sample_value;
            idx_lcp--;
        } else {
            if (sample_bucket_lcp > 0) {
                uint64_t tmp = buffer[idx_lcp] + lcp_sample_value;
                sample_bucket_lcp--;
                lcp.get_inter_sampled_values(sample_bucket_lcp, buffer.data());
                lcp_sample_value = lcp.sample(sample_bucket_lcp);
                lcp_len = tmp - (buffer.back() + lcp_sample_value);
                idx_lcp = lcp.get_sample_dens() - 2;
            } else { // LCP[0] = 0;
                lcp_len = 0;
            }
        }
        lcp_values.push_back(lcp_len);
    } while (lcp_len > 0);

    // Gap decoded rule_pos values
    vector<uint64_t> rule_pos_values;
    // Determine the bucket in which the rule_num falls
    uint64_t sample_bucket_rule_pos = (rule_num / rule_pos.get_sample_dens());
    //  Determine the rule's rule pos idx within the bucket
    uint64_t idx_rule_pos = rule_num % rule_pos.get_sample_dens();
    // Extract Gap encoded rule pos values
    rule_pos.get_inter_sampled_values(sample_bucket_rule_pos, buffer.data());
    // Rule pos sample from the bucket
    uint64_t rule_pos_sample_value = rule_pos.sample(sample_bucket_rule_pos);

    // Extract the necessary rule pos values
    for (uint64_t i = 0; i < lcp_values.size(); i++) {
        uint64_t pos;
        if (idx_rule_pos > 0) {
            pos = buffer[idx_rule_pos - 1] + rule_pos_sample_value;
            idx_rule_pos--;
        } else {
            if (sample_bucket_rule_pos > 0) {
                sample_bucket_rule_pos--;
                rule_pos.get_inter_sampled_values(sample_bucket_rule_pos,
                                                  buffer.data());
                rule_pos_sample_value = rule_pos.sample(sample_bucket_rule_pos);
                pos = buffer.back() + rule_pos_sample_value;
                idx_rule_pos = rule_pos.get_sample_dens() - 2;
            } else { // Rule[0] starts at pos 0
                pos = 0;
            }
        }
        rule_pos_values.push_back(pos);
    }

    // Previous rule id
    uint64_t prev_rule = rule_num - 1;
    // Current LCP length
    uint64_t cur_lcp_len;

    lcp_len = cur_lcp_len = lcp_values[0];
    // Next character to be extracted (right-to-left)
    int64_t idx = k + lcp_len - 1;

    // While there is LCP symbols to be extracted
    for (uint64_t j = 1; j < lcp_values.size(); j++) {
        uint64_t prev_lcp_len = lcp_values[j];
        /* If the previous LCP length is lesser than the current LCP length
         * we can extract symbols from the previous rule */
        if (prev_lcp_len < cur_lcp_len) {
            int64_t i = cur_lcp_len - 1;
            uint64_t prev_rule_pos = rule_pos_values[j];
            while (i >= (int64_t)prev_lcp_len && idx >= (int64_t)k) {
                // Copy the symbols to the extracted text
                extracted_text[idx] = rule[prev_rule_pos + i - prev_lcp_len];
                i--;
                idx--;
            }
            cur_lcp_len = prev_lcp_len;
        }
        prev_rule--;
    }
    k += lcp_len;
}

void gcis_gap_codec::extract_rule_suffix(uint64_t rule_num,
                                         sdsl::int_vector<> &extracted_text,
                                         uint64_t &k) {

    uint64_t rule_suffix_len = get_rule_length(rule_num);
    uint64_t rule_pos = get_rule_pos(rule_num);
    for (int64_t i = 0; i < rule_suffix_len; i++) {
        extracted_text[k++] = rule[rule_pos + i];
    }
}

void gcis_gap_codec::extract_rule(uint64_t rule_num,
                                  sdsl::int_vector<> &extracted_text,
                                  uint64_t &k) {

    extract_lcp(rule_num, extracted_text, k);
    extract_rule_suffix(rule_num, extracted_text, k);
}

void gcis_gap_codec::serialize(std::ostream &o) {
    o.write((char *)&alphabet_size, sizeof(alphabet_size));
    o.write((char *)&string_size, sizeof(string_size));
    lcp.serialize(o);
    rule_pos.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
    fully_decoded_rule_len.serialize(o);
    uint64_t n = reduced_string_ps.size();
    o.write((char *)&n, sizeof(n));
    o.write((char *)reduced_string_ps.data(), n * sizeof(uint64_t));
    o.write((char *)&fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}
void gcis_gap_codec::load(std::istream &i) {
    i.read((char *)&alphabet_size, sizeof(alphabet_size));
    i.read((char *)&string_size, sizeof(string_size));
    lcp.load(i);
    rule_pos.load(i);
    rule.load(i);
    tail.load(i);
    fully_decoded_rule_len.load(i);
    uint64_t n;
    i.read((char *)&n, sizeof(n));
    reduced_string_ps.resize(n);
    i.read((char *)reduced_string_ps.data(), sizeof(uint64_t) * n);
    i.read((char *)&fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}