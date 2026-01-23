#include "gcis.hpp"
#include "gcis_s8b_codec.hpp"
#include "simple8b.hpp"

class gcis_s8b_repair_hybrid_codec_level {
  public:
    sdsl::int_vector<> rule; // concatenated rules
    vector<uint_t> rule_pos; // position of concatenated rules

  public:
    /**
     * @brief Expands dictionary rule into a int_vector r_string.
     * Updates automatically the current index for r_string
     *
     * @param rule_num rule to be decompressed
     * @param r_string destination of the rule's expansion
     * @param l index for r_string
     */
    void expand_rule(uint_t rule_num, sdsl::int_vector<> &r_string, uint_t &l) {
        uint_t rule_start = rule_pos[rule_num];
        uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
        for (uint_t i = 0; i < rule_length; i++) {
            r_string[l] = rule[rule_start + i];
            l++;
        }
    }

    /**
     * @brief Expands dictionary rule into a C-string.
     * Updates automatically the current index for the final string.
     *
     * @param rule_num rule to be decompressed
     * @param s destination of the rule's expansion
     * @param l index for s
     */
    void expand_rule(uint_t rule_num, char *s, uint_t &l) {
        uint_t rule_start = rule_pos[rule_num];
        uint_t rule_length = rule_pos[rule_num + 1] - rule_pos[rule_num];
        for (uint_t i = 0; i < rule_length; i++) {
            s[l] = rule[rule_start + i];
            l++;
        }
    }
};

class gcis_s8b_repair_hybrid_codec : public gcis_s8b_codec {
  public:
    map<uint, tuple<uint, uint>> repair_dict;

  public:
    /**
     * @brief Stores the object into a stream
     *
     * @param o output stream
     */
    void serialize(std::ostream &o) {
        gcis_s8b_codec::serialize(o);
        uint repair_dict_size = repair_dict.size();
        o.write((char *)&repair_dict_size, sizeof(repair_dict_size));
        for (const auto &[k, v] : repair_dict) {
            const auto [p1, p2] = v;
            o.write((char *)&k, sizeof(k));
            o.write((char *)&p1, sizeof(p1));
            o.write((char *)&p2, sizeof(p2));
        }
    }
    /**
     * @brief Loads the object from a file
     *
     * @param i input stream
     */
    void load(std::istream &i) {
        gcis_s8b_codec::load(i);
        uint sz;
        i.read((char *)&sz, sizeof(sz));
        for (uint j = 0; j < sz; j++) {
            uint k, p1, p2;
            i.read((char *)&k, sizeof(k));
            i.read((char *)&p1, sizeof(p1));
            i.read((char *)&p2, sizeof(p2));
            repair_dict.insert({k, {p1, p2}});
        }
    }

    /**
     * @brief Decompresses a dictionary level entirely, making decompression
     * faster.
     *
     * @return gcis_s8b_repair_hybrid_codec_level returns the dictionary level
     * decompressed allowing fast acess.
     */
    gcis_s8b_repair_hybrid_codec_level decompress() {
        gcis_s8b_repair_hybrid_codec_level gd;
        uint64_t number_of_rules = rule_suffix_length.size();
        uint64_t total_lcp_length = 0;
        uint64_t total_rule_suffix_length = 0;
        // Compute the total LCP length
        lcp.reset();
        rule_suffix_length.reset();
        for (uint64_t i = 0; i < lcp.size(); i++) {
            total_lcp_length += lcp.get_next();
        }
        // Compute the total rule suffix length and the number of rules
        for (uint64_t i = 0; i < rule_suffix_length.size(); i++) {
            total_rule_suffix_length += rule_suffix_length.get_next();
        }

        // Resize data structures
        uint64_t total_length = total_lcp_length + total_rule_suffix_length;
        gd.rule.width(sdsl::bits::hi(repair_dict.size() + alphabet_size - 1) +
                      1);
        gd.rule.resize(total_length);
        uint64_t rule_start = 0;
        uint64_t prev_rule_start = 0;
        uint64_t start = 0;
        lcp.reset();
        rule_suffix_length.reset();

        for (uint64_t i = 0; i < repair_dict.size(); i++) {
            const auto [p1, p2] = repair_dict[i];
            gd.rule[rule_start] = p1;
            gd.rule[rule_start + 1] = p2;
            prev_rule_start = rule_start;
            rule_start += 2;
        }
        for (uint64_t i = repair_dict.size(); i < number_of_rules; i++) {
            int64_t k;

            uint64_t lcp_length = lcp.get_next();
            uint64_t rule_length = rule_suffix_length.get_next();

            total_length = lcp_length + rule_length;

            // Copy the contents of the previous rule by LCP length chars
            uint64_t j = 0;
            while (j < lcp_length) {
                gd.rule[rule_start + j] = gd.rule[prev_rule_start + j];
                j++;
            }

            // Copy the remaining suffix rule
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
};

class gcis_repair_hybrid_2
    : public gcis_abstract<gcis_s8b_repair_hybrid_codec> {
  public:
  public:
    void encode(char *s, int_t n) override {
        uint_t *SA = new uint_t[n];
        int_t K = 256;
        int_t *s2 = new int_t[n];
        for (int_t i = 0; i < n; i++) {
            s2[i] = s[i];
        }
        int cs = sizeof(int_t);
        int level = 0;
        gc_is(s2, SA, n, K, cs, level);
        delete[] SA;
    }

    pair<char *, int_t> decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char *str = 0;
        for (int64_t i = g.size() - 1; i >= 0; i--) {
            gcis_s8b_repair_hybrid_codec_level gd =
                std::move(g[i].decompress());
            sdsl::int_vector<> next_r_string;
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
            next_r_string.resize(g[i].string_size);
            uint_t l = 0;
            if (i == 0) {
                // Convert the reduced string in the original text
                str = new char[g[i].string_size];
                for (auto t : g[i].tail) {
                    str[l++] = (char)t;
                }
                for (uint_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], str, l);
                }
            } else {
                // Convert the reduced string in the previous reduced string
                cout << "Tail size = " << g[i].tail.size() << endl;
                cout << "String size = " << g[i].string_size << endl;
                for (uint_t j = 0; j < g[i].tail.size(); j++) {
                    next_r_string[l++] = g[i].tail[j];
                }
                for (uint_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], next_r_string, l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return make_pair(str, g[0].string_size);
    }

  private:
    void gc_is(int_t *s, uint_t *SA, int_t n, int_t K, int cs,
               int level) override {
        int_t i, j;
#ifdef MEM_MONITOR
        mm.event("GC-IS Level " + to_string(level));
#endif

#ifdef REPORT
        uint_t total_lcp = 0;
        uint_t total_rule_suffix_length = 0;
        uint_t run_length_potential = 0;
        uint_t total_rule_len = 0;
        uint_t discarded_rules_n = 0;
        uint_t discarded_rules_len = 0;
#endif

        g.push_back(gcis_s8b_repair_hybrid_codec());
        auto &codec = g.back();

        map<tuple<uint, uint>, uint> freq;
        for (i = 0; i < n - 1; i++) {
            tuple<uint, uint> t = {chr(i), chr(i + 1)};
            freq[t]++;
        }
        vector<pair<uint, tuple<uint, uint>>> freq_arr;
        for (const auto &[k, v] : freq) {
            freq_arr.push_back({v, k});
        }
        sort(freq_arr.begin(), freq_arr.end(),
             [](const auto &a, const auto &b) -> bool {
                 return a.first > b.first;
             });

        /** Get most frequent pair and insert it into repair's dictionary **/
        const auto [p1, p2] = get<1>(freq_arr[0]);
        codec.repair_dict.insert({0, {p1, p2}});
        freq_arr.clear();
        freq.clear();

        cout << "n K " << n << ' ' << K << '\n';

        /** Alphabet change **/

        for (i = 0, j = 0; i < n;) {
            if (i < n - 1 && chr(i) == p1 && chr(i + 1) == p2) {
                s[j++] = 0;
                i += 2;
            } else {
                s[j++] = s[i++] + 1;
            }
        }
        n = j;
        K = K + 1;
        cout << "n K " << n << ' ' << K << '\n';

        unsigned char *t =
            new unsigned char[n / 8 + 1]; // LS-type array in bits
        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        //  tset(n - 2, 0);
        tset(n - 1, 0); // the last symbol is L-type

        for (i = n - 2; i >= 0; i--) {
            tset(i, (chr(i) < chr(i + 1) ||
                     (chr(i) == chr(i + 1) && tget(i + 1) == 1))
                        ? 1
                        : 0);
        }

        int_t *bkt = new int_t[K]; // bucket counters

        int_t first = n - 1;

        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets

        for (i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }

        for (i = n - 2; i >= 0; i--) {
            if (isLMS(i)) {
                // cout << chr(i) << " at position " << i << " is LMS\n";
                // cout << "Inserting at bucket position " <<
                // bkt[rank[chr(i)]]
                //      << '\n';
                SA[bkt[chr(i)]--] = i;
                first = i;
            }
        }

        // SA[0] = n - 1; // set the single sentinel LMS-substring

        // Induce L-Type suffixes by using LMS-Type and L-Type suffixes
        induceSAl(t, SA, s, bkt, n, K, cs, level);

        // Induce S-Type suffixes by using L-Type and S-Type suffixes
        induceSAs(t, SA, s, bkt, n, K, cs, level);

        delete[] bkt;

        // compact all the sorted substrings into the first n1 items of s
        // 2*n1 must be not larger than n (proveable)
        // n1 contains the end of the lms positions
        int_t n1 = 0;
        for (i = 0; i < n; i++) {
            if (isLMS(SA[i])) {
                SA[n1++] = SA[i];
            }
        }

        // Init the name array buffer
        // SA[0,n1-1] = LMS starting positions
        // SA[n1,n-1] = Name for each LMS substring
        for (i = n1; i < n; i++) {
            SA[i] = EMPTY;
        }

        // find the lexicographic names of all LMS-substrings by comparing
        // the consecutive ones
        int_t name = -1;
        int_t prev = -1;

        int_t prev_len = -1;
        int_t cur_len = -1;

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        // Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            cur_len = 1;
            while (pos + cur_len < n && !isLMS(pos + cur_len))
                cur_len++;
            bool diff = false;
            int_t d;
            if (prev == -1 || prev_len != cur_len)
                diff = true;

            for (d = 0; d < min(cur_len, prev_len); d++) {
                if (chr(pos + d) != chr(prev + d)) {
                    diff = true;
                    break;
                }
            }

            // The consecutive LMS-substrings differs
            if (diff) {

                g[level].lcp.encode(d);
                g[level].rule_suffix_length.encode(cur_len - d);
                g[level].rule.resize(g[level].rule.size() + cur_len - d);

#ifdef REPORT
                total_rule_len += cur_len;
                total_lcp += d;
                total_rule_suffix_length += cur_len - d;
#endif

                for (j = 0; j < cur_len - d && j + pos + d < n; j++) {
#ifdef REPORT
                    if (j + pos + d + 1 < n &&
                        chr(j + pos + d) == chr(j + pos + d + 1)) {
                        run_length_potential++;
                    }
#endif
                    g[level].rule[rule_index] = (uint_t)chr(j + pos + d);
                    rule_index++;
                }
                name++;
                prev = pos;
                prev_len = cur_len;
            }
#ifdef REPORT
            else {
                size_t len = 1;
                if (pos != n - 1)
                    while (!isLMS(pos + len))
                        len++;
                discarded_rules_len += len;
                discarded_rules_n++;
            }
#endif
            pos = (pos % 2 == 0) ? pos / 2 : (pos - 1) / 2;
            SA[n1 + pos] = name;
        }

        sdsl::util::bit_compress(g[level].rule);
        g[level].lcp.encode();
        g[level].rule_suffix_length.encode();

        for (i = n - 1, j = n - 1; i >= n1; i--) {
            if (SA[i] != EMPTY) {
                SA[j--] = SA[i];
            }
        }

        // s1 is done now
        uint_t *SA1 = SA, *s1 = SA + n - n1;

        // Copy the first elements (not part of a LMS substring)
        g[level].tail.resize(first);
        for (j = 0; j < first; j++) {
            g[level].tail[j] =
                (uint64_t)(cs == sizeof(char) ? ((char *)s)[j] : s[j]);
        }
        sdsl::util::bit_compress(g[level].tail);

        // stage 2: solve the reduced problem
        // recurse if names are not yet unique

#ifdef REPORT
        gcis::util::print_report("Level ", level, "\n");
        gcis::util::print_report("Alphabet Size = ", K, "\n");
        gcis::util::print_report("String Size = ", n, "\n");
        gcis::util::print_report("Number of Rules = ", name + 1, "\n");
        gcis::util::print_report("Average Rule Length = ",
                                 (double)total_rule_len / (name + 1), "\n");
        gcis::util::print_report(
            "Number of Discarded Rules = ", discarded_rules_n, "\n");
        gcis::util::print_report(
            "Average Discarded Rules Length = ",
            (double)discarded_rules_len / discarded_rules_n, "\n");
        gcis::util::print_report(
            "Average LCP = ", (double)total_lcp / (name + 1), "\n");
        gcis::util::print_report("Average Rule Suffix Length = ",
                                 (double)total_rule_suffix_length / (name + 1),
                                 "\n");
        gcis::util::print_report(
            "Dictionary Level Size (bytes) =", g[level].size_in_bytes(), "\n");
        gcis::util::print_report("LCP Size (bits) = ", g[level].lcp.size(),
                                 "\n");
        gcis::util::print_report(
            "Rule Suffix Length (total) = ", g[level].rule.size(), "\n");
        gcis::util::print_report("Rule Suffix Width (bits per symbol) = ",
                                 (int_t)g[level].rule.width(), "\n");
        gcis::util::print_report("Tail Length = ", g[level].tail.size(), "\n");
        gcis::util::print_report("Tail Width (bits per symbol) = ",
                                 (int_t)g[level].tail.width(), "\n");
        gcis::util::print_report(
            "Run Length Potential (total) = ", run_length_potential, "\n");
        gcis::util::print_report("Avg Run Length per Rule Suffix = ",
                                 (double)run_length_potential / (name + 1),
                                 "\n");
#endif

        // bool premature_stop =
        //     evaluate_premature_stop(n, K, n1, name + 1, level);
        bool premature_stop = false;
        g[level].string_size = n;
        g[level].alphabet_size = K;

        if (name + 1 < n1 && !premature_stop) {
            g[level].string_size = n;
            g[level].alphabet_size = K;
            gc_is((int_t *)s1, SA1, n1, name + 1, sizeof(int_t), level + 1);
        } else { // generate the suffix array of s1 directly
            if (premature_stop) {
#ifdef REPORT
                gcis::util::print_report("Premature Stop employed at level ",
                                         level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s[j];
                }
                g.pop_back();
            } else {
                reduced_string.resize(n1);
                for (j = 0; j < n1; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s1[j];
                }
            }
            sdsl::util::bit_compress(reduced_string);

#ifdef REPORT
            gcis::util::print_report(
                "Reduced String Length = ", (int_t)reduced_string.size(), "\n");
            gcis::util::print_report(
                "Reduced String Width (bits per symbol) = ",
                (int_t)reduced_string.width(), "\n");
#endif
        }
        delete[] t;
    }
};