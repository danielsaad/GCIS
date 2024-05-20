//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_S8B_HPP
#define GC_IS_GCIS_S8B_HPP

#include "gcis.hpp"
#include "gcis_s8b_codec.hpp"
#include "util.hpp"
#include <print>
template <>
class gcis_dictionary<gcis_s8b_codec> : public gcis_abstract<gcis_s8b_codec> {
  public:
    // char *decode() override {
    //     sdsl::int_vector<> r_string = reduced_string;
    //     char *str = 0;
    //     for (int64_t i = g.size() - 1; i >= 0; i--) {
    //         sdsl::int_vector<> next_r_string;
    //         gcis_s8b_codec_level gd = std::move(g[i].decompress());
    //         gd.rule_delim_sel.set_vector(&gd.rule_delim);
    //         next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
    //         next_r_string.resize(g[i].string_size);
    //         uint64_t l = 0;
    //         if (i == 0) {
    //             // Convert the reduced string in the original text
    //             str = new char[g[i].string_size];
    //             for (auto t : g[i].tail) {
    //                 str[l++] = (char)t;
    //             }
    //             for (uint64_t j = 0; j < r_string.size(); j++) {
    //                 gd.expand_rule(r_string[j], str, l);
    //             }
    //         } else {
    //             // Convert the reduced string in the previous reduced string
    //             for (uint64_t j = 0; j < g[i].tail.size(); j++) {
    //                 next_r_string[l++] = g[i].tail[j];
    //             }
    //             for (uint64_t j = 0; j < r_string.size(); j++) {
    //                 gd.expand_rule(r_string[j], next_r_string, l);
    //             }
    //             r_string = std::move(next_r_string);
    //         }
    //     }
    //     return str;
    // }

    pair<char *, int_t> decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char *str = nullptr;
        for (int64_t i = g.size() - 1; i >= 0; i--) {
            sdsl::int_vector<> next_r_string;
            gcis_s8b_pointers_codec_level gd = std::move(g[i].decompress());
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
            next_r_string.resize(g[i].string_size);
            uint_t l = 0;
            if (i == 0) {
                // Convert the reduced string in the original text
                str = new char[g[i].string_size];
                for (auto t : g[i].tail) {
                    str[l++] = (char)t;
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], str, l);
                }
            } else {
                // Convert the reduced string in the previous reduced string
                for (uint64_t j = 0; j < g[i].tail.size(); j++) {
                    next_r_string[l++] = g[i].tail[j];
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], next_r_string, l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return make_pair(str, g[0].string_size);
    }

  private:
    void gc_is(int_t *s, uint_t *SA, int_t n, int_t K, int cs, int level) {
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

        int_t first = n;

        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets

        for (i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }

        for (i = n - 2; i >= 0; i--) {
            if (isLMS(i)) {
                // cout << chr(i) << " at position " << i << " is LMS\n";
                // cout << "Inserting at bucket position " << bkt[rank[chr(i)]]
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

        // find the lexicographic names of all LMS-substrings by comparing the
        // consecutive ones
        int_t name = -1;
        int_t prev = -1;

        int_t prev_len = -1;
        int_t cur_len = -1;

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        g.push_back(gcis_s8b_codec());
        // Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            cur_len = 1;
            while (pos + cur_len < n && !isLMS(pos + cur_len))
                cur_len++;
            bool diff = prev == -1 || prev_len != cur_len;
            int_t d;

            for (d = 0; d < min(cur_len, prev_len); d++) {
                if (chr(pos + d) != chr(prev + d)) {
                    diff = true;
                    break;
                }
            }

            // The consecutive LMS-substrings differs
            if (diff) {

                g[level].lcp.encode(d);
                g[level].rule_suffix_length.encode(cur_len - d - 1);
                g[level].rule.resize(g[level].rule.size() + cur_len - d - 1);

#ifdef REPORT
                total_rule_len += cur_len;
                total_lcp += d;
                total_rule_suffix_length += cur_len - d;
#endif

                for (j = 0; j < cur_len - d - 1 && j + pos + d < n; j++) {
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

class gcis_s8b_pointers : public gcis_dictionary<gcis_s8b_codec> {
  public:
    pair<char *, int_t> decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char *str = 0;
        for (int64_t i = g.size() - 1; i >= 0; i--) {
            gcis_s8b_pointers_codec_level gd = std::move(g[i].decompress());
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
                    // std::println("Number of rules = {}, Rule = {}, Level = {} "
                                //  "l = {}, next_r_string.size() = {}",
                                //  r_string.size(), j, i, l,
                                //  next_r_string.size());
                    gd.expand_rule(r_string[j], next_r_string, l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return make_pair(str, g[0].string_size);
    }

    pair<char *, int_t> decode_saca(uint_t **sa) override {

        sdsl::int_vector<> r_string = reduced_string;
        unsigned char *str;
        uint_t n = g[0].string_size;
        uint_t *SA = new uint_t[n];

        int_t *s = (int_t *)SA + n / 2;

        int cs = sizeof(int_t);
        if (g.size()) {

            for (int64_t level = g.size() - 1; level >= 0; level--) {
#if TIME
                auto start = timer::now();
#endif

                n = g[level].string_size;

                uint_t n1 = r_string.size();
                uint_t *SA1 = SA, *s1 = SA + n - n1;

                // copy to s1[1]
                if (level == g.size() - 1)
                    for (uint_t i = 0; i < n1; i++)
                        SA1[r_string[i]] = i;
                else
                    for (uint_t i = 0; i < n1; i++)
                        s1[i] = SA[i];

#if DEBUG
                cout << endl
                     << "level = " << level
                     << "\t string_size = " << g[level].string_size << "\n"
                     << "alphabet_size = " << g[level].alphabet_size << endl;
                cout << "n = " << n << "\nn1 = " << n1 << endl;
                cout << "\n####" << endl;
                cout << "s1 = ";
                for (uint_t i = 0; i < n1; i++) {
                    cout << s1[i] << " ";
                    cout << endl;
                }
#endif

#if TIME
                auto expand = timer::now();
#endif

                gcis_s8b_pointers_codec_level gd =
                    std::move(g[level].decompress());
                sdsl::int_vector<> next_r_string(
                    g[level].string_size, 0,
                    sdsl::bits::hi(g[level].alphabet_size - 1) + 1);
                uint_t l = 0;

                int_t K = g[level].alphabet_size; // alphabet

                int_t *bkt = new int_t[K]; // bucket
                int_t *cnt = new int_t[K]; // counters

                init_buckets(cnt, K);

                if (level == 0) {

                    // delete[] s;
                    // Convert the reduced string in the original text
                    str = new unsigned char[g[level].string_size];
                    for (uint64_t j = 0; j < g[level].tail.size(); j++) {
                        str[l] = g[level].tail[j];
                        cnt[str[l++]]++; // count frequencies
                    }
                    for (uint64_t j = 0; j < r_string.size(); j++) {
                        gd.expand_rule_bkt(r_string[j], str, l, cnt);
                    }
                    n = g[level].string_size;
                    // Classify the type of each character
                    uint_t cur_t, succ_t;
                    uint_t j = n1 - 1;
                    // s1[j--] = n - 1;
                    succ_t = 0; // s[n-2] must be L-type
                    for (uint_t i = n - 1; i > 0; i--) {
                        cur_t = (str[i - 1] < str[i] ||
                                 (str[i - 1] == str[i] && succ_t == 1))
                                    ? 1
                                    : 0;
                        if (cur_t == 0 && succ_t == 1)
                            s1[j--] = i;
                        succ_t = cur_t;
                    }
                } else {

                    init_buckets(bkt, K);

                    // Convert the reduced string in the previous reduced string
                    for (uint_t j = 0; j < g[level].tail.size(); j++) {
                        next_r_string[l++] = g[level].tail[j];
                        cnt[g[level].tail[j]]++; // count frequencies
                    }
                    for (uint_t j = 0; j < r_string.size(); j++) {
                        gd.expand_rule_bkt(r_string[j], next_r_string, l, cnt);
                    }
                    r_string = std::move(next_r_string);

                    // n=r_string.size();
                    n = g[level].string_size;
                    // copy to s[1]
                    for (uint_t i = 0; i < n; i++)
                        s[i] = r_string[i];

                    // Classify the type of each character
                    uint_t cur_t, succ_t;
                    uint_t j = n1 - 1;
                    // s1[j--] = n - 1;
                    succ_t = 0; // s[n-1] must be L-type
                    for (uint_t i = n - 1; i > 0; i--) {
                        cur_t =
                            (r_string[i - 1] < r_string[i] ||
                             (r_string[i - 1] == r_string[i] && succ_t == 1))
                                ? 1
                                : 0;
                        if (cur_t == 0 && succ_t == 1)
                            s1[j--] = i;
                        succ_t = cur_t;
                    }
                }

#if TIME
                auto end = timer::now();
                cout << "expand: "
                     << (double)duration_cast<seconds>(end - expand).count()
                     << " seconds" << endl;
#endif

#if DEGUB
                cout << "n = " << n << "\nn1 = " << n1 << endl;
                cout << "level = " << level
                     << "\t string_size = " << g[level].string_size << "\n"
                     << "alphabet_size = " << g[level].alphabet_size << endl;
                cout << "SA1: ";
                for (uint_t i = 0; i < n1; i++)
                    cout << SA[i] << ", ";
                cout << endl;
#endif

                // stage 3: induce the result for the original problem
                get_buckets(cnt, bkt, K, true);

#if TIME
                auto begin = timer::now();
#endif

#if TIME
                end = timer::now();
                cout << "classify: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

#if DEGUB
                for (int_t i = 0; i < n; i++) {
                    cout << tget(i) << " ";
                }
                cout << endl;
#endif

                int_t j = 0;
                for (int_t i = 0; i < n1; i++) {
                    SA1[i] = s1[SA1[i]]; // get index in s1
                }
                for (int_t i = n1; i < n; i++) {
                    SA[i] = EMPTY; // init SA[n1..n-1]
                }

                if (level) {
                    for (int_t i = n1 - 1; i >= 0; i--) {
                        j = SA[i];
                        SA[i] = EMPTY;
                        SA[bkt[chr(j)]--] = j;
                    }
                } else {
                    for (int_t i = n1 - 1; i >= 0; i--) {
                        j = SA[i];
                        SA[i] = EMPTY;
                        SA[bkt[str[j]]--] = j;
                    }
                    // SA[0] = n - 1;
                }

#if TIME
                end = timer::now();
                cout << "position: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

                if (level)
                    induceSAl(SA, s, cnt, bkt, n, K, cs, level);
                else
                    induceSAl(SA, (int_t *)str, cnt, bkt, n, K,
                              sizeof(unsigned char), level);

#if TIME
                end = timer::now();
                cout << "induce L: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

                if (level)
                    induceSAs(SA, s, cnt, bkt, n, K, cs, level);
                else
                    induceSAs(SA, (int_t *)str, cnt, bkt, n, K,
                              sizeof(unsigned char), level);

#if TIME
                end = timer::now();
                cout << "induce S: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
#endif

#if DEGUB
                cout << "SA: ";
                for (uint_t i = 0; i < n; i++) {
                    cout << SA[i] << ", ";
                }
                cout << endl;
#endif
                delete[] bkt;
                delete[] cnt;

#if TIME
                auto stop = timer::now();
                cout << "time: "
                     << (double)duration_cast<seconds>(stop - start).count()
                     << " seconds" << endl;
#endif
            }
        } else {
            str = new unsigned char[reduced_string.size()];
            for (uint64_t i = 0; i < reduced_string.size(); i++)
                str[i] = reduced_string[i];
        }

        *sa = SA;
        return make_pair((char *)str, g[0].string_size);
    } // end decode_sac
};

class gcis_lyndon : public gcis_s8b_pointers {
  public:
    pair<char *, int_t> decode_lyndon(int_t **lyn) {
        sdsl::int_vector<> r_string = reduced_string;
        unsigned char *str;
        uint_t n = g[0].string_size;
        uint_t *SA = new uint_t[n];
        // allocate space for LA
        int_t *LA = new int_t[n];
        // Initialize LA
        for (int i = 0; i < n; i++)
            LA[i] = i + 1;

        int_t *s = (int_t *)SA + n / 2;

        int cs = sizeof(int_t);
        if (g.size()) {

            for (int64_t level = g.size() - 1; level >= 0; level--) {
#if TIME
                auto start = timer::now();
#endif

                n = g[level].string_size;

                uint_t n1 = r_string.size();
                uint_t *SA1 = SA, *s1 = SA + n - n1;

                // copy to s1[1]
                if (level == g.size() - 1)
                    for (uint_t i = 0; i < n1; i++)
                        SA1[r_string[i]] = i;
                else
                    for (uint_t i = 0; i < n1; i++)
                        s1[i] = SA[i];

#if DEBUG
                cout << endl
                     << "level = " << level
                     << "\t string_size = " << g[level].string_size << "\n"
                     << "alphabet_size = " << g[level].alphabet_size << endl;
                cout << "n = " << n << "\nn1 = " << n1 << endl;
                cout << "\n####" << endl;
                cout << "s1 = ";
                for (uint_t i = 0; i < n1; i++) {
                    cout << s1[i] << " ";
                    cout << endl;
                }
#endif

#if TIME
                auto expand = timer::now();
#endif

                gcis_s8b_pointers_codec_level gd =
                    std::move(g[level].decompress());
                sdsl::int_vector<> next_r_string(
                    g[level].string_size, 0,
                    sdsl::bits::hi(g[level].alphabet_size - 1) + 1);
                uint_t l = 0;

                int_t K = g[level].alphabet_size; // alphabet

                int_t *bkt = new int_t[K]; // bucket
                int_t *cnt = new int_t[K]; // counters

                init_buckets(cnt, K);

                if (level == 0) {

                    // delete[] s;
                    // Convert the reduced string in the original text
                    str = new unsigned char[g[level].string_size];
                    for (uint64_t j = 0; j < g[level].tail.size(); j++) {
                        str[l] = g[level].tail[j];
                        cnt[str[l++]]++; // count frequencies
                    }
                    for (uint64_t j = 0; j < r_string.size(); j++) {
                        gd.expand_rule_bkt(r_string[j], str, l, cnt);
                    }
                    n = g[level].string_size;
                    // Classify the type of each character
                    uint_t cur_t, succ_t;
                    uint_t j = n1 - 1;
                    // s1[j--] = n - 1;
                    succ_t = 0; // s[n-1] must be L-type
                    for (uint_t i = n - 1; i > 0; i--) {
                        cur_t = (str[i - 1] < str[i] ||
                                 (str[i - 1] == str[i] && succ_t == 1))
                                    ? 1
                                    : 0;
                        if (cur_t == 0 && succ_t == 1)
                            s1[j--] = i;
                        succ_t = cur_t;
                    }
                } else {

                    init_buckets(bkt, K);

                    // Convert the reduced string in the previous reduced string
                    for (uint_t j = 0; j < g[level].tail.size(); j++) {
                        next_r_string[l++] = g[level].tail[j];
                        cnt[g[level].tail[j]]++; // count frequencies
                    }
                    // std::println("Level = {}", level);
                    // std::println("r_string_size() = {}", r_string.size());
                    for (uint_t j = 0; j < r_string.size(); j++) {
                        gd.expand_rule_bkt(r_string[j], next_r_string, l, cnt);
                    }
                    r_string = std::move(next_r_string);

                    // n=r_string.size();
                    n = g[level].string_size;
                    // copy to s[1]
                    for (uint_t i = 0; i < n; i++)
                        s[i] = r_string[i];

                    // Classify the type of each character
                    uint_t cur_t, succ_t;
                    uint_t j = n1 - 1;
                    // s1[j--] = n - 1;
                    succ_t = 0; // s[n-1] must be L-type
                    for (uint_t i = n - 1; i > 0; i--) {
                        cur_t =
                            (r_string[i - 1] < r_string[i] ||
                             (r_string[i - 1] == r_string[i] && succ_t == 1))
                                ? 1
                                : 0;
                        if (cur_t == 0 && succ_t == 1)
                            s1[j--] = i;
                        succ_t = cur_t;
                    }
                }

#if TIME
                auto end = timer::now();
                cout << "expand: "
                     << (double)duration_cast<seconds>(end - expand).count()
                     << " seconds" << endl;
#endif

#if DEGUB
                cout << "n = " << n << "\nn1 = " << n1 << endl;
                cout << "level = " << level
                     << "\t string_size = " << g[level].string_size << "\n"
                     << "alphabet_size = " << g[level].alphabet_size << endl;
                cout << "SA1: ";
                for (uint_t i = 0; i < n1; i++)
                    cout << SA[i] << ", ";
                cout << endl;
#endif

                // stage 3: induce the result for the original problem
                get_buckets(cnt, bkt, K, true);

#if TIME
                auto begin = timer::now();
#endif

#if TIME
                end = timer::now();
                cout << "classify: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

#if DEGUB
                for (int_t i = 0; i < n; i++) {
                    cout << tget(i) << " ";
                }
                cout << endl;
#endif

                int_t j = 0;
                for (int_t i = 0; i < n1; i++) {
                    SA1[i] = s1[SA1[i]]; // get index in s1
                }
                for (int_t i = n1; i < n; i++) {
                    SA[i] = EMPTY; // init SA[n1..n-1]
                }

                if (level) {
                    for (int_t i = n1 - 1; i >= 0; i--) {
                        j = SA[i];
                        SA[i] = EMPTY;
                        SA[bkt[chr(j)]--] = j;
                    }
                } else {
                    for (int_t i = n1 - 1; i >= 0; i--) {
                        j = SA[i];
                        SA[i] = EMPTY;
                        SA[bkt[str[j]]--] = j;
                    }
                    // SA[0] = n - 1;
                }

#if TIME
                end = timer::now();
                cout << "position: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

                if (level)
                    induceSAl(SA, s, cnt, bkt, n, K, cs, level);
                else
                    induceSAl(SA, (int_t *)str, cnt, bkt, n, K,
                              sizeof(unsigned char), level);

#if TIME
                end = timer::now();
                cout << "induce L: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

                if (level)
                    induceSAs(SA, s, cnt, bkt, n, K, cs, level);
                else
                    induceSAs_lyndon(SA, LA, (int_t *)str, cnt, bkt, n, K,
                                     sizeof(unsigned char), level);

#if TIME
                end = timer::now();
                cout << "induce S: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
#endif

#if DEGUB
                cout << "SA: ";
                for (uint_t i = 0; i < n; i++) {
                    cout << SA[i] << ", ";
                }
                cout << endl;
#endif
                delete[] bkt;
                delete[] cnt;

#if TIME
                auto stop = timer::now();
                cout << "time: "
                     << (double)duration_cast<seconds>(stop - start).count()
                     << " seconds" << endl;
#endif
                // std::println("Final suffix array\n");
                // if (level == 0) {
                //     for (int i = 0; i < n; i++) {
                //         std::println("SA[{}] = {}", i, SA[i]);
                //     }
                // }
            }
        } else {
            str = new unsigned char[reduced_string.size()];
            for (uint64_t i = 0; i < reduced_string.size(); i++)
                str[i] = reduced_string[i];
        }
        *lyn = LA;
        return make_pair((char *)str, g[0].string_size);
    } // end decode_lyndon

  private:
    void induceSAs_lyndon(uint_t *SA, int_t *LA, int_t *s, int_t *cnt,
                          int_t *bkt, uint_t n, int_t K, int_t cs, int level) {
        get_buckets(cnt, bkt, K, true);
        int_t i, j;
        for (i = n - 1; i > 0; i--) {
            int_t next, prev;
            j = SA[i];
            // std::println("SA[{}] = {}", i, SA[i]);
            if (SA[i] > 0) {
                if (LA[j - 1] < j)
                    prev = LA[j - 1];
                else
                    prev = j - 1;
                next = LA[j];
                LA[next - 1] = prev;
                // std::println("(SA[i] > 0) Computed LA[{}] = {}", next - 1,
                //              prev);
                if (prev >= 0) {
                    LA[prev] = next;
                    // std::println("(SA[i] > 0) Computed LA[{}] = {}", prev,
                    //              next);
                }

                j = SA[i] - 1;
                if (chr(j) <= chr(j + 1) && bkt[chr(j)] < i) { // Inducing
                                                               // S-type
                    // std::println("Inducing S-suffix {} at pos {}", j,
                                //  bkt[chr(j)]);

                    SA[bkt[chr(j)]] = j;
                    bkt[chr(j)]--;
                }
            } else {
                prev = j - 1;
                next = LA[j];
                LA[next - 1] = prev;
                // std::println("Computed LA[{}] = {}", next - 1, prev);
            }
        }
        LA[n - 1] = n;
        // for (int j = 0; j < n; j++) {
        //     std::println("Intermediary LA[{}] = {}", j, LA[j]);
        // }
        for (j = 0; j < n; j++) {
            if (LA[j] < j || LA[j] >= n) {
                LA[j] = 1;
                // std::println("Final computing LA[{}] = 1", j);
            } else {
                LA[j] = LA[j] - j;
                // std::println("Final computing LA[{}] = {}", j, LA[j]);
            }
        }
    }

    void gc_is(int_t *s, uint_t *SA, int_t n, int_t K, int cs, int level) {
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

        int_t first = n;

        std::println("Inducing GCIS Lyndon");
        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets

        for (i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }

        for (i = n - 2; i >= 0; i--) {
            if (isLMS(i)) {
                // cout << chr(i) << " at position " << i << " is LMS\n";
                // cout << "Inserting at bucket position " << bkt[rank[chr(i)]]
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

        // find the lexicographic names of all LMS-substrings by comparing the
        // consecutive ones
        int_t name = -1;
        int_t prev = -1;

        int_t prev_len = 0;
        int_t cur_len = -1;

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        g.push_back(gcis_s8b_codec());
        // Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            cur_len = 1;
            while (pos + cur_len < n && !isLMS(pos + cur_len))
                cur_len++;
            cur_len++;
            bool diff = prev == -1 || cur_len != prev_len;
            int_t d;
            //  d equals to the LCP between two consecutive LMS-substrings
            for (d = 0; d < min(prev_len, cur_len)-1; d++) {
                // If is first suffix in LMS order (sentinel), or one of the
                // suffixes reached the last position of T, or the
                // characters of T differs or the type os suffixes differ.
                if (pos + d == n || prev + d == n ||
                    chr(pos + d) != chr(prev + d)) {
                    diff = true;
                    break;
                }
            }

            if (diff) {
                std::println("d = {} cur_len = {}", d, cur_len);
                if (cur_len - d - 1 <= 0)
                    std::println("{} - 1 - {} = {}", cur_len, d,
                                 cur_len - d - 1);
                g[level].lcp.encode(d);
                g[level].rule_suffix_length.encode(cur_len - d - 1);
                g[level].rule.resize(g[level].rule.size() + cur_len - d - 1);

#ifdef REPORT
                total_rule_len += cur_len - 1;
                total_lcp += d;
                total_rule_suffix_length += cur_len - d - 1;
#endif
                // std::print("Rule {}: ", name + 1);
                // std::print("LCP = {}, ", d);
                for (j = 0; j < cur_len -d - 1 && j + pos + d < n; j++) {
#ifdef REPORT
                    if (j + pos + d + 1 < n &&
                        chr(j + pos + d) == chr(j + pos + d + 1)) {
                        run_length_potential++;
                    }
#endif

                    // std::print("{} ", (uint_t)chr(j + pos + d));
                    g[level].rule[rule_index] = (uint_t)chr(j + pos + d);
                    rule_index++;
                }
                // std::println("");
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
            // std::println("LMS string with name {} at pos {}", name, pos);
            // for (int k = 0; k < cur_len - 1; k++) {
            //     if (K == 256)
            //         std::print("{} ", (char)chr(pos + k));
            //     else
            //         std::print("{} ", (int)chr(pos + k));
            // }
            // std::print("\n");

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

#endif // GC_IS_GCIS_S8B_HPP
