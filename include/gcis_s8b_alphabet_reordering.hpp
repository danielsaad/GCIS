#include "gcis_s8b.hpp"
#include "gcis_statistics.hpp"

class gcis_s8b_pointers_alpha : public gcis_s8b_pointers {
  public:
    gcis_statistics stats;

  public:
    using gcis_s8b_pointers::get_buckets;
    using gcis_s8b_pointers::induceSAl;
    using gcis_s8b_pointers::induceSAs;
    
    void get_buckets(int_t *s, int_t *bkt, int_t n, int_t K, int cs, bool end,
                     map<int, int> &rank) {
        int_t i, sum = 0;

        // clear all buckets
        init_buckets(bkt, K);
        // compute the size of each bucket
        for (i = 0; i < n; i++) {
            bkt[rank[chr(i)]]++;
        }
        // Mark the end of each bucket
        if (end) {
            for (i = 0; i < K; i++) {
                sum += bkt[i];
                bkt[i] = sum - 1;
            }
        } else {
            for (i = 0; i < K; i++) {
                sum += bkt[i];
                bkt[i] = sum - bkt[i];
            }
        }
        // cout << "Get buckets\n";
        // for (int i = 0; i < n; i++) {
        //     cout << "bkt[" << i << "] = " << bkt[i] << '\n';
        // }
    }

    void induceSAs(unsigned char *t, uint_t *SA, int_t *s, int_t *bkt, int_t n,
                   int_t K, int cs, int level, map<int, int> &rank) {
        int_t i, j;
        // cout << "inducing SAs" << endl;
        get_buckets(s, bkt, n, K, cs, true, rank); // find ends of buckets
        for (i = n - 1; i >= 0; i--) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0 && tget(j)) {
                    // cout << chr(j) << " at position " << j
                    //      << " induced at position " << bkt[rank[chr(j)]] <<
                    //      endl;
                    SA[bkt[rank[chr(j)]]--] = j;
                }
            }
        }
    }

    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the
    // L-Suffixes
    void induceSAl(unsigned char *t, uint_t *SA, int_t *s, int_t *bkt, int_t n,
                   int_t K, int cs, int level, map<int, int> &rank) {
        int_t i, j;
        // find heads of buckets
        get_buckets(s, bkt, n, K, cs, false, rank);

        // cout << "inducing SAl" << endl;
        // cout << chr(n - 1) << " at position " << n - 1
        //      << " induced at position " << bkt[rank[chr(n - 1)]] << endl;
        SA[bkt[rank[chr(n - 1)]]++] = n - 1;
        //  if(level==0) bkt[0]++;
        for (i = 0; i < n; i++) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0 && !tget(j)) {
                    // cout << chr(j) << " at position " << j
                    //      << " induced at position " << bkt[rank[chr(j)]]
                    //      << endl;
                    SA[bkt[rank[chr(j)]]++] = j;
                }
            }
        }
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
        stats.lvl_stats.push_back(gcis_level_statistics());
        auto &lvl_stats = stats.lvl_stats.back();
#endif

        unsigned char *t =
            new unsigned char[n / 8 + 1]; // LS-type array in bits
        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        //  tset(n - 2, 0);
        tset(n - 1, 0); // the last symbol is L-type
        map<tuple<int, int>, int> freq;
        for (i = n - 1; i > 0; i--) {
            tuple<int, int> t = {chr(i), chr(i - 1)};
            freq[t]++;
        }
        vector<pair<int, tuple<int, int>>> freq_arr;
        for (const auto &[k, v] : freq) {
            freq_arr.push_back({v, k});
        }
        sort(freq_arr.begin(), freq_arr.end(),
             [](const pair<int, tuple<int, int>> &a,
                const pair<int, tuple<int, int>> &b) -> bool {
                 return a.first > b.first;
             });
        for (const auto &[f, t] : freq_arr) {
            // cout << get<0>(t) << ' ' << get<1>(t) << " = " << f << endl;
        }
        map<int, int> rank;
        int r = 0;
        for (const auto &[f, t] : freq_arr) {
            if (rank.find(get<0>(t)) == rank.end()) {
                rank[get<0>(t)] = r++;
            }
        }
        for (const auto &[k, v] : rank) {
            // cout << k << " rank: " << v << '\n';
        }
        freq_arr.clear();
        freq.clear();

        for (i = n - 2; i >= 0; i--) {
            tset(i, (rank[chr(i)] < rank[chr(i + 1)] ||
                     (rank[chr(i)] == rank[chr(i + 1)] && tget(i + 1) == 1))
                        ? 1
                        : 0);
        }
        int_t *bkt = new int_t[K]; // bucket counters

        int_t first = n - 1;

        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true,rank); // find ends of buckets

        for (i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }

        for (i = n - 2; i >= 0; i--) {
            if (isLMS(i)) {
                // cout << chr(i) << " at position " << i << " is LMS\n";
                // cout << "Inserting at bucket position " << bkt[rank[chr(i)]]
                //      << '\n';
                SA[bkt[rank[chr(i)]]--] = i;
                first = i;
            }
        }

        // SA[0] = n - 1; // set the single sentinel LMS-substring

        // Induce L-Type suffixes by using LMS-Type and L-Type suffixes
        induceSAl(t, SA, s, bkt, n, K, cs, level,rank);

        // Induce S-Type suffixes by using L-Type and S-Type suffixes
        induceSAs(t, SA, s, bkt, n, K, cs, level,rank);

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
            bool diff = false;
            int_t d;
            // d equals to the LCP between two consecutive LMS-substrings
            // for (d = 0; d < n; d++) {
            //     // If is first suffix in LMS order (sentinel), or one of the
            //     // suffixes reached the last position of T, or the
            //     // characters of T differs or the type os suffixes differ.
            //     if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||
            //         chr(pos + d) != chr(prev + d) ||
            //         (isLMS(pos + d) ^ (isLMS(prev + d)))) {
            //         diff = true;
            //         break;
            //     }
            //     // The comparison has reached the end of at least one
            //     // LMS-substring
            //     if (d > 0 && (isLMS(pos + d) || isLMS(prev + d))) {
            //         break;
            //     }
            // }
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
                lvl_stats.total_rule_len += cur_len;
                lvl_stats.total_lcp_len += d;
                lvl_stats.total_rule_suffix_len += cur_len - d;
                total_rule_len += cur_len;
                total_lcp += d;
                total_rule_suffix_length += cur_len - d;
#endif

                for (j = 0; j < cur_len - d && j + pos + d < n; j++) {
#ifdef REPORT
                    if (j + pos + d + 1 < n &&
                        chr(j + pos + d) == chr(j + pos + d + 1)) {
                        run_length_potential++;
                        lvl_stats.total_rl_potential++;
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
                lvl_stats.total_discarded_rules_len += len;
                lvl_stats.discarded_rules_n++;
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
        stats.level_n = stats.lvl_stats.size();
        lvl_stats.alphabet_size = K;
        lvl_stats.string_size = n;
        lvl_stats.rules_n = name + 1;
        lvl_stats.avg_rule_len =
            (double)lvl_stats.total_rule_len / lvl_stats.rules_n;
        lvl_stats.avg_discarded_rules_len =
            (double)lvl_stats.total_discarded_rules_len /
            lvl_stats.discarded_rules_n;
        lvl_stats.avg_rule_len =
            (double)lvl_stats.total_rule_len / lvl_stats.rules_n;
        lvl_stats.dictionary_size = g[level].size_in_bytes();
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
