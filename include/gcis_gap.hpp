#ifndef GCIS_GAP_HPP
#define GCIS_GAP_HPP

#include "gcis.hpp"
#include "gcis_gap_codec.hpp"

void print_text(sdsl::int_vector<> &v, size_t size) {
    for (size_t i = 0; i < size; i++) {
        cout << v[i] << " ";
    }
    cout << endl;
}
template <>
class gcis_dictionary<gcis_gap_codec> : public gcis_abstract<gcis_gap_codec> {

  public:
  private:
    std::vector<uint32_t> partial_sum;

  public:
    /**
     * @brief Serialize a gcis_gap_codec object
     *
     * @param o The ostream object in which the information will be stored.
     */
    void serialize(std::ostream &o) override {
        gcis_abstract::serialize(o);
        uint64_t n = partial_sum.size();
        o.write((char *)&n, sizeof(n));
        o.write((char *)partial_sum.data(), sizeof(uint32_t) * n);
    }

    void load(std::istream &i) override {
        gcis_abstract::load(i);
        uint64_t n;
        i.read((char *)&n, sizeof(n));
        partial_sum.resize(n);
        i.read((char *)partial_sum.data(), sizeof(uint32_t) * n);
    }

    /**
     * Extracts any valid substring T[l,r] from the text
     * @param l Beggining of such substring
     * @param r End of such substring
     * @return Returns the extracted substring
     */
    sdsl::int_vector<> extract(uint64_t l, uint64_t r) {
        uint64_t size =
            g.size() ? 4 * (g.back().fully_decoded_tail_len + (r - l + 1))
                     : (r - l + 1);
        sdsl::int_vector<> extracted_text(size);
        sdsl::int_vector<> tmp_text(size);
        extract(l, r, extracted_text, tmp_text);
        return extracted_text;
    }

    /**
     * @brief Extracts several valid substrings of the form T[l,r]
     * from the text.
     *
     * @param query A vector containing [l,r] pairs.
     */
    void extract_batch(vector<pair<int, int>> &query) {
        int l, r;
        std::tie(l, r) = query[0];
        uint64_t query_length = 50000;
        uint64_t size = query_length;
//            g.size() ?  (g.back().fully_decoded_tail_len + (query_length))
  //                   : (query_length);
        sdsl::int_vector<> extracted_text(size);
        sdsl::int_vector<> tmp_text(size);
        for (auto p : query) {
            cout << "Extracting"
                 << "[" << p.first << "," << p.second << "]" << endl;
            extract(p.first, p.second, extracted_text, tmp_text);
            for (uint64_t i = p.first; i <= p.second; i++) {
                cout << (char)extracted_text[i - p.first];
            }
            cout << endl;
        }
    }

    char *decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char *str;
        if (g.size()) {
            for (int64_t i = g.size() - 1; i >= 0; i--) {
                sdsl::int_vector<> next_r_string;
                gcis_gap_codec_level gd = std::move(g[i].decompress());
                next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
                next_r_string.resize(g[i].string_size);
                uint64_t l = 0;
                if (i == 0) {
                    // Convert the reduced string in the original text
                    str = new char[g[i].string_size];
                    for (uint64_t j = 0; j < g[i].tail.size(); j++) {
                        str[l++] = g[i].tail[j];
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
        } else {
            str = new char[reduced_string.size()];
            for (uint64_t i = 0; i < reduced_string.size(); i++) {
                str[i] = reduced_string[i];
            }
        }
        return str;
    }

    //     unsigned char *decode_saca(uint_t **sa) {

    //         sdsl::int_vector<> r_string = reduced_string;
    //         unsigned char *str;
    //         uint_t n = g[0].string_size;
    //         uint_t *SA = new uint_t[n];

    //         int_t *s = (int_t *)SA + n / 2;
    //         unsigned char *t =
    //             new unsigned char[n / 8 + 1]; // LS-type array in bits

    //         int cs = sizeof(int_t);
    //         if (g.size()) {

    //             for (int64_t level = g.size() - 1; level >= 0; level--) {

    // #if TIME
    //                 auto start = timer::now();
    // #endif

    //                 n = g[level].string_size;

    //                 uint_t n1 = r_string.size();
    //                 uint_t *SA1 = SA, *s1 = SA + n - n1;

    //                 // copy to s1[1]
    //                 if (level == g.size() - 1)
    //                     for (uint_t i = 0; i < n1; i++)
    //                         SA1[r_string[i]] = i;
    //                 else
    //                     for (uint_t i = 0; i < n1; i++)
    //                         s1[i] = SA[i];

    // #if DEBUG
    //                 cout << "\n##\nlevel = " << level << endl;
    //                 cout << "n = " << n << "\nn1 = " << n1 << endl;
    //                 cout << "\n####" << endl;
    //                 cout << "s1 = ";
    //                 for (uint_t i = 0; i < n1; i++) {
    //                     cout << s1[i] << " ";
    //                     cout << endl;
    //                 }
    // #endif

    // #if TIME
    //                 auto expand = timer::now();
    // #endif

    //                 sdsl::int_vector<> next_r_string;
    //                 gcis_gap_codec_level gd =
    //                     std::move(g[level].decompress());
    //                 next_r_string.width(sdsl::bits::hi(g[level].alphabet_size
    //                 - 1) +
    //                                     1);
    //                 next_r_string.resize(g[level].string_size);
    //                 uint64_t l = 0;

    //                 int_t K = g[level].alphabet_size; // alphabet

    //                 int *bkt = (int *)malloc(sizeof(int) * K);
    //                 init_buckets(bkt, K);

    //                 if (level == 0) {

    //                     //    delete[] s;
    //                     // Convert the reduced string in the original text
    //                     str = new unsigned char[g[level].string_size];
    //                     for (uint64_t j = 0; j < g[level].tail.size(); j++) {
    //                         str[l++] = g[level].tail[j];
    //                         bkt[g[level].tail[j]]++; // count frequencies
    //                     }
    //                     for (uint64_t j = 0; j < r_string.size(); j++) {
    //                         gd.expand_rule_bkt(r_string[j], str, l, bkt);
    //                     }
    //                     n = g[level].string_size;
    //                     // Classify the type of each character
    //                     tset(n - 2, 0);
    //                     tset(n - 1, 1); // the sentinel must be in s1,
    //                     important!!! for (int_t i = n - 3; i >= 0; i--) {
    //                         tset(i, (str[i] < str[i + 1] ||
    //                                  (str[i] == str[i + 1] && tget(i + 1) ==
    //                                  1))
    //                                     ? 1
    //                                     : 0);
    //                     }
    //                 } else {
    //                     // Convert the reduced string in the previous reduced
    //                     string for (uint64_t j = 0; j < g[level].tail.size();
    //                     j++) {
    //                         next_r_string[l++] = g[level].tail[j];
    //                         bkt[g[level].tail[j]]++; // count frequencies
    //                     }
    //                     for (uint64_t j = 0; j < r_string.size(); j++) {
    //                         gd.expand_rule_bkt(r_string[j], next_r_string, l,
    //                         bkt);
    //                     }
    //                     r_string = std::move(next_r_string);

    //                     // n=r_string.size();
    //                     n = g[level].string_size;
    //                     // copy to s[1]
    //                     for (uint_t i = 0; i < n; i++)
    //                         s[i] = r_string[i];

    //                     // Classify the type of each character
    //                     tset(n - 2, 0);
    //                     tset(n - 1, 1); // the sentinel must be in s1,
    //                     important!!! for (int_t i = n - 3; i >= 0; i--) {
    //                         tset(i, (r_string[i] < r_string[i + 1] ||
    //                                  (r_string[i] == r_string[i + 1] &&
    //                                   tget(i + 1) == 1))
    //                                     ? 1
    //                                     : 0);
    //                     }
    //                 }

    // #if TIME
    //                 auto end = timer::now();
    //                 cout << "expand: "
    //                      << (double)duration_cast<seconds>(end -
    //                      expand).count()
    //                      << " seconds" << endl;
    // #endif

    // #if DEGUB
    //                 cout << "n = " << n << "\nn1 = " << n1 << endl;
    //                 cout << "level = " << level
    //                      << "\t string_size = " << g[level].string_size <<
    //                      "\n"
    //                      << "alphabet_size = " << g[level].alphabet_size <<
    //                      endl;
    //                 cout << "SA1: ";
    //                 for (uint_t i = 0; i < n1; i++)
    //                     cout << SA[i] << ", ";
    //                 cout << endl;
    // #endif

    //                 // stage 3: induce the result for the original problem
    //                 get_buckets_end(bkt, true, K);

    // #if TIME
    //                 auto begin = timer::now();
    // #endif

    // #if TIME
    //                 end = timer::now();
    //                 cout << "classify: "
    //                      << (double)duration_cast<seconds>(end -
    //                      begin).count()
    //                      << " seconds" << endl;
    //                 begin = timer::now();
    // #endif

    // #if DEGUB
    //                 for (int_t i = 0; i < n; i++) {
    //                     cout << tget(i) << " ";
    //                 }
    //                 cout << endl;
    // #endif

    //                 int_t j = 0;
    //                 for (int_t i = 1; i < n; i++) {
    //                     if (isLMS(i)) {
    //                         s1[j++] = i; // get p1
    //                     }
    //                 }

    //                 for (int_t i = 0; i < n1; i++) {
    //                     SA1[i] = s1[SA1[i]]; // get index in s1
    //                 }
    //                 for (int_t i = n1; i < n; i++) {
    //                     SA[i] = EMPTY; // init SA[n1..n-1]
    //                 }

    //                 if (level) {
    //                     for (int_t i = n1 - 1; i >= 0; i--) {
    //                         j = SA[i];
    //                         SA[i] = EMPTY;
    //                         SA[bkt[chr(j)]--] = j;
    //                     }
    //                 } else {
    //                     for (int_t i = n1 - 1; i >= 0; i--) {
    //                         j = SA[i];
    //                         SA[i] = EMPTY;
    //                         if (i == 0)
    //                             SA[0] = n - 1;
    //                         else
    //                             SA[bkt[str[j]]--] = j;
    //                     }
    //                 }

    // #if TIME
    //                 end = timer::now();
    //                 cout << "position: "
    //                      << (double)duration_cast<seconds>(end -
    //                      begin).count()
    //                      << " seconds" << endl;
    //                 begin = timer::now();
    // #endif

    //                 if (level)
    //                     induceSAl(t, SA, s, bkt, n, K, cs, level);
    //                 else
    //                     induceSAl(t, SA, (int_t *)str, bkt, n, K,
    //                               sizeof(unsigned char), level);

    // #if TIME
    //                 end = timer::now();
    //                 cout << "induce L: "
    //                      << (double)duration_cast<seconds>(end -
    //                      begin).count()
    //                      << " seconds" << endl;
    //                 begin = timer::now();
    // #endif

    //                 if (level)
    //                     induceSAs(t, SA, s, bkt, n, K, cs, level);
    //                 else
    //                     induceSAs(t, SA, (int_t *)str, bkt, n, K,
    //                               sizeof(unsigned char), level);

    // #if TIME
    //                 end = timer::now();
    //                 cout << "induce S: "
    //                      << (double)duration_cast<seconds>(end -
    //                      begin).count()
    //                      << " seconds" << endl;
    // #endif

    // #if DEGUB
    //                 cout << "SA: ";
    //                 for (uint_t i = 0; i < n; i++) {
    //                     cout << SA[i] << ", ";
    //                 }
    //                 cout << endl;
    // #endif
    //                 free(bkt);

    // #if TIME
    //                 auto stop = timer::now();
    //                 cout << "time: "
    //                      << (double)duration_cast<seconds>(stop -
    //                      start).count()
    //                      << " seconds" << endl;
    // #endif
    //             }
    //         } else {
    //             str = new unsigned char[reduced_string.size()];
    //             for (uint64_t i = 0; i < reduced_string.size(); i++)
    //                 str[i] = reduced_string[i];
    //         }

    //         delete[] t;
    //         *sa = SA;
    //         return str;
    //     }

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

        unsigned char *t =
            new unsigned char[n / 8 + 1]; // LS-type array in bits
        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        //  tset(n - 2, 0);
        tset(n - 1, 1); // the sentinel must be in s1, important!!!
        for (i = n - 2; i >= 0; i--)
            tset(i, (chr(i) < chr(i + 1) ||
                     (chr(i) == chr(i + 1) && tget(i + 1) == 1))
                        ? 1
                        : 0);

        int_t *bkt = new int_t[K]; // bucket counters

        size_t first = n - 1;

        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets

        for (i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }

        for (i = n - 2; i >= 0; i--) {
            if (isLMS(i)) {
                SA[bkt[chr(i)]--] = i;
                first = i;
            }
        }

        SA[0] = n - 1; // set the single sentinel LMS-substring

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

        std::vector<uint64_t> fdrlen;

        uint_t rule_index = 0;
        vector<uint32_t> lcp;
        vector<uint32_t> rule_pos;
        uint64_t last_lcp = 0;
        uint64_t last_rule_pos = 0;
        g.push_back(gcis_gap_codec());
        // Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            bool diff = false;
            int_t d;
            // d equals to the LCP between two consecutive LMS-substrings
            for (d = 0; d < n; d++) {
                // If is first suffix in LMS order (sentinel), or one of the
                // suffixes reached the last position of T, or the
                // characters of T differs or the type os suffixes differ.
                if (prev == -1 ||
                    (chr(pos + d) != chr(prev + d) &&
                     (d == 0 || (!isLMS(pos + d) && !isLMS(prev + d)))) ||
                    (isLMS(pos + d) ^ (isLMS(prev + d)))) {
                    diff = true;
                    break;
                }
                // The comparison has reached the end of at least one
                // LMS-substring
                if (d > 0 && (isLMS(pos + d) || isLMS(prev + d))) {
                    break;
                }
            }

            // The consecutive LMS-substrings differs
            if (diff) {

                // TODO: can we replace len =1 for len	=d?
                // Get the length of the current lms-substring
                size_t len = 1;
                if (pos != n - 1)
                    while (!isLMS(pos + len))
                        len++;

                // Get the length of the previous LMS-Substring
                size_t len2 = 1;
                if (prev != n - 1)
                    while (!isLMS(prev + len2))
                        len2++;

                // Resizes Rule array, LCP array and rule delimiter bitvector
                uint64_t old_lcp_size, old_rule_delim_size;

                // Put a limit on how far Front Coding goes back
                // TODO: change magic number
                if (name % 32 == 0) {
                    d = 0;
                }

                lcp.push_back(d + last_lcp);
                last_lcp += d;

                rule_pos.push_back((len - d) + last_rule_pos);
                last_rule_pos += len - d;

                // cout << rule_pos.back() << endl;

                g[level].rule.resize(g[level].rule.size() + len - d);

#ifdef REPORT
                total_rule_len += len;
                total_lcp += d;
                total_rule_suffix_length += len - d;
#endif

                // Copy the symbols into the delimited rule positions
                for (j = 0; j < len - d && j + pos + d < n; j++) {
#ifdef REPORT
                    if (j + pos + d - 1 < n &&
                        chr(j + pos + d) == chr(j + pos + d + 1)) {
                        run_length_potential++;
                    }
#endif
                    g[level].rule[rule_index] = (uint_t)chr(j + pos + d);
                    rule_index++;
                }
                // Insert the fully decode rule length
                if (level == 0) {
                    // The symbols are terminal L(x) = 1, for every x
                    fdrlen.push_back(len);
                } else {
                    // The symbols are not necessarily terminal.
                    uint64_t sum =
                        g[level - 1].fully_decoded_rule_len[chr(pos)];
                    for (uint64_t i = 1; i + pos < n && !isLMS(pos + i); i++) {
                        sum +=
                            g[level - 1].fully_decoded_rule_len[chr(pos + i)];
                    }
                    fdrlen.push_back(sum);
                }
                name++;
                prev = pos;
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
        g[level].lcp =
            std::move(sdsl::enc_vector<sdsl::coder::elias_delta>(lcp));
        g[level].rule_pos =
            std::move(sdsl::enc_vector<sdsl::coder::elias_delta>(rule_pos));
        g[level].fully_decoded_rule_len = std::move(sdsl::dac_vector<>(fdrlen));

        lcp.clear();
        rule_pos.clear();

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
        // Compress the tail
        sdsl::util::bit_compress(g[level].tail);
        // Determine the accumulated tail length
        if (level > 0) {
            g[level].fully_decoded_tail_len =
                g[level - 1].fully_decoded_tail_len;
            for (auto t : g[level].tail) {
                g[level].fully_decoded_tail_len +=
                    g[level - 1].fully_decoded_rule_len[t];
            }
        } else {
            g[level].fully_decoded_tail_len = g[level].tail.size();
        }

        // stage 2: solve the reduced problem
        // recurse if names are not yet unique

#ifdef REPORT
        print_report("Level ", level, "\n");
        print_report("Alphabet Size = ", K, "\n");
        print_report("String Size = ", n, "\n");
        print_report("Number of Rules = ", name + 1, "\n");
        print_report("Average Rule Length = ",
                     (double)total_rule_len / (name + 1), "\n");
        print_report("Number of Discarded Rules = ", discarded_rules_n, "\n");
        print_report("Average Discarded Rules Length = ",
                     (double)discarded_rules_len / discarded_rules_n, "\n");
        print_report("Average LCP = ", (double)total_lcp / (name + 1), "\n");
        print_report("Average Rule Suffix Length = ",
                     (double)total_rule_suffix_length / (name + 1), "\n");
        print_report(
            "Dictionary Level Size (bytes) =", g[level].size_in_bytes(), "\n");
        print_report("LCP Size (bits) = ", g[level].lcp.size(), "\n");
        print_report("Rule Suffix Length (total) = ", g[level].rule.size(),
                     "\n");
        print_report("Rule Suffix Width (bits per symbol) = ",
                     (int_t)g[level].rule.width(), "\n");
        print_report("Tail Length = ", g[level].tail.size(), "\n");
        print_report("Tail Width (bits per symbol) = ",
                     (int_t)g[level].tail.width(), "\n");
        print_report("Run Length Potential (total) = ", run_length_potential,
                     "\n");
        print_report("Avg Run Length per Rule Suffix = ",
                     (double)run_length_potential / (name + 1), "\n");
#endif

        // TODO: reenable premature_stop comparison
        bool premature_stop = false;
        //	  bool premature_stop =
        // evaluate_premature_stop(n,K,n1,name+1,level);
        g[level].string_size = n;
        g[level].alphabet_size = K;
        if (name + 1 < n1 && !premature_stop) {
            // cout << "rules = ";
            // for (uint64_t i = 0; i < n1; i++) {
            //     cout << (int)s1[i];
            // }
            // cout << endl;
            gc_is((int_t *)s1, SA1, n1, name + 1, sizeof(int_t), level + 1);
        } else {
            // generate the suffix array of s1 directly
            if (premature_stop) {
                // The encoding algorithm is stopped prematurely
#ifdef REPORT
                print_report("Premature Stop employed at level ", level, "\n");
#endif
                reduced_string.resize(n);
                // cout << "Reduced string = ";
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] =
                        (uint64_t)(cs == sizeof(char) ? ((char *)s)[j] : s[j]);
                    // cout << reduced_string[j];
                }
                // cout << endl;
                // The last level of computation is discarded
                g.pop_back();
            } else {
                // cout << "Reduced string = ";
                reduced_string.resize(n1);
                for (j = 0; j < n1; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s1[j];
                    // cout << reduced_string[j];
                }
                // cout << endl;
            }
            sdsl::util::bit_compress(reduced_string);
            partial_sum.resize(reduced_string.size());
            partial_sum[0] = 0;
            for (uint64_t i = 1; i < reduced_string.size(); i++) {
                partial_sum[i] =
                    partial_sum[i - 1] +
                    g.back().fully_decoded_rule_len[reduced_string[i - 1]];
                // cout << "Partial sum = " << partial_sum[i] << "\n";
            }

#ifdef REPORT
            print_report(
                "Reduced String Length = ", (int_t)reduced_string.size(), "\n");
            print_report("Reduced String Width (bits per symbol) = ",
                         (int_t)reduced_string.width(), "\n");
#endif
        }
        delete[] t;
    }

    /**
     *
     * @param partial_sum The partial sum data-structure
     * @param sz the position we want to find
     * @return the leftmost index rk such that partial_sum[rk]>=sz
     */
    uint64_t bsearch_upperbound(vector<uint32_t> &partial_sum, uint64_t sz) {
        if (partial_sum.back() <= sz) {
            return partial_sum.size() - 1;
        }
        int64_t l = 0, r = partial_sum.size() - 1;
        while (l < r) {
            int64_t mid = l + (r - l) / 2;
            if (partial_sum[mid] +
                    g.back().fully_decoded_rule_len[reduced_string[mid]] >
                sz) {
                r = mid;
            } else {
                l = mid + 1;
            }
        }
        return l;
    }

    /**
     *
     * @param partial_sum The partial sum data-structure
     * @param sz the position we want to find
     * @return the rightmost index lk such that partial_sum[lk]<=sz
     */
    uint64_t bsearch_lowerbound(vector<uint32_t> &partial_sum, uint64_t sz) {
        if (partial_sum[0] >= sz) {
            return 0;
        }
        int64_t l = 0, r = partial_sum.size() - 1;
        while (l < r) {
            int64_t mid = l + (r - l + 1) / 2;
            if (partial_sum[mid] <= sz) {
                l = mid;
            } else {
                r = mid - 1;
            }
        }
        return l;
    }

    /**
     *
     * @param codec The dictionary level we are employing
     * @param extracted_text The string to be extracted
     * @param sz The index we want to find
     * @param text_r The tracked position of the original text representing the
     * end of extracted_text
     * @return The leftmost index such that text_r sum(extracted_text,i) <= sz
     */
    uint64_t sequential_upperbound(gcis_gap_codec &codec,
                                   sdsl::int_vector<> &extracted_text,
                                   int64_t extracted_text_size, int64_t sz,
                                   int64_t &text_r) {
        int64_t index = 0;
        int64_t rule_length = 0;
        for (index = extracted_text_size - 1; index >= 0; index--) {
            rule_length = codec.fully_decoded_rule_len[extracted_text[index]];
            if (text_r - rule_length <= sz) {
                break;
            }
            text_r -= rule_length;
        }
        return index;
    }

    /**
     *
     * @param codec The dictionary level we are employing
     * @param extracted_text The string to be extracted
     * @param sz The index we want to found
     * @param text_l The tracked position of the original text representing the
     * beggining of extracted_text
     * @return The rightmost index such that text_l + sum(extracted_text,i) >=
     * sz
     */
    uint64_t sequential_lowerbound(gcis_gap_codec &codec,
                                   sdsl::int_vector<> &extracted_text,
                                   int64_t extracted_text_size, int64_t sz,
                                   int64_t &text_l) {
        int64_t index = 0;
        int64_t rule_length = 0;
        for (index = 0; index < extracted_text_size; index++) {
            rule_length = codec.fully_decoded_rule_len[extracted_text[index]];
            if (text_l + rule_length > sz) {
                break;
            }
            text_l += rule_length;
        }
        return index;
    }

    /**
     * @brief Extract the substring T[l,r] from the text.
     *
     * @param l Beggining of the substring
     * @param r End of the substring
     * @param extracted_text Extracted substring buffer
     * @param tmp_text Temporary Buffer
     */
    void extract(int64_t l, int64_t r, sdsl::int_vector<> &extracted_text,
                 sdsl::int_vector<> &tmp_text) {
        // Stores the interval being tracked in the text
        int64_t text_l;
        int64_t text_r;
        // Stores the interval being tracked in the level
        uint64_t lk, rk;
        text_l = 0;
        text_r = g.size() > 0 ? g[0].string_size : reduced_string.size();
        uint64_t extracted_idx = 0;

        /**
         * The dictionary is only the original string
         * The extraction is done in a straight-forward fashion
         */
        if (g.size() == 0) {
            for (uint64_t j = 0; j < reduced_string.size(); j++) {
                extracted_text[j] = reduced_string[j];
            }
            return;
        }

        /**
         * Else, the extraction should proceed from the reduced string to the
         * decompressed text
         */

        if (r < g.back().fully_decoded_tail_len) {
            // The string lies on the tail. Copy all the tail.
            text_l = 0;
            text_r = g.back().fully_decoded_tail_len - 1;
            for (auto v : g.back().tail) {
                tmp_text[extracted_idx++] = v;
            }
        } else if (l < g.back().fully_decoded_tail_len) {
            // A prefix of the string lies on the tail.
            // Copy the tail
            text_l = 0;
            for (auto v : g.back().tail) {
                cout << v;
                tmp_text[extracted_idx++] = v;
            }
            cout << endl;
            // Find the leftmost index which covers r
            rk = bsearch_upperbound(partial_sum,
                                    r - g.back().fully_decoded_tail_len);
            text_r = g.back().fully_decoded_tail_len + partial_sum[rk] +
                     g.back().fully_decoded_rule_len[reduced_string[rk]];
            // Decompress the rules located at reduced_string[0..rk];
            // print_text(reduced_string, rk + 1);
            for (uint64_t i = 0; i <= rk; i++) {
                g.back().extract_rule(reduced_string[i], tmp_text,
                                      extracted_idx);
            }
        } else {
            // The string does not occur in the tail
            // Find the rightmost index which covers l
            lk = bsearch_lowerbound(partial_sum,
                                    l - g.back().fully_decoded_tail_len);
            text_l = g.back().fully_decoded_tail_len + partial_sum[lk];
            // Find the leftmost index which covers r
            rk = bsearch_upperbound(partial_sum,
                                    r - g.back().fully_decoded_tail_len);
            text_r = g.back().fully_decoded_tail_len + partial_sum[rk] +
                     g.back().fully_decoded_rule_len[reduced_string[rk]];

            // Decompress the rules located at reduced_string[0..rk];
            for (uint64_t i = lk; i <= rk; i++) {
                g.back().extract_rule(reduced_string[i], tmp_text,
                                      extracted_idx);
            }
        }
        int64_t level = g.size() - 2;
        // Extract the reduced string part
        while (level >= 0) {
            uint64_t extracted_text_len = extracted_idx;
            extracted_idx = 0;
            std::swap(extracted_text, tmp_text);
            // The extracted string lies on the tail
            if (r < g[level].fully_decoded_tail_len) {
                text_l = 0;
                text_r = g[level].fully_decoded_tail_len - 1;
                // Copy the tail
                if (text_r > tmp_text.size()) {
                    tmp_text.resize(text_r);
                    extracted_text.resize(text_r);
                }
                for (auto v : g[level].tail) {
                    tmp_text[extracted_idx++] = v;
                }
            } else if (l < g[level].fully_decoded_tail_len) {
                // A prefix of the string lies on the tail
                text_l = 0;
                // Copy the tail
                for (auto v : g[level].tail) {
                    tmp_text[extracted_idx++] = v;
                }
                rk = sequential_upperbound(g[level], extracted_text,
                                           extracted_text_len, r, text_r);
                // print_text(extracted_text, rk + 1);
                for (uint64_t i = 0; i <= rk; i++) {
                    g[level].extract_rule(extracted_text[i], tmp_text,
                                          extracted_idx);
                }
            } else {
                text_l =
                    std::max<int64_t>(text_l, g[level].fully_decoded_tail_len);
                lk = sequential_lowerbound(g[level], extracted_text,
                                           extracted_text_len, l, text_l);
                rk = sequential_upperbound(g[level], extracted_text,
                                           extracted_text_len, r, text_r);
                // print_text(extracted_text, rk - lk + 1);
                for (uint64_t i = lk; i <= rk; i++) {
                    g[level].extract_rule(extracted_text[i], tmp_text,
                                          extracted_idx);
                }
                // print_text(tmp_text,extracted_idx);
            }
            level--;
        }

        for (uint64_t i = 0, offset = l - text_l; i < r - l + 1; i++) {
            extracted_text[i] = tmp_text[i + offset];
        }
        // print_text(extracted_text, r - l + 1);
    }
};

#endif