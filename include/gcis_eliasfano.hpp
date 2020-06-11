//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_ELIASFANO_HPP
#define GC_IS_GCIS_ELIASFANO_HPP

#include "gcis.hpp"
#include "gcis_eliasfano_codec.hpp"
#include <iostream>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

#define TIME 0

template <>
class gcis_dictionary<gcis_eliasfano_codec>
    : public gcis_abstract<gcis_eliasfano_codec> {

  public:
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
        //            g.size() ?  (g.back().fully_decoded_tail_len +
        //            (query_length))
        //                   : (query_length);
        sdsl::int_vector<> extracted_text(size);
        sdsl::int_vector<> tmp_text(size);
        auto first = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::high_resolution_clock::now();
        for (auto p : query) {
            // cout << "Extracting"
            //      << "[" << p.first << "," << p.second << "]" << endl;
            auto t0 = std::chrono::high_resolution_clock::now();
            extract(p.first, p.second, extracted_text, tmp_text);
            auto t1 = std::chrono::high_resolution_clock::now();
            total_time += t1 - t0;
            for (uint64_t i = p.first; i <= p.second; i++) {
                cout << (unsigned char)extracted_text[i - p.first];
            }
            cout << endl;
        }
        std::chrono::duration<double> elapsed = total_time - first;
        cout << "Batch Extraction Total time(s): " << elapsed.count() << endl;
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

    pair<char *, int_t> decode() override {
        vector<uint_t> r_string(reduced_string.size());
        for (uint_t i = 0; i < reduced_string.size(); i++) {
            r_string[i] = reduced_string[i];
        }
        char *str;
        if (g.size()) {
            for (int64_t i = g.size() - 1; i >= 0; i--) {
                vector<uint_t> next_r_string;
                gcis_eliasfano_pointers_codec_level gd =
                    std::move(g[i].decompress());
                next_r_string.resize(g[i].string_size);
                uint_t l = 0;
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
        return make_pair(str, g[0].string_size);
    }

    pair<char *, int_t> decode_saca(uint_t **sa) override {

        std::vector<uint_t> r_string;
        r_string.resize(reduced_string.size());
        for (int i = 0; i < reduced_string.size(); i++) {
            r_string[i] = reduced_string[i];
        }
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

                vector<uint_t> next_r_string;
                gcis_eliasfano_pointers_codec_level gd =
                    std::move(g[level].decompress());
                next_r_string.resize(g[level].string_size);
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
                        if (cur_t == 0 && succ_t == 1) {
                            // cout << "s1[" << j << "] = " << i << endl;
                            s1[j--] = i;
                        }
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
    } // end decode_saca

    pair<char *,int_t> decode_saca_lcp(uint_t **sa, int_t **lcp) override {

        vector<uint_t> r_string;
        r_string.resize(reduced_string.size());
        for (uint_t i = 0; i < reduced_string.size(); i++)
            r_string[i] = reduced_string[i];
        unsigned char *str;
        uint_t n = g[0].string_size;
        uint_t *SA = new uint_t[n];
        int_t *LCP = new int_t[n];

        uint_t i;
        for (i = 0; i < n; i++)
            SA[i] = LCP[i] = 0;

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

                vector<uint_t> next_r_string;
                gcis_eliasfano_pointers_codec_level gd =
                    std::move(g[level].decompress());
                next_r_string.resize(g[level].string_size);
                uint_t l = 0;

                int_t K = g[level].alphabet_size; // alphabet

                int_t *bkt = new int_t[K]; // bucket
                int_t *cnt = new int_t[K]; // counters

                init_buckets(cnt, K);

                if (level == 0) {

                    // delete[] s;
                    // Convert the reduced string in the original text
                    str = new unsigned char[g[level].string_size];
                    for (uint_t j = 0; j < g[level].tail.size(); j++) {
                        str[l++] = g[level].tail[j];
                        cnt[g[level].tail[j]]++; // count frequencies
                    }
                    for (uint_t j = 0; j < r_string.size(); j++) {
                        gd.expand_rule_bkt(r_string[j], str, l, cnt);
                    }
                    n = g[level].string_size;
                    // str[n - 1] = 0;
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
                    for (uint64_t j = 0; j < g[level].tail.size(); j++) {
                        next_r_string[l] = g[level].tail[j];
                        cnt[next_r_string[l++]]++; // count frequencies
                    }
                    for (uint64_t j = 0; j < r_string.size(); j++) {
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
                    succ_t = 0; // s[n-2] must be L-type
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

                if (level == 0) {
                    uint_t *RA = s1;
                    int_t *PLCP = LCP + n - n1; // PHI is stored in PLCP array
                    // compute the LCP of consecutive LMS-suffixes
                    compute_lcp_phi_sparse_sais((int_t *)str, SA1, RA, LCP,
                                                PLCP, n1, sizeof(char),n);
                }

                int_t j = 0;
                for (int_t i = 0; i < n1; i++) {
                    SA1[i] = s1[SA1[i]]; // get index in s1
                }
                for (int_t i = n1; i < n; i++) {
                    SA[i] = EMPTY; // init SA[n1..n-1]
                }

                if (level == 0) {
                    for (i = n1; i < n; i++)
                        LCP[i] = 0;
                }

                if (level) {
                    for (int_t i = n1 - 1; i >= 0; i--) {
                        j = SA[i];
                        SA[i] = EMPTY;
                        SA[bkt[chr(j)]--] = j;
                    }
                } else {
                    int_t l;
                    for (int_t i = n1 - 1; i > 0; i--) {
                        j = SA[i];
                        SA[i] = U_MAX;
                        l = LCP[i];
                        LCP[i] = 0;

                        SA[bkt[str[j]]] = j;
                        LCP[bkt[str[j]]--] = l;
                    }
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
                    // induceSAl(SA, (int_t *)str, cnt, bkt, n, K,
                    //          sizeof(char), level);
                    induceSAl_LCP(SA, LCP, (int_t *)str, cnt, bkt, n, K,
                                  sizeof(unsigned char), level);

#if DEGUB
                if (level == 0) {
                    for (i = 0; i < n; i++)
                        printf("%d\t", SA[i]);
                    printf("\n\n");
                    for (i = 0; i < n; i++)
                        printf("%d\t", LCP[i]);
                    printf("\n");
                }
#endif

#if TIME
                end = timer::now();
                cout << "induce L: "
                     << (double)duration_cast<seconds>(end - begin).count()
                     << " seconds" << endl;
                begin = timer::now();
#endif

                if (level)
                    induceSAs(SA, s, cnt, bkt, n, K, cs, level);
                else {
                    // induceSAs(SA, (int_t *)str, cnt, bkt, n, K,
                    //          sizeof(char), level);
                    induceSAs_LCP(SA, LCP, (int_t *)str, cnt, bkt, n, K,
                                  sizeof(unsigned char), level);
                    SA[0] = n - 1;
                }

#if DEGUB
                if (level == 0) {
                    for (i = 0; i < n; i++)
                        printf("%d\t", SA[i]);
                    printf("\n\n");
                    for (i = 0; i < n; i++)
                        printf("%d\t", LCP[i]);
                    printf("\n");
                }
#endif
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
        *lcp = LCP;
        return make_pair((char*) str,g[0].string_size);
    } // end decode_saca

  private:
    std::vector<uint32_t> partial_sum;

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
        tset(n - 1, 0); // The last symbol is L-type
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
        int_t cur_len = -1;
        int_t prev_len = -1;

        std::vector<uint64_t> fdrlen;

        uint_t rule_index = 0;
        sdsl::bit_vector lcp;
        sdsl::bit_vector rule_delim;
        g.push_back(gcis_eliasfano_codec());
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

                // Resizes Rule array, LCP array and rule delimiter
                // bitvector
                uint64_t old_lcp_size, old_rule_delim_size;

                old_lcp_size = lcp.size();
                old_rule_delim_size = rule_delim.size();

                // Put a limit on how far Front Coding goes back
                // TODO: change magic number
                if (name % 32 == 0) {
                    d = 0;
                }

                g[level].rule.resize(g[level].rule.size() + cur_len - d);
                lcp.resize(lcp.size() + d + 1);
                rule_delim.resize(rule_delim.size() + cur_len - d + 1);

                for (uint64_t i = old_lcp_size; i < lcp.size(); i++) {
                    lcp[i] = 0;
                }
                for (uint64_t i = old_rule_delim_size; i < rule_delim.size();
                     i++) {
                    rule_delim[i] = 0;
                }

#ifdef REPORT
                total_rule_len += cur_len;
                total_lcp += d;
                total_rule_suffix_length += cur_len - d;
#endif

                // Encode LCP value in unary
                lcp[lcp.size() - 1] = 1;
                // Encode rule length in unary
                rule_delim[rule_delim.size() - 1] = 1;
                // Copy the symbols into the delimited rule positions
                for (j = 0; j < cur_len - d && j + pos + d < n; j++) {
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
                    fdrlen.push_back(cur_len);
                } else {
                    // The symbols are not necessarly terminal.
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
        g[level].lcp.encode(lcp);
        g[level].rule_suffix_length.encode(rule_delim);
        g[level].fully_decoded_rule_len = sdsl::dac_vector_dp<>(fdrlen);

        sdsl::util::clear(lcp);
        sdsl::util::clear(rule_delim);

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
        gcis::util::print_report("Level ", level, "\n");
        gcis::util::print_report("Alphabet Size = ", K, "\n");
        gcis::util::print_report("String Size = ", n, "\n");
        gcis::util::print_report("Number of Rules = ", name + 1, "\n");
        gcis::util::print_report("Total rules length  = ", total_rule_len,
                                 "\n");
        gcis::util::print_report("Average Rule Length = ",
                                 (double)total_rule_len / (name + 1), "\n");
        gcis::util::print_report(
            "Number of Discarded Rules = ", discarded_rules_n, "\n");
        gcis::util::print_report("Average Discarded Rules Length = ",
                                 discarded_rules_n > 0
                                     ? (double)discarded_rules_len /
                                           discarded_rules_n
                                     : 0,
                                 "\n");
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

        // TODO: reenable premature_stop comparison
        bool premature_stop = false;
        //	  bool premature_stop =
        // evaluate_premature_stop(n,K,n1,name+1,level);
        g[level].string_size = n;
        g[level].alphabet_size = K;
        if (name + 1 < n1 && !premature_stop) {
            gc_is((int_t *)s1, SA1, n1, name + 1, sizeof(int_t), level + 1);
        } else {
            // generate the suffix array of s1 directly
            if (premature_stop) {
                // The encoding algorithm is stopped prematurely
#ifdef REPORT
                gcis::util::print_report("Premature Stop employed at level ",
                                         level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] =
                        (uint64_t)(cs == sizeof(char) ? ((char *)s)[j] : s[j]);
                }
                // The last level of computation is discarded
                g.pop_back();
            } else {
                reduced_string.resize(n1);
                for (j = 0; j < n1; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s1[j];
                }
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
            gcis::util::print_report(
                "Reduced String Length = ", (int_t)reduced_string.size(), "\n");
            gcis::util::print_report(
                "Reduced String Width (bits per symbol) = ",
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
     * @param text_r The tracked position of the original text representing
     * the end of extracted_text
     * @return The leftmost index such that text_l r sum(extracted_text,i)
     * <= sz
     */
    uint64_t sequential_upperbound(gcis_eliasfano_codec &codec,
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
     * @param text_l The tracked position of the original text representing
     * the beggining of extracted_text
     * @return The rightmost index such that text_l + sum(extracted_text,i)
     * > sz
     */
    uint64_t sequential_lowerbound(gcis_eliasfano_codec &codec,
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

    void extract(int64_t l, int64_t r, sdsl::int_vector<> &extracted_text,
                 sdsl::int_vector<> &tmp_text) {
        //	  // Stores the interval being tracked in the text
        int64_t text_l;
        int64_t text_r;
        // Stores the interval being tracked in the level
        uint64_t lk, rk;
        text_l = 0;
        text_r = g.size() > 0 ? g[0].string_size : reduced_string.size();
        uint64_t extracted_idx = 0;

        /**
         * The dictionary is only the original string
         * The extraction is done in a straight-foward fashion
         */
        if (g.size() == 0) {
            for (uint64_t j = 0; j < reduced_string.size(); j++) {
                extracted_text[j] = reduced_string[j];
            }
            return;
        }

        /**
         * Else, the extraction should proceed from the reduced string to
         * the decompressed text
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
                tmp_text[extracted_idx++] = v;
            }
            // Find the leftmost index which covers r
            rk = bsearch_upperbound(partial_sum,
                                    r - g.back().fully_decoded_tail_len);
            text_r = g.back().fully_decoded_tail_len + partial_sum[rk] +
                     g.back().fully_decoded_rule_len[reduced_string[rk]];
            // Decompress the rules located at reduced_string[0..rk];
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
                for (uint64_t i = lk; i <= rk; i++) {
                    g[level].extract_rule(extracted_text[i], tmp_text,
                                          extracted_idx);
                }
            }
            level--;
        }

        for (uint64_t i = 0, offset = l - text_l; i < r - l + 1; i++) {
            extracted_text[i] = tmp_text[i + offset];
        }
    }
};

#endif // GC_IS_GCIS_ELIASFANO_HPP
