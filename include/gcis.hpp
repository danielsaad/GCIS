#ifndef GC_IS_HPP
#define GC_IS_HPP

#include "gcis_eliasfano_codec.hpp"
#include "gcis_s8b_codec.hpp"
#include "gcis_unary_codec.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "util.hpp"
#include <cstdint>
#include <cstring>
#define chr(i) (cs == sizeof(int_t) ? ((int_t *)s)[i] : ((unsigned char *)s)[i])

#define false 0
#define true 1

#define DEBUG 0
#define DEPTH                                                                  \
    1 // compute time and size of reduced problem for each recursion call
#define PHASES 0 // compute time for each phase

#ifdef m64
const int_t EMPTY = 0xffffffffffffffff;
#else
const int EMPTY = 0xffffffff;
#endif

unsigned char mask[] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

#define tget(i) ((t[(i) / 8] & mask[(i) % 8]) ? 1 : 0)
#define tset(i, b)                                                             \
    t[(i) / 8] =                                                               \
        (b) ? (mask[(i) % 8] | t[(i) / 8]) : ((~mask[(i) % 8]) & t[(i) / 8])

#define isLMS(i) (i > 0 && tget(i) && !tget(i - 1))

template <class codec_t> class gcis_abstract {
  public:
    std::vector<codec_t> g;
    sdsl::int_vector<> reduced_string;

  public:
    virtual uint64_t size_in_bytes() {
        uint64_t total_bytes = 0;
        for (uint64_t i = 0; i < g.size(); i++) {
            total_bytes += g[i].size_in_bytes();
        }
        total_bytes += sdsl::size_in_bytes(reduced_string);
        return total_bytes;
    }

    void encode(char *s) {
        int_t n = strlen(s) + 1;
        uint_t *SA = new uint_t[n];
        int_t K = 256;
        int cs = sizeof(char);
        int level = 0;

        gc_is((int_t *)s, SA, n, K, cs, level);

        delete[] SA;
    }

    virtual char *decode() = 0;

    virtual void serialize(std::ostream &o) {
        reduced_string.serialize(o);
        uint64_t size = g.size();
        o.write((char *)&size, sizeof(uint64_t));
        for (uint64_t i = 0; i < g.size(); i++) {
            g[i].serialize(o);
        }
    }

    virtual void load(std::istream &i) {
        uint64_t size;
        reduced_string.load(i);
        i.read((char *)&size, sizeof(uint64_t));
        g.resize(size);
        for (uint64_t j = 0; j < size; j++) {
            g[j].load(i);
        }
    }

    bool suffix_array_check(uint_t *SA, unsigned char *s, size_t len, int cs,
                            unsigned char sentinel) {
        int_t i, j, k;

        for (i = 0; i < len - 1; i++) {
            size_t min = SA[i + 1] < SA[i] ? (len - SA[i]) : (len - SA[i + 1]);
            if (!sleq(s, SA[i], SA[i + 1], min, cs, sentinel)) {

                printf("#%d) %d, %d&\n", i, SA[i], SA[i + 1]);

                for (j = SA[i], k = SA[i + 1]; (j < SA[i] + 5); j++, k++)
                    printf("%d | %d\n", chr(j), chr(k));
                printf("\n");

                return 0;
            }
        }

        unsigned char *tmp =
            (unsigned char *)malloc(len * sizeof(unsigned char));

        for (i = 0; i < len; i++)
            tmp[i] = 0;

        for (i = 0; i < len; i++)
            tmp[SA[i]] = 1;

        for (i = 0; i < len; i++) {
            if (!tmp[i]) {
                free(tmp);
                fprintf(stderr, "Error: not a permutation\n");
                return 0;
            }
        }

        free(tmp);

        return 1;
    }

    int suffix_array_write(uint_t *SA, uint_t n, char *c_file,
                           const char *ext) {

        FILE *f_out;
        char *c_out =
            (char *)malloc((strlen(c_file) + strlen(ext) + 3) * sizeof(char));

        sprintf(c_out, "%s.%s", c_file, ext);
        f_out = fopen(c_out, "wb");

        fwrite(SA, sizeof(uint_t), n, f_out);

        fclose(f_out);
        free(c_out);

        return 1;
    }

    uint_t *saca(char *s, uint_t *SA, int_t n) {

        int_t K = 256;
        int cs = sizeof(char);
        int level = 0;

        SAIS((int_t *)s, SA, n, K, cs, level);

        return SA;
    }

  protected:
    bool evaluate_premature_stop(int_t n0, int_t alphabet_size_s0, int_t n1,
                                 int_t alphabet_size_s1, int_t level) {

        uint64_t number_of_bits_s0 = sdsl::bits::hi(alphabet_size_s0) + 1;
        uint64_t number_of_bits_s1 = sdsl::bits::hi(alphabet_size_s1) + 1;
        uint64_t dict_size_bits = g[level].size_in_bytes() * 8;
        if (n0 * number_of_bits_s0 < n1 * number_of_bits_s1 + dict_size_bits)
            return true;
        return false;
    }

    //!
    //! \param s Sequence containing \0 in the end
    //! \param SA i-th level Suffix Array
    //! \param n size of sequence s
    //! \param K  Alphabet Size
    //! \param cs Symbol Length
    //! \param level current level of recursion.

    virtual void gc_is(int_t *s, uint_t *SA, int_t n, int_t K, int cs,
                       int level) {
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

        // Compact all the sorted substrings into the first n1 items of s
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

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        //        g.push_back(gcis_unary_codec());
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
                if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||
                    chr(pos + d) != chr(prev + d) ||
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

#ifdef REPORT
                total_rule_len += len;
                total_lcp += d;
                total_rule_suffix_length += len - d;
#endif

                for (j = 0; j < len - d && j + pos + d < n; j++) {
#ifdef REPORT
                    if (j + pos + d - 1 < n &&
                        chr(j + pos + d) == chr(j + pos + d + 1)) {
                        run_length_potential++;
                    }
#endif
                }
                // Since the adjacent LMS substrings differ, we must assign
                // a new name
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
            // Insert the name in upper half of the SA array
            // just after the LMS substring positions
            SA[n1 + pos] = name;
        }

        // Compact the n1 renamed substrings in the end of SA
        for (i = n - 1, j = n - 1; i >= n1; i--) {
            if (SA[i] != EMPTY) {
                SA[j--] = SA[i];
            }
        }

        // We use the same space of SA to store SA1 and s1
        // SA1 contains the space necessary to suffix sort the renamed string
        // s1 contains the renamed string
        uint_t *SA1 = SA;
        uint_t *s1 = SA + n - n1;

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
        //        print_report("LCP Size (bits) = ",g[level].lcp.size(),"\n");
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

        bool premature_stop =
            evaluate_premature_stop(n, K, n1, name + 1, level);
        if (name + 1 < n1 && !premature_stop) {
            gc_is((int_t *)s1, SA1, n1, name + 1, sizeof(int_t), level + 1);
        } else { // generate the suffix array of s1 directly
            if (premature_stop) {
#ifdef REPORT
                print_report("Premature Stop employed at level ", level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s[j];
                }
                // Discard the last computed level
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
            print_report(
                "Reduced String Length = ", (int_t)reduced_string.size(), "\n");
            print_report("Reduced String Width (bits per symbol) = ",
                         (int_t)reduced_string.width(), "\n");
#endif
        }
        // delete bitvector t
        delete[] t;
    }

    sdsl::int_vector<> gd_is(int level) {
        if (level == g.size()) {
            return (reduced_string);
        } else {
            // Picks the reduced string from the next level
            sdsl::int_vector<> r_string = gd_is(level + 1);
            // TODO: Compress next_r_string beforehand
            sdsl::int_vector<> next_r_string;
            // Decode the reduced string using the level grammar
            next_r_string.resize(g[level].tail.size());

            // Copy the tail to the beggining of the new reduced string
            for (uint64_t i = 0; i < g[level].tail.size(); i++) {
                next_r_string[i] = g[level].tail[i];
            }

            // Expand the rules to each position of the new reduced string
            for (uint64_t i = 0; i < r_string.size(); i++) {
                uint64_t rule = r_string[i];
                g[level].expand_rule(rule, next_r_string);
            }
            // returns the new reduced string to the level i-1
            sdsl::util::bit_compress(next_r_string);
            return next_r_string;
        }
    }

    // compute SA for the S-Type suffixes by inducing the L-Type suffixes and
    // the S-Type suffixes
    void induceSAs(unsigned char *t, uint_t *SA, int_t *s, int_t *bkt, int_t n,
                   int_t K, int cs, int level) {
        int_t i, j;
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets
        for (i = n - 1; i >= 0; i--) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0 && tget(j)) {
                    SA[bkt[chr(j)]--] = j;
                }
            }
        }
    }

    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the
    // L-Suffixes
    void induceSAl(unsigned char *t, uint_t *SA, int_t *s, int_t *bkt, int_t n,
                   int_t K, int cs, int level) {
        int_t i, j;
        // find heads of buckets
        get_buckets(s, bkt, n, K, cs, false);
        //  if(level==0) bkt[0]++;
        for (i = 0; i < n; i++) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0 && !tget(j)) {
                    SA[bkt[chr(j)]++] = j;
                }
            }
        }
    }

    // compute SA for the S-Type suffixes by inducing the L-Type suffixes and
    // the S-Type suffixes
    void induceSAs(uint_t *SA, int_t *s, int_t *cnt, int_t *bkt, int_t n,
                   int_t K, int cs, int level) {
        int_t i, j;
        get_buckets(cnt, bkt, K, true);
        for (i = n - 1; i >= 0; i--) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0)
                    if (chr(j) <= chr(j + 1) && bkt[chr(j)] < i) {
                        SA[bkt[chr(j)]--] = j;
                    }
            }
        }
    }

    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the
    // L-Suffixes
    void induceSAl(uint_t *SA, int_t *s, int_t *cnt, int_t *bkt, int_t n,
                   int_t K, int cs, int level) {
        int_t i, j;
        // find heads of buckets
        get_buckets(cnt, bkt, K, false);
        //  if(level==0) bkt[0]++;
        for (i = 0; i < n; i++) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0)
                    if (chr(j) >= chr(SA[i])) {
                        SA[bkt[chr(j)]++] = j;
                    }
            }
        }
    }
    // Init buckets
    void init_buckets(int_t *bkt, int_t K) {
        int_t i;
        // clear all buckets
        for (i = 0; i < K; i++) {
            bkt[i] = 0;
        }
    }

    // Count frequencies
    void get_counts(int_t *s, int_t *bkt, int_t n, int_t K, int cs) {
        init_buckets(bkt, K);           // clear all buckets
        for (int_t i = 0; i < n; i++) { // compute the size of each bucket
            bkt[chr(i)]++;
        }
    }

    // Compute the head or end of each bucket
    void get_buckets(int_t *tmp, int_t *bkt, int_t K, bool end) {
        int_t sum = 0;
        // compute the size of each bucket
        if (end) {
            for (int_t i = 0; i < K; i++) {
                sum += tmp[i];
                bkt[i] = sum - 1;
            }
        } else {
            for (int_t i = 0; i < K; i++) {
                sum += tmp[i];
                bkt[i] = sum - tmp[i];
            }
        }
    }

    // Compute the head or end of each bucket
    void get_buckets(int_t *s, int_t *bkt, int_t n, int_t K, int cs, bool end) {
        int_t i, sum = 0;

        // clear all buckets
        init_buckets(bkt, K);
        // compute the size of each bucket
        for (i = 0; i < n; i++) {
            bkt[chr(i)]++;
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
    }

    bool sleq(unsigned char *s, int_t a, int_t b, size_t len, int cs,
              unsigned char sentinel) {

        size_t i;

        for (i = 0; i < len; i++) {

            if (chr(a) < chr(b))
                return 1;
            else if (chr(a) > chr(b))
                return 0;

            else if (chr(a) == sentinel &&
                     chr(b) == sentinel) { // $_i < $_j iff i < j
                if (a < b)
                    return 1;
                else
                    return 0;
            }
            a++;
            b++;
        }

        return 1;
    }

    void SAIS(int_t *s, uint_t *SA, int_t n, int_t K, int cs, int level) {
        int i, j;

        unsigned char *t =
            (unsigned char *)malloc(n / 8 + 1); // LS-type array in bits

        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        tset(n - 2, 0);
        tset(n - 1, 1); // the sentinel must be in s1, important!!!
        for (i = n - 3; i >= 0; i--)
            tset(i, (chr(i) < chr(i + 1) ||
                     (chr(i) == chr(i + 1) && tget(i + 1) == 1))
                        ? 1
                        : 0);

        int *bkt = (int *)malloc(sizeof(int) * K); // bucket counters

        // sort all the S-substrings
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets
        for (i = 0; i < n; i++)
            SA[i] = EMPTY;
        for (i = n - 2; i >= 0; i--)
            if (isLMS(i))
                SA[bkt[chr(i)]--] = i;
        SA[0] = n - 1; // set the single sentinel LMS-substring

        induceSAl(t, SA, s, bkt, n, K, cs, level);
        induceSAs(t, SA, s, bkt, n, K, cs, level);

        free(bkt);

        // compact all the sorted substrings into the first n1 items of s
        // 2*n1 must be not larger than n (proveable)
        int n1 = 0;
        for (i = 0; i < n; i++)
            if (isLMS(SA[i]))
                SA[n1++] = SA[i];

        // Init the name array buffer
        for (i = n1; i < n; i++)
            SA[i] = EMPTY;
        // find the lexicographic names of all substrings
        int name = 0, prev = -1;
        for (i = 0; i < n1; i++) {
            int pos = SA[i];
            bool diff = false;
            for (int d = 0; d < n; d++)
                if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||
                    chr(pos + d) != chr(prev + d) ||
                    tget(pos + d) != tget(prev + d)) {
                    diff = true;
                    break;
                } else if (d > 0 && (isLMS(pos + d) || isLMS(prev + d)))
                    break;

            if (diff) {
                name++;
                prev = pos;
            }
            pos = (pos % 2 == 0) ? pos / 2 : (pos - 1) / 2;
            SA[n1 + pos] = name - 1;
        }
        for (i = n - 1, j = n - 1; i >= n1; i--)
            if (SA[i] != EMPTY)
                SA[j--] = SA[i];

        // s1 is done now
        uint_t *SA1 = SA, *s1 = SA + n - n1;

        // stage 2: solve the reduced problem

        // recurse if names are not yet unique
        if (name < n1) {
            SAIS((int_t *)s1, SA1, n1, name, sizeof(int), level + 1);
        } else // generate the suffix array of s1 directly
            for (i = 0; i < n1; i++)
                SA1[s1[i]] = i;

        // stage 3: induce the result for the original problem

        bkt = (int *)malloc(sizeof(int) * K); // bucket counters

        // put all left-most S characters into their buckets
        get_buckets(s, bkt, n, K, cs, true); // find ends of buckets
        j = 0;
        for (i = 1; i < n; i++)
            if (isLMS(i))
                s1[j++] = i; // get p1
        for (i = 0; i < n1; i++)
            SA1[i] = s1[SA1[i]]; // get index in s1
        for (i = n1; i < n; i++)
            SA[i] = EMPTY; // init SA[n1..n-1]
        for (i = n1 - 1; i >= 0; i--) {
            j = SA[i];
            SA[i] = EMPTY;
            if (level == 0 && i == 0)
                SA[0] = n - 1;
            else
                SA[bkt[chr(j)]--] = j;
        }

        induceSAl(t, SA, s, bkt, n, K, cs, level);
        induceSAs(t, SA, s, bkt, n, K, cs, level);

        free(bkt);
        free(t);
    }
};

//
//
// template<>
// void gc_is_dictionary<gcis_eliasfano_codec>::extract(uint64_t l,
//                                                      uint64_t r,
//                                                      uint64_t level,
//                                                      sdsl::int_vector<>&
//                                                      extracted_text,
//                                                      sdsl::int_vector<>&
//                                                      tmp_text){
//
//     assert(l<=r);
//     // Base case: the reduced string is explictly stored
//     // just extract the interval.
//     if(level==g.size()){
//         for(uint64_t i=l;i<=r;i++){
//             extracted_text[i-l] = reduced_string[i];
//         }
//         return;
//     }
//
//     // Base case: the extracted string lies on the tail.
//     // just extractt it.
//     if(r<g[level].tail.size()){
//         g[level].extract_rules(l,r,extracted_text,tmp_text);
//         return;
//     }
//     uint64_t l2,r2;
//     // Compute the extraction interval in the level+1
//     l2 = g[level].lms_rnk1(std::max(g[level].tail.size(),l) + 1)-1;
//     r2 = g[level].lms_rnk1(r+1)-1;
//     // Recursively extract the text from [l2,r2] in level+1
//     extract(l2,r2,level+1,extracted_text,tmp_text);
//     std::swap(extracted_text,tmp_text);
//     // Use the extracted text from level+1 to extract the text from level
//     g[level].extract_rules(l,r,extracted_text,tmp_text);
//
//     if(r<g[level].tail.size()){
//         g[level].extract_rules(l,r,extracted_text,tmp_text);
//         return;
//     }
//     uint64_t l2,r2;
//     // Compute the extraction interval in the level+1
//     l2 = g[level].lms_rnk1(std::max(g[level].tail.size(),l) + 1)-1;
//     r2 = g[level].lms_rnk1(r+1)-1;
//     // Recursively extract the text from [l2,r2] in level+1
//     extract(l2,r2,level+1,extracted_text,tmp_text);
//     std::swap(extracted_text,tmp_text);
//     // Use the extracted text from level+1 to extract the text from level
//     g[level].extract_rules(l,r,extracted_text,tmp_text);
// }

// template<>
// sdsl::int_vector<> gc_is_dictionary<gcis_eliasfano_codec>::extract(uint64_t
// l,
//                                              uint64_t r){
//     sdsl::int_vector<> extracted_text(r-l+1);
//     sdsl::int_vector<> tmp_text(r-l+1);
//     extract(l,r,0,extracted_text,tmp_text);
//     return extracted_text;
// }

template <typename T> class gcis_dictionary : public gcis_abstract<T> {};

#endif
