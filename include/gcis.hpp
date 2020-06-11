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

#define tget(i) ((t[(i) >> 3] & mask[(i)&0x7]) ? 1 : 0)
#define tset(i, b)                                                             \
    t[(i) >> 3] =                                                              \
        (b) ? (mask[(i)&0x7] | t[(i) >> 3]) : ((~mask[(i)&0x7]) & t[(i) >> 3])

#define isLMS(i) (i > 0 && tget(i) && !tget(i - 1))

#define max(a, b) ((a) > (b) ? (a) : (b))

/**/
#define RMQ 2          // variants = (1, trivial) (2, using Gog's stack)
#define BINARY 0       // binary search on stack operations
#define STACK_SIZE 895 // to use 10Kb of working space

typedef struct _pair {
    uint_t idx;
    int_t lcp;
} t_pair;

int compare(const void *a, const void *b) {
    if (*(const uint_t *)a < *(const uint_t *)b)
        return -1;
    if (*(const uint_t *)a > *(const uint_t *)b)
        return 1;
    return 0;
}

void stack_push(t_pair *STACK, int_t *top, uint_t idx, int_t lcp) {
    STACK[*top].idx = idx;
    STACK[*top].lcp = lcp;
    (*top)++;
}
/**/

class gcis_interface {
  public:
    virtual void encode(char *s, int_t n) = 0;
    virtual pair<char *, int_t> decode() = 0;
    virtual void extract_batch(vector<pair<int, int>> &v_query) = 0;
    virtual pair<char *, int_t> decode_saca(uint_t **SA) = 0;
    virtual pair<char *, int_t> decode_saca_lcp(uint_t **SA, int_t **LCP) = 0;
    virtual uint64_t size_in_bytes() = 0;
    virtual void serialize(std::ostream &o) = 0;
    virtual void load(std::istream &i) = 0;
};

template <class codec_t> class gcis_abstract : public gcis_interface {
  public:
    std::vector<codec_t> g;
    sdsl::int_vector<> reduced_string;

  public:
    void extract_batch(vector<pair<int, int>> &v_query) {
        throw(gcis::util::NotImplementedException("extract_batch"));
    }
    pair<char *, int_t> decode_saca(uint_t **SA) {
        throw(gcis::util::NotImplementedException("decode_saca"));
    }
    pair<char *, int_t> decode_saca_lcp(uint_t **SA, int_t **LCP) {
        throw(gcis::util::NotImplementedException("decode_saca_lcp"));
    }
    virtual uint64_t size_in_bytes() {
        uint64_t total_bytes = 0;
        for (uint64_t i = 0; i < g.size(); i++) {
            total_bytes += g[i].size_in_bytes();
        }
        total_bytes += sdsl::size_in_bytes(reduced_string);
        return total_bytes;
    }

    void encode(char *s, int_t n) {
        uint_t *SA = new uint_t[n];
        int_t K = 256;
        int cs = sizeof(char);
        int level = 0;

        gc_is((int_t *)s, SA, n, K, cs, level);

        delete[] SA;
    }

    virtual pair<char *, int_t> decode() = 0;

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

    bool lcp_array_check(uint_t *SA, int_t *LCP, unsigned char *s, size_t len,
                         int cs, unsigned char sentinel) {
        int r = 1;
        uint_t i, j, k;
        int_t h;
        double mean = 0.0;

        for (i = 1; i < len; i++) {

            j = SA[i - 1];
            k = SA[i];
            for (h = 0; j + h < len && k + h < len; h++)
                if (s[j + h] != s[k + h])
                    break;

            if (LCP[i] != h) {
                cout << i << "\t" << h << "\t" << LCP[i] << " (wrong) " << endl;
                r = 0;
                // return 0;
            }
            mean += (double)LCP[i] / (double)len;
        }

        printf("LCP_mean = %lf\n", mean);

        return r;
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
        //        gcis::util::print_report("LCP Size (bits) =
        //        ",g[level].lcp.size(),"\n");
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

        bool premature_stop =
            evaluate_premature_stop(n, K, n1, name + 1, level);
        if (name + 1 < n1 && !premature_stop) {
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
            gcis::util::print_report(
                "Reduced String Length = ", (int_t)reduced_string.size(), "\n");
            gcis::util::print_report(
                "Reduced String Width (bits per symbol) = ",
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

        SA[bkt[chr(n - 1)]++] = n - 1;
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

    void compute_lcp_phi_sparse_sais(int_t *s, uint_t *SA1, uint_t *RA,
                                     int_t *LCP, int_t *PLCP, uint_t n1, int cs,
                                     int_t n) {

        uint_t i;

        /*
        PLCP[SA1[0]]=0;//PLCP* (lms) is stored in PLCP array
        for(i=1; i<n1; i++)
          PLCP[SA1[i]] = LCP[i];
        */

        LCP[SA1[0]] = 0; // PHI is stored in LCP array
        for (i = 1; i < n1; i++)
            LCP[SA1[i]] = SA1[i - 1]; // RA[SA1[i-1]];

        int_t l = 0; // q=0;
        for (i = 0; i < n1 - 1; i++) {

            if (l < 0)
                l = 0;
            // l=max(0,l);

            while (RA[i] + l < n && RA[LCP[i]] + l < n &&
                   chr(RA[i] + l) == chr(RA[LCP[i]] + l))
                l++;
            PLCP[i] = l;

            if (LCP[i] == n1 - 1)
                l -= RA[i + 1] - RA[i];
            else
                l -= max(RA[i + 1] - RA[i],
                         RA[LCP[i] + 1] -
                             RA[LCP[i]]); // LCP[i] stores the distance of i-th
                                          // suffix to its successor
        }

        LCP[0] = 0;
        for (i = 1; i < n1; i++)
            LCP[i] = PLCP[SA1[i]];
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

    // compute SA for the S-Type suffixes by inducing the L-Type suffixes and
    // the S-Type suffixes
    void induceSAs_LCP(uint_t *SA, int_t *LCP, int_t *s, int_t *cnt, int_t *bkt,
                       int_t n, int_t K, int cs, int level) {
        int_t i, j;
        get_buckets(cnt, bkt, K, true);

#if RMQ == 1
        int_t *M = (int_t *)malloc(sizeof(int_t) * K);
        for (i = 0; i < K; i++)
            M[i] = I_MAX;
#elif RMQ == 2
        uint_t *last_occ = (uint_t *)malloc(K * sizeof(uint_t));
        uint_t *tmp = (uint_t *)malloc(K * sizeof(uint_t));
        t_pair *STACK = (t_pair *)malloc((STACK_SIZE + 1) * sizeof(t_pair));
        int_t top = 0;
        // init
        stack_push(STACK, &top, n, -1);
        for (i = 0; i < K; i++)
            last_occ[i] = n - 1;
#endif

        for (i = n - 1; i >= 0; i--) {
            // if (SA[i] != EMPTY) {
            if (SA[i] > 0) {
                j = SA[i] - 1;
                if (j >= 0 && chr(j) <= chr(j + 1) && bkt[chr(j)] < i) {
                    SA[bkt[chr(j)]] = j;
#if RMQ == 1
                    if (LCP[bkt[chr(j)] + 1] >= 0)
                        LCP[bkt[chr(j)] + 1] = M[chr(j)] + 1;
#elif RMQ == 2
                    int_t min = I_MAX, end = top - 1;

                    int_t last = last_occ[chr(j)];
// search (can be binary)
#if BINARY == 1
                    int_t a = 0, b = top - 1;
                    int_t m = (b - a) / 2;
                    while (a <= b) {
                        if (STACK[m].idx == last) {
                            break;
                        }
                        if (STACK[m].idx < last)
                            b = m - 1;
                        else
                            a = m + 1;
                        m = a + (b - a) / 2;
                    }
                    end = m - 1;
#else
                    while (STACK[end].idx <= last)
                        end--;
#endif

                    min = STACK[(end + 1)].lcp;
                    last_occ[chr(j)] = i;

                    if (LCP[bkt[chr(j)] + 1] >= 0)
                        LCP[bkt[chr(j)] + 1] = min + 1;

#endif

#if RMQ == 1
                    if (LCP[bkt[chr(j)]] > 0)
                        LCP[bkt[chr(j)]] = I_MAX;
                    M[chr(j)] = I_MAX;
#endif

                    bkt[chr(j)]--;

                    if (SA[bkt[chr(j)]] != U_MAX) { // L/S-seam
                        int_t l = 0;
                        while (SA[bkt[chr(j)] + 1] + l < n &&
                               SA[bkt[chr(j)]] + l < n &&
                               chr(SA[bkt[chr(j)] + 1] + l) ==
                                   chr(SA[bkt[chr(j)]] + l))
                            l++;
                        LCP[bkt[chr(j)] + 1] = l;
                    }
                }
            }
            if (LCP[i] < 0)
                LCP[i] = 0;

#if RMQ == 1
            int_t k;
            for (k = 0; k < K; k++)
                if (M[k] > LCP[i])
                    M[k] = LCP[i];
#elif RMQ == 2

            int_t lcp = max(0, LCP[i]);

            while (STACK[(top)-1].lcp >= lcp)
                (top)--;
            stack_push(STACK, &top, i, lcp);

            if (top >= STACK_SIZE) {

                int_t j;
                memcpy(tmp, last_occ, K * sizeof(uint_t));
                qsort(tmp, K, sizeof(uint_t), compare);

                int_t curr = 0, end = 1;
                STACK[top].idx = U_MAX;

                for (j = K - 1; j >= 0; j--) {
                    if (STACK[end - 1].idx > tmp[j]) {
                        while (STACK[curr].idx > tmp[j])
                            curr++;
                        STACK[end].idx = STACK[curr].idx;
                        STACK[end].lcp = STACK[curr].lcp;
                        end++;
                    }
                }

                if (end >= STACK_SIZE) {
                    fprintf(stderr, "ERROR: induceSAl0_LCP\n");
                    exit(1);
                }
                top = end;
            }
#endif
        } // for
        LCP[0] = 0;

// variant 1
#if RMQ == 1
        free(M);
#elif RMQ == 2
        free(STACK);
        free(last_occ);
        free(tmp);
#endif
    }

    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the
    // L-Suffixes
    void induceSAl(uint_t *SA, int_t *s, int_t *cnt, int_t *bkt, int_t n,
                   int_t K, int cs, int level) {
        int_t i, j;
        // find heads of buckets
        get_buckets(cnt, bkt, K, false);
        SA[bkt[chr(n - 1)]++] = n - 1;
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
    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the
    // L-Suffixes
    void induceSAl_LCP(uint_t *SA, int_t *LCP, int_t *s, int_t *cnt, int_t *bkt,
                       int_t n, int_t K, int cs, int level) {
        int_t i, j;

        for (i = 0; i < K; i++)
            if (bkt[i] + 1 < n)
                if (SA[bkt[i] + 1] != U_MAX)
                    LCP[bkt[i] + 1] = I_MIN;

        // find heads of buckets
        get_buckets(cnt, bkt, K, false);

        for (i = 0; i < K; i++)
            if (bkt[i] < n)
                LCP[bkt[i]] = -2;

#if RMQ == 1
        int_t *M = (int_t *)malloc(sizeof(int_t) * K);
        for (i = 0; i < K; i++) {
            M[i] = I_MAX;
        }
#elif RMQ == 2
        uint_t *last_occ = (uint_t *)malloc(K * sizeof(uint_t));
        uint_t *tmp = (uint_t *)malloc(K * sizeof(uint_t));
        t_pair *STACK = (t_pair *)malloc((STACK_SIZE + 1) * sizeof(t_pair));
        int_t top = 0;
        // init
        stack_push(STACK, &top, 0, -1);
        for (i = 0; i < K; i++)
            last_occ[i] = 0;
#endif

        //  bkt[0]++;
        for (i = 0; i < n; i++) {
            if (SA[i] != U_MAX) {

                if (LCP[i] == I_MIN) { // is a L/S-seam position
                    int_t l = 0;
                    if (SA[bkt[chr(SA[i])] - 1] < n - 1)
                        while (
                            SA[i] + l < n && SA[bkt[char(SA[i])] - 1] + l < n &&
                            chr(SA[i] + l) == chr(SA[bkt[chr(SA[i])] - 1] + l))
                            ++l;
                    LCP[i] = l;
                }
#if RMQ == 1
                uint_t k;
                for (k = 0; k < K; k++)
                    if (M[k] > LCP[i])
                        M[k] = max(0, LCP[i]);
#elif RMQ == 2
                int_t min_lcp = 0;
                uint_t last;
                if (!SA[i])
                    last = 0;
                else {
                    last = last_occ[chr(SA[i] - 1)];
                    last_occ[chr(SA[i] - 1)] = i + 1;
                }
                int_t lcp = max(0, LCP[i]);
#if BINARY == 1
                int_t a = 0, b = top - 1;
                int_t m = (b - a) / 2;
                while (a <= b) {
                    if (STACK[m].lcp == lcp) {
                        break;
                    }
                    if (STACK[m].lcp > lcp)
                        b = m - 1;
                    else
                        a = m + 1;
                    m = a + (b - a) / 2;
                }
                top = m;
#else
                while (STACK[(top)-1].lcp >= lcp)
                    (top)--;
#endif
                stack_push(STACK, &top, i + 1, lcp);
                j = top - 1;
#if BINARY == 1
                a = 0, b = top - 1;
                m = (b - a) / 2;
                while (a <= b) {
                    if (STACK[m].idx == last) {
                        m++;
                        break;
                    }
                    if (STACK[m].idx > last)
                        b = m - 1;
                    else
                        a = m + 1;
                    m = a + (b - a) / 2;
                }
                j = m - 1;
#else
                while (STACK[j].idx > last)
                    j--;
#endif
                min_lcp = STACK[(j + 1)].lcp;
#endif

                // j = SA[i] - 1;
                if (SA[i] > 0) {
                    j = SA[i] - 1;
                    if (chr(j) >= chr(j + 1)) {
                        SA[bkt[chr(j)]] = j;
#if RMQ == 1
                        LCP[bkt[chr(j)]] += M[chr(j)] + 1;
                        M[chr(j)] = I_MAX;
#elif RMQ == 2
                        LCP[bkt[chr(j)]] += min_lcp + 1;
#endif
                        bkt[chr(j)]++;
                    }
                    if (bkt[chr(SA[i])] - 1 < i) { // if is LMS-type
                        SA[i] = U_MAX;
                    }
                }
#if RMQ == 2
                if (top >= STACK_SIZE) { // if stack is full
                    int_t j;
                    memcpy(tmp, last_occ, K * sizeof(uint_t));
                    qsort(tmp, K, sizeof(uint_t), compare);
                    int_t curr = 1, end = 1;
                    STACK[top].idx = U_MAX;

                    for (j = 0; j < K; j++) {
                        if (STACK[end - 1].idx < tmp[j] + 1) {
// search (can be binary)
#if BINARY == 1
                            int_t a = curr - 1, b = top - 1;
                            int_t m = (b - a) / 2, last = tmp[j] + 1;
                            while (a <= b) {
                                if (STACK[m].idx == last) {
                                    break;
                                }
                                if (STACK[m].idx > last)
                                    b = m - 1;
                                else
                                    a = m + 1;
                                m = a + (b - a) / 2;
                            }
                            curr = m;
#else
                            while (STACK[curr].idx < tmp[j] + 1)
                                curr++;
#endif

                            if (curr < top) {
                                STACK[end].idx = STACK[curr].idx;
                                STACK[end].lcp = STACK[curr].lcp;
                                end++;
                                curr++;
                            }
                        }
                    }
                    if (end >= STACK_SIZE) {
                        fprintf(stderr, "ERROR: induceSAl0_LCP\n");
                        exit(1);
                    }
                    top = end;
                }
#endif
            }
        }
#if RMQ == 1
        free(M);
#elif RMQ == 2
        free(STACK);
        free(last_occ);
        free(tmp);
#endif
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
