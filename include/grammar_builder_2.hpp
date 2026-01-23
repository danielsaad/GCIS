#include "util.hpp"

#include "sais.h"
#include "sais_yuta.hpp"

#ifndef UCHAR_SIZE
#define UCHAR_SIZE 256
#endif
#ifndef MINBUCKETSIZE
#define MINBUCKETSIZE 256
#endif

#define sais_index_type int_t
#define sais_bool_type int_t
#define SAIS_LMSSORT2_LIMIT 0x3fffffff

#define SAIS_MYMALLOC(_num, _type) ((_type *)malloc((_num) * sizeof(_type)))
#define SAIS_MYFREE(_ptr, _num, _type) free((_ptr))
#define chr(_a)                                                                \
    (cs == sizeof(sais_index_type) ? ((sais_index_type *)T)[(_a)]              \
                                   : ((unsigned char *)T)[(_a)])

namespace gcis {

template <class grammar_t> class grammar_builder_2 {
  public:
    grammar_builder_2() = default;

    grammar_t build(char *s) {
        int n = strlen(s);
        sais_index_type *SA = new uint_t[n];

        m_ghelper.pre_process();

        gcis(m_text, SA, 0, n, UCHAR_SIZE, sizeof(unsigned char), 0);

        m_ghelper.post_process();

        delete[] SA;
        return m_ghelper.get_grammar();
    }

  private:
    typename grammar_t::grammar_builder_t m_ghelper;
    const void *m_text;
    sais_index_type m_text_size;
    sais_index_type* *m_SA;
    sais_index_type m_alphabet_size;
    int_t m_cs;
    sais_index_type m_level = 0;
    unsigned char *m_t_ptr;

    /**
     * @brief Induces the suffixes of type L and S, thus, sorting the
     * LMS substrings.
     * After this step, all LMS-substrings are in sorted order.
     * The suffix array at positions S[0,n1-1] will contain the starting
     * position of such LMS substrings after this.
     * @return The number of LMS substrings
     */
    int lms_sort() {
        int_t i;
        // stage 1: reduce the problem by at least half

        // Classify the type of each character
        tset(m_t_ptr, m_text_size - 1,
             1); // the sentinel must be in s1, important!!!
        for (i = m_text_size - 2; i >= 0; i--)
            tset(m_t_ptr, i,
                 (chr(m_text, m_cs, i) < chr(m_text, m_cs, i + 1) ||
                  (chr(m_text, m_cs, i) == chr(m_text, m_cs, i + 1) &&
                   tget(m_t_ptr, i + 1) == 1))
                     ? 1
                     : 0);

        int_t *bkt = new int_t[m_alphabet_size]; // bucket counters

        size_t first = m_text_size - 1;

        // sort all the S-substrings
        // find ends of buckets
        get_buckets(m_text, bkt, m_text_size, m_alphabet_size, m_cs, true);

        for (i = 0; i < (int_t)m_text_size; i++) {
            m_SA[i] = EMPTY;
        }

        for (i = m_text_size - 2; i >= 0; i--) {
            if (isLMS(m_t_ptr, i)) {
                m_SA[bkt[chr(m_text, m_cs, i)]--] = i;
                first = i;
            }
        }

        m_SA[0] = m_text_size - 1; // set the single sentinel LMS-substring

        // Induce L-Type suffixes by using LMS-Type and L-Type suffixes
        induceSAl(m_t_ptr, m_SA, m_text, bkt, m_text_size, m_alphabet_size,
                  m_cs, m_level);

        // Induce S-Type suffixes by using L-Type and S-Type suffixes
        induceSAs(m_t_ptr, m_SA, m_text, bkt, m_text_size, m_alphabet_size,
                  m_cs, m_level);

        delete[] bkt;

        // compact all the sorted substrings into the first n1 items of s
        // 2*n1 must be not larger than n (proveable)
        // n1 contains the end of the lms positions
        int_t n1 = 0;
        for (i = 0; i < (int_t)m_text_size; i++) {
            if (isLMS(m_t_ptr, m_SA[i])) {
                m_SA[n1++] = m_SA[i];
            }
        }

        // Init the name array buffer
        // SA[0,n1-1] = LMS starting positions
        // SA[n1,n-1] = Name for each LMS substring
        for (i = n1; i < (int_t)m_text_size; i++) {
            m_SA[i] = EMPTY;
        }
        return n1;
    }

    /**
     * @brief Compare LMS substrings starting at positions i and j
     * respectively
     *
     * @param i Starting position of the first LMS substrings
     * @param j Starting position of the second LMS
     * @return  a pair (b,l). b indicates wether the LMS substrings are
     * different and l indicates the LCP shared between the LMS substrings.
     */
    std::pair<bool, int_t> compare_lms_substrings(uint_t i, uint_t j) {
        bool diff = false;
        int_t lcp = 0;
        // d equals to the LCP between two consecutive LMS-substrings
        for (lcp = 0; lcp < (int_t)m_text_size; lcp++) {
            // If is first suffix in LMS order (sentinel), or one of the
            // suffixes reached the last position of T, or the
            // characters of T differs or the type os suffixes differ.
            if (j == -1 ||
                (chr(m_text, m_cs, i + lcp) != chr(m_text, m_cs, j + lcp) &&
                 (lcp == 0 ||
                  (!isLMS(m_t_ptr, i + lcp) && !isLMS(m_t_ptr, j + lcp)))) ||
                (isLMS(m_t_ptr, i + lcp) ^ (isLMS(m_t_ptr, j + lcp)))) {
                diff = true;
                break;
            }
            // The comparison has reached the end of at least one
            // LMS-substring
            // TODO: is this really needed?
            if (lcp > 0 &&
                (isLMS(m_t_ptr, i + lcp) || isLMS(m_t_ptr, j + lcp))) {
                break;
            }
        }
        return std::make_pair(diff, lcp);
    }
    /**
     * @brief Rename the rules.
     * After this step,  the rules name will be contained at the end
     * of SA. The unique LMS substrings starting position will be recorded
     * at the start of SA.
     * @param n1, the number of lms substrings.
     * @return the LCP between unique LMS substrings
     */
    std::vector<int_t> rename_lms(int_t n1) {
        // find the lexicographic names of all LMS-substrings by comparing the
        // consecutive ones
        int_t name = -1;
        int_t prev = -1;
        std::vector<int_t> lms_substrings_lcp;
        uint_t rule_index = 0;
        // Iterate over all suffixes in the LMS sorted array
        for (int_t i = 0; i < n1; i++) {
            int_t pos = m_SA[i];
            auto diff = compare_lms_substrings(pos, prev);

            // The consecutive LMS-substrings differs
            if (diff.first) {
                lms_substrings_lcp.push_back(diff.second);
                m_SA[rule_index++] = pos;
                name++;
                prev = pos;
            }
            pos = (pos % 2 == 0) ? pos / 2 : (pos - 1) / 2;
            m_SA[n1 + pos] = name;
        }
        // Compact every LMS substring name to the end of SA
        // After this SA[n-n1,n-1] will have the names of the LMS substrings
        // or in other words, the reduced text.
        for (int_t i = m_text_size - 1, j = m_text_size - 1; i >= n1; i--) {
            if (m_SA[i] != EMPTY) {
                m_SA[j--] = m_SA[i];
            }
        }
        return lms_substrings_lcp;
    }

    // gcis(m_text,SA,0,n,UCHAR_SIZE,sizeof(unsigned char),0);

    void gcis(const void *s, sais_index_type *SA, sais_index_type fs,
              sais_index_type n, sais_index_type k, int cs) {
        assert((T != NULL) && (SA != NULL));
        assert((0 <= fs) && (0 < n) && (1 <= k));

        m_text = s;
        m_SA = SA;
        m_text_size = n;
        m_alphabet_size = K;
        m_cs = cs;
        m_t_ptr = t;

        if (k <= MINBUCKETSIZE) {
            if ((C = SAIS_MYMALLOC(k, sais_index_type)) == NULL) {
                return -2;
            }
            if (k <= fs) {
                B = SA + (n + fs - k);
                flags = 1;
            } else {
                if ((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) {
                    SAIS_MYFREE(C, k, sais_index_type);
                    return -2;
                }
                flags = 3;
            }
        } else if (k <= fs) {
            C = SA + (n + fs - k);
            if (k <= (fs - k)) {
                B = C - k;
                flags = 0;
            } else if (k <= (MINBUCKETSIZE * 4)) {
                if ((B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) {
                    return -2;
                }
                flags = 2;
            } else {
                B = C;
                flags = 8;
            }
        } else {
            if ((C = B = SAIS_MYMALLOC(k, sais_index_type)) == NULL) {
                return -2;
            }
            flags = 4 | 8;
        }
        if ((n <= SAIS_LMSSORT2_LIMIT) && (2 <= (n / k))) {
            if (flags & 1) {
                flags |= ((k * 2) <= (fs - k)) ? 32 : 16;
            } else if ((flags == 0) && ((k * 2) <= (fs - k * 2))) {
                flags |= 32;
            }
        }

        m_ghelper.set_gcis_structures(s, SA, n, K, cs, level, m_t_ptr);
        int_t n1 = lms_sort();
        // stage 2: solve the reduced problem
        // recurse if names are not yet unique
        auto lcp = rename_lms(n1);
        m_ghelper.update_grammar(lcp);

        // s1 is done now
        uint_t *SA1 = SA, *s1 = SA + n - n1;
        uint_t name_n = lcp.size();

        if (name_n < n1) {
            // There is repeated non-terminals, recurse.
            gc_is((int_t *)s1, SA1, n1, name_n, sizeof(int), level + 1);

        } else {
            // all non-terminals are unique, compute the reduced string.
            m_text_size = n1;
            m_alphabet_size = n1;
            m_SA = SA;
            m_cs = cs;
            m_t_ptr = t;
            m_ghelper.set_gcis_structures((int_t *)s1, SA1, n1, name_n,
                                          sizeof(int), level, m_t_ptr);
            m_ghelper.compute_reduced_string();
        }
        delete[] t;
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
                if (j >= 0 && !tget(t, j)) {
                    SA[bkt[chr(s, cs, j)]++] = j;
                }
            }
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
                if (j >= 0 && tget(t, j)) {
                    SA[bkt[chr(s, cs, j)]--] = j;
                }
            }
        }
    }

    // Compute the head or end of each bucket
    void get_buckets(int_t *tmp, int_t *bkt, uint_t K, bool end) {
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
            bkt[chr(s, cs, i)]++;
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
            bkt[chr(s, cs, i)]++;
        }
    }
};

} // namespace gcis