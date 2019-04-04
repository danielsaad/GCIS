//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_S8B_HPP
#define GC_IS_GCIS_S8B_HPP

#include "gcis.hpp"
#include "gcis_s8b_codec.hpp"
#include "util.hpp"
template <>
class gcis_dictionary<gcis_s8b_codec> : public gcis_abstract<gcis_s8b_codec> {
public:

    char* decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char* str= 0;
        for(int64_t i=g.size()-1 ;i>=0 ; i--){
            sdsl::int_vector<> next_r_string;
            gcis_s8b_codec_level gd = std::move(g[i].decompress());
            gd.rule_delim_sel.set_vector(&gd.rule_delim);
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size-1) +1 );
            next_r_string.resize(g[i].string_size);
            uint64_t l = 0;
            if(i==0) {
                // Convert the reduced string in the original text
                str = new char[g[i].string_size];
                for(auto t: g[i].tail){
                    str[l++] = (char) t;
                }
                for(uint64_t j=0;j < r_string.size();j++){
                    gd.expand_rule(r_string[j],str,l);
                }
            }
            else{
                // Convert the reduced string in the previous reduced string
                for(uint64_t j=0;j<g[i].tail.size();j++){
                    next_r_string[l++] = g[i].tail[j];
                }
                for(uint64_t j=0;j<r_string.size();j++){
                    gd.expand_rule(r_string[j],next_r_string,l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return str;
    }

private:
    void gc_is(int_t *s,
               uint_t *SA, int_t n,
               int_t K,
               int cs,
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

        unsigned char *t = new unsigned char[n / 8 + 1]; // LS-type array in bits
        // stage 1: reduce the problem by at least 1/2

        // Classify the type of each character
        //  tset(n - 2, 0);
        tset(n - 1, 1); // the sentinel must be in s1, important!!!
        for (i = n - 2; i >= 0; i--)
            tset(i, (chr(i) < chr(i + 1) || (chr(i) == chr(i + 1) && tget(i + 1) == 1)) ? 1 : 0);

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

        //Induce L-Type suffixes by using LMS-Type and L-Type suffixes
        induceSAl(t, SA, s, bkt, n, K, cs, level);

        //Induce S-Type suffixes by using L-Type and S-Type suffixes
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

        // find the lexicographic names of all LMS-substrings by comparing the consecutive ones
        int_t name = -1;
        int_t prev = -1;

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        g.push_back(gcis_s8b_codec());
        //Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            bool diff = false;
            int_t d;
            // d equals to the LCP between two consecutive LMS-substrings
            for (d = 0; d < n; d++) {
                //If is first suffix in LMS order (sentinel), or one of the suffixes reached the last position of T, or the
                // characters of T differs or the type os suffixes differ.
                if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||  chr(pos + d) != chr(prev + d)
                    || ( isLMS(pos+d) ^ (isLMS(prev+d))) )  {
                    diff = true;
                    break;
                }
                //The comparison has reached the end of at least one LMS-substring
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

                g[level].lcp.encode(d);
                g[level].rule_suffix_length.encode(len-d);
                g[level].rule.resize(g[level].rule.size() + len - d);

#ifdef REPORT
                total_rule_len+=len;
                total_lcp +=d;
                total_rule_suffix_length += len-d;
#endif

                for (j = 0; j < len - d  && j + pos + d < n; j++) {
#ifdef REPORT
                    if(j+pos+d-1<n && chr(j+pos+d) == chr(j+pos+d+1)){
                        run_length_potential++;
                    }
#endif
                    g[level].rule[rule_index] = (uint_t) chr(j + pos + d);
                    rule_index++;
                }
                name++;
                prev = pos;
            }
#ifdef REPORT
            else{
                size_t len = 1;
                if (pos != n - 1)
                    while (!isLMS(pos + len))
                        len++;
                discarded_rules_len+=len;
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
            g[level].tail[j] = (uint64_t) (cs == sizeof(char) ? ((char*)s)[j] :
                                           s[j]);
        }
        sdsl::util::bit_compress(g[level].tail);

        // stage 2: solve the reduced problem
        // recurse if names are not yet unique

#ifdef REPORT
        print_report("Level ",level,"\n");
        print_report("Alphabet Size = ",K, "\n");
        print_report("String Size = ",n,"\n");
        print_report("Number of Rules = ", name+1, "\n");
        print_report("Average Rule Length = ",(double) total_rule_len/(name+1),"\n");
        print_report("Number of Discarded Rules = ",discarded_rules_n,"\n");
        print_report("Average Discarded Rules Length = ",(double) discarded_rules_len/discarded_rules_n,"\n");
        print_report("Average LCP = ",(double) total_lcp/(name+1),"\n");
        print_report("Average Rule Suffix Length = ",(double) total_rule_suffix_length/(name+1),"\n");
        print_report("Dictionary Level Size (bytes) =",g[level].size_in_bytes(),"\n");
        print_report("LCP Size (bits) = ",g[level].lcp.size(),"\n");
        print_report("Rule Suffix Length (total) = ",g[level].rule.size(),"\n");
        print_report("Rule Suffix Width (bits per symbol) = ",(int_t) g[level].rule.width(),"\n");
        print_report("Tail Length = ",g[level].tail.size(),"\n");
        print_report("Tail Width (bits per symbol) = ",(int_t) g[level].tail.width(),"\n");
        print_report("Run Length Potential (total) = ",run_length_potential,"\n");
        print_report("Avg Run Length per Rule Suffix = ",(double) run_length_potential/(name+1),"\n");
#endif

        bool premature_stop = evaluate_premature_stop(n,K,n1,name+1,level);
        if (name+1 < n1 && !premature_stop) {
            g[level].string_size = n;
            g[level].alphabet_size = K;
            gc_is((int_t *) s1, SA1, n1, name+1, sizeof(int_t), level + 1);
        } else { // generate the suffix array of s1 directly
            if(premature_stop){
#ifdef REPORT
                print_report("Premature Stop employed at level ",level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s[j];
                }
                g.pop_back();
            }
            else{
                reduced_string.resize(n1);
                for (j = 0; j < n1; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s1[j];
                }
            }
            sdsl::util::bit_compress(reduced_string);

#ifdef REPORT
            print_report("Reduced String Length = ",(int_t) reduced_string.size(),"\n");
            print_report("Reduced String Width (bits per symbol) = ",(int_t) reduced_string.width(),"\n");
#endif
        }
        delete[] t;
    }

};


#endif //GC_IS_GCIS_S8B_HPP
