//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_ELIASFANO_NO_LCP_HPP
#define GC_IS_GCIS_ELIASFANO_NO_LCP_HPP

#include <iostream>
#include "gcis.hpp"
#include "gcis_eliasfano_codec_no_lcp.hpp"


template<>
class gcis_dictionary<gcis_eliasfano_codec_no_lcp> : public gcis_abstract<gcis_eliasfano_codec_no_lcp> {

public:

public:

    void serialize(std::ostream& o) override {
        gcis_abstract::serialize(o);
        uint32_t n = partial_sum.size();
        o.write((char*)&n,sizeof(n));
        o.write((char*)partial_sum.data(),sizeof(uint32_t) * n);
    }

    void load(std::istream& i) override {
        gcis_abstract::load(i);
        uint32_t n;
        i.read((char*)&n,sizeof(n));
        partial_sum.resize(n);
        i.read((char*)partial_sum.data(),sizeof(uint32_t) * n);
    }


    /**
     * Extracts any valid substring T[l,r] from the text
     * @param l Beggining of such substring
     * @param r End of such substring
     * @return Returns the extracted substring
     */
    sdsl::int_vector<> extract(uint64_t l,uint64_t r){
        uint64_t size = g.size() ? 2* (g.back().fully_decoded_tail_len + (r-l+1)) : (r-l+1);
        sdsl::int_vector<> extracted_text(size);
        sdsl::int_vector<> tmp_text(size);
        extract(l,r,extracted_text,tmp_text);
        return extracted_text;
    }


    pair<char*,int_t> decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char* str;
        if(g.size()) {
            for (int64_t i = g.size() - 1; i >= 0; i--) {
                sdsl::int_vector<> next_r_string;
                gcis_eliasfano_codec_no_lcp_level gd = std::move(g[i].decompress());
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
        }
        else{
            str = new char[reduced_string.size()];
            for(uint64_t i=0 ; i< reduced_string.size();i++){
                str[i] = reduced_string[i];
            }
        }
        return make_pair(str,g[0].string_size);
    }


private:
    std::vector<uint32_t> partial_sum;

private:



    void gc_is(int_t *s,
               uint_t *SA,
               int_t n,
               int_t K,
               int cs,
               int level) override {
        int_t i, j;

#       ifdef MEM_MONITOR
        mm.event("GC-IS Level " + to_string(level));
#       endif

#ifdef REPORT
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

        std::vector<uint64_t> fdrlen;

        uint_t rule_index = 0;
        sdsl::bit_vector rule_delim;
        g.push_back(gcis_eliasfano_codec_no_lcp());
        //Iterate over all suffixes in the LMS sorted array
        for (i = 0; i < n1; i++) {
            int_t pos = SA[i];
            bool diff = false;
            int_t d;
            // d equals to the LCP between two consecutive LMS-substrings
            for (d = 0; d < n; d++) {
                //If is first suffix in LMS order (sentinel), or one of the suffixes reached the last position of T, or the
                // characters of T differs or the type os suffixes differ.
                if (prev == -1 ||  (chr(pos + d) != chr(prev + d) && (d==0 ||(!isLMS(pos+d)  && !isLMS(prev+d)) ))
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

                // TODO: can we replace len =1 for len  =d?
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
                uint64_t old_rule_delim_size;

                old_rule_delim_size = rule_delim.size();

  
                g[level].rule.resize(g[level].rule.size() + len);
                rule_delim.resize(rule_delim.size() + len +1);

                for(uint64_t i=old_rule_delim_size;i<rule_delim.size();i++){
                    rule_delim[i] = 0;
                }

#ifdef REPORT
                total_rule_len+=len;
                total_rule_suffix_length += len;
#endif
                // Encode rule length in unary
                rule_delim[rule_delim.size() - 1] = 1;
                // Copy the symbols into the delimited rule positions
                for (j = 0; j < len   && j + pos  < n; j++) {
#ifdef REPORT
                    if(j+pos-1<n && chr(j+pos) == chr(j+pos+1)){
                        run_length_potential++;
                    }
#endif
                    g[level].rule[rule_index] = (uint_t) chr(j + pos);
                    rule_index++;
                }
                // Insert the fully decode rule length
                if(level==0){
                    // The symbols are terminal L(x) = 1, for every x
                    fdrlen.push_back(len);
                }
                else{
                    // The symbols are not necessarily terminal.            
                    uint64_t sum = g[level-1].fully_decoded_rule_len[chr(pos)];
                    for(uint64_t i=1;  i+pos < n && !isLMS(pos+i);i++){
                        sum+=g[level-1].fully_decoded_rule_len[chr(pos+i)];
                    }
                    fdrlen.push_back(sum);
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
        g[level].rule_suffix_length.encode(rule_delim);
        g[level].fully_decoded_rule_len = sdsl::dac_vector_dp<>(fdrlen);

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
            g[level].tail[j] = (uint64_t) (cs == sizeof(char) ? ((char*)s)[j] :
                                           s[j]);
        }
        // Compress the tail
        sdsl::util::bit_compress(g[level].tail);
        // Determine the accumulated tail length
        if(level>0) {
            g[level].fully_decoded_tail_len = g[level-1].fully_decoded_tail_len;
            for (auto t: g[level].tail) {
                g[level].fully_decoded_tail_len += g[level-1].fully_decoded_rule_len[t];
            }
        }
        else {
            g[level].fully_decoded_tail_len = g[level].tail.size();
        }


        // stage 2: solve the reduced problem
        // recurse if names are not yet unique

#ifdef REPORT
        gcis::util::print_report("Level ",level,"\n");
        gcis::util::print_report("Alphabet Size = ",K, "\n");
        gcis::util::print_report("String Size = ",n,"\n");
        gcis::util::print_report("Number of Rules = ", name+1, "\n");
        gcis::util::print_report("Average Rule Length = ",(double) total_rule_len/(name+1),"\n");
        gcis::util::print_report("Number of Discarded Rules = ",discarded_rules_n,"\n");
        gcis::util::print_report("Average Discarded Rules Length = ",(double) discarded_rules_len/discarded_rules_n,"\n");
        gcis::util::print_report("Average Rule Suffix Length = ",(double) total_rule_suffix_length/(name+1),"\n");
        gcis::util::print_report("Dictionary Level Size (bytes) =",g[level].size_in_bytes(),"\n");
        gcis::util::print_report("Rule Suffix Length (total) = ",g[level].rule.size(),"\n");
        gcis::util::print_report("Rule Suffix Width (bits per symbol) = ",(int_t) g[level].rule.width(),"\n");
        gcis::util::print_report("Tail Length = ",g[level].tail.size(),"\n");
        gcis::util::print_report("Tail Width (bits per symbol) = ",(int_t) g[level].tail.width(),"\n");
        gcis::util::print_report("Run Length Potential (total) = ",run_length_potential,"\n");
        gcis::util::print_report("Avg Run Length per Rule Suffix = ",(double) run_length_potential/(name+1),"\n");
#endif

        // TODO: reenable premature_stop comparison
        // bool premature_stop = false;
       bool premature_stop = evaluate_premature_stop(n,K,n1,name+1,level);
        g[level].string_size = n;
        g[level].alphabet_size = K;
        if (name+1 < n1 && !premature_stop) {
            gc_is((int_t *) s1, SA1, n1, name+1, sizeof(int_t), level + 1);
        }
        else {
            // generate the suffix array of s1 directly
            if(premature_stop){
                // The encoding algorithm is stopped prematurely
#ifdef REPORT
                gcis::util::print_report("Premature Stop employed at level ", level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = (uint64_t) (cs == sizeof(char) ? ((char*)s)[j] :
                                                    s[j]);
                }
                // The last level of computation is discarded
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
            partial_sum.resize(reduced_string.size());
            partial_sum[0] = 0;
            for(uint64_t i = 1;i<reduced_string.size();i++){
                partial_sum[i] = partial_sum[i-1] + g.back().fully_decoded_rule_len[reduced_string[i-1]];
                // cout << "Partial sum = " << partial_sum[i] << "\n";
            }

#ifdef REPORT
            gcis::util::print_report("Reduced String Length = ",(int_t) reduced_string.size(),"\n");
            gcis::util::print_report("Reduced String Width (bits per symbol) = ",(int_t) reduced_string.width(),"\n");
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
    uint64_t bsearch_upperbound(vector<uint32_t>& partial_sum, uint64_t sz) {
        if(partial_sum.back()<=sz){
            return partial_sum.size()-1;
        }
        int64_t l=0,r=partial_sum.size()-1;
        while(l<r){
            int64_t mid = l+ (r-l)/2;
            if(partial_sum[mid] + g.back().fully_decoded_rule_len[reduced_string[mid]] > sz){
                r = mid;
            }
            else{
                l = mid+1;
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
    uint64_t bsearch_lowerbound(vector<uint32_t>& partial_sum, uint64_t sz) {
        if(partial_sum[0] >= sz){
            return 0;
        }
        int64_t l=0,r=partial_sum.size()-1;
        while(l<r){
            int64_t mid = l+ (r-l+1)/2;
            if(partial_sum[mid] <= sz){
                l = mid;
            }
            else{
                r = mid-1;
            }
        }
        return l;
    }

    /**
     *
     * @param codec The dictionary level we are employing
     * @param extracted_text The string to be extracted
     * @param sz The index we want to find
     * @param text_r The tracked position of the original text representing the end of extracted_text
     * @return The leftmost index such that text_l r sum(extracted_text,i) <= sz
     */
    uint64_t sequential_upperbound(gcis_eliasfano_codec_no_lcp &codec,
                                   sdsl::int_vector<>& extracted_text,
                                   int64_t extracted_text_size,
                                   int64_t sz, int64_t& text_r) {
        int64_t index = 0;
        int64_t rule_length = 0;
        for(index = extracted_text_size - 1 ; index>=0; index--){
            rule_length = codec.fully_decoded_rule_len[extracted_text[index]];
            if(text_r -rule_length <= sz) {
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
     * @param text_l The tracked position of the original text representing the beggining of extracted_text
     * @return The rightmost index such that text_l + sum(extracted_text,i) <= sz
    */
    uint64_t sequential_lowerbound(gcis_eliasfano_codec_no_lcp &codec,
                                   sdsl::int_vector<>& extracted_text,
                                   int64_t extracted_text_size,
                                   int64_t sz, int64_t& text_l) {
        int64_t index = 0;
        int64_t rule_length = 0;
        for(index = 0 ; index<extracted_text_size; index++){
            rule_length = codec.fully_decoded_rule_len[extracted_text[index]];
            if(text_l +rule_length > sz) {
                break;
            }
            text_l += rule_length;
        }
        return index;
    }

    void extract(int64_t l, int64_t r, sdsl::int_vector<>& extracted_text, sdsl::int_vector<>& tmp_text){
//        // Stores the interval being tracked in the text
        int64_t text_l;
        int64_t text_r;
        // Stores the interval being tracked in the level
        uint64_t lk,rk;
        text_l= 0;
        text_r = g.size()>0 ? g[0].string_size : reduced_string.size();
        uint64_t extracted_idx = 0;

        /**
         * The dictionary is only the original string
         * The extraction is done in a straight-foward fashion
         */
        if(g.size()==0){
            for(uint64_t j=0;j<reduced_string.size();j++){
                extracted_text[j] = reduced_string[j];
            }
            return;
        }

        /**
         * Else, the extraction should proceed from the reduced string to the
         * decompressed text
         */

        if(r < g.back().fully_decoded_tail_len) {
            // The string lies on the tail. Copy all the tail.
            text_l = 0;
            text_r = g.back().fully_decoded_tail_len -1;
            for (auto v:g.back().tail) {
                tmp_text[extracted_idx++] = v;
            }
        }
        else if(l<g.back().fully_decoded_tail_len){
            // A prefix of the string lies on the tail.
            // Copy the tail
            text_l = 0;
            for (auto v:g.back().tail) {
                tmp_text[extracted_idx++] = v;
            }
            // Find the leftmost index which covers r
            rk = bsearch_upperbound(partial_sum,r-g.back().fully_decoded_tail_len);
            text_r = g.back().fully_decoded_tail_len + partial_sum[rk] + g.back().fully_decoded_rule_len[reduced_string[rk]];
            // Decompress the rules located at reduced_string[0..rk];
            for(uint64_t i=0;i<=rk;i++){
                g.back().extract_rule(reduced_string[i],tmp_text,extracted_idx);
            }
        }
        else{
            // The string does not occur in the tail
            // Find the rightmost index which covers l
            lk = bsearch_lowerbound(partial_sum,l-g.back().fully_decoded_tail_len);
            text_l = g.back().fully_decoded_tail_len + partial_sum[lk];
            // Find the leftmost index which covers r
            rk = bsearch_upperbound(partial_sum,r-g.back().fully_decoded_tail_len);
            text_r = g.back().fully_decoded_tail_len + partial_sum[rk] + g.back().fully_decoded_rule_len[reduced_string[rk]];
            // Decompress the rules located at reduced_string[0..rk];
            for(uint64_t i=lk;i<=rk;i++){
                g.back().extract_rule(reduced_string[i],tmp_text,extracted_idx);
            }
        }
        int64_t level = g.size()-2;
        // Extract the reduced string part
        while(level>=0){
            uint64_t extracted_text_len = extracted_idx;
            extracted_idx = 0;
            std::swap(extracted_text,tmp_text);
            // The extracted string lies on the tail
            if(r<g[level].fully_decoded_tail_len) {
                text_l = 0;
                text_r = g[level].fully_decoded_tail_len-1;
                // Copy the tail
                for (auto v:g[level].tail) {
                    tmp_text[extracted_idx++] = v;
                }
            }
            else if(l<g[level].fully_decoded_tail_len){
                // A prefix of the string lies on the tail
                text_l = 0;
                // Copy the tail
                for (auto v:g[level].tail) {
                    tmp_text[extracted_idx++] = v;
                }
                rk = sequential_upperbound(g[level],extracted_text,extracted_text_len,r,text_r);
                for(uint64_t i=0;i<=rk;i++){
                    g[level].extract_rule(extracted_text[i],tmp_text,extracted_idx);
                }
            }
            else{
                text_l = std::max<int64_t>(text_l,g[level].fully_decoded_tail_len);
                lk = sequential_lowerbound(g[level],extracted_text,extracted_text_len,l,text_l);
                rk = sequential_upperbound(g[level],extracted_text,extracted_text_len,r,text_r);
                for(uint64_t i=lk;i<=rk;i++){
                    g[level].extract_rule(extracted_text[i],tmp_text,extracted_idx);
                }
            }
            level --;
        }

        for(uint64_t i = 0,offset=l-text_l;i<r-l+1;i++){
            extracted_text[i] = tmp_text[i+offset];
        }
    }

};


#endif //GC_IS_GCIS_ELIASFANO_NO_LCP_HPP
