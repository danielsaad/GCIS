//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_UNARY_HPP
#define GC_IS_GCIS_UNARY_HPP

#include "gcis.hpp"
#include "gcis_unary_codec.hpp"


template <>
class gcis_dictionary<gcis_unary_codec> : gcis_abstract<gcis_unary_codec>{
public:

    uint64_t size_in_bytes(){
        uint64_t total_bytes=0;
        for(uint64_t i=0;i<g.size();i++){
            total_bytes+=g[i].size_in_bytes();
        }
        total_bytes+= sdsl::size_in_bytes(reduced_string);
        return total_bytes;
    }

    void encode(char* s){
        int_t n = strlen(s)+1;
        uint_t* SA = new uint_t[n];
        int_t K = 256;
        int cs = sizeof(char);
        int level = 0;

        gc_is((int_t*) s, SA,n,K,cs,level);
        delete[] SA;
    }

    char* decode_inline(){
        int level=0;
        sdsl::int_vector<> final_string = gd_is(level);
        char* s = new char[final_string.size()+1];
        for(uint64_t i=0;i<final_string.size();i++){
            s[i]=final_string[i];
        }
        return s;
    }
    void extract(uint64_t l,
                 uint64_t r,
                 uint64_t level,
                 sdsl::int_vector<>& extracted_text,
                 sdsl::int_vector<>& tmp_text){}

    char* decode(){
        sdsl::int_vector<> r_string = reduced_string;
        char* str;
        for(int64_t i=g.size()-1 ;i>=0 ; i--){
            sdsl::int_vector<> next_r_string;
            gcis_unary_codec_level gd = std::move(g[i].decompress());
            gd.rule_delim_sel.set_vector(&gd.rule_delim);
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size-1) +1 );
            next_r_string.resize(g[i].string_size);
            uint64_t l = 0;
            if(i==0) {
                // Convert the reduced string in the original text
                str = new char[g[i].string_size];
                for(uint64_t j=0;j<g[i].tail.size();j++){
                    str[l++] = g[i].tail[j];
                }
                for(uint64_t j=0;j<r_string.size();j++){
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


    void serialize(std::ostream &o) {
        reduced_string.serialize(o);
        uint64_t size = g.size();
        o.write((char*) &size,sizeof(uint64_t));
        for(uint64_t i=0;i<g.size();i++){
            g[i].serialize(o);
        }
    }

    void load(std::istream &i) {
        uint64_t size;
        reduced_string.load(i);
        i.read((char*) &size,sizeof(uint64_t));
        g.resize(size);
        for(uint64_t j=0;j<size;j++){
            g[j].load(i);
        }
    }

private:

    bool evaluate_premature_stop(int_t n0,int_t alphabet_size_s0,
                                 int_t n1,int_t alphabet_size_s1,int_t level){

        uint64_t number_of_bits_s0 = sdsl::bits::hi(alphabet_size_s0)+1;
        uint64_t number_of_bits_s1 = sdsl::bits::hi(alphabet_size_s1)+1;
        uint64_t  dict_size_bits = g[level].size_in_bytes()*8;
        if(n0*number_of_bits_s0 < n1 *number_of_bits_s1 + dict_size_bits )
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

    void gc_is(int_t *s,
               uint_t *SA,
               int_t n,
               int_t K,
               int cs,
               int level) override{
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
        getBuckets(s, bkt, n, K, cs, true); // find ends of buckets

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

        // find the lexicographic names of all LMS-substrings by comparing the consecutive ones
        int_t name = -1;
        int_t prev = -1;

        int_t last_set_lcp_bit = -1;
        uint_t rule_index = 0;
        g.push_back(gcis_unary_codec());
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

                // Resizes Rule array, LCP array and rule delimiter bitvector
                uint64_t old_lcp_size, old_rule_delim_size;

                old_lcp_size = g[level].lcp.size();
                old_rule_delim_size = g[level].rule_delim.size();

                g[level].rule.resize(g[level].rule.size() + len - d);
                g[level].lcp.resize(g[level].lcp.size()+d+1);
                g[level].rule_delim.resize(g[level].rule_delim.size() + len - d+1);

                for(uint64_t i=old_lcp_size;i<g[level].lcp.size();i++){
                    g[level].lcp[i] = 0;
                }
                for(uint64_t i=old_rule_delim_size;i<g[level].rule_delim.size();i++){
                    g[level].rule_delim[i] = 0;
                }

#ifdef REPORT
                total_rule_len+=len;
                total_lcp +=d;
                total_rule_suffix_length += len-d;
#endif

                // Encode LCP value in unary
                g[level].lcp[g[level].lcp.size()-1] = 1;
                // Encode rule length in unary
                g[level].rule_delim[g[level].rule_delim.size() - 1] = 1;
                // Copy the symbols into the delimited rule positions
                for (j = 0; j < len - d  && j + pos + d < n; j++) {
#ifdef REPORT
                    if(j+pos+d-1<n && chr(j+pos+d) == chr(j+pos+d+1)){
                        run_length_potential++;
                    }
#endif
                    g[level].rule[rule_index] = (uint_t) chr(j + pos + d);
                    rule_index++;
                }
                // Since the adjacent LMS substrings differ, we must assign
                // a new name
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
            // Insert the name in upper half of the SA array
            // just after the LMS substring positions
            SA[n1 + pos] = name;
        }

        // Compress the  part of the rules not encoded by the lcp information
        sdsl::util::bit_compress(g[level].rule);

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
        gcis::util::print_report("Level ",level,"\n");
        gcis::util::print_report("Alphabet Size = ",K, "\n");
        gcis::util::print_report("String Size = ",n,"\n");
        gcis::util::print_report("Number of Rules = ", name+1, "\n");
        gcis::util::print_report("Average Rule Length = ",(double) total_rule_len/(name+1),"\n");
        gcis::util::print_report("Number of Discarded Rules = ",discarded_rules_n,"\n");
        gcis::util::print_report("Average Discarded Rules Length = ",(double) discarded_rules_len/discarded_rules_n,"\n");
        gcis::util::print_report("Average LCP = ",(double) total_lcp/(name+1),"\n");
        gcis::util::print_report("Average Rule Suffix Length = ",(double) total_rule_suffix_length/(name+1),"\n");
        gcis::util::print_report("Dictionary Level Size (bytes) =",g[level].size_in_bytes(),"\n");
        gcis::util::print_report("LCP Size (bits) = ",g[level].lcp.size(),"\n");
        gcis::util::print_report("Rule Suffix Length (total) = ",g[level].rule.size(),"\n");
        gcis::util::print_report("Rule Suffix Width (bits per symbol) = ",(int_t) g[level].rule.width(),"\n");
        gcis::util::print_report("Tail Length = ",g[level].tail.size(),"\n");
        gcis::util::print_report("Tail Width (bits per symbol) = ",(int_t) g[level].tail.width(),"\n");
        gcis::util::print_report("Run Length Potential (total) = ",run_length_potential,"\n");
        gcis::util::print_report("Avg Run Length per Rule Suffix = ",(double) run_length_potential/(name+1),"\n");
#endif

        bool premature_stop = evaluate_premature_stop(n,K,n1,name+1,level);
        if (name+1 < n1 && !premature_stop) {
            g[level].string_size = n;
            g[level].alphabet_size = K;
            gc_is((int_t *) s1, SA1, n1, name+1, sizeof(int_t), level + 1);
        } else { // generate the suffix array of s1 directly
            if(premature_stop){
#ifdef REPORT
                gcis::util::print_report("Premature Stop employed at level ",level, "\n");
#endif
                reduced_string.resize(n);
                for (j = 0; j < n; j++) {
                    // Copy the reduced substring
                    reduced_string[j] = s[j];
                }
                // Discard the last computed level
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
            gcis::util::print_report("Reduced String Length = ",(int_t) reduced_string.size(),"\n");
            gcis::util::print_report("Reduced String Width (bits per symbol) = ",(int_t) reduced_string.width(),"\n");
#endif
        }
        // Initialize support data structures for LCP
        // and rules delimiters
        for(uint64_t i=0;i<g.size();i++){
            sdsl::util::init_support(g[i].lcp_sel,&(g[i].lcp));
            sdsl::util::init_support(g[i].rule_sel,&(g[i].rule_delim));
        }
        // delete bitvector t
        delete[] t;
    }


    sdsl::int_vector<> gd_is(int level) {
        if(level==g.size()){
            return (reduced_string);
        }
        else{
            //Picks the reduced string from the next level
            sdsl::int_vector<> r_string= gd_is(level+1);
            //TODO: Compress next_r_string beforehand
            sdsl::int_vector<> next_r_string;
            //Decode the reduced string using the level grammar
            next_r_string.resize(g[level].tail.size());

            //Copy the tail to the beggining of the new reduced string
            for(uint64_t i=0;i<g[level].tail.size();i++){
                next_r_string[i]=g[level].tail[i];
            }

            //Expand the rules to each position of the new reduced string
            for(uint64_t i=0;i<r_string.size();i++){
                uint64_t rule = r_string[i];
                g[level].expand_rule(rule,next_r_string);
            }
            //returns the new reduced string to the level i-1
            sdsl::util::bit_compress(next_r_string);
            return next_r_string;
        }
    }



    // compute SA for the S-Type suffixes by inducing the L-Type suffixes and the S-Type suffixes
    void induceSAs(unsigned char *t,
                   uint_t *SA,
                   int_t *s,
                   int_t *bkt,
                   int_t n,
                   int_t K,
                   int cs,
                   int level) {
        int_t i, j;
        getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
        for (i = n - 1; i >= 0; i--) {
            if (SA[i] != EMPTY) {
                j = SA[i] - 1;
                if (j >= 0 && tget(j)) {
                    SA[bkt[chr(j)]--] = j;
                }
            }
        }
    }

    // compute SA for the L-Type suffixes by inducing the LMS-Suffixes and the L-Suffixes
    void induceSAl(unsigned char *t,
                   uint_t *SA,
                   int_t *s,
                   int_t *bkt,
                   int_t n,
                   int_t K,
                   int cs,
                   int level) {
        int_t i, j;
        // find heads of buckets
        getBuckets(s, bkt, n, K, cs, false);
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


// compute the head or end of each bucket
    void getBuckets(int_t *s,
                    int_t *bkt,
                    int_t n,
                    int_t K,
                    int cs,
                    int end) {
        int_t i, sum = 0;

        // clear all buckets
        for (i = 0; i < K; i++) {
            bkt[i] = 0;
        }

        // compute the size of each bucket
        for (i = 0; i < n; i++) {
            bkt[chr(i)]++;
        }
        //Mark the end of each bucket
        for (i = 0; i < K; i++) {
            sum += bkt[i];
            bkt[i] = end ? sum - 1 : sum - bkt[i];
        }
    }

};





#endif //GC_IS_GCIS_UNARY_HPP
