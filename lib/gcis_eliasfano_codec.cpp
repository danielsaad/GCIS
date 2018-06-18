#include "gcis_eliasfano_codec.hpp"


uint64_t gcis_eliasfano_codec::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += 2*sizeof(uint_t);
    total_bytes += lcp.size_in_bytes();
    total_bytes += rule_suffix_length.size_in_bytes();
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    return total_bytes;
}

void gcis_eliasfano_codec::serialize(std::ostream &o) {

    o.write((char*) &string_size, sizeof(string_size));
    o.write((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.serialize(o);
    rule_suffix_length.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
    fully_decoded_rule_len.serialize(o);
    o.write((char*) &fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}
void gcis_eliasfano_codec::load(std::istream &i) {
    i.read((char*) &string_size, sizeof(string_size));
    i.read((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.load(i);
    rule_suffix_length.load(i);
    rule.load(i);
    tail.load(i);
    fully_decoded_rule_len.load(i);
    i.read((char*) &fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}


gcis_eliasfano_codec_level gcis_eliasfano_codec::decompress(){
    gcis_eliasfano_codec_level gd;
    uint64_t number_of_rules = 0;
    for(uint64_t i=0;i<rule_suffix_length.size();i++){
        if(rule_suffix_length.access_bv(i)){
            number_of_rules++;
        }
    }
    uint64_t total_lcp_length = 0;
    uint64_t total_rule_suffix_length = 0;
    // Compute the total LCP length
    for(uint64_t i=0;i<number_of_rules;i++){
        total_lcp_length+=lcp[i];
    }
    // Compute the total rule suffix length and the number of rules
    for(uint64_t i=0;i<number_of_rules;i++){
        total_rule_suffix_length += rule_suffix_length[i];
    }

    // Resize data structures
    uint64_t total_length = total_lcp_length + total_rule_suffix_length;
    sdsl::bit_vector rule_delim(total_length+1);
    gd.rule.resize(total_length);
    gd.rule.width(sdsl::bits::hi(alphabet_size-1)+1);
    sdsl::util::set_to_value(rule_delim,0);


    uint64_t rule_start = 0;
    uint64_t prev_rule_start;
    uint64_t start = 0;

    for(uint64_t i=0;i<number_of_rules;i++){

        uint64_t lcp_length = lcp[i];
        uint64_t rule_length = rule_suffix_length[i];
        total_length = lcp_length+rule_length;

        // Copy the contents of the previous rule by LCP length chars
        uint64_t j = 0;
        while(j<lcp_length){
            gd.rule[rule_start+j] = gd.rule[prev_rule_start+j];
            j++;
        }

        // Copy the remaining suffix rule
        while(j<total_length){
            assert(start < rule.size());
            gd.rule[rule_start+j] = rule[start++];
            j++;
        }

        rule_delim[rule_start]=1;
        prev_rule_start = rule_start;
        rule_start+=total_length;
    }
    rule_delim[rule_delim.size()-1]=1;
    gd.rule_delim.encode(rule_delim);

    return gd;
}

void gcis_eliasfano_codec_level::expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        r_string[l] = rule[rule_start+i];
        l++;
    }
}

void gcis_eliasfano_codec_level::expand_rule(uint64_t rule_num, char* s,uint64_t &l) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        s[l] = rule[rule_start+i];
        l++;
    }
}

void gcis_eliasfano_codec_level::expand_rule_bkt(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l, int_t *bkt) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        r_string[l++] = rule[rule_start+i];
	bkt[rule[rule_start+i]]++;
    }
}

void gcis_eliasfano_codec_level::expand_rule_bkt(uint64_t rule_num, char* s,uint64_t &l, int_t *bkt) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        s[l++] = rule[rule_start+i];
	bkt[rule[rule_start+i]]++;
    }
}

uint64_t gcis_eliasfano_codec::get_lcp(uint64_t i) {
    return lcp[i];
}

uint64_t gcis_eliasfano_codec::get_rule_pos(uint64_t i) {
    return i == 0 ? 0 : rule_suffix_length.pos(i-1) -(i-1);
}


uint64_t gcis_eliasfano_codec::get_rule_length(uint64_t i){
    return rule_suffix_length[i];

}

 // Extract the LCP part L from a rule_num regarding L[l,r]
 // If l>r, does not extract
 void gcis_eliasfano_codec::extract_lcp(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
     if(l>r) return;
     uint64_t lcp_len = get_lcp(rule_num);
     uint64_t cur_lcp_len = lcp_len;
     // Index of LCP element
     int64_t idx = k+ lcp_len-1;
     uint64_t prev_rule = rule_num - 1;
     while(idx >= (int64_t) k+l){
         uint64_t prev_lcp_len = get_lcp(prev_rule);
         if(prev_lcp_len<cur_lcp_len){
             int64_t i = cur_lcp_len -1;
             uint64_t prev_rule_pos = get_rule_pos(prev_rule);
             while(i >= (int64_t)prev_lcp_len && idx >= (int64_t) k+l){
                 // We copy the value to the extracted text iff it falls in the interval
                 if(idx<= k+r){
                     extracted_text[idx-l] = rule[prev_rule_pos+i-prev_lcp_len];
                 }
                 i--;
                 idx--;
             }
             cur_lcp_len = prev_lcp_len;
         }
         prev_rule--;
     }
     k+=(r-l+1);
 }

 // Extract the Rule Suffix part S from a rule_num regarding S[l,r]
 // If l>r, does not extract
 void gcis_eliasfano_codec::extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
     uint64_t rule_pos = get_rule_pos(rule_num);
     for(int64_t i = l;i<= r;i++){
         extracted_text[k++] = rule[rule_pos+i];
     }
 }
//
 // Extract the Rule R regarding R[l,r]
 void gcis_eliasfano_codec::extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
     int64_t lcp_len = get_lcp(rule_num);
     int64_t rule_suffix_len = get_rule_length(rule_num);
     int64_t l_lcp,l_rule_suffix,r_lcp,r_rule_suffix;

     if(lcp_len>l) {
         //There is something in the LCP to extract

         //LCP starts at l, discarding previous LMS symbols
         l_lcp = l;
         //LCP ends prematurely at position r or copies all lcp symbols from L[l,r]
         r_lcp = std::min<int64_t>(lcp_len-1, r);
         extract_lcp(rule_num, l_lcp, r_lcp, extracted_text, k);
         if(r_lcp==lcp_len-1){
             l_rule_suffix = 0;
             r_rule_suffix = r - lcp_len;
             extract_rule_suffix(rule_num,l_rule_suffix,r_rule_suffix,extracted_text,k);
         }
     }
     else{
         // There is no LCP to extract. Some rule symbols must be discarded as well
         l_rule_suffix = l-lcp_len;
         r_rule_suffix = r - lcp_len;
         extract_rule_suffix(rule_num, l_rule_suffix, r_rule_suffix, extracted_text, k);
     }

 }


/**
 * Extract a rule R entirely
 * @param rule_num The number of the rule to be extracted
 * @param extracted_text the container for the extracted rule
 * @param k The index which tracks in which position of the container the extracted rule should be placed
 */
void gcis_eliasfano_codec::extract_rule(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k){
    extract_lcp(rule_num, extracted_text, k);
    extract_rule_suffix(rule_num, extracted_text, k);

}


/**
 *
 * @param rule_num
 * @param extracted_text
 * @param k
 */
void gcis_eliasfano_codec::extract_lcp(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k){
    uint64_t lcp_len = get_lcp(rule_num);
    uint64_t cur_lcp_len = lcp_len;
    // Index of LCP element
    int64_t idx = k+ lcp_len-1;
    uint64_t prev_rule = rule_num - 1;
    while(idx >= (int64_t) k){
        uint64_t prev_lcp_len = get_lcp(prev_rule);
        if(prev_lcp_len<cur_lcp_len){
            int64_t i = cur_lcp_len -1;
            uint64_t prev_rule_pos = get_rule_pos(prev_rule);
            while(i >= (int64_t) prev_lcp_len && idx >= (int64_t) k){
                // We copy the value to the extracted text iff it falls in the interval
                if(idx<= k+lcp_len-1){
                    extracted_text[idx] = rule[prev_rule_pos+i-prev_lcp_len];
                }
                i--;
                idx--;
            }
            cur_lcp_len = prev_lcp_len;
        }
        prev_rule--;
    }
    k += lcp_len;
}


void gcis_eliasfano_codec::extract_rule_suffix(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k){
    int64_t rule_suffix_len = get_rule_length(rule_num);
    uint64_t rule_pos = get_rule_pos(rule_num);
    for(int64_t i = 0;i< rule_suffix_len;i++){
        extracted_text[k++] = rule[rule_pos+i];
    }
}
