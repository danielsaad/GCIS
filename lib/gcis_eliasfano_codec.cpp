#include "gcis_eliasfano_codec.hpp"


uint64_t gcis_eliasfano_codec::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += 2*sizeof(uint_t);
    total_bytes += lcp.size_in_bytes();
    total_bytes += rule_suffix_length.size_in_bytes();
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    // total_bytes += sdsl::size_in_bytes(lms_bv);
    // total_bytes += sdsl::size_in_bytes(lms_rnk1);
    // total_bytes += sdsl::size_in_bytes(lms_sel1);
    return total_bytes;
}

void gcis_eliasfano_codec::serialize(std::ostream &o) {

    o.write((char*) &string_size, sizeof(string_size));
    o.write((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.serialize(o);
    rule_suffix_length.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
    // lms_bv.serialize(o);
    // lms_rnk1.serialize(o);
    // lms_sel1.serialize(o);
}
void gcis_eliasfano_codec::load(std::istream &i) {
    i.read((char*) &string_size, sizeof(string_size));
    i.read((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.load(i);
    rule_suffix_length.load(i);
    rule.load(i);
    tail.load(i);
    // lms_bv.load(i);
    // lms_rnk1.load(i);
    // lms_sel1.load(i);
    // lms_rnk1.set_vector(&lms_bv);
    // lms_sel1.set_vector(&lms_bv);
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

uint64_t gcis_eliasfano_codec::get_lcp(uint64_t i) {
    return lcp[i];
}

uint64_t gcis_eliasfano_codec::get_rule_pos(uint64_t i) {
    return i == 0 ? 0 : rule_suffix_length.pos(i-1) -(i-1);
}

uint64_t gcis_eliasfano_codec::get_rule_length(uint64_t i){
    return rule_suffix_length[i];

}

// // Extract the LCP part L from a rule_num regarding L[l,r]
// // If l>r, does not extract
// void gcis_eliasfano_codec::extract_lcp(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
//     if(l>r) return;
//     uint64_t lcp_len = get_lcp(rule_num);
//     uint64_t cur_lcp_len = lcp_len;
//     // Index of LCP element
//     int64_t idx = k+ lcp_len-1;
//     uint64_t prev_rule = rule_num - 1;
//     while(idx >= (int64_t) k+l){
//         uint64_t prev_lcp_len = get_lcp(prev_rule);
//         if(prev_lcp_len<cur_lcp_len){
//             int64_t i = cur_lcp_len -1;
//             uint64_t prev_rule_pos = get_rule_pos(prev_rule);
//             while(i >= (int64_t)prev_lcp_len && idx >= (int64_t) k+l){
//                 // We copy the value to the extracted text iff it falls in the interval
//                 if(idx<= k+r){
//                     extracted_text[idx-l] = rule[prev_rule_pos+i-prev_lcp_len];
//                 }
//                 i--;
//                 idx--;
//             }
//             cur_lcp_len = prev_lcp_len;
//         }
//         prev_rule--;
//     }
//     k+=(r-l+1);
// }
//
// // Extract the Rule Suffix part S from a rule_num regarding S[l,r]
// // If l>r, does not extract
// void gcis_eliasfano_codec::extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
//     uint64_t rule_pos = get_rule_pos(rule_num);
//     for(int64_t i = l;i<= r;i++){
//         extracted_text[k++] = rule[rule_pos+i];
//     }
// }
//
// // Extract the Rule R regarding R[l,r]
// void gcis_eliasfano_codec::extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
//     int64_t lcp_len = get_lcp(rule_num);
//     int64_t rule_suffix_len = get_rule_length(rule_num);
//     int64_t l_lcp,l_rule_suffix,r_lcp,r_rule_suffix;
//
//     if(lcp_len>l) {
//         //There is something in the LCP to extract
//
//         //LCP starts at l, discarding previous LMS symbols
//         l_lcp = l;
//         //LCP ends prematurely at position r or copies all lcp symbols from L[l,r]
//         r_lcp = std::min<int64_t>(lcp_len-1, r);
//         extract_lcp(rule_num, l_lcp, r_lcp, extracted_text, k);
//         if(r_lcp==lcp_len-1){
//             l_rule_suffix = 0;
//             r_rule_suffix = r - lcp_len;
//             extract_rule_suffix(rule_num,l_rule_suffix,r_rule_suffix,extracted_text,k);
//         }
//     }
//     else{
//         // There is no LCP to extract. Some rule symbols must be discarded as well
//         l_rule_suffix = l-lcp_len;
//         r_rule_suffix = r - lcp_len;
//         extract_rule_suffix(rule_num, l_rule_suffix, r_rule_suffix, extracted_text, k);
//     }
//
// }
//
// void gcis_eliasfano_codec::extract_rules(uint64_t l,
//                                          uint64_t r,
//                                          sdsl::int_vector<>& extracted_text,
//                                          sdsl::int_vector<>& tmp_text) {
//     // The subinterval to be extracted falls in the tail
//     uint64_t k = 0;
//     if(l<tail.size()) {
//         for (uint64_t i = l; i < std::min(r+1, tail.size()); i++) {
//             extracted_text[k++] = tail[i];
//         }
//         if (r < tail.size()) {
//             return;
//         }
//         l = tail.size();
//     }
//     //extract rules from lms substrings in the interval[l,r]
//     uint64_t leftmost_rule = lms_rnk1(l+1)-1;
//     uint64_t rightmost_rule = lms_rnk1(r+1)-1;
//     // Beggining of the leftmost LMS substring
//     uint64_t left_pos = lms_sel1(leftmost_rule+1);
//     // Beginning of the rightmost LMS substring
//     uint64_t right_pos = lms_sel1(rightmost_rule+1);
//     // Extract a suffix of the leftmost lms substring
//     if(leftmost_rule!=rightmost_rule) {
//         uint64_t leftmost_rule_len = get_lcp(tmp_text[0]) + get_rule_length(tmp_text[0]);
//         extract_rule(tmp_text[0],l-left_pos,leftmost_rule_len-1,extracted_text,k);
//         for (uint64_t i = leftmost_rule + 1; i < rightmost_rule; i++) {
//             // Extract the entire lms substring
//             uint64_t rule_len = get_lcp(tmp_text[i-leftmost_rule]) + get_rule_length(tmp_text[i - leftmost_rule]);
//             extract_rule(tmp_text[i - leftmost_rule],0,rule_len-1, extracted_text, k);
//         }
//         extract_rule(tmp_text[rightmost_rule-leftmost_rule],0,r-right_pos,extracted_text,k);
//     }
//     else{
//         extract_rule(tmp_text[0],l-left_pos,r-right_pos, extracted_text, k);
//     }


// }
