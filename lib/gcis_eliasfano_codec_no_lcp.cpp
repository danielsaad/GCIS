#include "gcis_eliasfano_codec_no_lcp.hpp"


uint64_t gcis_eliasfano_codec_no_lcp::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += 2*sizeof(uint_t);
    total_bytes += rule_suffix_length.size_in_bytes();
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    return total_bytes;
}

void gcis_eliasfano_codec_no_lcp::serialize(std::ostream &o) {

    o.write((char*) &string_size, sizeof(string_size));
    o.write((char*) &alphabet_size, sizeof(alphabet_size));
    rule_suffix_length.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
    fully_decoded_rule_len.serialize(o);
    o.write((char*) &fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}
void gcis_eliasfano_codec_no_lcp::load(std::istream &i) {
    i.read((char*) &string_size, sizeof(string_size));
    i.read((char*) &alphabet_size, sizeof(alphabet_size));
    rule_suffix_length.load(i);
    rule.load(i);
    tail.load(i);
    fully_decoded_rule_len.load(i);
    i.read((char*) &fully_decoded_tail_len, sizeof(fully_decoded_tail_len));
}


gcis_eliasfano_codec_no_lcp_level gcis_eliasfano_codec_no_lcp::decompress(){
    gcis_eliasfano_codec_no_lcp_level gd;
    uint64_t number_of_rules = 0;
    for(uint64_t i=0;i<rule_suffix_length.size();i++){
        if(rule_suffix_length.access_bv(i)){
            number_of_rules++;
        }
    }
    uint64_t total_rule_suffix_length = 0;
    // Compute the total rule suffix length and the number of rules
    for(uint64_t i=0;i<number_of_rules;i++){
        total_rule_suffix_length += rule_suffix_length[i];
    }

    // Resize data structures
    uint64_t total_length = total_rule_suffix_length;
    sdsl::bit_vector rule_delim(total_length+1);
    gd.rule.resize(total_length);
    gd.rule.width(sdsl::bits::hi(alphabet_size-1)+1);
    sdsl::util::set_to_value(rule_delim,0);


    uint64_t rule_start = 0;
    uint64_t prev_rule_start;
    uint64_t start = 0;

    for(uint64_t i=0;i<number_of_rules;i++){

        uint64_t rule_length = rule_suffix_length[i];
        total_length = rule_length;

        uint64_t j = 0;

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



void gcis_eliasfano_codec_no_lcp_level::expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        r_string[l] = rule[rule_start+i];
        l++;
    }
}

void gcis_eliasfano_codec_no_lcp_level::expand_rule(uint64_t rule_num, char* s,uint64_t &l) {
    uint64_t rule_start = rule_delim.pos(rule_num);
    uint64_t rule_length = rule_delim.pos(rule_num+1) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        s[l] = rule[rule_start+i];
        l++;
    }
}

uint64_t gcis_eliasfano_codec_no_lcp::get_rule_pos(uint64_t i) {
    return i == 0 ? 0 : rule_suffix_length.pos(i-1) -(i-1);
}


uint64_t gcis_eliasfano_codec_no_lcp::get_rule_length(uint64_t i){
    return rule_suffix_length[i];

}


 // Extract the Rule Suffix part S from a rule_num regarding S[l,r]
 // If l>r, does not extract
 void gcis_eliasfano_codec_no_lcp::extract_rule_suffix(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
     uint64_t rule_pos = get_rule_pos(rule_num);
     for(int64_t i = l;i<= r;i++){
         extracted_text[k++] = rule[rule_pos+i];
     }
 }
//
 // Extract the Rule R regarding R[l,r]
 void gcis_eliasfano_codec_no_lcp::extract_rule(uint64_t rule_num,int64_t l,int64_t r,sdsl::int_vector<>& extracted_text,uint64_t& k){
    int64_t rule_suffix_len = get_rule_length(rule_num);
    extract_rule_suffix(rule_num, l,r, extracted_text, k);
 }


/**
 * Extract a rule R entirely
 * @param rule_num The number of the rule to be extracted
 * @param extracted_text the container for the extracted rule
 * @param k The index which tracks in which position of the container the extracted rule should be placed
 */
void gcis_eliasfano_codec_no_lcp::extract_rule(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k){
    extract_rule_suffix(rule_num, extracted_text, k);
}



void gcis_eliasfano_codec_no_lcp::extract_rule_suffix(uint64_t rule_num,sdsl::int_vector<>& extracted_text,uint64_t& k){
    int64_t rule_suffix_len = get_rule_length(rule_num);
    uint64_t rule_pos = get_rule_pos(rule_num);
    for(int64_t i = 0;i< rule_suffix_len;i++){
        extracted_text[k++] = rule[rule_pos+i];
    }
}