#include "gcis_s8b_codec.hpp"

uint64_t gcis_s8b_codec::size_in_bytes() {
    uint64_t total_bytes = 0;
    total_bytes += 2*sizeof(uint_t);
    total_bytes += lcp.size_in_bytes();
    total_bytes += rule_suffix_length.size_in_bytes();
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    return total_bytes;
}

void gcis_s8b_codec::serialize(std::ostream &o) {
    o.write((char*) &string_size, sizeof(string_size));
    o.write((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.serialize(o);
    rule_suffix_length.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
}
void gcis_s8b_codec::load(std::istream &i) {
    i.read((char*) &string_size, sizeof(string_size));
    i.read((char*) &alphabet_size, sizeof(alphabet_size));
    lcp.load(i);
    rule_suffix_length.load(i);
    rule.load(i);
    tail.load(i);
}


gcis_s8b_codec_level gcis_s8b_codec::decompress(){
    gcis_s8b_codec_level gd;
    uint64_t number_of_rules = rule_suffix_length.size();
    uint64_t total_lcp_length = 0;
    uint64_t total_rule_suffix_length = 0;
    // Compute the total LCP length
    lcp.reset();
    rule_suffix_length.reset();
    for(uint64_t i=0;i<lcp.size();i++){
        total_lcp_length+=lcp.get_next();
    }
    // Compute the total rule suffix length and the number of rules
    for(uint64_t i=0;i<rule_suffix_length.size();i++){
        total_rule_suffix_length += rule_suffix_length.get_next();
    }

    // Resize data structures
    uint64_t total_length = total_lcp_length + total_rule_suffix_length;
    gd.rule.width(sdsl::bits::hi(alphabet_size-1)+1);
    gd.rule.resize(total_length);
    gd.rule_delim.resize(total_length+1);
    sdsl::util::set_to_value(gd.rule_delim,0);


    uint64_t rule_start = 0;
    uint64_t prev_rule_start;
    int64_t last_set_lcp = -1;
    int64_t last_set_rule = -1;
    uint64_t start = 0;
    lcp.reset();
    rule_suffix_length.reset();


    for(uint64_t i=0;i<number_of_rules;i++){
        int64_t k;

        uint64_t lcp_length = lcp.get_next();
        uint64_t rule_length = rule_suffix_length.get_next();

        total_length = lcp_length+rule_length;

        // Copy the contents of the previous rule by LCP length chars
        uint64_t j = 0;
        while(j<lcp_length){
            gd.rule[rule_start+j] = gd.rule[prev_rule_start+j];
            j++;
        }

        // Copy the remaining suffix rule
        while(j<total_length){
            gd.rule[rule_start+j] = rule[start++];
            j++;
        }

        gd.rule_delim[rule_start]=1;
        prev_rule_start = rule_start;
        rule_start+=total_length;
    }
    gd.rule_delim[gd.rule_delim.size()-1]=1;

//    sdsl::util::bit_compress(gd.rule);
    sdsl::util::init_support(gd.rule_delim_sel,&gd.rule_delim);
    return gd;
}



void gcis_s8b_codec_level::expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num+1);
    uint64_t rule_length = rule_delim_sel(rule_num+2) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        r_string[l] = rule[rule_start+i];
        l++;
    }
}

void gcis_s8b_codec_level::expand_rule(uint64_t rule_num, char* s,uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num+1);
    uint64_t rule_length = rule_delim_sel(rule_num+2) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        s[l] = rule[rule_start+i];
        l++;
    }
}
