#include "gcis_unary_codec.hpp"



void gcis_unary_codec::serialize(std::ostream &o) {
    o.write((char*) &alphabet_size,sizeof(alphabet_size));
    o.write((char*) &string_size,sizeof(string_size));
    lcp.serialize(o);
    rule_delim.serialize(o);
    rule.serialize(o);
    tail.serialize(o);
    rule_sel.serialize(o);
    lcp_sel.serialize(o);
}

void gcis_unary_codec::load(std::istream& i){
    i.read((char*) &alphabet_size,sizeof(alphabet_size));
    i.read((char*) &string_size,sizeof(string_size));
    lcp.load(i);
    rule_delim.load(i);
    rule.load(i);
    tail.load(i);
    rule_sel.load(i);
    lcp_sel.load(i);
    lcp_sel.set_vector(&lcp);
    rule_sel.set_vector(&rule_delim);
}



uint64_t gcis_unary_codec::size_in_bytes() {
    uint64_t total_bytes = sdsl::size_in_bytes(lcp);
    total_bytes += sdsl::size_in_bytes(rule_delim);
    total_bytes += sdsl::size_in_bytes(lcp_sel);
    total_bytes += sdsl::size_in_bytes(rule_sel);
    total_bytes += sdsl::size_in_bytes(rule);
    total_bytes += sdsl::size_in_bytes(tail);
    total_bytes += sizeof(alphabet_size);
    total_bytes += sizeof(string_size);
    return total_bytes;
}

inline uint64_t gcis_unary_codec::get_lcp_length(uint64_t rule_num) {
    if(rule_num==0)
        return 0;
    return lcp_sel(rule_num+1)-lcp_sel(rule_num)-1;
}

inline uint64_t gcis_unary_codec::get_s_length(uint64_t rule_num) {
    if(rule_num==0)
        return 1;
    return(rule_sel(rule_num+1)-rule_sel(rule_num)-1);
}

inline uint64_t gcis_unary_codec::get_s_start(uint64_t rule_num) {
    if(rule_num==0)
        return 0;
    return(rule_sel(rule_num)-rule_num+1);
}
//TODO: we dont need to perform select to recover the LCP entries.
gcis_unary_codec_level gcis_unary_codec::decompress(){
    gcis_unary_codec_level gd;
    uint64_t number_of_rules = 0;
    uint64_t lcp_length = 0;
    uint64_t rule_suffix_length = 0;
    // Compute the total LCP length
    for(uint64_t i=0;i<lcp.size();i++){
        if(lcp[i]==0)
            lcp_length++;
    }
    // Compute the total rule suffix length and the number of rules
    for(uint64_t i=0;i<rule_delim.size();i++){
        if(rule_delim[i]==0){
            rule_suffix_length++;
        }
        else{
            number_of_rules++;
        }
    }

    // Resize data structures
    uint64_t total_length = lcp_length + rule_suffix_length;
    gd.rule.width(sdsl::bits::hi(alphabet_size-1)+1);
    gd.rule.resize(total_length);
    gd.rule_delim.resize(total_length+1);
    sdsl::util::set_to_value(gd.rule_delim,0);


    uint64_t rule_start = 0;
    uint64_t prev_rule_start;
    int64_t last_set_lcp = -1;
    int64_t last_set_rule = -1;
    uint64_t start = 0;

    for(uint64_t i=0;i<number_of_rules;i++){
        int64_t k;
        // Computes LCP length of the ith rule
        for(k=last_set_lcp+1;lcp[k]==0 && k<lcp.size();k++);

        uint64_t lcp_length = k-last_set_lcp-1;
        last_set_lcp = k;

        // Computes ith Rule Suffix Length
        for(k=last_set_rule+1;rule_delim[k]==0 && k<rule_delim.size();k++);
        uint64_t rule_length = k - last_set_rule -1;
        last_set_rule = k;

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

void gcis_unary_codec_level::expand_rule(uint64_t rule_num, sdsl::int_vector<> &r_string, uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num+1);
    uint64_t rule_length = rule_delim_sel(rule_num+2) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        r_string[l] = rule[rule_start+i];
        l++;
    }
}

void gcis_unary_codec_level::expand_rule(uint64_t rule_num, char* s,uint64_t &l) {
    uint64_t rule_start = rule_delim_sel(rule_num+1);
    uint64_t rule_length = rule_delim_sel(rule_num+2) - rule_start;
    for(uint64_t i=0;i<rule_length;i++){
        s[l] = rule[rule_start+i];
        l++;
    }
}





void gcis_unary_codec::expand_rule(uint64_t rule_num,
                                sdsl::int_vector<>& r_string){
    uint64_t lcp_length = get_lcp_length(rule_num);
    uint64_t s_length = get_s_length(rule_num);
    uint64_t rule_length  = lcp_length + s_length;
    uint64_t rule_start = get_s_start(rule_num);
    uint64_t l = r_string.size();
    r_string.resize(r_string.size()+rule_length);

    if(lcp_length>0) {
        expand_lcp(rule_num - 1,lcp_length,r_string,l);
    }
    for(uint64_t i=lcp_length;i<rule_length;i++){
        r_string[l+i] = rule[rule_start++];
    }
}

void gcis_unary_codec::expand_lcp(uint64_t rule_num,
                               uint64_t lcp_length,
                               sdsl::int_vector<>& r_string,
                               uint64_t l){
    uint64_t new_lcp_length = get_lcp_length(rule_num);
    if(new_lcp_length>0){
        expand_lcp(rule_num-1,std::min(lcp_length,new_lcp_length),r_string,l);
    }
    if(lcp_length>new_lcp_length){
        uint64_t rule_str_start = get_s_start(rule_num);
        for(uint64_t i= new_lcp_length;i<lcp_length;i++){
            r_string[l+i] = rule[rule_str_start];
            rule_str_start++;
        }
    }

}
