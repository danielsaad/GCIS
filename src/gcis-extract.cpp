#include <iostream>
#include "gc-is.hpp"
int main(int argc,char* argv[]){
    std::ifstream input (argv[1]);
    std::ifstream query(argv[2]);
    gc_is_dictionary<gcis_eliasfano_codec> d;
    d.load(input);
    uint64_t l,r;
    while(query >> l >> r){
        d.extract(l,r);
    }
    return 0;
}