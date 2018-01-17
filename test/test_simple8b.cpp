#include <vector>
#include <cstdint>
#include <random>
#include <ostream>
#include <fstream>
#include "gtest/gtest.h"
#include "simple8b.hpp"



TEST(simple8b_exaustive, random){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, 1<<10);
    std::vector<uint64_t> v;
    simple8b_codec s8;
    for(uint64_t i=0;i<(1<<20);i++){
        uint64_t n = dis(gen);
        v.push_back(n);
        s8.encode(n);
    }
    s8.encode();
    for(uint64_t i=0;i<s8.size();i++){
        EXPECT_EQ(v[i],s8.get_next());
    }
//    s8.reset();
//    for(uint64_t i=0;i<s8.size();i++){
//        EXPECT_EQ(v[i],s8.get_next());
//    }
}


TEST(simplet8b_serialize_load,random){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, 1<<10);
    std::vector<uint64_t> v;
    simple8b_codec s8;
    std::ofstream o;
    o.open("teste.bin",std::ofstream::binary);
    for(uint64_t i=0;i<(1<<20);i++){
        uint64_t n = dis(gen);
        v.push_back(n);
        s8.encode(n);
    }
    s8.encode();
    s8.serialize(o);
    simple8b_codec s8_prime;
    std::ifstream i;
    i.open("teste.bin",std::ifstream::binary);
    s8_prime.load(i);
    for(uint64_t i=0;i<(1<<20);i++) {
        EXPECT_EQ(v[i], s8_prime.get_next());
    }
}