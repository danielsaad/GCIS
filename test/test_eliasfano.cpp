#include <vector>
#include <cstdint>
#include <random>
#include <ostream>
#include <fstream>
#include "gtest/gtest.h"
#include "eliasfano.hpp"



TEST(eliasfano_simple, int_vector_access){
    sdsl::int_vector<> v = {1,0,0,1,2,0,0,4,4,0};
    sdsl::int_vector<> v2 = v;
    eliasfano_codec ef;
    ef.encode(v);
    for(uint64_t i=0;i<v.size();i++){
        EXPECT_EQ(ef[i],v2[i]);
    }
}


TEST(eliasfano_simple, bitvector_access){
    sdsl::bit_vector bv = {1,0,0,1,1,0,1,0,0,0,1,1,0,0,1};
    eliasfano_codec ef;
    ef.encode(bv);
    EXPECT_EQ(ef.size(),bv.size());
    for(uint64_t i=0;i<bv.size();i++){
        EXPECT_EQ(ef.access_bv(i),bv[i]);
    }
}

TEST(eliasfano_simple, bitvector_sel){
    sdsl::bit_vector bv = {1,0,0,1,1,0,1,0,0,0,1,1,0,0,1};
    eliasfano_codec ef;
    ef.encode(bv);
    EXPECT_EQ(ef.size(),bv.size());
    uint64_t k = 0;
    for(uint64_t i=0;i<bv.size();i++){
        if(bv[i]) {
            EXPECT_EQ(ef.pos(k),i);
            k++;
        }
    }
}

TEST(eliasfano_simple, bitvector_value){
    sdsl::bit_vector bv = {1,0,0,1,1,0,1,0,0,0,1,1,0,0,1};
    sdsl::int_vector<> v = {0,2,0,1,3,0,2};
    eliasfano_codec ef;
    ef.encode(bv);
    EXPECT_EQ(ef.size(),bv.size());
    for(uint64_t i=0;i<v.size();i++){
        EXPECT_EQ(ef[i],v[i]);
    }
}



TEST(eliasfano_serialize,bitvector_access){
    std::ofstream o("/tmp/teste.dat");
    sdsl::bit_vector bv = {1,0,0,1,1,0,0,1,0,1,0,1,0,1};
    eliasfano_codec ef;
    ef.encode(bv);
    ef.serialize(o);
    eliasfano_codec ef2;
    o.close();
    std::ifstream i("/tmp/teste.dat");
    ef2.load(i);
    for(uint64_t i=0;i<bv.size();i++){
        EXPECT_EQ(ef2.access_bv(i),bv[i]);
    }
    i.close();
}

TEST(eliasfano_serialize,bitvector_sel){
    std::ofstream o("/tmp/teste.dat");
    sdsl::bit_vector bv = {1,0,0,1,1,0,0,1,0,1,0,1,0,1};
    eliasfano_codec ef;
    ef.encode(bv);
    ef.serialize(o);
    eliasfano_codec ef2;
    o.close();
    std::ifstream i("/tmp/teste.dat");
    ef2.load(i);
    uint64_t k = 0;
    for(uint64_t i=0;i<bv.size();i++){
        if(bv[i]) {
            EXPECT_EQ(ef2.pos(k),i);
            k++;
        }
    }
    i.close();
}
