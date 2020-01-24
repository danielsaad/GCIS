
#include <ostream>
#include <fstream>
#include <sdsl/lcp_bitcompressed.hpp>
#include "gtest/gtest.h"
#include "gcis_index_bs.hpp"

using namespace gcis_index_private;

gcis_index_bs<> GCIS_index_toy;


std::string T = "AGCTTTTCATTCTGACTGCAACAGCTTTTCATTCTGACTGCAAC";

std::map<char, uint32_t> sigma = {
        {'A',0},
        {'C',1},
        {'G',10},
        {'T',13},
};

std::vector<gcis_index_grid<>::lpoint> points = {
        {{ 6,1}, 43},
        {{ 7,2}, 32},
        {{ 1,3}, 35},
        {{13,4}, 26},
        {{13,5}, 42},
        {{ 3,6}, 37},
        {{ 10,7}, 39},
        {{ 9,8}, 17},
        {{14,9}, 16},
        {{11,10}, 31},
        {{ 1,11}, 36},
        {{15,12}, 22},
        {{15,13}, 41},
        {{ 1,14}, 28},
        {{12,15}, 10},
        {{12,16}, 38},
        {{ 1,17}, 9},
        {{14,18}, 25},
        {{14,19}, 30},
        {{14,20}, 21},
        {{14,21}, 15},
        {{ 2,22}, 24},
        {{ 2,23}, 29},
        {{ 1,24}, 20},
        {{14,25}, 14},
        {{14,26}, 13},
        {{ 2,27}, 12}
};
std::map<uint , std::string > suffixes = {
        {1,"AAC"},
        {2,"AC"},
        {3,"ACTGC"},
        {4,"ACTGCAAC"},
        {5,"AGCTTTTC"},
        {6,"C"},
        {7,"C"},
        {8,"C"},
        {9,"CTGACTGC"},
        {10,"CTGACTGCAAC"},
        {11,"CTGC"},
        {12,"CTTTTC"},
        {13,"CTTTTC"},
        {14,"G"},
        {15,"G"},
        {16,"GC"},
        {17,"T"},
        {18,"TC"},
        {19,"TG"},
        {20,"TGC"},
        {21,"TT"},
        {22,"TTC"},
        {23,"TTTC"},
        {24,"TTTTC"}
};
std::map<uint32_t ,std::string > non_t_exp = {
        /*$A****/{0,"A"},
        /*$C****/{1,"C"},
        /*X1****/{2,"AAC"},
        /*Y2****/{3,"ATTCTGACTGCAAC"},
        /*S*****/{4,"AGCTTTTCATTCTGACTGCAACAGCTTTTCATTCTGACTGCAAC"},
        /*X2****/{5,"ACTGC"},
        /*Y3****/{6,"ATTCTGACTGC"},
        /*X6****/{7,"CTTTTC"},
        /*tail2*/{8,"AGCTTTTC"},
        /*Y1****/{9,"AACAGCTTTTC"},
        /*$G****/{10,"G"},
        /*X3****/{11,"AG"},
        /*X5****/{12,"CTG"},
        /*$T****/{13,"T"},
        /*X4****/{14,"ATT"}

};

std::map<uint32_t ,uint > rule_pre = {
        /*$A****/{0,2},
        /*$C****/{1,3},
        /*S*****/{2,33},
        /*X1****/{3,39},
        /*Y2****/{4,1},
        /*X2****/{5,26},
        /*Y3****/{6,17},
        /*X6****/{7,10},
        /*tail2*/{8,6},
        /*Y1****/{9,32},
        /*$G****/{10,4},
        /*X3****/{11,7},
        /*X5****/{12,22},
        /*$T****/{13,5},
        /*X4****/{14,18}
};



sdsl::bit_vector dfuds = {
        1,1,0, //DFUDS INIT
        1,1,1,1,1,1,1,1,0,//4
        0,0,0,0,1,1,0,//AGCT 8
        1,1,0,//11
        0,0,//AG
        1,1,1,1,1,1,0,//7
        0,0,0,0,0,0,//CTTTTC
        1,1,1,0,//6
        1,1,1,0,//14
        0,0,0,//ATT
        1,1,1,0,//12
        0,0,0,//CTG
        1,1,1,1,1,0,//5
        0,0,0,0,0,//ACTGC
        1,1,1,0,//9
        1,1,1,0,//2
        0,0,0,//AAC
        0,0,//11 7
        1,1,1,1,0,//3
        0,0,0,0// 14 12 5 2
};
std::map<uint ,long> off ={
        {3,-4},//4
        {12,-4},{13,-3},{14,-2},{15,-1},{16,0}, //AGCT 8
        {19,0},//11
        {22,0},{23,1},//AG
        {24,2},//7
        {31,2},{32,3},{33,4},{34,5},{35,6},{36,7},//CTTTTC
        {37,8},//6
        {41,8},//14
        {45,8},{46,9},{47,10}, //ATT
        {48,11},//12
        {52,11},{53,12},{54,13},//CTG
        {55,14},//5
        {61,14},{62,15},{63,16},{64,17},{65,18},//ACTGC
        {66,19},//9
        {70,19},//2
        {74,19},{75,20},{76,21},//AAC
        {77,22},{78,24},//11 7
        {79,30},//3
        {84,30},{85,33},{86,36},{87,41}
};

///////////////////////////S//////////////t2///X3////////X6///////////////////Y3///X4///////////X5///////////X2////////////////Y1//X1////////////X3//X6//Y2//X4//X5//X2//X1
sdsl::int_vector<> X   = {  4, 0, 1,10,13,  8, 11, 0,10,  7, 1,13,13,13,13, 1,  6, 14, 0,13,13, 12, 1,13,10,  5, 0, 1,13,10, 1,  9,  2, 0, 0, 1, 11,  7,  3, 14, 12,  5,  2};
sdsl::bit_vector focc  = {  1, 1, 1, 1, 1,  1,  1, 0, 0,  1, 0, 0, 0, 0, 0, 0,  1,  1, 0, 0, 0,  1, 0, 0, 0,  1, 0, 0, 0, 0, 0,  1,  1, 0, 0, 0,  0,  0,  1,  0,  0,  0,  0};
sdsl::int_vector<> pi  = {  4, 0, 1,10,13,  8, 11,        7,                    6, 14,          12,           5,                 9,  2,                   3};
sdsl::bit_vector t     = {                         1, 1,     1, 1, 1, 1, 1, 1,         1, 1, 1,     1, 1, 1,     1, 1, 1, 1, 1,         1, 1, 1,  0,  0,      0,  0,  0,  0};
std::string st         = "AGCTTTTCATTCTGACTGCAAC";
sdsl::int_vector<>wt   = {11,7,14,12,5,2};
//////////////////////////////////////////////1///////////////////2///////////////////3///////////////////4//////////////
//////////////////////////0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7
//////////////////////////A,C,G,T,A,G,C,T,T,T,T,C,A,T,T,C,T,G,A,C,T,G,C,A,A,C,A,G,C,T,T,T,T,C,A,T,T,C,T,G,A,C,T,G,C,A,A,C
sdsl::bit_vector l     = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0};

void load_toy_example(){


    sdsl::inv_perm_support<> inv_pi (&pi);
    sdsl::wt_gmr<> gmrwt;
    sdsl::construct_im(gmrwt,wt);

    gcis_index_grid<> m_grid (points, 15, 24);
    GCIS_index_toy.set_grid(m_grid);
    GCIS_index_toy.set_bvfocc(sdsl::sd_vector<>(focc));
    GCIS_index_toy.set_tree(dfuds);
    GCIS_index_toy.set_bvt(sdsl::sd_vector<>(t));
    GCIS_index_toy.set_vt(st);
    GCIS_index_toy.set_nt(gmrwt);
    GCIS_index_toy.set_pi(pi);
    GCIS_index_toy.set_l(l);
    GCIS_index_toy.set_sigma(4);

}


TEST(dumbSuite, dumb){
    ASSERT_TRUE(true);
}

//
//TEST(gcisIndexBSToyExample,lowerBound)
//{
//    for (int i = 2; i <= 100 ; ++i) {
//        uint64_t l = 1, r = 100;
//        bool f = GCIS_index_toy.bsearch_lowerBound(l,r,[&i](const uint64_t& X) -> int{
//            return (i == X)?0:( i > X )?1:-1;
//        });
//        ASSERT_EQ(l,i);
//        ASSERT_FALSE(!f);
//     }
//
//}
//TEST(gcisIndexBSToyExample,upperBound)
//{
//    for (int i = 2; i <= 100 ; ++i) {
//        uint64_t l = 1, r = 100;
//        bool f = GCIS_index_toy.bsearch_upperBound(l,r,[&i](const uint64_t& X) -> int{
//            return (i == X)?0:( i > X )?1:-1;
//        });
//        ASSERT_EQ(r,i);
//        ASSERT_FALSE(!f);
//    }
//
//}
//
//TEST(gcisIndexToyExampleBS,findRange){
//
//    load_toy_example();
//
//    for (const auto &sfx : suffixes) {
////        std::cout<<sfx.second<<std::endl;
//        std::vector< gcis_index<>::_node > nodes;
//        GCIS_index_toy.find_range(sfx.second,nodes);
//    }
//
//}
TEST( sdsltest, lcp_build){
    std::string a = "abracadabra";
    sdsl::lcp_bitcompressed<> LCP;
    sdsl::construct_im(LCP,a.c_str(),1);

}
TEST(gcisIndexToyExampleBS,locate){


    load_toy_example();
    srand(time(nullptr));
    for (uint i = 0; i < 100000 ; ++i) {

        uint l = rand()%T.size();
        uint r = rand()%T.size();

        if( l > r ) std::swap(l,r);
        uint m = r - l + 1 ;
        if( m < 5  ) r+=5;
        if(r >= T.size() || m == 1) continue;

        std::string s; s.resize(m);
        std::copy(T.begin()+l,T.begin()+l+m, s.begin());
        std::vector< gcis_index<>::len_type > occ;
//        s = "AG";
        GCIS_index_toy.locate(s,occ);
        std::sort(occ.begin(),occ.end());

        size_t pos = T.find(s.c_str(),0);
        std::vector<gcis_index<>::len_type > test_occ;
        while(pos!= std::string::npos){
            test_occ.push_back(pos);
            pos = T.find(s.c_str(),pos+1);
        }

        for (int i = 0 ; i < occ.size(); ++i) occ[i]-=4;

//        std::cout<<s<<std::endl;

        ASSERT_EQ(test_occ,occ);

    }



}

