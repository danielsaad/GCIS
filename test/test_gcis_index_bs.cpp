//
// Created by ale on 14-01-20.
//

#include "gcis_index.hpp"
#include "gtest/gtest.h"
#include <fstream>
#include <ostream>

using namespace gcis_index_private;

struct dumbGrammar;
std::string T = "AGCTTTTCATTCTGACTGCAACAGCTTTTCATTCTGACTGCAAC";

std::map<char, uint32_t> sigma = {
    {'A', 0},
    {'C', 1},
    {'G', 10},
    {'T', 13},
};
std::map<uint32_t, std::string> non_t_exp = {
    /*$A****/ {0, "A"},
    /*$C****/ {1, "C"},
    /*X1****/ {2, "AAC"},
    /*Y2****/ {3, "ATTCTGACTGCAAC"},
    /*S*****/ {4, "AGCTTTTCATTCTGACTGCAACAGCTTTTCATTCTGACTGCAAC"},
    /*X2****/ {5, "ACTGC"},
    /*Y3****/ {6, "ATTCTGACTGC"},
    /*X6****/ {7, "CTTTTC"},
    /*tail2*/ {8, "AGCTTTTC"},
    /*Y1****/ {9, "AACAGCTTTTC"},
    /*$G****/ {10, "G"},
    /*X3****/ {11, "AG"},
    /*X5****/ {12, "CTG"},
    /*$T****/ {13, "T"},
    /*X4****/ {14, "ATT"}

};

std::map<uint32_t, uint> rule_pre = {
    /*$A****/ {0, 2},
    /*$C****/ {1, 3},
    /*S*****/ {2, 33},
    /*X1****/ {3, 39},
    /*Y2****/ {4, 1},
    /*X2****/ {5, 26},
    /*Y3****/ {6, 17},
    /*X6****/ {7, 10},
    /*tail2*/ {8, 6},
    /*Y1****/ {9, 32},
    /*$G****/ {10, 4},
    /*X3****/ {11, 7},
    /*X5****/ {12, 22},
    /*$T****/ {13, 5},
    /*X4****/ {14, 18}};

sdsl::bit_vector dfuds = {
    1, 1, 0,                   // DFUDS INIT
    1, 1, 1, 1, 1, 1, 1, 1, 0, // 4
    0, 0, 0, 0, 1, 1, 0,       // AGCT 8
    1, 1, 0,                   // 11
    0, 0,                      // AG
    1, 1, 1, 1, 1, 1, 0,       // 7
    0, 0, 0, 0, 0, 0,          // CTTTTC
    1, 1, 1, 0,                // 6
    1, 1, 1, 0,                // 14
    0, 0, 0,                   // ATT
    1, 1, 1, 0,                // 12
    0, 0, 0,                   // CTG
    1, 1, 1, 1, 1, 0,          // 5
    0, 0, 0, 0, 0,             // ACTGC
    1, 1, 1, 0,                // 9
    1, 1, 1, 0,                // 2
    0, 0, 0,                   // AAC
    0, 0,                      // 11 7
    1, 1, 1, 1, 0,             // 3
    0, 0, 0, 0                 // 14 12 5 2
};
std::map<uint, long> off = {
    {3, -4},                                                   // 4
    {12, -4}, {13, -3}, {14, -2}, {15, -1}, {16, 0},           // AGCT 8
    {19, 0},                                                   // 11
    {22, 0},  {23, 1},                                         // AG
    {24, 2},                                                   // 7
    {31, 2},  {32, 3},  {33, 4},  {34, 5},  {35, 6},  {36, 7}, // CTTTTC
    {37, 8},                                                   // 6
    {41, 8},                                                   // 14
    {45, 8},  {46, 9},  {47, 10},                              // ATT
    {48, 11},                                                  // 12
    {52, 11}, {53, 12}, {54, 13},                              // CTG
    {55, 14},                                                  // 5
    {61, 14}, {62, 15}, {63, 16}, {64, 17}, {65, 18},          // ACTGC
    {66, 19},                                                  // 9
    {70, 19},                                                  // 2
    {74, 19}, {75, 20}, {76, 21},                              // AAC
    {77, 22}, {78, 24},                                        // 11 7
    {79, 30},                                                  // 3
    {84, 30}, {85, 33}, {86, 36}, {87, 41}};

///////////////////////////S//////////////t2///X3////////X6///////////////////Y3///X4///////////X5///////////X2////////////////Y1//X1////////////X3//X6//Y2//X4//X5//X2//X1
sdsl::int_vector<> X = {4, 0, 1,  10, 13, 8,  11, 0, 10, 7,  1,  13, 13, 13, 13,
                        1, 6, 14, 0,  13, 13, 12, 1, 13, 10, 5,  0,  1,  13, 10,
                        1, 9, 2,  0,  0,  1,  11, 7, 3,  14, 12, 5,  2};
sdsl::bit_vector focc = {1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0,
                         0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                         0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
sdsl::int_vector<> pi = {4, 0, 1, 10, 13, 8, 11, 7, 6, 14, 12, 5, 9, 2, 3};
sdsl::bit_vector t = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
std::string st = "AGCTTTTCATTCTGACTGCAAC";
sdsl::int_vector<> wt = {11, 7, 14, 12, 5, 2};
//////////////////////////////////////////////1///////////////////2///////////////////3///////////////////4//////////////
//////////////////////////0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7
//////////////////////////A,C,G,T,A,G,C,T,T,T,T,C,A,T,T,C,T,G,A,C,T,G,C,A,A,C,A,G,C,T,T,T,T,C,A,T,T,C,T,G,A,C,T,G,C,A,A,C
sdsl::bit_vector l = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0};

gcis_index<> GCIS_index_toy;

void load_toy_example() {

    sdsl::inv_perm_support<> inv_pi(&pi);
    sdsl::wt_gmr<> gmrwt;
    sdsl::construct_im(gmrwt, wt);

    GCIS_index_toy.set_bvfocc(sdsl::sd_vector<>(focc));
    GCIS_index_toy.set_tree(dfuds);
    GCIS_index_toy.set_bvt(sdsl::sd_vector<>(t));
    GCIS_index_toy.set_vt(st);
    GCIS_index_toy.set_nt(gmrwt);
    GCIS_index_toy.set_pi(pi);
    GCIS_index_toy.set_l(l);
    GCIS_index_toy.set_sigma(4);
}

TEST(dumbSuite, dumb) { ASSERT_TRUE(true); }

TEST(gcisIndexToyExample, map_preorder) {
    load_toy_example();
    for (int i = 0; i < X.size(); ++i) {
        //        std::cout<<i<<std::endl;
        auto c = GCIS_index_toy.map_node(i + 1);
        ASSERT_EQ(c > 14 ? sigma[c] : c, X[i]);
    }
}
TEST(gcisIndexToyExample, map_rule) {
    load_toy_example();
    for (uint32_t i = 0; i < 14; ++i) {
        //        std::cout<<i<<std::endl;
        ASSERT_EQ(GCIS_index_toy.map_rule(i), rule_pre[i]);
    }
}
// TEST(gcisIndexToyExample,offset){
//    load_toy_example();
//    for (const auto &item : off) {
//        std::cout<<item.first<<std::endl;
//        ASSERT_EQ(GCIS_index_toy.offset_node(item.first),item.second + 4);
//    }
//}
//
TEST(gcisIndexToyExample, display) {
    load_toy_example();
    srand(time(nullptr));

    for (int i = 0; i < 1000000; ++i) {
        uint p = rand() % T.size();
        uint m = rand() % (T.size() - p);

        std::string str, s;
        s.resize(m);
        std::copy(T.begin() + p, T.begin() + p + m, s.begin());
        str.reserve(m);
        //        std::cout<<p<<","<<m<<" s:"<<s<<std::endl;
        GCIS_index_toy.display(p, m, str);
        ASSERT_EQ(str, s);
    }
}

TEST(gcisIndexToyExample, cmp_prefix_rule) {
    load_toy_example();
    srand(time(nullptr));

    for (int i = 0; i < 1000000; ++i) {
        uint p = rand() % T.size();
        uint m = rand() % (T.size() - p);
        uint rule = rand() % 15;
        while (rule == 0 || rule == 1 || rule == 13 || rule == 10 || rule == 4)
            rule = rand() % 15;

        std::string str, s;
        s.resize(m);
        std::copy(T.begin() + p, T.begin() + p + m, s.begin());
        //        std::cout<<p<<","<<m<<" s:"<<s<<std::endl;
        //        std::cout<<"rule:"<<rule<<std::endl;
        long long j = 0;

        int t = GCIS_index_toy.cmp_prefix_rule(rule, s, j);

        if (t == 0) {
            size_t f = s.find_first_of(non_t_exp[rule].c_str());
            ASSERT_EQ(f, 0);
        } else {
            int t2 = std::strcmp(s.c_str(), non_t_exp[rule].c_str());

            t2 = (t2 > 0) ? 1 : ((t2 < 0) ? -1 : 0);
            ASSERT_EQ(t2, t);
        }
    }
}
TEST(gcisIndexToyExample, cmp_suffix_rule) {
    load_toy_example();
    srand(time(nullptr));

    for (int i = 0; i < 1000000; ++i) {
        uint p = rand() % T.size();
        uint m = rand() % (T.size() - p);
        uint rule = rand() % 15;

        //        std::cout<<p<<","<<m<<","<<rule<<std::endl;
        if (!m)
            continue;
        if (rule == 0 || rule == 1 || rule == 13 || rule == 10 || rule == 4)
            continue;

        std::string str, s;
        s.resize(m);
        std::copy(T.begin() + p, T.begin() + p + m, s.begin());
        //        std::cout<<p<<","<<m<<" s:"<<s<<std::endl;
        //        std::cout<<"rule:"<<rule<<std::endl;
        uint j = s.size() - 1;

        int t = GCIS_index_toy.cmp_suffix_rule(rule, s, j);

        std::string ss = s;
        std::reverse(ss.begin(), ss.end());
        std::string sr = non_t_exp[rule];
        std::reverse(sr.begin(), sr.end());

        if (t == 0) {

            size_t f = ss.find_first_of(sr.c_str());
            ASSERT_EQ(f, 0);

        } else {
            int t2 = std::strcmp(ss.c_str(), sr.c_str());

            t2 = (t2 > 0) ? 1 : ((t2 < 0) ? -1 : 0);
            ASSERT_EQ(t2, t);
        }
    }
}
