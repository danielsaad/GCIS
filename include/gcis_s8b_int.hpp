//
// Created by danielsaad on 1/19/18.
//

#ifndef GC_IS_GCIS_S8B_INT_HPP
#define GC_IS_GCIS_S8B_INT_HPP

#include "gcis.hpp"
#include "gcis_s8b_codec.hpp"
#include "util.hpp"

class gcis_s8b_int : public gcis_abstract<gcis_s8b_codec> {
  public:
    void encode(uint_t *s, int_t n) {
        uint_t *SA = new uint_t[n];
        std::cout << "Computing alphabet size..." << std::endl;
        set<int> alphabet(s, s + n);
        int_t K = alphabet.size();
        std::cout << "Alphabet size = " << K << std::endl;
        int cs = sizeof(int);
        int level = 0;
        std::cout << "Compressing..." << std::endl;
        gc_is((int_t *)s, SA, n, K, cs, level);
        delete[] SA;
    }

    pair<uint_t *, int_t> decode_int() {
        sdsl::int_vector<> r_string = reduced_string;
        uint_t *str = nullptr;
        for (int64_t i = g.size() - 1; i >= 0; i--) {
            sdsl::int_vector<> next_r_string;
            gcis_s8b_pointers_codec_level gd = std::move(g[i].decompress());
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
            next_r_string.resize(g[i].string_size);
            uint_t l = 0;
            if (i == 0) {
                // convert the reduced string in the original text
                str = new uint_t[g[i].string_size];
                for (auto t : g[i].tail) {
                    str[l++] = (uint_t)t;
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], str, l);
                }
            } else {
                // convert the reduced string in the previous reduced string
                for (uint64_t j = 0; j < g[i].tail.size(); j++) {
                    next_r_string[l++] = g[i].tail[j];
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], next_r_string, l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return make_pair(str, g[0].string_size);
    }

    pair<char *, int_t> decode() override {
        sdsl::int_vector<> r_string = reduced_string;
        char *str = nullptr;
        for (int64_t i = g.size() - 1; i >= 0; i--) {
            sdsl::int_vector<> next_r_string;
            gcis_s8b_pointers_codec_level gd = std::move(g[i].decompress());
            next_r_string.width(sdsl::bits::hi(g[i].alphabet_size - 1) + 1);
            next_r_string.resize(g[i].string_size);
            uint_t l = 0;
            if (i == 0) {
                // convert the reduced string in the original text
                str = new char[g[i].string_size];
                for (auto t : g[i].tail) {
                    str[l++] = (char)t;
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], str, l);
                }
            } else {
                // convert the reduced string in the previous reduced string
                for (uint64_t j = 0; j < g[i].tail.size(); j++) {
                    next_r_string[l++] = g[i].tail[j];
                }
                for (uint64_t j = 0; j < r_string.size(); j++) {
                    gd.expand_rule(r_string[j], next_r_string, l);
                }
                r_string = std::move(next_r_string);
            }
        }
        return make_pair(str, g[0].string_size);
    }
};

#endif // GC_IS_GCIS_S8B_INT_HPP
