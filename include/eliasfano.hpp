//
// Created by danielsaad on 23/10/17.
//

#ifndef GC_IS_ELIASFANO_HPP_HPP
#define GC_IS_ELIASFANO_HPP_HPP

#include "sdsl/sd_vector.hpp"
#include "sdsl/int_vector.hpp"


class eliasfano_codec{
public:
    eliasfano_codec() = default;
    eliasfano_codec(eliasfano_codec& rhs) = default;
    eliasfano_codec(eliasfano_codec&& rhs);
    void encode(sdsl::int_vector<>& v);
    void encode(sdsl::bit_vector& v);
    uint64_t operator[](uint64_t i);
    uint64_t pos(uint64_t i);
    uint64_t access_bv(uint64_t i);
    uint64_t get_next();

    uint64_t size();

    uint64_t size_in_bytes();

    void serialize(std::ostream& o);

    void load(std::istream& i);

private:
    uint64_t m_size = 0;
    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::select_1_type m_sel;
};


#endif //GC_IS_ELIASFANO_HPP_HPP
