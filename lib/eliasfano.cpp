#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"
#include "eliasfano.hpp"



eliasfano_codec::eliasfano_codec(eliasfano_codec &&rhs) :
m_size(rhs.m_size),
sdv(std::move(rhs.sdv)),
m_sel(std::move(rhs.m_sel))
{
    sdsl::util::init_support(m_sel,&sdv);
}


void eliasfano_codec::encode(sdsl::bit_vector&v){
    sdv = sdsl::sd_vector<>(v);
    sdsl::util::init_support(m_sel,&sdv);
    m_size = v.size();
}

void eliasfano_codec::encode(sdsl::int_vector<> &v) {
    for(uint64_t i=1;i<v.size();i++){
        v[i]+= v[i-1]+1;
    }
    sdv = sdsl::sd_vector<>(v.begin(),v.end());
    sdsl::util::init_support(m_sel,&sdv);
    m_size = v.size();
}

uint64_t eliasfano_codec::size_in_bytes() {
    uint64_t total_bytes = sdsl::size_in_bytes(sdv);
    total_bytes += sdsl::size_in_bytes(m_sel);
    return total_bytes;
}

uint64_t eliasfano_codec::operator[](uint64_t i) {
    return i==0 ? m_sel(1) : m_sel(i+1) - m_sel(i) -1;
}

uint64_t eliasfano_codec::get_next(){
   //TODO: implement sequencial scan in a efficient way
}

uint64_t eliasfano_codec::size(){
    return m_size;
}

void eliasfano_codec::load(std::istream &i) {
    i.read((char*) &m_size,sizeof(m_size));
    sdv.load(i);
    m_sel.load(i);
    m_sel.set_vector(&sdv);
}
void eliasfano_codec::serialize(std::ostream &o) {
    o.write((char*) &m_size,sizeof(m_size));
    sdv.serialize(o);
    m_sel.serialize(o);
}

uint64_t eliasfano_codec::pos(uint64_t i) {
    return m_sel(i+1);
}

uint64_t eliasfano_codec::access_bv(uint64_t i) {
    return sdv[i];
}
