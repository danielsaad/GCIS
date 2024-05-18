#include "simple8b.hpp"
#include "sdsl/bits.hpp"

class simple8b_selector {
public:
    uint32_t n_items;
    uint32_t n_bits;
    uint32_t n_waste;
};

static simple8b_selector s8b_selector[16] = {
        {240, 0, 60},
        {120, 0, 60},
        {60, 1, 0},
        {30, 2, 0},
        {20,  3, 0},{15, 4, 0},{12, 5, 0},
        {10, 6, 0},
        {8, 7, 4},
        {7, 8, 4},
        {6,10,0},
        {5,12,0},
        {4,15,0},
        {3,20,0},
        {2,30,0},
        {1,60,0}
};


// Encode a entire array in simple8b
void simple8b_codec::encode(uint64_t* v,uint64_t n){
    //TODO: implement
    m_size+=n;
    return;
}

// Encode a single integer in simple8b
void simple8b_codec::encode(uint64_t k){
    m_buf[m_buf_i++] = k;
    if(m_buf_i == BUF_SIZE) {
        m_buf_i = flush(BUF_SIZE);
    }
    m_size++;
}

// Finish the encoding process
void simple8b_codec::encode(){
    pack();
}



uint64_t simple8b_codec::get_next(){
        if(m_buf_i == m_buf_size){
            m_buf_size = decode_to_buffer();
            m_buf_i = 0;
        }
        return m_buf[m_buf_i++];
}

uint64_t simple8b_codec::size_in_bytes() {
    uint64_t total_bytes = m_v.size() * sizeof(uint64_t);
    // total_bytes += (5+BUF_SIZE) * sizeof(uint64_t);
    return total_bytes;
}


void simple8b_codec::serialize(std::ostream& o){
    uint64_t size = m_v.size();
    o.write((char*) &m_size,sizeof(uint64_t));
    o.write((char*) &size,sizeof(uint64_t));
    o.write((char*) m_v.data(),sizeof(uint64_t)*size);
}



void simple8b_codec::load(std::istream& i){
    uint64_t size;
    i.read((char*) &m_size,sizeof(uint64_t));
    i.read((char*) &size,sizeof(uint64_t));
    m_v.resize(size);
    i.read((char*) m_v.data(),sizeof(uint64_t)*size);
}

uint64_t simple8b_codec::size(){
    return m_size;
}

void simple8b_codec::reset(){
    m_buf_i = 0;
    m_buf_size = 0;
    m_cur_word = 0;
}

// Flushes the content in buffer to the final vector
// Returns the number of itens remaining in the buffer
uint64_t simple8b_codec::flush(uint64_t buf_size){
    uint64_t index = pack(buf_size);
    /* We reached the last element of the buffer but we do not know if this is the correct selector. We need to look to the other elements */
    if(index < BUF_SIZE){
        // copy elements to the buffer start
        for(uint64_t k=0;k<BUF_SIZE -index;k++){
            m_buf[k] = m_buf[index + k];
        }
    }
    return BUF_SIZE-index;
}

uint64_t simple8b_codec::decode_to_buffer(){
    uint64_t index = 0;
    while(true && m_cur_word<m_v.size()){
        uint64_t word = m_v[m_cur_word];
        uint64_t n_items = s8b_selector[word & 0xf].n_items;
        uint64_t n_bits =  s8b_selector[word & 0xf].n_bits;
        word >>= 4;
        if(index+n_items>BUF_SIZE)
            break;
        for(uint64_t i=0;i<n_items;i++){
            m_buf[index++] = word & sdsl::bits::lo_set[n_bits];
            word >>= n_bits;
        }
        m_cur_word++;
    }
    return index;
}

// Pack a series of integers.
// Returns the number of packed integers
uint64_t simple8b_codec::pack(uint64_t buf_size){
    uint64_t index = 0;
    while(index<buf_size){
        for(uint64_t i=0;i<16;i++){
            uint64_t word = i;
            uint64_t shift = 4;
            uint64_t n_items = 0;
            for(uint64_t j=index;j<buf_size;j++){
                if(n_items == s8b_selector[i].n_items)
                    break;
                if(m_buf[j] > ((1ULL << s8b_selector[i].n_bits) - 1))
                    break;
                word |=  m_buf[j] << shift;
                shift += s8b_selector[i].n_bits;
                n_items++;
            }

            // We have a filled a word with selector 'i'
            if(n_items==s8b_selector[i].n_items){
                m_v.push_back(word);
                index += n_items;
                break;
            }
        }
    }
    return index;
}

uint64_t simple8b_codec::pack(){
    uint64_t index = 0;
    while(index<m_buf_i){
        for(uint64_t i=0;i<16;i++){
            uint64_t word = i;
            uint64_t shift = 4;
            uint64_t n_items = 0;
            for(uint64_t j=index;j<m_buf_i;j++){
                if(n_items == s8b_selector[i].n_items)
                    break;
                if(m_buf[j] > ((1ULL << s8b_selector[i].n_bits) - 1))
                    break;
                word |=  m_buf[j] << shift;
                shift += s8b_selector[i].n_bits;
                n_items++;
            }

            // We have a filled a word with selector 'i'
            if((n_items==s8b_selector[i].n_items) || (index + n_items == m_buf_i) ){
                m_v.push_back(word);
                index += n_items;
                break;
            }
        }
    }
    m_buf_i = 0;
    return index;
}
