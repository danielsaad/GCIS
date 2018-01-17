#include <cstdint>
#include <vector>
#include <limits>
#include <fstream>

#ifndef SIMPLE8B_HPP_DEFINED
#define SIMPLE8B_HPP_DEFINED


class simple8b_codec{
public:

    // Encode a entire array in simple8b
    void encode(uint64_t* v,uint64_t n);

    // Encode a single integer in simple8b
    void encode(uint64_t k);

    // Finish the encoding process
    void encode();

    uint64_t get_next();

    uint64_t size();

    void reset();

    uint64_t size_in_bytes();

    void serialize(std::ostream& o);

    void load(std::istream& i);

private:

    std::vector<uint64_t> m_v;

    static const uint64_t BUF_SIZE = 4800;
    // cache to store integers to be coded and decoded
    uint64_t m_buf[BUF_SIZE];
    uint64_t m_buf_i = 0;
    // Number of integers coded
    uint64_t m_size = 0;
    // Number of integers decoded
    uint64_t m_buf_size = 0;
    // Current word to be extracted from m_v
    uint64_t m_cur_word = 0;



    // Flushes the content in buffer to the final vector
    // Returns the number of itens remaining in the buffer
    uint64_t flush(uint64_t buf_size);

    uint64_t decode_to_buffer();

    // Pack a series of integers.
    // Returns the number of packed integers
    uint64_t pack(uint64_t buf_size);
    uint64_t pack();
};

#endif //SIMPLE8B_HPP_DEFINED
