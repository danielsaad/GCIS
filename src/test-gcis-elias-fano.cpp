#include "gcis_eliasfano_index.hpp"
#include <fstream>
void load_string_from_file(char *&str, const char *filename) {
    std::ifstream f(filename, std::ios::binary);
    f.seekg(0, std::ios::end);
    uint64_t size = f.tellg();
    f.seekg(0, std::ios::beg);
    str = new char[size + 1];
    f.read(str, size);
    str[size] = 0;
    f.close();
};

int main() {
    char* str = "AGCTTTTCATTCTGACTGCAACAGCTTTTCATTCTGACTGCAAC";
    gcis::grammar_builder<gcis::elias_fano_grammar> builder;
    auto g = builder.build(str);
    std::cout << "Grammar has been built!\n";
}
