#include "gcis_eliasfano_index.hpp"
#include "index_builder.hpp"
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

int main(int argc, char *argv[]) {
    std::ofstream output(argv[3], std::ios::binary);
    char *str;
    if (strcmp(argv[1], "-c") == 0) {
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        g.serialize(output);
        delete[] str;
    } else if (strcmp(argv[1], "-d") == 0) {
        ifstream grammar_file(argv[2], std::ifstream::in);
        gcis::elias_fano_grammar g;
        g.load(grammar_file);
        str = g.decode();
        output.write(str, strlen(str));
        delete[] str;
    } else if (strcmp(argv[1], "-i") == 0) {
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        gcis::index_basics<gcis::elias_fano_grammar> index(g, str);
    }
    output.close();
}
