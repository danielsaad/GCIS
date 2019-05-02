#include "sais_nong.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

void load_string_from_file(char *&str, char *filename) {
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

    if (argc != 3) {
        std::cerr << "Usage: ./sais-nong <input_file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    auto start = timer::now();
    char *str;
    load_string_from_file(str, argv[1]);

    size_t n = strlen(str) + 1;
    int32_t *SA = new int32_t[n];
    int32_t k = 256;

    std::cout << "Building SA with SAIS NONG." << std::endl;
    SA_IS((unsigned char *)str, SA, n, 255, sizeof(char), 0);
    auto stop = timer::now();

    cout << "time: " << (double)duration_cast<seconds>(stop - start).count()
         << " seconds" << endl;

    string ouf_basename = argv[2];
    std::ofstream output(ouf_basename + ".sa", std::ios::binary);
    n--;
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)&SA[1], sizeof(int32_t) * n);

    output.close();
    return 0;
}