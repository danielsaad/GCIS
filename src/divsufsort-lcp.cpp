#include "divsufsort-lcp.h"
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
        std::cerr << "Usage: ./divsufsort <input_file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    auto start = timer::now();
    char *str;
    load_string_from_file(str, argv[1]);
    size_t n = strlen(str);
    saidx_t *SA = new saidx_t[n];
    saidx_t *LCP = new saidx_t[n];

    std::cout << "Building SA + LCP with Divsufsort." << std::endl;
    divsuflcpsort((sauchar_t *)str, SA, LCP, n);
    auto stop = timer::now();

    string output_file_basename(argv[2]);
    std::ofstream output(output_file_basename + ".sa", std::ios::binary);
    std::ofstream output_lcp(output_file_basename + ".lcp", std::ios::binary);
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)SA, sizeof(saidx_t) * n);
    output_lcp.write((const char *)&n, sizeof(n));
    output_lcp.write((const char *)LCP, sizeof(saidx_t) * n);
    output.close();
    output_lcp.close();
    return 0;
}