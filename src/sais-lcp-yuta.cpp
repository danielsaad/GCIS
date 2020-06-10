#include "sais-lcp.h"
#include "util.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

void load_string_from_file(char *&str, char *filename, int_t &n) {
    std::ifstream f(filename, std::ios::binary);
    f.seekg(0, std::ios::end);
    n = f.tellg();
    f.seekg(0, std::ios::beg);
    str = new char[n];
    f.read(str, n);
    f.close();
};

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: ./sais-yuta <input_file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    char *str = nullptr;
    int_t n = 0;
    load_string_from_file(str, argv[1], n);

    int *SA = new int_t[n];
    int *LCP = new int_t[n];

    std::cout << "Building SA+LCP with SAIS-LCP-YUTA." << std::endl;
    auto start = timer::now();
    sais_lcp((unsigned char *)str, SA, LCP, n);
    auto stop = timer::now();

    cout << "input:\t" << strlen(str) << " bytes" << endl;

    cout << "time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << " seconds" << endl;

    string ouf_basename(argv[2]);
    std::ofstream output(ouf_basename + ".sa", std::ios::binary);

    output.write((const char *)&n, sizeof(n));
    output.write((const char *)SA, sizeof(n) * n);

    std::ofstream output_lcp(ouf_basename + ".lcp", std::ios::binary);
    output_lcp.write((const char *)&n, sizeof(n));
    output_lcp.write((const char *)LCP, sizeof(int) * n);

    output.close();
    output_lcp.close();

    return 0;
}
