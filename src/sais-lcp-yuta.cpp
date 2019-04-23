#include "sais-lcp.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <chrono>

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
        std::cerr << "Usage: ./sais-yuta <input_file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    char *str;
    load_string_from_file(str, argv[1]);

    size_t n = strlen(str) + 1;
    int *SA = new int[n];
    int* LCP = new int[n-1];

    std::cout << "Building SA+LCP with SAIS-LCP-YUTA." << std::endl;
    auto start = timer::now();
    sais_lcp((unsigned char*) str,SA,LCP,n);
    auto stop = timer::now();

    cout << "input:\t" << strlen(str) << " bytes" << endl;

    std::ofstream output(argv[2], std::ios::binary);
    // output.write((const char*) &n,sizeof(n));
    // output.write((const char*)SA,sizeof(sa_int32_t)*n);
    output.close();

    cout << "time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << " seconds" << endl;
    return 0;
}
