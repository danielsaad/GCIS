#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include "divsufsort.h"

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


    std::cout << "Building SA with Divsufsort." << std::endl;
    divsufsort((sauchar_t*) str, SA, n);
    auto stop = timer::now();

    std::ofstream output(argv[2], std::ios::binary);

    // We don't store the data for performance reasons
    // output.write((const char*) &n,sizeof(n));
    // output.write((const char*)SA,sizeof(int32_t)*n);
    cout << "time: " << (double)duration_cast<milliseconds>(stop - start).count()/1000.0
         << " seconds" << endl;
    output.close();
    return 0;
}