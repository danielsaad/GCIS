#include <cstring>
#include <fstream>
#include <iostream>
#include "gcis.hpp"
#include "gcis_eliasfano.hpp"
#include "sais.h"
#include "sais_nong.hpp"

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: ./decode-sais-yuta <compressed-file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Decode compressed file
    auto start = timer::now();
    gcis_dictionary<gcis_eliasfano_codec> d;
    std::ifstream compressed_file(argv[1], std::ios::binary);
    std::ofstream output_file(argv[2], std::ios::binary);
    cout << "Loading " << argv[1] << endl;
    d.load(compressed_file);
    cout << "Decompressing " << argv[1] << endl;
    char *str = d.decode();
    auto stop = timer::now();

    cout << "Decompression Time: " << (double)duration_cast<milliseconds>(stop - start).count()/1000.0
         << setprecision(2) << fixed << " seconds" << endl;

    // Compute the Suffix Array
    start = timer::now();
    size_t n = strlen(str) + 1;
    cout << "Computing the Suffix Array with SAIS-NONG" << endl;
    int32_t *SA = new int32_t[n];
    SA_IS((unsigned char *)str, SA, n, 255, sizeof(char), 0);
    stop = timer::now();

    // We just want to compute the suffix array, we do not need to store it.
    // output.write((const char*)SA,sizeof(sa_int32_t)*n);
    output_file.close();

    cout << "Suffix Array Construction Time: " << (double)duration_cast<milliseconds>(stop - start).count()/1000.0
         << setprecision(2) << fixed << " seconds" << endl;
    return 0;
}
