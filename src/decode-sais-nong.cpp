#include "gcis.hpp"
#include "gcis_eliasfano.hpp"
#include "sais.h"
#include "sais_nong.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

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
    int_t n;
    char *str;
    tie(str,n) = d.decode();
    auto stop = timer::now();

    cout << "Decompression Time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << setprecision(2) << fixed << " seconds" << endl;

    // Compute the Suffix Array
    start = timer::now();
    cout << "Computing the Suffix Array with SAIS-NONG" << endl;
    int32_t *SA = new int32_t[n];
    SA_IS((unsigned char *)str, SA, n, 255, sizeof(char), 0);
    stop = timer::now();
    cout << "Suffix Array Construction Time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << setprecision(2) << fixed << " seconds" << endl;

    string ouf_basename = argv[2];
    std::ofstream output(ouf_basename + ".sa", std::ios::binary);
    n--;
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)&SA[1], sizeof(int32_t) * n);
    return 0;
}
