#include "gcis.hpp"
#include "gcis_eliasfano.hpp"
#include "sais.h"
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

    auto start = timer::now();
    gcis_dictionary<gcis_eliasfano_codec> d;
    std::ifstream compressed_file(argv[1], std::ios::binary);
    cout << "Loading " << argv[2] << endl;
    d.load(compressed_file);
    cout << "Decompressing " << argv[2] << endl;
    int_t n;
    char *str;
    tie(str,n) = d.decode();
    auto stop = timer::now();
    cout << "Decompression Time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << setprecision(2) << fixed << " seconds" << endl;

    start = timer::now();
    cout << "Computing the Suffix Array with SAIS-YUTA" << endl;
    sa_int32_t *SA = new sa_int32_t[n];
    sa_int32_t k = 256;
    sais_u8((sa_uint8_t *)str, SA, n, k);
    cout << "input:\t" << strlen(str) << " bytes" << endl;
    stop = timer::now();

    cout << "Suffix Array Construction Time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << " seconds" << endl;

    string ouf_basename(argv[2]);
    std::ofstream output(ouf_basename + ".sa", std::ios::binary);
    n--;
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)&SA[1], sizeof(sa_int32_t) * n);
    output.close();
    return 0;
}
