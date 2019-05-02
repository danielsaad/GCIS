#include <cstring>
#include <fstream>
#include <iostream>
#include "gcis.hpp"
#include "gcis_eliasfano.hpp"
#include "divsufsort.h"

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: ./decode-divsufsort <compressed-file> <output_file>"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    auto start = timer::now();
    gcis_dictionary<gcis_eliasfano_codec> d;
    std::ifstream compressed_file(argv[1], std::ios::binary);
    cout << "Loading " << argv[2] << endl;
    d.load(compressed_file);
    cout << "Decompressing " << argv[2] << endl;
    char *str = d.decode();
    auto stop = timer::now();
    cout << "Decompression Time: " << (double)duration_cast<milliseconds>(stop - start).count()/1000.0
         << setprecision(2) << fixed << " seconds" << endl;

    start = timer::now();
    cout << "Computing the Suffix Array with Divsufsort" << endl;
    size_t n = strlen(str);
    saidx_t* SA = new saidx_t[n];
    divsufsort((sauchar_t*)str,SA,n);

    cout << "input:\t" << strlen(str) << " bytes" << endl;

    stop = timer::now();

    cout << "Suffix Array Construction Time: " << (double)duration_cast<milliseconds>(stop - start).count()/1000.0
         << " seconds" << endl;

    string output_file_basename(argv[2]);
    std::ofstream output(output_file_basename + ".sa", std::ios::binary);
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)SA, sizeof(saidx_t) * n);
    output.close();

    return 0;
}
