#include "gcis_s8b_int.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

void load_int_string_from_file(uint_t *&str, char *filename, int_t &n) {
    std::ifstream f(filename, std::ios::binary);
    f.seekg(0, std::ios::end);
    n = f.tellg()/sizeof(uint_t);
    f.seekg(0, std::ios::beg);
    str = new uint_t[n];
    f.read((char *)str, n * sizeof(uint_t));
    f.close();
};

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: \n"
                  << argv[0] << " -c <file_to_be_encoded> <output>\n"
                  << argv[0] << " -d <file_to_be_decoded> <output>\n";
        exit(EXIT_FAILURE);
    }

    // Dictionary type
    char *mode = argv[1];
    gcis_s8b_int d;
    if (strcmp(mode, "-c") == 0) {
        int_t n;
        uint_t *str;
        load_int_string_from_file(str, argv[2], n);
        std::ofstream output(argv[3], std::ios::binary);
        auto start = timer::now();
        d.encode(str, n);
        auto stop = timer::now();

        cout << "input:\t" << n << " bytes" << endl;
        cout << "output:\t" << d.size_in_bytes() << " bytes" << endl;
        cout << "time: " << (double)duration_cast<seconds>(stop - start).count()
             << " seconds" << endl;

        d.serialize(output);
        output.close();
        delete[] str;
    } else if (strcmp(mode, "-d") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3], std::ios::binary);

        d.load(input);

        auto start = timer::now();
        uint_t *str;
        int_t n;
        tie(str, n) = d.decode_int();
        auto stop = timer::now();

        cout << "input:\t" << d.size_in_bytes() << " bytes" << endl;
        cout << "output:\t" << n * sizeof(uint_t) << " bytes" << endl;
        cout << "time: "
             << (double)duration_cast<milliseconds>(stop - start).count() /
                    1000.0
             << setprecision(2) << fixed << " seconds" << endl;

        output.write((char *)str, n * sizeof(uint_t));
        input.close();
        output.close();
    } else {
        std::cerr << "Invalid mode. Use -c for encoding and -d for decoding.\n";
        exit(EXIT_FAILURE);
    }

    return 0;
}