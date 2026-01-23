#include "gcis_s8b_lyndon.hpp"
#include "sais.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>

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
    std::cout.setf(std::ios_base::unitbuf);

#ifdef MEM_MONITOR
    mm.event("GC-IS Init");
#endif

    if (argc != 4 and argc != 3) {
        std::cerr << "Usage: \n"
                  << argv[0] << " -c <file_to_be_encoded> <output> to compress using GCIS\n"
                  << argv[0] << " -d <file_to_be_decoded> <output> to decompress using GCIS\n"
                  << argv[0] << "-dd <file_to_be_decode> to decompress using GCIS and discard the result (dummy decompress)\n"
                  << argv[0] << " -l <file_to_be_decoded> to compute the LA during GCIS decompression\n"
                  << argv[0] << " -lse <file_to_be_decoded> to compute the LA during GCIS decompression (semi-external version, the SA is streamed from disk)\n"
                  << argv[0] << " -lite <file_to_be_decoded> to compute the LA during GCIS decompression (lite version, does not require the original string)\n"
                  << argv[0] << " -litese <file_to_be_decoded> to compute the LA during GCIS decompression (lite semi-external version, does not require the original string and the SA is streamed from disk)\n"
                  << argv[0] << " -sa <file_to_be_decoded> to compute the SA during GCIS decompression\n";
        exit(EXIT_FAILURE);
    }

    // Dictionary type
    char *mode = argv[1];
    if (strcmp(mode, "-c") == 0) {
        int_t n;
        char *str;
        load_string_from_file(str, argv[2], n);
        std::ofstream output(argv[3], std::ios::binary);

#ifdef MEM_MONITOR
        mm.event("GC-IS Compress");
#endif

        auto start = timer::now();
        auto g = gcis_lyndon();
        g.encode(str, n);
        auto stop = timer::now();

#ifdef MEM_MONITOR
        mm.event("GC-IS Save");
#endif

        cout << "input:\t" << n << " bytes" << endl;
        cout << "output:\t" << g.size_in_bytes() << " bytes" << endl;
        cout << "time: " << (double)duration_cast<seconds>(stop - start).count()
             << " seconds" << endl;

        g.serialize(output);
        output.close();
        delete[] str;
    } else if (strcmp(mode, "-d") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3], std::ios::binary);

#ifdef MEM_MONITOR
        mm.event("GC-IS Load");
#endif

        auto g = gcis_lyndon();
        g.load(input);

#ifdef MEM_MONITOR
        mm.event("GC-IS Decompress");
#endif

        auto start = timer::now();
        char *str;
        int_t n;
        tie(str, n) = g.decode();
        auto stop = timer::now();

        cout << "input:\t" << g.size_in_bytes() << " bytes" << endl;
        cout << "output:\t" << n << " bytes" << endl;
        cout << "time: "
             << (double)duration_cast<milliseconds>(stop - start).count() /
                    1000.0
             << setprecision(2) << fixed << " seconds" << endl;

        output.write(str, n);
        input.close();
        output.close();
    } else if (strcmp(mode, "-dd") == 0) {
        std::ifstream input(argv[2]);

#ifdef MEM_MONITOR
        mm.event("GC-IS Load");
#endif

        auto g = gcis_lyndon();
        g.load(input);

#ifdef MEM_MONITOR
        mm.event("GC-IS Decompress");
#endif

        auto start = timer::now();
        char *str;
        int_t n;
        tie(str, n) = g.decode();
        auto stop = timer::now();

        cout << "input:\t" << g.size_in_bytes() << " bytes" << endl;
        cout << "output:\t" << n << " bytes" << endl;
        cout << "time: "
             << (double)duration_cast<milliseconds>(stop - start).count() /
                    1000.0
             << setprecision(2) << fixed << " seconds" << endl;
        input.close();
        delete[] str;
    } else if (strcmp(mode, "-l") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3]);
        auto g = gcis_lyndon();
        g.load(input);
        int_t *LA = nullptr;
        auto [str, n] = g.decode_lyndon(&LA);
        output.write((const char *)LA, sizeof(int_t) * n);
        delete[] LA;
        delete[] str;
    } else if (strcmp(mode, "-lse") == 0) {
        auto glse = gcis_lyndon_semi_external();
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3]);
        glse.load(input);
        int_t *LA = nullptr;
        auto [str, n] = glse.decode_lyndon(&LA);
        output.write((const char *)LA, sizeof(int_t) * n);
        delete[] LA;
        delete[] str;
    } else if (strcmp(mode, "-lite") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3]);
        auto glite = gcis_lyndon_lite();
        glite.load(input);
        int_t *LA = nullptr;
        uint_t n = glite.decode_lyndon(&LA);
        output.write((const char *)LA, sizeof(int_t) * n);
        delete[] LA;
    } else if (strcmp(mode, "-litese") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3]);
        auto glitese = gcis_lyndon_lite_semi_external();
        glitese.load(input);
        int_t *LA = nullptr;
        uint_t n = glitese.decode_lyndon(&LA);
        output.write((const char *)LA, sizeof(int_t) * n);
        delete[] LA;
    } else if (strcmp(mode, "-sa") == 0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3]);
        auto g = gcis_lyndon();
        g.load(input);
        uint_t *SA = nullptr;
        auto [t, n] = g.decode_saca(&SA);
        output.write((const char *)SA, sizeof(uint_t) * n);
        delete[] SA;
        delete[] t;
    }

#ifdef MEM_MONITOR
    mm.event("GC-IS Finish");
#endif

    return 0;
}