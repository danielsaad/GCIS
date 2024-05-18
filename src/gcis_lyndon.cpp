#include "../external/malloc_count/malloc_count.h"
#include "gcis_s8b.hpp"
#include "sais.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <print>

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

#ifdef MEM_MONITOR
    mm.event("GC-IS Init");
#endif

    if (argc != 4 and argc != 3) {
        std::cerr << "Usage: \n"
                  << argv[0] << " -c <file_to_be_encoded> <output>\n"
                  << argv[0] << " -d <file_to_be_decoded> <output>\n"
                  << argv[0] << " -l <file_to_be_decoded> \n"
                  << argv[0] << " -sa <file_to_be_decoded>\n";
        exit(EXIT_FAILURE);
    }

    // Dictionary type
    gcis_lyndon g;
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
    } else if (strcmp(mode, "-l") == 0) {
        std::ifstream input(argv[2]);
        g.load(input);
        uint_t *LA = nullptr;
        auto [str, n] = g.decode_lyndon(&LA);
        for (int i = 0; i < n; i++) {
            std::println("LA[{}] = {}", i, LA[i]);
        }
    } else if (strcmp(mode, "-sa") == 0) {
        std::ifstream input(argv[2]);
        g.load(input);
        uint_t *SA = nullptr;
        auto [t, n] = g.decode_saca(&SA);
        for (int i = 0; i < n; i++) {
            std::println("SA[{}] = {}", i, SA[i]);
        }
    }

#ifdef MEM_MONITOR
    mm.event("GC-IS Finish");
#endif

    return 0;
}
