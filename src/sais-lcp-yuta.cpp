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

    size_t n = strlen(str)+1;
    int *SA = new int[n];
    int* LCP = new int[n];

    std::cout << "Building SA+LCP with SAIS-LCP-YUTA." << std::endl;
    auto start = timer::now();
    sais_lcp((unsigned char*) str,SA,LCP,n);
    auto stop = timer::now();

    cout << "input:\t" << strlen(str) << " bytes" << endl;

    cout << "time: "
         << (double)duration_cast<milliseconds>(stop - start).count() / 1000.0
         << " seconds" << endl;


    string ouf_basename(argv[2]);
    std::ofstream output(ouf_basename+".sa", std::ios::binary);

    n--;
    output.write((const char*) &n,sizeof(n));
    output.write((const char*) &SA[1],sizeof(int)*n);

    std::ofstream output_lcp(ouf_basename+".lcp", std::ios::binary);
    output_lcp.write((const char*) &n,sizeof(n));
    output_lcp.write((const char*)&LCP[1],sizeof(int)*n);


    output.close();
    output_lcp.close();

    return 0;
}
