#include <iostream>
#include <fstream>
#include <cstring>
#include "sais.h"


void load_string_from_file(char*& str,char* filename){
        std::ifstream f(filename,std::ios::binary);
        f.seekg (0, std::ios::end);
        uint64_t size = f.tellg();
        f.seekg(0,std::ios::beg);
        str = new char[size+1];
        f.read(str,size);
        str[size]='\0';
        f.close();
};

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The output array of suffixes.
 * @param n The length of the given string.
 * @param k The alphabet size.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
// SAIS_API
// sa_int32_t
// sais_u8(const sa_uint8_t *T, sa_int32_t *SA, sa_int32_t n, sa_int32_t k);

int main(int argc, char* argv[]){

    if(argc!=3){
        std::cerr << "Usage: ./sais -c <input_file>" << std::endl;
        exit(EXIT_FAILURE);
    }
    char* str;
    load_string_from_file(str, argv[2]);

    sa_int32_t n = strlen(str);
    sa_int32_t* SA = new sa_int32_t[n];    
    sa_int32_t k = 256;

    std::cout << "Building SA with SAIS." << std::endl;
    sais_u8((sa_uint8_t*) str,SA,n,k);
    return 0;
}
