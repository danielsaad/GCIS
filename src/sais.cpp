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


int main(int argc, char* argv[]){

    if(argc!=4) {
        std::cerr << "Usage: ./sais -c <input_file> <output_file>" << std::endl;
        exit(EXIT_FAILURE);
    }
    char* str;
    load_string_from_file(str, argv[2]);

    size_t n = strlen(str)+1;
    sa_int32_t* SA = new sa_int32_t[n];    
    sa_int32_t k = 256;

    std::cout << "Building SA with SAIS." << std::endl;
    sais_u8((sa_uint8_t*) str,SA,n,k);
    std::ofstream output(argv[3], std::ios::binary);
    output.write((const char*) &n,sizeof(n));
    output.write((const char*)SA,sizeof(sa_int32_t)*n);
    output.close();
    return 0;
}
