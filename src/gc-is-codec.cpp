#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include "gcis.hpp"
#include "gcis_s8b.hpp"
#include "gcis_unary.hpp"
#include "gcis_eliasfano.hpp"


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


    #ifdef MEM_MONITOR
    mm.event("GC-IS Init");
    #endif


    if(argc!=4){
        std::cerr << "Usage: ./gc-is-codec -c <file_to_be_encoded> <output>\n"  <<
                  "./gc-is-codec -d <file_to_be_decoded> <output>\n" <<
                  "./gc-is-codec -e <encoded_file> <query file>\n";

        exit(EXIT_FAILURE);
    }
    char* mode = argv[1];
    if(strcmp(mode,"-c")==0) {
        char* str;
        load_string_from_file(str, argv[2]);
        std::ofstream output(argv[3],std::ios::binary);
        gcis_dictionary<gcis_eliasfano_codec> d;

        #ifdef MEM_MONITOR
        mm.event("GC-IS Compress");
        #endif

        d.encode(str);

        #ifdef MEM_MONITOR
        mm.event("GC-IS Save");
        #endif

        d.serialize(output);
        delete[] str;
    }
    else if(strcmp(mode,"-d")==0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3], std::ios::binary);
        gcis_dictionary<gcis_eliasfano_codec> d;

        #ifdef MEM_MONITOR
                mm.event("GC-IS Load");
        #endif

        d.load(input);

        #ifdef MEM_MONITOR
                mm.event("GC-IS Decompress");
        #endif

        char *str = d.decode();
        output.write(str, strlen(str));
    }
    else if(strcmp(mode,"-e")==0){
        std::ifstream input (argv[2],std::ios::binary);
        std::ifstream query(argv[3]);
        gcis_dictionary<gcis_eliasfano_codec> d;

        #ifdef MEM_MONITOR
                mm.event("GC-IS Load");
        #endif

        d.load(input);

        #ifdef MEM_MONITOR
                mm.event("GC-IS Extract");
        #endif
        uint64_t l,r;
        while(query >> l >> r){
            sdsl::int_vector<> str = d.extract(l,r);
            cout << l << " " << r << " = ";
            cout.flush();
            for(uint64_t i=0;i<r-l+1;i++){
                cout << (char) str[i];
            }
            cout << endl;
            cout.flush();
        }
    }
    else{
        std::cerr << "Invalid mode: use -c for compression, -d for decompression or -e for extraction.\n";
        exit(EXIT_FAILURE);
    }

    #ifdef MEM_MONITOR
    mm.event("GC-iS Finish");
    #endif

    return 0;
}
