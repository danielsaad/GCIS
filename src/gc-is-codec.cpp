#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include "gcis.hpp"
#include "gcis_s8b.hpp"
#include "gcis_unary.hpp"
#include "gcis_eliasfano.hpp"
#include "gcis_eliasfano_no_lcp.hpp"
#include "../external/malloc_count/malloc_count.h"
#include "sais.h"

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

void load_string_from_file(char*& str,char* filename){
        std::ifstream f(filename,std::ios::binary);
        f.seekg (0, std::ios::end);
        uint64_t size = f.tellg();
        f.seekg(0,std::ios::beg);
        str = new char[size+1];
        f.read(str,size);
        str[size]=0;
        f.close();
};

int main(int argc, char* argv[]){


    #ifdef MEM_MONITOR
    mm.event("GC-IS Init");
    #endif


    if(argc!=4){
        std::cerr << "Usage: ./gc-is-codec -c <file_to_be_encoded> <output>\n"  <<
                  "./gc-is-codec -d <file_to_be_decoded> <output>\n" <<
                  "./gc-is-codec -s <file_to_be_decoded> <output>\n" <<
                  "./gc-is-codec -a <input_file> <output>\n" <<
                  "./gc-is-codec -e <encoded_file> <query file>\n";

        exit(EXIT_FAILURE);
    }

    // Dictionary type

    gcis_dictionary<gcis_eliasfano_codec> d;
    char* mode = argv[1];
    if(strcmp(mode,"-c")==0) {
        char* str;
        load_string_from_file(str, argv[2]);
        std::ofstream output(argv[3],std::ios::binary);

        #ifdef MEM_MONITOR
        mm.event("GC-IS Compress");
        #endif

				auto start = timer::now();
        d.encode(str);
				auto stop = timer::now();

        #ifdef MEM_MONITOR
        mm.event("GC-IS Save");
        #endif

				cout<<"input:\t"<<strlen(str)<<" bytes"<<endl;
				cout<<"output:\t"<<d.size_in_bytes()<<" bytes"<<endl;

        d.serialize(output);
        delete[] str;
        output.close();
				cout << "time: " << (double)duration_cast<seconds>(stop-start).count()<<" seconds" << endl;
			}
			else if(strcmp(mode,"-d")==0) {

        std::ifstream input(argv[2]);
        std::ofstream output(argv[3], std::ios::binary);

        #ifdef MEM_MONITOR
                mm.event("GC-IS Load");
        #endif

        d.load(input);

        #ifdef MEM_MONITOR
                mm.event("GC-IS Decompress");
        #endif

				auto start = timer::now();
        char *str = d.decode();
				auto stop = timer::now();

        output.write(str, strlen(str));
        input.close();
        output.close();

				cout<<"input:\t"<<d.size_in_bytes()<<" bytes"<<endl;
				cout<<"output:\t"<<strlen(str)<<" bytes"<<endl;
				cout << "time: " << (double)duration_cast<seconds>(stop-start).count()<<" seconds" << endl;
    }
    else if(strcmp(mode,"-s")==0) {
        std::ifstream input(argv[2]);
        std::ofstream output(argv[3], std::ios::binary);

        #ifdef MEM_MONITOR
        mm.event("GC-IS/SACA Load");
        #endif

        d.load(input);

        #ifdef MEM_MONITOR
        mm.event("GC-IS/SACA Decompress");
        #endif

				uint_t *SA;
				auto start = timer::now();
        char *str = d.decode_saca(&SA);
				auto stop = timer::now();

				size_t n = strlen(str)+1;

				#if CHECK
					if(!d.suffix_array_check(SA, (unsigned char*)str, (uint_t) n, sizeof(char), 0)) std::cout << "isNotSorted!!\n";
					else std::cout << "isSorted!!\n";
				#endif

        // output.write(str, n-1);
        d.suffix_array_write(SA, n, argv[3], "sa");
        output.write((const char*) &n,sizeof(n));
        //output.write((const char*)SA,sizeof(uint_t)*n);
        input.close();
        output.close();

				cout<<"input:\t"<<d.size_in_bytes()<<" bytes"<<endl;
				cout<<"output:\t"<<n-1<<" bytes"<<endl;
				cout<<"SA:\t"<<n*sizeof(uint_t)<<" bytes"<<endl;
							
				cout << "time: " << (double)duration_cast<seconds>(stop-start).count()<<" seconds" << endl;
				delete[] SA;
        delete[] str;
    }
    else if(strcmp(mode,"-A")==0) {

				char* str;
    		load_string_from_file(str, argv[2]);

    		size_t n = strlen(str)+1;
    		uint_t* SA = new uint_t[n];    

    		std::cout << "Building SA with SAIS." << std::endl;
				auto start = timer::now();
    		//sais_u8((sa_uint8_t*) str,SA,n,k);
    		d.saca(str, SA, n);
				auto stop = timer::now();

				#if CHECK
					if(!d.suffix_array_check(SA, (unsigned char*)str, (uint_t) n, sizeof(char), 0)) std::cout << "isNotSorted!!\n";
					else std::cout << "isSorted!!\n";
				#endif
    		//std::ofstream output(argv[3], std::ios::binary);
    		//output.write((const char*) &n,sizeof(n));
    		//output.write((const char*)SA,sizeof(sa_int32_t)*n);
        d.suffix_array_write((uint_t*)SA, n, argv[2], "SA");
    		//output.close();

				cout<<"input:\t"<<d.size_in_bytes()<<" bytes"<<endl;
				cout<<"SA:\t"<<n*sizeof(uint_t)<<" bytes"<<endl;
							
				cout << "time: " << (double)duration_cast<seconds>(stop-start).count()<<" seconds" << endl;
		}
    else if(strcmp(mode,"-e")==0){
        std::ifstream input (argv[2],std::ios::binary);
        std::ifstream query(argv[3]);

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
    mm.event("GC-IS Finish");
    #endif

    return 0;
}
