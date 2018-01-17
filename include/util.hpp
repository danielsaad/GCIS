//
// Created by danielsaad on 4/4/17.
//

#ifndef GC_IS_UTIL_HPP
#define GC_IS_UTIL_HPP


#include <iostream>
#include <fstream>
#include <sdsl/int_vector.hpp>
#include <cstdint>

using namespace std;

#ifdef m64
    typedef int64_t  int_t;
    typedef uint64_t uint_t;
#else
    typedef int32_t  int_t;
    typedef uint32_t uint_t;
#endif


#ifdef MYDEBUG
extern ofstream dbg_file;
#endif

#ifdef REPORT
extern ofstream report_file;
#endif

#ifdef MEM_MONITOR
#include "mem_monitor.hpp"
extern mem_monitor mm;
#endif

#ifdef REPORT

template <typename T>
void print_report(T value){
    report_file << value;
}


template<typename T, typename... TArgs>
void print_report(T value, TArgs... args) {
    report_file << value << " ";
    print_report(args...);
}

#endif  //REPORT


#ifdef MYDEBUG

template<typename T, typename... TArgs>
void pdbg2(T value, TArgs... args) {
    dbg_file  << value << " ";
    pdbg2(args...);
}

template<typename T>
void pdbg2(T value){
    dbg_file << value << endl;
}


template<typename T>
void pdbg(T value) {
    dbg_file << value << endl;
}

template<typename T, typename... TArgs>
void pdbg(T value, TArgs... args) {
    dbg_file << value << " ";
    pdbg(args...);
}


template<uint8_t x>
ostream& operator<<(ostream& os,const sdsl::int_vector<x>& v){
    for(uint64_t i=0;i<v.size();i++){
        dbg_file << v[i] << " ";
    }
    dbg_file << "\nSize = " << v.size() << " Width = " << (int) v.width();
    return dbg_file;
}


template<typename T>
ostream& operator<<(ostream& os,const std::vector<T>& v){
    for(uint64_t i=0;i<v.size();i++){
        dbg_file << v[i] << " ";
    }
    return dbg_file;
}
#endif //MY_DEBUGs


#endif //GC_IS_UTIL_HPP
