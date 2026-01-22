//
// Created by danielsaad on 4/11/17.
//
#include "util.hpp"
#include <filesystem>
#include <fstream>

namespace gcis {
namespace util {
std::string make_temp_filename() {
    std::string base = std::filesystem::temp_directory_path();
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_int_distribution<int> dist(0, 0xFFFFFF);
    return base + "/tmp_" + std::to_string(dist(rng));
}

} // namespace util
} // namespace gcis
#ifdef MYDEBUG
ofstream dbg_file("dbg.log");
#endif

#ifdef REPORT
ofstream report_file("report.log");
#endif

#ifdef MEM_MONITOR
mem_monitor mm("mem-mon-out.csv");
#endif