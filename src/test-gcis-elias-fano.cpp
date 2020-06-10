#include "gcis_eliasfano_index.hpp"
#include "gcis_index_bs.hpp"
#include "index_builder.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <fstream>

// #define TEST
// #define PROFILER

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#ifdef MEM_MONITOR
#include "mem_monitor.hpp"
#include "util.hpp"
#endif

const int PATTERN_LEN = 10;

using namespace sdsl;
using std::pair;
using std::vector;
using std::chrono::steady_clock;
using mapwtcfg_t =
    sdsl::wt_gmr<sdsl::int_vector<>,
                 sdsl::inv_multi_perm_support<8UL, sdsl::int_vector<>>>;

bool test_display(const gcis_index_private::gcis_index_bs<> &G,
                  const unsigned char *T, std::ifstream &queries_file) {
    int number_of_queries;
    size_t substring_sz;

    queries_file >> number_of_queries >> substring_sz;
    cout << "Number of substrings: " << number_of_queries << endl;
    cout << "Substring length: " << substring_sz << endl;
    vector<pair<uint64_t, uint64_t>> queries(number_of_queries);
    for (int i = 0; i < number_of_queries; i++) {
        queries_file >> queries[i].first >> queries[i].second;
    }
    uint64_t total_time_ns = 0;
    for (size_t i = 0; i < queries.size(); i++) {
        std::string str;
        str.resize(substring_sz);
        auto t1 = std::chrono::steady_clock::now();
        G.display_L(queries[i].first, substring_sz, str);
        auto t2 = std::chrono::steady_clock::now();
        total_time_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)
                .count();

#ifdef TEST
        std::string s;
        s.resize(substring_sz);
        std::copy(T + queries[i].first, T + queries[i].first + substring_sz,
                  s.begin());
        if (str != s) {
            std::cout << "s:" << s << "\n display:" << str << std::endl;
            return false;
        }
#endif
    }

    cout << "Extract took: " << total_time_ns / 1e6 << " ms." << endl;
    cout << "Extract took: " << total_time_ns / (1e3 * number_of_queries)
         << " us per substring" << endl;
    cout << "Extract took: "
         << total_time_ns / (1e3 * number_of_queries * substring_sz)
         << " us per symbol" << endl;
    return true;
}

csa_wt<wt_huff<rrr_vector<>>, 32, 32> build_fm_index(const char *file) {
    csa_wt<wt_huff<rrr_vector<>>, 32, 32> fmidx;
    construct(fmidx, file, 1);
    return fmidx;
}

bool test_locate(const gcis_index_private::gcis_index_bs<> &G,
                 std::ifstream &queries_file, const char *text_file) {

#ifdef TEST
    auto fmidx = build_fm_index(text_file);
#endif

    vector<string> queries;
    int number_of_queries;
    int pattern_len;
    queries_file.read((char *)&number_of_queries, sizeof(number_of_queries));
    queries_file.read((char *)&pattern_len, sizeof(pattern_len));
    cout << "Number of patterns: " << number_of_queries << endl;
    cout << "Pattern length: " << pattern_len << endl;
    while (queries_file) {
        string pattern;
        pattern.resize(pattern_len);
        queries_file.read((char *)pattern.data(), sizeof(char) * pattern_len);
        queries.emplace_back(pattern);
    }

    uint64_t total_time_ns = 0;
    uint64_t number_of_occ = 0;
    for (auto &pattern : queries) {
        vector<gcis_index_private::gcis_index<>::len_type> occ;
        auto t1 = std::chrono::steady_clock::now();
        G.locate(pattern, occ);
        auto t2 = std::chrono::steady_clock::now();
        total_time_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)
                .count();
        std::sort(occ.begin(), occ.end());
        number_of_occ += occ.size();
#ifdef TEST
        auto test_occ_ctnr =
            sdsl::locate(fmidx, pattern.begin(), pattern.end());
        std::sort(test_occ_ctnr.begin(), test_occ_ctnr.end());
        vector<gcis_index_private::gcis_index<>::len_type> test_occ;
        for (uint64_t i = 0; i < test_occ_ctnr.size(); i++) {
            test_occ.push_back(test_occ_ctnr[i]);
        }
        if (test_occ != occ) {
            return false;
        }
#endif
    }
    cout << "Total of occ: " << number_of_occ << endl;
    cout << "Locate took: " << total_time_ns / 1e6 << " ms.\n";
    cout << "Locate took: " << total_time_ns / (1e3 * number_of_queries)
         << " us per pattern." << endl;
    cout << "Locate took: " << total_time_ns / (1e3 * number_of_occ)
         << " us/occ." << endl;
    return true;
}

void print_usage(char *argv[]) {
    cout << "Invalid option." << endl;
    cout << "Usage:" << endl;
    cout << "\t" << argv[0]
         << " -i <input_text> <output_index> (Builds the index)" << endl;
    cout << "\t" << argv[0]
         << " -x <input_text> <input_index> <queries file> (Perform  "
            "extraction of substrings)"
         << endl;
    cout << "\t" << argv[0]
         << " -p <input_text> <input_index> <queries file> (Perform  pattern "
            "matching)"
         << endl;
    cout << "\t" << argv[0] << " -info <input_index> (Prints index information)"
         << endl;
}

void load_string_from_file(char *&str, const char *filename,uint_t& n) {
    std::ifstream f(filename, std::ios::binary);
    f.seekg(0, std::ios::end);
    uint_t size = f.tellg();
    f.seekg(0, std::ios::beg);
    str = new char[size];
    f.read((char*) str, size);
    n = size;
    f.close();
};

int main(int argc, char *argv[]) {
#ifdef PROFILER
    ProfilerStart("gperf.out");
#endif
    char *str;
    uint_t n;
    if (argc < 2) {
        print_usage(argv);
        return 0;
    }
    if (strcmp(argv[1], "-c") == 0) {
        std::ofstream output(argv[3], std::ios::binary);
        load_string_from_file(str, argv[2],n);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str,n);
        g.serialize(output);
        output.close();
        delete[] str;
    } else if (strcmp(argv[1], "-d") == 0) {
        std::ofstream output(argv[3], std::ios::binary);
        ifstream grammar_file(argv[2], std::ifstream::in);
        gcis::elias_fano_grammar g;
        g.load(grammar_file);
        tie(str,n) = g.decode();
        output.write(str, n);
        output.close();
        delete[] str;
    } else if (strcmp(argv[1], "-i") == 0) {
        std::cout << "INDEX CASE\n";
        //        sleep(5);
        load_string_from_file(str, argv[2],n);

#ifdef MEM_MONITOR
        mm.event("GCIS Grammar Build");
#endif

        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str,n);

#ifdef MEM_MONITOR
        mm.event("GCIS Index Basics");
#endif
        gcis::index_basics<gcis::elias_fano_grammar, mapwtcfg_t> index(g, str);
        std::cout << "gcis::index_basics<gcis::elias_fano_grammar>\n";
        //        sleep(5);
        gcis_index_private::gcis_index_bs<sdsl::sd_vector<>, sdsl::sd_vector<>,
                                          mapwtcfg_t>
            gcisIndexBs;
        std::cout << "default gcis_index_private::gcis_index_bs<>\n";
        //        sleep(5);
        gcisIndexBs.set_bvfocc(index.m_focc);
        gcisIndexBs.set_vt(index.m_str);
        gcisIndexBs.set_bvt(index.m_t);
        gcisIndexBs.set_pi(index.m_pi);
        gcisIndexBs.set_nt(index.m_wt);
        gcisIndexBs.set_tree(index.m_bv_dfuds);
        gcisIndexBs.set_l(index.m_l);
        std::cout << "loading index structures\n";
        gcisIndexBs.print();
#ifdef MEM_MONITOR
        mm.event("GCIS Index Points");
#endif

        std::vector<gcis_index_private::gcis_index_grid<>::lpoint> points(
            index.m_grid_points.size());
        for (uint i = 0; i < points.size(); ++i)
            points[i] = {{index.m_grid_points[i].prev_rule + 1, i + 1},
                         index.m_grid_points[i].id};
#ifdef MEM_MONITOR
        mm.event("GCIS Index Grid");
#endif

        gcis_index_private::gcis_index_grid<> _grid(points, index.m_pi.size(),
                                                    index.m_grid_points.size());
        std::cout << "gcis_index_private::gcis_index_grid<> _grid" << std::endl;
        gcisIndexBs.set_grid(_grid);

        std::cout << "size in bytes:" << gcisIndexBs.size_in_bytes()
                  << std::endl;
        std::ofstream output(argv[3], std::ios::binary);
#ifdef MEM_MONITOR
        mm.event("GCIS Index Serialize");
#endif

        gcisIndexBs.serialize(output);
        output.close();
    } else if (strcmp(argv[1], "-x") == 0) {
        // argv[2] = text
        // argv[3] = index
        // argv[4] = queries file
        cout << "Extract mode" << endl;
        ifstream index_file(argv[3], std::ifstream::in);
        ifstream queries_file(argv[4], std::ifstream::in);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;
        gcisIndexBs.load(index_file);
        load_string_from_file(str, argv[2],n);
        if (!test_display(gcisIndexBs, (const unsigned char *)str,
                          queries_file)) {
            std::cout << "TEST DISPLAY DOES NOT PASS\n";
            return 0;
        }
        std::cout << "TEST DISPLAY PASSED\n";
    } else if (strcmp(argv[1], "-p") == 0) {
        // argv[2] = tex
        // argv[3] = index
        // argv[4] = queries file
        ifstream index_file(argv[3], std::ifstream::in);
        ifstream queries_file(argv[4], std::ifstream::in);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;
        gcisIndexBs.load(index_file);
        cout << "Locate mode" << endl;
        if (!test_locate(gcisIndexBs, queries_file, argv[2])) {
            std::cout << "TEST LOCATE DOES NOT PASS\n";
            return 0;
        }
        std::cout << "TEST LOCATE PASSED\n";
    } else if (strcmp(argv[1], "-info") == 0) {
        // Displays index info
        cout << "Printing index info" << endl;
        ifstream index_file(argv[2], std::ifstream::in);
        gcis_index_private::gcis_index_bs<> index;
        index.load(index_file);
        index.print_size_in_bytes();
    } else {       
        print_usage(argv);
    }
#ifdef PROFILER
    ProfilerStop();
#endif
}