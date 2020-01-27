#include "gcis_eliasfano_index.hpp"
#include "gcis_index_bs.hpp"
#include "index_builder.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <fstream>

const int PATTERN_LEN = 10;

using namespace sdsl;

bool test_display(const gcis_index_private::gcis_index_bs<> &G,
                  std::string &T) {
    // srand(time(nullptr));
    size_t N = T.size();

    for (int i = 0; i < 1000000; ++i) {
        uint p = 9;  // rand()%N;
        uint m = 16; /// rand()%(N - p);

        std::string str, s;
        s.resize(m);
        std::copy(T.begin() + p, T.begin() + p + m, s.begin());
        str.reserve(m);
        G.display(p, m, str);
        if (str != s) {
            std::cout << "s:" << s << "\n display:" << str << std::endl;
            //         G.display(p,m,str);
            return false;
        }
    }

    return true;
}

csa_wt<wt_huff<rrr_vector<>>, 32, 32> build_fm_index(const char *file) {
    csa_wt<wt_huff<rrr_vector<>>, 32, 32> fmidx;
    construct(fmidx, file, 1);
        return fmidx;
}

bool test_locate(const gcis_index_private::gcis_index_bs<> &G, std::string &T,
                 const char *text_file) {

    auto fmidx = build_fm_index(text_file);

    srand(time(nullptr));
    size_t N = T.size();
    for (uint i = 0; i < 1000; ++i) {

        // uint l = rand() % N;
        // uint r = std::min<int>(N - 1, l + PATTERN_LEN);

        uint_t r = std::max<int>(rand() % N, PATTERN_LEN - 1);
        uint_t l = r - PATTERN_LEN + 1;
        // uint l = 494988;
        // uint r = 494997;

        cout << "[l,r] = "
             << "[" << l << "," << r << "]" << endl;
        std::string s;
        s.resize(r - l + 1);
        std::copy(T.begin() + l, T.begin() + r + 1, s.begin());
        cout << "Pattern = " << s << endl;
        std::vector<gcis_index_private::gcis_index<>::len_type> occ;
        cout << "GCIS locating" << endl;
        G.locate(s, occ);
        std::sort(occ.begin(), occ.end());
        std::cout << "GCIS finished" << endl;
        cout << "FM-Index locating" << endl;
        auto test_occ_ctnr = sdsl::locate(fmidx, s.begin(), s.end());
        std::cout << "FM-Index finished" << endl;
        std::sort(test_occ_ctnr.begin(), test_occ_ctnr.end());
        vector<gcis_index_private::gcis_index<>::len_type> test_occ;
        for (uint64_t i = 0; i < test_occ_ctnr.size(); i++) {
            test_occ.push_back(test_occ_ctnr[i]);
        }

        cout << "Number of occurences: " << occ.size() << endl;
        cout << "Number of occurences (FM): " << test_occ.size() << endl;

        if (test_occ != occ) {
            return false;
        }
    }

    return true;
}

void load_string_from_file(char *&str, const char *filename) {
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
    char *str;
    if (strcmp(argv[1], "-c") == 0) {
        std::ofstream output(argv[3], std::ios::binary);
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        g.serialize(output);
        output.close();
        delete[] str;
    } else if (strcmp(argv[1], "-d") == 0) {
        std::ofstream output(argv[3], std::ios::binary);
        ifstream grammar_file(argv[2], std::ifstream::in);
        gcis::elias_fano_grammar g;
        g.load(grammar_file);
        str = g.decode();
        output.write(str, strlen(str));
        output.close();
        delete[] str;
    } else if (strcmp(argv[1], "-i") == 0) {
        std::cout << "INDEX CASE\n";
        //        sleep(5);
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        gcis::index_basics<gcis::elias_fano_grammar> index(g, str);
        std::cout << "gcis::index_basics<gcis::elias_fano_grammar>\n";
        //        sleep(5);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;
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

        std::vector<gcis_index_private::gcis_index_grid<>::lpoint> points(
            index.m_grid_points.size());
        for (uint i = 0; i < points.size(); ++i)
            points[i] = {{index.m_grid_points[i].prev_rule + 1, i + 1},
                         index.m_grid_points[i].id};
        gcis_index_private::gcis_index_grid<> _grid(points, index.m_pi.size(),
                                                    index.m_grid_points.size());
        std::cout << "gcis_index_private::gcis_index_grid<> _grid" << std::endl;
        gcisIndexBs.set_grid(_grid);

        std::cout << "size in bytes:" << gcisIndexBs.size_in_bytes()
                  << std::endl;
        std::ofstream output(argv[3], std::ios::binary);
        gcisIndexBs.serialize(output);
        output.close();
    } else if (strcmp(argv[1], "-p") == 0) {
        ifstream index_file(argv[3], std::ifstream::in);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;
        gcisIndexBs.load(index_file);
        load_string_from_file(str, argv[2]);
        std::string ss = str;
        //        if(!test_display(gcisIndexBs,ss)){
        //            std::cout<<"TEST DISPLAY DOES NOT PASS\n";
        //            return 0;
        //        }
        //        std::cout<<"TEST DISPLAY PASSED\n";

        if (!test_locate(gcisIndexBs, ss, argv[2])) {
            std::cout << "TEST LOCATE DOES NOT PASS\n";
            return 0;
        }
        std::cout << "TEST LOCATE PASSED\n";
    }
}