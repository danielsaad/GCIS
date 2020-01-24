#include "gcis_eliasfano_index.hpp"
#include "index_builder.hpp"
#include "gcis_index_bs.hpp"

#include <fstream>
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
    std::ofstream output(argv[3], std::ios::binary);
    char *str;
    if (strcmp(argv[1], "-c") == 0) {
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        g.serialize(output);
        delete[] str;
    } else if (strcmp(argv[1], "-d") == 0) {
        ifstream grammar_file(argv[2], std::ifstream::in);
        gcis::elias_fano_grammar g;
        g.load(grammar_file);
        str = g.decode();
        output.write(str, strlen(str));
        delete[] str;
    } else if (strcmp(argv[1], "-i") == 0) {
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        gcis::index_basics<gcis::elias_fano_grammar> index(g, str);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;

        gcisIndexBs.set_bvfocc(index.m_focc);
        gcisIndexBs.set_vt(index.m_str);
        gcisIndexBs.set_bvt(index.m_t);
        gcisIndexBs.set_pi(index.m_pi);
        gcisIndexBs.set_nt(index.m_wt);
        gcisIndexBs.set_tree(index.m_bv_dfuds);
        gcisIndexBs.set_l(index.m_l);

        std::vector<gcis_index_private::gcis_index_grid<>::lpoint> points(index.m_grid_points.size());
        for (uint i = 0; i < points.size(); ++i)
            points[i] = {{index.m_grid_points[i].prev_rule+1,i+1},index.m_grid_points[i].id};
        gcis_index_private::gcis_index_grid<> _grid(points, index.m_pi.size(), index.m_grid_points.size());
        gcisIndexBs.set_grid(_grid);

        std::cout<<"size in bytes:"<<gcisIndexBs.size_in_bytes()<<std::endl;
    }

    output.close();
}
