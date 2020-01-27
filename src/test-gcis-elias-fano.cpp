#include "gcis_eliasfano_index.hpp"
#include "index_builder.hpp"
#include "gcis_index_bs.hpp"
#include <sdsl/suffix_arrays.hpp>

#include <fstream>


bool test_display(const gcis_index_private::gcis_index_bs<>& G, std::string& T)
{
    srand(time(nullptr));
    size_t N = T.size();

    for (int i = 0; i < 1000000; ++i)
    {
        uint p = 9;//rand()%N;
        uint m = 16;///rand()%(N - p);

        std::string str,s;
        s.resize(m);
        std::copy(T.begin()+p,T.begin()+p+m,s.begin());
        str.reserve(m);
        G.display(p,m,str);
        if(str != s){
            std::cout<<"s:"<<s<<"\n display:"<<str<<std::endl;
   //         G.display(p,m,str);
            return false;
        }

    }

    return true;
}


bool test_locate(const gcis_index_private::gcis_index_bs<>& G, std::string& T){





    srand(time(nullptr));
    size_t N = T.size();
    for (uint i = 0; i < 100000 ; ++i) {

        uint l = rand()%N;
        uint r = rand()%N;

        if( l > r ) std::swap(l,r);
        uint m = r - l + 1 ;
        if( m < 5  ) r+=5;
        if(r >= T.size() || m == 1) continue;

        std::string s; s.resize(m);
        std::copy(T.begin()+l,T.begin()+l+m, s.begin());
        std::vector< gcis_index_private::gcis_index<>::len_type > occ;
        G.locate(s,occ);
        std::sort(occ.begin(),occ.end());

        size_t pos = T.find(s.c_str(),0);
        std::vector<gcis_index_private::gcis_index<>::len_type > test_occ;

        while(pos!= std::string::npos)
        {
            test_occ.push_back(pos);
            pos = T.find(s.c_str(),pos+1);
        }


//        std::cout<<s<<std::endl;

        if(test_occ != occ) {
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
        std::cout<<"INDEX CASE\n";
//        sleep(5);
        load_string_from_file(str, argv[2]);
        gcis::grammar_builder<gcis::elias_fano_grammar> builder;
        auto g = builder.build(str);
        gcis::index_basics<gcis::elias_fano_grammar> index(g, str);
        std::cout<<"gcis::index_basics<gcis::elias_fano_grammar>\n";
//        sleep(5);
        gcis_index_private::gcis_index_bs<> gcisIndexBs;
        std::cout<<"default gcis_index_private::gcis_index_bs<>\n";
//        sleep(5);
        gcisIndexBs.set_bvfocc(index.m_focc);
        gcisIndexBs.set_vt(index.m_str);
        gcisIndexBs.set_bvt(index.m_t);
        gcisIndexBs.set_pi(index.m_pi);
        gcisIndexBs.set_nt(index.m_wt);
        gcisIndexBs.set_tree(index.m_bv_dfuds);
        gcisIndexBs.set_l(index.m_l);
        std::cout<<"loading index structures\n";
        gcisIndexBs.print();

        std::vector<gcis_index_private::gcis_index_grid<>::lpoint> points(index.m_grid_points.size());
        for (uint i = 0; i < points.size(); ++i)
            points[i] = {{index.m_grid_points[i].prev_rule+1,i+1},index.m_grid_points[i].id};
        gcis_index_private::gcis_index_grid<> _grid(points, index.m_pi.size(), index.m_grid_points.size());
        std::cout<<"gcis_index_private::gcis_index_grid<> _grid"<<std::endl;
        gcisIndexBs.set_grid(_grid);

        std::cout<<"size in bytes:"<<gcisIndexBs.size_in_bytes()<<std::endl;
        gcisIndexBs.serialize(output);
        std::string ss = str;
//        if(!test_display(gcisIndexBs,ss)){
//            std::cout<<"TEST DISPLAY DOES NOT PASS\n";
//            return 0;
//        }
//        std::cout<<"TEST DISPLAY PASSED\n";

        if(!test_locate(gcisIndexBs,ss)){
            std::cout<<"TEST LOCATE DOES NOT PASS\n";
            return 0;
        }
        std::cout<<"TEST LOCATE PASSED\n";


    }

    output.close();
}