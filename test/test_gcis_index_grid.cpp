//
// Created by ale on 10-01-20.
//

#include <ostream>
#include <fstream>
#include "gtest/gtest.h"
#include "gcis_index_grid.hpp"

#define GRID_MAX_ROW 100
#define GRID_MAX_COLS 100
#define GRID_MAX_LAB 100
#define GRID_MAX_2D_QUERY 10000


typedef std::pair<gcis_index_private::gcis_index_grid<>::point,gcis_index_private::gcis_index_grid<>::point> query;

void generate_random_points(std::vector<gcis_index_private::gcis_index_grid<>::lpoint>& );
void generate_random_2d_queries(std::vector<query>& );
bool evaluate_point(const gcis_index_private::gcis_index_grid<>::point&, const query&);

TEST(dumbSuite, dumb){
    ASSERT_TRUE(true);
}

TEST(gcisGrid, build){

    std::srand(std::time(nullptr));
    std::vector<gcis_index_private::gcis_index_grid<>::lpoint> v(GRID_MAX_COLS);


    for (uint i = 1; i <= GRID_MAX_COLS; ++i)
    {
        uint row = std::rand()%GRID_MAX_ROW+1;
        uint tag = std::rand()%GRID_MAX_LAB+1;
        gcis_index_private::gcis_index_grid<>::lpoint p = std::make_pair(std::make_pair(row,i),tag);
        v[i-1] = p;
    }

    gcis_index_private::gcis_index_grid<> G(v, GRID_MAX_ROW, GRID_MAX_COLS);

    ASSERT_TRUE(true);
}

TEST(gcisGrid, firstColumnLabel){

    std::vector<gcis_index_private::gcis_index_grid<>::lpoint> v;
    generate_random_points(v);
    gcis_index_private::gcis_index_grid<> G(v, GRID_MAX_ROW, GRID_MAX_COLS);

    for (uint i = 1; i <= GRID_MAX_COLS; ++i){
#ifdef DEBUG
        std::cout<<i<<"->"<<v[i-1].second<<std::endl;
#endif
        EXPECT_EQ(v[i-1].second , G.first_column_point(i));
    }

}

TEST(gcisGrid, labels_search_2d)
{

    std::vector<gcis_index_private::gcis_index_grid<>::lpoint> v;
    generate_random_points(v);
    gcis_index_private::gcis_index_grid<> G(v, GRID_MAX_ROW, GRID_MAX_COLS);

    std::vector<query> q;
    generate_random_2d_queries(q);
    std::vector<std::vector<gcis_index_private::gcis_index_grid<>::ltype>> l(GRID_MAX_2D_QUERY);

    for (int i = 0; i < GRID_MAX_2D_QUERY ; ++i)
    {
        for (int j = 0; j < GRID_MAX_COLS; ++j)
            if(evaluate_point(v[j].first,q[i]))  l[i].push_back(v[j].second);

        std::sort(l[i].begin(),l[i].end());
    }

    for (int i = 0; i < GRID_MAX_2D_QUERY ; ++i)
    {

#ifdef DEBUG
        std::cout<<"q["<<i+1<<"] = (["<<q[i].first.first<<" "<<q[i].first.second<<"],["<<q[i].second.first<<" "<<q[i].second.second<<"])"<<std::endl;
#endif
        std::vector<gcis_index_private::gcis_index_grid<>::ltype >  occ;

        G.labels_search_2d(q[i].first,q[i].second,occ);
        std::sort(occ.begin(),occ.end());

        ASSERT_EQ(occ,l[i]);
    }
}


void generate_random_points(std::vector<gcis_index_private::gcis_index_grid<>::lpoint>& v){

    std::srand(std::time(nullptr));

    v.resize(GRID_MAX_COLS);

#ifdef DEBUG
    std::cout<<"Generating random points...."<<std::endl;
#endif
    for (uint i = 1; i <= GRID_MAX_COLS; ++i)
    {
        uint row = std::rand()%GRID_MAX_ROW+1;
        uint tag = std::rand()%GRID_MAX_LAB+1;
#ifdef DEBUG
        std::cout<<"<"<<row<<","<<i<<">["<<tag<<"]"<<std::endl;
#endif
        gcis_index_private::gcis_index_grid<>::lpoint p = std::make_pair(std::make_pair(row,i),tag);
        v[i-1] = p;
    }


}
void generate_random_2d_queries(std::vector<query>& q)
{
    std::srand(std::time(nullptr));

#ifdef DEBUG
    std::cout<<"generate_random_2d_queries........."<<std::endl;
#endif
    q.resize(GRID_MAX_2D_QUERY);

    for (int i = 0; i < GRID_MAX_2D_QUERY; ++i)
    {
         gcis_index_private::gcis_index_grid<>::point lup = std::make_pair(std::rand()%GRID_MAX_ROW + 1,std::rand()%GRID_MAX_COLS + 1);
         gcis_index_private::gcis_index_grid<>::point rbt = std::make_pair(std::rand()%GRID_MAX_ROW + 1,std::rand()%GRID_MAX_COLS + 1);

         if(lup.first > rbt.first) std::swap(lup.first,rbt.first);
         if(lup.second > rbt.second) std::swap(lup.second,rbt.second);

         q[i] = query(lup,rbt);

#ifdef DEBUG
         std::cout<<"q["<<i+1<<"] = (["<<lup.first<<" "<<lup.second<<"],["<<rbt.first<<" "<<rbt.second<<"])"<<std::endl;
#endif

    }
}
bool evaluate_point(const gcis_index_private::gcis_index_grid<>::point& p, const query& bbox){
    return (p.first  >= bbox.first.first && p.first <= bbox.second.first &&
            p.second >= bbox.first.second && p.second <= bbox.second.second ) ;
};
