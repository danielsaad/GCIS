//
// Created by ale on 13-01-20.
//

#ifndef GCIS_GCIS_INDEX_BS_HPP
#define GCIS_GCIS_INDEX_BS_HPP

#include <string>
#include <vector>
#include "gcis_index.hpp"
#include <sdsl/rmq_succinct_sada.hpp>

namespace gcis_index_private {

    template<typename, typename, typename, typename, typename>
    class gcis_index_bs;

    /**
    * IS-based compressed index binary search strategy
    *
    **/
    template<
            typename t_mapfbv = sdsl::sd_vector<>,
            typename t_maptbv = sdsl::sd_vector<>,
            typename t_mapwt  = sdsl::wt_gmr<>,
            typename t_gridbv = sdsl::rrr_vector<>,
            typename t_gridwt = sdsl::wt_int<>
    >
    class gcis_index_bs : public gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt> {

        typedef gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt> abs_gcis;
        typedef typename abs_gcis::len_type len_type;
        typedef typename gcis_index_grid<>::point point;
        typedef typename abs_gcis::_node _node;
//
    public:
        gcis_index_bs() = default;

        virtual ~gcis_index_bs() = default;

        /** Return the positions of the text where the pattern occurs */
        void locate(const std::string &s, std::vector<len_type> &) const;

    public:
        /**Implement an strategy to find the range of valid points of the grid*/
        void find_primary_occ(const std::string &s, std::vector<_node> &bfs_nodes) const;


        bool find_partition_range(const std::string &s, const uint & i, point& X, point&Y, bool ) const ;

        /**binary search lower bound*/
        template<typename C>
        bool bsearch_lowerBound(uint64_t &l, uint64_t r, const C &) const;

        /**binary search upper bound*/
        template<typename C>
        bool bsearch_upperBound(uint64_t l, uint64_t &r, const C &) const;
    };

    /**
      *
      * @tparam _gcis
      * @tparam C function to compare prefix / suffix rules
      * @param l left limit
      * @param r right limit
      * @param f_cmp function object to compare prefix / suffix rules
      * @param patt pattern string to search
      * @param n number of characters to compare
      * @return
      */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    template<typename C>
    bool gcis_index_bs<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::bsearch_lowerBound(uint64_t &l, uint64_t r,
                                                                                            const C &f_cmp
    ) const {
        bool found = false;

        while (l < r) {

            uint64_t mid = (l + r) / 2;
            int c = f_cmp(mid);

            if (c < 0) /** the rule is greater or equal than the pattern */
            {
                r = mid - 1;
            } else {
                if (c > 0)
                    l = mid + 1;
                else {
                    r = mid;
                    found = true;
                }
            }
        }


        if (!found && l == r && f_cmp(l) == 0) return true;

        return found;
    }

    /**
     *
     * @tparam _gcis
     * @tparam C function to compare prefix / suffix rules
     * @param l left limit
     * @param r right limit
     * @param f_cmp function object to compare prefix / suffix rules
     * @param patt pattern string to search
     * @param n number of characters to compare
     * @return
     */

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    template<typename C>
    bool gcis_index_bs<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::bsearch_upperBound(uint64_t l, uint64_t &r,
                                                                                            const C &f_cmp) const {
        bool found = false;

        while (l < r) {

            uint64_t mid = ceil((l + r) / 2.0);
            int c = f_cmp(mid);

            if (c < 0) /** the rule is greater or equal than the pattern */
            {
                r = mid - 1;
            } else {
                if (c > 0)
                    l = mid + 1;
                else {
                    l = mid;
                    found = true;
                }
            }

        }

        if (!found && l == r && f_cmp(r) == 0) return true;

        return found;
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool gcis_index_bs<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::find_partition_range(const std::string &s,
                                                                                              const uint &i,
                                                                                              gcis_index_bs::point &X,
                                                                                              gcis_index_bs::point &Y,
                                                                                              bool ch) const {

            uint64_t r_1 = 256, r_2 = this->n_rules-1;

            if(ch){
                 r_1 = 1; r_2 = 255;
            }

            /**
             * Binary search over rules
             * */
            bool stop = false;

            if (!bsearch_lowerBound(r_1, r_2, [this, &s, &i](const uint32_t &mid) {
                len_type j = i-1;
                int r = this->cmp_suffix_rule(mid , s, j);
//                int r = 1;
                if (r == 0 && j < 0) return 0; //string drain
                if (r == 0) return 1;// string not drain
                return r;//string differents
            })) return false;

            if (!stop && !bsearch_upperBound(r_1, r_2, [this, &s,&i](const uint32_t &mid) {
                len_type j = i-1;
                int r = this->cmp_suffix_rule(mid , s, j);
                if (r == 0 && j < 0) return 0; //string drain
                if (r == 0) return 1;// string not drain
                return r;//string differents
            })) return false;


            uint64_t c_1 = 1, c_2 = this->n_suffixes;
            /**
             * Binary search over suffix rules
             *
             * */
            if (!stop && !bsearch_lowerBound(c_1, c_2, [this, &s,&i](const uint32_t &mid) {
                auto sfx = this->_grid.first_column_point(mid)+1;
                len_type j = i;
                int r = this->cmp_suffix_grammar(sfx, s, j);
                if (r == 0 && j < s.size()) return 1;
                if (j == s.size()) return 0;
                return r;

            })) return false;

            if (!stop && !bsearch_upperBound(c_1, c_2, [this, &s,&i](const uint32_t &mid) {
                auto sfx = this->_grid.first_column_point(mid)+1;
                len_type j = i;
                int r = this->cmp_suffix_grammar(sfx, s, j);
                if (r == 0 && j < s.size()) return 1;
                if (j == s.size()) return 0;
                return r;

            })) return false;

            X = std::make_pair(r_1+1, c_1);
            Y = std::make_pair(r_2+1, c_2);

            return true;
    }



    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index_bs<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::find_primary_occ(const std::string &s,
                                                                                    std::vector<_node> &nodes) const
    {
        size_t n_s = s.size();

        for (long i = 1; i < n_s; ++i)
        {
                point X,Y;
                if(find_partition_range(s,i,X,Y,(i == 1))){
                    /**
                     * search points in the grid
                     * */
                    std::vector<uint32_t> labels;
                    this->_grid.labels_search_2d(X, Y, labels);
                    for (const auto &item : labels)
                        nodes.push_back(std::make_pair(item + 1, -1*i));
                }
        }

        point X,Y;

        if(find_partition_range(s,1,X,Y,false)){
            /**
             * search points in the grid
             * */
            std::vector<uint32_t> labels;
            this->_grid.labels_search_2d(X, Y, labels);
            for (const auto &item : labels)
                nodes.push_back(std::make_pair(item + 1, -1));
        }



    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index_bs<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::locate(const std::string &s,
                                                                                std::vector<len_type> &occ) const
    {
        /**
        * Find occ that cross boundaries grammar partition
        * */
        std::vector<_node> primary_occ;
        find_primary_occ(s, primary_occ);
        /**
         * Find all secondary occurences
         * */
        for (const auto &v : primary_occ)
            this->find_secondary_occ(v, occ);
    }


};


#endif //GCIS_GCIS_INDEX_BS_HPP
