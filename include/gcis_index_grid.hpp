//
// Created by ale on 10-01-20.
//

#ifndef GCIS_GCIS_INDEX_GRID_HPP
#define GCIS_GCIS_INDEX_GRID_HPP

#include <bits/forward_list.h>
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/bit_vectors.hpp"

//#define DEBUG 1

namespace gcis_index_private{


    template <typename ,typename,typename > class gcis_index_grid;

    /**
    * IS-based compressed index
    **/
    template <
            typename wavelet_tree  = sdsl::wt_int<>,
            typename XB_bit_vector = sdsl::rrr_vector<>,
            typename int_sequence = sdsl::int_vector<>

    > class gcis_index_grid {

        public:

            typedef uint32_t dtype; // type for dimensions
            typedef uint32_t ltype; // type for labels
            typedef std::pair<std::pair<dtype , dtype >, ltype > lpoint;
            typedef std::pair<dtype , dtype > point;

            /** Default constructor */
            gcis_index_grid();
            /** Copy constructor */
            gcis_index_grid( const gcis_index_grid& );
            /** Assignament operator */
            gcis_index_grid& operator = ( const gcis_index_grid& );
            /** Destructor */
            virtual ~gcis_index_grid();
            /**
             * Constructor receive a vector of points with labels
             * build the grid and return the vector sorted by rows and then by columns
             *
             * */
            gcis_index_grid( std::vector<lpoint> , const dtype& ,const dtype&);
            /** Serialize the structure to a file */
            virtual void serialize(std::ofstream & ) const;
            /** Load the structure from a file */
            virtual void load( std::ifstream & );
            /** Return the total size of the structure in bytes */
            virtual size_t size_in_bytes() const;

            /**
             * Return the label associated to the first point in the column.
             * */
            ltype first_column_point(const dtype& ) const;
            /**
             * return a list with all the labels in the 2d search range
             * */
            void  labels_search_2d(const point&, const point&, std::vector<ltype>& ) const;

            uint32_t n_cols() const{ return sb.size();}

    private:
            /** maps from virtual row to the position in XB bitvector */
            dtype  map(const dtype &) const;

//            dtype  unmap(const dtype  &) const;

            /**Build rank and select support structures for the grid */
            void build_rank_select_structures();
            /** Build the structures from a vector of labeled points */
            void build( std::vector<lpoint>&, const dtype& ,const dtype& );


            /**
             *  Attributes
             * */

            wavelet_tree sb;

            int_sequence sl;

            XB_bit_vector xb;
            typename XB_bit_vector::rank_1_type xb_rank1;
            typename XB_bit_vector::select_1_type xb_sel1;
            typename XB_bit_vector::select_0_type xb_sel0;
    };





    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::gcis_index_grid()
    {

    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::gcis_index_grid(const gcis_index_grid & B)
    {
        sb = B.sb;
        xb = B.xb;
        sl = B.sl;
        build_rank_select_structures();

    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>& gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::operator=(const gcis_index_grid & B)
    {
        sb = B.sb;
        xb = B.xb;
        sl = B.sl;
        build_rank_select_structures();
        return *this;
    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::~gcis_index_grid()
    {

    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::gcis_index_grid(std::vector<lpoint> points, const dtype& nrows,const dtype& ncols)
    {
        build(points, nrows, ncols);
    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    void gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::serialize(std::ofstream & fout) const
    {
        sdsl::serialize(sb,fout);

        sdsl::serialize(xb,fout);
        sdsl::serialize(xb_rank1,fout);
        sdsl::serialize(xb_sel0,fout);
        sdsl::serialize(xb_sel1,fout);

        sdsl::serialize(sl,fout);
    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    void gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::load(std::ifstream &fin)
    {
        sdsl::load(sb,fin);
        sdsl::load(xb,fin);

        sdsl::load(xb_rank1,fin);
        sdsl::load(xb_sel0 ,fin);
        sdsl::load(xb_sel1 ,fin);

        sdsl::load(sl,fin);

        xb_sel0.set_vector(&xb);
        xb_sel1.set_vector(&xb);
        xb_rank1.set_vector(&xb);

    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    size_t gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::size_in_bytes() const {
        return sdsl::size_in_bytes(sb) +
               sdsl::size_in_bytes(xb) +
               sdsl::size_in_bytes(xb_rank1) +
               sdsl::size_in_bytes(xb_sel0) +
               sdsl::size_in_bytes(xb_sel1) +
               sdsl::size_in_bytes(sl);
    }


    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    typename gcis_index_grid<wavelet_tree,XB_bit_vector,int_sequence>::ltype gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::first_column_point(
            const gcis_index_grid::dtype & col) const
    {
        return sl[sb.select(1, col)];
    }


    /**
     *  Recibe 2 points: X(Xrow , Xcol) and Y(Yrow , Ycol) identifying left upper corner and right bottom corner.
     *
     * */
    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    void gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::labels_search_2d(const gcis_index_grid::point & X,
                                                                                      const gcis_index_grid::point & Y,
                                                                                      std::vector<ltype> & labels) const
    {
        long p1 = 0, p2 = 0;
        /**
         * map the rows to the real position in the wavelet tree sequence
         * */
        p1 = map(X.first);
        p2 = map(Y.first + 1);  --p2;

        if (p1 > p2 || p2 < 0) return;

//        auto res = _wt_sb.range_search_2d2(p1, p2, X.second, Y.second);
        /** WT 2d orthogonal range search */
        auto res = sb.range_search_2d((dtype)p1, (dtype)p2, X.second, Y.second);
        /** Store just the labels of the points */
        for (auto point : res.second)
            labels.emplace_back(sl[sb.select(1,point.second)]);
    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    typename gcis_index_grid<wavelet_tree,XB_bit_vector,int_sequence>::dtype
    gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::map(const gcis_index_grid::dtype & row) const
    {
        long r = row;
        r = xb_sel1(r) - r + 1;
        return r;
    }

    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    void gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::build(std::vector<lpoint> & points, const dtype& max_rows,const dtype& max_cols )
    {


#ifdef DEBUG
        std::cout<<"gcis_index_grid::build"<<std::endl;
#endif
        dtype n_points = points.size();
//        dtype max_cols = 0, max_rows = 0;

        /**Find  max col and max row */
//        for (const auto &item : points) {
//            if( item.first.first > max_rows ) max_rows = item.first.first;
//            if( item.first.second > max_cols ) max_cols = item.first.second;
//        }


#ifdef DEBUG
        std::cout<<"n_points:"<<n_points<<std::endl;
        std::cout<<"max_cols:"<<max_cols<<std::endl;
        std::cout<<"max_rows:"<<max_rows<<std::endl;
#endif
        /**
         * Sort points by rows (rules) then by cols(suffix)
         * */
        auto begin = points.begin();
        auto end = points.end();

        sort(begin, end, [](const lpoint &a, const lpoint &b) -> bool {

            if ((a.first.first) < (b.first.first)) return true;

            if ((a.first.first) > (b.first.first)) return false;

            return (a.first.second) < (b.first.second);
        });

#ifdef DEBUG
        std::cout<<"sorted points:"<<n_points<<std::endl;
        for(uint i  = 0; i < 10; i++)
            std::cout<<"<"<<points[i].first.first<<","<<points[i].first.second<<">["<<points[i].second<<"]"<<std::endl;
#endif
        /**
         * Building bit_vector XB
         * XB store the row representation in unary
         * if |row[i]| = 4 => XB store 00001
         * */

        {

            std::vector<size_t> card_rows(max_rows, 0);
            /**
             * Computing the cardinal of every column and every row
             * */
            for (auto i = points.begin(); i != points.end(); ++i)
                card_rows[i->first.first - 1]++;

            sdsl::bit_vector _XB(n_points + max_rows + 1, 0);

            _XB[0] = true;
            /**
             * Put a 0 in XB for each element in the ith row
             * then add 1
             * */
            size_t p = 1;

            for (size_t j = 0; j < max_rows; ++j)
            {
                _XB[p + card_rows[j]] = true;
                p += card_rows[j] + 1;
            }

            xb = XB_bit_vector(_XB);

        }


        /**
         * Build a wavelet_tree on sb ( index of the columns not empty ) and plain representation for sl(labels)
         *
         * */

        {
            std::ofstream sb_file("sb_file", std::ios::binary);

            sdsl::int_vector<> _sl(n_points, 0);
            sdsl::int_vector<> _sb(n_points, 0);
            size_t j = 0;
            for (auto i = points.begin(); i != points.end(); ++i) {
                _sl[j] = i->second;
                _sb[j] = i->first.second;
                ++j;
            }

            sl = int_sequence(_sl);
            sdsl::util::bit_compress(sl);

            sdsl::util::bit_compress(_sb);
            sdsl::serialize(_sb, sb_file);
            sb_file.close();

            std::string id_sb = sdsl::util::basename("sb_file") + "_";
            sdsl::cache_config file_conf_sb(false, "./", id_sb);
            sdsl::construct(sb, "sb_file", file_conf_sb, 0);

        }


        build_rank_select_structures();


#ifdef DEBUG
        std::cout<<"xb:";
        for(uint i = 0 ; i < xb.size(); ++i)
            std::cout<<xb[i];
        std::cout<<std::endl;
        std::cout<<"sb:";
        for(uint i = 0 ; i < sb.size(); ++i)
            std::cout<<sb[i]<<" ";
        std::cout<<std::endl;
        std::cout<<"sl:";
        for(uint i = 0 ; i < sl.size(); ++i)
            std::cout<<sl[i]<<" ";
        std::cout<<std::endl;
#endif
    }


    template<typename wavelet_tree, typename XB_bit_vector, typename int_sequence>
    void gcis_index_grid<wavelet_tree, XB_bit_vector, int_sequence>::build_rank_select_structures() {

        xb_sel1 = typename XB_bit_vector::select_1_type(&xb);
        xb_sel0 = typename  XB_bit_vector::select_0_type(&xb);
        xb_rank1 =  typename XB_bit_vector::rank_1_type(&xb);

    }

};






#endif //GCIS_GCIS_INDEX_GRID_HPP
