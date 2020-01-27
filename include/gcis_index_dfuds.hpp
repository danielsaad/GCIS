//
// Created by ale on 14-01-20.
//

#ifndef GCIS_GCIS_INDEX_DFUDS_HPP
#define GCIS_GCIS_INDEX_DFUDS_HPP


#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/bp_support_sada.hpp>

namespace gcis_index_private{





    template<
            typename bv = sdsl::bit_vector,
            typename bp_sup = sdsl::bp_support_sada<>,
            typename rank_00_sup = sdsl::rank_support_v<00,2>,
            typename select_00_sup = sdsl::select_support_mcl<00,2>
    >
    class gcis_index_dfuds {

        typedef sdsl::bit_vector::select_0_type bv_select_0;


    public:

        typedef unsigned long ul;

        gcis_index_dfuds(){};

        gcis_index_dfuds(const gcis_index_dfuds& T)
        {
            bit_vector = T.bit_vector;
            bps         = bp_sup(&bit_vector);
            rank_00     = rank_00_sup (&bit_vector);
            select_00   = select_00_sup(&bit_vector);
            select_0    = bv_select_0(&bit_vector);
        };

        virtual ~gcis_index_dfuds(){};

        void print(){
            uint count = 0, pre = 0;
            for (uint i = 0; i < bit_vector.size() ; ++i) {
                if(bit_vector[i] == 0){
                    std::cout<<bit_vector[i]<<"["<<count<<"]"<<"<"<<i<<","<<++pre<<">"<<std::endl;
                    count = 0;
                }else{
                    count++;
                    std::cout<<bit_vector[i];
                }


            }
            std::cout<<std::endl;

        }

        void set_tree( const sdsl::bit_vector& b)
        {
            bit_vector = b;
            bps         = bp_sup(&bit_vector);
            rank_00     = rank_00_sup (&bit_vector);
            select_00   = select_00_sup(&bit_vector);
            select_0    = bv_select_0(&bit_vector);
        }

        void build( const sdsl::bit_vector &v)
        {
            auto _bv = sdsl::bit_vector(v.size()+3);
            _bv[0]=1;_bv[1]=1;_bv[2]=0;

            for (int i = 0; i < v.size(); ++i) {
                _bv[3+i] = v[i];
            }

            bit_vector  = bv(_bv);
            bps         = bp_sup(&bit_vector);
            rank_00     = rank_00_sup (&bit_vector);
            select_00   = select_00_sup(&bit_vector);
            select_0    = bv_select_0(&bit_vector);
        }

        /*
         * tree operations
         * */

        ul root()const{ return 3ul;}

        ul pre_order( const ul & v ) const {  return  v-bps.rank(v)+bit_vector[v] ; }

        ul pre_order_select( const ul & i) const { return select_0(i)+1; }

        ul operator[](const ul& i)const { return select_0(i)+1;}

        ul parent( const ul & v) const
        {

            if(v == 3) return 0;

            return pred0( bps.find_open(v-1) )+1;
        }

        ul ichild( const ul & v, const ul & t) const
        {
            if(!bit_vector[v]) return 0;

            size_t close = bps.find_close(succ0(v)-t);

            if( close == bps.size() ) return 0;

            return close + 1;

        }

        ul fchild( const ul & v) const{ return succ0(v)+1;};
//
//    ul lchild( const ul &, const ul &) const;

        ul nsibling( const ul & v) const {  return bps.fwd_excess(v-1,-1)+1; }

        ul children (const ul & v) const
        {
            if(!bit_vector[v])
                return  0;

            size_t zeros = v - bps.rank(v)+1 ;

            return select_0(zeros+1) - v;
        }

        ul childrank (const ul & v) const
        {
            if(v == 3) return 0;

            size_t bp_parent = bps.find_open(v-1);

            return succ0(bp_parent) - bp_parent;
        }

        ul is_leaf (const ul & v) const { return !bit_vector[v]; }

        bool is_ancestor (const ul & u, const ul & v) const
        {
            return (u < v) && (v < bps.fwd_excess(u-1,-1));
        }

        ul lca( ul & u,  ul & v) const {
            if(u > v) std::swap(u,v);
            return bps.rmq( succ0(u), v - 1) + 1;
//            ul _lca = u;
//
//            while(!is_ancestor(_lca,v)){
//                _lca = parent(_lca);
//            }
//
//            return _lca;
        }

        ul lcu(const ul & u, const ul & v) const {
            /*
             * preorder(u) < preorder(v)
             * */

            ul _lca = u;
            ul p_lca = parent( _lca);

            while(p_lca != 3 && !is_ancestor( p_lca , v)){
                _lca = p_lca;
                p_lca = parent(_lca);
            }

            return _lca;
        }


        /*
         *
         * leaf operations
         *
         * */

        ul leaf_rank (const ul & v) const {
            return rank_00(v+1);
        }

        ul leaf_select (const ul & i) const { return select_00(i);}

        ul leaf_num (const ul & v) const { return leaf_rank( bps.fwd_excess(v-1,-1) +1)-leaf_rank(v);}

        ul last_leaf (const ul & v) const { return leaf_rank(v) + leaf_num(v) - 1;}

        ul next_leaf(const ul & i) const { return select_00(i+1);}



        void load(std::ifstream &f)
        {
            sdsl::load(bit_vector,f);
            bps         = bp_sup(&bit_vector);
            rank_00     = rank_00_sup (&bit_vector);
            select_00   = select_00_sup(&bit_vector);
            select_0    = bv_select_0(&bit_vector);
        }

        void serialize(std::ofstream & f) const
        {
            sdsl::serialize(bit_vector,f);
        }


        gcis_index_dfuds& operator= (const gcis_index_dfuds& T){
            bit_vector = T.bit_vector;
            bps         = bp_sup(&bit_vector);
            rank_00     = rank_00_sup (&bit_vector);
            select_00   = select_00_sup(&bit_vector);
            select_0    = bv_select_0(&bit_vector);
        }

        size_t size_in_bytes()const{
            return sdsl::size_in_bytes(bit_vector)+
                   sdsl::size_in_bytes(bps)+
                   sdsl::size_in_bytes(rank_00)+
                   sdsl::size_in_bytes(select_0)+
                   sdsl::size_in_bytes(select_00);

        }

    protected:

        bv bit_vector;
        bv_select_0 select_0;


        bp_sup bps;
        rank_00_sup rank_00;
        select_00_sup select_00;

        ul pred0(const ul & i) const
        {
            if(bit_vector[i]==0)
                return i;
            return select_0(i - bps.rank(i)+1);
        }

        ul succ0(const ul& i) const
        {
            size_t zeros = i - bps.rank(i)+1;
            return select_0(zeros+bit_vector[i]);
        }

    };



}



#endif //GCIS_GCIS_INDEX_DFUDS_HPP
