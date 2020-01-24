//
// Created by ale on 10-01-20.
//

#ifndef GCIS_GCIS_INDEX_HPP
#define GCIS_GCIS_INDEX_HPP

#include <string>
#include <vector>
#include <sdsl/inv_perm_support.hpp>
#include "gcis_index_grid.hpp"
#include "gcis_index_dfuds.hpp"

#define TOYEXAMPLE 1
#define BUFFER_SIZE 1000;
#define MSIGMA 256;

namespace gcis_index_private{


    template < typename , typename , typename, typename, typename  > class gcis_index;

    /**
    * IS-based compressed index
    *
    **/
    template <
                typename t_mapfbv = sdsl::sd_vector<> ,
                typename t_maptbv = sdsl::sd_vector<>,
                typename t_mapwt  = sdsl::wt_gmr<>,
                typename t_gridbv = sdsl::rrr_vector<>,
                typename t_gridwt = sdsl::wt_int<>
                >

    class gcis_index {


        public:

            typedef long long len_type; // type for len of the text
            typedef uint64_t rule_type; // type for rules
            typedef gcis_index_dfuds<>::ul size_tree;
            typedef gcis_index_grid<t_gridwt,t_gridbv,sdsl::int_vector<>> t_grid; // typdef for grid
            typedef std::pair<size_tree,len_type> _node; // type to track secondary occ



            gcis_index();
            gcis_index( const gcis_index& );
//            gcis_index& operator= (const gcis_index &);
            virtual ~gcis_index();

            virtual void serialize(std::ofstream& ) const;
            virtual void load(std::ifstream& );

            void set_grid(const t_grid& );
            void set_tree(const sdsl::bit_vector & );
            void set_bvfocc( const t_mapfbv& );
            void set_bvt( const t_maptbv & );
            void set_vt( const std::string& );
            void set_pi( const sdsl::int_vector<>& );
            void set_nt( const t_mapwt &);
            void set_l(const sdsl::bit_vector & );
            void set_sigma( const uint& );

            /** Return the positions of the text where the pattern occurs */
//            virtual void locate(const std::string & s, std::vector<len_type>& )const  = 0;
            /** Return the substring T[i...i+m)*/
            virtual void display(const len_type& , const len_type& , std::string& ) const ;


        protected:

            /** grid for primary occurrences */
            t_grid _grid;
            /** bitvector for first time nodes occ */
            t_mapfbv bvfocc;
            typename t_mapfbv::select_1_type bvfocc_sel1;
            typename t_mapfbv::select_0_type bvfocc_sel0;
            typename t_mapfbv::rank_1_type bvfocc_rank1;
            /** permutation for first occ node rules associated to 1 in bvfocc*/
            sdsl::int_vector<> pi;
            sdsl::inv_perm_support<> inv_pi;
            /** bitvector for characters nodes in seccondary nodes associated to 0 in bvfocc*/
            t_maptbv bvt;
            typename t_maptbv::rank_1_type bvt_rank1;
            typename t_maptbv::select_0_type bvt_sel0;
            /** sequence for secondary occ nodes labels that aren't characters associated to 0 in bvt*/
            t_mapwt wtnt;
            /** plain array for sequence of characters associated to 1 in bvt*/
            std::string vt;
            /** Tree loud topology */
            gcis_index_dfuds<> dfuds_tree;
            /** offset of leaves in the text*/
            sdsl::sd_vector<> bvl;
            sdsl::sd_vector<>::rank_1_type bvl_rank1;
            sdsl::sd_vector<>::select_1_type bvl_select1;

            uint sigma;
            uint32_t n_rules;
            uint32_t n_suffixes;

        public:

            /**Compare a rule suffix with an string left partition*/
            virtual int cmp_suffix_rule(const rule_type& , const std::string &,  len_type& j) const;
            virtual int cmp_suffix_rule_pre(const size_tree & , const std::string &,  len_type j) const;
            /**Compare a rule suffix with an string left partition*/
            virtual int cmp_prefix_rule(const rule_type& , const std::string &,  len_type  &  ) const;

            /**Compare a grammar suffix with an string right partition*/
            virtual int cmp_suffix_grammar(size_tree , const std::string &,  len_type & j) const;

//            virtual void dfs_expand_rule(const rule_type& X, std::string& s, len_type& off, const len_type& i, const len_type& m, bool suffix) const;

            virtual bool expand_prefix(const size_tree& node, std::string& s,const len_type& m, const size_tree & off) const;
            /**Implement an strategy to find the range of valid points of the grid*/
//            virtual void find_range(const std::string & s, std::vector<len_type>& ) const  = 0;
            /**Algorithm to recursively track secondary occ in the grammar*/
            virtual void find_secondary_occ(const _node &, std::vector<len_type>& ) const ;


            /** given a rule return the node of its first occ */
            size_tree map_rule(const rule_type& ) const;
            /** given a node return the label rule*/
            rule_type map_node(const size_tree& ) const;

            /** compute the offset of a node in the text*/
            len_type offset_node( const size_tree& )const;
            /** true if the  node is a terminal node X -> "ACGTTT" expand just characters nodes*/
            bool is_t_node(const size_tree& )const;
            /** true if the node is a character node*/
            bool is_ch_node(const size_tree&)const;
            bool is_ch_node(const size_tree&, size_tree&)const;
            /** take a terminal node and append the block of it and its siblings to the string*/
            void append_block( const size_tree&, const size_tree&, std::string&) const;
            /** take a character node and compare the block of it and its siblings to the string*/
            int cmp_block( len_type pos,  uint l,const std::string& str,  len_type & i, const bool& dir) const;

    };

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::gcis_index(const gcis_index &) {

    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::~gcis_index() = default;

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::serialize(std::ofstream &f) const {
        _grid.serialize(f);
        sdsl::serialize(bvfocc,f);
        sdsl::serialize(pi,f);
        sdsl::serialize(bvt,f);
        sdsl::serialize(wtnt,f);
        dfuds_tree.serialize(f);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::load(std::ifstream & f) {

        _grid.load(f);
        sdsl::load(bvfocc,f);
        sdsl::load(pi,f);
        sdsl::load(bvt,f);
        sdsl::load(wtnt,f);
        dfuds_tree.load(f);
    }

    /***
     * @brief This function return the preorder of the node where the first occ of rule appears,
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @param X rule id
     * @return preorder of the node
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<>::size_tree
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::map_rule(const gcis_index::rule_type & X) const {

        unsigned long p = inv_pi[X];
        auto pre = bvfocc_sel1.select(p+1);
        return pre+1;
    }
    /**
     * @brief This function return the rule id of the node with preorder p.
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @return
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<>::rule_type
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::map_node(const gcis_index::size_tree & p) const {

        if(bvfocc[p-1]){
            //primary node
            auto r = bvfocc_rank1.rank(p);
            return pi[r-1];
        }
        //secondary node
        auto r0 = p - bvfocc_rank1(p);
        if(bvt[r0-1]){
            //character symbol
            auto sr1 = bvt_rank1.rank(r0);
            return vt[sr1-1];

        }
        //non-terminal
        auto sr0 = r0 - bvt_rank1.rank(r0);
        return wtnt[sr0-1];
    }
    /**
     * @brief extract a substring T[i...i+m]
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @param i initial position
     * @param m len of the substring to extract
     * @param str result
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::display(const gcis_index::len_type & ii,
                                                                              const gcis_index::len_type & m,
                                                                              std::string & str) const
    {

        len_type i = ii + sigma;
        /**found leaf containing the position rank [1...i+1) */
        len_type leaf_boundary = bvl_rank1.rank(i+1);
        size_tree leaf_pos = bvl_select1(leaf_boundary);
        len_type off = i-leaf_pos;
        /**found the associated leaf node */
        size_tree node = dfuds_tree.leaf_select(leaf_boundary);
        // in case that is a single character

        std::string sblock;
        sblock.reserve(m);
        long strs = sblock.size();
        size_tree pos = 0,parent = 0,ch = 0, j = 0 ;
        // extraction process bottom up
        parent = dfuds_tree.parent(node);
        ch = dfuds_tree.children(parent);
        j = dfuds_tree.childrank(node);

        if(is_ch_node(node,pos)){

            append_block(pos-1,ch-j+1,sblock);
            strs = sblock.size();

            do{
                node = dfuds_tree.parent(node);
                parent = dfuds_tree.parent(node);
                ch = dfuds_tree.children(parent);
                j = dfuds_tree.childrank(node);

            }while(j == ch);
            node = dfuds_tree.nsibling(node);
            ++j;
            /// cambiar aqui por extraer hasta el final;;

        }
        while(  strs - off < m ){

            while(j <= ch && expand_prefix(node,sblock,m,off)){
                node = dfuds_tree.nsibling(node); ++j;
            }
            if( j > ch ){
                node = parent;
                parent = dfuds_tree.parent(node);
                ch = dfuds_tree.children(parent);
                j = dfuds_tree.childrank(node) + 1;
                if( j <= ch) node = dfuds_tree.nsibling(node);
            }else{
                strs = sblock.size();
            }
        }
        /** remove innecesary characters */
        str.resize(m);
        std::copy(sblock.begin() + off,sblock.begin() + off + m , str.begin());

    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<>::len_type
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::offset_node(const gcis_index::size_tree & node) const {
        long long lrank = dfuds_tree.leaf_rank(node-1)+1;
//        if(dfuds_tree.is_leaf(node))
//            --lrank;
        return bvl_select1.select(lrank);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_prefix(const gcis_index::size_tree &node,
                                                                                    std::string &s,
                                                                                    const gcis_index::len_type &m,
                                                                                    const size_tree& off) const{
        /** if the node is not a leaf */
        if(!dfuds_tree.is_leaf(node))
        {
            /** check if its first child is a character node */
            size_tree fchild = dfuds_tree.fchild(node);
            size_tree l = dfuds_tree.children(node);
            // if its a leaf and a character node then al the children are character and we extract the block
            size_tree r = 0;
            if(dfuds_tree.is_leaf(fchild) && is_ch_node(fchild,r)){
                // append the corresponding block of characters
                size_tree pos = bvt_rank1.rank(r);
                append_block(pos-1,l,s);
                return s.size()-off < m; // continue expanding
            }
            // if it is no a character node -> expand in dfs
            uint i = 1;
            while( i <= l )
            {
                if(!expand_prefix(fchild,s,m,off))
                    return false;
                fchild = dfuds_tree.ichild(node,++i);
            }
            return true;
        }
        // leaf case !!!!! WE DO NOT TREAT CHARACTER LEAVES
        // as we never arrive to character leaves the leaf must to be a non-terminal leaf
        // we need to jump to its first occ
        size_tree pre = dfuds_tree.pre_order(node);
        rule_type X = map_node(pre);
        size_tree preorder_focc= map_rule(X);

        return expand_prefix(dfuds_tree.pre_order_select(preorder_focc),s,m,off);
    }


    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_suffix_rule_pre(const size_tree & pre,
                                                                                     const std::string & str,
                                                                                     len_type j) const{
        rule_type X = map_node(pre);
        size_tree node = dfuds_tree.pre_order_select(pre);
        //case leaf
        if(  dfuds_tree.is_leaf(node)  )  //if it is a leaf then is a symbol
        {
            if(str[j--] == X) return 0;
            if(str[j--] > X) return 1;

            return -1;

        }
        // lambda function for recursively cmp
        std::function< int (const size_tree&, const size_tree&) > dfs_cmp_suffix;

        dfs_cmp_suffix  = [this, &dfs_cmp_suffix, &str, &j ](const size_tree& node,const size_tree& pre_node){
            // !!!!!!WE NEVER REACH A CHARACTER NODE
            // case leaf
            if(dfuds_tree.is_leaf(node)) // as we never reach a character node if we reach a leaf it has to be non-terminal leaf
            {   //jump to it definition
                rule_type Y = map_node(pre_node);
                size_tree pre = map_rule(Y);
                size_tree def_node = dfuds_tree.pre_order_select(pre);
                int r = dfs_cmp_suffix(def_node,pre);
                return r;
            }

            // check if it is a terminal node X->"AJSKLDJALKSJD"

            size_tree ch = dfuds_tree.children(node);
            size_tree lchild = dfuds_tree.ichild(node,ch);
            size_tree pre_lchild =  dfuds_tree.pre_order(lchild);

            if(dfuds_tree.is_leaf(lchild)){
                auto r0 = pre_lchild - bvfocc_rank1(pre_lchild);
                if( bvt[r0-1] ){ // is a terminal node (just compare characters) and return
                    auto r1 = bvt_rank1.rank(r0);
                    return cmp_block(r1-1,(uint) ch,str,j,true);
                }
            }
            //

            // compare first child, is already precomputed
            int r = dfs_cmp_suffix(lchild,pre_lchild);
            if( r != 0) return r; // return if missmatch

            uint i = ch-1;
            // compare the rest of children
            while(i > 0)
            {
                lchild = dfuds_tree.ichild(node,i);
                pre_lchild = dfuds_tree.pre_order(lchild);
                r = dfs_cmp_suffix(lchild,pre_lchild);
                if( r != 0) return r;
                --i;
            }
            return 0;
        };


        return dfs_cmp_suffix(node,pre);

    }
        /**
         * @brief compare a suffix of the rule with the left part of pattern partition in reverse order cmp(X,str,j) rev(X) ?
         * rev(str[1...j])
         * @tparam t_mapfbv
         * @tparam t_maptbv
         * @tparam t_mapwt
         * @tparam t_gridbv
         * @tparam t_gridwt
         * @param X rule to compare
         * @param str pattern to compare
         * @param j suffix pos by reference return the number of ch drains
         * @return 0 if the left partition is a suffix of the rule
         */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_suffix_rule(const gcis_index::rule_type & X,
                                                                                      const std::string & str,
                                                                                      len_type& j) const{


        // compute the preorder node and node pos
        size_tree pre = map_rule(X);
        size_tree node = dfuds_tree.pre_order_select(pre);
        //case leaf
        if(  dfuds_tree.is_leaf(node)  )  //if it is a leaf then is a symbol
        {
#ifdef TOYEXAMPLE
            char Y = X;
            switch(Y){
                case 0: Y = 65;
                    break;
                case 1: Y = 67;
                    break;
                case 10: Y = 71;
                    break;
                case 13: Y = 84;
                    break;
            }


            if(str[j] < Y) return -1;
            if(str[j] > Y) return 1;
            --j;
            return 0;
#endif

            if(str[j] < Y) return -1;
            if(str[j] > Y) return 1;
            --j;
            return 0;

        }
        // lambda function for recursively cmp
        std::function< int (const size_tree&, const size_tree&) > dfs_cmp_suffix;

        dfs_cmp_suffix  = [this, &dfs_cmp_suffix, &str, &j ](const size_tree& node,const size_tree& pre_node){
            // !!!!!!WE NEVER REACH A CHARACTER NODE
            // case leaf
            if(dfuds_tree.is_leaf(node)) // as we never reach a character node if we reach a leaf it has to be non-terminal leaf
            {   //jump to it definition
                rule_type Y = map_node(pre_node);
                size_tree pre = map_rule(Y);
                size_tree def_node = dfuds_tree.pre_order_select(pre);
                int r = dfs_cmp_suffix(def_node,pre);
                return r;
            }

            // check if it is a terminal node X->"AJSKLDJALKSJD"

            size_tree ch = dfuds_tree.children(node);
            size_tree lchild = dfuds_tree.ichild(node,ch);
            size_tree pre_lchild =  dfuds_tree.pre_order(lchild);

            if(dfuds_tree.is_leaf(lchild)){
                auto r0 = pre_lchild - bvfocc_rank1(pre_lchild);
                if( bvt[r0-1] ){ // is a terminal node (just compare characters) and return
                    auto r1 = bvt_rank1.rank(r0);
                    return cmp_block(r1-1,(uint) ch,str,j,true);
                }
            }
            //

            // compare first child, is already precomputed
            int r = dfs_cmp_suffix(lchild,pre_lchild);
            if( r != 0) return r; // return if missmatch

            uint i = ch-1;
            // compare the rest of children
            while(j >= 0 && i > 0 )
            {
                lchild = dfuds_tree.ichild(node,i);
                pre_lchild = dfuds_tree.pre_order(lchild);
                r = dfs_cmp_suffix(lchild,pre_lchild);
                if( r != 0) return r;
                --i;
            }
            return 0;
        };


        return dfs_cmp_suffix(node,pre);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_t_node(const gcis_index::size_tree & node ) const {
        // THE NODE MUST NOT BE A LEAF
        assert(!dfuds_tree.is_leaf(node));
        size_tree fchild = dfuds_tree.fchild(node);
        return is_ch_node(fchild);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_ch_node(const gcis_index::size_tree & node) const {
        assert(dfuds_tree.is_leaf(node));
        size_tree pre = dfuds_tree.pre_order(node);
        auto r0 = pre - bvfocc_rank1(pre);
        return bvt[r0-1];
    }
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_ch_node(const gcis_index::size_tree & node, size_tree& idx) const {
        assert(dfuds_tree.is_leaf(node));
        size_tree pre = dfuds_tree.pre_order(node);
        idx = pre - bvfocc_rank1(pre);
        return bvt[idx-1];
    }

    /**
     * @brief this function performs pattern matching in the passed direction
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @param pos in the character string
     * @param l number of characters to analyse
     * @param str string to compare
     * @param i point partition string
     * @param dir direction to compare
     * @return 0 if the left/right partition is suffix/prefix of the block
     *         1 if the left/right partition is greater
     *        -1 if is left
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_block( gcis_index::len_type pos,
                                                                                uint l,
                                                                               const std::string &str,
                                                                               len_type & i,
                                                                               const bool &dir) const {
        if(dir){
            assert(l <= pos+1);
            while(l && i >= 0){
                if(str[i] < vt[pos]) return -1;
                if(str[i] > vt[pos]) return 1;
                --pos;
                --i;
                --l;
            }

            //if( i < 0 && l != 0 ) return -1;

            return 0;
        }
        assert( pos + l <= vt.size());

        uint ssz = str.size();
        while(l && i < ssz){
            if(str[i] < vt[pos]) return -1;
            if(str[i] > vt[pos]) return 1;
            ++pos;
            ++i;
            --l;
        }
        //if( i == ssz && l != 0 ) return -1;

        return 0;
    }
    /**
     * @brief same that cmp_suffix_rule but with prefix heheheh
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @param X
     * @param str
     * @param j by reference return how many character where drain
     * @return 0 if the right partition is prefix of the rule
     *          -1 if it is less
     *          1 if it is greater
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_prefix_rule(const gcis_index::rule_type & X,
                                                                                     const std::string & str,
                                                                                     len_type  &j) const{


        // compute the preorder node and node pos
        size_tree pre = map_rule(X);
        size_tree node = dfuds_tree.pre_order_select(pre);
        //case leaf
        if(  dfuds_tree.is_leaf(node)  )  //if it is a leaf then is a symbol
        {
#ifdef TOYEXAMPLE
            char Y = X;
            switch(Y){
                case 0: Y = 65;
                    break;
                case 1: Y = 67;
                    break;
                case 10: Y = 71;
                    break;
                case 13: Y = 84;
                    break;
            }


            if(str[j] < Y) return -1;
            if(str[j] > Y) return 1;
            j++;
            return 0;
#endif

            if(str[j] < Y) return -1;
            if(str[j] > Y) return 1;
            j++;
            return 0;

        }
        // lambda function for recursively cmp
        std::function< int (const size_tree&, const size_tree&) > dfs_cmp_prefix;

        dfs_cmp_prefix  = [this, &dfs_cmp_prefix, &str, &j ](const size_tree& node,const size_tree& pre_node)->int{
            // !!!!!!WE NEVER REACH A CHARACTER NODE
            // case leaf
            if(dfuds_tree.is_leaf(node)) // as we never reach a character node if we reach a leaf it has to be non-terminal leaf
            {   //jump to it definition
                rule_type Y = map_node(pre_node);
                size_tree pre = map_rule(Y);
                size_tree def_node = dfuds_tree.pre_order_select(pre);
                uint r = dfs_cmp_prefix(def_node,pre);
                return r;
            }

            // check if it is a terminal node X->"AJSKLDJALKSJD"

            size_tree ch = dfuds_tree.children(node);
            size_tree fchild = dfuds_tree.fchild(node);
            size_tree pre_fchild =  pre_node+1;
            // ask if fchild is a leaf
            if(dfuds_tree.is_leaf(fchild)){
                auto r0 = pre_fchild - bvfocc_rank1(pre_fchild);
                if( bvt[r0-1] ){ // is a terminal node (just compare characters) and return
                    auto r1 = bvt_rank1.rank(r0);
                    int r = cmp_block(r1-1,(uint) ch,str,j,false);
                    return r;
                }
                //

            }

            // compare first child is already precomputed
            uint r = dfs_cmp_prefix(fchild,pre_fchild);
            if( r != 0) return r; // return if missmatch

            uint i = 2;
            // compare the rest of children
            while( j < str.size() && i <= ch)
            {
                fchild = dfuds_tree.ichild(node,i);
                pre_fchild = dfuds_tree.pre_order(fchild);
                r = dfs_cmp_prefix(fchild,pre_fchild);
                if( r != 0) return r;
                ++i;
            }
            return 0;
        };


        return dfs_cmp_prefix(node,pre);
    }

    /**
     * @brief compare a suffix of a rule with the right partition of the string
     * @tparam t_mapfbv
     * @tparam t_maptbv
     * @tparam t_mapwt
     * @tparam t_gridbv
     * @tparam t_gridwt
     * @param pre_node
     * @param str
     * @param j
     * @return 0 if the right partition is a prefix of suffix starting in the node
     *          -1 if is less
     *          1 if is greater
     */
    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_suffix_grammar( gcis_index::size_tree pre_node,
                                                                                        const std::string & str,
                                                                                         len_type  &j) const{

        size_tree node = dfuds_tree.pre_order_select(pre_node);
        size_tree parent = dfuds_tree.parent(node);
        size_tree ch = dfuds_tree.children(parent);
        size_tree rch = dfuds_tree.childrank(node);
        rule_type Y = map_node(pre_node);

#ifdef TOYEXAMPLE

        switch(Y){
            case 65: Y = 0;
                break;
            case 67: Y = 1;
                break;
            case 71: Y = 10;
                break;
            case 84: Y = 13;
                break;
        }
#endif
        int r = cmp_prefix_rule(Y,str,j);
//        if(j == str.size()) return 0;
        if(r != 0) return r;
//        if(r == 0 && j < str.size()) return 1;

        for (size_tree i = rch+1; j > 0 && i <= ch ; ++i) {

            node = dfuds_tree.ichild(parent,i);
            pre_node = dfuds_tree.pre_order(node);
            Y = map_node(pre_node);

#ifdef TOYEXAMPLE

            switch(Y){
                case 65: Y = 0;
                    break;
                case 67: Y = 1;
                    break;
                case 71: Y = 10;
                    break;
                case 84: Y = 13;
                    break;
            }
#endif
            r = cmp_prefix_rule(Y,str,j);
//            if(j == str.size()) return 0;
            if(r != 0) return r;
//            if(r == 0 && j < str.size()) return 1;
        }

        if(r == 0 && j < str.size()) return 1;

        return 0;
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::find_secondary_occ(const _node & v,
                                                                                         std::vector<len_type> & occ) const{

        struct secc_node{
            size_tree pre;
            size_tree pos;
            rule_type X;
            long long off;
            secc_node(const size_tree & _pre, const size_tree & _pos,const rule_type & _x , const len_type & _off):pre(_pre),pos(_pos),X(_x), off(_off){ }
        };


        std::deque<secc_node> S;

        {
            size_tree node = dfuds_tree.pre_order_select(v.first);
            size_tree parent = dfuds_tree.parent(node);
            size_tree pre_parent = dfuds_tree.pre_order(parent);
            // this rule has to be a non-terminal so we search it in the golynski sequence
            rule_type X = map_node(pre_parent);
            long long off =  v.second +  offset_node(node) - offset_node(parent);

            size_t nocc = this->wtnt.rank(wtnt.size(),X);
            secc_node ts(pre_parent, parent, X, off);
            S.emplace_back(ts);

            for (int i = 1; i <= nocc ; ++i) {

                auto ps = wtnt.select(i,X)+1;
                auto bvt0 = bvt_sel0.select(ps)+1;
                auto bvfocc0 = bvfocc_sel0(bvt0);
                size_tree tnode = dfuds_tree.pre_order_select(bvfocc0+1);
                secc_node t(bvfocc0+1, tnode, X, off);
                S.emplace_back(t);
            }
        }

        while(!S.empty())
        {
            if(S.front().pre == 1) occ.push_back((len_type) S.front().off);
            else
                {
                    size_tree parent = dfuds_tree.parent(S.front().pos);
                    size_tree pre = dfuds_tree.pre_order(parent);
                    rule_type X = map_node(pre);
                    long long off =  S.front().off +  offset_node(S.front().pos) - offset_node(parent);
                    secc_node ts(pre,parent,X,off);
                    S.push_back(ts);
                    size_tree nocc = this->wtnt.rank(wtnt.size(),X);
                    for (int i = 1; i <= nocc ; ++i) {

                        auto ps = wtnt.select(i,X) + 1;
                        auto bvt0 = bvt_sel0.select(ps)+1;
                        auto bvfocc0 = bvfocc_sel0(bvt0);
                        size_tree tnode = dfuds_tree.pre_order_select(bvfocc0+1);
                        secc_node t(bvfocc0+1, tnode, X, off);
                        S.emplace_back(t);
                    }
                }
            S.pop_front();
        }
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_tree(const sdsl::bit_vector & b) {
         dfuds_tree.set_tree(b);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_bvfocc(const t_mapfbv & b ) {

        bvfocc = b;
        bvfocc_sel1  = typename t_mapfbv::select_1_type(&bvfocc) ;
        bvfocc_sel0  = typename t_mapfbv::select_0_type(&bvfocc) ;
        bvfocc_rank1 = typename t_mapfbv::rank_1_type(&bvfocc) ;


    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_bvt(const t_maptbv & b) {
        bvt = b;
        bvt_rank1 = typename t_maptbv::rank_1_type (&bvt);
        bvt_sel0  = typename t_maptbv::select_0_type (&bvt);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_vt(const std::string & str) {
        vt = str;
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_pi(const sdsl::int_vector<> & _pi) {
        pi = _pi;
        inv_pi = sdsl::inv_perm_support<>(&pi);
        n_rules = pi.size();
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::gcis_index() = default;

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_grid(const gcis_index::t_grid & g) {
        _grid = g;
        n_suffixes = _grid.n_cols();
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_nt(const t_mapwt & nt) {  wtnt = nt; }


    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::append_block(const gcis_index::size_tree & pos,
                                                                                   const gcis_index::size_tree & l,
                                                                                   std::string & str
                                                                                   ) const
   {
        std::string s;
        s.resize(l);
        std::copy(vt.begin()+pos,vt.begin()+pos+l, s.begin());
        str += s;

    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_l(const sdsl::bit_vector & _l) {
        bvl = _l;
        bvl_rank1 = sdsl::sd_vector<>::rank_1_type (&bvl);
        bvl_select1 = sdsl::sd_vector<>::select_1_type (&bvl);
    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_sigma(const uint &s) {
        sigma = s;
    }




//    /**
//     *
//     * @tparam _gcis
//     * @param v is a pair< node bfs-order, initial offset occ >
//     * @param occ is the vector for storing occ
//     */
//    template<typename _gcis>
//    void gcis_index<_gcis>::find_secondary_occ(const gcis_index::_node & v, std::vector<len_type> & occ) const
//    {
//        std::deque<_node> S;
//
//        {
//            typename _gcis::rule X = sequence_tree[v.first];
//            t_wt::size_type n_occ = sequence_tree.rank(sequence_tree.size(),v.first);
//
//            /**
//             * Initialize the deque with all the occ of v
//             * */
//            S.push_back(v);
//            for (int i = 1; i < n_occ ; ++i)
//                S.emplace_back(sequence_tree.select(i+1,v.first)+1,v.second);
//
//        }
//
//        while(!S.empty())
//        {
//            /** base case we reach the root of the tree */
//            if(S.front().first == 1) occ.push_back(S.front().second);
//            else
//                {
//                    /**
//                    * find the parent of the current node
//                    * */
//                    size_tree vp = loud_tree.node_pos(S.front().first);
//                    size_tree parent_v;
//                    size_tree ch_rank = loud_tree.child_rank(vp,parent_v);
//                    t_wt::size_type label = sequence_tree[loud_tree.node_id(parent_v)];
//                    /** Compute offset and insert the node to a future processing*/
//                    auto off = S.front().second + _grammar->offset(label,ch_rank);
//                    S.emplace_back(label,(len_type)off);
//
//                    /** Compute total occ of the node */
//                    t_wt::size_type n_occ = sequence_tree.rank( sequence_tree.size(),label);
//
//                    for (int i = 1; i < n_occ ; ++i)
//                        S.emplace_back(sequence_tree.select(i+1,label)+1,off);
//                }
//
//            /**Remove the processed node*/
//            S.pop_front();
//        }
//
//    }

}







#endif //GCIS_GCIS_INDEX_HPP
