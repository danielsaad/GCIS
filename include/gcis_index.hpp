//
// Created by ale on 10-01-20.
//

#ifndef GCIS_GCIS_INDEX_HPP
#define GCIS_GCIS_INDEX_HPP

#include "gcis_index_dfuds.hpp"
#include "gcis_index_grid.hpp"
#include <sdsl/inv_perm_support.hpp>
#include <string>
#include <vector>

// TODO: Remove
using std::cout;
using std::endl;

//#define TOYEXAMPLE 1

namespace gcis_index_private {

template <typename, typename, typename, typename, typename> class gcis_index;

/**
 * IS-based compressed index
 *
 **/
template <
    typename t_mapfbv = sdsl::sd_vector<>,
    typename t_maptbv = sdsl::sd_vector<>, typename t_mapwt = sdsl::wt_gmr<>,
    typename t_gridbv = sdsl::rrr_vector<>, typename t_gridwt = sdsl::wt_int<>>

class gcis_index {

  public:
    typedef long long len_type; // type for len of the text
    typedef uint64_t rule_type; // type for rules
    typedef gcis_index_dfuds<>::ul size_tree;
    typedef gcis_index_grid<t_gridwt, t_gridbv, sdsl::int_vector<>>
        t_grid;                                   // typdef for grid
    typedef std::pair<size_tree, len_type> _node; // type to track secondary occ

    gcis_index();
    gcis_index(const gcis_index &);
    //            gcis_index& operator= (const gcis_index &);
    virtual ~gcis_index();

    virtual void serialize(std::ofstream &) const;
    virtual void load(std::ifstream &);
    virtual size_t size_in_bytes() const;

    void set_grid(const t_grid &);
    void set_tree(const sdsl::bit_vector &);
    void set_bvfocc(const t_mapfbv &);
    void set_bvt(const t_maptbv &);
    void set_vt(const std::string &);
    void set_pi(const sdsl::int_vector<> &);
    void set_nt(const t_mapwt &);
    void set_l(const sdsl::sd_vector<> &);
    void set_sigma(const uint &);

    /** Return the positions of the text where the pattern occurs */
    //            virtual void locate(const std::string & s,
    //            std::vector<len_type>& )const  = 0;
    /** Return the substring T[i...i+m)*/
    virtual void display(const len_type &, const len_type &,
                         std::string &) const;

    virtual void displayII(const len_type &, const len_type &,
                         std::string &) const;

    virtual void display_L(const std::size_t &i, const std::size_t &len, std::string &str) const{
        str.resize(len);
        size_t p = 0;
        expand_interval(dfuds_tree.root(), std::make_pair(i, i+len-1), str, p);
    }

    virtual bool extract_prefix_with_gap
    (
            const size_tree& node,
            std::string& str,
            long& m_off,
            const uint &gap
    ) const;

    virtual void print_size_in_bytes() const {

        std::cout << "Printing size in bytes" << std::endl;

        std::cout << "Mapping Structures:" << std::endl;
        std::cout << "bvfocc:" << sdsl::size_in_bytes(bvfocc) << std::endl;
        std::cout << "bvfocc_rank1:" << sdsl::size_in_bytes(bvfocc_rank1)
                  << std::endl;
        std::cout << "bvfocc_sel0:" << sdsl::size_in_bytes(bvfocc_sel0)
                  << std::endl;
        std::cout << "bvfocc_sel1:" << sdsl::size_in_bytes(bvfocc_sel1)
                  << std::endl;

        std::cout << "-----------" << std::endl;

        std::cout << "pi:" << sdsl::size_in_bytes(pi) << std::endl;
        auto t_pi = pi;
        sdsl::util::bit_compress(t_pi);
        std::cout << "pi(bit_compress):" << sdsl::size_in_bytes(t_pi)
                  << std::endl;
        std::cout << "inv_pi_support:" << sdsl::size_in_bytes(inv_pi)
                  << std::endl;

        std::cout << "-----------" << std::endl;

        std::cout << "bvt:" << sdsl::size_in_bytes(bvt) << std::endl;
        std::cout << "bvt_rank1:" << sdsl::size_in_bytes(bvt_rank1)
                  << std::endl;
        std::cout << "bvt_sel0:" << sdsl::size_in_bytes(bvt_sel0) << std::endl;

        std::cout << "-----------" << std::endl;

        std::cout << "wt:" << sdsl::size_in_bytes(wtnt) << std::endl;

        std::cout << "-----------" << std::endl;
        std::cout << "str:" << vt.size() << std::endl;

        std::cout << "-----------" << std::endl;
        std::cout << "Grid:" << std::endl;
        _grid.print_size_in_bytes();

        std::cout << "-----------" << std::endl;
        std::cout << "Tree:" << std::endl;
        dfuds_tree.print_size_in_bytes();
    }

    virtual void print_size_in_bytes(std::fstream& f) const {

        sdsl::int_vector<> T(vt.size(),0);
        std::copy(vt.begin(),vt.end(),T.begin());
        sdsl::util::bit_compress(T);

        f << "Printing size in bytes, BYTES,MB\n" ;

        f << "bv-f-occ,"        <<std::to_string(  sdsl::size_in_bytes(bvfocc)          ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvfocc)        )/1000000.0  ) << std::endl;
        f << "bv-f-occ-rank-1," <<std::to_string(  sdsl::size_in_bytes(bvfocc_rank1)    ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvfocc_rank1)  )/1000000.0  ) << std::endl;
        f << "bv-f-occ-sel-0,"  <<std::to_string(  sdsl::size_in_bytes(bvfocc_sel0)     ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvfocc_sel0)   )/1000000.0  ) << std::endl;
        f << "bv-f-occ-sel-1,"  <<std::to_string(  sdsl::size_in_bytes(bvfocc_sel1)     ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvfocc_sel1)   )/1000000.0  ) << std::endl;
        f << "wt-s-occ,"        <<std::to_string(  sdsl::size_in_bytes(wtnt)            ) << ","<<std::to_string(  double(sdsl::size_in_bytes(wtnt)          )/1000000.0  ) << std::endl;
        f << "vt-str,"          <<std::to_string(  vt.size()                            ) << ","<<std::to_string(  double(vt.size()                          )/1000000.0  ) << std::endl;
        f << "vt-int,"          <<std::to_string(  sdsl::size_in_bytes(T)               ) << ","<<std::to_string(  double(sdsl::size_in_bytes(T)             )/1000000.0  ) << std::endl;
        f << "bvt,"             <<std::to_string(  sdsl::size_in_bytes(bvt)             ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvt)           )/1000000.0  ) << std::endl;
        f << "bvt-rank-1,"      <<std::to_string(  sdsl::size_in_bytes(bvt_rank1)       ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvt_rank1)     )/1000000.0  ) << std::endl;
        f << "bvt-sel-0,"       <<std::to_string(  sdsl::size_in_bytes(bvt_sel0)        ) << ","<<std::to_string(  double(sdsl::size_in_bytes(bvt_sel0)      )/1000000.0  ) << std::endl;
        f << "bvt-sel-0,"       <<std::to_string(  sdsl::size_in_bytes(pi)              ) << ","<<std::to_string(  double(sdsl::size_in_bytes(pi)            )/1000000.0  ) << std::endl;
        f << "bvt-sel-0,"       <<std::to_string(  sdsl::size_in_bytes(inv_pi)          ) << ","<<std::to_string(  double(sdsl::size_in_bytes(inv_pi)        )/1000000.0  ) << std::endl;
        f << "grid,"            <<std::to_string(  _grid.size_in_bytes()                ) << ","<<std::to_string(  double(_grid.size_in_bytes()              )/1000000.0  ) << std::endl;
        f << "grid-sb,"         <<std::to_string(  sdsl::size_in_bytes(_grid.sb)        ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.sb)      )/1000000.0  ) << std::endl;
        f << "grid-sl,"         <<std::to_string(  sdsl::size_in_bytes(_grid.sl)        ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.sl)      )/1000000.0  ) << std::endl;
        f << "grid-xb,"         <<std::to_string(  sdsl::size_in_bytes(_grid.xb)        ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.xb)      )/1000000.0  ) << std::endl;
        f << "grid-xb_rank1,"   <<std::to_string(  sdsl::size_in_bytes(_grid.xb_rank1)  ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.xb_rank1))/1000000.0  ) << std::endl;
        f << "grid-xb_sel1,"    <<std::to_string(  sdsl::size_in_bytes(_grid.xb_sel1)   ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.xb_sel1) )/1000000.0  ) << std::endl;
        f << "grid-xb_sel0,"    <<std::to_string(  sdsl::size_in_bytes(_grid.xb_sel0)   ) << ","<<std::to_string(  double(sdsl::size_in_bytes(_grid.xb_sel0) )/1000000.0  ) << std::endl;
        f << "tree,"            <<std::to_string(  dfuds_tree.size_in_bytes()           ) << ","<<std::to_string(  double(dfuds_tree.size_in_bytes()         )/1000000.0  ) << std::endl;



        std::cout << "Mapping Structures:" << std::endl;
        std::cout << "bvfocc:" << sdsl::size_in_bytes(bvfocc) << std::endl;
        std::cout << "bvfocc_rank1:" << sdsl::size_in_bytes(bvfocc_rank1)
                  << std::endl;
        std::cout << "bvfocc_sel0:" << sdsl::size_in_bytes(bvfocc_sel0)
                  << std::endl;
        std::cout << "bvfocc_sel1:" << sdsl::size_in_bytes(bvfocc_sel1)
                  << std::endl;

        std::cout << "-----------" << std::endl;

        std::cout << "pi:" << sdsl::size_in_bytes(pi) << std::endl;
        auto t_pi = pi;
        sdsl::util::bit_compress(t_pi);
        std::cout << "pi(bit_compress):" << sdsl::size_in_bytes(t_pi)
                  << std::endl;
        std::cout << "inv_pi_support:" << sdsl::size_in_bytes(inv_pi)
                  << std::endl;




        std::cout << "-----------" << std::endl;

        std::cout << "bvt:" << sdsl::size_in_bytes(bvt) << std::endl;
        std::cout << "bvt_rank1:" << sdsl::size_in_bytes(bvt_rank1)
                  << std::endl;
        std::cout << "bvt_sel0:" << sdsl::size_in_bytes(bvt_sel0) << std::endl;

        std::cout << "-----------" << std::endl;

        std::cout << "wt:" << sdsl::size_in_bytes(wtnt) << std::endl;

        std::cout << "-----------" << std::endl;
        std::cout << "str:" << vt.size() << std::endl;
        std::cout << "str(int-vector):" << sdsl::size_in_bytes(T) << std::endl;

        std::cout << "-----------" << std::endl;
        std::cout << "Grid:" << std::endl;
        _grid.print_size_in_bytes();

        std::cout << "-----------" << std::endl;
        std::cout << "Tree:" << std::endl;
        dfuds_tree.print_size_in_bytes();
        std::cout <<std::endl;
    }

    virtual void print() const {

        std::cout << "Printing mapping" << std::endl;
        std::cout << "bvfocc:";
        for (uint i = 0; i < bvfocc.size(); i++)
            std::cout << bvfocc[i] << " ";
        std::cout << "\npi:";
        for (int i = 0; i < pi.size(); i++)
            std::cout << pi[i] << " ";
        std::cout << "\nIpi:";
        for (int i = 0; i < pi.size(); i++)
            std::cout << inv_pi[i] << " ";
        std::cout << "\nbvt:";
        for (uint i = 0; i < bvt.size(); i++)
            std::cout << bvt[i] << " ";
        std::cout << "\nwt:";
        for (uint i = 0; i < wtnt.size(); i++)
            std::cout << wtnt[i] << " ";
        std::cout << "\ns:" << vt << std::endl;
        std::cout << "\nL:" ;
        for (uint i = 0; i < bvl.size(); i++)
            std::cout << bvl[i] << " ";
        std::cout << "\ndfuds:";
        dfuds_tree.print();
    }

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
    /** bitvector for characters nodes in seccondary nodes associated to 0 in
     * bvfocc*/
    t_maptbv bvt;
    typename t_maptbv::rank_1_type bvt_rank1;
    typename t_maptbv::select_0_type bvt_sel0;
    /** sequence for secondary occ nodes labels that aren't characters
     * associated to 0 in bvt*/
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


    void get_parent_ch_rank(const size_tree& node,size_tree& padre,size_tree & ch_r)const{

        size_tree bp_p = dfuds_tree.find_open(node-1);

        ch_r = (size_tree)(dfuds_tree.succ0(bp_p) - bp_p);

        padre = dfuds_tree.pred0(bp_p)+1;
    }
    /**Compare a rule suffix with an string left partition*/
    virtual int cmp_suffix_rule(const rule_type &, const std::string &,
                                len_type &j) const;
    virtual int cmp_suffix_rule_pre(const size_tree &, const std::string &,
                                    len_type j) const;
    /**Compare a rule suffix with an string left partition*/
    virtual int cmp_prefix_rule(const rule_type &, const std::string &,
                                len_type &) const;

    /**Compare a grammar suffix with an string right partition*/
    virtual int cmp_suffix_grammar(size_tree, const std::string &,
                                   len_type &j) const;

    //            virtual void dfs_expand_rule(const rule_type& X, std::string&
    //            s, len_type& off, const len_type& i, const len_type& m, bool
    //            suffix) const;

    virtual bool expand_prefix(const size_tree &node, std::string &s,
                               const len_type &m, const size_tree &off) const;
    virtual bool expand_prefixII(const size_tree &node, std::string &s,
                               const len_type &m,  long &off) const;


    virtual bool expand_suffix_gimp(const size_tree &node, std::string &s,
                                    const size_t &m, size_t &off) const;

    virtual bool expand_prefix_gimp(const size_tree &node, std::string &s,
                                    const size_t &m, size_t &off) const;
    virtual bool expand_all_rule(const size_tree &node, std::string &s, size_t &off) const;

    /**Implement an strategy to find the range of valid points of the grid*/
    //            virtual void find_range(const std::string & s,
    //            std::vector<len_type>& ) const  = 0;
    /**Algorithm to recursively track secondary occ in the grammar*/
    virtual void find_secondary_occ(const _node &,
                                    std::vector<len_type> &) const;

    /** given a rule return the node of its first occ */
    size_tree map_rule(const rule_type &) const;
    /** given a node return the label rule*/
    rule_type map_node(const size_tree &) const;

    /** compute the offset of a node in the text*/
    len_type offset_node(const size_tree &) const;
    /** true if the  node is a terminal node X -> "ACGTTT" expand just
     * characters nodes*/
    bool is_t_node(const size_tree &) const;
    /** true if the node is a character node*/
    bool is_ch_node(const size_tree &) const;
    bool is_ch_node(const size_tree &, size_tree &) const;
    /** take a terminal node and append the block of it and its siblings to the
     * string*/
    void append_block(const size_tree &, const size_tree &,
                      std::string &) const;
    /** take a terminal node and append the block(i ... j) to the
     * string str[off..off+j-i]*/

    void append_block(const size_t& off,const size_tree &, const size_tree &,
                      std::string &) const;

    void append_suffix_block(const size_t& off,const size_tree &, const size_tree &,
                      std::string &) const;
    /** take a character node and compare the block of it and its siblings to
     * the string*/
    int cmp_block(len_type pos, uint l, const std::string &str, len_type &i,
                  const bool &dir) const;

    bool eval_leaf( size_t&llb, size_t&pos,std::string &s,const size_tree& node ,const long &len)const{
        size_tree r = 0;
        if(is_ch_node(node,r)){
            size_t off_str = bvt_rank1(r)-1;

            if(len < 0){

                size_tree parent = 0,ch_r = 0;
                get_parent_ch_rank(node,parent,ch_r);
                auto ch = dfuds_tree.succ0(parent) - parent;
                size_tree cp = ch-ch_r+1;
                cp = cp > (s.size()- pos)?(s.size()- pos):cp;
//
//                if(cp < 0){
//                    std::cout<<"(ch-ch_r+1) < 0\n";
//                    sleep(10);
//                }
//                if(off_str + cp> vt.size() || s.size() < pos + cp){
//                    std::cout<<"OUT OF RANGE\n";
//                    std::cout<<"pos:"<<pos<<std::endl;
//                    std::cout<<"(ch-ch_r+1):"<<(ch-ch_r+1)<<std::endl;
//                    std::cout<<"ch:"<<ch<<std::endl;
//                    std::cout<<"ch_r:"<<ch_r<<std::endl;
//                    std::cout<<"len:"<<len<<std::endl;
//
////                    sleep(10);
//                }
                std::copy(vt.begin()+off_str,vt.begin()+(off_str+ cp),s.begin()+pos);
                pos += cp;
                llb += cp-1;
                return true;
            }
            std::copy(vt.begin()+ off_str,vt.begin()+off_str+len,s.begin()+pos);
            pos += len;
            llb += len-1;
            return true;

        }
        return false;

    }

    void expand_interval(const size_t &node, const std::pair<size_t, size_t> &range, std::string &s, std::size_t &pos)const
    {
        /*precondition range always belongs to [node.leftmost-leaf, node.rigthmost-leaf]*/

        /*base case*/

        /*recursion case*/



        auto llb = bvl_rank1(range.first+1); //left leaf begin
        auto pllb = bvl_select1(llb); // position of left leaf begin

        auto lle = bvl_rank1(range.second+1); //left leaf begin
        auto plle = bvl_select1(lle); // position of left leaf begin

//        auto f_eval_leaf = [this,&llb,&pos,&s](const size_tree& node ,const long &len){
//
//            size_tree r = 0;
//            if(is_ch_node(node,r)){
//                size_tree parent = 0,ch_r = 0;
//                get_parent_ch_rank(node,parent,ch_r);
//                size_tree ch = dfuds_tree.succ0(parent)-parent;
//                size_t off_str = bvt_rank1(r)-1;
//                if(len < 0){
//
//                    std::copy(vt.begin()+off_str,vt.begin()+(off_str+(ch-ch_r+1)),s.begin()+pos);
//                    pos += ch-ch_r+1;
//                    llb += ch-ch_r;
//                    return true;
//                }
//                std::copy(vt.begin()+ off_str,vt.begin()+off_str+len,s.begin()+pos);
//                pos += len;
//                llb += len-1;
//                return true;
//
//            }
//            return false;
//        };

        if(llb == lle){
            size_tree leaf_node = dfuds_tree.leaf_select(llb+256);

            if(!eval_leaf(llb, pos,s,leaf_node,s.size()-pos)){

                /***jump to first occ***/
                size_tree pre = dfuds_tree.pre_order(leaf_node);
                rule_type X = map_node(pre);
                size_tree preorder_focc = map_rule(X);
                size_tree _node = dfuds_tree.pre_order_select(preorder_focc);
                size_tree node_lr = dfuds_tree.leaf_rank(_node);
                size_t p_node_lr = bvl_select1(node_lr-256);
                auto _off1 = range.first-pllb;
                auto _off2 = range.second-pllb;
                auto l =  p_node_lr+_off1;
                auto r =  p_node_lr+_off2;
                return expand_interval(_node,std::make_pair(l,r),s,pos);

            }
        }




        /* left range*/
        {
            /* Need to jump to llb definition */
            size_tree left_leaf_node = dfuds_tree.leaf_select(llb+256);
            if(!eval_leaf(llb, pos,s,left_leaf_node,-1)){

                expand_suffix_gimp(left_leaf_node,s,bvl_select1(llb + 1)-range.first,pos);

            }

            if(pos  == s.size()) return;

        }

        /* middle range*/
        while(pos < s.size() && ++llb < lle){
            /*expand all the phrases in the middle*/
            size_tree middle_leaf_node = dfuds_tree.leaf_select(llb+256);
            if(!eval_leaf(llb, pos,s,middle_leaf_node,-1)){
                expand_all_rule(middle_leaf_node,s,pos);
                if(pos  >= s.size()) return;
            }
        }

        if(pos  == s.size()) return;
        /* right range*/
        if(pos < s.size())
        {
            /* Need to jump to llb definition */
            /* compute label X of lle */
            size_tree right_leaf_node = dfuds_tree.leaf_select(lle+256);
            if(!eval_leaf(llb, pos,s,right_leaf_node,s.size()-pos)){
                expand_prefix_gimp(right_leaf_node,s,s.size()-pos,pos);
            }
        }

    }
};
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::gcis_index(
    const gcis_index &) {}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::~gcis_index() {

}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::serialize(
    std::ofstream &f) const {

    _grid.serialize(f);

    sdsl::serialize(bvfocc, f);
    sdsl::serialize(bvfocc_sel0, f);
    sdsl::serialize(bvfocc_rank1, f);
    sdsl::serialize(bvfocc_sel1, f);
    sdsl::serialize(pi, f);
    sdsl::serialize(inv_pi, f);
    sdsl::serialize(bvt, f);
    sdsl::serialize(bvt_rank1, f);
    sdsl::serialize(bvt_sel0, f);
    sdsl::serialize(wtnt, f);

    size_t sz = vt.size();
    f.write((const char *)&sz, sizeof(sz));
    f.write((const char *)vt.data(), sizeof(char) * vt.size());

    sdsl::serialize(bvl, f);
    sdsl::serialize(bvl_select1, f);
    sdsl::serialize(bvl_rank1, f);
    dfuds_tree.serialize(f);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::load(
    std::ifstream &f) {

    _grid.load(f);
    sdsl::load(bvfocc, f);
    sdsl::load(bvfocc_sel0, f);
    sdsl::load(bvfocc_rank1, f);
    sdsl::load(bvfocc_sel1, f);
    sdsl::load(pi, f);
    sdsl::load(inv_pi, f);
    sdsl::load(bvt, f);
    sdsl::load(bvt_rank1, f);
    sdsl::load(bvt_sel0, f);
    sdsl::load(wtnt, f);

    size_t sz;
    f.read((char *)&sz, sizeof(sz));
    vt.resize(sz);
    f.read((char *)vt.data(), sizeof(char) * vt.size());

    sdsl::load(bvl, f);
    sdsl::load(bvl_select1, f);
    sdsl::load(bvl_rank1, f);

    bvl_select1.set_vector(&bvl);
    bvl_rank1.set_vector(&bvl);

    bvt_rank1.set_vector(&bvt);
    bvt_sel0.set_vector(&bvt);

    inv_pi.set_vector(&pi);

    bvfocc_sel0.set_vector(&bvfocc);
    bvfocc_sel1.set_vector(&bvfocc);
    bvfocc_rank1.set_vector(&bvfocc);

    dfuds_tree.load(f);

    this->n_suffixes = this->_grid.n_cols();
    this->n_rules = this->pi.size();
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
size_t
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::size_in_bytes()
    const {
    return _grid.size_in_bytes() + sdsl::size_in_bytes(bvfocc) +
           sdsl::size_in_bytes(pi) + sdsl::size_in_bytes(bvt) +
           sdsl::size_in_bytes(wtnt) + dfuds_tree.size_in_bytes() +
           sizeof(n_rules) + sizeof(n_suffixes);
}
/***
 * @brief This function return the preorder of the node where the first occ of
 * rule appears,
 * @tparam t_mapfbv
 * @tparam t_maptbv
 * @tparam t_mapwt
 * @tparam t_gridbv
 * @tparam t_gridwt
 * @param X rule id
 * @return preorder of the node
 */
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<>::size_tree
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::map_rule(
    const gcis_index::rule_type &X) const {
    using std::cout;
    using std::endl;
    // cout << "X = " << X << endl;
    unsigned long p = inv_pi[X];
    //        std::cout<<pi.size()<<std::endl;
    //        for(int i  = 0; i < pi.size(); ++i)
    //            std::cout<<"elem["<<i<<"]:"<<pi[i]<<std::endl;
    auto pre = bvfocc_sel1.select(p + 1);
    return pre + 1;
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
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<>::rule_type
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::map_node(
    const gcis_index::size_tree &p) const {

    if (bvfocc[p - 1]) {
        // cout<< "primary node" << endl;
        auto r = bvfocc_rank1.rank(p);
        return pi[r - 1];
    }
    // cout<< "secondary node" << endl;
    // secondary node
    auto r0 = p - bvfocc_rank1(p);
    if (bvt[r0 - 1]) {
        // character symbol
        // cout<< "character symbol" << endl;

        auto sr1 = bvt_rank1.rank(r0);
        return (unsigned char)vt[sr1 - 1];
    }
    // non-terminal
    //  cout<< "non-terminal symbol" << endl;

    auto sr0 = r0 - bvt_rank1.rank(r0);
    //        for (int i = 0; i < wtnt.size() ; ++i) {
    //            std::cout<<wtnt[i]<<" ";
    //        }
    //        std::cout<<std::endl;
    return wtnt[sr0 - 1];
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

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
        typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::displayII(
        const gcis_index::len_type &ii, const gcis_index::len_type &m,
        std::string &str) const {

        str.resize(m);
        long m_off = 0;
        extract_prefix_with_gap(dfuds_tree.root(),str,m_off,ii);


}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::display(
    const gcis_index::len_type &ii, const gcis_index::len_type &m,
    std::string &str) const {

    len_type i = ii + 256;
    /**found leaf containing the position rank [1...i+1)
     * the 256 terminals symbols are not included in the L bitvector
     * */
    len_type leaf_boundary = bvl_rank1.rank(ii + 1);
    size_tree leaf_pos = bvl_select1(leaf_boundary) + 256;
    len_type off = i - leaf_pos;
    /**found the associated leaf node */
    size_tree node = dfuds_tree.leaf_select(leaf_boundary + 256);
    // in case that is a single character
    std::string sblock;
    sblock.reserve(m);
    long strs = sblock.size();
    size_tree pos = 0, parent = 0, ch = 0, j = 0;
    // extraction process bottom up
    parent = dfuds_tree.parent(node);
    ch = dfuds_tree.children(parent);
    j = dfuds_tree.childrank(node);

    if (is_ch_node(node, pos)) {

        append_block(bvt_rank1(pos) - 1, ch - j + 1, sblock);
        strs = sblock.size();


        /**
         * go up in the tree while the current node will the last child
         * */
        do {
            node = parent;
            parent = dfuds_tree.parent(node);
            ch = dfuds_tree.children(parent);
            j = dfuds_tree.childrank(node);

        } while (j == ch);

        node = dfuds_tree.nsibling(node);
        ++j;


    }

    while (strs - off < m) {

        while (j <= ch && expand_prefix(node, sblock, m, off)) {
            node = dfuds_tree.nsibling(node);
            ++j;
        }
        if (j > ch) {
            node = parent;
            parent = dfuds_tree.parent(node);
            ch = dfuds_tree.children(parent);
            j = dfuds_tree.childrank(node) + 1;
            if (j <= ch)
                node = dfuds_tree.nsibling(node);
        } else {
            strs = sblock.size();
        }
    }
    //    if(dfuds_tree.is_leaf(node)){
    //        expand_prefix(node, sblock, m, off);
    //    }

    /** remove innecesary characters */
    str.resize(m);
    std::copy(sblock.begin() + off, sblock.begin() + off + m, str.begin());
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<>::len_type
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::offset_node(
    const gcis_index::size_tree &node) const {

    if (node == 3)
        return 0;
//todo arreglar problemas con las hojas
    long long lrank = dfuds_tree.leaf_rank(node);
    return bvl_select1.select(lrank - 256);
}



template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
        typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::extract_prefix_with_gap(
        const size_tree& node,
        std::string& str,
        long &  m_off,
        const uint & gap) const
{



    if (!dfuds_tree.is_leaf(node)) {

        size_tree leaf_r = 0;

        if(node == 3)
            leaf_r = 257;
        else
            leaf_r = dfuds_tree.leaf_rank(node);

//        auto leafnode = dfuds_tree.leaf_select( leaf_r );
        auto node_off = bvl_select1(leaf_r-256);
        auto pos = node_off+gap;
        auto c_leaf = bvl_rank1( pos + 1 );
        size_tree nnode = dfuds_tree.leaf_select( c_leaf + 256);

        size_tree r = 0;
        size_tree parent = 0;
//        size_tree pppp = dfuds_tree.parent(nnode);
//        size_tree chhh = dfuds_tree.childrank(nnode);

        if(is_ch_node(nnode, r)) {
            // si es un nodo caracter

            size_tree j = 0;
            get_parent_ch_rank(nnode,parent,j);

            size_tree l = dfuds_tree.children(parent);

            auto rest = str.size() - m_off;
            auto rest_node = l - j+1;
            auto min = rest_node < rest ? rest_node :rest;
            size_t str_off = bvt_rank1(r)-1;
            append_block(m_off,str_off, min, str);
            m_off += min;

            if (m_off >= str.size()) return false;

            nnode = parent;

        }else{
            // si es un nodo segunda occ
            //calculamos
            size_t  bvls = bvl_select1(c_leaf);
            auto ngap = (node_off+gap) - bvls;

            if(ngap == 0){
                if(!expand_prefixII(nnode,str,str.size(),m_off))
                    return false;
            }else{
                if(!extract_prefix_with_gap(nnode,str,m_off,ngap))
                    return false;
            }

        }

        // seguir extrayendo recursivamente al nodo inicial
        /**
        * go up in the tree while the current node will the last child
        * */

        /**
         * para saber si es ancestro
         *
         * */



        size_tree fwdExcessNode = dfuds_tree.fwd_excess(node-1,-1);
        bool expandir = true;
        nnode = dfuds_tree.nsibling(nnode);
        /**
         * aun estamos en el subarbol de node
         * */

        while( m_off < str.size() && node < nnode && nnode <= fwdExcessNode && expandir){
            expandir = expand_prefixII(nnode,str,str.size(),m_off);
            nnode = dfuds_tree.nsibling(nnode);
        }

        if(m_off >= str.size()) return false;

        return true;

    }

    // leaf case !!!!! WE DO NOT TREAT CHARACTER LEAVES
    // as we never arrive to character leaves the leaf must to be a non-terminal
    assert(!is_ch_node(node));
    // leaf we need to jump to its first occ
    size_tree pre = dfuds_tree.pre_order(node);
    rule_type X = map_node(pre);
    size_tree preorder_focc = map_rule(X);

    return extract_prefix_with_gap(dfuds_tree.pre_order_select(preorder_focc), str, m_off,gap);

}


template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
        typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_prefixII(

        const gcis_index::size_tree &node,
        std::string &s,//is string to store
        const gcis_index::len_type & m, //is the number of character to extract
        long & off// current position in str to copy,
        ) const {
        /** if the node is not a leaf */
        if (!dfuds_tree.is_leaf(node)) {
            /** check if its first child is a character node */
            size_tree fchild = dfuds_tree.fchild(node);
            size_tree l = dfuds_tree.children(node);
            // if its a leaf and a character node then all children are character
            // and we extract the block
            size_tree r = 0;

            if (dfuds_tree.is_leaf(fchild) && is_ch_node(fchild, r)) {
                // append the corresponding block of characters
                //            size_tree pos = bvt_rank1.rank(r);
                auto min = l < (m-off)?l:(m-off);
                append_block(off,bvt_rank1(r)- 1, min, s);
                off += min;
                return (off < m ); // continue expanding
            }
            // if it is no a character node -> expand in dfs
            uint i = 1;
            while (i <= l) {
                if (!expand_prefixII(fchild, s, m, off))
                    return false;
                fchild = dfuds_tree.nsibling(fchild);
                ++i;
            }
            return true;
        }

        // leaf case !!!!! WE DO NOT TREAT CHARACTER LEAVES
        // as we never arrive to character leaves the leaf must to be a non-terminal
        assert(!is_ch_node(node));
        // leaf we need to jump to its first occ
        size_tree pre = dfuds_tree.pre_order(node);
        rule_type X = map_node(pre);
        size_tree preorder_focc = map_rule(X);

        return expand_prefixII(dfuds_tree.pre_order_select(preorder_focc), s, m, off);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_prefix(
    const gcis_index::size_tree &node, std::string &s,
    const gcis_index::len_type &m, const size_tree &off) const {
    /** if the node is not a leaf */


    if (!dfuds_tree.is_leaf(node)) {
        /** check if its first child is a character node */
        size_tree fchild = dfuds_tree.fchild(node);
        size_tree l = dfuds_tree.children(node);
        // if its a leaf and a character node then all children are character
        // and we extract the block
        size_tree r = 0;
        if (dfuds_tree.is_leaf(fchild) && is_ch_node(fchild, r)) {
            // append the corresponding block of characters
            //            size_tree pos = bvt_rank1.rank(r);
            append_block(bvt_rank1(r)- 1, l, s);
            return long(s.size() - off) < m; // continue expanding
        }
        // if it is no a character node -> expand in dfs
        uint i = 1;
        while (i <= l) {
            if (!expand_prefix(fchild, s, m, off))
                return false;
            fchild = dfuds_tree.ichild(node, ++i);
        }
        return true;
    }

    // leaf case !!!!! WE DO NOT TREAT CHARACTER LEAVES
    // as we never arrive to character leaves the leaf must to be a non-terminal
    assert(!is_ch_node(node));
    // leaf we need to jump to its first occ
    size_tree pre = dfuds_tree.pre_order(node);
    rule_type X = map_node(pre);
    size_tree preorder_focc = map_rule(X);

    return expand_prefix(dfuds_tree.pre_order_select(preorder_focc), s, m, off);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv,
               t_gridwt>::cmp_suffix_rule_pre(const size_tree &pre,
                                              const std::string &str,
                                              len_type j) const {
    rule_type X = map_node(pre);
    size_tree node = dfuds_tree.pre_order_select(pre);
    // case leaf
    if (dfuds_tree.is_leaf(node)) // if it is a leaf then is a symbol
    {
        if ((unsigned char)str[j] == X)
            return 0;
        if ((unsigned char)str[j] > X)
            return 1;
        j--;
        return -1;
    }
    // lambda function for recursively cmp
    std::function<int(const size_tree &, const size_tree &)> dfs_cmp_suffix;

    dfs_cmp_suffix = [this, &dfs_cmp_suffix, &str,
                      &j](const size_tree &node, const size_tree &pre_node) {
        // !!!!!!WE NEVER REACH A CHARACTER NODE
        // case leaf
        if (dfuds_tree.is_leaf(
                node)) // as we never reach a character node if we reach a leaf
                       // it has to be non-terminal leaf
        {              // jump to it definition
            rule_type Y = map_node(pre_node);
            size_tree pre = map_rule(Y);
            size_tree def_node = dfuds_tree.pre_order_select(pre);
            int r = dfs_cmp_suffix(def_node, pre);
            return r;
        }

        // check if it is a terminal node X->"AJSKLDJALKSJD"

        size_tree ch = dfuds_tree.children(node);
        size_tree lchild = dfuds_tree.ichild(node, ch);
        size_tree pre_lchild = dfuds_tree.pre_order(lchild);

        if (dfuds_tree.is_leaf(lchild)) {
            auto r0 = pre_lchild - bvfocc_rank1(pre_lchild);
            if (bvt[r0 - 1]) { // is a terminal node (just compare characters)
                               // and return
                auto r1 = bvt_rank1.rank(r0);
                return cmp_block(r1 - 1, (uint)ch, str, j, true);
            }
        }
        //

        // compare first child, is already precomputed
        int r = dfs_cmp_suffix(lchild, pre_lchild);
        if (r != 0)
            return r; // return if missmatch

        uint i = ch - 1;
        // compare the rest of children
        while (i > 0) {
            lchild = dfuds_tree.ichild(node, i);
            pre_lchild = dfuds_tree.pre_order(lchild);
            r = dfs_cmp_suffix(lchild, pre_lchild);
            if (r != 0)
                return r;
            --i;
        }
        return 0;
    };

    return dfs_cmp_suffix(node, pre);
}
/**
 * @brief compare a suffix of the rule with the left part of pattern partition
 * in reverse order cmp(X,str,j) rev(X) ? rev(str[1...j])
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
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv,
               t_gridwt>::cmp_suffix_rule(const gcis_index::rule_type &X,
                                          const std::string &str,
                                          len_type &j) const {

    if (X < 256) {
        unsigned char m = (unsigned char)X;
        if ((unsigned char)str[j] < X) {
            if ((unsigned char)str[j] >= m)
                std::cout << "UNSIGNED CHAR PROBLEM\n";
            return -1;
        }

        if ((unsigned char)str[j] > X) {
            if ((unsigned char)str[j] <= m)
                std::cout << "UNSIGNED CHAR PROBLEM\n";
            return 1;
        }

        if ((unsigned char)str[j] != m)
            std::cout << "UNSIGNED CHAR PROBLEM\n";

        --j;
        return 0;
    }
    // compute the preorder node and node pos
    size_tree pre = map_rule(X);
    size_tree node = dfuds_tree.pre_order_select(pre);

    std::function<int(const size_tree &, const size_tree &)> dfs_cmp_suffix;

    dfs_cmp_suffix = [this, &dfs_cmp_suffix, &str,
                      &j](const size_tree &node, const size_tree &pre_node) {
        // !!!!!!WE NEVER REACH A CHARACTER NODE
        // case leaf
        if (dfuds_tree.is_leaf(
                node)) // as we never reach a character node if we reach a leaf
                       // it has to be non-terminal leaf
        {              // jump to it definition
            rule_type Y = map_node(pre_node);
            size_tree pre = map_rule(Y);
            size_tree def_node = dfuds_tree.pre_order_select(pre);
            int r = dfs_cmp_suffix(def_node, pre);
            return r;
        }

        // check if it is a terminal node X->"AJSKLDJALKSJD"

        size_tree ch = dfuds_tree.children(node);
        size_tree lchild = dfuds_tree.ichild(node, ch);
        size_tree pre_lchild = dfuds_tree.pre_order(lchild);

        if (dfuds_tree.is_leaf(lchild)) {
            auto r0 = pre_lchild - bvfocc_rank1(pre_lchild);
            if (bvt[r0 - 1]) { // is a terminal node (just compare characters)
                               // and return
                auto r1 = bvt_rank1.rank(r0);
                return cmp_block(r1 - 1, (uint)ch, str, j, true);
            }
        }
        //

        // compare first child, is already precomputed
        int r = dfs_cmp_suffix(lchild, pre_lchild);
        if (r != 0)
            return r; // return if missmatch

        size_tree i = ch - 1;
        // compare the rest of children
        while (j >= 0 && i > 0) {
            lchild = dfuds_tree.ichild(node, i);
            pre_lchild = dfuds_tree.pre_order(lchild);
            r = dfs_cmp_suffix(lchild, pre_lchild);
            if (r != 0)
                return r;
            --i;
        }
        return 0;
    };

    return dfs_cmp_suffix(node, pre);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_t_node(
    const gcis_index::size_tree &node) const {
    // THE NODE MUST NOT BE A LEAF
    assert(!dfuds_tree.is_leaf(node));
    size_tree fchild = dfuds_tree.fchild(node);
    return is_ch_node(fchild);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_ch_node(
    const gcis_index::size_tree &node) const {
    assert(dfuds_tree.is_leaf(node));
    size_tree pre = dfuds_tree.pre_order(node);
    auto r0 = pre - bvfocc_rank1(pre);
    return bvt[r0 - 1];
}
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
bool gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::is_ch_node(
    const gcis_index::size_tree &node, size_tree &idx) const {

    assert(dfuds_tree.is_leaf(node));

    size_tree pre = dfuds_tree.pre_order(node);

    idx = pre - bvfocc_rank1(pre);
//    indx_str = bvt_rank1(idx);
    return bvt[idx - 1];
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
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::cmp_block(
    gcis_index::len_type pos, uint l, const std::string &str, len_type &i,
    const bool &dir) const {

    if (dir) {

        //    std::string sss; sss.resize(l);
        //    std::copy(vt.begin()+pos-l , vt.begin()+pos,sss.begin());
        //    std::cout<<"comparing  rule"<<std::endl;
        //    std::cout<<"cmp_block("<<sss<<","<< "str["<<i<<"]
        //    "<<str<<")"<<std::endl;

        assert(l <= pos + 1);
        while (l && i >= 0) {
            if ((unsigned char)str[i] < (unsigned char)vt[pos])
                return -1;
            if ((unsigned char)str[i] > (unsigned char)vt[pos])
                return 1;
            --pos;
            --i;
            --l;
        }

        // if( i < 0 && l != 0 ) return -1;

        return 0;
    }

    assert(pos + l <= vt.size());

    //    std::string sss; sss.resize(l);
    //    std::copy(vt.begin()+pos , vt.begin()+pos+l,sss.begin());
    //    std::cout<<"comparing  suffix"<<std::endl;
    //    std::cout<<"cmp_block("<<sss<<","<< "str["<<i<<"]
    //    "<<str<<")"<<std::endl;

    uint ssz = str.size();
    while (l && i < ssz) {
        if ((unsigned char)str[i] < (unsigned char)vt[pos])
            return -1;
        if ((unsigned char)str[i] > (unsigned char)vt[pos])
            return 1;
        ++pos;
        ++i;
        --l;
    }
    // if( i == ssz && l != 0 ) return -1;

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
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv,
               t_gridwt>::cmp_prefix_rule(const gcis_index::rule_type &X,
                                          const std::string &str,
                                          len_type &j) const {

    if (X < 256) {
        unsigned char m = (unsigned char)X;
        if ((unsigned char)str[j] < X) {
            return -1;
        }

        if ((unsigned char)str[j] > X) {
            return 1;
        }

        j++;
        return 0;
    }

    // compute the preorder node and node pos
    size_tree pre = map_rule(X);
    size_tree node = dfuds_tree.pre_order_select(pre);
    // lambda function for recursively cmp
    std::function<int(const size_tree &, const size_tree &)> dfs_cmp_prefix;

    dfs_cmp_prefix = [this, &dfs_cmp_prefix, &str,
                      &j](const size_tree &node,
                          const size_tree &pre_node) -> int {
        // !!!!!!WE NEVER REACH A CHARACTER NODE
        // case leaf
        if (dfuds_tree.is_leaf(
                node)) // as we never reach a character node if we reach a leaf
                       // it has to be non-terminal leaf
        {              // jump to it definition
            rule_type Y = map_node(pre_node);
            size_tree pre = map_rule(Y);
            size_tree def_node = dfuds_tree.pre_order_select(pre);
            uint r = dfs_cmp_prefix(def_node, pre);
            return r;
        }

        // check if it is a terminal node X->"AJSKLDJALKSJD"

        size_tree ch = dfuds_tree.children(node);
        size_tree fchild = dfuds_tree.fchild(node);
        size_tree pre_fchild = pre_node + 1;
        // ask if fchild is a leaf
        if (dfuds_tree.is_leaf(fchild)) {
            auto r0 = pre_fchild - bvfocc_rank1(pre_fchild);
            if (bvt[r0 - 1]) { // is a terminal node (just compare characters)
                               // and return
                auto r1 = bvt_rank1.rank(r0);
                int r = cmp_block(r1 - 1, (uint)ch, str, j, false);
                return r;
            }
            //
        }

        // compare first child is already precomputed
        uint r = dfs_cmp_prefix(fchild, pre_fchild);
        if (r != 0)
            return r; // return if missmatch

        uint i = 2;
        // compare the rest of children
        while (j < str.size() && i <= ch) {
            fchild = dfuds_tree.ichild(node, i);
            pre_fchild = dfuds_tree.pre_order(fchild);
            r = dfs_cmp_prefix(fchild, pre_fchild);
            if (r != 0)
                return r;
            ++i;
        }
        return 0;
    };

    return dfs_cmp_prefix(node, pre);
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
template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
int gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv,
               t_gridwt>::cmp_suffix_grammar(gcis_index::size_tree pre_node,
                                             const std::string &str,
                                             len_type &j) const {

    size_tree node = dfuds_tree.pre_order_select(pre_node);
    size_tree parent = dfuds_tree.parent(node);
    size_tree ch = dfuds_tree.children(parent);
    size_tree rch = dfuds_tree.childrank(node);
    rule_type Y = map_node(pre_node);

    // cout << "Y = " << Y << endl;
    // cout << "pre_node = " << pre_node << endl;

#ifdef TOYEXAMPLE

    switch (Y) {
    case 65:
        Y = 0;
        break;
    case 67:
        Y = 1;
        break;
    case 71:
        Y = 10;
        break;
    case 84:
        Y = 13;
        break;
    }
#endif
    int r = cmp_prefix_rule(Y, str, j);
    //        if(j == str.size()) return 0;
    if (r != 0)
        return r;
    //        if(r == 0 && j < str.size()) return 1;

    for (size_tree i = rch + 1; j < str.size() && i <= ch; ++i) {

        node = dfuds_tree.ichild(parent, i);
        pre_node = dfuds_tree.pre_order(node);
        Y = map_node(pre_node);

        // cout << "Y2 = " << Y << endl;
        // cout << "pre_node2 = " << pre_node << endl;

#ifdef TOYEXAMPLE

        switch (Y) {
        case 65:
            Y = 0;
            break;
        case 67:
            Y = 1;
            break;
        case 71:
            Y = 10;
            break;
        case 84:
            Y = 13;
            break;
        }
#endif
        r = cmp_prefix_rule(Y, str, j);
        //            if(j == str.size()) return 0;
        if (r != 0)
            return r;
        //            if(r == 0 && j < str.size()) return 1;
    }

    if (r == 0 && j < str.size())
        return 1;

    return 0;
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv,
                t_gridwt>::find_secondary_occ(const _node &v,
                                              std::vector<len_type> &occ)
    const {

    struct secc_node {
        size_tree pre;
        size_tree pos;
        rule_type X;
        long long off;
        secc_node(const size_tree &_pre, const size_tree &_pos,
                  const rule_type &_x, const len_type &_off)
            : pre(_pre), pos(_pos), X(_x), off(_off) {}
    };

    std::deque<secc_node> S;

    {
        size_tree node = dfuds_tree.pre_order_select(v.first);
        size_tree parent = dfuds_tree.parent(node);
        size_tree pre_parent = dfuds_tree.pre_order(parent);
        // this rule has to be a non-terminal so we search it in the golynski
        // sequence
        rule_type X = map_node(pre_parent);
        long long off = v.second + offset_node(node) - offset_node(parent);

        size_t nocc = this->wtnt.rank(wtnt.size(), X);
        secc_node ts(pre_parent, parent, X, off);
        S.emplace_back(ts);

        for (int i = 1; i <= nocc; ++i) {

            auto ps = wtnt.select(i, X) + 1;
            auto bvt0 = bvt_sel0.select(ps) + 1;
            auto bvfocc0 = bvfocc_sel0(bvt0);
            size_tree tnode = dfuds_tree.pre_order_select(bvfocc0 + 1);
            secc_node t(bvfocc0 + 1, tnode, X, off);
            S.emplace_back(t);
        }
    }

    while (!S.empty()) {
        if (S.front().pre == 1)
            occ.push_back((len_type)S.front().off);
        else {
            size_tree parent = dfuds_tree.parent(S.front().pos);
            size_tree pre = dfuds_tree.pre_order(parent);
            rule_type X = map_node(pre);
            long long off = S.front().off + offset_node(S.front().pos) -
                            offset_node(parent);
            secc_node ts(pre, parent, X, off);
            S.push_back(ts);
            size_tree nocc = this->wtnt.rank(wtnt.size(), X);
            for (int i = 1; i <= nocc; ++i) {

                auto ps = wtnt.select(i, X) + 1;
                auto bvt0 = bvt_sel0.select(ps) + 1;
                auto bvfocc0 = bvfocc_sel0(bvt0);
                size_tree tnode = dfuds_tree.pre_order_select(bvfocc0 + 1);
                secc_node t(bvfocc0 + 1, tnode, X, off);
                S.emplace_back(t);
            }
        }
        S.pop_front();
    }
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_tree(
    const sdsl::bit_vector &b) {
    dfuds_tree.set_tree(b);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_bvfocc(
    const t_mapfbv &b) {

    bvfocc = std::move(b);
    bvfocc_sel1 = typename t_mapfbv::select_1_type(&bvfocc);
    bvfocc_sel0 = typename t_mapfbv::select_0_type(&bvfocc);
    bvfocc_rank1 = typename t_mapfbv::rank_1_type(&bvfocc);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_bvt(
    const t_maptbv &b) {
    bvt = std::move(b);
    bvt_rank1 = typename t_maptbv::rank_1_type(&bvt);
    bvt_sel0 = typename t_maptbv::select_0_type(&bvt);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_vt(
    const std::string &str) {
    vt = str;
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_pi(
    const sdsl::int_vector<> &_pi) {
    pi = std::move(_pi);
    inv_pi = sdsl::inv_perm_support<>(&pi);
    n_rules = pi.size();
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::gcis_index() =
    default;

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_grid(
    const gcis_index::t_grid &g) {
    _grid = g;
    n_suffixes = _grid.n_cols();
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_nt(
    const t_mapwt &nt) {
    wtnt = std::move(nt);
    // sdsl::util::assign(wtnt, nt);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::append_block(
    const gcis_index::size_tree &pos, const gcis_index::size_tree &l,
    std::string &str) const {

    uint b = str.size();
    str.resize(str.size() + l);
    std::copy(vt.begin() + pos, vt.begin() + pos + l, str.begin() + b);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
        typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::append_block(
        const size_t & off, const gcis_index::size_tree &pos, const gcis_index::size_tree &l,
        std::string &str) const {
        std::copy(vt.begin() + pos, vt.begin() + pos + l, str.begin() + off);
    }

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_l(
    const sdsl::sd_vector<> &_l) {
    bvl = std::move(_l);
    bvl_rank1 = sdsl::sd_vector<>::rank_1_type(&bvl);
    bvl_select1 = sdsl::sd_vector<>::select_1_type(&bvl);
}

template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
          typename t_gridbv, typename t_gridwt>
void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::set_sigma(
    const uint &s) {
    sigma = s;
}

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_suffix_gimp(const gcis_index::size_tree &node,
                                                                                    std::string &s,
                                                                                    const size_t &m,
                                                                                    size_t &pos) const {

        size_t off = 0;
        std::deque<std::pair<size_tree ,std::pair<size_tree ,size_tree >>> Q;

        size_tree _node = node;
        if(dfuds_tree.is_leaf(_node)){
            size_tree pre = dfuds_tree.pre_order(node);
            rule_type X = map_node(pre);
            size_tree preorder_focc = map_rule(X);
            _node = dfuds_tree.pre_order_select(preorder_focc);
        }


        auto first_leaf = dfuds_tree.leaf_rank(_node);
        auto current_leaf= first_leaf+dfuds_tree.leaf_num(_node)-1;


        while(current_leaf >= first_leaf )
        {
            // current leaf is a Terminal rule?

            auto _c_node = dfuds_tree.leaf_select(current_leaf);
            size_tree r = 0;
            if(is_ch_node(_c_node,r))
            {
                auto parent = dfuds_tree.parent(_c_node);
                auto ch = dfuds_tree.children(parent);

                size_t off_str = bvt_rank1(r)-1;

                if( ch + off  > m ){
                    uint cp = m-off;
                    // extraemos solo los m-off restantes. retornamos false

                    std::copy(vt.begin()+(off_str-cp+1),vt.begin()+off_str+1,s.begin()+m-cp-off);
                    pos +=cp;

                    return false;

                }else{
                    //sino extraemos ch elements e incrementamos off en ch
                    // y decrementamos current leaf en ch

                    std::copy(vt.begin()+(off_str-ch+1),vt.begin()+off_str+1,s.begin()+m-ch-off);
                    off+=ch;
                    pos+=ch;
                    current_leaf-=ch;
                }

                if(off == m) return false;

                while(current_leaf < first_leaf && !Q.empty()){

                    current_leaf = Q.front().second.first-1;

                    first_leaf = Q.front().second.second;

                    Q.pop_front();
                }
            }
            else
            {

                //Save actua state;
                Q.emplace_front(std::make_pair(_node,std::make_pair(current_leaf,first_leaf)));

                // leaf we need to jump to its first occ
                size_tree pre = dfuds_tree.pre_order(_c_node);
                rule_type X = map_node(pre);
                size_tree preorder_focc = map_rule(X);
                size_tree _node = dfuds_tree.pre_order_select(preorder_focc);
                first_leaf= dfuds_tree.leaf_rank(_node);
                current_leaf = first_leaf+dfuds_tree.leaf_num(_node)-1;
            }
        }

        return true;

    }


    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_prefix_gimp(const gcis_index::size_tree &node,
                                                                                    std::string &s,
                                                                                    const size_t  &m,
                                                                                    size_t &pos) const {

        size_t off = 0;
        size_tree _node = node;
        if(dfuds_tree.is_leaf(_node)){
            size_tree pre = dfuds_tree.pre_order(node);
            rule_type X = map_node(pre);
            size_tree preorder_focc = map_rule(X);
            _node = dfuds_tree.pre_order_select(preorder_focc);
        }

        std::deque<std::pair<size_tree ,std::pair<size_tree ,size_tree > > > Q;

        auto current_leaf = dfuds_tree.leaf_rank(_node);
        auto last_leaf= current_leaf +dfuds_tree.leaf_num(_node)-1;


        while(current_leaf <= last_leaf )
        {
            // current leaf is a Terminal rule?

            auto _c_node = dfuds_tree.leaf_select(current_leaf);
            size_tree r = 0;
            if(is_ch_node(_c_node,r))
            {
                auto parent = dfuds_tree.parent(_c_node);
                auto ch = dfuds_tree.children(parent);
                size_t off_str = bvt_rank1(r)-1;

                if( ch + off  > m ){
                    // extraemos solo los m-off restantes. retornamos false

//                    append_block(off,r-1,m-off,s);
                    std::copy(vt.begin()+off_str,vt.begin()+off_str+m-off,s.begin()+pos);
                    pos+=m-off;
                    return false;

                }else{
                    //sino extraemos ch elements e incrementamos off en ch
                    // y decrementamos current leaf en ch
//                    append_block(off,r-1,ch,s);

                    std::copy(vt.begin()+off_str,vt.begin()+off_str+ch,s.begin()+pos);
                    off+=ch;
                    pos+=ch;
                    current_leaf += ch;
                }
                if(off == m) return false;
                while(current_leaf > last_leaf && !Q.empty()){

                    current_leaf = Q.front().second.first+1;

                    last_leaf  = Q.front().second.second;

                    Q.pop_front();
                }
            }
            else
            {

                //Save actua state;
                Q.emplace_front(std::make_pair(_node,std::make_pair(current_leaf,last_leaf)));

                // leaf we need to jump to its first occ
                size_tree pre = dfuds_tree.pre_order(_c_node);
                rule_type X = map_node(pre);
                size_tree preorder_focc = map_rule(X);
                size_tree _node = dfuds_tree.pre_order_select(preorder_focc);

                current_leaf = dfuds_tree.leaf_rank(_node);
                last_leaf= current_leaf +dfuds_tree.leaf_num(_node)-1;

            }
        }

        return true;

    }

    template<typename t_mapfbv, typename t_maptbv, typename t_mapwt, typename t_gridbv, typename t_gridwt>
    bool
    gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::expand_all_rule(const gcis_index::size_tree &node,
                                                                                    std::string &s,
                                                                                    size_t &pos) const {


        size_t off =0;

        std::deque<std::pair<size_tree ,std::pair<size_tree ,size_tree >>> Q;


        size_tree _node = node;
        if(dfuds_tree.is_leaf(_node)){
            size_tree pre = dfuds_tree.pre_order(node);
            rule_type X = map_node(pre);
            size_tree preorder_focc = map_rule(X);
            _node = dfuds_tree.pre_order_select(preorder_focc);
        }


        auto current_leaf = dfuds_tree.leaf_rank(_node);
        auto last_leaf= current_leaf +dfuds_tree.leaf_num(_node)-1;




        while(current_leaf <= last_leaf )
        {
            // current leaf is a Terminal rule?

            auto _c_node = dfuds_tree.leaf_select(current_leaf);
            size_tree r = 0;
            if(is_ch_node(_c_node,r))
            {
                auto parent = dfuds_tree.parent(_c_node);
                auto ch = dfuds_tree.children(parent);
                size_t off_str = bvt_rank1(r)-1;

                std::copy(vt.begin()+off_str,vt.begin()+off_str+ch,s.begin()+pos);
//                append_block(off,r-1,ch,s);
                off+=ch;
                pos+=ch;
                current_leaf+=ch;

                if(pos == s.size()) return false;


                while(current_leaf > last_leaf && !Q.empty()){

                    current_leaf = Q.front().second.first+1;

                    last_leaf  = Q.front().second.second;

                    Q.pop_front();
                }
            }
            else
            {

                //Save actua state;
                Q.emplace_front(std::make_pair(_node,std::make_pair(current_leaf,last_leaf)));

                // leaf we need to jump to its first occ
                size_tree pre = dfuds_tree.pre_order(_c_node);
                rule_type X = map_node(pre);
                size_tree preorder_focc = map_rule(X);
                _node = dfuds_tree.pre_order_select(preorder_focc);

                current_leaf = dfuds_tree.leaf_rank(_node);
                last_leaf= current_leaf +dfuds_tree.leaf_num(_node)-1;

            }
        }

        return true;

    }

    template <typename t_mapfbv, typename t_maptbv, typename t_mapwt,
            typename t_gridbv, typename t_gridwt>
    void gcis_index<t_mapfbv, t_maptbv, t_mapwt, t_gridbv, t_gridwt>::append_suffix_block(
            const size_t & off, const gcis_index::size_tree &pos, const gcis_index::size_tree &l,
            std::string &str) const {
        std::copy(vt.begin() + pos - l, vt.begin() + pos, str.begin() + off);
    }

} // namespace gcis_index_private

#endif // GCIS_GCIS_INDEX_HPP
