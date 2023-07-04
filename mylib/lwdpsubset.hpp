#ifndef SIMFRONTIER_LWDPSUBSET_HPP
#define SIMFRONTIER_LWDPSUBSET_HPP

#include "common.hpp"
#include "graph.hpp"
#include "../tdzdd/DdStructure.hpp"
#include "../tdzdd/dd/Node.hpp"

#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <utility>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <cstdio>

class StateSubset{
public:
  std::array<int8_t, 16> comp;
  std::array<int8_t, 16> vset;
  uint64_t bcode;
  int8_t cnum;
  
  StateSubset() {};
  StateSubset(uint64_t _bcode): bcode(_bcode) {};
          
  bool operator==(const StateSubset& rhs) const {
    return bcode == rhs.bcode && !memcmp(comp.data(), rhs.comp.data(), 16) && !memcmp(vset.data(), rhs.vset.data(), 16);
  }
};

namespace std{
  template <>
  struct hash<StateSubset>{
    public:
    uint64_t operator()(const StateSubset& s) const{
      const uint64_t *pc = reinterpret_cast<const uint64_t*>(s.comp.data());
      const uint64_t *pv = reinterpret_cast<const uint64_t*>(s.vset.data());
      return s.bcode * 17687306105512151009ULL + *pc * 9088778391368939145ULL + *(pc+1) * 98696038758792049ULL + *pv * 314159257ULL + *(pv+1);
    }
  };
}

template <typename VALTYPE>
class DPBlockSubset{
public:
  StateSubset s;
  VALTYPE p, r;
  std::vector<VALTYPE> q;
  int8_t cnum;
  int64_t lo;
  int64_t hi;
  std::vector<int8_t> vlo;
  std::vector<int8_t> vhi;
  
  DPBlockSubset() {};
  DPBlockSubset(const StateSubset& _s, int8_t _siz, VALTYPE ZEROELEM): s(_s), p(ZEROELEM), r(ZEROELEM), q(_siz, ZEROELEM), cnum(_siz), vlo(_siz), vhi(_siz) {};
  DPBlockSubset(StateSubset&& _s, int8_t _siz, VALTYPE ZEROELEM): s(_s), p(ZEROELEM), r(ZEROELEM), q(_siz, ZEROELEM), cnum(_siz), vlo(_siz), vhi(_siz) {};
};

template <typename VALTYPE>
class LWDPSubset{
public:
  std::vector<VALTYPE> wp;
  std::vector<VALTYPE> wn;
  
  int vsetcnt;
  std::vector<int> vtovset;
  std::vector<int> vsetlast;
  bool markedempty;
  
  std::vector<std::vector<DPBlockSubset<VALTYPE>>> dp;
  std::vector<std::unordered_map<StateSubset, size_t>> maps;
  
  LWDPSubset(){};
  
  void solve(Graph& G, const tdzdd::DdStructure<2>& base, const std::vector<VALTYPE>& _wp, const std::vector<VALTYPE>& _wn, const std::vector<int>& _vsets, std::vector<VALTYPE>& res, VALTYPE IDENTELEM, VALTYPE ZEROELEM);
  void innerSolve(const Graph& G, const tdzdd::DdStructure<2>& base, std::vector<VALTYPE>& res, VALTYPE IDENTELEM, VALTYPE ZEROELEM);
  void allClear(){
    dp.clear();
  }
};

template <typename VALTYPE>
void LWDPSubset<VALTYPE>::solve(Graph& G, const tdzdd::DdStructure<2>& base, const std::vector<VALTYPE>& _wp, const std::vector<VALTYPE>& _wn, const std::vector<int>& _vsets, std::vector<VALTYPE>& res, VALTYPE IDENTELEM, VALTYPE ZEROELEM){
  allClear();
  
  G.buildFrontiers();
  fprintf(stderr, "max_frontier_width = %d\n", G.maxFWidth());
  
  wp.resize(_wp.size());
  std::copy(_wp.begin(), _wp.end(), wp.begin());
  wn.resize(_wn.size());
  std::copy(_wn.begin(), _wn.end(), wn.begin());
  
  int n = G.numV();
  int m = G.numE();
  
  vtovset.resize(n+1);
  vsetcnt = 1;
  markedempty = true;
  for(const int v : _vsets){
    if(!v){
      ++vsetcnt;
    }else{
      vtovset[v] = vsetcnt;
      if(vsetcnt == 1) markedempty = false;
    }
  }
  --vsetcnt;
  
  vsetlast.resize(vsetcnt+1);
  for(int i=0; i<m; ++i){
    vsetlast[vtovset[G.e[i].first ]] = i;
    vsetlast[vtovset[G.e[i].second]] = i;
  }
  vsetlast[0] = -1;
  
  dp.resize(m+1);
  maps.resize(m+1);
  
  innerSolve(G, base, res, IDENTELEM, ZEROELEM);
}

template <typename VALTYPE>
void LWDPSubset<VALTYPE>::innerSolve(const Graph& G, const tdzdd::DdStructure<2>& base, std::vector<VALTYPE>& res, VALTYPE IDENTELEM, VALTYPE ZEROELEM){
  int n = G.numV();
  int m = G.numE();
  std::vector<int64_t> ssizes(m+1);
  
  maps.resize(m+1);
  
  // root node
  tdzdd::NodeId broot;
  if(base.getRoot(broot)){
    StateSubset root(broot.code());
    root.cnum = 0;
    root.comp.fill(-1);
    root.vset.fill(0);
    maps[0].emplace(root, 0);
    dp[0].emplace_back(root, 0, ZEROELEM);
    dp[0][0].p = IDENTELEM;
    ssizes[0] = 1;
  }
  
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    size_t kk = now_fro.size();
    size_t tt = med_fro.size();
    size_t ll = next_fro.size();
    maps[i+1].reserve(ssizes[i] * 2);
    dp[i+1].reserve(ssizes[i] * 2);
    ssizes[i+1] = 0;
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const StateSubset& now_state = ent.first;
      const auto& now_comp = now_state.comp;
      const auto& now_vset = now_state.vset;
      uint64_t now_bcode = now_state.bcode;
      tdzdd::NodeId now_bnode(now_bcode);
      int tdzddlevel = now_bnode.row();
      size_t blevel = static_cast<size_t>(m - tdzddlevel);
      int8_t cc = now_state.cnum;
      int8_t cc_old = cc;
      // generate intermediate state
      StateSubset med_state;
      auto& med_comp = med_state.comp;
      auto& med_vset = med_state.vset;
      med_comp.fill(-1);
      med_vset.fill(0);
      memcpy(med_comp.data(), now_comp.data(), tt);
      memcpy(med_vset.data(), now_vset.data(), cc);
      for(const auto& pos : now_ent){
        med_comp[pos] = cc;
        med_vset[cc++] = vtovset[med_fro[pos]];
      }
      int8_t minusvsetcnt = std::max(1, -*std::min_element(med_vset.begin(), med_vset.begin() + cc));
      {
        std::vector<int8_t> reminus(vsetcnt+1, -1);
        reminus[1] = 1;
        for(int8_t c=0; c<cc; ++c){
          if(med_vset[c] > 0 && i >= vsetlast[med_vset[c]]){
            if(reminus[med_vset[c]] < 0) reminus[med_vset[c]] = ++minusvsetcnt;
            med_vset[c] = -reminus[med_vset[c]];
          }
        }
      }
      std::vector<int8_t> minuscnt(minusvsetcnt+1);
      for(size_t c=0; c<cc; ++c){
        if(med_vset[c] < 0) ++minuscnt[-med_vset[c]];
      }
      
      // lo_state processing
      tdzdd::NodeId lo_bnode(now_bcode);
      if(i < blevel || base.getChild(lo_bnode, tdzddlevel, 0)){
        // generate lo_state
        StateSubset lo_state(lo_bnode.code());
        bool prune = false;
        auto& lo_comp = lo_state.comp;
        auto& lo_vset = lo_state.vset;
        lo_comp.fill(-1);
        lo_vset.fill(0);
        memcpy(lo_comp.data(), med_comp.data(), tt);
        std::vector<int8_t> erased;
        erased.reserve(2);
        for(const auto& pos : now_lve){
          lo_comp[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : lo_comp){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        for(const auto& pos : now_lve){
          if(renum[med_comp[pos]] < 0) erased.emplace_back(med_comp[pos]);
        }
        if(erased.size() >= 2 && erased[0] == erased[1]) erased.pop_back();
        for(const auto& c : erased){
          if(med_vset[c] != 0){
            if(med_vset[c] > 0) prune = true;
            else if(minuscnt[-med_vset[c]] >= 2) prune = true;
          }
        }
        for(int8_t c=0; c<cc; ++c){
          if(renum[c] >= 0) lo_vset[renum[c]] = med_vset[c];
        }
        {
          std::vector<int8_t> minusrenum(minusvsetcnt+1, -1);
          int8_t newminusvsetcnt = 1;
          minusrenum[1] = 1;
          for(int8_t c=0; c<cc_new; ++c){
            if(lo_vset[c] >= 0) continue;
            if(minusrenum[-lo_vset[c]] < 0) minusrenum[-lo_vset[c]] = ++newminusvsetcnt;
            lo_vset[c] = -minusrenum[-lo_vset[c]];
          }
        }
        lo_state.cnum = cc_new;
        
        // find or generate id
        int64_t lo_id;
        if(prune){
          lo_id = -1;
        }else{
          auto it = maps[i+1].find(lo_state);
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            lo_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(lo_state, ssizes[i+1]);
            lo_id = ssizes[i+1]++;
            dp[i+1].emplace_back(lo_state, cc_new, ZEROELEM);
          }
        }
        dp[i][now_id].lo = lo_id;
        
        if(!prune){
          // vlo equals renum
          std::copy(renum.begin(), renum.begin() + cc_old, dp[i][now_id].vlo.begin());
          // special treatment for erased components
          for(const auto& c : erased){
            if(c >= cc_old) continue;
            if(med_vset[c] == -1 || (med_vset[c] == 0 && markedempty)) dp[i][now_id].vlo[c] = -2;
          }
        }
      }else{
        dp[i][now_id].lo = -1;
      }
      
      // hi_state processing
      int8_t cat_to   = med_comp[now_vpos.first];
      int8_t cat_from = med_comp[now_vpos.second];
      tdzdd::NodeId hi_bnode(now_bcode);
      if(i == blevel && base.getChild(hi_bnode, tdzddlevel, 1)){
        // modify med_state
        bool prune = false;
        if(cat_from != cat_to){
          if(med_vset[cat_from] != 0 && med_vset[cat_to] != 0){
            if(med_vset[cat_from] != med_vset[cat_to]) prune = true;
            else if(med_vset[cat_from] < 0) --minuscnt[-med_vset[cat_from]];
          }else if(med_vset[cat_from] != 0){
            med_vset[cat_to] = med_vset[cat_from];
          }
          for(auto&& val : med_comp){
            if(val == cat_from) val = cat_to;
          }
        }
        // generate hi_state
        StateSubset hi_state(hi_bnode.code());
        auto& hi_comp = hi_state.comp;
        auto& hi_vset = hi_state.vset;
        hi_comp.fill(-1);
        hi_vset.fill(0);
        memcpy(hi_comp.data(), med_comp.data(), tt);
        std::vector<int8_t> erased;
        erased.reserve(2);
        for(const auto& pos : now_lve){
          hi_comp[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : hi_comp){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        for(const auto& pos : now_lve){
          if(renum[med_comp[pos]] < 0) erased.emplace_back(med_comp[pos]);
        }
        if(erased.size() >= 2 && erased[0] == erased[1]) erased.pop_back();
        for(const auto& c : erased){
          if(med_vset[c] != 0){
            if(med_vset[c] > 0) prune = true;
            else if(minuscnt[-med_vset[c]] >= 2) prune = true;
          }
        }
        for(int8_t c=0; c<cc; ++c){
          if(renum[c] >= 0) hi_vset[renum[c]] = med_vset[c];
        }
        renum[cat_from] = renum[cat_to];
        {
          std::vector<int8_t> minusrenum(minusvsetcnt+1, -1);
          int8_t newminusvsetcnt = 1;
          minusrenum[1] = 1;
          for(int8_t c=0; c<cc_new; ++c){
            if(hi_vset[c] >= 0) continue;
            if(minusrenum[-hi_vset[c]] < 0) minusrenum[-hi_vset[c]] = ++newminusvsetcnt;
            hi_vset[c] = -minusrenum[-hi_vset[c]];
          }
        }
        hi_state.cnum = cc_new;
        
        // find or generate id
        int64_t hi_id;
        if(prune){
          hi_id = -1;
        }else{
          auto it = maps[i+1].find(hi_state);
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            hi_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(hi_state, ssizes[i+1]);
            hi_id = ssizes[i+1]++;
            dp[i+1].emplace_back(hi_state, cc_new, ZEROELEM);
          }
        }
        dp[i][now_id].hi = hi_id;
        
        if(!prune){
          // vhi equals renum
          std::copy(renum.begin(), renum.begin() + cc_old, dp[i][now_id].vhi.begin());
          // special treatment for erased components
          for(const auto& c : erased){
            if(med_vset[c] == -1 || (med_vset[c] == 0 && markedempty)){
              if(c < cc_old)dp[i][now_id].vhi[c] = -2;
              if(c == cat_to && cat_from < cc_old) dp[i][now_id].vhi[cat_from] = -2;
            }
          }
        }
      }else{
        dp[i][now_id].hi = -1;
      }
    }
    dp[i+1].shrink_to_fit();
  }
  
  // dp computation
  // top-down computation of p
  for(size_t i=0; i<m; ++i){
    for(const auto& ent : dp[i]){
      if(ent.lo >= 0) dp[i+1][ent.lo].p += wn[i] * ent.p;
      if(ent.hi >= 0) dp[i+1][ent.hi].p += wp[i] * ent.p;
    }
  }
  
  // bottom-up computation of r and q
  if(ssizes[m] == 0){
    res.resize(n+1);
    std::fill(res.begin(), res.end(), ZEROELEM);
    return;
  }
  dp[m][0].r = IDENTELEM;
  for(size_t i=m-1; i>=1; --i){
    for(auto&& ent : dp[i]){
      // computation of r
      if(ent.lo >= 0) ent.r  = wn[i] * dp[i+1][ent.lo].r;
      if(ent.hi >= 0) ent.r += wp[i] * dp[i+1][ent.hi].r;
      // computation of q
      for(int8_t c=0; c<ent.cnum; ++c){
        if(ent.lo >= 0) ent.q[c]  = wn[i] * (ent.vlo[c] >= 0 ? dp[i+1][ent.lo].q[ent.vlo[c]] : (ent.vlo[c] == -2 ? dp[i+1][ent.lo].r : ZEROELEM));
        if(ent.hi >= 0) ent.q[c] += wp[i] * (ent.vhi[c] >= 0 ? dp[i+1][ent.hi].q[ent.vhi[c]] : (ent.vhi[c] == -2 ? dp[i+1][ent.hi].r : ZEROELEM));
      }
    }
  }
  
  // levelwise computation
  res.resize(n+1);
  std::fill(res.begin(), res.end(), ZEROELEM);
  std::vector<int> computed(n+1);
  computed[0] = 1;
  // vertices present in vset
  for(int i=1; i<=n; ++i){
    if(vtovset[i]){
      if(vtovset[i] == 1) res[i] = dp[m][0].p;
      computed[i] = 1;
    }
  }
  // vertices with degree >= 2
  for(size_t i=0; i<m; ++i){
    std::vector<size_t> compute_pos;
    const auto& now_fro = G.fros[i];
    for(size_t pos=0; pos<now_fro.size(); ++pos){
      int now_v = now_fro[pos];
      if(now_v >= 1 && !computed[now_v]){
        compute_pos.emplace_back(pos);
        computed[now_v] = 1;
        res[now_v] = ZEROELEM;
      }
    }
    for(const auto& ent : dp[i]){
      for(size_t pos : compute_pos){
        res[now_fro[pos]] += ent.p * ent.q[ent.s.comp[pos]];
      }
    }
  }
  // vertices with degree = 1
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    
    size_t pos, pos2;
    int tv;
    bool flg = false;
    if((now_vpos.first >= now_fro.size() || now_fro[now_vpos.first] <= 0) && (now_vpos.first >= next_fro.size() || next_fro[now_vpos.first] <= 0)){
      flg = true;
      pos = now_vpos.first;
      pos2 = now_vpos.second;
      tv = med_fro[pos];
    }
    if((now_vpos.second >= now_fro.size() || now_fro[now_vpos.second] <= 0) && (now_vpos.second >= next_fro.size() || next_fro[now_vpos.second] <= 0)){
      flg = true;
      pos = now_vpos.second;
      pos2 = now_vpos.first;
      tv = med_fro[pos];
    }
    if(flg){
      if(computed[tv]) continue;
      res[tv] = ZEROELEM;
      if(pos2 < next_fro.size() && next_fro[pos2] > 0){
        for(const auto& ent : dp[i]){
          if(ent.hi >= 0) res[tv] += ent.p * wp[i] * dp[i+1][ent.hi].q[dp[i+1][ent.hi].s.comp[pos2]];
          if(markedempty && ent.lo >= 0) res[tv] += ent.p * wn[i] * dp[i+1][ent.lo].r;
        }
      }else{
        for(const auto& ent : dp[i]){
          int8_t c = ent.s.comp[pos2];
          if(ent.hi >= 0) res[tv] += ent.p * wp[i] * (ent.vhi[c] >= 0 ? dp[i+1][ent.hi].q[ent.vhi[c]] : (ent.vhi[c] == -2 ? dp[i+1][ent.hi].r : ZEROELEM));
          if(markedempty && ent.lo >= 0) res[tv] += ent.p * wn[i] * dp[i+1][ent.lo].r;
        }
      }
    }
  }

  size_t totalsize = 0;
  for(size_t i=0; i<=m; ++i) totalsize += ssizes[i];
  fprintf(stderr, "#(states) = %zu\n", totalsize);
}

#endif // SIMFRONTIER_LWDPSUBSET_HPP
