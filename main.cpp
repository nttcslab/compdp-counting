#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/lwdpsubset.hpp"

#include "tdzdd/DdSpec.hpp"
#include "tdzdd/DdSpecOp.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/util/Graph.hpp"
#include "tdzdd/spec/FrontierBasedSearch.hpp"
#include "tdzdd/util/IntSubset.hpp"
#include "tdzdd/spec/DegreeConstraint.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <iostream>
#include <algorithm>

#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [mode] [graph_file] [point_file]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 3){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  
  int mode = atoi(argv[1]);

  if(!G.readfromFile(argv[2])){
    fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[2]);
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  std::vector<int> vs;
  int tmp;
  {
    FILE *fp;
    if((fp = fopen(argv[3], "r")) == NULL){
      fprintf(stderr, "ERROR: reading point file file %s failed.\n", argv[3]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    while(fscanf(fp, "%d", &tmp) != EOF){
      vs.emplace_back(tmp);
    }
    fclose(fp);
  }
  
  auto cstart = std::chrono::system_clock::now();
  
  std::vector<int> vsets;
  
  n = G.numV();
  m = G.numE();
  
  tdzdd::Graph tG;
  for(const auto& edg : G.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  tG.update();
  
  std::vector<mp::uint512_t> wn(m, 1);
  std::vector<mp::uint512_t> wp(m, 1);
  std::vector<mp::uint512_t> res;
  LWDPSubset<mp::uint512_t> LWSolver;
  
  if(mode == 0){       // path
    int src = vs[0];
    int dst = vs[1];
    vsets.emplace_back(src);
    vsets.emplace_back(dst);
    vsets.emplace_back(0);
    
    tdzdd::FrontierBasedSearch fbs(tG, -1, true, true);
    tdzdd::IntRange OneOnly(1, 1, 1);
    tdzdd::IntRange ZeroOrTwo(0, 2, 2);
    tdzdd::DegreeConstraint dc(tG, &ZeroOrTwo, true);
    dc.setConstraint(std::to_string(src), &OneOnly);
    dc.setConstraint(std::to_string(dst), &OneOnly);
    
    tdzdd::DdStructure<2> basedd(tdzdd::zddIntersection(dc, fbs));
    basedd.zddReduce();
    LWSolver.solve(G, basedd, wp, wn, vsets, res, 1, 0);
  }else if(mode == 1){ // cycle
    int src = vs[0];
    vsets.emplace_back(src);
    vsets.emplace_back(0);
    
    tG.setColor(std::to_string(src), 1);
    tG.update();
    tdzdd::FrontierBasedSearch fbs(tG, 0, false, true);
    tdzdd::IntRange ZeroOrTwo(0, 2, 2);
    tdzdd::DegreeConstraint dc(tG, &ZeroOrTwo, true);
    
    tdzdd::DdStructure<2> basedd(tdzdd::zddIntersection(dc, fbs));
    basedd.zddReduce();
    LWSolver.solve(G, basedd, wp, wn, vsets, res, 1, 0);
  }else if(mode == 2){ // Steiner tree
    vsets.resize(vs.size()+1);
    std::copy(vs.begin(), vs.end(), vsets.begin());
    
    tdzdd::FrontierBasedSearch fbs(tG, 1, true, true);
    
    tdzdd::DdStructure<2> basedd(fbs);
    basedd.zddReduce();
    LWSolver.solve(G, basedd, wp, wn, vsets, res, 1, 0);
  }else{               // rooted spanning tree
    int cnt = static_cast<int>(vs.size());
    vsets.resize(cnt*2+1);
    for(int i=0; i<cnt; ++i){
      vsets[i*2+1] = vs[i];
    }
    
    for(int i=0; i<cnt; ++i){
      tG.setColor(std::to_string(vsets[i*2+1]), i+1);
    }
    tG.update();
    tdzdd::IntRange MoreThanOne(1);
    tdzdd::DegreeConstraint dc(tG, &MoreThanOne, true);
    tdzdd::FrontierBasedSearch fbs(tG, 1, true, true);
    
    tdzdd::DdStructure<2> basedd(tdzdd::zddIntersection(dc, fbs));
    basedd.zddReduce();
    LWSolver.solve(G, basedd, wp, wn, vsets, res, 1, 0);
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  for(int i=1; i<=n; ++i){
    std::cout << res[i] << std::endl;
  }
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}