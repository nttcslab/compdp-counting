#include "mylib/common.hpp"
#include "mylib/graph.hpp"

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
#include <numeric>

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
  
  n = G.numV();
  m = G.numE();
  
  std::vector<std::string> res(n+1);
  std::vector<size_t> ddsize;
  
  tdzdd::Graph tG;
  for(const auto& edg : G.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  tG.update();
  
  if(mode == 0){       // path
    int src = vs[0];
    int dst = vs[1];
    
    tdzdd::IntRange OneOnly(1, 1, 1);
    tdzdd::IntRange ZeroOrTwo(0, 2, 2);
    tdzdd::DegreeConstraint dc(tG, &ZeroOrTwo, true);
    dc.setConstraint(std::to_string(src), &OneOnly);
    dc.setConstraint(std::to_string(dst), &OneOnly);
    tdzdd::DdStructure<2> basedd(dc);
    basedd.zddReduce();
    
    for(int i=1; i<=n; ++i){
      tG.clearColors();
      tG.setColor(std::to_string(src), 1);
      tG.setColor(std::to_string(dst), 1);
      if(i != src && i != dst) tG.setColor(std::to_string(i), 1);
      tG.update();
      tdzdd::FrontierBasedSearch fbs(tG, -1, true, true);
      tdzdd::DdStructure<2> compdd(tdzdd::zddIntersection(basedd, fbs));
      ddsize.emplace_back(compdd.size());
      compdd.zddReduce();
      res[i] = compdd.zddCardinality();
    }
  }else if(mode == 1){ // cycle
    int src = vs[0];
    
    tdzdd::IntRange ZeroOrTwo(0, 2, 2);
    tdzdd::DegreeConstraint dc(tG, &ZeroOrTwo, true);
    tdzdd::DdStructure<2> basedd(dc);
    basedd.zddReduce();
    
    for(int i=1; i<=n; ++i){
      tG.clearColors();
      tG.setColor(std::to_string(src), 1);
      if(i != src) tG.setColor(std::to_string(i), 1);
      tG.update();
      tdzdd::FrontierBasedSearch fbs(tG, 0, false, true);
      tdzdd::DdStructure<2> compdd(tdzdd::zddIntersection(basedd, fbs));
      ddsize.emplace_back(compdd.size());
      compdd.zddReduce();
      res[i] = compdd.zddCardinality();
    }
  }else if(mode == 2){ // Steiner tree
    for(int i=1; i<=n; ++i){
      tG.clearColors();
      bool exist = false;
      for(int v : vs){
        tG.setColor(std::to_string(v), 1);
        if(v == i) exist = true;
      }
      if(!exist) tG.setColor(std::to_string(i), 1);
      tG.update();
      tdzdd::FrontierBasedSearch fbs(tG, 0, true, true);
      tdzdd::DdStructure<2> compdd(fbs);
      ddsize.emplace_back(compdd.size());
      compdd.zddReduce();
      res[i] = compdd.zddCardinality();
    }
  }else{               // rooted spanning forest
    int numc = static_cast<int>(vs.size());
    tdzdd::IntRange MoreThanOne(1);
    tdzdd::DegreeConstraint dc(tG, &MoreThanOne, true);
    tdzdd::DdStructure<2> basedd(dc);
    basedd.zddReduce();
    
    for(int i=1; i<=n; ++i){
      tG.clearColors();
      bool exist = false;
      for(int k=0; k<numc; ++k){
        tG.setColor(std::to_string(vs[k]), k+1);
        if(vs[k] == i) exist = true;
      }
      if(exist){
        res[i] = "0";
        continue;
      }
      tG.setColor(std::to_string(i), numc+1);
      tG.update();
      tdzdd::FrontierBasedSearch fbs(tG, 0, true, true);
      tdzdd::DdStructure<2> compdd(tdzdd::zddIntersection(basedd, fbs));
      ddsize.emplace_back(compdd.size());
      compdd.zddReduce();
      res[i] = compdd.zddCardinality();
    }
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  for(int i=1; i<=n; ++i){
    std::cout << res[i] << std::endl;
  }
  
  fprintf(stderr, "max size   : %zu\n", *std::max_element(ddsize.begin(), ddsize.end()));
  fprintf(stderr, "total size : %zu\n", std::accumulate(ddsize.begin(), ddsize.end(), static_cast<size_t>(0)));
  fprintf(stderr, "calc time  : %.6lf ms\n", ctime);
  
  return 0;
}