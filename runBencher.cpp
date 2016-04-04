//Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
//Code must not be used, distributed, without written consent by the authors
#include "include/global.hpp"
#include "include/log.hpp"
#include "include/worker.hpp"
#include "extra.hpp"

#include <string>
#include <stdlib.h>
#include <unordered_map>
#include <fstream>
#include <vector>

#define GEN_BFS_TASK(X, WIDTH) \
   X(width==WIDTH) { \
      bfsRunner = new CustomRunner<WIDTH>(); \
   } \


std::vector<std::unordered_map<int, int> > vid2bid;
std::vector<std::vector<int> > bid2vid;

bool build_source_vector(const std::string& sources, const size_t numThreads, const int numBFSs){
   bid2vid.resize(numThreads);
   std::ifstream infile(sources);
   if(infile.is_open()){
      int a, b = 0;
      size_t t = 0;
      while (infile >> a && t < numThreads){
         bid2vid[t].push_back(a);
         b ++;
         if(b == numBFSs){
        	 t ++;
        	 b = 0;
         }
      }

      infile.close();

      if(t < numThreads)
         return false;

      return true;
   } else {
      return false;
   }
}

void build_vid2bid(const size_t numThreads, const Query4::PersonSubgraph& subgraph){
	vid2bid.resize(numThreads);
	for(size_t i = 0; i < numThreads; i++){
		for(size_t j = 0; j < bid2vid[i].size(); j ++){
			int externalId = subgraph.mapExternalNodeId(bid2vid[i][j]);
			vid2bid[i].insert(std::make_pair(externalId, j));
		}
	}
}


int main(int argc, char** argv) {
   if(argc < 5) {
      FATAL_ERROR("Not enough parameters");
   }

   std::string graph = argv[1];
   std::string sources = argv[2];
   size_t numThreads = std::stoi(argv[3]);
   int numBFSs = std::stoi(argv[4]);
   std::string outprefix = "";
   if(argc > 5)
	   outprefix = argv[5];

   const int width = numBFSs / 64;

   // build source map
   if(numBFSs % 64 != 0){
      FATAL_ERROR("Number of concurrent BFSs should be multiples of 64");
   }
   if(numBFSs > 1024){
      FATAL_ERROR("support maximum 1024 concurrent BFSs per thread");
   }
   if(!build_source_vector(sources, numThreads, numBFSs)){
      FATAL_ERROR("No enough sources");
   }

   LOG_PRINT("[Main] Using "<< numThreads <<" threads, each running " << numBFSs << " concurrent BFSs");

   VirtualRunner *bfsRunner = NULL;
   GEN_BFS_TASK(if, 1)
   GEN_BFS_TASK(else if, 2)
   GEN_BFS_TASK(else if, 3)
   GEN_BFS_TASK(else if, 4)
   GEN_BFS_TASK(else if, 5)
   GEN_BFS_TASK(else if, 6)
   GEN_BFS_TASK(else if, 7)
   GEN_BFS_TASK(else if, 8)
   GEN_BFS_TASK(else if, 9)
   GEN_BFS_TASK(else if, 10)
   GEN_BFS_TASK(else if, 11)
   GEN_BFS_TASK(else if, 12)
   GEN_BFS_TASK(else if, 13)
   GEN_BFS_TASK(else if, 14)
   GEN_BFS_TASK(else if, 15)
   GEN_BFS_TASK(else if, 16)

   if(bfsRunner != NULL){
	   auto personGraph = Graph<Query4::PersonId>::loadFromPath(graph);

	   /*
	   for(int i = 0; i < numThreads; i++){
		   const std::unordered_map<int, int>& v2b = vid2bid[i];
		   const std::vector<int>& b2v = bid2vid[i];
		   for(std::unordered_map<int, int>::const_iterator it = v2b.begin(); it != v2b.end(); it++){
			   printf("%d  %d - %d\n", i, it->first, b2v[it->second]);
		   }
	   }
	   */

	   tschrono::Time start = tschrono::now();
	   // Allocate additional worker threads
	   Workers workers(numThreads-1);

	   bfsRunner -> run(numThreads, workers, personGraph, outprefix);

	   workers.close();

	   tschrono::Time runtime = tschrono::now() - start;
	   LOG_PRINT("finish in " << runtime << " ms");

   }
   return 0;
}
