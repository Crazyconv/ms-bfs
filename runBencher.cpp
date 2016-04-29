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

#define GEN_BFS_RUNNER(X, WIDTH) \
   X(width==WIDTH) { \
      bfsRunner = new CustomRunner<WIDTH>(); \
   } \


/* std::vector<std::unordered_map<int, int> > vid2bid; */
// each thread has a mapping from bit position to source id
std::vector<std::vector<uint64_t> > bid2vid;
std::vector<uint64_t> exe_time;
std::vector<uint64_t> write_time;

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

/*
void build_vid2bid(const size_t numThreads, const Query4::PersonSubgraph& subgraph){
	vid2bid.resize(numThreads);
	for(size_t i = 0; i < numThreads; i++){
		for(size_t j = 0; j < bid2vid[i].size(); j ++){
			int externalId = subgraph.mapExternalNodeId(bid2vid[i][j]);
			vid2bid[i].insert(std::make_pair(externalId, j));
		}
	}
}
*/


int main(int argc, char** argv) {
   if(argc < 5) {
      FATAL_ERROR("Not enough parameters");
   }

   std::string graph = argv[1];
   std::string sources = argv[2];
   size_t numThreads = std::stoi(argv[3]);
   int numBFSs = std::stoi(argv[4]);
   std::string stat_file = "";
   std::string outprefix = "";
   if(argc > 5)
   	   stat_file = argv[5];
   if(argc > 6)
	   outprefix = argv[6];

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
   /* We need a constant (here is 1 ~ 16) for the template variable in CustomRunner<WIDTH> */
   GEN_BFS_RUNNER(if, 1)
   GEN_BFS_RUNNER(else if, 2)
   GEN_BFS_RUNNER(else if, 3)
   GEN_BFS_RUNNER(else if, 4)
   GEN_BFS_RUNNER(else if, 5)
   GEN_BFS_RUNNER(else if, 6)
   GEN_BFS_RUNNER(else if, 7)
   GEN_BFS_RUNNER(else if, 8)
   GEN_BFS_RUNNER(else if, 9)
   GEN_BFS_RUNNER(else if, 10)
   GEN_BFS_RUNNER(else if, 11)
   GEN_BFS_RUNNER(else if, 12)
   GEN_BFS_RUNNER(else if, 13)
   GEN_BFS_RUNNER(else if, 14)
   GEN_BFS_RUNNER(else if, 15)
   GEN_BFS_RUNNER(else if, 16)

   if(bfsRunner != NULL){
	   tschrono::Time start = tschrono::now();

	   // load graph
	   auto personGraph = Graph<Query4::PersonId>::loadFromPath(graph);

	   tschrono::Time load_time = tschrono::now() - start;

	   // Allocate additional worker threads
	   Workers workers(numThreads-1);
	   exe_time.resize(numThreads);
	   write_time.resize(numThreads);

	   // perform BFS
	   bfsRunner -> run(numThreads, workers, personGraph, outprefix);

	   workers.close();

	   if(stat_file.length() > 0){

			std::ofstream outfile(stat_file);
			if(outfile.is_open()){
				outfile << numThreads * numBFSs << "," << numThreads << "," << numBFSs << ","
						<< load_time/1000.0 << "," << (*(std::max_element(exe_time.begin(), exe_time.end())))/1000.0 << ","
						<< (*(std::max_element(write_time.begin(), write_time.end())))/1000.0 << std::endl;
				outfile.close();
			} else {
				LOG_PRINT("Unable to open file " << stat_file);
			}
	   }

   }
   return 0;
}
