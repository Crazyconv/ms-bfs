#ifndef global_hpp
#define global_hpp

#include <unordered_map>
#include <vector>
#include "graph.hpp"

extern std::vector<std::unordered_map<uint64_t, int> > vid2bid;
extern std::vector<std::vector<uint64_t> > bid2vid;
extern std::vector<uint64_t> exe_time;
extern std::vector<uint64_t> write_time;


namespace Query4 {
   typedef uint64_t PersonId; // Type for person ids
   typedef Graph<PersonId> PersonSubgraph;
}

#endif
