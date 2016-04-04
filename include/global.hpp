#ifndef global_hpp
#define global_hpp

#include <unordered_map>
#include <vector>
#include "graph.hpp"

extern std::vector<std::unordered_map<int, int> > vid2bid;
extern std::vector<std::vector<int> > bid2vid;

namespace Query4 {
   typedef uint32_t PersonId; // Type for person ids
   typedef Graph<PersonId> PersonSubgraph;
}

#endif
