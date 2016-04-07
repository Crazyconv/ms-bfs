#ifndef EXTRA_HPP
#define EXTRA_HPP

#include "include/global.hpp"
#include "include/bfs/batchdistance.hpp"
#include "include/bfs/bitops.hpp"
#include "include/graph.hpp"
#include "include/worker.hpp"
#include "include/scheduler.hpp"
#include "include/log.hpp"

#include <array>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <fstream>


#define SORTED_NEIGHBOR_PROCESSING
#define DO_PREFETCH

// data structure for bit field
template<uint64_t width>
struct BatchBits {
    static const size_t TYPE_BITS_COUNT = 64;
    uint64_t data[width];

    BatchBits() : data() {
    }

    bool isAllZero() const {
        for (unsigned i = 0; i < width; ++i) {
            if(BitBaseOp<uint64_t>::notZero(data[i])) {
                return false;
            }
        }
        return true;
    }

    size_t count() const {
        size_t sum=0;
        for (unsigned i = 0; i < width; ++i) {
            sum+=BitBaseOp<uint64_t>::popCount(data[i]);
        }
        return sum;
    }

    void negate() {
        for (unsigned i = 0; i < width; ++i) {
            data[i] = ~data[i];
        }
    }

    void setBit(const size_t ix) {
        auto field = ix/TYPE_BITS_COUNT;
        auto field_bit = ix-(field*TYPE_BITS_COUNT);
        data[field] |= BitBaseOp<uint64_t>::getSetMask(field_bit);
    }
};

static const unsigned int PREFETCH=38;

// a BFSTask is executed by a thread, the execution is defined in void operator()()
template<uint64_t width>
struct BFSTask{
    static const size_t TYPE_BITS = 64;
    static const size_t WIDTH = width;
    static const size_t BATCH_BITS_COUNT = TYPE_BITS * WIDTH;
    typedef BatchBits<width> Bitset;

    size_t index;                                         /* thread id */
    const Query4::PersonSubgraph& subgraph;               /* graph data */
    /* const std::unordered_map<int, int>& v2b; */
    const std::vector<int>& b2v;                          /* bit pos to source id mapping */
    const int subgraphSize;                               /* graph size */
    const std::string& output;                            /* output file name prefix */

    int totalReachable[BATCH_BITS_COUNT];                 /* #reachable vertices from a souuce */

public:
    BFSTask(size_t index, const Query4::PersonSubgraph& subgraph, const std::string& output): index(index), subgraph(subgraph),
		/*v2b(vid2bid[index]),*/ b2v(bid2vid[index]), subgraphSize(subgraph.size()), output(output){

    	// initialize totalReachable
        for(size_t i = 0; i < BATCH_BITS_COUNT; i++){
        	int externalId = exId(b2v[i]);
        	/* personComponents[i]: id of the component that vertex i belongs to
        	 * componentSize[j]: size of component j
        	 */
        	totalReachable[i] = subgraph.componentSizes[subgraph.personComponents[externalId]] - 1;
        }
    }
    inline int exId(int vid){
    	/* mapping from internal id (real id, read from graph file) to external id (consecutive id) */
    	return subgraph.mapExternalNodeId(vid);
    }
    void operator()(){

        // Initialize visit lists
        std::array<Bitset*,2> visitLists;
        for(int a=0; a<2; a++) {
           const auto ret=posix_memalign(reinterpret_cast<void**>(&(visitLists[a])),64,sizeof(Bitset)*subgraphSize);
           if(unlikely(ret!=0)) {
              throw -1;
           }
           new(visitLists[a]) Bitset[subgraphSize]();
        }

        // Initialize seen vector
        Query4::PersonId minPerson = std::numeric_limits<Query4::PersonId>::max();
        Bitset* seen;
        const auto ret=posix_memalign(reinterpret_cast<void**>(&seen),64,sizeof(Bitset)*subgraphSize);
        if(unlikely(ret!=0)) {
           throw -1;
        }
        new(seen) Bitset[subgraphSize]();

        // Initialize distance vector
        /* dist[subgraphSize][BATCH_BITS_COUNT] */
        int ** dist = new int*[subgraphSize];
        for(int i = 0; i < subgraphSize; i++){
            dist[i] = new int[BATCH_BITS_COUNT];
            memset(dist[i],-1,BATCH_BITS_COUNT*sizeof(int));
        }

        // Initialize reachable vector
        int * reachable = new int[BATCH_BITS_COUNT];

        // Initialize active queries
        for(int i = 0; i < BATCH_BITS_COUNT; i++){
        	int externalId = exId(b2v[i]);
            seen[externalId].setBit(i);
            visitLists[0][externalId].setBit(i);
            minPerson = std::min(minPerson, (Query4::PersonId)externalId);
            dist[externalId][i] = 0;
        }

        // Initialize iteration workstate
        /* bit 1 means that BFS has not complete */
        Bitset processQuery;
        processQuery.negate();

        /* some data structure to check #vertices discovered by each BFS */
        /* don't really understand the internal logic */
        uint32_t queriesToProcess=BATCH_BITS_COUNT;
        alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT];
        memset(numDistDiscovered,0,BATCH_BITS_COUNT*sizeof(uint32_t));
        Query4::BatchDistance<uint64_t, width> batchDist(numDistDiscovered);

        size_t curToVisitQueue = 0;
        uint32_t nextDistance = 1;

        Query4::PersonId startPerson=minPerson;

        // Run iterations
        do {
            Bitset* const toVisit = visitLists[curToVisitQueue];
            Bitset* const nextToVisit = visitLists[1-curToVisitQueue];

            runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, dist, batchDist, processQuery, nextDistance);

			#ifdef DEBUG
			uint64_t newReached = 0;
			uint64_t numNotZero = 0;
			#endif /*DEBUG*/

			for(uint32_t pos=0; pos<BATCH_BITS_COUNT; pos++) {
			    /* update processQuery
			     * if a bfs has discovered all vertices reachable, change the corresponding bit to 1
			     * and decrease queriesToProcess
			     * */
                updateProcessQuery(processQuery, pos, numDistDiscovered[pos], reachable, queriesToProcess);

				#ifdef DEBUG
				if(numDistDiscovered[pos]>0) {
				   newReached += numDistDiscovered[pos];
				   numNotZero++;
				}
				#endif /*DEBUG*/
            }

			/* if no more to process, break from the loop */
            if(queriesToProcess==0) {
                break;
            }
            nextDistance++;

			#ifdef DEBUG
            assert(newReached!=0);
			#endif /*DEBUG*/

            // Reset iteration state
            memset(numDistDiscovered,0,BATCH_BITS_COUNT*sizeof(uint32_t));

            // Swap queues
            startPerson = 0;
            curToVisitQueue = 1-curToVisitQueue;
        } while(true);

        if(output.length() > 0){
        	std::ofstream outfile(output + "_" + std::to_string(index) + ".txt");
        	if(outfile.is_open()){
    			for(size_t i = 0; i < BATCH_BITS_COUNT; i++){
    				for(int j = 0; j < subgraphSize; j++){
    					outfile << b2v[i] << " " << subgraph.mapInternalNodeId(j) << " " << dist[j][i] << std::endl;
//    					printf("%d %d %d\n", b2v[i], subgraph.mapInternalNodeId(j), dist[j][i]);
    				}
    			}
        	    outfile.close();
			} else {
				LOG_PRINT("Unable to open file " << output << "_" << index << ".txt");
			}
        }

        free(seen);
        free(visitLists[0]);
        free(visitLists[1]);
        for(int i = 0; i < subgraphSize; i++){
            delete[] dist[i];
        }
        delete[] dist;
        delete[] reachable;
    }


    void updateProcessQuery(Bitset& processQuery, const uint32_t pos, const uint32_t numDiscovered,
        int* reachable, uint32_t& queriesToProcess) {

        auto field = pos/Bitset::TYPE_BITS_COUNT;
        auto field_bit = pos-(field*Bitset::TYPE_BITS_COUNT);

        if(BitBaseOp<uint64_t>::notZero(processQuery.data[field] & BitBaseOp<uint64_t>::getSetMask(field_bit))) {

            if(totalReachable[pos] == reachable[pos] || numDiscovered==0) {
                processQuery.data[field] = BitBaseOp<uint64_t>::andNot(processQuery.data[field], BitBaseOp<uint64_t>::getSetMask(field_bit));
                queriesToProcess--;
            }
        } else {
            assert(numDiscovered == 0);
        }
    }


    Bitset createVisitList(const Bitset& visitList, const Bitset& processQuery) {
        Bitset validVisit;

        bool nonZero=false;
        for(unsigned i=0; i<width; i++) {
            if(BitBaseOp<uint64_t>::notZero(visitList.data[i])) {
                nonZero=true;
                break;
            }
        }
        if(nonZero) {
            for(unsigned i=0; i<width; i++) {
                validVisit.data[i] = visitList.data[i] & processQuery.data[i];
            }
        }
        return validVisit;
    }

    void updateDist(unsigned field, uint64_t newVisits, int* dist, uint32_t nextDistance){
    	uint64_t lone = 1;
        for(uint64_t i = 0; i < 64 && newVisits != 0; i++){
            if((newVisits & (lone << i)) != 0){
                dist[field*64+i] = nextDistance;
                newVisits ^= (lone << i);
            }
        }
    }

	#ifdef SORTED_NEIGHBOR_PROCESSING
    void runBatchRound(const Query4::PersonSubgraph& subgraph, const Query4::PersonId startPerson, const Query4::PersonId limit,
            Bitset* visitList, Bitset* nextVisitList, Bitset* __restrict__ seen, int** dist,
    		Query4::BatchDistance<uint64_t, width>& batchDist, const Bitset processQuery, uint32_t nextDistance) {

		#ifdef DO_PREFETCH
		const int p2=std::min(PREFETCH, (unsigned int)(limit-startPerson));
		for(int a=1; a<p2; a++) {
		   __builtin_prefetch(visitList + a,0);
		}
		#endif /*DO_PREFETCH*/

		/* for each vertex */
		for (Query4::PersonId curPerson = startPerson; curPerson<limit; ++curPerson) {
		   auto curVisit = visitList[curPerson];

		   #ifdef DO_PREFETCH
		   if(curPerson+PREFETCH < limit) {
			  __builtin_prefetch(visitList + curPerson + PREFETCH,0);
		   }
		   #endif

		   bool zero=true;
		   for(int i=0; i<width; i++) {
			  if(BitBaseOp<uint64_t>::notZero(curVisit.data[i])) {
				 zero=false;
				 break;
			  }
		   }
		   if(zero) {
			  continue;
		   }

		   /* neighbors */
		   const auto& curFriends=*subgraph.retrieve(curPerson);
		   auto friendsBounds = curFriends.bounds();

		   #ifdef DO_PREFETCH
		   const int p=std::min(PREFETCH, (unsigned int)(friendsBounds.second-friendsBounds.first));
		   for(int a=1; a<p; a++) {
			  __builtin_prefetch(nextVisitList + *(friendsBounds.first+a),1);
		   }
		   #endif /*DO_PREFETCH*/

		   for(int i=0; i<width; i++) {
			  curVisit.data[i] &= processQuery.data[i];
		   }

		   /* for each neighbor, update nextVisitList */
		   while(friendsBounds.first != friendsBounds.second) {
			  #ifdef DO_PREFETCH
			  if(friendsBounds.first+PREFETCH < friendsBounds.second) {
				 __builtin_prefetch(nextVisitList + *(friendsBounds.first+PREFETCH),1);
			  }
			  #endif /*DO_PREFETCH*/

			  for(int i=0; i<width; i++) {
				 nextVisitList[*friendsBounds.first].data[i] |= curVisit.data[i];
			  }
			  ++friendsBounds.first;
		   }

		   /* reset visitList */
		   for(int i=0; i<width; i++) {
			  visitList[curPerson].data[i] = BitBaseOp<uint64_t>::zero();
		   }
		}

		/* udpate seen and dist */
		for (Query4::PersonId curPerson = 0; curPerson<limit; ++curPerson) {
		   for(unsigned i=0; i<width; i++) {
			  const uint64_t nextVisit = nextVisitList[curPerson].data[i];
			  if(BitBaseOp<uint64_t>::notZero(nextVisit)) {
				 const uint64_t newVisits = BitBaseOp<uint64_t>::andNot(nextVisit, seen[curPerson].data[i]);
				 nextVisitList[curPerson].data[i] = newVisits;
				 if(BitBaseOp<uint64_t>::notZero(newVisits)) {
					 updateDist(i, newVisits, dist[curPerson], nextDistance);
					 seen[curPerson].data[i] |= newVisits;
					 batchDist.updateDiscovered(newVisits, i);
				 }
			  }
		   }
		}

		batchDist.finalize();
    }

	#else /*SORTED_NEIGHBOR_PROCESSING*/

    void runBatchRound(const Query4::PersonSubgraph& subgraph, const Query4::PersonId startPerson, const Query4::PersonId limit,
        Bitset* visitList, Bitset* nextVisitList, Bitset* __restrict__ seen, int** dist,
		Query4::BatchDistance<uint64_t, width>& batchDist, const Bitset processQuery, uint32_t nextDistance) {

        for (Query4::PersonId curPerson = startPerson; curPerson<limit; ++curPerson) {
            Bitset validVisit = createVisitList(visitList[curPerson], processQuery);

            /* Skip persons with empty visit list */
            if(validVisit.count() == 0) { continue; }

            const auto& curFriends=*subgraph.retrieve(curPerson);

            auto friendsBounds = curFriends.bounds();

			#ifdef DO_PREFETCH
			__builtin_prefetch(seen + *(friendsBounds.first+1),0);
			__builtin_prefetch(nextVisitList + *(friendsBounds.first+1),1);
			__builtin_prefetch(seen + *(friendsBounds.first+2),0);
			__builtin_prefetch(nextVisitList + *(friendsBounds.first+2),1);
			#endif /*DO_PREFETCH*/

			/* for each neighbor, update visitNext, seen and dist */
            while(friendsBounds.first != friendsBounds.second) {

				#ifdef DO_PREFETCH
				if(friendsBounds.first+3 < friendsBounds.second) {
				   __builtin_prefetch(seen + *(friendsBounds.first+3),0);
				   __builtin_prefetch(nextVisitList + *(friendsBounds.first+3),1);
				}
				#endif /*DO_PREFETCH*/

                for(unsigned i=0; i<width; i++) {
                    uint64_t newVisits = BitBaseOp<uint64_t>::andNot(validVisit.data[i], seen[*friendsBounds.first].data[i]);
                    if(BitBaseOp<uint64_t>::notZero(newVisits)) {
                        updateDist(i, newVisits, dist[*friendsBounds.first], nextDistance);
                        seen[*friendsBounds.first].data[i] |= validVisit.data[i];
                        nextVisitList[*friendsBounds.first].data[i] |= newVisits;

                        // Normal case until uint64_t
                        batchDist.updateDiscovered(newVisits, i);
                    }
                }
                ++friendsBounds.first;
            }
        }
        batchDist.finalize();
    }

	#endif /*SORTED_NEIGHBOR_PROCESSING*/

};

struct VirtualRunner {
	VirtualRunner(){}
	virtual ~VirtualRunner(){}
	virtual void run(const size_t& numThreads, Workers & workers, const Query4::PersonSubgraph& subgraph, const std::string& output) = 0;
};

template<uint64_t width>
struct CustomRunner : public VirtualRunner{
	CustomRunner(){}
    ~CustomRunner() { }

    /* initialize task and execute each task in one thread */
    /* don't understand the details of the Scheduler, Executor, etc */
    void run(const size_t& numThreads, Workers & workers, const Query4::PersonSubgraph& subgraph, const std::string& output){
    	   TaskGroup tasks;
    	   for(size_t i = 0; i < numThreads; i++){
    		   BFSTask<width> bfsTask(i, subgraph, output);
    		   tasks.schedule(LambdaRunner::createLambdaTask(bfsTask));
    	   }

    	   Scheduler scheduler;
    	   scheduler.schedule(tasks.close());
    	   scheduler.setCloseOnEmpty();

    	   workers.assist(scheduler);

    	   // Always run one executor on the main thread
    	   Executor executor(scheduler,0, false);
    	   executor.run();

    	   scheduler.waitAllFinished();
    }
};

#endif /* EXTRA_HPP_ */
