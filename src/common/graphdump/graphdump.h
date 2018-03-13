#ifndef EXT_GRAPH_DUMP_H
#define EXT_GRAPH_DUMP_H

#include <tbb/concurrent_queue.h>
#include "junctionapi/junctionapi.h"

void run_graph_dump(int argc, char * argv[], std::string file_name, tbb::concurrent_queue<TwoPaCo::JunctionPosition> *queue, std::atomic<bool> * done);

#endif
