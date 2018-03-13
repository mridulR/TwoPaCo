#ifndef _EDGE_RESULT_
#define _EDGE_RESULT_

 #include "junctionapi/junctionapi.h"

namespace TwoPaCo
{

 struct EdgeResult
    {
        uint32_t pieceId;
        std::vector<JunctionPosition> junction;
    };
}
#endif
