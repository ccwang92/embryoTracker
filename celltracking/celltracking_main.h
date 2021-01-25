#ifndef CELLTRACKINGMAIN_H
#define CELLTRACKINGMAIN_H
#include "../cellsegmentation/cellsegment_main.h"


class cellTrackingMain
{
public:
    cellTrackingMain(cellSegmentMain &cellSegment);
    ~cellTrackingMain(){};
protected:
    cellCensus movieInfo;
};

#endif // CELLTRACKINGMAIN_H
