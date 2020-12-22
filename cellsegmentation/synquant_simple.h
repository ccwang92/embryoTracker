#ifndef SYNQUANTSIMPLE_H
#define SYNQUANTSIMPLE_H
#include "img_basic_proc_declare.h"
#include "cellsegment_main.h"

#define UNDEFINED (byte)0
#define NOT_OBJECT (byte)1
#define OBJECT (byte)2
#define NOT_USED_AS_N (byte)0
#define USED_AS_N_ONCE (byte)1
#define USED_AS_N_MORE (byte)2

class synQuantSimple{
public:
    synQuantSimple(Mat *_srcVolume, segParameter p4segVol, odStatsParameter p4odStats);

    void componentTree3d(segParameter p4segVol, odStatsParameter p4odStats);
    void objLabel(size_t minSize ,size_t maxSize);
    size_t findNode(size_t e);
    size_t mergeNodes(size_t e1,size_t e2);
public:
    Mat zMap, idMap;
    vector<float> zscore_list;
    vector<size_t> intensity_levels;
    float max_exist_zscore;
    size_t maxZ_intensity_level;
protected:
    Mat *srcVolumeUint8;
    vector<size_t> sortedIndex, parNode;
    vector<size_t> voxSum, voxSumN, areas, areasN;
    vector<byte> usedN;
    vector<byte> outputArray;
    vector<float> diffN;
    size_t width, height, zSlice;
    vector<vector<size_t>> BxCor;
};

#endif // SYNQUANTSIMPLE_H
