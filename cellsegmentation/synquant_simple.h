#ifndef SYNQUANTSIMPLE_H
#define SYNQUANTSIMPLE_H

#include "cellsegment_main.h"
#include "img_basic_proc_declare.h"
#include <algorithm>

#define UNDEFINED (byte)0
#define NOT_OBJECT (byte)1
#define OBJECT (byte)2
#define NOT_USED_AS_N (byte)0
#define USED_AS_N_ONCE (byte)1
#define USED_AS_N_MORE (byte)2

using namespace cv;
using namespace std;
class synQuantSimple
{
public:
    synQuantSimple(Mat *_srcVolume, float _src_var, segParameter &p4segVol, odStatsParameter &p4odStats);
    synQuantSimple(singleCellSeed &seed, odStatsParameter &p4odStats);
    //~synQuantSimple(){} //TODO: delete the pointers
public:
    void processVoxLabel(size_t j);
    void componentTree3d(segParameter p4segVol, odStatsParameter p4odStats);
    void componentTree3d4seed(singleCellSeed &seed, odStatsParameter &p4odStats);
    float zscoreCal(float t0, size_t M/*in*/, size_t N/*nei*/);

    size_t findNode(size_t e);
    size_t mergeNodes(size_t e1,size_t e2);


    void objLabel(size_t minSize ,size_t maxSize); //label objects simply using size constraints
    void objLabel_zscore(float zscore_thres); // TODO: label objects using a zscore threshold
    void objLabel_descending();//label object from the largest zscore one by one;

    void fdr_control();
    // zscore for comparing fg and bg neighbors
    float debiasedFgBgCompare(unsigned debiasMethod);
//    template <typename T> T debiasedFgBgCompare(vector<T> const & fg, vector<T> const & bg, vector<T> const & neglectVals,
//                                                unsigned debiasMethod);
public:
    Mat *zMap, *idMap;
    int cell_num;
    vector<float> zscore_list;
    vector<float> valid_zscore;
    vector<size_t> valid_zscore_idx;
    vector<size_t> intensity_levels;
    float max_exist_zscore;
    size_t maxZ_intensity_level;
protected:
    Mat *srcVolumeUint8;
    float src_var;
    unsigned char* imArray;//point to the same address as srcVolume
    vector<size_t> sortedIndex, parNode;
    vector<size_t> voxSum, voxSumN, areas, areasN;
    vector<byte> usedN;
    vector<byte> outputArray;
    vector<float> diffN;
    size_t width, height, zSlice;
    vector<vector<size_t>> BxCor;
};

#endif // SYNQUANTSIMPLE_H
