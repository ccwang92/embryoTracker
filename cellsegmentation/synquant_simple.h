#ifndef SYNQUANTSIMPLE_H
#define SYNQUANTSIMPLE_H

#include "img_basic_proc.h"
//#include "vol_basic_proc.hpp"
#include "types_define.h"
#include <algorithm>
#include <set>

#define UNDEFINED (byte)0
#define NOT_OBJECT (byte)1
#define OBJECT (byte)2
#define NOT_USED_AS_N (byte)0
#define USED_AS_N_ONCE (byte)1
#define USED_AS_N_MORE (byte)2


class synQuantSimple
{
public:
    synQuantSimple(cv::Mat *_srcVolume, float _src_var, segParameter &p4segVol, odStatsParameter &p4odStats);
    synQuantSimple(singleCellSeed &seed);
    ~synQuantSimple(){
        if(idMap){
            delete idMap;
        }
        if(zMap) {
            delete zMap;
        }
    }
public:
    void processVoxLabel(std::size_t j);
    void componentTree3d(segParameter p4segVol, odStatsParameter p4odStats);
    void cellTerritoryExtractFromSeed(singleCellSeed &seed, odStatsParameter &p4odStats, size_t minSize = 100);
    float zscoreCal(float t0, std::size_t M/*in*/, std::size_t N/*nei*/); // directly use non_overlap_gaussian

    std::size_t findNode(std::size_t e);
    std::size_t mergeNodes(std::size_t e1,std::size_t e2);


    void objLabel(std::size_t minSize ,std::size_t maxSize); //label objects simply using size constraints
    void objLabel_zscore(float zscore_thres); // TODO: label objects using a zscore threshold
    void objLabel_descending();//label object from the largest zscore one by one;

    void fdr_control();
    // zscore for comparing fg and bg neighbors using the boundary pixels among fg and bg
    float debiasedFgBgBandCompare(cv::Mat *cur_reg, cv::Mat *validNei, singleCellSeed *seed, odStatsParameter p4odStats);
//    void refineCellTerritoryWithSeedRegion(singleCellSeed &seed, segParameter &p4segVol);
//    void refineCellsTerritoriesWithSeedRegions(singleCellSeed &seed, segParameter &p4segVol);
//    void cellShrinkTest(singleCellSeed &seed, segParameter &p4segVol);
//    void fgGapRemoval(singleCellSeed &seed, segParameter &p4segVol);
//    void gapBasedRegionSegment(singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats);
//    void gapTest2SplitCellTerritory(cv::Mat* seeds_Map /*CV_32S*/, int n, singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats);
//    void removeOtherSeedsInfgMap(singleCellSeed &seed, segParameter &p4segVol);
public:
    cv::Mat *zMap, *idMap, fgMap, fgMapGapRemoved;
    int cell_num;
    std::vector<float> zscore_list;
    std::vector<float> valid_zscore;
    std::vector<std::size_t> valid_zscore_idx;
    std::vector<std::size_t> intensity_levels;
    float max_exist_zscore;
    std::size_t maxZ_intensity_level;
protected:
    cv::Mat *srcVolumeUint8;//point to the same address as srcVolume
    float src_var;
    unsigned char* imArray;//point to the same address as srcVolume
    std::vector<std::size_t> sortedIndex, parNode;
    std::vector<std::size_t> voxSum, voxSumN, areas, areasN;
    std::vector<byte> usedN;
    std::vector<byte> outputArray;
    std::vector<float> diffN;
    std::size_t width, height, zSlice;
    std::vector<std::vector<std::size_t>> BxCor;
};

#endif // SYNQUANTSIMPLE_H
