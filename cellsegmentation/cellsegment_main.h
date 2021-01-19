#ifndef CELLSEGMENT_MAIN_H
#define CELLSEGMENT_MAIN_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include "synquant_simple.h"
//#include "vol_basic_proc.hpp"
//#define VOLUME_WH_MAX 100000

class cellSegmentMain
{
public:
    cellSegmentMain(void *data_grayim4d, int _data_type, long buffSize[5]/*(x,y,z,c,t)*/);
    ~cellSegmentMain(){delete data_rows_cols_slices;};
    void cellSegmentSingleFrame(cv::Mat *data_grayim3d, std::size_t curr_frame);
    void regionWiseAnalysis4d(cv::Mat *data_grayim3d, cv::Mat *dataVolFloat, cv::Mat * volStblizedFloat,
                              cv::Mat *idMap /*int*/, int seed_num, int frame, std::vector<int> test_ids);
    void cropSeed(int seed_id, std::vector<std::size_t> idx_yxz, cv::Mat *data_grayim3d, cv::Mat *data_stbized,
                  cv::Mat *idMap, int frame, singleCellSeed &seed, segParameter p4segVol);
    // the key function to get cell territory
    int refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol);
    // functions called in refineSeed2Region
    void refineCellTerritoryWithSeedRegion(synQuantSimple& cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol);
    void refineCellsTerritoriesWithSeedRegions(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol);
    void cellShrinkTest(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol);
    void fgGapRemoval(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol);
    void gapBasedRegionSegment(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats);
    void gapTest2SplitCellTerritory(synQuantSimple &cellSegFromSynQuant, cv::Mat* seeds_Map /*CV_32S*/, int n, singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats);
    void removeOtherSeedsInfgMap(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol);
protected:
    string debug_folder;
    string default_name;
    int data_type;
    int * data_rows_cols_slices;
    long time_points = 0;
    long curr_time_point;
    segParameter p4segVol;
    odStatsParameter p4odStats;
    std::vector<cv::Mat> cell_label_maps;
    std::vector<cv::Mat> threshold_maps;
    std::vector<cv::Mat> principalCurv2d;
    std::vector<cv::Mat> principalCurv3d;
    std::vector<cv::Mat> varMaps;
    std::vector<cv::Mat> stblizedVarMaps;
    std::vector<std::vector<float>> varTrends;
    std::vector<std::vector<float>> stblizedVarTrends;
    std::vector<float> variances;
    std::vector<std::size_t> number_cells;
    std::vector<std::vector<std::size_t>> voxIdxList; // cell voxIdx list
public:
    void init_parameter(){
        debug_folder = "/home/ccw/Desktop/embryo_res_folder/";
        default_name = debug_folder + "test.tiff";

        p4segVol.min_intensity = 0.0;
        p4segVol.fdr = .05;
        p4segVol.min_cell_sz = 100;
        p4segVol.max_cell_sz = 3000;
        p4segVol.min_fill = 0.0001;
        p4segVol.max_WHRatio = 100;
        p4segVol.min_seed_size = 10;
        p4segVol.graph_cost_design[0] = ARITHMETIC_AVERAGE; //default 1, GEOMETRIC_AVERAGE = 2;
        p4segVol.graph_cost_design[1] = 2;
        p4segVol.growConnectInTest = 4;
        p4segVol.growConnectInRefine = 6;
        p4segVol.edgeConnect = 48;
        p4segVol.neiMap = 26;
        p4segVol.connect4fgGapRemoval = 26;
        p4segVol.shift_yxz[0] = 20;
        p4segVol.shift_yxz[1] = 20;
        p4segVol.shift_yxz[2] = 4;
        p4segVol.shrink_flag = true;
        p4segVol.shrink_scale_yxz[0] = 2;
        p4segVol.shrink_scale_yxz[1] = 2;
        p4segVol.shrink_scale_yxz[2] = 1;
        p4segVol.fgBoundaryHandle = LEAVEALONEFIRST;
        p4segVol.gapTestMinMaxRadius[0] = 2;
        p4segVol.gapTestMinMaxRadius[1] = 4;
        p4segVol.growSeedInTracking = false;

        p4odStats.gap4varTrendEst = 2;
        p4odStats.gap4fgbgCompare = 0; // the fg and bg should be adjacent
        p4odStats.roundNum4fgbgCompare = 3; // fg/bg both are width-3 band
        p4odStats.varAtRatio = 0.80;
        p4odStats.fgSignificanceTestWay = KSEC;
        p4odStats.minGapWithOtherCell_yxz[0] = 3;
        p4odStats.minGapWithOtherCell_yxz[1] = 3;
        p4odStats.minGapWithOtherCell_yxz[2] = 1;
        p4odStats.connectInSeedRefine = 6;
        p4odStats.gapTestMethod = GAP_LOCALORDERSTATS;
        p4odStats.gapTestSkippedBandWidth = 2;
        p4odStats.gapTestThreshold = 0.01;
    }
    void reset_shift(){
        p4segVol.shift_yxz[0] = 20;
        p4segVol.shift_yxz[1] = 20;
        p4segVol.shift_yxz[2] = 4;
    }
};

#endif // CELLSEGMENT_MAIN_H
