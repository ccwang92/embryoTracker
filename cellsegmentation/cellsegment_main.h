#ifndef CELLSEGMENT_MAIN_H
#define CELLSEGMENT_MAIN_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include "img_basic_proc_declare.h"

//#define VOLUME_WH_MAX 100000
using namespace cv;
using namespace std;

struct segParameter {
    union
    {
        // for synQuant
        float min_intensity;
        size_t min_cell_sz;
        size_t max_cell_sz;
        float min_fill;
        size_t max_WHRatio;
        float noise_estimate_ratio;
        float fdr;
        float min_zscore;

        // for seed refine
        size_t min_seed_size;
        int graph_cost_design[2];
        int growConnectInTest, growConnectInRefine;
        int edgeConnect, neiMap;
        int connect4fgGapRemoval;
        int shift_yxz[3];
        bool shrink_flag;
        int shrink_scale_yxz[3];
        int fgBoundaryHandle;// leaveAloneFirst, compete or repeat:
//            if the fg is touching the boundary of the intial foreground, we can choose to
//            (1) leaveAloneFirst: if after gap testing, the region is still large, we
//            enlarge the foreground and re-do our pipeline
//            (2) repeat: enlarge the foreground and re-do pipeline immediately, not
//            waiting for the next gap testing step
//            (3) compete: using the boundary as another region to compete with seed
//            region the split the fg. This method is default when do tracking
    // more possible parameters
//        % If we scaled the data
//        q.minSize = q.minSize/scaling;% cell size
//        q.minSeedSize = q.minSeedSize/scaling; % seed size
//        q.shift = floor(q.shift/scaling);
//        if tracking_flag
//            % should we remove small regions?
//            q.shortestTrack = 0; % the shortest path we want, cells in other tracks will be removed
//            q.removeSamllRegion = false; % though we have threshold, we did not
//            q.fgBoundaryHandle = 'leaveAloneFirst';
//            q.growSeedTimes = 2;
//            q.growSeedinTracking = true;
//            q.multi_frames_flag = false; % we did not consider multiple frames. Otherwise
//            % there may be results dominated by other bright cells
//            q.multiSeedProcess = true; % did we add multiple cells simultaneously or one by one
//            q.splitRegionInTracking = true;

//            q.updateCellsAdjMissingCell = false;% when add missing cell, do we need to update other regions nearby
//            q.sqrtDistance = false; % euclidian distance or squared euclidian distance
//        end
        int gapTestMinMaxRadius[2];
        bool growSeedInTracking;
    };
};

struct odStatsParameter {
    union
    {
        int gap4varTrendEst;
        int gap4fgbgCompare; // default is 0
        int roundNum4fgbgCompare;
        float varAtRatio;
        int fgSignificanceTestWay;
        int minGapWithOtherCell_yxz[3];
        int connectInSeedRefine;
        int gapTestMethod;
        int gapTestSkippedBandWidth;
        float gapTestThreshold;
    };
};
struct singleCellSeed{
    int id;
    vector<size_t> idx_yxz, idx_yxz_cropped;
    vector<int> y, x, z;
    Range crop_range_yxz[3];
    Mat seedMap;
    Mat eigMap2d, eigMap3d;
    Mat score2d, score3d;
    Mat scoreMap;
    Mat gap2dMap, gap3dMap;
    Mat varMap;
    Mat stblizedVarMap;
    Mat volUint8;
    Mat volStblizedFloat;
    Mat idMap; //idComp
    //Mat fgMapGapRemoved; //newIdComp
    Mat fgMap;
    Mat otherIdMap;
    Mat outputIdMap;
    int outCell_num;
    int bestFgThreshold = -1;
};

class cellSegmentMain
{
public:
    cellSegmentMain(unsigned char *data_grayim4d, int _data_type, long buffSize[5]/*(x,y,z,c,t)*/);
    ~cellSegmentMain(){delete data_rows_cols_slices;};
    void cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame);
    void regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
                                               Mat *eigMap3d, Mat *varMap, Mat * stblizedVarMap, vector<int> test_ids);
    void cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d, Mat *data_stbized, Mat *idMap, Mat *eigMap2d,
                                   Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, singleCellSeed &seed, segParameter p4segVol);
    int refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol);
protected:
    string debug_folder;
    string default_name;
    int data_type;
    int * data_rows_cols_slices;
    long time_points = 0;
    long curr_time_point;
    segParameter p4segVol;
    odStatsParameter p4odStats;
    vector<Mat> cell_label_maps;
    vector<Mat> threshold_maps;
    vector<Mat> principalCurv2d;
    vector<Mat> principalCurv3d;
    vector<Mat> varMaps;
    vector<Mat> stblizedVarMaps;
    vector<vector<float>> varTrends;
    vector<vector<float>> stblizedVarTrends;
    vector<float> variances;
    vector<size_t> number_cells;
    vector<vector<size_t>> voxIdxList; // cell voxIdx list

    void init_parameter(){
        string debug_folder = "/home/ccw/Desktop/embryo_res_folder/";
        string default_name = debug_folder + "test.tiff";
        p4segVol.min_intensity = 0.0;
        p4segVol.fdr = .05;
        p4segVol.min_cell_sz =100;
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
        p4segVol.shrink_scale_yxz[0] = 4;
        p4segVol.shrink_scale_yxz[1] = 4;
        p4segVol.shrink_scale_yxz[2] = 4;
        p4segVol.fgBoundaryHandle = LEAVEALONEFIRST;
        p4segVol.gapTestMinMaxRadius[0] = 2;
        p4segVol.gapTestMinMaxRadius[1] = 4;
        p4segVol.growSeedInTracking = false;

        p4odStats.gap4varTrendEst = 2;
        p4odStats.gap4fgbgCompare = 0;
        p4odStats.roundNum4fgbgCompare = 3;
        p4odStats.varAtRatio = 0.95;
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
