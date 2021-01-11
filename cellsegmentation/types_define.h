#ifndef TYPES_DEFINE_H
#define TYPES_DEFINE_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <string>
#include <vector>

struct segParameter {
    // for synQuant
    float min_intensity;
    std::size_t min_cell_sz;
    std::size_t max_cell_sz;
    float min_fill;
    std::size_t max_WHRatio;
    float noise_estiMate_ratio;
    float fdr;
    float min_zscore;

    // for seed refine
    std::size_t min_seed_size;
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

struct odStatsParameter {
//    union
//    {
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
//    };
};
struct singleCellSeed{
    int id;
    std::vector<std::size_t> idx_yxz, idx_yxz_cropped;
    std::vector<int> y, x, z;
    cv::Range crop_range_yxz[3];
    cv::Mat seedMap;
    cv::Mat eigMap2d, eigMap3d;
    cv::Mat score2d, score3d;
    cv::Mat scoreMap;
    cv::Mat gap2dMap, gap3dMap;
    cv::Mat varMap;
    cv::Mat stblizedVarMap;
    cv::Mat volUint8;
    cv::Mat volStblizedFloat;
    cv::Mat idMap; //idComp
    //cv::Mat fgMapGapRemoved; //newIdComp
    cv::Mat validSearchAreaMap;
    cv::Mat otherIdMap;
    cv::Mat outputIdMap;
    int outCell_num;
    int bestFgThreshold = -1;
};

#endif // TYPES_DEFINE_H
