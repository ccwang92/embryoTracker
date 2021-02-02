#ifndef TYPES_DEFINE_H
#define TYPES_DEFINE_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <string>
#include <vector>

struct segParameter {
    /** *******  for synQuant ***************/
    float min_intensity;
    std::size_t min_cell_sz;
    std::size_t max_cell_sz;
    float min_fill;
    std::size_t max_WHRatio;
    float noise_estiMate_ratio;
    float fdr;
    float min_zscore;

    /** ******* for seed refine ***************/
    std::size_t min_seed_size;
    int graph_cost_design[2];
    int growConnectInTest, growConnectInRefine;
    int edgeConnect, neiMap;
    int connect4fgGapRemoval;
    int shift_yxz[3];
    bool shrink_flag;
    int shrink_scale_yxz[3];
    int gapTestMinMaxRadius[2];
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
    /** ******* for tracking ***************/
    long validTrackLength; // the shortest path we want, cells in other tracks will be removed
    bool removeSamllRegion; // though we have threshold, we did not
    int growSeedTimes;
    bool growSeedInTracking;
    bool multi_frames_flag; // we did not consider multiple frames. Otherwise there may be results dominated by other bright cells
    bool multiSeedProcess;  // did we add multiple cells simultaneously or one by one
    bool splitRegionInTracking;

    bool updateCellsAdjMissingCell; //  when add missing cell, do we need to update other regions nearby
    bool sqrtDistance; // euclidian distance or squared euclidian distance


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

struct trackParameter{
    bool cycle_track;//circulation framework to solve tracking problem
    int stdCombined;//we are cal difference between two locations (transitCostInitial_cell line 58)
    int maxDistXYZ[3]; // max moving distance for a cell in adjacent frames
    int k; // maximum number of jump allowed
    int transitionFactor;// the weight of transition cost
    int validPre;// check 2 previous point to get the mean
    int validPost;// check 2 post point to get the mean
    long cellNum;
    bool timeJump;//Consier the jump x5 = x4+v*(t5-t4) , t5-t4 may not be 1
    float initEnter;// initial enter/exit cost, force every node link to its neighbor
    float realEnter;// 12.3546 is from chi2inv(1-0.01/particleNum) / 2
    float c_en;// cost of appearance and disappearance in the scene
    float c_ex;//
    float observationCost;// make sure detections are all included
    float jumpCost[3];// how much we should punish jump frames
    int varEstMethod; // median and independent
    int costCalMethod; //chi1Square:use 1df chi-squre, fisher: 2df chi-square,zscore: use z-score
    int validtrackLength4var;// tracks with smaller length will not be used to cal variance
    float truncatedGaussian ;// remove 20// of extreme edges in variance estiamtion (truncated gaussian)
    int varPrior;// use gamma distribution as variance prior, use longest 100 tracks to estimate gamma parameters
    int priorType; // prior distribution: gamma/weiAvg/Gauss/ scaled inverse chi-square  (sInvX2)
    bool directionWise;// cal pvalue direction-wise and use fisher method (or chi-square) to combine
    float dependencyFactor;// the particles are not independent, based on their correlation, we add a correction factor for variance estimation.
    int splitMergeHandle;// noJumpAll: all linking has no jump; noJump: no jump for spliting and merge; none: consider jump
    int maxIter;// maximum iteration number
    bool removeSingleRegion;// if a cell does not participate any trace, remove it.
    bool detect_missing_head_tail;// add missing cell, do we consider head/tail node
    bool applyDrift2allCoordinate;// correct drifting and change all y-, x- and z- coordinate

    // parameters invovle segmentation
    bool blindmergeNewCell;// for a new detected region, blindly merge it to existing one if there is touch
    bool simple_merge;// for two regions, if we decide to merge them, simply add their pixels together
    bool par_kid_consistency_check;// when split a region in regionRefresh.m, check first if its two parents and two kids are consistent
    bool reSplitTest;// for a new detected cell, if it is adjacent to another cell, re-split these two using their union voxels
    bool stableNodeTest;//
    bool useOldArcProps;//
    bool considerBrokenCellOnly;// for linking allowing split/merge, does not consider nodes that has another good linkage already
    bool addCellMissingPart;// if a cell missed part, we want to detect it, otherwise we can remove seeds that highly overlapped with an existing cell
    bool splitMergeCost;// if cost of a+b->c is 20, then cost of a->c and b->c are set to 10 if true; otherwise both 20
};

// structures to save the detection infor
struct nodeRelation{
    size_t node_id;
    float dist_p2k, dist_k2p; // current to neighbor in next frames or in counter direction
    long overlap_size;
    float link_cost;
};
//struct directFamily{ // for each node, it has at most two parents or two kids
//    simpleCell parents[2];
//    int parent_num = 0;
//    simpleCell kids[2];
//    int kid_num = 0;
//};
struct nodeInfo{
    size_t node_id;
    //directFamily family_members; // neighboring relationship, at most two kids or parents
    size_t parents[2];
    int parent_num = 0;
    size_t kids[2];
    int kid_num = 0;
    size_t nodeId2trackId, nodeLocInTrack;
    std::vector<nodeRelation> neighbors; // candidate kids
};
//struct nodeInfoInTrack{
//    size_t nodeId2trackId, nodeLocInTrack;
//};
struct singleCellCensus{
    nodeInfo node_info; //
    float xCoord, yCoord, zCoord;
    int frames;
    int labelInMap; //read id in label_maps
    std::vector<size_t> voxIdx;
    std::vector<int> vox_x, vox_y, vox_z;
//    std::vector<nodeRelation> neighbors; // candidate kids
//    nodeInfoInTrack particle2track; //particle2track in matlab
    //directFamily family_members; // neighboring relationship, at most two kids or parents
};


struct allCellsCensus{
    std::vector<float> xCoord, yCoord, zCoord;
    std::vector<int> frames;
    std::vector<int> labelInMap; //read id in label_maps
    std::vector<std::vector<size_t>> voxIdx;
    std::vector<std::vector<int>> vox_x, vox_y, vox_z;
//    std::vector<std::vector<nodeRelation>> neighbors; // candidate kids
//    std::vector<singleCellCensus> cells;
    std::vector<nodeInfo> nodes;
    std::vector<std::vector<size_t>> tracks; //a set of node_ids
//    std::vector<nodeInfoInTrack> particle2track; //particle2track in matlab
//    std::vector<directFamily> parents, kids; // neighboring relationship, at most two kids or parents
};

#endif // TYPES_DEFINE_H
