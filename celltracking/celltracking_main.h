#ifndef CELLTRACKINGMAIN_H
#define CELLTRACKINGMAIN_H
#include "../cellsegmentation/cellsegment_main.h"


class cellTrackingMain
{
public:
    cellTrackingMain(cellSegmentMain &cellSegment);
    ~cellTrackingMain(){};
private:
    void cellInfoAccumuate(cellSegmentMain &cellSegment);
    void initTransitionCost(cellSegmentMain &cellSegment);

    // distance calculation
    float voxelwise_avg_distance(size_t cell_curr, size_t cell_nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(size_t joint_cells_curr[], size_t joint_cells_nei[], float &c2n, float &n2c);
    float voxelwise_avg_distance(vector<size_t> &joint_cells_curr, vector<size_t> &joint_cells_nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(size_t cell_curr, vector<size_t> &joint_cells_nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(vector<size_t> &joint_cells_curr, size_t cell_nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(combinedCellsCensus &curr, combinedCellsCensus &nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(size_t curr, combinedCellsCensus &nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(combinedCellsCensus &curr, size_t nei, float &c2n, float &n2c);
    float voxelwise_avg_distance(vector<size_t> &curr_voxIdx, vector<int> &curr_vox_x, vector<int> &curr_vox_y,
                                 vector<int> &curr_vox_z, vector<int> &curr_range_xyz, vector<int> &curr_start_coord_xyz, int curr_frame,
                                 vector<size_t> &nei_voxIdx, vector<int> &nei_vox_x, vector<int> &nei_vox_y,
                                 vector<int> &nei_vox_z, vector<int> &nei_range_xyz, vector<int> &nei_start_coord_xyz, int nei_frame,
                                 float &c2n, float &n2c);
    float voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
                                 vector<size_t> &nei_voxIdx, int nei_frame,
                                 MatSize data3d_sz, float &c2n, float &n2c);
    float voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
                                 vector<size_t> &nei_voxIdx, int nei_frame,
                                 int data3d_sz[3], float &c2n, float &n2c);
    void updatePreNeighborInfo();
    size_t cellOverlapSize(size_t c0, size_t c1, cellSegmentMain &cellSegment);
    void reCalculateCellsDistances(vector<float> &nn_dist);
    void calCell2neighborDistance(vector<float> &nn_dist);
    float distance2cost(float distance, float alpha, float beta, float punish);
    void extractNeighborIds(vector<Mat> &cell_label_maps, size_t node_idx, vector<size_t> & nei_idxs);
    void extractPreNeighborIds(vector<Mat> &cell_label_maps, size_t cell_idx, vector<size_t> &nei_idxs);
    void mccTracker_one2one();
    void mccTracker_splitMerge(vector<splitMergeNodeInfo> &split_merge_node_info);
    void track2parentKid();
    void refreshTracks(); //remove the empty tracks
    void mergeOvTracks(); //
    void node2trackUpt();
    void updateJumpCost();
    // functions to update cost given new gamma fitting results
    void driftCorrection();
    void updateGammaParam(vector<float> &nn_dist);
    void updateArcCost();
    void getArcCostOne2OneTrack(size_t track_id, vector<float> &arc_costs);
    void stableSegmentFixed();
    void movieInfoUpate();

    bool isValid(size_t node_idx, cellSegmentMain *cellSegment = nullptr);
    bool isBestNeighbor(size_t n1, size_t n2, float &cost);
    bool findBestNeighbor(size_t n1, size_t &out_n, float &cost, int target_frame);
    void infinitifyCellRelation(size_t n1, size_t n2);
    void nullifyCellOrNode(size_t node_idx, cellSegmentMain *cellSegment = nullptr); //
    void nullifyCellOrNode(size_t node_idx[], cellSegmentMain *cellSegment = nullptr); //
    void nullifyCellOrNode(vector<size_t> node_idx, cellSegmentMain *cellSegment = nullptr); //
    void appendNewCellOrNode(cellSegmentMain &cellSegment, simpleNodeInfo &newCell, size_t append_loc_idx);
private: // remaining for split/merge module
    void split_merge_module(cellSegmentMain &cellSegment);

    void handleMergeSplitRegions();
    void regionRefresh(cellSegmentMain &cellSegment, vector<simpleNodeInfo> &newCells, vector<size_t> &uptCell_idxs);

    void detectPeerRegions(vector<splitMergeNodeInfo> &split_merge_node_info,
                           unordered_map<long, long> &node_id2split_merge_node_id);
    void combineCellsIntoOneRegion(vector<size_t> &cell_idxes, combinedCellsCensus &out_region_info);
    float bestPeerCandidate(size_t node_id, vector<size_t> &out_bestPeer, bool parent_flag);
    void peerRegionVerify(size_t node_id, float cost_good2go, bool parents_test,
                          vector<splitMergeNodeInfo> &split_merge_node_info,
                          unordered_map<long, long> &node_id2split_merge_node_id);
    size_t sizeCumulate(size_t curr_cell, size_t familiy_members[2]);


    bool exist_in_pairs(vector<pair<size_t[2], int>> &pairs, size_t id);
    int parentsKidsConsistency(size_t node_id);
    int handleInconsistentParentKid(cellSegmentMain &cellSegment, size_t node_id);

    // TODO: it seems we have already avoid conflict decisions
    void conflict_decision_handle(vector<tuple<size_t, int, float>> &merged_split_peers);
    bool seedsRefine_intensity(cellSegmentMain &cellSegment, vector<size_t> &root_idx, int root_frame,
                               vector<vector<size_t>> &seeds_idx, int seed_frame,
                               vector<vector<size_t>> &ref_seeds_idx);
    bool seedsRefine_gap(Mat1b &possibleGaps, vector<vector<size_t>> &seeds_idx, Mat1i &outLabelMap);


    bool binary_seedsMap_create(Mat1b &fgMap, Mat1b *possbileGap3d, Mat1b *possbileGap2d,
                                vector<vector<size_t>> seeds_idx, Mat1i &seeds_map, size_t minSz);
    bool bisectRegion_gapGuided(cellSegmentMain &cellSegment, size_t reg2split_idx,
                                vector<size_t> &reg2split, int reg2split_frame,
                                vector<vector<size_t>> &reg4seeds, bool usePriorGapMap, vector<vector<size_t>> &splitRegs);
    bool bisectRegion_bruteforce(cellSegmentMain &cellSegment, size_t reg2split_idx,
                                 vector<size_t> &reg2split, int reg2split_frame,
                                 vector<vector<size_t>> &reg4seeds, int reg4seeds_frame,
                                 bool usePriorGapMap, vector<vector<size_t>> &splitRegs);
    bool bisectValidTest(cellSegmentMain &cellSegment, size_t reg2split_idx, vector<size_t> reg2split,
                         int reg2split_frame, vector<vector<size_t>> reg4seeds, int reg4seeds_frame,
                         bool gapBasedSplit, bool usePriorGapMap, vector<vector<size_t>> &splitRegs,
                         float *reg4seeds2splitRes_costs);
    int regionSplitMergeJudge(cellSegmentMain &cellSegment, size_t curr_node_id, bool one2multiple_flag, float &pvalue);
    bool testCellsInOneTrackAdjacentOrNot(cellSegmentMain &cellSegment, vector<unordered_set<size_t>> left_or_right_cells);
    bool mergeValidTest(size_t curr_node_id, size_t seedRegs4split[2]);


    bool handleNonSplitReg_link2oneSeed(size_t node_idx, size_t seeds[2]);
    bool separateRegion(cellSegmentMain &cellSegment, size_t node_idx, size_t seeds[2], bool oneSeedIsGood, vector<simpleNodeInfo> &newCells);

    bool mergedRegionGrow(cellSegmentMain &cellSegment, size_t seeds[2], vector<simpleNodeInfo> &newCells);
private:    // TODO: missing cell module
    void missing_cell_module(cellSegmentMain &cellSegment);

private:
    void movieInfo_update(cellSegmentMain &cellSegment, vector<simpleNodeInfo> &newCells, vector<size_t> &uptCell_idxs);
private:
    allCellsCensus movieInfo;
    trackParameter p4tracking;
    vector<size_t> cumulative_cell_nums;
    vector<Mat1b> validGapMaps;
    //friend class cellSegmentMain;
public:
    void init_parameter(segParameter &p4seg, long cell_num){
        p4tracking.cycle_track = true; // true: circulation framework to solve tracking problem
        p4tracking.stdCombined = 1; // we are cal difference between two locations (transitCostInitial_cell line 58)
        p4tracking.maxDistXYZ[0] = 25; // max moving distance for a cell in adjacent frames
        p4tracking.maxDistXYZ[1] = 25; // max moving distance for a cell in adjacent frames
        p4tracking.maxDistXYZ[2] = 5; // max moving distance for a cell in adjacent frames
        p4tracking.k = 3; // maximum number of jump allowed
        p4tracking.cellNum = cell_num; // cell number
        p4tracking.transitionFactor = 1;// the weight of transition cost
        p4tracking.validPre = 4; // check 2 previous point to get the mean
        p4tracking.validPost = 4; // check 2 post point to get the mean
        //p4tracking.maxEdgeNum = 4; at most check maxEdgeNum following or previous edges
        p4tracking.timeJump = false; // Consier the jump x5 = x4+v*(t5-t4) , t5-t4 may not be 1
        p4tracking.initEnter = 100; // initial enter/exit cost, force every node link to its neighbor
        p4tracking.realEnter = chi2inv(1-0.01/cell_num, 1)/2; // 12.3546 is from chi2inv(1-0.01/particleNum) / 2
        p4tracking.c_en = p4tracking.realEnter;// cost of appearance and disappearance in the scene
        p4tracking.c_ex = p4tracking.c_en;
        p4tracking.observationCost = -(p4tracking.c_en+p4tracking.c_ex); // make sure detections are all included
        p4tracking.jumpCost.resize(p4tracking.k);
        p4tracking.jumpCost[0] = -1;
        p4tracking.jumpCost[1] = -1;
        p4tracking.jumpCost[2] = -1;
        p4tracking.varEstMethod = 1; // 0 for median and 1 for independent
        p4tracking.costCalMethod= ccm_CHI1SQUARE; // chi1Square:use 1df chi-squre, fisher: 2df chi-square,zscore: use z-score
        p4tracking.validtrackLength4var = 5;// tracks with smaller length will not be used to cal variance
        p4tracking.truncatedGaussian = 0.2; // remove 20// of extreme edges in variance estiamtion (truncated gaussian)
        p4tracking.varPrior = 100;// use gamma distribution as variance prior, use longest 100 tracks to estimate gamma parameters
        p4tracking.priorType = pri_GAUSS; // prior distribution: gamma/weiAvg/Gauss/ scaled inverse chi-square  (sInvX2)
        p4tracking.directionWise = false;// cal pvalue direction-wise and use fisher method (or chi-square) to combine
        p4tracking.dependencyFactor = 1.4;// the particles are not independent, based on their correlation, we add a correction factor for variance estimation.
        p4tracking.splitMergeHandle = NOJUMPALL;// noJumpAll: all linking has no jump; noJump: no jump for spliting and merge; none: consider jump
        p4tracking.maxIter = 7; // maximum iteration number
        p4tracking.removeSingleRegion = true; // if a cell does not participate any trace, remove it.
        p4tracking.detect_missing_head_tail = true; // add missing cell, do we consider head/tail node
        p4tracking.applyDrift2allCoordinate = false; // correct drifting and change all y-, x- and z- coordinate


        p4tracking.blindmergeNewCell = false; // for a new detected region, blindly merge it to existing one if there is touch
        p4tracking.simple_merge = true;// for two regions, if we decide to merge them, simply add their pixels together
        p4tracking.par_kid_consistency_check = true; // when split a region in regionRefresh.m, check first if its two parents and two kids are consistent
        p4tracking.reSplitTest = true; // for a new detected cell, if it is adjacent to another cell, re-split these two using their union voxels
        p4tracking.stableNodeTest =true;
        p4tracking.useOldArcProps = false;
        p4tracking.considerBrokenCellOnly = true; // for linking allowing split/merge, does not consider nodes that has another good linkage already
        p4tracking.addCellMissingPart = false; // if a cell missed part, we want to detect it, otherwise we can remove seeds that highly overlapped with an existing cell
        p4tracking.splitMergeCost = true;// if cost of a+b->c is 20, then cost of a->c and b->c are set to 10 if true; otherwise both 20
        p4tracking.min_stable_node_cluster_sz = 5;

        ////update cost of p4seg
        p4seg.validTrackLength = 0; // the shortest path we want, cells in other tracks will be removed
        p4seg.removeSamllRegion = false; // though we have threshold, we did not
        p4seg.fgBoundaryHandle = LEAVEALONEFIRST;
        p4seg.growSeedTimes = 2;
        p4seg.growSeedInTracking = true;
        p4seg.multi_frames_flag = false; // we did not consider multiple frames. Otherwise
        // there may be results dominated by other bright cells
        p4seg.multiSeedProcess = true; // did we add multiple cells simultaneously or one by one
        p4seg.splitRegionInTracking = true;

        p4seg.updateCellsAdjMissingCell = false;// when add missing cell, do we need to update other regions nearby
        p4seg.sqrtDistance = false; // euclidian distance or squared euclidian distance

    }
};

#endif // CELLTRACKINGMAIN_H
