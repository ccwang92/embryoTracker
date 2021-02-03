#ifndef CELLTRACKINGMAIN_H
#define CELLTRACKINGMAIN_H
#include "../cellsegmentation/cellsegment_main.h"


class cellTrackingMain
{
public:
    cellTrackingMain(cellSegmentMain &cellSegment);
    void cellInfoAccumuate(cellSegmentMain &cellSegment);
    void initTransitionCost(cellSegmentMain &cellSegment);
    float voxelwise_avg_distance(size_t cell_curr, size_t cell_nei, float &c2n, float &n2c);
    void extractNeighborIds(vector<Mat> &cell_label_maps, size_t node_idx, vector<size_t> & nei_idxs);
    void cellInfo2graph();
    ~cellTrackingMain(){};

private:
    allCellsCensus movieInfo;
    trackParameter p4tracking;
    vector<size_t> cumulative_cell_nums;
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
        //p4tracking.realEnter = chi2inv(1-0.01/cell_num, 1)/2; // 12.3546 is from chi2inv(1-0.01/particleNum) / 2
        p4tracking.c_en = p4tracking.realEnter;// cost of appearance and disappearance in the scene
        p4tracking.c_ex = p4tracking.c_en;
        p4tracking.observationCost = -(p4tracking.c_en+p4tracking.c_ex); // make sure detections are all included
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
