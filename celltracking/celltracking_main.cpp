#include "celltracking_main.h"

cellTrackingMain::cellTrackingMain(cellSegmentMain &cellSegment)
{
    /////////////////////////////////////////////////
    //   step 1. initialize cell info              //
    /////////////////////////////////////////////////
    cellInfoAccumuate(cellSegment);
    init_parameter(cellSegment.p4segVol, movieInfo.frames.size());
    /////////////////////////////////////////////////
    //   step 2. try a naive tracking first        //
    /////////////////////////////////////////////////
    initTransitionCost(cellSegment);
    cellInfo2graph();
    //circulationBasedTracker();
    /////////////////////////////////////////////////
    //   step 3. main loop for cell tracking       //
    /////////////////////////////////////////////////

}

//struct cellCensus{
//    std::vector<float> xCoord, yCoord, zCoord;
//    std::vector<int> frames;
//    std::vector<int> labelInMap; //read id in label_maps
//    std::vector<std::vector<size_t>> voxIdx;
//    std::vector<std::vector<int>> vox_x, vox_y, vox_z;
//    std::vector<std::vector<size_t>> tracks; //a set of node_ids
//    std::vector<size_t> nodeID2trackID; //particle2track in matlab
//    std::vector<size_t> nodeLocInTrack; //particle2track in matlab
//    std::vector<std::vector<size_t>> parents, kids; // neighboring relationship, at most two kids or parents
//};
void cellTrackingMain:: cellInfoAccumuate(cellSegmentMain &cellSegment){
    extractVoxIdxList(cellSegment.cell_label_maps, movieInfo.voxIdx, cellSegment.number_cells);
    FOREACH_i(movieInfo.voxIdx){
        vec_ind2sub(movieInfo.voxIdx[i], movieInfo.vox_y[i],
                    movieInfo.vox_x[i], movieInfo.vox_z[i], cellSegment.cell_label_maps[0].size);

        movieInfo.xCoord[i] = vec_mean(movieInfo.vox_x[i]);
        movieInfo.yCoord[i] = vec_mean(movieInfo.vox_y[i]);
        movieInfo.zCoord[i] = vec_mean(movieInfo.vox_z[i]);
    }
    movieInfo.nodes.reserve(movieInfo.xCoord.size());
    long increment = 0;
    FOREACH_i(cellSegment.number_cells){
        increment += i>0? cellSegment.number_cells[i]:0;
        for(int j = 0; j < cellSegment.number_cells[i]; j++){
            movieInfo.frames[j + increment] = i;
            movieInfo.labelInMap[j + increment] = j;
            movieInfo.nodes[j + increment].node_id = j + increment;
        }
    }

    cumulative_cell_nums.resize(cellSegment.number_cells.size());
    cumulative_cell_nums[0] = cellSegment.number_cells[0];
    for (int i = 1; i< cumulative_cell_nums.size(); i++){
        cumulative_cell_nums[i] = (cumulative_cell_nums[i-1] + cellSegment.number_cells[i]);
    }
}
/**
 * @brief initTransitionCost: initialize transition cost
 */
void cellTrackingMain::initTransitionCost(cellSegmentMain &cellSegment){
    FOREACH_i(movieInfo.nodes){
        vector<size_t> nei_ids;
        extractNeighborIds(cellSegment.cell_label_maps, i, nei_ids);
        FOREACH_j(nei_ids){
            // cal the distance between two nodes
        }
    }

}
/**
 * @brief extractNeighborIds: for a given node idx, get its possible neighbors
 * @param cell_label_maps
 * @param node_idx
 * @param nei_idxs
 */
void cellTrackingMain::extractNeighborIds(vector<Mat> &cell_label_maps, size_t cell_idx, vector<size_t> &nei_idxs){

    for (size_t i = movieInfo.frames[cell_idx] + 1; i < cell_label_maps.size(); i++){
        unordered_set<int> nei_labels;
        FOREACH_j(movieInfo.voxIdx[cell_idx]){
            if (cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j])>0){
                nei_labels.insert(cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j]));
            }
        }
        for(auto it : nei_labels){
            nei_idxs.push_back(it + cumulative_cell_nums[i - 1]);
        }
    }
}
/**
 * @brief cellInfo2graph: using the movieInfo to build the graph for CINDA
 */
void cellTrackingMain::cellInfo2graph(){

}
