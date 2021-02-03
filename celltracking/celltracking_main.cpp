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
    movieInfo.xCoord.resize(movieInfo.voxIdx.size());
    movieInfo.yCoord.resize(movieInfo.voxIdx.size());
    movieInfo.zCoord.resize(movieInfo.voxIdx.size());
    movieInfo.range_xyz.resize(movieInfo.voxIdx.size());
    movieInfo.start_coord_xyz.resize(movieInfo.voxIdx.size());
    FOREACH_i(movieInfo.voxIdx){
        vec_ind2sub(movieInfo.voxIdx[i], movieInfo.vox_y[i],
                    movieInfo.vox_x[i], movieInfo.vox_z[i], cellSegment.cell_label_maps[0].size);
        movieInfo.xCoord[i] = vec_mean(movieInfo.vox_x[i]);
        movieInfo.yCoord[i] = vec_mean(movieInfo.vox_y[i]);
        movieInfo.zCoord[i] = vec_mean(movieInfo.vox_z[i]);

        movieInfo.start_coord_xyz[i].resize(3);
        size_t dummy;
        movieInfo.start_coord_xyz[i][0] = vec_min(movieInfo.vox_x[i], dummy);
        movieInfo.start_coord_xyz[i][1] = vec_min(movieInfo.vox_y[i], dummy);
        movieInfo.start_coord_xyz[i][2] = vec_min(movieInfo.vox_z[i], dummy);
        movieInfo.range_xyz[i].resize(3);
        movieInfo.range_xyz[i][0] = vec_max(movieInfo.vox_x[i], dummy) - movieInfo.range_xyz[i][0] + 1;
        movieInfo.range_xyz[i][1] = vec_max(movieInfo.vox_y[i], dummy) - movieInfo.range_xyz[i][1] + 1;
        movieInfo.range_xyz[i][2] = vec_max(movieInfo.vox_z[i], dummy) - movieInfo.range_xyz[i][2] + 1;
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

    movieInfo.frame_shift.resize(cellSegment.number_cells.size());
    FOREACH_i(movieInfo.frame_shift){
        movieInfo.frame_shift[i] = {0, 0, 0};
    }
    cumulative_cell_nums.resize(cellSegment.number_cells.size());
    cumulative_cell_nums[0] = cellSegment.number_cells[0];
    for (int i = 1; i< cumulative_cell_nums.size(); i++){
        cumulative_cell_nums[i] = (cumulative_cell_nums[i-1] + cellSegment.number_cells[i]);
    }

}
/**
 * @brief voxelwise_avg_distance: use matlab's matrix indexing method
 * @param cell_curr
 * @param cell_nei
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(size_t cell_curr, size_t cell_nei, float &c2n, float &n2c){
    bool *ref_cell = new bool[movieInfo.voxIdx[cell_curr].size()];
    memset(ref_cell, false, movieInfo.voxIdx[cell_curr].size());
    size_t idx;
    size_t frame_sz = movieInfo.range_xyz[cell_curr][0] * movieInfo.range_xyz[cell_curr][1];
    for(size_t i = 0; i < movieInfo.voxIdx[cell_curr].size(); i++){
        idx = (movieInfo.vox_z[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][2]) * frame_sz +
               (movieInfo.vox_x[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][0]) * movieInfo.range_xyz[cell_curr][1] +
                movieInfo.vox_y[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][1];
        ref_cell[idx] = true;
    }
    bool *mov_cell = new bool[movieInfo.voxIdx[cell_nei].size()];
    memset(mov_cell, false, movieInfo.voxIdx[cell_nei].size());
    frame_sz = movieInfo.range_xyz[cell_nei][0] * movieInfo.range_xyz[cell_nei][1];
    for(size_t i = 0; i < movieInfo.voxIdx[cell_nei].size(); i++){
        idx = (movieInfo.vox_z[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][2]) * frame_sz +
               (movieInfo.vox_x[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][0]) * movieInfo.range_xyz[cell_nei][1] +
                movieInfo.vox_y[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][1];
        mov_cell[idx] = true;
    }
    vector<float> dist(2);
    vector<float> shift_xyz(3);
    shift_xyz = vec_Minus(movieInfo.frame_shift[movieInfo.frames[cell_curr]], movieInfo.frame_shift[movieInfo.frames[cell_nei]]);
    distanceTransRegion2Region(ref_cell, movieInfo.range_xyz[cell_curr],
                               mov_cell, movieInfo.range_xyz[cell_nei],
                               shift_xyz, dist);
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
