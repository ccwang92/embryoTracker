#include "celltracking_main.h"
#include "CINDA/src_c/cinda_funcs.c"

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
    /////////////////////////////////////////////////
    //   step 3. main loop for cell tracking       //
    /////////////////////////////////////////////////
    int loop_cnt = 1;
    while (loop_cnt <= p4tracking.maxIter){
        // MODULE ONE: split/merge test from current tracking results
        for(int dummy = 0; dummy < 3; dummy ++){
            split_merge_module(cellSegment);
        }
        // MODULE TWO: retrieve missing cells
        missing_cell_module(cellSegment);
        loop_cnt ++;
    }
}

void cellTrackingMain::cellInfoAccumuate(cellSegmentMain &cellSegment){
    vector<vector<size_t>> curr_cell_idx;
    extractVoxIdxList(cellSegment.cell_label_maps, curr_cell_idx, cellSegment.number_cells);
    /** KEY: we generate an overly large voxIdx in movieInfo */
    movieInfo.voxIdx.resize(2*curr_cell_idx.size()); // we assume each cell at most can be split once
    long long increment = 0;
    for(int i : cellSegment.number_cells){
        for(size_t j=increment; j<(increment+i); j++){
            movieInfo.voxIdx[j] = curr_cell_idx[i];
        }
        for(size_t j=increment+i; j<(increment+2*i); j++){
            movieInfo.voxIdx[j].reserve(0);
        }
        increment += i*2;
    }
    cumulative_cell_nums.resize(cellSegment.number_cells.size());
    cumulative_cell_nums[0] = cellSegment.number_cells[0]*2;
    for (int i = 1; i< cumulative_cell_nums.size(); i++){
        cumulative_cell_nums[i] = (cumulative_cell_nums[i-1] + 2*cellSegment.number_cells[i]);
    }
    movieInfo.xCoord.resize(movieInfo.voxIdx.size());
    movieInfo.yCoord.resize(movieInfo.voxIdx.size());
    movieInfo.zCoord.resize(movieInfo.voxIdx.size());
    movieInfo.range_xyz.resize(movieInfo.voxIdx.size());
    movieInfo.start_coord_xyz.resize(movieInfo.voxIdx.size());
    FOREACH_i(movieInfo.voxIdx){
        if(movieInfo.voxIdx[i].size()>0){
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
        }else{
            movieInfo.vox_y[i].reserve(0);
            movieInfo.vox_x[i].reserve(0);
            movieInfo.vox_z[i].reserve(0);
            movieInfo.xCoord[i] = -INFINITY;
            movieInfo.yCoord[i] = -INFINITY;
            movieInfo.zCoord[i] = -INFINITY;
            movieInfo.start_coord_xyz[i] = {-1,-1,-1};
            movieInfo.range_xyz[i] = {0,0,0};
        }
    }
    movieInfo.nodes.reserve(movieInfo.xCoord.size());
    increment = 0;
    FOREACH_i(cellSegment.number_cells){
        increment += i>0? 2*cellSegment.number_cells[i]:0;
        for(int j = 0; j < 2*cellSegment.number_cells[i]; j++){
            movieInfo.frames[j + increment] = i;
            movieInfo.labelInMap[j + increment] = 0;
            movieInfo.nodes[j + increment].node_id = j + increment;
        }
    }

    movieInfo.frame_shift_xyz.resize(cellSegment.number_cells.size());
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i] = {0, 0, 0};
    }

    movieInfo.tracks.resize(0);

    validGapMaps.resize(cellSegment.number_cells.size());
    for(int i = 0; i < cellSegment.number_cells.size(); i++){
        validGapMaps = Mat(cellSegment.cell_label_maps[i].dims, cellSegment.cell_label_maps[i].size,
                           CV_8U, Scalar(255));
    }
}
/**
 * @brief voxelwise_avg_distance: use matlab's matrix indexing method
 * @param cell_curr: id of current cell
 * @param cell_nei: id of neighbor cell
 * @param c2n: distance from cell_curr to cell_nei
 * @param n2c: distance from cell_nei to cell_curr
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(size_t cell_curr, size_t cell_nei, float &c2n, float &n2c){
//    bool *ref_cell = new bool[movieInfo.voxIdx[cell_curr].size()];
//    memset(ref_cell, false, movieInfo.voxIdx[cell_curr].size());
//    size_t idx;
//    size_t frame_sz = movieInfo.range_xyz[cell_curr][0] * movieInfo.range_xyz[cell_curr][1];
//    for(size_t i = 0; i < movieInfo.voxIdx[cell_curr].size(); i++){
//        idx = (movieInfo.vox_z[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][2]) * frame_sz +
//               (movieInfo.vox_x[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][0]) * movieInfo.range_xyz[cell_curr][1] +
//                movieInfo.vox_y[cell_curr][i] - movieInfo.start_coord_xyz[cell_curr][1];
//        ref_cell[idx] = true;
//    }
//    bool *mov_cell = new bool[movieInfo.voxIdx[cell_nei].size()];
//    memset(mov_cell, false, movieInfo.voxIdx[cell_nei].size());
//    frame_sz = movieInfo.range_xyz[cell_nei][0] * movieInfo.range_xyz[cell_nei][1];
//    for(size_t i = 0; i < movieInfo.voxIdx[cell_nei].size(); i++){
//        idx = (movieInfo.vox_z[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][2]) * frame_sz +
//               (movieInfo.vox_x[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][0]) * movieInfo.range_xyz[cell_nei][1] +
//                movieInfo.vox_y[cell_nei][i] - movieInfo.start_coord_xyz[cell_nei][1];
//        mov_cell[idx] = true;
//    }
//    vector<float> dist(2);
//    vector<double> shift_xyz(3); // the true shift from ref cell to mov cell (:nei_start-ref_start-drift)
//    shift_xyz = vec_Minus(movieInfo.frame_shift_xyz[movieInfo.frames[cell_curr]],
//            movieInfo.frame_shift_xyz[movieInfo.frames[cell_nei]]); //this is -drift
//    shift_xyz[0] += vec_min(movieInfo.vox_x[cell_nei]) - vec_min(movieInfo.vox_x[cell_curr]);
//    shift_xyz[1] += vec_min(movieInfo.vox_y[cell_nei]) - vec_min(movieInfo.vox_y[cell_curr]);
//    shift_xyz[2] += vec_min(movieInfo.vox_z[cell_nei]) - vec_min(movieInfo.vox_z[cell_curr]);
//    float max_dist = distanceTransRegion2Region(ref_cell, movieInfo.range_xyz[cell_curr],
//                               mov_cell, movieInfo.range_xyz[cell_nei],
//                               shift_xyz, dist);
//    c2n = dist[0];
//    n2c = dist[1];
//    return max_dist;

    return voxelwise_avg_distance(movieInfo.voxIdx[cell_curr], movieInfo.vox_x[cell_curr], movieInfo.vox_y[cell_curr],
                                  movieInfo.vox_z[cell_curr], movieInfo.range_xyz[cell_curr], movieInfo.start_coord_xyz[cell_curr],
                                  movieInfo.frames[cell_curr], movieInfo.voxIdx[cell_nei], movieInfo.vox_x[cell_nei], movieInfo.vox_y[cell_nei],
                                  movieInfo.vox_z[cell_nei], movieInfo.range_xyz[cell_nei], movieInfo.start_coord_xyz[cell_nei],
                                  movieInfo.frames[cell_nei], c2n, n2c);
}
float cellTrackingMain::voxelwise_avg_distance(size_t joint_cells_curr[], size_t joint_cells_nei[], float &c2n, float &n2c){
    vector<size_t> curr(joint_cells_curr, joint_cells_curr + sizeof(joint_cells_curr)/sizeof(joint_cells_curr[0]));
    vector<size_t> nei(joint_cells_nei, joint_cells_nei + sizeof(joint_cells_nei)/sizeof(joint_cells_nei[0]));

    return voxelwise_avg_distance(curr, nei, c2n, n2c);
}
/**
 * @brief voxelwise_avg_distance: the distance between two regions, each region contains >=1 cells
 * @param joint_cells_curr
 * @param joint_cells_nei
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &joint_cells_curr, vector<size_t> &joint_cells_nei, float &c2n, float &n2c){
    combinedCellsCensus curr_joint_reg, nei_joint_reg;
    combineCellsIntoOneRegion(joint_cells_curr, curr_joint_reg);
    combineCellsIntoOneRegion(joint_cells_nei, nei_joint_reg);
    return voxelwise_avg_distance(curr_joint_reg, nei_joint_reg, c2n, n2c);
}

float cellTrackingMain::voxelwise_avg_distance(size_t cell_curr, vector<size_t> &joint_cells_nei, float &c2n, float &n2c){
    combinedCellsCensus nei_joint_reg;
    combineCellsIntoOneRegion(joint_cells_nei, nei_joint_reg);
    return voxelwise_avg_distance(cell_curr, nei_joint_reg, c2n, n2c);
}
float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &joint_cells_curr, size_t cell_nei, float &c2n, float &n2c){
    return voxelwise_avg_distance(cell_nei, joint_cells_curr, n2c, c2n);
}
/**
 * @brief voxelwise_avg_distance
 * @param curr
 * @param nei
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(combinedCellsCensus &curr, combinedCellsCensus &nei, float &c2n, float &n2c){
//    bool *ref_cell = new bool[curr.voxIdx.size()];
//    memset(ref_cell, false, curr.voxIdx.size());
//    size_t idx;
//    size_t frame_sz = curr.range_xyz[0] * curr.range_xyz[1];
//    for(size_t i = 0; i < curr.voxIdx.size(); i++){
//        idx = (curr.vox_z[i] - curr.start_coord_xyz[2]) * frame_sz +
//               (curr.vox_x[i] - curr.start_coord_xyz[0]) * curr.range_xyz[1] +
//                curr.vox_y[i] - curr.start_coord_xyz[1];
//        ref_cell[idx] = true;
//    }
//    bool *mov_cell = new bool[nei.voxIdx.size()];
//    memset(mov_cell, false, nei.voxIdx.size());
//    frame_sz = nei.range_xyz[0] * nei.range_xyz[1];
//    for(size_t i = 0; i < nei.voxIdx.size(); i++){
//        idx = (nei.vox_z[i] - nei.start_coord_xyz[2]) * frame_sz +
//               (nei.vox_x[i] - nei.start_coord_xyz[0]) * nei.range_xyz[1] +
//                nei.vox_y[i] - nei.start_coord_xyz[1];
//        mov_cell[idx] = true;
//    }
//    vector<float> dist(2);
//    vector<double> shift_xyz(3); // the true shift from ref cell to mov cell (:nei_start-ref_start-drift)
//    shift_xyz = vec_Minus(movieInfo.frame_shift_xyz[curr.frame],
//            movieInfo.frame_shift_xyz[nei.frame]); //this is -drift
//    shift_xyz[0] += vec_min(nei.vox_x) - vec_min(curr.vox_x);
//    shift_xyz[1] += vec_min(nei.vox_y) - vec_min(curr.vox_y);
//    shift_xyz[2] += vec_min(nei.vox_z) - vec_min(curr.vox_z);
//    float max_dist = distanceTransRegion2Region(ref_cell, curr.range_xyz,
//                               mov_cell, nei.range_xyz,
//                               shift_xyz, dist);
//    c2n = dist[0];
//    n2c = dist[1];
//    return max_dist;

    return voxelwise_avg_distance(curr.voxIdx, curr.vox_x, curr.vox_y, curr.vox_z, curr.range_xyz,
                                  curr.start_coord_xyz, curr.frame,
                                  nei.voxIdx, nei.vox_x, nei.vox_y, nei.vox_z, nei.range_xyz,
                                  nei.start_coord_xyz, nei.frame,
                                  c2n, n2c);
}
/**
 * @brief voxelwise_avg_distance
 * @param curr
 * @param nei
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(size_t curr, combinedCellsCensus &nei, float &c2n, float &n2c){
//    voxelwise_avg_distance(vector<size_t> &curr_voxIdx, vector<int> &curr_vox_x, vector<int> &curr_vox_y,
//                                 vector<int> &curr_vox_z, vector<int> &curr_range_xyz, vector<int> &curr_start_coord_xyz, int curr_frame,
//                                 vector<size_t> &nei_voxIdx, vector<int> &nei_vox_x, vector<int> &nei_vox_y,
//                                 vector<int> &nei_vox_z, vector<int> &nei_range_xyz, vector<int> &nei_start_coord_xyz, int nei_frame,
//                                 float &c2n, float &n2c)
    return voxelwise_avg_distance(movieInfo.voxIdx[curr], movieInfo.vox_x[curr], movieInfo.vox_y[curr],
                                  movieInfo.vox_z[curr], movieInfo.range_xyz[curr], movieInfo.start_coord_xyz[curr],
                                  movieInfo.frames[curr], nei.voxIdx, nei.vox_x, nei.vox_y, nei.vox_z, nei.range_xyz,
                                  nei.start_coord_xyz, nei.frame, c2n, n2c);
}
/**
 * @brief voxelwise_avg_distance
 * @param curr
 * @param nei
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(combinedCellsCensus &curr, size_t nei, float &c2n, float &n2c){
    return voxelwise_avg_distance(nei, curr, n2c, c2n);
}
/**
 * @brief voxelwise_avg_distance: the final function to calculate the vox-wise distance between two regions
 * @param curr_voxIdx
 * @param curr_vox_x
 * @param curr_vox_y
 * @param curr_vox_z
 * @param curr_range_xyz
 * @param curr_start_coord_xyz
 * @param curr_frame
 * @param nei_voxIdx
 * @param nei_vox_x
 * @param nei_vox_y
 * @param nei_vox_z
 * @param nei_range_xyz
 * @param nei_start_coord_xyz
 * @param nei_frame
 * @param c2n
 * @param n2c
 * @return
 */
float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &curr_voxIdx, vector<int> &curr_vox_x, vector<int> &curr_vox_y,
                             vector<int> &curr_vox_z, vector<int> &curr_range_xyz, vector<int> &curr_start_coord_xyz, int curr_frame,
                             vector<size_t> &nei_voxIdx, vector<int> &nei_vox_x, vector<int> &nei_vox_y,
                             vector<int> &nei_vox_z, vector<int> &nei_range_xyz, vector<int> &nei_start_coord_xyz, int nei_frame,
                             float &c2n, float &n2c){
    if(!(curr_voxIdx.size() > 0 && nei_voxIdx.size() > 0)){
        c2n = INFINITY;
        n2c = c2n;
        return INFINITY;
    }
    bool *ref_cell = new bool[curr_voxIdx.size()];
    memset(ref_cell, false, curr_voxIdx.size());
    size_t idx;
    size_t frame_sz = curr_range_xyz[0] * curr_range_xyz[1];
    for(size_t i = 0; i < curr_voxIdx.size(); i++){
        idx = (curr_vox_z[i] - curr_start_coord_xyz[2]) * frame_sz +
               (curr_vox_x[i] - curr_start_coord_xyz[0]) * curr_range_xyz[1] +
                curr_vox_y[i] - curr_start_coord_xyz[1];
        ref_cell[idx] = true;
    }
    bool *mov_cell = new bool[nei_voxIdx.size()];
    memset(mov_cell, false, nei_voxIdx.size());
    frame_sz = nei_range_xyz[0] * nei_range_xyz[1];
    for(size_t i = 0; i < nei_voxIdx.size(); i++){
        idx = (nei_vox_z[i] - nei_start_coord_xyz[2]) * frame_sz +
               (nei_vox_x[i] - nei_start_coord_xyz[0]) * nei_range_xyz[1] +
                nei_vox_y[i] - nei_start_coord_xyz[1];
        mov_cell[idx] = true;
    }
    vector<float> dist(2);
    vector<double> shift_xyz(3); // the true shift from ref cell to mov cell (:nei_start-ref_start-drift)
    shift_xyz = vec_Minus(movieInfo.frame_shift_xyz[curr_frame],
            movieInfo.frame_shift_xyz[nei_frame]); //this is -drift
    shift_xyz[0] += vec_min(nei_vox_x) - vec_min(curr_vox_x);
    shift_xyz[1] += vec_min(nei_vox_y) - vec_min(curr_vox_y);
    shift_xyz[2] += vec_min(nei_vox_z) - vec_min(curr_vox_z);
    float max_dist = distanceTransRegion2Region(ref_cell, curr_range_xyz,
                               mov_cell, nei_range_xyz,
                               shift_xyz, dist);
    c2n = dist[0];
    n2c = dist[1];
    return max_dist;
}

float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
                             vector<size_t> &nei_voxIdx, int nei_frame,
                             int data3d_sz[3], float &c2n, float &n2c){
    if(!(curr_voxIdx.size() > 0 && nei_voxIdx.size() > 0)){
        c2n = INFINITY;
        n2c = c2n;
        return INFINITY;
    }
    vector<int> curr_vox_x, curr_vox_y, curr_vox_z;
    vector<int> curr_range_xyz(3), curr_start_coord_xyz(3);
    vector<int> nei_vox_x, nei_vox_y, nei_vox_z;
    vector<int> nei_range_xyz(3), nei_start_coord_xyz(3);

    vec_ind2sub(curr_voxIdx, curr_vox_y, curr_vox_x, curr_vox_z, data3d_sz);
    curr_start_coord_xyz[0] = vec_min(curr_vox_x);
    curr_start_coord_xyz[1] = vec_min(curr_vox_y);
    curr_start_coord_xyz[2] = vec_min(curr_vox_z);
    curr_range_xyz[0] = vec_max(curr_vox_x) - curr_start_coord_xyz[0] + 1;
    curr_range_xyz[1] = vec_max(curr_vox_y) - curr_start_coord_xyz[1] + 1;
    curr_range_xyz[2] = vec_max(curr_vox_z) - curr_start_coord_xyz[2] + 1;

    vec_ind2sub(nei_voxIdx, nei_vox_y, nei_vox_x, nei_vox_z, data3d_sz);
    nei_start_coord_xyz[0] = vec_min(nei_vox_x);
    nei_start_coord_xyz[1] = vec_min(nei_vox_y);
    nei_start_coord_xyz[2] = vec_min(nei_vox_z);
    nei_range_xyz[0] = vec_max(nei_vox_x) - nei_start_coord_xyz[0] + 1;
    nei_range_xyz[1] = vec_max(nei_vox_y) - nei_start_coord_xyz[1] + 1;
    nei_range_xyz[2] = vec_max(nei_vox_z) - nei_start_coord_xyz[2] + 1;


    bool *ref_cell = new bool[curr_voxIdx.size()];
    memset(ref_cell, false, curr_voxIdx.size());
    size_t idx;
    size_t frame_sz = curr_range_xyz[0] * curr_range_xyz[1];
    for(size_t i = 0; i < curr_voxIdx.size(); i++){
        idx = (curr_vox_z[i] - curr_start_coord_xyz[2]) * frame_sz +
               (curr_vox_x[i] - curr_start_coord_xyz[0]) * curr_range_xyz[1] +
                curr_vox_y[i] - curr_start_coord_xyz[1];
        ref_cell[idx] = true;
    }
    bool *mov_cell = new bool[nei_voxIdx.size()];
    memset(mov_cell, false, nei_voxIdx.size());
    frame_sz = nei_range_xyz[0] * nei_range_xyz[1];
    for(size_t i = 0; i < nei_voxIdx.size(); i++){
        idx = (nei_vox_z[i] - nei_start_coord_xyz[2]) * frame_sz +
               (nei_vox_x[i] - nei_start_coord_xyz[0]) * nei_range_xyz[1] +
                nei_vox_y[i] - nei_start_coord_xyz[1];
        mov_cell[idx] = true;
    }
    vector<float> dist(2);
    vector<double> shift_xyz(3); // the true shift from ref cell to mov cell (:nei_start-ref_start-drift)
    shift_xyz = vec_Minus(movieInfo.frame_shift_xyz[curr_frame],
            movieInfo.frame_shift_xyz[nei_frame]); //this is -drift
    shift_xyz[0] += vec_min(nei_vox_x) - vec_min(curr_vox_x);
    shift_xyz[1] += vec_min(nei_vox_y) - vec_min(curr_vox_y);
    shift_xyz[2] += vec_min(nei_vox_z) - vec_min(curr_vox_z);
    float max_dist = distanceTransRegion2Region(ref_cell, curr_range_xyz,
                               mov_cell, nei_range_xyz,
                               shift_xyz, dist);
    c2n = dist[0];
    n2c = dist[1];
    return max_dist;
}

float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
                             vector<size_t> &nei_voxIdx, int nei_frame,
                             MatSize data3d_sz, float &c2n, float &n2c){
    int data_sz[3] = {data3d_sz[0], data3d_sz[1], data3d_sz[2]};
    return voxelwise_avg_distance(curr_voxIdx, curr_frame, nei_voxIdx, nei_frame, data_sz, c2n, n2c);
}
/**
 * @brief combineCellsIntoOneRegion: make a region by joining several cells
 * @param cell_idxes
 * @param out_region_info
 */
void cellTrackingMain::combineCellsIntoOneRegion(vector<size_t> &cell_idxes, combinedCellsCensus &out_region_info){
    out_region_info.frame = movieInfo.frames[cell_idxes[0]];
    out_region_info.start_coord_xyz.resize(3);
    fill(out_region_info.range_xyz.begin(), out_region_info.range_xyz.end(), -1);

    for(size_t n : cell_idxes){
        if(out_region_info.voxIdx.size() == 0) continue;
        out_region_info.voxIdx.insert(out_region_info.voxIdx.end(), movieInfo.voxIdx[n].begin(), movieInfo.voxIdx[n].end());
        out_region_info.vox_x.insert(out_region_info.vox_x.end(), movieInfo.vox_x[n].begin(), movieInfo.vox_x[n].end());
        out_region_info.vox_y.insert(out_region_info.vox_y.end(), movieInfo.vox_y[n].begin(), movieInfo.vox_y[n].end());
        out_region_info.vox_z.insert(out_region_info.vox_z.end(), movieInfo.vox_z[n].begin(), movieInfo.vox_z[n].end());

        if(out_region_info.range_xyz[0] == -1){
            out_region_info.range_xyz = movieInfo.range_xyz[n];
            out_region_info.start_coord_xyz =  movieInfo.start_coord_xyz[n];
        }else{
            if(movieInfo.start_coord_xyz[n][0] < out_region_info.start_coord_xyz[0]){
                out_region_info.range_xyz[0] += (out_region_info.start_coord_xyz[0] - movieInfo.start_coord_xyz[n][0]);
                out_region_info.start_coord_xyz[0] = movieInfo.start_coord_xyz[n][0];
            }
            if(movieInfo.start_coord_xyz[n][1] < out_region_info.start_coord_xyz[1]){
                out_region_info.range_xyz[1] += (out_region_info.start_coord_xyz[1] - movieInfo.start_coord_xyz[n][1]);
                out_region_info.start_coord_xyz[1] = movieInfo.start_coord_xyz[n][1];
            }
            if(movieInfo.start_coord_xyz[n][2] < out_region_info.start_coord_xyz[2]){
                out_region_info.range_xyz[2] += (out_region_info.start_coord_xyz[2] - movieInfo.start_coord_xyz[n][2]);
                out_region_info.start_coord_xyz[2] = movieInfo.start_coord_xyz[n][2];
            }
        }
    }
}

/**
 * @brief updatePreNeighborInfo: update the information of candidate parents of a node
 */
void cellTrackingMain::updatePreNeighborInfo(){
    //// update the preNeighbors information
    for (nodeInfo n : movieInfo.nodes){
        for (nodeRelation neighbor : n.neighbors){
            nodeRelation tmp;
            tmp.node_id = n.node_id;
            tmp.link_cost = neighbor.link_cost;
            tmp.overlap_size = neighbor.overlap_size;
            movieInfo.nodes[neighbor.node_id].preNeighbors.push_back(tmp); //pre-nei has no dist_c2n or dist_n2c
        }
    }
}
size_t cellTrackingMain::cellOverlapSize(size_t c0, size_t c1, cellSegmentMain &cellSegment){
    size_t overlap_sz = 0;
    int curr_frame = movieInfo.frames[c0];
    int curr_cell_label = movieInfo.labelInMap[c0];
    for(size_t k : movieInfo.voxIdx[c1]){
        if(cellSegment.cell_label_maps[curr_frame].at<int>(k) == curr_cell_label){
            overlap_sz ++;
        }
    }
    return overlap_sz;
}
/**
 * @brief distance2cost
 * @param distance
 * @param alpha
 * @param beta
 * @param punish
 * @return
 */
float cellTrackingMain::distance2cost(float distance, float alpha, float beta, float punish = 1){
    float p = gammacdf(distance, alpha, beta);
    float cost = normInv(p * punish / 2);
    return cost*cost;
}
/**
 * @brief calCell2neighborDistance: get the cell-cell voxelwise distance
 * @return
 */
void cellTrackingMain::calCell2neighborDistance(vector<float> &nn_dist){
    float max_dist;
    nn_dist.resize(movieInfo.nodes.size());
    size_t nn_dist_cnt = 0;
    if (movieInfo.tracks.empty()){ // if there is no track info, use all node with nearest neighbor
        FOREACH_i(movieInfo.nodes){
            nn_dist[nn_dist_cnt] = INFINITY;
            FOREACH_j(movieInfo.nodes[i].neighbors){
                // cal the distance between two nodes
                max_dist = voxelwise_avg_distance(i, movieInfo.nodes[i].neighbors[j].node_id,
                                                  movieInfo.nodes[i].neighbors[j].dist_c2n,
                                                movieInfo.nodes[i].neighbors[j].dist_n2c);
                if (max_dist < nn_dist[nn_dist_cnt]){
                    nn_dist[nn_dist_cnt] = max_dist;
                }
            }
            nn_dist_cnt ++;
        }
    }else{ // if there are tracks, use long tracks
        FOREACH_i(movieInfo.tracks){
            if(movieInfo.tracks[i].size() < p4tracking.validtrackLength4var){
                continue;
            }
            for(size_t j=0; j<movieInfo.tracks[i].size()-1; j++){
                nn_dist[nn_dist_cnt] = movieInfo.nodes[movieInfo.tracks[i][j]].kid_cost[0];
                nn_dist_cnt ++;
            }
        }
    }
    nn_dist.resize(nn_dist_cnt);
}
/**
 * @brief initTransitionCost: initialize transition cost
 */
void cellTrackingMain::initTransitionCost(cellSegmentMain &cellSegment){
    // discover the neighbors for each cell
    movieInfo.overall_neighbor_num = 0;
    FOREACH_i(movieInfo.nodes){
        movieInfo.nodes[i].detect_confidence = p4tracking.observationCost;
        movieInfo.nodes[i].in_cost = p4tracking.c_en;
        movieInfo.nodes[i].out_cost = p4tracking.c_ex;
        vector<size_t> nei_ids(0);
        extractNeighborIds(cellSegment.cell_label_maps, i, nei_ids);
        movieInfo.nodes[i].neighbors.reserve(nei_ids.size());
        movieInfo.nodes[i].preNeighbors.reserve(0); //initialized preNeighbor
        if (nei_ids.size() < 1){
            continue;
        }
        movieInfo.overall_neighbor_num += nei_ids.size();
        FOREACH_j(nei_ids){
            movieInfo.nodes[i].neighbors[j].overlap_size = cellOverlapSize(i, nei_ids[j], cellSegment);
        }
    }
    // shift initialize as all 0s
    movieInfo.frame_shift_xyz.resize(cellSegment.number_cells.size());
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i] = {0, 0, 0};
    }
    // calculate the distances between each cell and its neighbors
    updateGammaParam();
    updateArcCost(); // calculate the link cost from given new gamma parameter
    // initialize all node as not stable
    FOREACH_i(movieInfo.nodes){
        movieInfo.nodes[i].stable_status = NOT_STABLE;
    }
    // start a naive linking
    mccTracker_one2one(); // with no split/merge vectors
    updateGammaParam();
    updateArcCost();
    stableSegmentFixed();
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
 * @brief mccTracker_one2one: build the graph without split and merge settings, but allows jump
 */
void cellTrackingMain::mccTracker_one2one(){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      We directly use the interface I wrote for python with the following function name               //
    // long long* pyCS2(long *msz, double *mtail, double *mhead, double *mlow, double *macap, double *mcost)//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long n = 2*movieInfo.nodes.size() + 1; // # vertices in graph
    long long m = movieInfo.overall_neighbor_num + 3*movieInfo.nodes.size(); // # arcs in graph (upper bound)
    double *mtail = new double[m];
    double *mhead = new double[m];
    double *mlow = new double[m]; memset(mlow, 0, m);
    double *macap = new double[m]; memset(macap, 1, m);
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 0;
    int rand_max = 100;
    float linkCostUpBound = INFINITY / 1e7;
    for(nodeInfo node : movieInfo.nodes){
        // in arc
        mtail[arc_cnt] = src_id;
        mhead[arc_cnt] =  2 * node.node_id+1;
        mcost[arc_cnt] = round(node.in_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // observation arc
        mtail[arc_cnt] = 2 * node.node_id+1;
        mhead[arc_cnt] =  2 * node.node_id + 2;
        mcost[arc_cnt] = round(p4tracking.observationCost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // out arc
        mhead[arc_cnt] = src_id;
        mtail[arc_cnt] =  2 * node.node_id + 2;
        mcost[arc_cnt] = round(node.out_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        linkCostUpBound = node.out_cost;
        // link with neighbors
        for(nodeRelation neighbor : node.neighbors){
            if(neighbor.link_cost < linkCostUpBound){
                mtail[arc_cnt] = 2 * node.node_id + 2;
                mhead[arc_cnt] =  2 * neighbor.node_id + 1;
                mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;
            }

        }
    }
    //assert(arc_cnt == m);
    long *msz = new long[3]; // should be long long, but for our data, long is enough
    msz[0] = 12;
    msz[1] = n;
    msz[2] = arc_cnt;

    // call pyCS2: note, everything needs to be very precise
    long long *track_vec = pyCS2(msz, mtail, mhead, mlow, macap, mcost);

    // update movieInfo.tracks
    vector<long long> cost;
    vector<size_t> curr_track;
    for(long long i = 1; i < track_vec[0] + 1; i++){
        if(track_vec[i] > 0){
            curr_track.push_back(track_vec[i]);
        }
        else{
            cost.push_back(track_vec[i]);
            vector<size_t> new_track;
            for(size_t j = 0; j < curr_track.size(); j+=2){
                new_track.push_back(size_t(curr_track[j] / 2));
            }
            if (new_track.size() > 1){
                movieInfo.tracks.push_back(new_track);
            }
        }
    }
    // update parent and kids given tracks
    track2parentKid();

    // update the jumpCost in p4tracking if no split/merge allowed
    updateJumpCost();
    // by the way, also update stable segmentations
    stableSegmentFixed();
    if (false){
        //TODO: re-link jump-over cells, this will only be called in the last iteration

    }

    delete[] mtail;
    delete[] mhead;
    delete[] mlow;
    delete[] macap;
    delete[] mcost;
    free(track_vec);
}
/**
 * @brief mccTracker_splitMerge: tracking using CINDA allows for split/merge, but no jump
 * @param split_merge_node_info
 * @param node_id2split_merge_node_id
 */
void cellTrackingMain::mccTracker_splitMerge(vector<splitMergeNodeInfo> &split_merge_node_info){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      We directly use the interface I wrote for python with the following function name               //
    // long long* pyCS2(long *msz, double *mtail, double *mhead, double *mlow, double *macap, double *mcost)//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long m_append = split_merge_node_info.size();
    // build the graph without split and merge settings
    long long n = 2*movieInfo.nodes.size() + 1; // # vertices in graph
    long long m = movieInfo.overall_neighbor_num + 3*movieInfo.nodes.size() + m_append; // # arcs in graph (upper bound)
    double *mtail = new double[m];
    double *mhead = new double[m];
    double *mlow = new double[m]; memset(mlow, 0, m);
    double *macap = new double[m]; memset(macap, 1, m);
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 0;
    int rand_max = 100;
    bool link_adj_frs = false;
    if(p4tracking.splitMergeHandle == NOJUMPALL){
        link_adj_frs = true;
    }
    bool *visited = new bool[movieInfo.nodes.size()];
    memset(visited, false, movieInfo.nodes.size());

    float linkCostUpBound = INFINITY / 1e7;
    // add the arcs of splitting and merging
    for(splitMergeNodeInfo sp_mg_info : split_merge_node_info){
        if(sp_mg_info.node_id >= 0){
            if(sp_mg_info.parent_flag){
                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.family_nodes[0] + 1;
                mcost[arc_cnt] = round(sp_mg_info.link_costs[0] * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.family_nodes[1] + 1;
                mcost[arc_cnt] = round(sp_mg_info.link_costs[1] * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = src_id;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 2;
                mcost[arc_cnt] = round(sp_mg_info.src_link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                visited[sp_mg_info.node_id] = true;
            }else{
                mtail[arc_cnt] = 2 * sp_mg_info.family_nodes[0] + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 1;
                mcost[arc_cnt] = round(sp_mg_info.link_costs[0] * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.family_nodes[1] + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 1;
                mcost[arc_cnt] = round(sp_mg_info.link_costs[1] * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 1;
                mhead[arc_cnt] =  src_id;
                mcost[arc_cnt] = round(sp_mg_info.src_link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                visited[sp_mg_info.family_nodes[0]] = true;
                visited[sp_mg_info.family_nodes[1]] = true;
            }
        }
    }
    for(nodeInfo node : movieInfo.nodes){
        if(visited[node.node_id]){
            continue;
        }
        // in arc
        mtail[arc_cnt] = src_id;
        mhead[arc_cnt] =  2 * node.node_id + 1;
        mcost[arc_cnt] = round(node.in_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // observation arc
        mtail[arc_cnt] = 2 * node.node_id+1;
        mhead[arc_cnt] =  2 * node.node_id + 2;
        mcost[arc_cnt] = round(p4tracking.observationCost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // out arc
        mhead[arc_cnt] = src_id;
        mtail[arc_cnt] =  2 * node.node_id+2;
        mcost[arc_cnt] = round(node.out_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // link with neighbors
        if(link_adj_frs){
            int curr_frame = movieInfo.frames[node.node_id];
            for(nodeRelation neighbor : node.neighbors){
                if(movieInfo.frames[neighbor.node_id]-1 == curr_frame){
                    if(neighbor.link_cost < linkCostUpBound){
                        mtail[arc_cnt] = 2 * node.node_id+2;
                        mhead[arc_cnt] =  2 * neighbor.node_id+1;
                        mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                        arc_cnt ++;
                    }
                }
            }
        }else{
            for(nodeRelation neighbor : node.neighbors){
                if(neighbor.link_cost < linkCostUpBound){
                    mtail[arc_cnt] = 2 * node.node_id+1;
                    mhead[arc_cnt] =  2 * neighbor.node_id+2;
                    mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                    arc_cnt ++;
                }
            }
        }
    }
    //assert(arc_cnt == m);
    long *msz = new long[3]; // should be long long, but for our data, long is enough
    msz[0] = 12;
    msz[1] = n;
    //msz[2] = m;
    msz[2] = arc_cnt; // m is the upper bound

    // call pyCS2: note, everything needs to be very precise
    long long *track_vec = pyCS2(msz, mtail, mhead, mlow, macap, mcost);

    // update movieInfo.tracks
    vector<long long> cost;
    vector<size_t> curr_track;
    for(long long i = 1; i < track_vec[0] + 1; i++){
        if(track_vec[i] > 0){
            curr_track.push_back(track_vec[i]);
        }
        else{
            cost.push_back(track_vec[i]);
            vector<size_t> new_track;
            for(size_t j = 0; j < curr_track.size(); j+=2){
                new_track.push_back(size_t(curr_track[j] / 2));
            }
            if (new_track.size() > 1){
                movieInfo.tracks.push_back(new_track);
            }
        }
    }
    // update parent and kids given tracks
    track2parentKid();
    if(m_append > 0){// merge tracks with overlapped nodes
        mergeOvTracks();
    }

    delete[] mtail;
    delete[] mhead;
    delete[] mlow;
    delete[] macap;
    delete[] mcost;
    delete[] visited;
    free(track_vec);
}
void cellTrackingMain::refreshTracks(){
    vector<std::vector<size_t>> new_tracks(movieInfo.tracks.size());
    long valid_num = 0;
    for(vector<size_t> track : movieInfo.tracks){
        if(track[0] >= 0){
            new_tracks[valid_num] = track;
        }
    }
    new_tracks.resize(valid_num);
    movieInfo.tracks.clear();
    movieInfo.tracks = new_tracks;
}
/**
 * @brief mergeOvTrack_node2trackUpt: merge overlapped tracks and update node id to track information
 */
void cellTrackingMain::mergeOvTracks(){
    // first reset all node to track ids as -1
    vector<long> tmp_nodeId2trackId(movieInfo.xCoord.size());
    fill(tmp_nodeId2trackId.begin(), tmp_nodeId2trackId.end(), -1);
    size_t curr_id;
    int overlapped_track_id;
    for(size_t t = 0; t < movieInfo.tracks.size(); t++){
        overlapped_track_id = -1;
        for (size_t i = 1; i < movieInfo.tracks[t].size(); i++){
            curr_id = movieInfo.tracks[t][i];
            if(tmp_nodeId2trackId[curr_id] > 0){
                overlapped_track_id = tmp_nodeId2trackId[curr_id];
                break;
            }else{
                tmp_nodeId2trackId[curr_id] = t;
            }
        }
        if(overlapped_track_id > 0){
            for (size_t i = 1; i < movieInfo.tracks[t].size(); i++){
                tmp_nodeId2trackId[movieInfo.tracks[t][i]] = overlapped_track_id;
            }
            vector<size_t> tmp_track = movieInfo.tracks[overlapped_track_id];
            movieInfo.tracks[overlapped_track_id].resize(tmp_track.size() + movieInfo.tracks[t].size());
            merge(tmp_track.begin(), tmp_track.end(), movieInfo.tracks[t].begin(), movieInfo.tracks[t].end(),
                  movieInfo.tracks[overlapped_track_id].begin());
            movieInfo.tracks[t][0] = -1;
        }
    }
    refreshTracks();
}
void cellTrackingMain::node2trackUpt(){
    FOREACH_i(movieInfo.tracks){
        FOREACH_j(movieInfo.tracks[i]){
            movieInfo.nodes[movieInfo.tracks[i][j]].nodeId2trackId = i;
            movieInfo.nodes[movieInfo.tracks[i][j]].nodeLocInTrack = j;
        }
    }
}
void cellTrackingMain::updateJumpCost(){
    size_t overall_arcs = 0;
    vector<size_t> jumpArcs(p4tracking.k);
    fill(jumpArcs.begin(), jumpArcs.end(), 0);
    for(vector<size_t> track : movieInfo.tracks){
        overall_arcs += track.size() - 1;
        for (size_t i = 1; i < track.size(); i++){
            jumpArcs[movieInfo.frames[track[i]] - movieInfo.frames[track[i-1]] - 1] ++;
        }
    }
    FOREACH_i(p4tracking.jumpCost){
        p4tracking.jumpCost[i] = (float)((double)jumpArcs[i] / overall_arcs);
    }
}
/**
 * @brief track2parentKid: !! this function should be called before we do track merge
 */
void cellTrackingMain::track2parentKid(){
    for(nodeInfo n : movieInfo.nodes){
        n.parent_num = 0;
        n.kid_num = 0;
    }
    for(vector<size_t> track : movieInfo.tracks){
        for (size_t i = 1; i < track.size(); i++){
            movieInfo.nodes[track[i-1]].kids[movieInfo.nodes[track[i-1]].kid_num] = track[i];
            for (int k = 0; k<movieInfo.nodes[track[i-1]].neighbors.size(); k++){
                if (track[i] == movieInfo.nodes[track[i-1]].neighbors[k].node_id){
                    movieInfo.nodes[track[i-1]].kid_cost[movieInfo.nodes[track[i-1]].kid_num] =
                            movieInfo.nodes[track[i-1]].neighbors[k].link_cost;
                    break;
                }
            }
            movieInfo.nodes[track[i-1]].kid_num ++;

            movieInfo.nodes[track[i]].parents[movieInfo.nodes[track[i]].parent_num] = track[i-1];
            for (int k = 0; k<movieInfo.nodes[track[i]].preNeighbors.size(); k++){
                if (track[i] == movieInfo.nodes[track[i]].preNeighbors[k].node_id){
                    movieInfo.nodes[track[i]].parent_cost[movieInfo.nodes[track[i]].parent_num] =
                            movieInfo.nodes[track[i]].preNeighbors[k].link_cost;
                    break;
                }
            }
            movieInfo.nodes[track[i]].parent_num ++;
        }
    }
}
// functions to update cost given new gamma fitting results
void cellTrackingMain::driftCorrection(){
    long num_frames = movieInfo.frame_shift_xyz.size();
    vector<long> samples4driftCumulate(num_frames);
    size_t pre_node_id, curr_node_id;
    int curr_frame;
    for(vector<size_t> track : movieInfo.tracks){
        if(track.size() >= p4tracking.validtrackLength4var){
            for (size_t i = 1; i < track.size(); i++){
                pre_node_id = track[i-1];
                curr_node_id = track[i];
                curr_frame = movieInfo.frames[curr_node_id];
                if(curr_frame - movieInfo.frames[pre_node_id] == 1){ // only consider the adjacent moving
                    movieInfo.frame_shift_xyz[curr_frame][0] += movieInfo.xCoord[curr_node_id] - movieInfo.xCoord[pre_node_id];
                    movieInfo.frame_shift_xyz[curr_frame][1] += movieInfo.yCoord[curr_node_id] - movieInfo.yCoord[pre_node_id];
                    movieInfo.frame_shift_xyz[curr_frame][2] += movieInfo.zCoord[curr_node_id] - movieInfo.zCoord[pre_node_id];
                    samples4driftCumulate[curr_frame] ++;
                }
            }
        }
    }
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i][0] /= samples4driftCumulate[curr_frame];
        movieInfo.frame_shift_xyz[i][1] /= samples4driftCumulate[curr_frame];
        movieInfo.frame_shift_xyz[i][2] /= samples4driftCumulate[curr_frame];
    }
}
void cellTrackingMain::updateGammaParam(){
    vector<float> nn_dist;
    calCell2neighborDistance(nn_dist);
    truncatedGammafit(nn_dist, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
}
void cellTrackingMain::updateArcCost(){
    for(nodeInfo n : movieInfo.nodes){
        for(nodeRelation &neighbor : n.neighbors){
            neighbor.link_cost = distance2cost(MAX(neighbor.dist_c2n, neighbor.dist_n2c),
                                               movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
        }
    }
}

void cellTrackingMain::getArcCostOne2OneTrack(size_t track_id, vector<float> &arc_costs){
    arc_costs.resize(movieInfo.tracks[track_id].size()-1);
    for(size_t i=0; i<movieInfo.tracks[track_id].size()-1; i++){
        arc_costs[i] = movieInfo.nodes[movieInfo.tracks[track_id][i]].kid_cost[0];
    }
}
/**
 * @brief stableArcFixed: some segmentation are stable, should not be split/merge any more
 * This function will only be called when we has one to one linking.
 */
void cellTrackingMain::stableSegmentFixed(){
    movieInfo.track_arcs_avg_mid_std.resize(movieInfo.tracks.size());
    float max_cost = pow(normInv(0.5*p4tracking.jumpCost[0]/2), 2);
    int min_valid_node_cluster_sz = 5;
    FOREACH_i(movieInfo.tracks){
        vector<float> arc_costs;
        getArcCostOne2OneTrack(i, arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][0] = vec_mean(arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][1] = vec_median(arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][2] = vec_stddev(arc_costs);

        size_t start_id = 0, end_id = 0;
        for(size_t j = 0; j<arc_costs.size(); j++){
            if (arc_costs[j] <= max_cost){
                end_id = j + 1;
            }else{
                if (end_id - start_id + 1 >= min_valid_node_cluster_sz){
                    for(size_t k = start_id + 1; k <= end_id - 1; k++){
                        movieInfo.nodes[movieInfo.tracks[i][k]].stable_status = STABLE_TRACK_MID;
                    }
                    if(movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status == NOT_STABLE ||
                            movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status == STABLE_TRACK_HEAD){
                        movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status = STABLE_TRACK_HEAD;
                    }else{
                        movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status = STABLE_TRACK_MID;
                    }
                    if(movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status == NOT_STABLE ||
                            movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status == STABLE_TRACK_END){
                        movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status = STABLE_TRACK_END;
                    }else{
                        movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status = STABLE_TRACK_MID;
                    }
                }
                start_id = j+1;
                end_id = j+1;
            }
        }

        if (end_id == arc_costs.size() && (end_id - start_id + 1) >= min_valid_node_cluster_sz){
            for(size_t k = start_id + 1; k <= end_id - 1; k++){
                movieInfo.nodes[movieInfo.tracks[i][k]].stable_status = STABLE_TRACK_MID;
            }
            if(movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status == NOT_STABLE ||
                    movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status == STABLE_TRACK_HEAD){
                movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status = STABLE_TRACK_HEAD;
            }else{
                movieInfo.nodes[movieInfo.tracks[i][start_id]].stable_status = STABLE_TRACK_MID;
            }
            if(movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status == NOT_STABLE ||
                    movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status == STABLE_TRACK_END){
                movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status = STABLE_TRACK_END;
            }else{
                movieInfo.nodes[movieInfo.tracks[i][end_id]].stable_status = STABLE_TRACK_MID;
            }
        }
    }
}
/**
 * @brief bestPeerCandidate: given a region id, find if there are two regions that can together fit into its territory
 * @param node_id
 * @param bestPeer
 * @param parent_flag
 */
float cellTrackingMain::bestPeerCandidate(size_t node_id, vector<size_t> &out_bestPeer, bool parent_flag){
    vector<size_t> peerCandidates(0);
    out_bestPeer.resize(0);
    if(parent_flag && movieInfo.nodes[node_id].neighbors.size() >= 2){
        float curr_cost;
        bool currBest;
        for(nodeRelation neighbor : movieInfo.nodes[node_id].neighbors){
            curr_cost = neighbor.dist_c2n;
            currBest = true;
            for(nodeRelation preNeighbor : movieInfo.nodes[neighbor.node_id].preNeighbors){
                if(curr_cost > preNeighbor.dist_c2n && preNeighbor.node_id != node_id){
                    currBest = false;
                    break;
                }
            }
            if(currBest) peerCandidates.push_back(neighbor.node_id);
        }
    }else if(movieInfo.nodes[node_id].preNeighbors.size() >= 2){
        float curr_cost;
        bool currBest;
        for(nodeRelation preNeighbor : movieInfo.nodes[node_id].preNeighbors){
            curr_cost = preNeighbor.dist_c2n;
            currBest = true;
            for(nodeRelation neighbor : movieInfo.nodes[preNeighbor.node_id].neighbors){
                if(curr_cost > neighbor.dist_c2n && neighbor.node_id != node_id){
                    currBest = false;
                    break;
                }
            }
            if(currBest) peerCandidates.push_back(preNeighbor.node_id);
        }
    }
    if(peerCandidates.size() < 2){
        return INFINITY;
    }else if(peerCandidates.size() == 2){
        out_bestPeer = peerCandidates;
        float dummy_c2n, dummy_n2c;
        return voxelwise_avg_distance(node_id, out_bestPeer, dummy_c2n, dummy_n2c);
    }else{ // more than two regions available
        float min_distance = INFINITY, curr_distance;
        float dummy_c2n, dummy_n2c;
        vector<size_t> curr_peer(2);
        FOREACH_i(peerCandidates){
            curr_peer[0] = peerCandidates[i];
            for(size_t j = i + 1; j < peerCandidates.size(); j ++){
                curr_peer[0] = peerCandidates[j];
                curr_distance = voxelwise_avg_distance(node_id, curr_peer, dummy_c2n, dummy_n2c);
                if(curr_distance < min_distance){
                    out_bestPeer = curr_peer;
                    min_distance = curr_distance;
                }
            }
        }
        return min_distance;
    }
}

/**
 * @brief peerRegionVerify: to verify if there are two regions in adjacent frame that jointly has better linking cost
 * @param node_id
 * @param cost_good2go
 * @param parents_test
 * @param split_merge_node_info
 */
void cellTrackingMain::peerRegionVerify(size_t node_id, float cost_good2go, bool parents_test,
                                        vector<splitMergeNodeInfo> &split_merge_node_info,
                                        unordered_map<long, long> &node_id2split_merge_node_id){
    float merged_cost, min_exist_cost = INFINITY;
    vector<size_t> bestPeers(2);
    merged_cost = bestPeerCandidate(node_id, bestPeers, parents_test);
    if(parents_test){
        for(nodeRelation neighbor : movieInfo.nodes[node_id].neighbors){
            if(neighbor.dist_c2n < min_exist_cost){
                min_exist_cost = neighbor.dist_c2n;
            }
        }
    }else{
        for(nodeRelation preNeighbor : movieInfo.nodes[node_id].preNeighbors){
            if(preNeighbor.dist_c2n < min_exist_cost){
                min_exist_cost = preNeighbor.dist_c2n;
            }
        }
    }
    splitMergeNodeInfo split_merge_peer;
    if(merged_cost < min_exist_cost && min_exist_cost > cost_good2go){
        if (p4tracking.splitMergeCost){
            merged_cost /= 2;
        }
        if(merged_cost < MIN(p4tracking.c_en, p4tracking.c_ex)){
            // add to the list
            split_merge_peer.family_nodes[0] = bestPeers[0];
            split_merge_peer.family_nodes[1] = bestPeers[1];
            split_merge_peer.link_costs[0] = merged_cost;
            split_merge_peer.link_costs[1] = merged_cost;
            split_merge_peer.node_id = node_id;
            split_merge_peer.parent_flag = parents_test;
            split_merge_peer.src_link_cost = 0.0;
            split_merge_peer.invalid = false;
            split_merge_node_info.push_back(split_merge_peer);
            // 'long' is enough for our data analysis, indeed all the size_t are redundantly large
            if(parents_test){
                node_id2split_merge_node_id.insert(make_pair<long, long>(node_id, split_merge_node_info.size()));
            }else{
                node_id2split_merge_node_id.insert(make_pair<long, long>(-node_id, split_merge_node_info.size()));
            }
        }
    }
}
/**
 * @brief detectPeerRegions: test if two regions should be merged based on their common kid or parent region
 * @param merge_node_idx
 * @param split_node_idx
 */
void cellTrackingMain::detectPeerRegions(vector<splitMergeNodeInfo> &split_merge_node_info,
                                         unordered_map<long, long> &node_id2split_merge_node_id){
    float cost_good2go; // the average cost of arcs in a track: only arc with larger cost will be test
    for(nodeInfo node : movieInfo.nodes){
        cost_good2go = movieInfo.track_arcs_avg_mid_std[node.nodeId2trackId][0];
        if (cost_good2go > MIN(p4tracking.c_en, p4tracking.c_ex)){
            cost_good2go = 0;
        }
        bool parents_test = true;
        bool kids_test = true;
        if(p4tracking.stableNodeTest){ // for stable node, no need to split/merge
            switch (node.stable_status){
                case STABLE_TRACK_HEAD:
                    kids_test = false;
                    break;
                case STABLE_TRACK_MID:
                    parents_test = false;
                    kids_test = false;
                    break;
                case STABLE_TRACK_END:
                    parents_test = false;
                    break;
                default:
                    break;
            }
        }

        if(parents_test){
            peerRegionVerify(node.node_id, cost_good2go, true, split_merge_node_info, node_id2split_merge_node_id);
        }
        if(kids_test){
            peerRegionVerify(node.node_id, cost_good2go, false, split_merge_node_info, node_id2split_merge_node_id);
        }
    }
}
/**
 * @brief sizeCumulate: calculate the overall size of three nodes
 * @param curr_cell
 * @param familiy_members
 * @return
 */
size_t cellTrackingMain::sizeCumulate(size_t curr_cell, size_t familiy_members[2]){
    return movieInfo.voxIdx[curr_cell].size() + movieInfo.voxIdx[familiy_members[0]].size() +
            movieInfo.voxIdx[familiy_members[1]].size();
}
/**
 * @brief handleMergeSplitRegions: build the network allowing split/merge
 */
void cellTrackingMain::handleMergeSplitRegions(){
    vector<splitMergeNodeInfo> split_merge_node_info;
    unordered_map<long, long> node_id2split_merge_node_id;
    detectPeerRegions(split_merge_node_info, node_id2split_merge_node_id);
    // 1. check these contradict conditions: e.g. a->b+c and a+d->b
    for(size_t ss = 0; ss < split_merge_node_info.size(); ss++){
        auto sp_mgInfo = split_merge_node_info[ss];
        if(sp_mgInfo.invalid) continue;
        vector<pair<size_t,size_t>> contradict_groups;
        if(sp_mgInfo.parent_flag){ // current is a->b+c, we are finding if there is d+e->b or f+g->c
            for(int i = 0; i < 2; i++){
                auto it = node_id2split_merge_node_id.find(-sp_mgInfo.family_nodes[i]);
                if (it != node_id2split_merge_node_id.end() &&
                        !split_merge_node_info[it->second].invalid){
                    contradict_groups.push_back({it->second,
                                                 sizeCumulate(sp_mgInfo.family_nodes[i],
                                                 split_merge_node_info[it->second].family_nodes)});
                }
            }
        }else{
            for(int i = 0; i < 2; i++){
                auto it = node_id2split_merge_node_id.find(sp_mgInfo.family_nodes[i]);
                if (it != node_id2split_merge_node_id.end() &&
                        !split_merge_node_info[it->second].invalid){
                    contradict_groups.push_back({it->second,
                                                 sizeCumulate(sp_mgInfo.family_nodes[i],
                                                 split_merge_node_info[it->second].family_nodes)});
                }
            }
        }
        if (contradict_groups.size() > 0){
            contradict_groups.push_back({ss,
                                            sizeCumulate(sp_mgInfo.node_id, sp_mgInfo.family_nodes)});
            int best_idx = 0;
            size_t max_group_sz = 0;
            for(int i=0; i<contradict_groups.size(); i++){
                auto it = contradict_groups[i];
                if (it.second > max_group_sz){
                    max_group_sz = it.second;
                    best_idx = i;
                }
            }
            for(int i=0; i<contradict_groups.size(); i++){
                if(i != best_idx){ // remove contradict ones
                    split_merge_node_info[ss].invalid = true; //
                    //node_id2split_merge_node_id.erease(-sp_mgInfo.family_nodes[i]);
                }
            }
        }
    }
    mccTracker_splitMerge(split_merge_node_info);

}
/**
 * @brief handleInconsistentParentKid: we assume the inconsistent parents/kids are not all wrong, either
 * parents or kids are correct. We split the current cell region to see which one is correct.
 *
 * @param node_id
 * @return
 */
int cellTrackingMain::handleInconsistentParentKid(cellSegmentMain &cellSegment, size_t node_id){
    // test if the two parents can be used to split the cell region reprented by node_id
    vector<vector<size_t>> reg4seeds(2);
    reg4seeds[0] = movieInfo.voxIdx[movieInfo.nodes[node_id].parents[0]];
    reg4seeds[1] = movieInfo.voxIdx[movieInfo.nodes[node_id].parents[1]];
    int reg4seeds_frame = movieInfo.frames[movieInfo.nodes[node_id].parents[0]];
    bool gapBasedSplit = true;
    vector<vector<size_t>> splitRegs(2);
    float reg4seeds2splitRes_costs[2];
//    bisectValidTest(cellSegmentMain &cellSegment, size_t reg2split_idx, vector<size_t> reg2split,
//                                           int reg2split_frame, vector<vector<size_t>> reg4seeds, int reg4seeds_frame,
//                                           bool gapBasedSplit, vector<vector<size_t>> &splitRegs,
//                                           float *reg4seeds2splitRes_costs)
    bool valid_p = bisectValidTest(cellSegment, node_id, movieInfo.voxIdx[node_id], movieInfo.frames[node_id],
                        reg4seeds, reg4seeds_frame, gapBasedSplit, true,
                         splitRegs, reg4seeds2splitRes_costs);
    if(!valid_p){
        gapBasedSplit = false;
        valid_p = bisectValidTest(cellSegment, node_id, movieInfo.voxIdx[node_id], movieInfo.frames[node_id],
                        reg4seeds, reg4seeds_frame, gapBasedSplit, true,
                         splitRegs, reg4seeds2splitRes_costs);
    }

    // test if the two kids can be used to split the cell region reprented by node_id
//    vector<vector<size_t>> reg4seeds(2);
    reg4seeds[0] = movieInfo.voxIdx[movieInfo.nodes[node_id].kids[0]];
    reg4seeds[1] = movieInfo.voxIdx[movieInfo.nodes[node_id].kids[1]];
    reg4seeds_frame = movieInfo.frames[movieInfo.nodes[node_id].kids[0]];
    gapBasedSplit = true;
    bool valid_k = bisectValidTest(cellSegment, node_id, movieInfo.voxIdx[node_id], movieInfo.frames[node_id],
                        reg4seeds, reg4seeds_frame, gapBasedSplit, true,
                         splitRegs, reg4seeds2splitRes_costs);
    if(!valid_k){
        gapBasedSplit = false;
        valid_k = bisectValidTest(cellSegment, node_id, movieInfo.voxIdx[node_id], movieInfo.frames[node_id],
                        reg4seeds, reg4seeds_frame, gapBasedSplit, true,
                         splitRegs, reg4seeds2splitRes_costs);
    }
    if(valid_p && valid_k){
        return MERGE_BOTH_PARENTS_KIDS;
    }else if(valid_p){
        return SPLIT_BY_PARENTS;
    }else if(valid_k){
        return SPLIT_BY_KIDS;
    }else{
        return MERGE_SPLIT_NOT_SURE;
    }
    
}
int cellTrackingMain::parentsKidsConsistency(size_t node_id){
    if(movieInfo.nodes[node_id].parent_num != 2 || movieInfo.nodes[node_id].kid_num != 2 ||
            !p4tracking.par_kid_consistency_check){
        return CONSISTENCY_NOT_SURE;
    }
    float dummy_c2n, dummy_n2c;
    float merged_distance = voxelwise_avg_distance(movieInfo.nodes[node_id].parents, movieInfo.nodes[node_id].kids,
                                                   dummy_c2n, dummy_n2c);
    float merged_cost = distance2cost(merged_distance, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);

    if (merged_cost >= abs(p4tracking.observationCost)){
        return CONSISTENCY_NOT_SURE;
    }
    float map_2x2[4];
    memset(map_2x2, INFINITY, 4);
    int cnt = 0;
    for(size_t p : movieInfo.nodes[node_id].parents){
        for(size_t k : movieInfo.nodes[node_id].kids){
            for(nodeRelation nr : movieInfo.nodes[p].neighbors){
                if (nr.node_id == k){
                    map_2x2[cnt] = distance2cost(nr.dist_c2n, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
                }
                cnt++;
            }
        }
    }
    if (MAX(map_2x2[0], map_2x2[3]) < abs(p4tracking.observationCost) ||
            MAX(map_2x2[1], map_2x2[2]) < abs(p4tracking.observationCost) ){
        return PARENTS_KIDS_CONSISTENT;
    }else{
        return PARENTS_KIDS_NOT_CONSISTENT;
    }
}
// bool cellTrackingMain::exist_in_pairs(vector<size_t> &pairs, size_t id){
//     if (pairs.size() == 0) return false;

//     for(size_t p : pairs){
//         if (p == id){
//             return true;
//         }
//     }
//     return false;
// }
/**
 * @brief seed_map_gen: a map refine function since 3d map may be too liberal that removes some valid seeds;
 * We will check if there is a seed totally masked by 3d gap map from 3d
 * principal curvature.
 * @param fgMap
 * @param possibleGap3d
 * @param possibleGap2d
 * @param minSz
 * @param newFgMap: output
 */
void seed_map_gen(Mat1b &fgMap, Mat1b &possibleGap3d, Mat1b &possibleGap2d, Mat1b &newFgMap, size_t minSz = 0){
    // first test 3d pv
    Mat1b tmp0;
    bitwise_not(possibleGap3d, tmp0);
    bitwise_and(fgMap, tmp0, newFgMap);
    if(minSz > 0){
//        Mat1i tmp_label;
//        int numCC = connectedComponents3d(&newFgMap, tmp_label, 26);
//        removeSmallCC(tmp_label, numCC, minSz, false);
//        newFgMap = tmp_label > 0;
        bwareaopenMat(newFgMap, newFgMap, minSz, 26);
    }
    if(!isempty(newFgMap, CV_8U, 0)){ // use 2d pv to test if seeds are removed by 3d pv
        Mat1b newFgMap_2d;
        bitwise_not(possibleGap2d, tmp0);
        bitwise_and(fgMap, tmp0, newFgMap_2d);

        bitwise_or(newFgMap_2d, newFgMap, tmp0);
        Mat1i tmp_label;
        int numCC = connectedComponents3d(&tmp0, tmp_label, 26);

        vector<bool> seeds_covered_by_3d_map(numCC);
        fill(seeds_covered_by_3d_map.begin(), seeds_covered_by_3d_map.end(), false);
        FOREACH_i_MAT(tmp_label){
            if(tmp_label.at<int>(i) > 0){
                seeds_covered_by_3d_map[tmp_label.at<int>(i) - 1] = true;
            }
        }
        unordered_set<int> new_seeds;
        FOREACH_i(seeds_covered_by_3d_map){
            if(!seeds_covered_by_3d_map[i]){
                new_seeds.insert(i+1);
            }
        }
        if (new_seeds.size() > 0){
            FOREACH_i_MAT(tmp_label){
                if(tmp_label.at<int>(i) > 0){
                    auto it = new_seeds.find(tmp_label.at<int>(i));
                    if(it != new_seeds.end()){
                        newFgMap.at<unsigned char>(i) = 255;
                    }
                }
            }
        }
        if(minSz > 0){
            bwareaopenMat(newFgMap, newFgMap, minSz, 26);
        }
    }

}

/**
 * @brief seedsRefine_intensity:remove voxes in seeds that are too dim or bright
 * @param data3d
 * @param root_idx
 * @param root_frame
 * @param seeds_idx
 * @param seed_frame
 * @param ref_seeds_idx
 * @return
 */
bool seedsRefine_intensity(cellSegmentMain &cellSegment, vector<size_t> &root_idx, int root_frame,
                           vector<size_t> &seeds_idx, int seed_frame,
                           vector<size_t> &ref_seeds_idx){
    float low = 0.1, high = 0.9; // remove 20% voxels
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];

    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*seed_frame; // sub-matrix pointer
    Mat *frame_root = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind);

    unsigned char *ind2 = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*root_frame; // sub-matrix pointer
    Mat *frame_seed = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind2);
    ref_seeds_idx.resize(0);

    vector<float> vals = extractValsGivenIdx(frame_seed, seeds_idx, CV_8U);
    sort(vals.begin(), vals.end()); //add greater<>() as the 3rd paramter will make it descending
    float lb = vals[(int)round(vals.size()*low)];
    float ub = vals[(int)round(vals.size()*high)];
    auto tmp = intersection(seeds_idx, root_idx);
    for(size_t j : tmp){
        unsigned char val = frame_root->at<unsigned char>(j);
        if(val > lb && val <ub){
            ref_seeds_idx.push_back(j);
        }
    }
    if(ref_seeds_idx.size() == 0){
        return false;
    }
    return true;

}
/**
 * @brief seedsRefine_intensity:remove voxes in seeds that are too dim or bright
 * @param data3d
 * @param root_idx
 * @param root_frame
 * @param seeds_idx
 * @param seed_frame
 * @param ref_seeds_idx
 * @return
 */
bool cellTrackingMain::seedsRefine_intensity(cellSegmentMain &cellSegment, vector<size_t> &root_idx, int root_frame,
                           vector<vector<size_t>> &seeds_idx, int seed_frame,
                           vector<vector<size_t>> &ref_seeds_idx){
    float low = 0.1, high = 0.9; // remove 20% voxels
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];

    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*seed_frame; // sub-matrix pointer
    Mat *frame_root = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind);

    unsigned char *ind2 = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*root_frame; // sub-matrix pointer
    Mat *frame_seed = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind2);
    ref_seeds_idx.resize(seeds_idx.size());
    for(int i = 0; i<seeds_idx.size(); i++){
        vector<float> vals = extractValsGivenIdx(frame_seed, seeds_idx[i], CV_8U);
        sort(vals.begin(), vals.end()); //add greater<>() as the 3rd paramter will make it descending
        float lb = vals[(int)round(vals.size()*low)];
        float ub = vals[(int)round(vals.size()*high)];
        auto tmp = intersection(seeds_idx[i], root_idx);
        for(size_t j : tmp){
            unsigned char val = frame_root->at<unsigned char>(j);
            if(val > lb && val <ub){
                ref_seeds_idx[i].push_back(j);
            }
        }
        if(ref_seeds_idx.size() == 0){
            return false;
        }
    }
    return true;
}
/**
 * @brief seedRefine_gap
 * @param cellSegment
 * @param root_idx
 * @param root_frame
 * @param seeds_idx
 * @param seed_frame
 * @param ref_seeds_idx
 * @return
 */
bool cellTrackingMain::seedsRefine_gap(Mat1b &possibleGaps, vector<vector<size_t>> &seeds_idx, Mat1i &outLabelMap){
    outLabelMap = Mat(possibleGaps.dims, possibleGaps.size, CV_32S, Scalar(0));
    Mat1b noGaps;
    bitwise_not(possibleGaps, noGaps);
    for(int i = 0; i< seeds_idx.size(); i++){
        Mat1b tmp = Mat(possibleGaps.dims, possibleGaps.size, CV_8U, Scalar(0));
        setValMat(tmp, CV_8U, seeds_idx[i], 255);
        bitwise_and(noGaps, tmp, tmp);
        if(isempty(tmp, CV_8U)){
            return false;
        }else{
            setValMat(outLabelMap, CV_32S, tmp, i+1);
        }
    }

}

/**
 * @brief binary_seedsMap_create: create the seeds map given fgMap and two seed regions
 * @param fgMap
 * @param possibleGap3d
 * @param possibleGap2d
 * @param seeds_idx
 * @param seeds_map
 * @param minSz
 * @return
 */
bool cellTrackingMain::binary_seedsMap_create(Mat1b &fgMap, Mat1b *possibleGap3d,
                                              Mat1b *possibleGap2d, vector<vector<size_t>> seeds_idx, Mat1i &seeds_map, size_t minSz = 0){
    Mat1i label_map;
    int numCC;
    if(possibleGap3d == nullptr){
        Mat1b tmp0,  tmp1;
        bitwise_not(*possibleGap2d, tmp0);
        bitwise_and(fgMap, tmp0, tmp1);
        numCC = connectedComponents3d(&tmp1, label_map, 6);
    }else if(possibleGap2d == nullptr){
        Mat1b tmp0,  tmp1;
        bitwise_not(*possibleGap3d, tmp0);
        bitwise_and(fgMap, tmp0, tmp1);
        numCC = connectedComponents3d(&tmp1, label_map, 6);
    }else{
        Mat1b tmp0;
        seed_map_gen(fgMap, *possibleGap3d, *possibleGap2d, tmp0);
        numCC = connectedComponents3d(&tmp0, label_map, 6);
    }

    if(numCC < 2){
        return false;
    }
    Mat1i tmp_seedMap = Mat(fgMap.dims, fgMap.size, CV_32S, Scalar(0));
    for(int i = 0; i<seeds_idx.size(); i++){
        for(size_t j:seeds_idx[i]){
            tmp_seedMap.at<int>(j) = i + 1;
        }
    }
    vector<vector<size_t>> voxIdxList(numCC);
    extractVoxIdxList(&label_map, voxIdxList, numCC);
    seeds_map = Mat(fgMap.dims, fgMap.size, CV_32S, Scalar(0));
    vector<bool> seeds_idx_used(seeds_idx.size());
    fill(seeds_idx_used.begin(), seeds_idx_used.end(), false);
    for(int i = 0 ; i < numCC; i++){
        vector<int> vals;
        for(size_t j : voxIdxList[i]){
            if(tmp_seedMap.at<int>(j) > 0){
                vals.push_back(tmp_seedMap.at<int>(j));
            }
        }
        int it = Mode(vals.begin(), vals.end());
        if(it != 0){
            setValMat(seeds_map, CV_32S, voxIdxList[i], (float)it);
            seeds_idx_used[it-1] = true;
        }
    }
    for(bool v : seeds_idx_used){
        if(!v){
            return false;
        }
    }
    return true;
}
bool cellTrackingMain::bisectRegion_gapGuided(cellSegmentMain &cellSegment, size_t reg2split_idx,
                            vector<size_t> &reg2split, int reg2split_frame,
                            vector<vector<size_t>> &reg4seeds, bool usePriorGapMap, vector<vector<size_t>> &splitRegs){
    if(reg4seeds.size() != 2 || reg4seeds[0].size() == 0 || reg4seeds[1].size() == 0){
        return false;
    }
    // build score map
    singleCellSeed seed;
    vector<size_t> overall_idx(reg2split);
    overall_idx.insert(overall_idx.end(), reg4seeds[0].begin(), reg4seeds[0].end());
    overall_idx.insert(overall_idx.end(), reg4seeds[1].begin(), reg4seeds[1].end());
    vec_unique(overall_idx);
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];

    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*reg2split_frame; // sub-matrix pointer
    Mat *single_frame = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind);

    cellSegment.cropSeed(-1, overall_idx, single_frame, nullptr, &cellSegment.cell_label_maps[reg2split_frame],
                         reg2split_frame, seed, cellSegment.p4segVol);



    Mat mask = seed.eigMap2d < 0;
    setValMat(seed.eigMap2d, CV_32F, &mask, 0.0);
    Mat1b possibleGap2d = seed.eigMap2d > 0;
    float max_val = getMaxValMat(seed.eigMap2d, CV_32F, reg2split);
    scale_vol(&seed.eigMap2d, CV_32F, &seed.score2d, 0.001, 1, 0, max_val);

    mask = seed.eigMap3d < 0;
    setValMat(seed.eigMap3d, CV_32F, &mask, 0.0);
    Mat1b possibleGap3d = seed.eigMap3d > 0;
    max_val = getMaxValMat(seed.eigMap3d, CV_32F, reg2split);
    scale_vol(&seed.eigMap3d, CV_32F, &seed.score3d, 0.001, 1, 0, max_val);
    if(usePriorGapMap){
        Mat subGapMap;
        subVolExtract(&validGapMaps[reg2split_frame], CV_8U, subGapMap, seed.crop_range_yxz);// deep copy
        bitwise_and(possibleGap2d, subGapMap, possibleGap2d);
        bitwise_and(possibleGap3d, subGapMap, possibleGap3d);
    }

    // start to build seed map
    Mat1b fgMap = seed.idMap == movieInfo.labelInMap[reg2split_idx];
    vector<vector<size_t>> seeds_idx(reg4seeds.size());
    int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    FOREACH_i(reg4seeds){
        coordinateTransfer(reg4seeds[i], single_frame->size, seeds_idx[i], crop_start, seed.idMap.size);
    }
    Mat1i seedsMap;
    bool success = binary_seedsMap_create(fgMap, &possibleGap3d, nullptr, seeds_idx, seedsMap);
    if(!success){
        success = binary_seedsMap_create(fgMap, nullptr, &possibleGap3d, seeds_idx, seedsMap);
        if(!success){
            return false;
        }
    }
    seed.scoreMap = seed.score2d + seed.score3d;
    Mat1i outLabelMap;
    regionGrow(&seedsMap, 2, outLabelMap, &seed.scoreMap, &fgMap,
               cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, false);

    splitRegs.resize(2);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[0], 1);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[1], 2);
    int crop_start_yxz[3] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start, -seed.crop_range_yxz[2].start};
    coordinateTransfer(splitRegs[0], outLabelMap.size, splitRegs[0], crop_start_yxz,single_frame->size);
    coordinateTransfer(splitRegs[1], outLabelMap.size, splitRegs[1], crop_start_yxz,single_frame->size);
}
bool cellTrackingMain::bisectRegion_bruteforce(cellSegmentMain &cellSegment, size_t reg2split_idx,
                                               vector<size_t> &reg2split, int reg2split_frame,
                                               vector<vector<size_t>> &reg4seeds, int reg4seeds_frame,
                                               bool usePriorGapMap, vector<vector<size_t>> &splitRegs){
    if(reg4seeds.size() != 2 || reg4seeds[0].size() == 0 || reg4seeds[1].size() == 0){
        return false;
    }
    // build score map
    singleCellSeed seed;
    vector<size_t> overall_idx(reg2split);
    overall_idx.insert(overall_idx.end(), reg4seeds[0].begin(), reg4seeds[0].end());
    overall_idx.insert(overall_idx.end(), reg4seeds[1].begin(), reg4seeds[1].end());
    vec_unique(overall_idx);
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];

    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*reg2split_frame; // sub-matrix pointer
    Mat *single_frame = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind);

    cellSegment.cropSeed(-1, overall_idx, single_frame, nullptr, &cellSegment.cell_label_maps[reg2split_frame],
                         reg2split_frame, seed, cellSegment.p4segVol);

    Mat mask = seed.eigMap2d < 0;
    setValMat(seed.eigMap2d, CV_32F, &mask, 0.0);
    Mat1b possibleGap2d = seed.eigMap2d > 0;
    float max_val = getMaxValMat(seed.eigMap2d, CV_32F, reg2split);
    scale_vol(&seed.eigMap2d, CV_32F, &seed.score2d, 0.001, 1, 0, max_val);

    mask = seed.eigMap3d < 0;
    setValMat(seed.eigMap3d, CV_32F, &mask, 0.0);
    Mat1b possibleGap3d = seed.eigMap3d > 0;
    max_val = getMaxValMat(seed.eigMap3d, CV_32F, reg2split);
    scale_vol(&seed.eigMap3d, CV_32F, &seed.score3d, 0.001, 1, 0, max_val);
    if(usePriorGapMap){
        Mat subGapMap;
        subVolExtract(&validGapMaps[reg2split_frame], CV_8U, subGapMap, seed.crop_range_yxz);// deep copy
        bitwise_and(possibleGap2d, subGapMap, possibleGap2d);
        bitwise_and(possibleGap3d, subGapMap, possibleGap3d);
    }

    // start to build seed map
    Mat1b fgMap = seed.idMap == movieInfo.labelInMap[reg2split_idx];

    /** ***** THIS the only difference between brute-force and gap-based method ******** **/
    vector<vector<size_t>> ref_seeds_idx(reg4seeds.size());
    bool valid_seeds_flag = true;
    valid_seeds_flag = seedsRefine_intensity(cellSegment, reg2split, reg2split_frame,
                                             reg4seeds, reg4seeds_frame, ref_seeds_idx);
    if (!valid_seeds_flag) return false;

    vector<vector<size_t>> seeds_idx_in_subVol(ref_seeds_idx.size());
    int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    FOREACH_i(reg4seeds){
        coordinateTransfer(ref_seeds_idx[i], single_frame->size, seeds_idx_in_subVol[i], crop_start, seed.idMap.size);
    }
    Mat1i seedsMap;
    valid_seeds_flag = seedsRefine_gap(possibleGap3d, seeds_idx_in_subVol, seedsMap);
    /** ************************************ Done ************************************* **/

    seed.scoreMap = seed.score2d + seed.score3d;
    Mat1i outLabelMap;
    regionGrow(&seedsMap, 2, outLabelMap, &seed.scoreMap, &fgMap,
               cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, false);

    splitRegs.resize(2);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[0], 1);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[1], 2);
    int crop_start_yxz[3] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start, -seed.crop_range_yxz[2].start};
    coordinateTransfer(splitRegs[0], outLabelMap.size, splitRegs[0], crop_start_yxz,single_frame->size);
    coordinateTransfer(splitRegs[1], outLabelMap.size, splitRegs[1], crop_start_yxz,single_frame->size);
}
bool cellTrackingMain::bisectValidTest(cellSegmentMain &cellSegment, size_t reg2split_idx, vector<size_t> reg2split,
                                       int reg2split_frame, vector<vector<size_t>> reg4seeds, int reg4seeds_frame,
                                       bool gapBasedSplit, bool usePriorGapMap, vector<vector<size_t>> &splitRegs,
                                       float *reg4seeds2splitRes_costs){
    splitRegs.resize(reg4seeds.size());

//    bool bisectRegion_gapGuided(cellSegmentMain &cellSegment, size_t reg2split_idx,
//                                vector<size_t> &reg2split, int reg2split_frame,
//                                vector<vector<size_t>> &reg4seeds, vector<vector<size_t>> &splitRegs);
//    bool bisectRegion_bruteforce(cellSegmentMain &cellSegment, size_t reg2split_idx,
//                                 vector<size_t> &reg2split, int reg2split_frame,
//                                 vector<vector<size_t>> &reg4seeds, int reg4seeds_frame,
//                                 vector<vector<size_t>> &splitRegs);
    bool valid_reg_found;
    if(gapBasedSplit){
        valid_reg_found = bisectRegion_gapGuided(cellSegment, reg2split_idx, reg2split,
                               reg2split_frame, reg4seeds, usePriorGapMap, splitRegs);
    }else{
        valid_reg_found = bisectRegion_bruteforce(cellSegment, reg2split_idx, reg2split,
                                                  reg2split_frame, reg4seeds, reg4seeds_frame, usePriorGapMap, splitRegs);
    }
    if(valid_reg_found){
//        float voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
//                                                                           vector<size_t> &nei_voxIdx, int nei_frame,
//                                                                           MatSize data3d_sz, float &c2n, float &n2c);
        float dummy_c2n, dummy_n2c;
        reg4seeds2splitRes_costs[0] = voxelwise_avg_distance(reg2split, reg2split_frame, splitRegs[0], reg4seeds_frame,
                cellSegment.cell_label_maps[0].size, dummy_c2n, dummy_n2c);
        reg4seeds2splitRes_costs[1] = voxelwise_avg_distance(reg2split, reg2split_frame, splitRegs[1], reg4seeds_frame,
                cellSegment.cell_label_maps[0].size, dummy_c2n, dummy_n2c);
        if(MAX(reg4seeds2splitRes_costs[0], reg4seeds2splitRes_costs[1]) < abs(p4tracking.observationCost)){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}
/**
 * @brief testCellsInOneTrackAdjacentOrNot: there can be several regions in the same frame belonging to the same track.
 * If they are spatially disjoint, we can tell that this tracking contains at least two cells.
 * @param left_or_right_cells
 */
bool cellTrackingMain::testCellsInOneTrackAdjacentOrNot(cellSegmentMain &cellSegment, vector<unordered_set<size_t>> left_or_right_cells){
    for(size_t i=1; i < left_or_right_cells.size(); i++){ // neglect the first, which is the current node for test
        if(left_or_right_cells[i].size() > 1){
            // if there is a stable node, return not adjacent(false)
            if (p4tracking.stableNodeTest){
                for(auto nid : left_or_right_cells[i]){
                    assert(movieInfo.nodes[nid].node_id >= 0);
                    if(movieInfo.nodes[nid].stable_status != NOT_STABLE){
                        return false;
                    }
                }
            }
            unordered_set<int> labelInMaps;
            for(auto it : left_or_right_cells[i]){
                labelInMaps.insert(movieInfo.labelInMap[it]);
            }
            size_t test_cnt = 0;
            for(auto it : left_or_right_cells[i]){
                unordered_set<int> tmp;
                assert(it = movieInfo.nodes[it].node_id);
                int frame = movieInfo.frames[it];
                adjacentRegions(cellSegment.cell_label_maps[frame], movieInfo.voxIdx[it],
                                movieInfo.labelInMap[it], tmp);
                tmp.insert(movieInfo.labelInMap[it]);
                if (!set_exist(labelInMaps, tmp)){
                    return false;
                }
                test_cnt ++;
                if((test_cnt+1)==left_or_right_cells[i].size()){
                    break;
                }
            }
        }else{ // once there is a frame that there is only one region, quit the test
            return true;
        }
    }
    return true;
}
/**
 * @brief mergeValidTest: for the linking reuslts a0->a1+b1->a2+b2->...->end, if such track is short,
 * highly likely we should merge a1 and b1.
 * @param curr_node_id
 * @param seedRegs4split
 * @return
 */
bool cellTrackingMain::mergeValidTest(size_t curr_node_id, size_t seedRegs4split[2]){
    int root_fr = movieInfo.frames[curr_node_id];
    int seed_fr = movieInfo.frames[seedRegs4split[0]];
    bool kid_flag = true;
    if (seed_fr < root_fr) kid_flag = false;

    int loopCnt = 0;
    vector<size_t> bases(2);
    vector<size_t> pre_bases(2);
    bases[0] = seedRegs4split[0];
    bases[1] = seedRegs4split[1];
    while (loopCnt < p4tracking.validtrackLength4var){ // if the length is already enough to be a valid track, break;
        loopCnt ++;
        pre_bases = bases;
        if(kid_flag){
            if(movieInfo.nodes[bases[0]].kid_num == 1 && movieInfo.nodes[bases[1]].kid_num == 1){
                bases[0] = movieInfo.nodes[bases[0]].kids[0];
                bases[1] = movieInfo.nodes[bases[1]].kids[0];
            }else if(movieInfo.nodes[bases[0]].kid_num == 0 && movieInfo.nodes[bases[1]].kid_num == 0){
                return true;
            }else if(movieInfo.nodes[bases[0]].kid_num == 1 && movieInfo.nodes[bases[1]].kid_num == 0){
                bases[0] = movieInfo.nodes[bases[0]].kids[0];
                bases[1] = movieInfo.nodes[bases[0]].kids[0];
            }else if(movieInfo.nodes[bases[0]].kid_num == 0 && movieInfo.nodes[bases[1]].kid_num == 1){
                bases[0] = movieInfo.nodes[bases[1]].kids[0];
                bases[1] = movieInfo.nodes[bases[1]].kids[0];
            }else{
                return false;
            }
        }else{
            if(movieInfo.nodes[bases[0]].parent_num == 1 && movieInfo.nodes[bases[1]].parent_num == 1){
                bases[0] = movieInfo.nodes[bases[0]].parents[0];
                bases[1] = movieInfo.nodes[bases[1]].parents[0];
            }else if(movieInfo.nodes[bases[0]].parent_num == 0 && movieInfo.nodes[bases[1]].parent_num == 0){
                return true;
            }else if(movieInfo.nodes[bases[0]].parent_num == 1 && movieInfo.nodes[bases[1]].parent_num == 0){
                bases[0] = movieInfo.nodes[bases[0]].parents[0];
                bases[1] = movieInfo.nodes[bases[0]].parents[0];
            }else if(movieInfo.nodes[bases[0]].parent_num == 0 && movieInfo.nodes[bases[1]].parent_num == 1){
                bases[0] = movieInfo.nodes[bases[1]].parents[0];
                bases[1] = movieInfo.nodes[bases[1]].parents[0];
            }else {
                return false;
            }
        }
        if(bases[0] == bases[1]){
            size_t root_node2 = bases[0];

            if(kid_flag && movieInfo.nodes[root_node2].parent_num > 1){
                bases[0] = movieInfo.nodes[root_node2].parents[0];
                bases[1] = movieInfo.nodes[root_node2].parents[1];
                if((bases[0] != pre_bases[0] && bases[0] != pre_bases[1]) ||
                        (bases[1] != pre_bases[0] && bases[1] != pre_bases[1])){
                    return false;
                }
            }else if(!kid_flag && movieInfo.nodes[root_node2].kid_num > 1){
                bases[0] = movieInfo.nodes[root_node2].kids[0];
                bases[1] = movieInfo.nodes[root_node2].kids[1];
                if((bases[0] != pre_bases[0] && bases[0] != pre_bases[1]) ||
                        (bases[1] != pre_bases[0] && bases[1] != pre_bases[1])){
                    return false;
                }
            }

            float dummy_c2n, dummy_n2c;
            float maxDistance = voxelwise_avg_distance(root_node2, pre_bases, dummy_c2n, dummy_n2c);
            if(distance2cost(maxDistance, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]) < p4tracking.c_ex){
                return true;
            }else{
                return false;
            }
        }
    }
    return false;

}
/** regionSplitMergeJudge: tell if a region represented by curr_node_id should be split or merge its parents/kids
 * 
**/
int cellTrackingMain::regionSplitMergeJudge(cellSegmentMain &cellSegment, size_t curr_node_id, bool one2multiple_flag, float &pvalue){
    int curr_frame = movieInfo.frames[curr_node_id];
    //build the previous tree
    vector<unordered_set<size_t>> left_cells;
    left_cells.push_back({curr_node_id});
    int frame_parents = curr_frame;
    while(true){
        unordered_set<size_t> tmp_parents;
        frame_parents--;
        for(size_t nid : left_cells[left_cells.size()-1]){
            assert(movieInfo.nodes[nid].node_id == nid);
            for(int i = 0 ; i < movieInfo.nodes[nid].parent_num; i++){
                assert(frame_parents == movieInfo.frames[movieInfo.nodes[nid].parents[i]]);
                tmp_parents.insert(movieInfo.nodes[nid].parents[i]);
            }
        }
        
        unordered_set<size_t> tmp_kids;
        int frame_kids = frame_parents + 1;
        for(size_t nid : tmp_parents){
            assert(movieInfo.nodes[nid].node_id == nid);
            for(int i = 0 ; i < movieInfo.nodes[nid].kid_num; i++){
                assert(frame_kids == movieInfo.frames[movieInfo.nodes[nid].kids[i]]);
                tmp_kids.insert(movieInfo.nodes[nid].kids[i]);
            }
        }
        
        if(set_equal(tmp_kids, left_cells[left_cells.size()-1])){
            left_cells.push_back(tmp_parents);
        }else{
            break;
        }
    }
    //build the following tree
    vector<unordered_set<size_t>> right_cells;
    right_cells.push_back({curr_node_id});
    int frame_kids = curr_frame;
    while(true){
        unordered_set<size_t> tmp_kids;
        frame_kids++;
        for(size_t nid : right_cells[right_cells.size()-1]){
            assert(movieInfo.nodes[nid].node_id == nid);
            for(int i = 0 ; i < movieInfo.nodes[nid].kid_num; i++){
                assert(frame_kids == movieInfo.frames[movieInfo.nodes[nid].kids[i]]);
                tmp_kids.insert(movieInfo.nodes[nid].kids[i]);
            }
        }
        unordered_set<size_t> tmp_parents;
        for(size_t nid : tmp_kids){
            assert(movieInfo.nodes[nid].node_id == nid);
            for(int i = 0 ; i < movieInfo.nodes[nid].parent_num; i++){
                assert((frame_kids - 1) == movieInfo.frames[movieInfo.nodes[nid].parents[i]]);
                tmp_parents.insert(movieInfo.nodes[nid].parents[i]);
            }
        }
        if(set_equal(tmp_parents, right_cells[right_cells.size()-1])){
            right_cells.push_back(tmp_kids);
        }else{
            break;
        }
    }

    // start to judge if the region should be merged or split or leave it to next iteration
    // we totally has 5 ways to determine if a one2multiple or multiple2one link is one region or should be split
    pvalue = 0.5;
    float pVal_thres = 0.05;
    /** way 1: the neighbors of the root node lead to consistent conclusion **/
    size_t mergeVSseg[2] = {1, 0}; // 1-vs-more
    float oneCellpValue = 0.5; // initialize as non-significant
    float multiCellpValue = 0.5; // initialize as non-significant
    for(int i = 0; i < MIN(MIN(left_cells.size()-1, right_cells.size()-1), 10); i++){
        if (i+1 < left_cells.size()){
            if (left_cells[i+1].size()>1){
                mergeVSseg[1]++;
            }
            else{
                mergeVSseg[0]++;
            }
        }
        if (i+1 < right_cells.size()){
            if (right_cells[i+1].size()>1){
                mergeVSseg[1]++;
            }else{
                mergeVSseg[0]++;
            }
        }
        oneCellpValue = binocdf(mergeVSseg[0], mergeVSseg[0] + mergeVSseg[1], 0.5);
        multiCellpValue = binocdf(mergeVSseg[1], mergeVSseg[0] + mergeVSseg[1], 0.5);
        if (MIN(multiCellpValue, oneCellpValue) < pVal_thres){
            break;
        }
    }
    if (oneCellpValue < pVal_thres){ // probability of one cell in this region is low
        pvalue = oneCellpValue;
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }
    if (multiCellpValue < pVal_thres){ // probability of multi-cell in this region is low
        pvalue = multiCellpValue;
        if (one2multiple_flag) return MERGE_KIDS;
        else return MERGE_PARENTS;
    }

    /** way 2: if the regions in the same frame are not spatially adjacent
     * or there is stable nodes, split them **/
    if(!testCellsInOneTrackAdjacentOrNot(cellSegment, left_cells) ||
            !testCellsInOneTrackAdjacentOrNot(cellSegment, right_cells)){
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_KIDS;
    }
    /** way 3: segmentation is not consistent in consecutive two frames **/
    // adjframesConsistency.m: I think this is over-fitting, should not add here
    // Its assumption is that if cost(a+b->c+d) < MIN(cost(a+b->c, a+b->d)), a and b should merge.
    /** way 4: the neighbors of the root with similar size node lead to consistent conclusion **/
    int max_num_frame_include = 10;
    vector<long> sz_diff_adj_frames(left_cells.size() + right_cells.size() - 2);
    long sz_diff_adj_frames_cnt = 0;
    vector<long> sz_diff2test_region_left(left_cells.size());
    for (size_t i=0; i< left_cells.size(); i++){
        sz_diff2test_region_left[i] = 0;
        for(auto it : left_cells[i]){
            sz_diff2test_region_left[i] += movieInfo.voxIdx[it].size();
        }
    }
    for (size_t i=1; i< left_cells.size(); i++){
        sz_diff_adj_frames[sz_diff_adj_frames_cnt] = sz_diff2test_region_left[i] -
                sz_diff2test_region_left[i-1];
        sz_diff_adj_frames_cnt++;
    }
    for (size_t i=1; i< left_cells.size(); i++){
        sz_diff2test_region_left[i] = sz_diff2test_region_left[i] -
                sz_diff2test_region_left[0];
    }
    sz_diff2test_region_left[0] = 0;

    long sz_diff2test_region_right[right_cells.size()];
    for (size_t i=0; i< right_cells.size(); i++){
        sz_diff2test_region_right[i] = 0;
        for(auto it : right_cells[i]){
            sz_diff2test_region_right[i] += movieInfo.voxIdx[it].size();
        }
    }
    for (size_t i=1; i< right_cells.size(); i++){
        sz_diff_adj_frames[sz_diff_adj_frames_cnt] = sz_diff2test_region_right[i] -
                sz_diff2test_region_right[i-1];
        sz_diff_adj_frames_cnt++;
    }
    for (size_t i=1; i< right_cells.size(); i++){
        sz_diff2test_region_right[i] = sz_diff2test_region_right[i] -
                sz_diff2test_region_right[0];
    }
    sz_diff2test_region_right[0] = 0;

    float sz_std = vec_stddev(sz_diff_adj_frames);
    mergeVSseg[0] = 1; // 1-vs-more
    mergeVSseg[1] = 0;
    oneCellpValue = 0.5; // initialize as non-significant
    multiCellpValue = 0.5; // initialize as non-significant
    for (size_t i=0; i<MIN(MAX(left_cells.size()-1, right_cells.size()-1), max_num_frame_include); i++){
        if (i+1 < left_cells.size()){
            if (abs(sz_diff2test_region_left[i+1]) < sz_std){
                if (left_cells[i+1].size() > 1){
                    mergeVSseg[1] ++;
                }else{
                    mergeVSseg[0] ++;
                }
            }
        }
        if (i+1 < right_cells.size()){
            if (abs(sz_diff2test_region_right[i+1]) < sz_std){
                if (right_cells[i+1].size()>1){
                    mergeVSseg[1] ++;
                }else{
                    mergeVSseg[0] ++;
                }
            }
        }
        oneCellpValue = binocdf(mergeVSseg[0], mergeVSseg[0] + mergeVSseg[1], 0.5);
        multiCellpValue = binocdf(mergeVSseg[1], mergeVSseg[0] + mergeVSseg[1], 0.5);
        if (MIN(multiCellpValue, oneCellpValue) < pVal_thres){
            break;
        }
    }
    if (oneCellpValue < pVal_thres){ // probability of one cell in this region is low
        pvalue = oneCellpValue;
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }
    if (multiCellpValue < pVal_thres){ // probability of multi-cell in this region is low
        pvalue = multiCellpValue;
        if (one2multiple_flag) return MERGE_KIDS;
        else return MERGE_PARENTS;
    }

    /** way 5: test if the region is a small piece of shit **/
    // Way 5.1: one side show solid evidence of spliting.
    int split_cnt_left = 0, split_cnt_right = 0;
    for(int i = 0; i<MIN(6, left_cells.size()); i++){
        if(left_cells[i].size() > 1){
            split_cnt_left ++;
        }
    }
    for(int i = 0; i<MIN(6, right_cells.size()); i++){
        if(right_cells[i].size() > 1){
            split_cnt_right ++;
        }
    }
    if(split_cnt_right > 4 || split_cnt_left > 4){
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }
    // Way 5.2: more than two ancestors or kids exist.
    vector<size_t> bases(2);
    bool split_again_flag = false, split_consistent_flag = false;
    if(one2multiple_flag){
        bases[0] = movieInfo.nodes[curr_node_id].kids[0];
        bases[1] = movieInfo.nodes[curr_node_id].kids[1];
        if(movieInfo.nodes[bases[0]].kid_num > 1 || movieInfo.nodes[bases[1]].kid_num > 1){
            split_again_flag = true;
        }else if(movieInfo.nodes[bases[0]].kid_num == movieInfo.nodes[bases[1]].kid_num){
            split_consistent_flag = true;
        }
    }else{
        bases[0] = movieInfo.nodes[curr_node_id].parents[0];
        bases[1] = movieInfo.nodes[curr_node_id].parents[1];
        if(movieInfo.nodes[bases[0]].parent_num > 1 || movieInfo.nodes[bases[1]].parent_num > 1){
            split_again_flag = true;
        }else if(movieInfo.nodes[bases[0]].parent_num == movieInfo.nodes[bases[1]].parent_num){
            split_consistent_flag = true;
        }
    }
    if(split_again_flag){
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }else if(split_consistent_flag){
        if (one2multiple_flag) return MERGE_KIDS;
        else return MERGE_PARENTS;
    }
    // Way 5.3: gap exists for spliting.
    vector<vector<size_t>> reg4seeds(2), dummy_splitRegs;
    float dummy_costs[2];
    reg4seeds[0] = movieInfo.voxIdx[bases[0]];
    reg4seeds[1] = movieInfo.voxIdx[bases[1]];
    bool sectTest = bisectValidTest(cellSegment, movieInfo.voxIdx[curr_node_id], movieInfo.frames[curr_node_id],
                                    reg4seeds, movieInfo.frames[bases[0]], false, dummy_splitRegs, *dummy_costs);
    if(sectTest){
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }

    return MERGE_SPLIT_NOT_SURE;
}
/**
 * @brief regionRefresh: based on the tracking results, check if the region should be split or merge
 * @param cellSegment
 */
void cellTrackingMain::regionRefresh(cellSegmentMain &cellSegment, vector<simpleNodeInfo> &newCells, vector<size_t> &uptCell_idxs){
    //vector<tuple<size_t, int, float>> all_merge_split_peers;
    for(vector<size_t> track : movieInfo.tracks){
        if(track.size() < 2 || track[0] < 0){
            continue;
        }
        vector<tuple<size_t, int, float>> merged_split_peers; //node_id, decision, confidence
        unordered_set<size_t> nodes4merge;
        unordered_map<size_t, int> p_k_InConsistenct;
        for(size_t curr_node_id : track){
            if(movieInfo.nodes[curr_node_id].kid_num > 1 && !set_exist(nodes4merge, curr_node_id)){
                if (parentsKidsConsistency(curr_node_id) == PARENTS_KIDS_NOT_CONSISTENT){
                    int curr_decision = handleInconsistentParentKid(cellSegment, curr_node_id);
                    p_k_InConsistenct.insert({curr_node_id, curr_decision});

                    if (curr_decision != MERGE_SPLIT_NOT_SURE){ // if not sure, leave it to future iterations
                        if (curr_decision == MERGE_BOTH_PARENTS_KIDS){
                            if(adjacentRegions(cellSegment.cell_label_maps[movieInfo.frames[movieInfo.nodes[curr_node_id].kids[0]]],
                            movieInfo.voxIdx[movieInfo.nodes[curr_node_id].kids[0]], movieInfo.nodes[curr_node_id].kids[1])){
                                merged_split_peers.push_back(make_tuple(curr_node_id, MERGE_KIDS, INFINITY));
                                nodes4merge.insert(movieInfo.nodes[curr_node_id].kids[0]);
                                nodes4merge.insert(movieInfo.nodes[curr_node_id].kids[1]);
                            }else{
                                merged_split_peers.push_back(make_tuple(curr_node_id, MERGE_KIDS, 0.0));
                            }
                        }
                        else{
                            merged_split_peers.push_back(make_tuple(curr_node_id, curr_decision, 1.0));
                        }
                    }
                }else{
                    float pvalue;
                    int curr_decision = regionSplitMergeJudge(cellSegment, curr_node_id, true, pvalue);
                    if (curr_decision != MERGE_SPLIT_NOT_SURE){
                        merged_split_peers.push_back(make_tuple(curr_node_id, curr_decision, pvalue));
                        if(curr_decision == MERGE_KIDS){
                            nodes4merge.insert(movieInfo.nodes[curr_node_id].kids[0]);
                            nodes4merge.insert(movieInfo.nodes[curr_node_id].kids[1]);
                        }
                    }
                }
                
            }
        }
        // for multiple to one link
        unordered_set<size_t> nodes4split;
        for(size_t n : nodes4merge){
            nodes4split.insert(n);
            for(int j = 0; j< movieInfo.nodes[n].kid_num; j++){
                nodes4split.insert(movieInfo.nodes[n].kids[j]);
            }
        }
        for(size_t curr_node_id : track){
            if(movieInfo.nodes[curr_node_id].parent_num > 1 && !set_exist(nodes4split, curr_node_id)){
                if (movieInfo.nodes[curr_node_id].kid_num > 1){
                    auto it0 = find_if(merged_split_peers.begin(), merged_split_peers.end(),
                                                 [&curr_node_id](const tuple<size_t, int, float>& e){ return get<0>(e) == curr_node_id;} );

                    auto it1 = find_if(p_k_InConsistenct.begin(), p_k_InConsistenct.end(),
                                                 [&curr_node_id](const pair<size_t, int>& e){ return e.first == curr_node_id;} );
                    if(it0 != merged_split_peers.end()){
                        if(it1 != p_k_InConsistenct.end()){
                            int curr_decision = it1->second;
                            if(curr_decision == MERGE_BOTH_PARENTS_KIDS){
                                if(adjacentRegions(cellSegment.cell_label_maps[movieInfo.frames[movieInfo.nodes[curr_node_id].kids[0]]],
                                                   movieInfo.voxIdx[movieInfo.nodes[curr_node_id].parents[0]], movieInfo.nodes[curr_node_id].parents[1])){
                                    merged_split_peers.push_back(make_tuple(curr_node_id, MERGE_PARENTS, 1));
                                }else{
                                    merged_split_peers.push_back(make_tuple(curr_node_id, MERGE_PARENTS, 0));
                                }
                            }
                            if(curr_decision != MERGE_SPLIT_NOT_SURE){
                                continue;
                            }
                        }else if(get<1>(*it0) == SPLIT_BY_KIDS){
                            nodes4split.insert(movieInfo.nodes[curr_node_id].kids[0]);
                            nodes4split.insert(movieInfo.nodes[curr_node_id].kids[1]);
                            continue;
                        }
                    }
                }
                bool parents_already_merged = false;
                for(int p = 0; p<movieInfo.nodes[curr_node_id].parent_num; p++ ){
                    size_t p_c = movieInfo.nodes[curr_node_id].parents[p];
                    auto it = find_if(merged_split_peers.begin(), merged_split_peers.end(),
                                                 [&p_c](const tuple<size_t, int, float>& e){ return get<0>(e) == p_c;} );

                    if(it != merged_split_peers.end() && get<1>(*it) == MERGE_KIDS){
                        parents_already_merged = true;
                        break;
                    }
                }
                if(parents_already_merged) continue;

                float pvalue;
                int curr_decision = regionSplitMergeJudge(cellSegment, curr_node_id, false, pvalue);
                if (curr_decision != MERGE_SPLIT_NOT_SURE){
                    merged_split_peers.push_back(make_tuple(curr_node_id, curr_decision, pvalue));
                    if(curr_decision == SPLIT_BY_PARENTS){
                        nodes4split.insert(movieInfo.nodes[curr_node_id].kids[0]);
                        nodes4split.insert(movieInfo.nodes[curr_node_id].kids[1]);
                    }
                }

            }
        }
        //all_merge_split_peers.insert(all_merge_split_peers.end(), merged_split_peers.begin(), merged_split_peers.end());

        // deal with inconsistent decisions: not needed

        FOREACH_i(merged_split_peers){
            //separateRegion(cellSegmentMain &cellSegment, size_t node_idx,
            //size_t seeds[2], vector<simpleNodeInfo> &newCells)
            size_t curr_node_idx = get<0>(merged_split_peers[i]);
            size_t seeds[2];
            bool merge_flag = false, split_flag = false;
            switch (get<1>(merged_split_peers[i])){
            case MERGE_PARENTS:
                merge_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].parents[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].parents[1];
                break;
            case MERGE_KIDS:
                merge_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].kids[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].kids[1];
                break;
            case SPLIT_BY_PARENTS:
                split_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].parents[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].parents[1];
                break;
            case SPLIT_BY_KIDS:
                split_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].kids[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].kids[1];
                break;
            default:
                break;
            }

            if(merge_flag){
                if(mergedRegionGrow(cellSegment, seeds, newCells)){
                    uptCell_idxs.push_back(seeds[0]);
                    uptCell_idxs.push_back(seeds[1]);
                    // if success, we need to check if gap should be removed
                    if(get<2>(merged_split_peers[i]) < 0.05){ // we are confident with the merging, so remove gaps between regions
                        int frame = movieInfo.frames[seeds[0]];
                        setValMat(cellSegment.cell_label_maps[frame], CV_8U, newCells[newCells.size()-1].voxIdx, 0);
                    }
                }
            }else if(split_flag){
                if(!separateRegion(cellSegment, curr_node_idx, seeds, newCells)){
                    //curr_node_idx (1) fails to be split and (2) none of its parents can link to it (cost > obzCost)
                    // such case, we merge two parents
                    if(mergedRegionGrow(cellSegment, seeds, newCells)){
                        uptCell_idxs.push_back(seeds[0]);
                        uptCell_idxs.push_back(seeds[1]);
                        // for the un-seperable region, we do not want to frequently check them; so remove the gaps between them
                        int frame = movieInfo.frames[seeds[0]];
                        setValMat(cellSegment.cell_label_maps[frame], CV_8U, newCells[newCells.size()-1].voxIdx, 0);
                    }
                }else{
                    uptCell_idxs.push_back(curr_node_idx);
                }
            }
        }
    }

}
/**
 * @brief findBestNeighbor: find the best neighbor (linking cost) of n1 cell in a given frame
 * @param n1
 * @param target_frame
 * @return
 */
bool cellTrackingMain::findBestNeighbor(size_t n1, size_t &out_best_nei, float &cost, int target_frame){
    cost = INFINITY;
    if (target_frame > movieInfo.frames[n1]){
        for(nodeRelation n : movieInfo.nodes[n1].neighbors){
            if(movieInfo.frames[n.node_id]==target_frame && n.link_cost < cost){
                cost = n.link_cost;
                out_best_nei = n.node_id;
            }
        }
    }else{
        for(nodeRelation n : movieInfo.nodes[n1].preNeighbors){
            if(movieInfo.frames[n.node_id]==target_frame && n.link_cost < cost){
                cost = n.link_cost;
                out_best_nei = n.node_id;
            }
        }
    }
    if(cost == INFINITY){
        return false;
    }else{
        return true;
    }
}
/**
 * @brief isBestNeighbor: is n1 and n2 are best neighbor (correspond to linking cost)
 * @param n1
 * @param n2
 * @param cost
 * @return
 */
bool cellTrackingMain::isBestNeighbor(size_t n1, size_t n2, float &cost){
    size_t best_nei;
    if(findBestNeighbor(n1, best_nei, cost, movieInfo.frames[n2])){
        if (n2 != best_nei){
            return false;
        }else{
            if(findBestNeighbor(n2, best_nei, cost, movieInfo.frames[n1])){
                if(n1 != best_nei){
                    return false;
                }else{
                    return true;
                }
            }else{
                return false;
            }
        }
    }else{
        return false;
    }

}
/**
 * @brief regNodeValid: we may removed some region already
 * @param node_idx
 * @return
 */
bool cellTrackingMain::isValid(size_t node_idx, cellSegmentMain *cellSegment){
    if(node_idx >= movieInfo.nodes.size()) return false;
    if(movieInfo.nodes[node_idx].node_id != node_idx) return false;
    if(cellSegment == nullptr) return true;

    int f = movieInfo.frames[node_idx];
    if(cellSegment->cell_label_maps[f].at<int>(movieInfo.voxIdx[node_idx][0]) != movieInfo.labelInMap[node_idx]){
        return false;
    }
    return true;

}
void cellTrackingMain::infinitifyCellRelation(size_t n1, size_t n2){
    if(n1 > movieInfo.frames.size() || n2 > movieInfo.frames.size()) return;
    int f1 = movieInfo.frames[n1];
    int f2 = movieInfo.frames[n2];
    if(f1==f2) return;

    if(f1 > f2){
        size_t tmp = n1;
        n1 = n2;
        n2 = tmp;
    }
    for(nodeRelation &n : movieInfo.nodes[n1].neighbors){
        if(n.node_id == n2){
            n.link_cost = INFINITY;
            n.dist_c2n = 1000;
            n.dist_n2c = 1000;
            break;
        }
    }
    for(nodeRelation &n : movieInfo.nodes[n2].preNeighbors){
        if(n.node_id == n1){
            n.link_cost = INFINITY;
            n.dist_c2n = 1000;
            n.dist_n2c = 1000;
            break;
        }
    }
}
//!!!!!!!!!!! TODO: have not finished the implementation !!!!!!!!!!!!!!!!!!!!**/
void cellTrackingMain::addOneNewCell(cellSegmentMain &cellSegment, vector<size_t> voxIdx, int frame){
//    nodeInfo n1;
//    n1.node_id = movieInfo.frames.size() + 1;
//    n1.stable_status = NOT_STABLE;
//    movieInfo.frames.push_back(frame);
//    movieInfo.voxIdx.push_back(voxIdx);

//    vector<int> y, x, z;
//    vec_ind2sub(voxIdx, y, x, z, cellSegment.cell_label_maps[0].size);
//    movieInfo.vox_y.push_back(y);
//    movieInfo.vox_x.push_back(x);
//    movieInfo.vox_z.push_back(z);
//    cellSegment.cell_label_maps

}
/**
 * @brief nullifyCellOrNode: nullify a cell or a node from movieInfo
 * NOTE: not necessary to remove from the array at current stage
 * @param node_idx
 */
void cellTrackingMain::nullifyCellOrNode(size_t node_idx){
    movieInfo.nodes[node_idx].node_id = -1;
}
void cellTrackingMain::nullifyCellOrNode(size_t node_idx[]){
    int num = sizeof(node_idx)/sizeof(node_idx[0]);
    for(int i=0; i<num; i++){
        nullifyCellOrNode(node_idx[i]);
    }
}
void cellTrackingMain::nullifyCellOrNode(vector<size_t> node_idx){
    for(size_t n : node_idx){
        nullifyCellOrNode(n);
    }
}
/**
 * @brief separateRegion
 * @param cellSegment
 * @param node_idx
 * @param seeds
 */
bool cellTrackingMain::separateRegion(cellSegmentMain &cellSegment, size_t node_idx,
                                      size_t seeds[2], vector<simpleNodeInfo> &newCells){
    //newCells.resize(0);
    // first check that the regions is valid
    if(isValid(node_idx, &cellSegment)){
        bool split_success = false;
        vector<vector<size_t>> splitRegs(2);
        float reg4seeds2splitRes_costs[2];
        bool gapBasedSplit = false, usePriorGapMap = true;
        if(p4tracking.stableNodeTest == false || movieInfo.nodes[node_idx].stable_status == NOT_STABLE){
            vector<vector<size_t>> reg4seeds(2);
            reg4seeds[0] = movieInfo.voxIdx[seeds[0]];
            reg4seeds[1] = movieInfo.voxIdx[seeds[1]];
            int reg4seeds_frame = movieInfo.frames[seeds[0]];
            gapBasedSplit = false;
            split_success = bisectValidTest(cellSegment, node_idx, movieInfo.voxIdx[node_idx], movieInfo.frames[node_idx],
                                reg4seeds, reg4seeds_frame, gapBasedSplit, usePriorGapMap,
                                 splitRegs, reg4seeds2splitRes_costs);
            if(!split_success){
                splitRegs.clear();
                gapBasedSplit = true;
                split_success = bisectValidTest(cellSegment, node_idx, movieInfo.voxIdx[node_idx], movieInfo.frames[node_idx],
                                    reg4seeds, reg4seeds_frame, gapBasedSplit, usePriorGapMap,
                                     splitRegs, reg4seeds2splitRes_costs);
            }
        }
        if(!split_success){ // fail to split: let's check if merge is a possible solution
            if(handleNonSplitReg_link2oneSeed(node_idx, seeds)){
                return true;
            }else{
                return false;
            }
        }else{ // split success
            // 1. nullify the current region
            // nullifyCellOrNode(node_idx); // leave it to movieIinfo_update
            // 2. save the two new regions in
            simpleNodeInfo newC1, newC2;
            newCells.push_back(newC1);
            newCells.push_back(newC2);
            newCells[newCells.size() - 2].frame = movieInfo.frames[node_idx];
            newCells[newCells.size() - 2].voxIdx = splitRegs[0];
            newCells[newCells.size() - 1].frame = movieInfo.frames[node_idx];
            newCells[newCells.size() - 1].voxIdx = splitRegs[1];
            return true;
        }
    }else{
        qFatal("This region or seeds has been nullified!");
    }
}
/**
 * @brief handleNonSplitRegion:
 * @param node_idx
 * @param seeds
 */
bool cellTrackingMain::handleNonSplitReg_link2oneSeed(size_t node_idx, size_t seeds[2]){
    float curr_link_cost;
    bool one_is_good = false;
    for(int i=0; i<2; i++){
        if(isBestNeighbor(node_idx, seeds[0], curr_link_cost)){
            if(curr_link_cost < abs(p4tracking.observationCost)){
                infinitifyCellRelation(node_idx, seeds[1-i]); // remove the other node
                one_is_good = true;
            }
            break;
        }
    }
    return one_is_good;
}
/**
 * @brief mergedRegionGrow:merge two cells
 * @param cellSegment
 * @param seeds
 */
bool cellTrackingMain::mergedRegionGrow(cellSegmentMain &cellSegment, size_t seeds[2], vector<simpleNodeInfo> &newCells){
    assert(movieInfo.frames[seeds[0]] == movieInfo.frames[seeds[1]]);
    int curr_frame = movieInfo.frames[seeds[0]];
    if(p4tracking.stableNodeTest && movieInfo.nodes[seeds[0]].stable_status == NOT_STABLE &&
            movieInfo.nodes[seeds[1]].stable_status == NOT_STABLE){
        if(isValid(seeds[0], &cellSegment) && isValid(seeds[1], &cellSegment)){
            simpleNodeInfo newC;
            newCells.push_back(newC);
            newCells[newCells.size() - 1].frame = curr_frame;
            newCells[newCells.size() - 1].voxIdx = movieInfo.voxIdx[seeds[0]];
            newCells[newCells.size() - 1].voxIdx.insert(newCells[newCells.size() - 1].voxIdx.end(),
                    movieInfo.voxIdx[seeds[0]].begin(), movieInfo.voxIdx[seeds[0]].end());
            //nullifyCellOrNode(seeds); //leave it to movieInfo_update
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}



void cellTrackingMain::split_merge_module(cellSegmentMain &cellSegment){
    // 1. build the graph allowing split and merge and conduct min-cost circulation based tracking
    handleMergeSplitRegions();
    // 2. for each one2two or two2one linking, choose to split or merge the correspond nodes
    vector<simpleNodeInfo> newCells;
    vector<size_t> uptCell_idxs;
    regionRefresh(cellSegment, newCells, uptCell_idxs);

    // 3. update movieInfo/cellSegment based on the newly detected cells or removed cells

}
/**
 * @brief movieInfo_update:update movieInfo/cellSegment based on the newly detected cells or removed cells
 * @param cellSegment
 * @param newCells
 * @param uptCell_idxs
 */
void cellTrackingMain::movieInfo_update(cellSegmentMain &cellSegment, vector<simpleNodeInfo> &newCells,
                                        vector<size_t> &uptCell_idxs){
    if(vec_unique(uptCell_idxs)){
        qFatal("There is duplicated updating! check first func: regionRefresh() first");
    }
    vector<unordered_set<size_t>> uptCell_framewise (cellSegment.number_cells.size());
    //vector<unordered_set<size_t>> newCells_framewise (cellSegment.number_cells.size());
    FOREACH_i(uptCell_idxs){
        int frame = movieInfo.frames[uptCell_idxs[i]];
        uptCell_framewise[frame].insert(uptCell_idxs[i]);
    }
//    FOREACH_i(newCells){
//    }

    if(p4tracking.stableNodeTest){

    }
    // TODO: there should be an update related to missing cell re-detect (named: node_tested_st_end_jump)


}
void cellTrackingMain::missing_cell_module(cellSegmentMain &cellSegment){

}
