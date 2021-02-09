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

    movieInfo.frame_shift_xyz.resize(cellSegment.number_cells.size());
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i] = {0, 0, 0};
    }
    cumulative_cell_nums.resize(cellSegment.number_cells.size());
    cumulative_cell_nums[0] = cellSegment.number_cells[0];
    for (int i = 1; i< cumulative_cell_nums.size(); i++){
        cumulative_cell_nums[i] = (cumulative_cell_nums[i-1] + cellSegment.number_cells[i]);
    }
    movieInfo.tracks.resize(0);
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
        vector<size_t> nei_ids;
        extractNeighborIds(cellSegment.cell_label_maps, i, nei_ids);
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
    mccTracker({}, {}); // with no split/merge vectors
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
 * @brief cellInfo2graph: using the movieInfo to build the graph for CINDA
 */
void cellTrackingMain::mccTracker(vector<pair<size_t, float>> merge_node_idx, vector<pair<size_t, float>> split_node_idx){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      We directly use the interface I wrote for python with the following function name               //
    // long long* pyCS2(long *msz, double *mtail, double *mhead, double *mlow, double *macap, double *mcost)//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long m_append = merge_node_idx.size() + split_node_idx.size();
    // build the graph without split and merge settings
    long long n = 2*movieInfo.nodes.size() + 1; // # vertices in graph
    long long m = movieInfo.overall_neighbor_num + 3*movieInfo.nodes.size() + m_append; // # arcs in graph (upper bound)
    double *mtail = new double[m];
    double *mhead = new double[m];
    double *mlow = new double[m]; memset(mlow, 0, m);
    double *macap = new double[m]; memset(macap, 1, m);
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 1;
    int rand_max = 100;
    bool link_adj_frs = false;
    if(p4tracking.splitMergeHandle == NOJUMPALL && m_append > 0){
        link_adj_frs = true;
    }
    for(nodeInfo node : movieInfo.nodes){
        // in arc
        mtail[arc_cnt] = src_id;
        mhead[arc_cnt] =  2 * node.node_id;
        mcost[arc_cnt] = round(node.in_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // out arc
        mhead[arc_cnt] = src_id;
        mtail[arc_cnt] =  2 * node.node_id;
        mcost[arc_cnt] = round(node.out_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // link with neighbors
        if(link_adj_frs){
            int curr_frame = movieInfo.frames[node.node_id];
            for(nodeRelation neighbor : node.neighbors){
                if(movieInfo.frames[neighbor.node_id]-1 == curr_frame){
                    mtail[arc_cnt] = 2 * node.node_id;
                    mhead[arc_cnt] =  2 * neighbor.node_id;
                    mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                    arc_cnt ++;
                }
            }
        }else{
            for(nodeRelation neighbor : node.neighbors){
                mtail[arc_cnt] = 2 * node.node_id;
                mhead[arc_cnt] =  2 * neighbor.node_id;
                mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
                arc_cnt ++;
            }
        }
    }
    // add the arcs of splitting and merging
    if(m_append > 0){
        for(pair<size_t, float> node_id_cost : merge_node_idx){
            mtail[arc_cnt] = 2 * node_id_cost.first;
            mhead[arc_cnt] =  src_id;
            mcost[arc_cnt] = round(node_id_cost.second * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
            arc_cnt ++;
        }

        for(pair<size_t, float> node_id_cost : split_node_idx){
            mtail[arc_cnt] = src_id;
            mhead[arc_cnt] = 2 * node_id_cost.first + 1;
            mcost[arc_cnt] = round(node_id_cost.second * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
            arc_cnt ++;
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
    // update the jumpCost in p4tracking if no split/merge allowed
    if (m_append == 0){ // by the way, also update stable segmentations
        updateJumpCost();
        stableSegmentFixed();
        if (false){
            //re-link jump-over cells
        }
    }

    delete[] mtail;
    delete[] mhead;
    delete[] mlow;
    delete[] macap;
    delete[] mcost;
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
float cellTrackingMain::bestPeerCandidate(size_t node_id, vector<size_t> &bestPeer, bool parent_flag){
    vector<size_t> peerCandidates(0);
    bestPeer.resize(0);
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
        bestPeer = peerCandidates;
        float dummy_c2n, dummy_n2c;
        return voxelwise_avg_distance(node_id, bestPeer, dummy_c2n, dummy_n2c);
    }else{ // more than two regions available
        float min_distance = INFINITY, curr_distance;

        FOREACH_i(peerCandidates){
            bestPeer[0] = peerCandidates
            for(size_t j = i + 1; j < peerCandidates.size(); j ++){
                curr_distance = voxelwise_avg_distance(())
            }
        }
    }
}
/**
 * @brief detectPeerRegions: test if two regions should be merged based on their common kid or parent region
 * @param merge_node_idx
 * @param split_node_idx
 */
void cellTrackingMain::detectPeerRegions(vector<splitMergeNodeInfo> &split_merge_node_info){
    for(nodeInfo node : movieInfo.nodes){
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

        }
        if(kids_test){

        }
    }
}
/**
 * @brief handleMergeSplitRegions: build the network allowing split/merge
 */
void cellTrackingMain::handleMergeSplitRegions(){
    vector<splitMergeNodeInfo> split_merge_node_info;
    detectPeerRegions(split_merge_node_info);
}
void cellTrackingMain::split_merge_module(cellSegmentMain &cellSegment){
    handleMergeSplitRegions();
}

void cellTrackingMain::missing_cell_module(cellSegmentMain &cellSegment){

}
