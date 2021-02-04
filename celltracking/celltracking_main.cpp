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
    float max_dist = distanceTransRegion2Region(ref_cell, movieInfo.range_xyz[cell_curr],
                               mov_cell, movieInfo.range_xyz[cell_nei],
                               shift_xyz, dist);
    c2n = dist[0];
    n2c = dist[1];
    return max_dist;
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

    // calculate the distances between each cell and its neighbors
    vector<float> nn_dist;
    calCell2neighborDistance(nn_dist);
    truncatedGammafit(nn_dist, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
    updateArcCost(); // calculate the link cost from given new gamma parameter
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
void cellTrackingMain::cellInfo2graph(vector<pair<size_t, float>> merge_node_idx = {}, vector<pair<size_t, float>> split_node_idx = {}){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      We directly use the interface I wrote for python with the following function name               //
    // long long* pyCS2(long *msz, double *mtail, double *mhead, double *mlow, double *macap, double *mcost)//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long m_append = merge_node_idx.size() + split_node_idx.size();
    // build the graph without split and merge settings
    long long n = 2*movieInfo.nodes.size() + 1; // # vertices in graph
    long long m = movieInfo.overall_neighbor_num + 3*movieInfo.nodes.size() + m_append; // # arcs in graph
    double *mtail = new double[m];
    double *mhead = new double[m];
    double *mlow = new double[m]; memset(mlow, 0, m);
    double *macap = new double[m]; memset(macap, 1, m);
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 1;
    int rand_max = 100;
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
        for(nodeRelation neighbor : node.neighbors){
            mtail[arc_cnt] = 2 * node.node_id;
            mhead[arc_cnt] =  2 * neighbor.node_id;
            mcost[arc_cnt] = round(neighbor.link_cost * 1e7 + rand() % rand_max); // add some randomness to make sure unique optimal solution
            arc_cnt ++;
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
    assert(arc_cnt == m);
    long *msz = new long[3]; // should be long long, but for our data, long is enough
    msz[0] = 12;
    msz[1] = n;
    msz[2] = m;

    // call pyCS2
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
            movieInfo.tracks.push_back(new_track);
        }
    }
}


// functions to update cost given new gamma fitting results
void cellTrackingMain::driftCorrection(){

}
void cellTrackingMain::updateGammaParam(){

}
void cellTrackingMain::updateArcCost(){
    for(nodeInfo n : movieInfo.nodes){
        for(nodeRelation &neighbor : n.neighbors){
            neighbor.link_cost = distance2cost(MAX(neighbor.dist_c2n, neighbor.dist_n2c),
                                               movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
        }
    }
}
void cellTrackingMain::stableArcFixed(){

}
