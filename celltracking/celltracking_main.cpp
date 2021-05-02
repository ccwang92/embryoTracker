#include "celltracking_main.h"
#include <stack>          // std::stack
#include <deque>
#include <fstream> // for file stream
#include <algorithm>    // std::fill
#include "CINDA/src_c/cinda_funcs.c"
using namespace std;
using namespace cv;

// inline functions for crop detections
inline size_t key(int batch, int time,int section, int label) {
    return (size_t) label << 32 | (unsigned int) ((time << 6) + (batch << 2) + section);
}
inline void dekey(size_t k, int &batch, int &time,int &section, int &label) {
    label = k >> 32;
    int ts = k & 0xFFFF;
    section = ts & 3;
    batch = (ts & 63) >> 2;
    time = ts >> 6;
}
// for fused detections
inline size_t newkey(int time,int label) {
    return (size_t) label << 32 | (unsigned int) time;
}
inline void denewkey(size_t k, int &time, int &label) {
    label = k >> 32;
    time = k & 0xFFFF;
}


cellTrackingMain::cellTrackingMain(cellSegmentMain &cellSegment, const QStringList &fileNames, bool _debugMode)
{
    debugMode = _debugMode;
    tracking_sucess =false;
    /////////////////////////////////////////////////
    //   step 1. initialize cell info              //
    /////////////////////////////////////////////////
    cellInfoAccumuate(cellSegment);
    init_parameter(cellSegment.p4segVol, movieInfo.frames.size()/2);
    /////////////////////////////////////////////////
    //   step 2. try a naive tracking first        //
    /////////////////////////////////////////////////
    initTransitionCost(cellSegment);
    /////////////////////////////////////////////////
    //   step 3. main loop for cell tracking       //
    /////////////////////////////////////////////////
    int loop_cnt = 1;
    while (loop_cnt <= p4tracking.maxIter){
        //missing_cell_module(cellSegment);
        // MODULE ONE: split/merge test from current tracking results
        for(int dummy = 0; dummy < 3; dummy ++){
            split_merge_module(cellSegment);
        }
        // MODULE TWO: retrieve missing cells
        missing_cell_module(cellSegment);
        loop_cnt ++;
    }
    //missing_cell_module(cellSegment);//split_merge_module(cellSegment);
    /////////////////////////////////////////////////
    //   step 4. merge broken tracks               //
    /////////////////////////////////////////////////
    mccTracker_one2one();
    //merge_broken_tracks();
    //mergeOvTracks2(); // wrap up trackes using parent-kid relation (may contain redundant computation)
    tracking_sucess = true;
    size_t total_cells = accumulate(cellSegment.number_cells.begin(), cellSegment.number_cells.end(), 0);
    qInfo("----------------%ld cells after iterative correction-------------------", total_cells);
    saveTrackAndFinalSegmentResults(cellSegment, fileNames);
}
cellTrackingMain::cellTrackingMain(vector<int> data_size_yxzt, const QStringList &fileNames){
    tracking_sucess =false;
    if(!loadTrackResults(data_size_yxzt, fileNames)){
        qDebug("Fail to load files");
    }
    bool get_res_from_txt = true, get_jumpCost_only = false;
    mccTracker_one2one(get_jumpCost_only, get_res_from_txt);
    tracking_sucess = true;
}
bool cellTrackingMain::saveTrackAndFinalSegmentResults(cellSegmentMain &cellSegment, const QStringList &fileNames){
    //save segment results
    for (int i=0; i<fileNames.size(); i++){
        QString fileName = fileNames.at(i);
        string fileNameNoExt = fileName.left(fileName.lastIndexOf('.')).toStdString();
        string label_file_name = fileNameNoExt + "_label_map_int32_final.bin";
        ofstream label_file(label_file_name, ios::binary);
        if (!label_file.is_open()) return false;
        // tmp is local variable, which will be released soon, so we need copyTo
        label_file.write((const char*)(cellSegment.cell_label_maps[i].data),
                         cellSegment.cell_label_maps[i].elemSize() * cellSegment.cell_label_maps[i].total());
        label_file.close();
    }
    // write neighbor infor
    QString fileName = fileNames.at(0);
    string fileNameNoExt = fileName.left(fileName.lastIndexOf('/')).toStdString();
    string movieInfo_txt_name = fileNameNoExt + "/movieInfo.txt";
    ofstream movieInfo_txt(movieInfo_txt_name);
    movieInfo_txt << "p cell num: " << movieInfo.nodes.size() << "\n";
    movieInfo_txt << "p gamma: " << movieInfo.ovGammaParam[0] << " " << movieInfo.ovGammaParam[1] << "\n";
    movieInfo_txt << "p jumpRatio: " << p4tracking.jumpCost[0] << " "<< p4tracking.jumpCost[1] << " "<< p4tracking.jumpCost[2] << "\n";
    movieInfo_txt << "p enexObz: " << p4tracking.c_en << " "<< p4tracking.c_ex << " "<< p4tracking.observationCost << "\n";
    movieInfo_txt << "c node info: id, frame, label_in_frame \n";
    FOREACH_i(movieInfo.nodes){
        movieInfo_txt <<"n "<< i <<" "<< movieInfo.frames[i] << " " << movieInfo.labelInMap[i] << "\n";
    }
    movieInfo_txt << "c link info: id1, id2, distance, link_cost \n";
    long long arc_cnt = 0;
    float linkCostUpBound = INFINITY;
    for(size_t i=0; i < movieInfo.nodes.size(); i++){
        nodeInfo &node = movieInfo.nodes[i];
        if(movieInfo.voxIdx[i].size() == 0) continue;
        // link with neighbors
        for(nodeRelation neighbor : node.neighbors){
            if(neighbor.link_cost < linkCostUpBound){
                movieInfo_txt <<"a "<< node.node_id <<" "<< neighbor.node_id << " "
                             << neighbor.dist_c2n << " " << neighbor.link_cost << "\n";
                arc_cnt ++;
            }
        }
    }
    movieInfo_txt.close();
    return true;
}
bool cellTrackingMain::loadTrackResults(vector<int> data_size_yxzt, const QStringList &fileNames){
    ///-------------------------read movieInfo information from txt file-----------------------//
    QString fileName = fileNames.at(0);
    QString movieInfo_txt_name = fileName.left(fileName.lastIndexOf('/')) + "/movieInfo.txt";
    QFile txt_file(movieInfo_txt_name);
    txt_file.open(QFile::Text | QFile::ReadOnly);
    QTextStream in(&txt_file);
    long cell_num;
    //QString tmp = in.readLine();
    sscanf(in.readLine().toStdString().c_str(), "p cell num: %ld", &cell_num);
    for(int i=0; i<2; i++){
        // skip the first 2 lines
        //            p gamma: 1.50283 8.1355
        //            p jumpRatio: 0.691873 0.162243 0.145884
        in.readLine();
    }
    movieInfo.nodes.resize(cell_num);
    movieInfo.frames.resize(cell_num);
    movieInfo.xCoord.resize(cell_num);
    movieInfo.yCoord.resize(cell_num);
    movieInfo.zCoord.resize(cell_num);

    double en_cost, ex_cost, boz_cost;
    sscanf(in.readLine().toStdString().c_str(), "p enexObz: %lf %lf %lf", &en_cost, &ex_cost, &boz_cost);
    p4tracking.c_en = en_cost;// cost of appearance and disappearance in the scene
    p4tracking.c_ex = ex_cost;
    p4tracking.observationCost = -(p4tracking.c_en+p4tracking.c_ex) + 0.00001; // make sure detections are all included

    FOREACH_i(movieInfo.nodes){
        movieInfo.nodes[i].node_id = i;
        movieInfo.nodes[i].in_cost = p4tracking.c_en;
        movieInfo.nodes[i].out_cost = p4tracking.c_ex;
    }
    unordered_map<size_t, size_t> frame_label2cell_id;
    movieInfo.overall_neighbor_num = 0;
    while(!in.atEnd()) {
        string line = in.readLine().toStdString();
        // start parsing each line
        switch (line[0]) {
        case 'c':                  /* skip lines with comments */
        case '\n':                 /* skip empty lines   */
        case 'p':
        case '\0':                 /* skip empty lines at the end of file */
            break;
        case 'n':{
            int id, frame, labelInMap;
            sscanf(line.c_str(), "%*c %d %d %d", &id, &frame, &labelInMap);
            if(labelInMap != 0){
                frame_label2cell_id[newkey(frame, labelInMap)] = id;
                movieInfo.frames[id] = frame;
            }
            break;
        }
        case 'a':{
            int tail = 0;
            int head = 0;
            double distance, link_cost;
            sscanf(line.c_str(), "%*c %d %d %lf %lf", &tail, &head, &distance, &link_cost);

            nodeRelation tmp;
            tmp.node_id = head;
            tmp.dist_c2n = distance;
            tmp.link_cost = link_cost;
            movieInfo.nodes[tail].neighbors.push_back(tmp);
            movieInfo.overall_neighbor_num ++;
            break;
        }
        default:
            break;
        }
    }
    txt_file.close();
    ///----------------------------------read node location information from .bin files--------------//
    for (int i=0; i<fileNames.size(); i++){
        QString fileName = fileNames.at(i);
        QString label_file_name = fileName.left(fileName.lastIndexOf('.')) + "_label_map_int32_final.bin";
        QFile label_file(label_file_name);
        if (!label_file.open(QIODevice::ReadOnly)){
            qFatal("lost files!");
        }
        Mat mat_cur_time;
        Mat(3, data_size_yxzt.data(), CV_32S, label_file.readAll().data()).copyTo(mat_cur_time);
        label_file.close();

        double max_id;
        minMaxIdx(mat_cur_time, nullptr, &max_id);
        vector<vector<size_t>> cell_voxIdx;
        extractVoxIdxList(&mat_cur_time, cell_voxIdx, (int)max_id);

        FOREACH_j(cell_voxIdx){
            if(cell_voxIdx[j].size() == 0) continue;
            vector<int> y, x, z;
            vec_ind2sub(cell_voxIdx[j], y, x, z, mat_cur_time.size);
            size_t tmp = frame_label2cell_id[newkey(i, (int)j+1)];
            //extractSumGivenIdx(mat_cur_time, vector<size_t> idx, int datatype)
//            if (movieInfo.xCoord[tmp] > 330) {
//                qDebug("check point");
//            }
            movieInfo.xCoord[tmp] = vec_mean(x);
            movieInfo.yCoord[tmp] = vec_mean(y);
            movieInfo.zCoord[tmp] = vec_mean(z);
        }
    }
    return true;
}

void cellTrackingMain::extractTraceLocations(vector<int> data_size_yxzt, int width){
    //vector<vector<unordered_set<size_t>>> trace_sets; // time x trackes
    trace_sets.resize(data_size_yxzt[3]);
    for(int t=0; t<data_size_yxzt[3]; t++){
        trace_sets[t].resize(movieInfo.tracks.size());
    }
    vector<int> yxz_sz = {data_size_yxzt[0], data_size_yxzt[1], data_size_yxzt[2]};
    for(size_t i=0; i<movieInfo.tracks.size(); i++){
        //qInfo("%ld", i);
        for(size_t j=0; j<movieInfo.tracks[i].size(); j++){ // start from the second point
            //qInfo("%ld-%ld", i, j);
            if(j==0){
                size_t curr = movieInfo.tracks[i][j];
                vector<float> end_yxz = {movieInfo.yCoord[curr], movieInfo.xCoord[curr], movieInfo.zCoord[curr]};
                //qInfo("%f-%f-%f", end_yxz[0], end_yxz[1], end_yxz[2]);
                traceExtract(end_yxz, end_yxz, yxz_sz, width, trace_sets[movieInfo.frames[curr]][i]);
            }else{
                size_t pre = movieInfo.tracks[i][j-1];
                size_t curr = movieInfo.tracks[i][j];
                vector<float> start_yxz = {movieInfo.yCoord[pre], movieInfo.xCoord[pre], movieInfo.zCoord[pre]};
                vector<float> end_yxz = {movieInfo.yCoord[curr], movieInfo.xCoord[curr], movieInfo.zCoord[curr]};
                vector<float> diff = vec_Minus(end_yxz, start_yxz);
                for(int e=0; e<3;e++){
                    diff[e] /= (movieInfo.frames[curr]-movieInfo.frames[pre]);
                }
                for(int t = movieInfo.frames[pre]+1; t <= movieInfo.frames[curr]; t++){
                    for(int e=0; e<3;e++){
                        end_yxz[e] = diff[e]*(t-movieInfo.frames[pre]) + start_yxz[e];
                    }
                    traceExtract(start_yxz, end_yxz, yxz_sz, width, trace_sets[t][i]);
                }
            }
        }
    }
}
/**
 * @brief merge_broken_tracks: for node-neighbor pair with cost of dist_c2n or dist_n2c smaller than obz_cost, link them
 */
void cellTrackingMain::merge_broken_tracks(){

    for (nodeInfo & nf : movieInfo.nodes){
        if(nf.kid_num > 0) continue;
        for(nodeRelation &nr : nf.neighbors){
            float best_cost = 0;
            if (movieInfo.nodes[nr.node_id].parent_num == 0 &&
                    isBestNeighbor(nf.node_id, nr.node_id, best_cost)){
                int f_diff = movieInfo.frames[nr.node_id] - movieInfo.frames[nf.node_id];
                best_cost = distance2cost(MIN(nr.dist_n2c, nr.dist_c2n), p4tracking.jumpCost[f_diff-1]);
                if(best_cost < abs(p4tracking.observationCost)){
                    nf.kids[nf.kid_num] = nr.node_id;
                    nf.kid_cost[nf.kid_num] = best_cost;
                    nf.kid_num ++;

                    movieInfo.nodes[nr.node_id].parents[0] = nf.node_id;
                    movieInfo.nodes[nr.node_id].parent_cost[0] = best_cost;
                    movieInfo.nodes[nr.node_id].parent_num ++;
                }
                break;
            }
        }
    }
}
void cellTrackingMain::cellInfoAccumuate(cellSegmentMain &cellSegment){
    vector<vector<size_t>> curr_cell_idx;
    extractVoxIdxList(cellSegment.cell_label_maps, curr_cell_idx, cellSegment.number_cells);
    /** KEY: we generate an overly large voxIdx in movieInfo */
    movieInfo.voxIdx.resize(2*curr_cell_idx.size()); // we assume each cell at most can be split once
    long long increment = 0;
    long long cell_cnt = 0;
    for(int i : cellSegment.number_cells){
        for(size_t j=increment; j<(increment+i); j++){
            movieInfo.voxIdx[j] = curr_cell_idx[cell_cnt++];
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
    movieInfo.range_xyz.resize(movieInfo.voxIdx.size(), vector<int>(3));
    movieInfo.start_coord_xyz.resize(movieInfo.voxIdx.size(), vector<int>(3));
    movieInfo.vox_x.resize(movieInfo.voxIdx.size());
    movieInfo.vox_y.resize(movieInfo.voxIdx.size());
    movieInfo.vox_z.resize(movieInfo.voxIdx.size());
    FOREACH_i(movieInfo.voxIdx){
        if(movieInfo.voxIdx[i].size()>0){
            vec_ind2sub(movieInfo.voxIdx[i], movieInfo.vox_y[i],
                        movieInfo.vox_x[i], movieInfo.vox_z[i], cellSegment.cell_label_maps[0].size);
            movieInfo.xCoord[i] = vec_mean(movieInfo.vox_x[i]);
            movieInfo.yCoord[i] = vec_mean(movieInfo.vox_y[i]);
            movieInfo.zCoord[i] = vec_mean(movieInfo.vox_z[i]);

            //movieInfo.start_coord_xyz[i].resize(3);
            //size_t dummy;
            movieInfo.start_coord_xyz[i][0] = vec_min(movieInfo.vox_x[i]);
            movieInfo.start_coord_xyz[i][1] = vec_min(movieInfo.vox_y[i]);
            movieInfo.start_coord_xyz[i][2] = vec_min(movieInfo.vox_z[i]);
            //movieInfo.range_xyz[i].resize(3);
            movieInfo.range_xyz[i][0] = vec_max(movieInfo.vox_x[i]) - movieInfo.start_coord_xyz[i][0] + 1;
            movieInfo.range_xyz[i][1] = vec_max(movieInfo.vox_y[i]) - movieInfo.start_coord_xyz[i][1] + 1;
            movieInfo.range_xyz[i][2] = vec_max(movieInfo.vox_z[i]) - movieInfo.start_coord_xyz[i][2] + 1;
        }else{
            movieInfo.vox_y[i].reserve(0);
            movieInfo.vox_x[i].reserve(0);
            movieInfo.vox_z[i].reserve(0);
            movieInfo.xCoord[i] = -INFINITY;
            movieInfo.yCoord[i] = -INFINITY;
            movieInfo.zCoord[i] = -INFINITY;
            //movieInfo.start_coord_xyz[i].resize(3);
            movieInfo.start_coord_xyz[i] = {-1,-1,-1};
            //movieInfo.range_xyz[i].resize(3);
            movieInfo.range_xyz[i] = {0,0,0};
        }
    }
    movieInfo.nodes.resize(movieInfo.xCoord.size());
    movieInfo.frames.resize(movieInfo.xCoord.size());
    movieInfo.labelInMap.resize(movieInfo.xCoord.size());
    increment = 0;
    FOREACH_i(cumulative_cell_nums){
        for(int j = increment; j < cumulative_cell_nums[i]; j++){
            movieInfo.frames[j] = i;
            movieInfo.labelInMap[j] = (j-increment) < cellSegment.number_cells[i] ? (j-increment+1):0;
            movieInfo.nodes[j].node_id = j;
        }
        increment = cumulative_cell_nums[i];
    }

    movieInfo.frame_shift_xyz.resize(cellSegment.number_cells.size(), vector<double>(3));
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i] = {0, 0, 0};
    }

    movieInfo.tracks.resize(0);
    validGapMaps.reserve(cellSegment.number_cells.size());
    for(int i = 0; i < cellSegment.number_cells.size(); i++){
        Mat1b tmp = Mat(cellSegment.cell_label_maps[i].dims, cellSegment.cell_label_maps[i].size,
                           CV_8U, Scalar(255));
        validGapMaps.push_back(tmp);
    }
    movieInfo.node_tested_st_end_jump.resize(movieInfo.xCoord.size());
    for(missingCellTestPassed &mct : movieInfo.node_tested_st_end_jump){
        mct.jump_tested = false;
        mct.track_head_tested = false;
        mct.track_tail_tested = false;
        mct.region_id_jumped_to = 0;
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
    bool *ref_cell = new bool[curr_range_xyz[0] * curr_range_xyz[1] * curr_range_xyz[2]];
    memset(ref_cell, false, curr_range_xyz[0] * curr_range_xyz[1] * curr_range_xyz[2]);
    size_t idx;
    size_t frame_sz = curr_range_xyz[0] * curr_range_xyz[1];
    for(size_t i = 0; i < curr_voxIdx.size(); i++){
        idx = (curr_vox_z[i] - curr_start_coord_xyz[2]) * frame_sz +
               (curr_vox_x[i] - curr_start_coord_xyz[0]) * curr_range_xyz[1] +
                curr_vox_y[i] - curr_start_coord_xyz[1];
        if(curr_range_xyz[0] * curr_range_xyz[1] * curr_range_xyz[2] <= idx || idx < 0){
            qFatal("Leaking memory");
        }
        ref_cell[idx] = true;
    }
    bool *mov_cell = new bool[nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2]];
    memset(mov_cell, false, nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2]);
    frame_sz = nei_range_xyz[0] * nei_range_xyz[1];
    for(size_t i = 0; i < nei_voxIdx.size(); i++){
        idx = (nei_vox_z[i] - nei_start_coord_xyz[2]) * frame_sz +
               (nei_vox_x[i] - nei_start_coord_xyz[0]) * nei_range_xyz[1] +
                nei_vox_y[i] - nei_start_coord_xyz[1];
        if(nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2] <= idx || idx < 0){
            qFatal("Leaking memory");
        }
        mov_cell[idx] = true;  // !!! NOTE: here may cause memory leaking
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

    delete[] ref_cell;
    delete[] mov_cell;
    return max_dist;
}
/** debug Done
 * @brief voxelwise_avg_distance
 * @param curr_voxIdx
 * @param curr_frame
 * @param nei_voxIdx
 * @param nei_frame
 * @param data3d_sz
 * @param c2n
 * @param n2c
 * @return
 */
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


    bool *ref_cell = new bool[curr_range_xyz[0] * curr_range_xyz[1] * curr_range_xyz[2]];
    memset(ref_cell, false, curr_range_xyz[0] * curr_range_xyz[1] * curr_range_xyz[2]);
    size_t idx;
    size_t frame_sz = curr_range_xyz[0] * curr_range_xyz[1];
    for(size_t i = 0; i < curr_voxIdx.size(); i++){
        idx = (curr_vox_z[i] - curr_start_coord_xyz[2]) * frame_sz +
               (curr_vox_x[i] - curr_start_coord_xyz[0]) * curr_range_xyz[1] +
                curr_vox_y[i] - curr_start_coord_xyz[1];

        ref_cell[idx] = true;
    }
    bool *mov_cell = new bool[nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2]];
    memset(mov_cell, false, nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2]);
    frame_sz = nei_range_xyz[0] * nei_range_xyz[1];
    for(size_t i = 0; i < nei_voxIdx.size(); i++){
        idx = (nei_vox_z[i] - nei_start_coord_xyz[2]) * frame_sz +
               (nei_vox_x[i] - nei_start_coord_xyz[0]) * nei_range_xyz[1] +
                nei_vox_y[i] - nei_start_coord_xyz[1];
        if(nei_range_xyz[0] * nei_range_xyz[1] * nei_range_xyz[2] <= idx || idx < 0){
            qFatal("Leaking memory");
        }
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
    delete[] mov_cell;
    delete[] ref_cell;
    return max_dist;
}

float cellTrackingMain::voxelwise_avg_distance(vector<size_t> &curr_voxIdx, int curr_frame,
                             vector<size_t> &nei_voxIdx, int nei_frame,
                             MatSize data3d_sz, float &c2n, float &n2c){
    int data_sz[3] = {data3d_sz[0], data3d_sz[1], data3d_sz[2]};
    return voxelwise_avg_distance(curr_voxIdx, curr_frame, nei_voxIdx, nei_frame, data_sz, c2n, n2c);
}
/** debug done
 * @brief combineCellsIntoOneRegion: make a region by joining several cells
 * @param cell_idxes
 * @param out_region_info
 */
void cellTrackingMain::combineCellsIntoOneRegion(vector<size_t> &cell_idxes, combinedCellsCensus &out_region_info){
    out_region_info.frame = movieInfo.frames[cell_idxes[0]];
    out_region_info.start_coord_xyz.resize(3);
    out_region_info.range_xyz.resize(3);
    fill(out_region_info.range_xyz.begin(), out_region_info.range_xyz.end(), -1);

    for(size_t n : cell_idxes){
        if(movieInfo.voxIdx[n].size() == 0) continue;
        out_region_info.voxIdx.insert(out_region_info.voxIdx.end(), movieInfo.voxIdx[n].begin(), movieInfo.voxIdx[n].end());
        out_region_info.vox_x.insert(out_region_info.vox_x.end(), movieInfo.vox_x[n].begin(), movieInfo.vox_x[n].end());
        out_region_info.vox_y.insert(out_region_info.vox_y.end(), movieInfo.vox_y[n].begin(), movieInfo.vox_y[n].end());
        out_region_info.vox_z.insert(out_region_info.vox_z.end(), movieInfo.vox_z[n].begin(), movieInfo.vox_z[n].end());

        if(out_region_info.range_xyz[0] == -1){
            out_region_info.range_xyz = movieInfo.range_xyz[n];
            out_region_info.start_coord_xyz =  movieInfo.start_coord_xyz[n];
        }else{
            for (int i = 0; i < 3; i++){
                if(movieInfo.start_coord_xyz[n][i] < out_region_info.start_coord_xyz[i]){
                    out_region_info.range_xyz[i] += (out_region_info.start_coord_xyz[i] - movieInfo.start_coord_xyz[n][i]);
                    out_region_info.start_coord_xyz[i] = movieInfo.start_coord_xyz[n][i];
                }
                if ((movieInfo.start_coord_xyz[n][i] + movieInfo.range_xyz[n][i]) >
                        (out_region_info.start_coord_xyz[i] + out_region_info.range_xyz[i])){
                    out_region_info.range_xyz[i] = (movieInfo.start_coord_xyz[n][i] + movieInfo.range_xyz[n][i] - 1) -
                            out_region_info.start_coord_xyz[i] + 1;
                }
            }
//            if(movieInfo.start_coord_xyz[n][1] < out_region_info.start_coord_xyz[1]){
//                out_region_info.range_xyz[1] += (out_region_info.start_coord_xyz[1] - movieInfo.start_coord_xyz[n][1]);
//                out_region_info.start_coord_xyz[1] = movieInfo.start_coord_xyz[n][1];
//            }
//            if ((movieInfo.start_coord_xyz[n][1] + movieInfo.range_xyz[n][1]) >
//                    (out_region_info.start_coord_xyz[1] + out_region_info.range_xyz[1])){
//                out_region_info.range_xyz[1] = (movieInfo.start_coord_xyz[n][1] + movieInfo.range_xyz[n][1]) -
//                        out_region_info.range_xyz[1] + 1;
//            }
//            if(movieInfo.start_coord_xyz[n][2] < out_region_info.start_coord_xyz[2]){
//                out_region_info.range_xyz[2] += (out_region_info.start_coord_xyz[2] - movieInfo.start_coord_xyz[n][2]);
//                out_region_info.start_coord_xyz[2] = movieInfo.start_coord_xyz[n][2];
//            }
        }
    }
}

/**
 * @brief initPreNeighborInfo: initialize the information of candidate parents of a node
 */
void cellTrackingMain::initPreNeighborInfo(){
    //// update the preNeighbors information
    for (nodeInfo n : movieInfo.nodes){
        for (nodeRelation neighbor : n.neighbors){
            nodeRelation tmp;
            tmp.node_id = n.node_id;
            tmp.dist_c2n = neighbor.dist_n2c;
            tmp.dist_n2c = neighbor.dist_c2n;
            tmp.link_cost = neighbor.link_cost;
            tmp.overlap_size = neighbor.overlap_size;
            movieInfo.nodes[neighbor.node_id].preNeighbors.push_back(tmp); //pre-nei has no dist_c2n or dist_n2c
        }
    }
}

/**
 * @brief updatePreNeighborInfo: update the information of candidate parents of a node
 */
void cellTrackingMain::updatePreNeighborInfo(bool link_cost_only = true){
    //// update the preNeighbors information
    if(link_cost_only){
        for (nodeInfo n : movieInfo.nodes){
            for (nodeRelation neighbor : n.neighbors){
                for(nodeRelation &preNei : movieInfo.nodes[neighbor.node_id].preNeighbors){
                    if (preNei.node_id == n.node_id){
                        preNei.link_cost = neighbor.link_cost;
                        break;
                    }
                }
            }
        }
    }else{
        for (nodeInfo n : movieInfo.nodes){
            for (nodeRelation neighbor : n.neighbors){
                for(nodeRelation &preNei : movieInfo.nodes[neighbor.node_id].preNeighbors){
                    if (preNei.node_id == n.node_id){
                        preNei.dist_c2n = neighbor.dist_n2c;
                        preNei.dist_n2c = neighbor.dist_c2n;
                        preNei.link_cost = neighbor.link_cost;
                        break;
                    }
                }
            }
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
        //qDebug("val: %d",cellSegment.cell_label_maps[curr_frame].at<int>(k));
    }
    return overlap_sz;
}
/** exactly the same as Matlab version
 * @brief distance2cost
 * @param distance
 * @param alpha
 * @param beta
 * @param punish
 * @return
 */
float cellTrackingMain::distance2cost(float distance, float punish=1){
    //assert(frame_diff>0 && frame_diff)
    if (distance > 100) return INFINITY;
    float alpha = movieInfo.ovGammaParam[0];
    float beta = movieInfo.ovGammaParam[1];
    //float punish = p4tracking.jumpCost[frame_diff-1];
    double p = gammacdf(distance, alpha, beta);
    double cost = normInv(p * punish / 2);
//    if(isinf(cost*cost)){
//        qDebug("check point");
//    }
    return (float)(cost*cost);
}

void cellTrackingMain::reCalculateCellsDistances(){
    FOREACH_i(movieInfo.nodes){
        FOREACH_j(movieInfo.nodes[i].neighbors){
            // cal the distance between two nodes
            voxelwise_avg_distance(i, movieInfo.nodes[i].neighbors[j].node_id,
                                              movieInfo.nodes[i].neighbors[j].dist_c2n,
                                            movieInfo.nodes[i].neighbors[j].dist_n2c);
//            if (movieInfo.nodes[i].neighbors[j].dist_c2n == 0 ||
//                    movieInfo.nodes[i].neighbors[j].dist_n2c == 0){
//                qDebug("cost not consistent");
//            }
        }
    }
}
/**
 * @brief calCellFootprintsDistance: get the cell-cell voxelwise distance
 * @return
 */
void cellTrackingMain::calCellFootprintsDistance(vector<float> &nn_dist){
    if (movieInfo.tracks.empty()){ // if there is no track info, use all node with nearest neighbor
        //reCalculateCellsDistances(nn_dist);
        float cur_dist;
        nn_dist.resize(movieInfo.nodes.size());
        size_t nn_dist_cnt = 0;
        FOREACH_i(movieInfo.nodes){
            nn_dist[nn_dist_cnt] = INFINITY;
            FOREACH_j(movieInfo.nodes[i].neighbors){
                // cal the distance between two nodes
                cur_dist = MAX(movieInfo.nodes[i].neighbors[j].dist_c2n,
                               movieInfo.nodes[i].neighbors[j].dist_n2c);
                if (cur_dist < nn_dist[nn_dist_cnt]){
                    nn_dist[nn_dist_cnt] = cur_dist;
                }
            }
            if(nn_dist[nn_dist_cnt] == 0){
                qDebug("we found two cell exactly overlapped, this should not happen");
            }
            if(!isinf(nn_dist[nn_dist_cnt])){
                nn_dist_cnt ++;
            }
        }
        nn_dist.resize(nn_dist_cnt);
    }else{ // if there are tracks, use long tracks
        float min_dist;
        nn_dist.resize(movieInfo.nodes.size());
        size_t nn_dist_cnt = 0;
        FOREACH_i(movieInfo.tracks){
            if(movieInfo.tracks[i].size() < p4tracking.validtrackLength4var){
                continue;
            }
            for(size_t j=0; j<movieInfo.tracks[i].size()-1; j++){
                for (nodeRelation nr : movieInfo.nodes[movieInfo.tracks[i][j]].neighbors){
                    if(nr.node_id == movieInfo.tracks[i][j+1]){
                        nn_dist[nn_dist_cnt] = MAX(nr.dist_c2n, nr.dist_n2c);
                        nn_dist_cnt ++;
                        break;
                    }
                }
            }
        }
        nn_dist.resize(nn_dist_cnt);
    }

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
        movieInfo.nodes[i].neighbors.resize(nei_ids.size());
        movieInfo.nodes[i].preNeighbors.reserve(0); //initialized preNeighbor
        if (nei_ids.size() < 1){
            continue;
        }
        movieInfo.overall_neighbor_num += nei_ids.size();
        FOREACH_j(nei_ids){
            movieInfo.nodes[i].neighbors[j].node_id = nei_ids[j];
            movieInfo.nodes[i].neighbors[j].overlap_size = cellOverlapSize(i, nei_ids[j], cellSegment);
        }
    }
    // shift initialize as all 0s
    movieInfo.frame_shift_xyz.resize(cellSegment.number_cells.size());
    FOREACH_i(movieInfo.frame_shift_xyz){
        movieInfo.frame_shift_xyz[i] = {0, 0, 0};
    }
    // calculate the distances between each cell and its neighbors
    reCalculateCellsDistances();
    updateGammaParam();
    bool update_preNei_cost = false; // preNei has not been initilized
    updateArcCost(update_preNei_cost); // calculate the link cost from given new gamma parameter
    // initialize all node as not stable
    FOREACH_i(movieInfo.nodes){
        movieInfo.nodes[i].stable_status = NOT_STABLE;
    }
    // start a naive linking
    bool tracking4jumpCost_and_driftCorrection_Only = true;
    mccTracker_one2one(tracking4jumpCost_and_driftCorrection_Only); // jumpCost are updated, so the distance will

    /**
     * These 4 functions will not be called again, which means the drift and cell-cell distance and Gamma
     * parameters will be fixed in the following test.
     */
    driftCorrection(); // update frame drifting
    reCalculateCellsDistances(); //this distance will not wholy update any more
    updateGammaParam();
    initPreNeighborInfo(); // preNei initilized
    // update arc cost and stable status with the parameters (these generally are conducted in mccTracker_one2one).
    updateArcCost();
    // we need to update parent and kids given tracks first
    node2trackUpt(true);
    track2parentKid();
    stableSegmentFixed();
}
/**
 * @brief extractNeighborIds: for a given node idx, get its possible neighbors
 * @param cell_label_maps
 * @param node_idx
 * @param nei_idxs
 */
void cellTrackingMain::extractNeighborIds(vector<Mat> &cell_label_maps, size_t cell_idx, vector<size_t> &nei_idxs){
    if(movieInfo.frames[cell_idx] == (int)cell_label_maps.size()) return;
    for (long i = movieInfo.frames[cell_idx] + 1; i <= MIN(cell_label_maps.size()-1, movieInfo.frames[cell_idx]+p4tracking.k); i++){
        unordered_set<int> nei_labels;
        FOREACH_j(movieInfo.voxIdx[cell_idx]){
            if (cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j])>0){
                nei_labels.insert(cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j]));
            }
        }
        for(auto it : nei_labels){
            nei_idxs.push_back(it - 1 + cumulative_cell_nums[i - 1]);
            vector<size_t> tmp = intersection(movieInfo.voxIdx[it - 1 + cumulative_cell_nums[i - 1]],
                    movieInfo.voxIdx[cell_idx]);
            //qDebug("tmp size %ld", tmp.size());
        }
    }
}

/**
 * @brief extractPreNeighborIds: for a given node idx, get its possible pre-neighbors
 * @param cell_label_maps
 * @param node_idx
 * @param nei_idxs
 */
void cellTrackingMain::extractPreNeighborIds(vector<Mat> &cell_label_maps, size_t cell_idx, vector<size_t> &nei_idxs){
    if(movieInfo.frames[cell_idx] == (int)0) return;
    for (long i = movieInfo.frames[cell_idx] - 1; i >= MAX(0, movieInfo.frames[cell_idx]-p4tracking.k); i--){
        unordered_set<int> nei_labels;
        FOREACH_j(movieInfo.voxIdx[cell_idx]){
            if (cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j])>0){
                nei_labels.insert(cell_label_maps[i].at<int>(movieInfo.voxIdx[cell_idx][j]));
            }
        }
        for(auto it : nei_labels){
            nei_idxs.push_back(it - 1 + cumulative_cell_nums[i - 1]);
        }
    }
}
/**
 * @brief mccTracker_one2one: build the graph without split and merge settings, but allows jump
 */
void cellTrackingMain::mccTracker_one2one(bool get_jumpCost_only, bool get_res_from_txt){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      We directly use the interface I wrote for python with the following function name               //
    // long long* pyCS2(long *msz, double *mtail, double *mhead, double *mlow, double *macap, double *mcost)//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long n = 2*movieInfo.nodes.size() + 1; // # vertices in graph
    long long m = movieInfo.overall_neighbor_num + 3*movieInfo.nodes.size(); // # arcs in graph (upper bound)
    double *mtail = new double[m];
    double *mhead = new double[m];
    double *mlow = new double[m]; memset(mlow, 0, sizeof(double) * m);
    double *macap = new double[m]; // cannot use memset since sizeof(double) = 8bytes
    for(long long i = 0; i<m; i++) macap[i] = 1;
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 0;
    int rand_max = 100;
    int node_valid_cnt = 0;
    float linkCostUpBound = INFINITY;
    for(size_t i=0; i < movieInfo.nodes.size(); i++){
        nodeInfo &node = movieInfo.nodes[i];
        if(!get_res_from_txt && movieInfo.voxIdx[i].size() == 0) continue;

        node_valid_cnt ++;
        // in arc
        mtail[arc_cnt] = src_id;
        mhead[arc_cnt] =  2 * node.node_id+1;
        mcost[arc_cnt] = round((double)node.in_cost * 1e7); //  + rand() % rand_max add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // observation arc
        mtail[arc_cnt] = 2 * node.node_id+1;
        mhead[arc_cnt] =  2 * node.node_id + 2;
        //mcost[arc_cnt] = round((double)p4tracking.observationCost * 1e7); // add some randomness to make sure unique optimal solution
        mcost[arc_cnt] = round((-((double)node.in_cost+(double)node.out_cost) + 0.00001)  * 1e7);
        arc_cnt ++;
        // out arc
        mhead[arc_cnt] = src_id;
        mtail[arc_cnt] =  2 * node.node_id + 2;
        mcost[arc_cnt] = round((double)node.out_cost * 1e7); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        //linkCostUpBound = node.out_cost;
        // link with neighbors
        for(nodeRelation neighbor : node.neighbors){
            if(neighbor.link_cost < linkCostUpBound){
                mtail[arc_cnt] = 2 * node.node_id + 2;
                mhead[arc_cnt] =  2 * neighbor.node_id + 1;
                mcost[arc_cnt] = round((double)neighbor.link_cost * 1e7); // add some randomness to make sure unique optimal solution
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
    vector<double> cost;
    double all_cost = 0;
    vector<size_t> curr_track;
    movieInfo.tracks.clear();
    for(long long i = 1; i < track_vec[0] + 1; i++){
        if(track_vec[i] > 0){
            curr_track.push_back(track_vec[i]);
        }
        else{
            cost.push_back(((double)track_vec[i])/1e7);
            all_cost += cost[cost.size()-1];
            vector<size_t> new_track;
            for(size_t j = 0; j < curr_track.size(); j+=2){
                if(curr_track[j] < 1 || (!get_res_from_txt &&
                        movieInfo.voxIdx[size_t((curr_track[j]-1) / 2)].size() == 0)){
                    qFatal("Node id wrong!");
                }
                new_track.push_back(size_t((curr_track[j]-1) / 2));
            }
            if (new_track.size() > 1){
                movieInfo.tracks.push_back(new_track);
            }
            curr_track.clear();
        }
    }
    if(!get_res_from_txt){
        // update the jumpCost in p4tracking if no split/merge allowed
        updateJumpCost();
        if(!get_jumpCost_only){ // means this time of tracking is just to estimate jump cost and/or drift
            node2trackUpt(true);
            // update parent and kids given tracks
            track2parentKid();
            updateArcCost();
            // by the way, also update stable segmentations
            stableSegmentFixed();
        }
    }
    delete[] msz;
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
    double *mlow = new double[m]; memset(mlow, 0, sizeof(double) * m);
    double *macap = new double[m]; // cannot use memset since sizeof(double) = 8bytes
    for(long long i = 0; i<m; i++) macap[i] = 1;
    double *mcost = new double[m];

    long long arc_cnt = 0;
    double src_id = 0;
    int rand_max = 100;
    bool link_adj_frs = false;
    if(p4tracking.splitMergeHandle == NOJUMPALL){
        link_adj_frs = true;
    }
    bool *visited = new bool[movieInfo.nodes.size()]; // ok to use memset
    memset(visited, false, sizeof(bool) * movieInfo.nodes.size());

    float linkCostUpBound = INFINITY;
    // add the arcs of splitting and merging
    for(splitMergeNodeInfo &sp_mg_info : split_merge_node_info){
        if(sp_mg_info.node_id >= 0 && sp_mg_info.invalid == false){
            if(sp_mg_info.parent_flag){
                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.family_nodes[0] + 1;
                mcost[arc_cnt] = round((double)sp_mg_info.link_costs[0] * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.family_nodes[1] + 1;
                mcost[arc_cnt] = round((double)sp_mg_info.link_costs[1] * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = src_id;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 2;
                mcost[arc_cnt] = round((double)sp_mg_info.src_link_cost * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                visited[sp_mg_info.node_id] = true;
            }else{
                mtail[arc_cnt] = 2 * sp_mg_info.family_nodes[0] + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 1;
                mcost[arc_cnt] = round((double)sp_mg_info.link_costs[0] * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.family_nodes[1] + 2;
                mhead[arc_cnt] =  2 * sp_mg_info.node_id + 1;
                mcost[arc_cnt] = round((double)sp_mg_info.link_costs[1] * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                mtail[arc_cnt] = 2 * sp_mg_info.node_id + 1;
                mhead[arc_cnt] =  src_id;
                mcost[arc_cnt] = round((double)sp_mg_info.src_link_cost * 1e7); // add some randomness to make sure unique optimal solution
                arc_cnt ++;

                visited[sp_mg_info.family_nodes[0]] = true;
                visited[sp_mg_info.family_nodes[1]] = true;
            }
        }
    }
    for(nodeInfo node : movieInfo.nodes){
        if(movieInfo.voxIdx[node.node_id].size() == 0) continue;
        // in arc
        mtail[arc_cnt] = src_id;
        mhead[arc_cnt] =  2 * node.node_id + 1;
        mcost[arc_cnt] = round((double)node.in_cost * 1e7); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // observation arc
        mtail[arc_cnt] = 2 * node.node_id+1;
        mhead[arc_cnt] =  2 * node.node_id + 2;
        mcost[arc_cnt] = round((double)p4tracking.observationCost * 1e7); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        // out arc
        mhead[arc_cnt] = src_id;
        mtail[arc_cnt] =  2 * node.node_id+2;
        mcost[arc_cnt] = round((double)node.out_cost * 1e7); // add some randomness to make sure unique optimal solution
        arc_cnt ++;
        if(visited[node.node_id]){ // for those have split/merge arcs, not more neighbor-linking needed
            continue;
        }
        // link with neighbors:!!! here needs remove arcs added in previous section
        if(link_adj_frs){
            int curr_frame = movieInfo.frames[node.node_id];
            for(nodeRelation neighbor : node.neighbors){
                if(movieInfo.frames[neighbor.node_id]-1 == curr_frame){
                    if(neighbor.link_cost < linkCostUpBound){
                        mtail[arc_cnt] = 2 * node.node_id+2;
                        mhead[arc_cnt] =  2 * neighbor.node_id+1;
                        mcost[arc_cnt] = round((double)neighbor.link_cost * 1e7); // add some randomness to make sure unique optimal solution
                        arc_cnt ++;
                    }
                }
            }
        }else{
            for(nodeRelation neighbor : node.neighbors){
                if(neighbor.link_cost < linkCostUpBound){
                    mtail[arc_cnt] = 2 * node.node_id+1;
                    mhead[arc_cnt] =  2 * neighbor.node_id+2;
                    mcost[arc_cnt] = round((double)neighbor.link_cost * 1e7); // add some randomness to make sure unique optimal solution
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
    vector<double> cost;
    double all_cost = 0;
    vector<size_t> curr_track;
    movieInfo.tracks.clear();
    for(long long i = 1; i < track_vec[0] + 1; i++){
        if(track_vec[i] > 0){
            curr_track.push_back(track_vec[i]);
        }
        else{
            cost.push_back(((double)track_vec[i])/1e7);
            all_cost += cost[cost.size()-1];
            vector<size_t> new_track;
            for(size_t j = 0; j < curr_track.size(); j++){
                if (new_track.size() == 0 ||
                        size_t((curr_track[j]-1) / 2) != new_track[new_track.size()-1]){
                    new_track.push_back(size_t((curr_track[j]-1) / 2));
                    if(movieInfo.voxIdx[size_t((curr_track[j]-1) / 2)].size() == 0){
                        qFatal("node id assignment is wrong!");
                    }
                }
            }
            if (new_track.size() > 1){
                movieInfo.tracks.push_back(new_track);
            }
            curr_track.clear();
        }
    }
    // update parent and kids given tracks
    track2parentKid();
    if(m_append > 0){// merge tracks with overlapped nodes
        mergeOvTracks2();
    }
    delete[] msz;
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
 * @brief mergeOvTracks2: merge overlapped tracks and update node id to track information
 */
void cellTrackingMain::mergeOvTracks2(){
    vector<std::vector<size_t>> new_tracks;
    long valid_num = 0;
    vector<bool> visited(movieInfo.nodes.size());
    fill(visited.begin(), visited.end(), false);
    FOREACH_i(movieInfo.nodes){
        if(visited[i] || (movieInfo.nodes[i].parent_num == 0
                          && movieInfo.nodes[i].kid_num==0)){
            continue;
        }
        vector<size_t> one_track;
        deque<size_t> n_ids;
        n_ids.push_back(i);
        while(!n_ids.empty()){
            size_t cur_id = n_ids.front();
            one_track.push_back(cur_id);
            for (int kk = 0; kk < movieInfo.nodes[cur_id].kid_num; kk++){
                if (!visited[movieInfo.nodes[cur_id].kids[kk]]){
                    n_ids.push_back(movieInfo.nodes[cur_id].kids[kk]);
                }
            }
            for (int kk = 0; kk < movieInfo.nodes[cur_id].parent_num; kk++){
                if (!visited[movieInfo.nodes[cur_id].parents[kk]]){
                    n_ids.push_back(movieInfo.nodes[cur_id].parents[kk]);
                }
            }
            visited[cur_id] = true;
            n_ids.pop_front();
        }
        vec_unique(one_track);
        new_tracks.push_back(one_track);
        valid_num ++;
    }
    new_tracks.resize(valid_num);
    movieInfo.tracks.clear();
    movieInfo.tracks = new_tracks;
    node2trackUpt(false);
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
    node2trackUpt(false);
}
void cellTrackingMain::node2trackUpt(bool one2one_track){ // default is false
    for(nodeInfo &nf : movieInfo.nodes){
        nf.nodeId2One2OneTrackId = -1;
        nf.nodeId2trackId = -1;
        nf.nodeLocInTrack = -1;
    }
    FOREACH_i(movieInfo.tracks){
        FOREACH_j(movieInfo.tracks[i]){
            movieInfo.nodes[movieInfo.tracks[i][j]].nodeId2trackId = i;
            movieInfo.nodes[movieInfo.tracks[i][j]].nodeLocInTrack = j;
        }
    }
    if(one2one_track){
        for(nodeInfo &nf : movieInfo.nodes){
            nf.nodeId2One2OneTrackId = nf.nodeId2trackId;
        }
    }
}
void cellTrackingMain::updateJumpCost(){
    size_t overall_arcs = 0;
    vector<size_t> jumpArcs(p4tracking.k);
    fill(jumpArcs.begin(), jumpArcs.end(), 0);
    for(vector<size_t> &track : movieInfo.tracks){
        if(track.size() < p4tracking.validtrackLength4var) continue;
        overall_arcs += track.size() - 1;
        for (size_t i = 1; i < track.size(); i++){
            jumpArcs[movieInfo.frames[track[i]] - movieInfo.frames[track[i-1]] - 1] ++;
        }
    }
    FOREACH_i(p4tracking.jumpCost){
        if(overall_arcs == 0){//there is no meaningful arcs
            p4tracking.jumpCost[i] = 1.0;
        }else{
            p4tracking.jumpCost[i] = (float)((double)jumpArcs[i] / overall_arcs);
        }
    }
}
/**
 * @brief track2parentKid: !! this function should be called before we do track merge
 */
void cellTrackingMain::track2parentKid(){
    for(nodeInfo &n : movieInfo.nodes){
        n.parent_num = 0;
        n.kid_num = 0;
    }
    for(size_t t_id=0; t_id<movieInfo.tracks.size(); t_id++){
        vector<size_t> &track = movieInfo.tracks[t_id];
        for (size_t i = 1; i < track.size(); i++){
//            if(track[i] == 259){// && track[i] == 190
//                qDebug("Possible mem leaking");
//            }
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
                if (track[i-1] == movieInfo.nodes[track[i]].preNeighbors[k].node_id){
                    movieInfo.nodes[track[i]].parent_cost[movieInfo.nodes[track[i]].parent_num] =
                            movieInfo.nodes[track[i]].preNeighbors[k].link_cost;
                    break;
                }
            }
            movieInfo.nodes[track[i]].parent_num ++;

            if(movieInfo.nodes[4].kid_num > 2){
                qFatal("Possible mem leaking");
            }
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
    for(int i=1;i < movieInfo.frame_shift_xyz.size(); i++){
        movieInfo.frame_shift_xyz[i][0] /= MAX(1, samples4driftCumulate[i]);
        movieInfo.frame_shift_xyz[i][0] += movieInfo.frame_shift_xyz[i-1][0];
        movieInfo.frame_shift_xyz[i][1] /= MAX(1, samples4driftCumulate[i]);
        movieInfo.frame_shift_xyz[i][1] += movieInfo.frame_shift_xyz[i-1][1];
        movieInfo.frame_shift_xyz[i][2] /= MAX(1, samples4driftCumulate[i]);
        movieInfo.frame_shift_xyz[i][2] += movieInfo.frame_shift_xyz[i-1][2];
    }
}
void cellTrackingMain::updateGammaParam(){
    vector<float> nn_dist;
    calCellFootprintsDistance(nn_dist); // the nearest neighbor
    if (nn_dist.size() > 200){
        truncatedGammafit(nn_dist, movieInfo.ovGammaParam[0], movieInfo.ovGammaParam[1]);
    }
}
void cellTrackingMain::updateArcCost(bool updatePreNei){
    for(nodeInfo &n : movieInfo.nodes){
        for(nodeRelation &neighbor : n.neighbors){
            neighbor.link_cost = distance2cost(MAX(neighbor.dist_c2n, neighbor.dist_n2c),
                                               p4tracking.jumpCost[movieInfo.frames[neighbor.node_id] - movieInfo.frames[n.node_id]-1]);
        }
    }
    if(updatePreNei){
        bool update_cost_only = true;
        updatePreNeighborInfo(update_cost_only);
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
    int min_valid_node_cluster_sz = p4tracking.min_stable_node_cluster_sz;
    FOREACH_i(movieInfo.tracks){
        if(movieInfo.tracks[i].size() < p4tracking.validtrackLength4var) continue;
        vector<float> arc_costs;
        getArcCostOne2OneTrack(i, arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][0] = vec_mean(arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][1] = vec_median(arc_costs);
        movieInfo.track_arcs_avg_mid_std[i][2] = vec_stddev(arc_costs);

        size_t start_id = 0, end_id = 0;
        for(size_t j = 0; j<arc_costs.size(); j++){
            if (arc_costs[j] <= max_cost &&
                    movieInfo.frames[movieInfo.tracks[i][j+1]] - movieInfo.frames[movieInfo.tracks[i][j]] == 1){
                end_id = j + 1;
            }else{
                if ((end_id - start_id) >= min_valid_node_cluster_sz){
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
/** debug done
 * @brief bestPeerCandidate: given a region id, find if there are two regions that can together fit into its territory
 * @param node_id
 * @param bestPeer
 * @param parent_flag
 */
float cellTrackingMain::bestPeerCandidate(size_t node_id, vector<size_t> &out_bestPeer, bool parent_flag){
    vector<size_t> peerCandidates(0);
    out_bestPeer.resize(0);
    if(parent_flag){ // test kids
        if (movieInfo.nodes[node_id].neighbors.size() >= 2){
            int frame2test = movieInfo.frames[node_id] + 1;
            float curr_cost;
            bool currBest;
            for(nodeRelation neighbor : movieInfo.nodes[node_id].neighbors){
                if (movieInfo.frames[neighbor.node_id] != frame2test) continue;
                curr_cost = neighbor.dist_n2c;
                currBest = true;
                for(nodeRelation preNeighbor : movieInfo.nodes[neighbor.node_id].preNeighbors){
                    if((curr_cost > preNeighbor.dist_c2n || preNeighbor.link_cost < abs(p4tracking.observationCost))
                            && preNeighbor.node_id != node_id &&
                            movieInfo.frames[node_id] == movieInfo.frames[preNeighbor.node_id]){
                        currBest = false;
                        break;
                    }
                }
                if(currBest) {
                    peerCandidates.push_back(neighbor.node_id);
                }
            }
        }
    }else{
        if(movieInfo.nodes[node_id].preNeighbors.size() >= 2){
            int frame2test = movieInfo.frames[node_id] - 1;
            float curr_cost;
            bool currBest;
            for(nodeRelation preNeighbor : movieInfo.nodes[node_id].preNeighbors){
                if (movieInfo.frames[preNeighbor.node_id] != frame2test) continue;
                curr_cost = preNeighbor.dist_n2c;
                currBest = true;
                for(nodeRelation neighbor : movieInfo.nodes[preNeighbor.node_id].neighbors){
                    if((curr_cost > neighbor.dist_c2n  || neighbor.link_cost < abs(p4tracking.observationCost))
                            && neighbor.node_id != node_id &&
                            movieInfo.frames[node_id] == movieInfo.frames[neighbor.node_id]){
                        currBest = false;
                        break;
                    }
                }
                if(currBest) {
                    peerCandidates.push_back(preNeighbor.node_id);
                }
            }
        }
    }
    if(peerCandidates.size() < 2){
        return INFINITY;
    }else if(peerCandidates.size() == 2){
        out_bestPeer = peerCandidates;
        float dummy_c2n, dummy_n2c;
        return distance2cost(voxelwise_avg_distance(node_id, out_bestPeer, dummy_c2n, dummy_n2c),
                             p4tracking.jumpCost[0]); // no jump will be allowed
    }else{ // more than two regions available
        float min_distance = INFINITY, curr_distance;
        float dummy_c2n, dummy_n2c;
        vector<size_t> curr_peer(2);
        FOREACH_i(peerCandidates){
            curr_peer[0] = peerCandidates[i];
            for(size_t j = i + 1; j < peerCandidates.size(); j ++){
                curr_peer[1] = peerCandidates[j];
                curr_distance = voxelwise_avg_distance(node_id, curr_peer, dummy_c2n, dummy_n2c);
                if(curr_distance < min_distance){
                    out_bestPeer = curr_peer;
                    min_distance = curr_distance;
                }
            }
        }
        return distance2cost(min_distance, p4tracking.jumpCost[0]);
    }
}

/**
 * @brief peerRegionVerify: to verify if there are two regions in adjacent frame that jointly has better linking cost
 * @param node_id
 * @param cost_good2go
 * @param parents_test
 * @param split_merge_node_info
 */
void cellTrackingMain::peerRegionVerify(size_t node_id, float cost_good2go, bool kids_test,
                                        vector<splitMergeNodeInfo> &split_merge_node_info,
                                        unordered_map<long, long> &node_id2split_merge_node_id){
    float merged_cost, min_exist_cost = INFINITY;
    vector<size_t> bestPeers(2);
    merged_cost = bestPeerCandidate(node_id, bestPeers, kids_test);
    if(isinf(merged_cost)) return;
    int cur_frame = movieInfo.frames[node_id];
    if(kids_test){// testing current node's kids
        for(nodeRelation &neighbor : movieInfo.nodes[node_id].neighbors){
            if(neighbor.link_cost < min_exist_cost &&
                    movieInfo.frames[neighbor.node_id] == cur_frame + 1){
                min_exist_cost = neighbor.link_cost;
            }
        }
    }else{
        for(nodeRelation &preNeighbor : movieInfo.nodes[node_id].preNeighbors){
            if(preNeighbor.link_cost < min_exist_cost &&
                    movieInfo.frames[preNeighbor.node_id] == cur_frame - 1){
                min_exist_cost = preNeighbor.link_cost;
            }
        }
    }
    splitMergeNodeInfo split_merge_peer;
    if(merged_cost < (min_exist_cost-0.001) && min_exist_cost > cost_good2go){
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
            split_merge_peer.parent_flag = kids_test;
            split_merge_peer.src_link_cost = 0.0;
            split_merge_peer.invalid = false;
            split_merge_node_info.push_back(split_merge_peer);
            // 'long' is enough for our data analysis, indeed all the size_t are redundantly large
            if(kids_test){
                node_id2split_merge_node_id.insert(make_pair<long, long>(node_id, split_merge_node_info.size()-1));
            }else{
                node_id2split_merge_node_id.insert(make_pair<long, long>(-node_id, split_merge_node_info.size()-1));
            }
        }
    }
}
/** debug done
 * @brief detectPeerRegions: test if two regions should be merged based on their common kid or parent region
 * @param merge_node_idx
 * @param split_node_idx
 */
void cellTrackingMain::detectPeerRegions(vector<splitMergeNodeInfo> &split_merge_node_info,
                                         unordered_map<long, long> &node_id2split_merge_node_id){
    float cost_good2go; // the average cost of arcs in a track: only arc with larger cost will be test
    vector<size_t> p_cnt, k_cnt, all_cnt;
    for(nodeInfo node : movieInfo.nodes){
        if(node.nodeId2One2OneTrackId >= 0 &&
                node.nodeId2One2OneTrackId >= movieInfo.track_arcs_avg_mid_std.size()){
            qFatal("Error on node2trackid");
        }
        if(node.nodeId2One2OneTrackId < 0){
            cost_good2go = 0; // for isolated node, does not consider
        }
        else {
            cost_good2go = MIN(movieInfo.track_arcs_avg_mid_std[node.nodeId2One2OneTrackId][0],
                            movieInfo.track_arcs_avg_mid_std[node.nodeId2One2OneTrackId][1]); // use median to replace mean
        }

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
//        if (node.node_id == 8275){
//            qDebug("check point");
//        }
//        size_t tmp1 = split_merge_node_info.size();
        if(parents_test){
            peerRegionVerify(node.node_id, cost_good2go, false, split_merge_node_info, node_id2split_merge_node_id);

        }
//        if (split_merge_node_info.size() > tmp1){
//            p_cnt.push_back(node.node_id);
//        }
//        size_t tmp2 = split_merge_node_info.size();
        if(kids_test){
            peerRegionVerify(node.node_id, cost_good2go, true, split_merge_node_info, node_id2split_merge_node_id);
        }
//        if (split_merge_node_info.size() > tmp2){
//            k_cnt.push_back(node.node_id);
//        }
//        if (split_merge_node_info.size() > MIN(tmp1, tmp2)){
//            all_cnt.push_back(node.node_id);
//        }
    }
    //qDebug("check point");
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
    unordered_map<long, long> node_id2split_merge_node_id; // save node_id and its index in split_merge_node_info
    detectPeerRegions(split_merge_node_info, node_id2split_merge_node_id);
    vector<size_t> invalid_cnt;
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
                if(i != best_idx){ // remove contradict ones with smaller size (only keep best_idx)
                    split_merge_node_info[contradict_groups[i].first].invalid = true; //
                    invalid_cnt.push_back(split_merge_node_info[contradict_groups[i].first].node_id);
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
    if(valid_p && valid_k){// both parents/kids can be used to split root, but these 4 regions does not fit
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
    float merged_cost = distance2cost(merged_distance, 1.0); // adjacent frame, no need to add punishment

    if (merged_cost >= abs(p4tracking.observationCost)){
        return CONSISTENCY_NOT_SURE;
    }
    float map_2x2[4] = {INFINITY, INFINITY, INFINITY, INFINITY};
    int cnt = 0;
    for(size_t p : movieInfo.nodes[node_id].parents){
        for(size_t k : movieInfo.nodes[node_id].kids){
            for(nodeRelation nr : movieInfo.nodes[p].neighbors){
                if (nr.node_id == k){
                    map_2x2[cnt] = distance2cost(nr.dist_c2n, 1.0);
                    break;
                }
            }
            cnt++;
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
    delete frame_root;
    delete frame_seed;
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

    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*root_frame; // sub-matrix pointer
    Mat *frame_root = new Mat(3, cellSegment.normalized_data4d.size, CV_8U, ind);

    unsigned char *ind2 = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*seed_frame; // sub-matrix pointer
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
            delete frame_root;
            delete frame_seed;
            return false;
        }
    }
    delete frame_root;
    delete frame_seed;
    return true;
}
/**
 * @brief seedRefine_gap: put seeds' location idx into a binary seed map
 * @param cellSegment
 * @param root_idx
 * @param root_frame
 * @param seeds_idx
 * @param seed_frame
 * @param ref_seeds_idx
 * @return
 */
bool cellTrackingMain::seedsRefine_gap(Mat1b &possibleGaps, vector<vector<size_t>> &seeds_idx, Mat1i &outLabelMap){
    assert(seeds_idx.size() > 0);
    outLabelMap = Mat(possibleGaps.dims, possibleGaps.size, CV_32S, Scalar(0));
    Mat1b noGaps;
    bitwise_not(possibleGaps, noGaps);
    for(int i = 0; i< seeds_idx.size(); i++){
        Mat1b tmp = Mat(possibleGaps.dims, possibleGaps.size, CV_8U, Scalar(0));
        setValMat(tmp, CV_8U, seeds_idx[i], 255);
        bitwise_and(noGaps, tmp, tmp);
        Mat1i tmp_label;
        int ncc = connectedComponents3d(&tmp, tmp_label, 26);
        if(ncc > 1){
            int lid = largestRegionIdExtract(&tmp_label, ncc);
            tmp = tmp_label == lid;
        }
        if(ncc==0 || isempty(tmp, CV_8U)){
            return false;
        }else{
            setValMat(outLabelMap, CV_32S, &tmp, (float)i+1);
        }
    }
    return true;
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
    //ccShowSliceLabelMat(label_map);
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
        if(vals.size() > 0){
            int it = Mode(vals.begin(), vals.end());
            if(it != 0){
                setValMat(seeds_map, CV_32S, voxIdxList[i], (float)it);
                seeds_idx_used[it-1] = true;
            }
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

    int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    vector<size_t> reg2split_in_seed;
    coordinateTransfer(reg2split, single_frame->size, reg2split_in_seed, crop_start, seed.idMap.size);

    Mat mask = seed.eigMap2d < 0;
    setValMat(seed.eigMap2d, CV_32F, &mask, 0.0);
    Mat1b possibleGap2d = seed.eigMap2d > 0;
    float max_val = getMaxValMat(seed.eigMap2d, CV_32F, reg2split_in_seed);
    scale_vol(&seed.eigMap2d, CV_32F, &seed.score2d, 0.001, 1, 0, max_val);

    mask = seed.eigMap3d < 0;
    setValMat(seed.eigMap3d, CV_32F, &mask, 0.0);
    Mat1b possibleGap3d = seed.eigMap3d > 0;
    max_val = getMaxValMat(seed.eigMap3d, CV_32F, reg2split_in_seed);
    scale_vol(&seed.eigMap3d, CV_32F, &seed.score3d, 0.001, 1, 0, max_val);
    if(usePriorGapMap){
        Mat subGapMap;
        subVolExtract(&validGapMaps[reg2split_frame], CV_8U, subGapMap, seed.crop_range_yxz);// deep copy
        bitwise_and(possibleGap2d, subGapMap, possibleGap2d);
        bitwise_and(possibleGap3d, subGapMap, possibleGap3d);
    }

    // start to build seed map
    Mat1b fgMap;
    if (reg2split_idx < movieInfo.labelInMap.size()){
        fgMap = seed.idMap == movieInfo.labelInMap[reg2split_idx];
    }else{ // input is not a valid region exist in the label map
        fgMap = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_8U);
//        int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
//        vector<size_t> tmp_idx;
//        coordinateTransfer(reg2split, single_frame->size, tmp_idx, crop_start, seed.idMap.size);
        setValMat(fgMap, CV_8U, reg2split_in_seed, 255);
    }

    vector<vector<size_t>> seeds_idx(reg4seeds.size());
    //int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    FOREACH_i(reg4seeds){
        coordinateTransfer(reg4seeds[i], single_frame->size, seeds_idx[i], crop_start, seed.idMap.size);
    }
    Mat1i seedsMap;
    bool success = binary_seedsMap_create(fgMap, &possibleGap3d, nullptr, seeds_idx, seedsMap);
    if(!success){
        success = binary_seedsMap_create(fgMap, nullptr, &possibleGap2d, seeds_idx, seedsMap);
        if(!success){
            delete single_frame;
            return false;
        }
    }
    seed.scoreMap = seed.score2d + seed.score3d;
    Mat1i outLabelMap;
    regionGrow(&seedsMap, 2, outLabelMap, &seed.scoreMap, &fgMap,
               cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, false);

    //ccShowSliceLabelMat(outLabelMap);
    splitRegs.resize(2);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[0], 1);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[1], 2);
    int crop_start_yxz[3] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start, -seed.crop_range_yxz[2].start};
    coordinateTransfer(splitRegs[0], outLabelMap.size, splitRegs[0], crop_start_yxz,single_frame->size);
    coordinateTransfer(splitRegs[1], outLabelMap.size, splitRegs[1], crop_start_yxz,single_frame->size);
    delete single_frame;

//    vector<size_t> tmp1 = intersection(reg4seeds[0], splitRegs[0]);
//    vector<size_t> tmp2 = intersection(reg4seeds[1], splitRegs[1]);
    return true;
}
/** debug done
 * @brief bisectRegion_bruteforce
 * @param cellSegment
 * @param reg2split_idx
 * @param reg2split
 * @param reg2split_frame
 * @param reg4seeds
 * @param reg4seeds_frame
 * @param usePriorGapMap
 * @param splitRegs
 * @return
 */
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

    int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    vector<size_t> reg2split_in_seed;
    coordinateTransfer(reg2split, single_frame->size, reg2split_in_seed, crop_start, seed.idMap.size);

    Mat mask = seed.eigMap2d < 0;
    setValMat(seed.eigMap2d, CV_32F, &mask, 0.0);
    Mat1b possibleGap2d = seed.eigMap2d > 0;
    float max_val = getMaxValMat(seed.eigMap2d, CV_32F, reg2split_in_seed);
    scale_vol(&seed.eigMap2d, CV_32F, &seed.score2d, 0.001, 1, 0, max_val);

    mask = seed.eigMap3d < 0;
    setValMat(seed.eigMap3d, CV_32F, &mask, 0.0);
    Mat1b possibleGap3d = seed.eigMap3d > 0;
    max_val = getMaxValMat(seed.eigMap3d, CV_32F, reg2split_in_seed);
    scale_vol(&seed.eigMap3d, CV_32F, &seed.score3d, 0.001, 1, 0, max_val);
    if(usePriorGapMap){
        Mat subGapMap;
        subVolExtract(&validGapMaps[reg2split_frame], CV_8U, subGapMap, seed.crop_range_yxz);// deep copy
        bitwise_and(possibleGap2d, subGapMap, possibleGap2d);
        bitwise_and(possibleGap3d, subGapMap, possibleGap3d);
    }

    // start to build seed map
    Mat1b fgMap;
    if (reg2split_idx < movieInfo.labelInMap.size()){
        fgMap = seed.idMap == movieInfo.labelInMap[reg2split_idx];
    }else{ // input is not a valid region exist in the label map
        fgMap = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_8U);
        //int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
        //vector<size_t> tmp_idx;
        //coordinateTransfer(reg2split, single_frame->size, tmp_idx, crop_start, seed.idMap.size);
        setValMat(fgMap, CV_8U, reg2split_in_seed, 255);
    }

    /** ***** THIS the only difference between brute-force and gap-based method ******** **/
    vector<vector<size_t>> ref_seeds_idx(reg4seeds.size());
    bool valid_seeds_flag = true;
    valid_seeds_flag = seedsRefine_intensity(cellSegment, reg2split, reg2split_frame,
                                             reg4seeds, reg4seeds_frame, ref_seeds_idx);
    if (!valid_seeds_flag) {
        delete single_frame;
        return false;
    }

    vector<vector<size_t>> seeds_idx_in_subVol(ref_seeds_idx.size());
    //int crop_start[3] = {seed.crop_range_yxz[0].start, seed.crop_range_yxz[1].start, seed.crop_range_yxz[2].start};
    FOREACH_i(reg4seeds){
        coordinateTransfer(ref_seeds_idx[i], single_frame->size, seeds_idx_in_subVol[i], crop_start, seed.idMap.size);
    }
    Mat1i seedsMap;
    valid_seeds_flag = seedsRefine_gap(possibleGap3d, seeds_idx_in_subVol, seedsMap);

    //ccShowSliceLabelMat(seedsMap);
    /** ************************************ Done ************************************* **/

    seed.scoreMap = seed.score2d + seed.score3d;
    Mat1i outLabelMap;
    regionGrow(&seedsMap, 2, outLabelMap, &seed.scoreMap, &fgMap,
               cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, false);

    //ccShowSliceLabelMat(outLabelMap);
    splitRegs.resize(2);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[0], 1);
    extractVoxIdxGivenId(&outLabelMap, splitRegs[1], 2);
    int crop_start_yxz[3] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start, -seed.crop_range_yxz[2].start};
    coordinateTransfer(splitRegs[0], outLabelMap.size, splitRegs[0], crop_start_yxz,single_frame->size);
    coordinateTransfer(splitRegs[1], outLabelMap.size, splitRegs[1], crop_start_yxz,single_frame->size);
    //vector<size_t> tmp1 = intersection(reg4seeds[0], splitRegs[0]);
    //vector<size_t> tmp2 = intersection(reg4seeds[1], splitRegs[1]);
    delete single_frame;
    return true;
}
/** debug done
 * @brief bisectValidTest
 * @param cellSegment
 * @param reg2split_idx
 * @param reg2split
 * @param reg2split_frame
 * @param reg4seeds
 * @param reg4seeds_frame
 * @param gapBasedSplit
 * @param usePriorGapMap
 * @param splitRegs
 * @param reg4seeds2splitRes_costs
 * @return
 */
bool cellTrackingMain::bisectValidTest(cellSegmentMain &cellSegment, size_t reg2split_idx, vector<size_t> reg2split,
                                       int reg2split_frame, vector<vector<size_t>> reg4seeds, int reg4seeds_frame,
                                       bool gapBasedSplit, bool usePriorGapMap, vector<vector<size_t>> &splitRegs,
                                       float *reg4seeds2splitRes_costs){
    splitRegs.clear();
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
        reg4seeds2splitRes_costs[0] = voxelwise_avg_distance(reg4seeds[0], reg4seeds_frame, splitRegs[0], reg2split_frame,
                cellSegment.cell_label_maps[0].size, dummy_c2n, dummy_n2c);
        reg4seeds2splitRes_costs[0] = distance2cost(reg4seeds2splitRes_costs[0], p4tracking.jumpCost[0]); // adjacent frames
        reg4seeds2splitRes_costs[1] = voxelwise_avg_distance(reg4seeds[1], reg4seeds_frame, splitRegs[1], reg2split_frame,
                cellSegment.cell_label_maps[0].size, dummy_c2n, dummy_n2c);
        reg4seeds2splitRes_costs[1] = distance2cost(reg4seeds2splitRes_costs[1], p4tracking.jumpCost[0]);
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
                if(it < 0 || it != movieInfo.nodes[it].node_id){
                    qDebug("check point: the node has been removed or changed");
                }
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

    int loopCnt = 1;
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

//            if(kid_flag && movieInfo.nodes[root_node2].parent_num > 1){
//                bases[0] = movieInfo.nodes[root_node2].parents[0];
//                bases[1] = movieInfo.nodes[root_node2].parents[1];
//                if((bases[0] != pre_bases[0] && bases[0] != pre_bases[1]) ||
//                        (bases[1] != pre_bases[0] && bases[1] != pre_bases[1])){
//                    return false;
//                }
//            }else if(!kid_flag && movieInfo.nodes[root_node2].kid_num > 1){
//                bases[0] = movieInfo.nodes[root_node2].kids[0];
//                bases[1] = movieInfo.nodes[root_node2].kids[1];
//                if((bases[0] != pre_bases[0] && bases[0] != pre_bases[1]) ||
//                        (bases[1] != pre_bases[0] && bases[1] != pre_bases[1])){
//                    return false;
//                }
//            }
            if(kid_flag){
                for(int pp = 0 ; pp < movieInfo.nodes[root_node2].parent_num; pp ++){
                    pre_bases.push_back(movieInfo.nodes[root_node2].parents[pp]);
                }
            }else if(!kid_flag && movieInfo.nodes[root_node2].kid_num > 1){
                for(int kk = 0 ; kk < movieInfo.nodes[root_node2].kid_num; kk ++){
                    pre_bases.push_back(movieInfo.nodes[root_node2].kids[kk]);
                }
            }
            vec_unique(pre_bases);
            float dummy_c2n, dummy_n2c;
            float maxDistance = voxelwise_avg_distance(root_node2, pre_bases, dummy_c2n, dummy_n2c);
            if(distance2cost(maxDistance, p4tracking.jumpCost[0]) < p4tracking.c_ex){
                return true;
            }else{
                return false;
            }
        }// else continue;
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
        else return SPLIT_BY_PARENTS;
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
    /** way 4.1 two split regions that consistent with its two kid or parent regions **/
    if(one2multiple_flag && right_cells[1].size() == 2 &&
            mergeValidTest(*right_cells[0].begin(), movieInfo.nodes[*right_cells[0].begin()].kids)){
        return MERGE_KIDS;
    }else if(!one2multiple_flag && left_cells[1].size() == 2 &&
             mergeValidTest(*left_cells[0].begin(), movieInfo.nodes[*left_cells[0].begin()].parents)){
        return MERGE_PARENTS;
    }
    if((one2multiple_flag && right_cells[1].size() != 2) ||
            (!one2multiple_flag && left_cells[1].size() != 2)){
        return MERGE_SPLIT_NOT_SURE;
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
    }
    // Way 5.3: gap exists for spliting.
    vector<vector<size_t>> reg4seeds(2), dummy_splitRegs;
    float dummy_costs[2];
    reg4seeds[0] = movieInfo.voxIdx[bases[0]];
    reg4seeds[1] = movieInfo.voxIdx[bases[1]];
    bool sectTest = bisectValidTest(cellSegment, curr_node_id, movieInfo.voxIdx[curr_node_id], movieInfo.frames[curr_node_id],
                                    reg4seeds, movieInfo.frames[bases[0]], true, true, dummy_splitRegs, dummy_costs);
    if(sectTest){
        if (one2multiple_flag) return SPLIT_BY_KIDS;
        else return SPLIT_BY_PARENTS;
    }
    // way 5.4: all other ways failed
    if(split_consistent_flag){
        if (one2multiple_flag) return MERGE_KIDS;
        else return MERGE_PARENTS;
    }
    return MERGE_SPLIT_NOT_SURE;
}
/**
 * @brief regionRefresh: based on the tracking results, check if the region should be split or merge
 * @param cellSegment
 */
void cellTrackingMain::regionRefresh(cellSegmentMain &cellSegment, vector<newFoundCellInfo> &newCells, vector<size_t> &uptCell_idxs){
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
//                if(curr_decision == SPLIT_BY_KIDS && movieInfo.nodes[curr_node_id].kid_num < 2){
//                    int tmp = regionSplitMergeJudge(cellSegment, curr_node_id, false, pvalue);
//                }
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
        unordered_set<size_t> splitted_reg_ids, merged_reg_ids;
        FOREACH_i(merged_split_peers){
            //separateRegion(cellSegmentMain &cellSegment, size_t node_idx,
            //size_t seeds[2], vector<newFoundCellInfo> &newCells)
            size_t curr_node_idx = get<0>(merged_split_peers[i]);
            if(splitted_reg_ids.find(curr_node_idx)!=splitted_reg_ids.end() ||
                    merged_reg_ids.find(curr_node_idx)!=splitted_reg_ids.end()){
                continue;
            }
            size_t seeds[2];
            bool merge_flag = false, split_flag = false;
            switch (get<1>(merged_split_peers[i])){
            case MERGE_PARENTS:
                merge_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].parents[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].parents[1];
                if(merged_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        merged_reg_ids.find(seeds[1])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[1])!=splitted_reg_ids.end()){
                    merge_flag = false;
                }
                break;
            case MERGE_KIDS:
                merge_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].kids[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].kids[1];
                if(merged_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        merged_reg_ids.find(seeds[1])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[1])!=splitted_reg_ids.end()){
                    merge_flag = false;
                }
                break;
            case SPLIT_BY_PARENTS:
                split_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].parents[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].parents[1];
                if(merged_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        merged_reg_ids.find(seeds[1])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[1])!=splitted_reg_ids.end()){
                    split_flag = false;
                }
                break;
            case SPLIT_BY_KIDS:
                split_flag = true;
                seeds[0] = movieInfo.nodes[curr_node_idx].kids[0];
                seeds[1] = movieInfo.nodes[curr_node_idx].kids[1];
                if(merged_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        merged_reg_ids.find(seeds[1])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[0])!=splitted_reg_ids.end() ||
                        splitted_reg_ids.find(seeds[1])!=splitted_reg_ids.end()){
                    split_flag = false;
                }
                break;
            default:
                break;
            }
//            if(seeds[0] == 3792 || seeds[1] == 3792){
//                qDebug("Possible error location");
//            }
            if(merge_flag){
                if(mergedRegionGrow(cellSegment, seeds, newCells)){
//                    if(find(uptCell_idxs.begin(), uptCell_idxs.end(), seeds[0]) != uptCell_idxs.end() ||
//                            find(uptCell_idxs.begin(), uptCell_idxs.end(), seeds[1]) != uptCell_idxs.end()){
//                        qDebug("duplicated update found!");
//                    }
                    uptCell_idxs.push_back(seeds[0]);
                    uptCell_idxs.push_back(seeds[1]);
                    merged_reg_ids.insert(seeds[0]);
                    merged_reg_ids.insert(seeds[1]);
                    // if success, we need to check if gap should be removed
                    if(get<2>(merged_split_peers[i]) < 0.05){ // we are confident with the merging, so remove gaps between regions
                        int frame = movieInfo.frames[seeds[0]];
                        setValMat(validGapMaps[frame], CV_8U, newCells[newCells.size()-1].voxIdx, 0);
                    }
                }
            }else if(split_flag){
                bool oneSeedIsGood;
                if(!separateRegion(cellSegment, curr_node_idx, seeds, oneSeedIsGood, newCells)){
                    //curr_node_idx (1) fails to be split and (2) none of its parents can link to it (cost > obzCost)
                    // such case, we merge two parents
                    if(!oneSeedIsGood){ // if there is one seed is good, we have linked the region to one of the seed
                        if(mergedRegionGrow(cellSegment, seeds, newCells)){
//                            if(find(uptCell_idxs.begin(), uptCell_idxs.end(), seeds[0]) != uptCell_idxs.end() ||
//                                    find(uptCell_idxs.begin(), uptCell_idxs.end(), seeds[1]) != uptCell_idxs.end()){
//                                qDebug("duplicated update found!");
//                            }
                            uptCell_idxs.push_back(seeds[0]);
                            uptCell_idxs.push_back(seeds[1]);
                            merged_reg_ids.insert(seeds[0]);
                            merged_reg_ids.insert(seeds[1]);
                            // for the un-seperable region, we do not want to frequently check them; so remove the gaps between them
                            int frame = movieInfo.frames[seeds[0]];
                            setValMat(validGapMaps[frame], CV_8U, newCells[newCells.size()-1].voxIdx, 0);
                        }
                    }
                }else{
//                    if(find(uptCell_idxs.begin(), uptCell_idxs.end(), curr_node_idx) != uptCell_idxs.end() ){
//                        qDebug("duplicated update found!");
//                    }
                    uptCell_idxs.push_back(curr_node_idx);
                    splitted_reg_ids.insert(curr_node_idx);
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
    if(movieInfo.voxIdx[node_idx].size() == 0) return false;
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
/**
 * @brief appendNewCellOrNode:append a new cell to the frame
 * @param cellSegment
 * @param newCells
 * @param append_loc_idx: the location to append in voxIdx
 */
void cellTrackingMain::appendNewCellOrNode(cellSegmentMain &cellSegment, newFoundCellInfo &newCell,
                                           size_t append_loc_idx){
    // 1. update movieInfo
    int frame = newCell.frame;
    int labelInLabelMap;
    if (frame>0) labelInLabelMap = append_loc_idx + 1 - cumulative_cell_nums[frame-1];
    else labelInLabelMap = append_loc_idx + 1;

//    if(frame == 12 && labelInLabelMap==6){
//        qDebug("check point");
//    }
    movieInfo.frames[append_loc_idx] = newCell.frame;
    movieInfo.voxIdx[append_loc_idx] = newCell.voxIdx;
    movieInfo.labelInMap[append_loc_idx] = labelInLabelMap;

    vector<int> y, x, z;
    vec_ind2sub(newCell.voxIdx, y, x, z, cellSegment.cell_label_maps[0].size);
    movieInfo.vox_y[append_loc_idx] = y;
    movieInfo.vox_x[append_loc_idx] = x;
    movieInfo.vox_z[append_loc_idx] = z;

    movieInfo.xCoord[append_loc_idx] = vec_mean(x);
    movieInfo.yCoord[append_loc_idx] = vec_mean(y);
    movieInfo.zCoord[append_loc_idx] = vec_mean(z);

    movieInfo.start_coord_xyz[append_loc_idx] = {vec_min(x), vec_min(y), vec_min(z)};

    movieInfo.range_xyz[append_loc_idx] = {vec_max(x)-movieInfo.start_coord_xyz[append_loc_idx][0] + 1,
                                           vec_max(y)-movieInfo.start_coord_xyz[append_loc_idx][1] + 1,
                                           vec_max(z)-movieInfo.start_coord_xyz[append_loc_idx][2] + 1};
    if(vec_max(movieInfo.range_xyz[append_loc_idx]) > 250){
        qDebug("possibly lose information");
    }
    nodeInfo n;
    n.node_id = append_loc_idx;
    n.in_cost = p4tracking.c_en;
    n.out_cost = p4tracking.c_ex;
    n.detect_confidence = p4tracking.observationCost;
    n.stable_status = NOT_STABLE;
    // update the new neighbors
    vector<size_t> nei_ids(0);
    extractNeighborIds(cellSegment.cell_label_maps, append_loc_idx, nei_ids);
    movieInfo.nodes[append_loc_idx].neighbors.resize(nei_ids.size());
    for(int i=0; i<nei_ids.size(); i++){
        movieInfo.nodes[append_loc_idx].neighbors[i].node_id = nei_ids[i];
        float max_dist = voxelwise_avg_distance(append_loc_idx, nei_ids[i],
                                                movieInfo.nodes[append_loc_idx].neighbors[i].dist_c2n,
                                                movieInfo.nodes[append_loc_idx].neighbors[i].dist_n2c);
        int frame_diff = movieInfo.frames[nei_ids[i]] - newCell.frame;
        movieInfo.nodes[append_loc_idx].neighbors[i].link_cost = distance2cost(max_dist, p4tracking.jumpCost[frame_diff-1]);

        movieInfo.nodes[append_loc_idx].neighbors[i].overlap_size = cellOverlapSize(append_loc_idx, nei_ids[i],
                                                                                    cellSegment);
        bool existing_flag = false;
        for(nodeRelation &preNei : movieInfo.nodes[nei_ids[i]].preNeighbors){
            if(preNei.node_id == append_loc_idx){
                preNei.dist_n2c = movieInfo.nodes[append_loc_idx].neighbors[i].dist_c2n;
                preNei.dist_c2n = movieInfo.nodes[append_loc_idx].neighbors[i].dist_n2c;
                preNei.link_cost = movieInfo.nodes[append_loc_idx].neighbors[i].link_cost;
                preNei.overlap_size = movieInfo.nodes[append_loc_idx].neighbors[i].overlap_size;
                existing_flag = true;
                break;
            }
        }
        if(!existing_flag){
            nodeRelation newPreNei;
            newPreNei.node_id = append_loc_idx;
            newPreNei.dist_n2c = movieInfo.nodes[append_loc_idx].neighbors[i].dist_c2n;
            newPreNei.dist_c2n = movieInfo.nodes[append_loc_idx].neighbors[i].dist_n2c;
            newPreNei.link_cost = movieInfo.nodes[append_loc_idx].neighbors[i].link_cost;
            newPreNei.overlap_size = movieInfo.nodes[append_loc_idx].neighbors[i].overlap_size;
            movieInfo.nodes[nei_ids[i]].preNeighbors.push_back(newPreNei);
        }
    }
    movieInfo.overall_neighbor_num += nei_ids.size();
    // update the preNeighbors
    vector<size_t> preNei_ids(0);
    extractPreNeighborIds(cellSegment.cell_label_maps, append_loc_idx, preNei_ids);
    movieInfo.nodes[append_loc_idx].preNeighbors.resize(preNei_ids.size());
    for(int i=0; i<preNei_ids.size(); i++){
        movieInfo.nodes[append_loc_idx].preNeighbors[i].node_id = preNei_ids[i];
        float max_dist = voxelwise_avg_distance(append_loc_idx, preNei_ids[i],
                                                movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_c2n,
                                                movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_n2c);
        int frame_diff = newCell.frame - movieInfo.frames[preNei_ids[i]];
        movieInfo.nodes[append_loc_idx].preNeighbors[i].link_cost = distance2cost(max_dist, p4tracking.jumpCost[frame_diff-1]);

        movieInfo.nodes[append_loc_idx].preNeighbors[i].overlap_size = cellOverlapSize(append_loc_idx, preNei_ids[i],
                                                                                    cellSegment);
        bool existing_flag = false;
        for(nodeRelation &nei : movieInfo.nodes[preNei_ids[i]].neighbors){
            if(nei.node_id == append_loc_idx){
                nei.dist_n2c = movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_c2n;
                nei.dist_c2n = movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_n2c;
                nei.link_cost = movieInfo.nodes[append_loc_idx].preNeighbors[i].link_cost;
                nei.overlap_size = movieInfo.nodes[append_loc_idx].preNeighbors[i].overlap_size;
                existing_flag = true;
                break;
            }
        }
        if(!existing_flag){
            nodeRelation newNei;
            newNei.node_id = append_loc_idx;
            newNei.dist_n2c = movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_c2n;
            newNei.dist_c2n = movieInfo.nodes[append_loc_idx].preNeighbors[i].dist_n2c;
            newNei.link_cost = movieInfo.nodes[append_loc_idx].preNeighbors[i].link_cost;
            newNei.overlap_size = movieInfo.nodes[append_loc_idx].preNeighbors[i].overlap_size;
            movieInfo.nodes[preNei_ids[i]].neighbors.push_back(newNei);
        }
    }
    // 2. update cellSegment
//    if((newCell.voxIdx[2] == 751285 && newCell.voxIdx.size()==300) || append_loc_idx == 6182){
//        qDebug("check point");
//    }
    setValMat(cellSegment.cell_label_maps[frame], CV_32S, newCell.voxIdx, labelInLabelMap);
    setValMat(cellSegment.threshold_maps[frame], CV_8U, newCell.voxIdx, newCell.threshold);
    if(labelInLabelMap > cellSegment.number_cells[frame]){
        assert(labelInLabelMap == cellSegment.number_cells[frame]+1);
        cellSegment.number_cells[frame] = labelInLabelMap;
    }
}
/**
 * @brief nullifyCellOrNode: nullify a cell or a node from movieInfo and cellSegment
 * @param node_idx
 */
void cellTrackingMain::nullifyCellOrNode(size_t node_idx, cellSegmentMain *cellSegment){
    // 1. remove from cell_label_maps (threshold_maps may be not needed)
    // NOTE: we do not need to decrease cell_num since the max cell label in the frame is not changed (not for sure)
    if(cellSegment != nullptr){
        int frame = movieInfo.frames[node_idx];
        //qDebug("%d val", cellSegment->cell_label_maps[frame].at<int>(movieInfo.voxIdx[node_idx][0]));
        setValMat(cellSegment->cell_label_maps[frame], CV_32S, movieInfo.voxIdx[node_idx], 0);
        //qDebug("%d val", cellSegment->cell_label_maps[frame].at<int>(movieInfo.voxIdx[node_idx][0]));
        setValMat(cellSegment->threshold_maps[frame], CV_8U, movieInfo.voxIdx[node_idx], 0);
    }
    // 2. nullify movieInfo
    movieInfo.voxIdx[node_idx].resize(0);
    movieInfo.vox_y[node_idx].resize(0);
    movieInfo.vox_x[node_idx].resize(0);
    movieInfo.vox_z[node_idx].resize(0);
    movieInfo.xCoord[node_idx] = -INFINITY;
    movieInfo.yCoord[node_idx] = -INFINITY;
    movieInfo.zCoord[node_idx] = -INFINITY;
    movieInfo.start_coord_xyz[node_idx] = {-1,-1,-1};
    movieInfo.range_xyz[node_idx] = {0,0,0};
    movieInfo.labelInMap[node_idx] = 0;
    movieInfo.node_tested_st_end_jump[node_idx].track_head_tested = false;
    movieInfo.node_tested_st_end_jump[node_idx].track_tail_tested = false;
    movieInfo.node_tested_st_end_jump[node_idx].jump_tested = false;

    for(nodeRelation n : movieInfo.nodes[node_idx].neighbors){
        for(nodeRelation &pn : movieInfo.nodes[n.node_id].preNeighbors){
            if(pn.node_id == node_idx){
                pn.dist_c2n = INFINITY;
                pn.dist_n2c = INFINITY;
                pn.link_cost = INFINITY;
                pn.overlap_size = 0;
                break;
            }
        }
    }
    movieInfo.overall_neighbor_num -= movieInfo.nodes[node_idx].neighbors.size();
    movieInfo.nodes[node_idx].neighbors.resize(0);

    for(nodeRelation pn : movieInfo.nodes[node_idx].preNeighbors){
        for(nodeRelation &n : movieInfo.nodes[pn.node_id].neighbors){
            if(n.node_id == node_idx){
                n.dist_c2n = INFINITY;
                n.dist_n2c = INFINITY;
                n.link_cost = INFINITY;
                n.overlap_size = 0;
                break;
            }
        }
    }
    movieInfo.nodes[node_idx].preNeighbors.resize(0);

}
void cellTrackingMain::nullifyCellOrNode(size_t node_idx[], cellSegmentMain *cellSegment){
    int num = sizeof(node_idx)/sizeof(node_idx[0]);
    for(int i=0; i<num; i++){
        nullifyCellOrNode(node_idx[i], cellSegment);
    }
}
void cellTrackingMain::nullifyCellOrNode(vector<size_t> node_idx, cellSegmentMain *cellSegment){
    for(size_t n : node_idx){
        nullifyCellOrNode(n, cellSegment);
    }
}
/**
 * @brief separateRegion
 * @param cellSegment
 * @param node_idx
 * @param seeds
 */
bool cellTrackingMain::separateRegion(cellSegmentMain &cellSegment, size_t node_idx,
                                      size_t seeds[2], bool &oneSeedIsGood, vector<newFoundCellInfo> &newCells){
    //newCells.resize(0);
    // first check that the regions is valid
    if(isValid(node_idx, &cellSegment)){
        oneSeedIsGood = false;
        bool split_success = false;
        vector<vector<size_t>> splitRegs(2);
        float reg4seeds2splitRes_costs[2];
        bool gapBasedSplit = false, usePriorGapMap = true;
        if(p4tracking.stableNodeTest == false || movieInfo.nodes[node_idx].stable_status == NOT_STABLE){
            vector<vector<size_t>> reg4seeds(2);
            reg4seeds[0] = movieInfo.voxIdx[seeds[0]];
            reg4seeds[1] = movieInfo.voxIdx[seeds[1]];
            int reg4seeds_frame = movieInfo.frames[seeds[0]];
            gapBasedSplit = true;
            split_success = bisectValidTest(cellSegment, node_idx, movieInfo.voxIdx[node_idx], movieInfo.frames[node_idx],
                                reg4seeds, reg4seeds_frame, gapBasedSplit, usePriorGapMap,
                                 splitRegs, reg4seeds2splitRes_costs);
            if(!split_success){
                splitRegs.clear();
                gapBasedSplit = false;
                split_success = bisectValidTest(cellSegment, node_idx, movieInfo.voxIdx[node_idx], movieInfo.frames[node_idx],
                                    reg4seeds, reg4seeds_frame, gapBasedSplit, usePriorGapMap,
                                     splitRegs, reg4seeds2splitRes_costs);
            }
        }
        if(!split_success){ // fail to split: let's check if merge is a possible solution
            if(handleNonSplitReg_link2oneSeed(node_idx, seeds)){
                oneSeedIsGood = true;
            }
            return false;
        }else{ // split success
            // 1. nullify the current region
            // nullifyCellOrNode(node_idx); // leave it to movieIinfo_update
            unsigned char thres_cur =
                    cellSegment.threshold_maps[movieInfo.frames[node_idx]].at<unsigned char>(movieInfo.voxIdx[node_idx][0]);
            assert(thres_cur > 0);
            // 2. save the two new regions in
            newFoundCellInfo newC1, newC2;
            newCells.push_back(newC1);
            newCells.push_back(newC2);
            newCells[newCells.size() - 2].frame = movieInfo.frames[node_idx];
            newCells[newCells.size() - 2].threshold = thres_cur;
            newCells[newCells.size() - 2].voxIdx = splitRegs[0];
            newCells[newCells.size() - 1].frame = movieInfo.frames[node_idx];
            newCells[newCells.size() - 1].threshold = thres_cur;
            newCells[newCells.size() - 1].voxIdx = splitRegs[1];
            return true;
        }
    }else{
        qDebug("This region or seeds has been nullified!");
        return false;
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
        if(isBestNeighbor(node_idx, seeds[i], curr_link_cost)){
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
bool cellTrackingMain::mergedRegionGrow(cellSegmentMain &cellSegment, size_t seeds[2], vector<newFoundCellInfo> &newCells){
    assert(movieInfo.frames[seeds[0]] == movieInfo.frames[seeds[1]]);
    int curr_frame = movieInfo.frames[seeds[0]];
    if(p4tracking.stableNodeTest && movieInfo.nodes[seeds[0]].stable_status == NOT_STABLE &&
            movieInfo.nodes[seeds[1]].stable_status == NOT_STABLE){
        if(isValid(seeds[0], &cellSegment) && isValid(seeds[1], &cellSegment)){
            newFoundCellInfo newC;
            newCells.push_back(newC);
            newCells[newCells.size() - 1].frame = curr_frame;
            newCells[newCells.size() - 1].voxIdx = movieInfo.voxIdx[seeds[0]];
            newCells[newCells.size() - 1].voxIdx.insert(newCells[newCells.size() - 1].voxIdx.end(),
                    movieInfo.voxIdx[seeds[1]].begin(), movieInfo.voxIdx[seeds[1]].end());

            float r1 = ((float)movieInfo.voxIdx[seeds[0]].size()) / newCells[newCells.size() - 1].voxIdx.size();
            float r2 = ((float)movieInfo.voxIdx[seeds[1]].size()) / newCells[newCells.size() - 1].voxIdx.size();
            newCells[newCells.size() - 1].threshold = (unsigned char)(
                        r1 * cellSegment.threshold_maps[curr_frame].at<unsigned char>(movieInfo.voxIdx[seeds[0]][0]) +
                    r2 * cellSegment.threshold_maps[curr_frame].at<unsigned char>(movieInfo.voxIdx[seeds[1]][0])
                    );

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
    vector<newFoundCellInfo> newCells;
    vector<size_t> uptCell_idxs;
    regionRefresh(cellSegment, newCells, uptCell_idxs);

    // 3. update movieInfo/cellSegment based on the newly detected cells or removed cells
    movieInfo_update(cellSegment, newCells, uptCell_idxs);
}
/**
 * @brief movieInfo_update:update movieInfo/cellSegment based on the newly detected cells or removed cells
 * @param cellSegment
 * @param newCells
 * @param uptCell_idxs
 */
void cellTrackingMain::movieInfo_update(cellSegmentMain &cellSegment, vector<newFoundCellInfo> &newCells,
                                        vector<size_t> &uptCell_idxs){
    if(vec_unique(uptCell_idxs)){
        qFatal("There is duplicated updating! check first func: regionRefresh() first");
    }
    vector<unordered_set<size_t>> uptCell_framewise (cellSegment.number_cells.size());
    //vector<unordered_set<size_t>> newCells_framewise (cellSegment.number_cells.size());
    FOREACH_i(uptCell_idxs){
        int frame = movieInfo.frames[uptCell_idxs[i]];
        qInfo("remove node_id %ld, with label %d, threshold %d", uptCell_idxs[i], movieInfo.labelInMap[uptCell_idxs[i]],
                cellSegment.threshold_maps[frame].at<unsigned char>(movieInfo.voxIdx[uptCell_idxs[i]][0]));

        uptCell_framewise[frame].insert(uptCell_idxs[i]);
        nullifyCellOrNode(uptCell_idxs[i], &cellSegment);
    }
    FOREACH_i(newCells){
        newFoundCellInfo &sNi = newCells[i];
        size_t new_idx;
        if(uptCell_framewise[sNi.frame].size() > 0){
            new_idx = *uptCell_framewise[sNi.frame].begin();
            uptCell_framewise[sNi.frame].erase(uptCell_framewise[sNi.frame].begin());
        }else{
            new_idx = cellSegment.number_cells[sNi.frame] + 1;
            new_idx -= 1; // save index starts from 0
            new_idx += sNi.frame>0 ? cumulative_cell_nums[sNi.frame-1]:0;
        }

        appendNewCellOrNode(cellSegment, sNi, new_idx);

        qInfo("add node_id %ld, with label %d, threshold %d", new_idx, movieInfo.labelInMap[new_idx],
                cellSegment.threshold_maps[sNi.frame].at<unsigned char>(movieInfo.voxIdx[new_idx][0]));
//        if(new_idx == 6182){
//            qDebug("check piont");
//        }
//        if(movieInfo.voxIdx[6182].size() == 300 && cellSegment.cell_label_maps[24].at<int>(751035) != 155){
//            qDebug("check piont");
//        }
    }
    newCells.clear();
    uptCell_idxs.clear();
}
void cellTrackingMain::missing_cell_module(cellSegmentMain &cellSegment){
    // step 1. link cells allowing jump but no split/merge
    mccTracker_one2one();

    // step 2. re-detect missed cells (manipulate existing cells if needed)
    vector<newFoundCellInfo> newCells;
    vector<size_t> uptCell_idxs;
    retrieve_missing_cells(cellSegment, newCells, uptCell_idxs);

    // step 3. update movieInfo/cellSegment based on the newly detected cells or removed cells
    // movieInfo_update(cellSegment, newCells, uptCell_idxs);//Has been moved to retrieve_missing_cells
}
/**
 * @brief retrieve_missing_cells:
 * @param cellSegment
 */
void cellTrackingMain::retrieve_missing_cells(cellSegmentMain &cellSegment, vector<newFoundCellInfo> &newCells,
                                              vector<size_t> &uptCell_idxs){
    // check the cell one by one
    for(size_t i=0; i<movieInfo.nodes.size(); i++){
        if(i % 1000 == 0){
            qDebug("--------------processing %ld / %ld nodes.--------------------", i, movieInfo.nodes.size());
        }
        if(!isValid(i, nullptr)) continue;
        int cur_frame = movieInfo.frames[i];
        if(movieInfo.nodes[i].kid_num > 0){
            assert(movieInfo.nodes[i].kid_num == 1);
            size_t kid_idx = movieInfo.nodes[i].kids[0];
            if(movieInfo.frames[kid_idx] - cur_frame > 1){
                if (!(movieInfo.node_tested_st_end_jump[i].jump_tested
                        && movieInfo.node_tested_st_end_jump[i].region_id_jumped_to == kid_idx)){
                    if(!deal_single_missing_case(cellSegment, newCells, uptCell_idxs, i, MISS_AS_JUMP)){
                        movieInfo.node_tested_st_end_jump[i].jump_tested = true;
                        movieInfo.node_tested_st_end_jump[i].region_id_jumped_to = kid_idx;
                    }
                }
            }
        }
        if(p4tracking.detect_missing_head_tail){
            // tail extending
            if(movieInfo.nodes[i].kid_num == 0 && movieInfo.nodes[i].parent_num > 0 &&
                    cur_frame < cellSegment.number_cells.size()-1 &&
                    !movieInfo.node_tested_st_end_jump[i].track_tail_tested){
                if(!deal_single_missing_case(cellSegment, newCells, uptCell_idxs, i, MISS_AT_TRACK_END)){
                    movieInfo.node_tested_st_end_jump[i].track_tail_tested = true;
                }
            }
            // head extending
            if(movieInfo.nodes[i].parent_num == 0 && movieInfo.nodes[i].kid_num > 0 &&
                    cur_frame > 0 && !movieInfo.node_tested_st_end_jump[i].track_head_tested){
                if(!deal_single_missing_case(cellSegment, newCells, uptCell_idxs, i, MISS_AT_TRACK_START)){
                    movieInfo.node_tested_st_end_jump[i].track_head_tested = true;
                }
            }
        }

        /** directly update the movieInfo here **/
        movieInfo_update(cellSegment, newCells, uptCell_idxs);

//        if(cellSegment.cell_label_maps[movieInfo.frames[112846]].at<int>(movieInfo.voxIdx[112846][0]) == 0){
//            qFatal("found a wrong place %ld", i);
//        }
    }
}
/**
 * @brief deal_single_missing_case: single case of missing retrieval
 * @param cellSegment
 * @param newCells
 * @param uptCell_idxs
 * @param missing_type
 */
bool cellTrackingMain::deal_single_missing_case(cellSegmentMain &cellSegment, vector<newFoundCellInfo> &newCells,
                                                vector<size_t> &uptCell_idxs, size_t cur_node_idx, int missing_type){
    vector<int> missing_frames;
    vector<size_t> seed_loc_idx;
    size_t parent_idx = cur_node_idx;
    size_t kid_idx = missing_type == MISS_AS_JUMP ? movieInfo.nodes[parent_idx].kids[0] : 0;
    if (missing_type == MISS_AT_TRACK_START){
        size_t tmp = parent_idx;
        parent_idx = kid_idx;
        kid_idx = tmp;
    }
    if(!extractSeedFromGivenCell(cellSegment, missing_type, parent_idx,
                                 kid_idx, missing_frames, seed_loc_idx)){
        return false;
    }
//    if(parent_idx == 28528 && kid_idx == 65841){
//        qDebug("check point");
//    }

    bool detected_missing = false;
    // start the big game for missing cell retrieve
    for(int frame : missing_frames){
        vector<vector<size_t>> valid_seeds_loc_idx(0);
        vector<pair<size_t, int>> seeds_missing_type (0);
        unordered_set<size_t> cell_idx_can_be_removed (0);
        vector<pair<size_t, float>> removed_cell_idx_threshold (0);
        //vector<nodeRelation> removed_cells (0);
        vector<int> seed_loc_existing_labels = extractValsGivenIdx_type<int>(&cellSegment.cell_label_maps[frame], seed_loc_idx, CV_32S);
        vector<size_t> valid_seed_loc_idx = vec_atrange(seed_loc_idx, seed_loc_existing_labels, 0, 0, false);
        size_t min_seed_sz = cellSegment.p4segVol.min_seed_size;
        if(!p4tracking.addCellMissingPart){
            min_seed_sz = max(min_seed_sz, (1+seed_loc_idx.size())/2);
            if(valid_seed_loc_idx.size() < min_seed_sz){
                valid_seed_loc_idx.clear();
                checkSeedCoveredByExistingCell(cellSegment, missing_type, parent_idx, kid_idx, frame, min_seed_sz,
                                               seed_loc_idx, seed_loc_existing_labels, cell_idx_can_be_removed,
                                               seeds_missing_type, valid_seeds_loc_idx);
                // NOTE: these removed cells may be too much, we may add them back later
                if(cell_idx_can_be_removed.size() > 0){
                    for(auto it : cell_idx_can_be_removed){// nullify the region in label maps
                        removed_cell_idx_threshold.push_back(make_pair(it,
                                                 cellSegment.threshold_maps[frame].at<unsigned char>(movieInfo.voxIdx[it][0])));
                        setValMat(cellSegment.cell_label_maps[frame], CV_32S, movieInfo.voxIdx[it], 0);
                        setValMat(cellSegment.threshold_maps[frame], CV_8U, movieInfo.voxIdx[it], 0);
                    }
                }
            }else{
                valid_seeds_loc_idx.push_back(valid_seed_loc_idx);
                seeds_missing_type.push_back(make_pair(cur_node_idx, missing_type));
            }
        }else if(valid_seed_loc_idx.size() >= min_seed_sz){
            valid_seeds_loc_idx.push_back(valid_seed_loc_idx);
            seeds_missing_type.push_back(make_pair(cur_node_idx, missing_type));
        }else{
            continue;
        }
        vector<size_t> seed_idx4fgRefine;
        for(int i=0; i<valid_seeds_loc_idx.size(); i++){
            seed_idx4fgRefine.insert(seed_idx4fgRefine.end(), valid_seeds_loc_idx[i].begin(),
                                     valid_seeds_loc_idx[i].end());
        }

        /** debug*/
        //        vector<int> tmp = extractValsGivenIdx_type<int>(&cellSegment.cell_label_maps[frame], seed_idx4fgRefine, CV_32S);
        //        vector<size_t> invalid_locidx = vec_atrange(seed_idx4fgRefine, tmp, 100000000, 0, true);
        //        if(invalid_locidx.size() > 0){
        //            qDebug("we set some cell with existing labels to 0");
        //        }
        //---debug end --//
        int seed_label4fgRefine = cellSegment.number_cells[frame] + 1;// make sure no existing cell has this label
        setValMat(cellSegment.cell_label_maps[frame], CV_32S, seed_idx4fgRefine, seed_label4fgRefine);
        //vector<newFoundCellInfo> newCells_currentRun;
        size_t parentKid_idx[] = {parent_idx, kid_idx};
        size_t init_uptCellSize = uptCell_idxs.size();
        if(!redetectCellinTrackingwithSeed(cellSegment, seed_idx4fgRefine, seed_label4fgRefine, frame,
                                           valid_seeds_loc_idx, parentKid_idx, missing_type, newCells, uptCell_idxs)){
            if(removed_cell_idx_threshold.size() > 0){
                for(auto it : removed_cell_idx_threshold){
                    setValMat(cellSegment.cell_label_maps[frame], CV_32S, movieInfo.voxIdx[it.first], movieInfo.labelInMap[it.first]);
                    setValMat(cellSegment.threshold_maps[frame], CV_8U, movieInfo.voxIdx[it.first], it.second);
                }
            }
//        }else{
//            if(newCells_currentRun.size() == 1){

//            }else{

//            }
        }else{
            if(uptCell_idxs.size() > init_uptCellSize){ // partially nullify the region
                for(size_t rdr = init_uptCellSize; rdr<uptCell_idxs.size(); rdr++){
                    setValMat(cellSegment.cell_label_maps[frame], CV_32S, movieInfo.voxIdx[uptCell_idxs[rdr]], 0);
                    setValMat(cellSegment.threshold_maps[frame], CV_8U, movieInfo.voxIdx[uptCell_idxs[rdr]], 0);
                }
            }
            if(cell_idx_can_be_removed.size() > 0){
//                for(auto cell_iid : cell_idx_can_be_removed){
//                    if(find(uptCell_idxs.begin(), uptCell_idxs.end(), cell_iid) != uptCell_idxs.end()){
//                        qDebug("duplicated update found!");
//                    }
//                }
                uptCell_idxs.insert(uptCell_idxs.end(), cell_idx_can_be_removed.begin(), cell_idx_can_be_removed.end());
            }
            detected_missing = true;
        }
        // Set the seed region back to 0s, NOTE:!!! this place may contradict with line-3457,
        // so we need to check, rather than using: setValMat(cellSegment.cell_label_maps[frame], CV_32S, seed_idx4fgRefine, 0);
        FOREACH_i(seed_idx4fgRefine){
            if(cellSegment.cell_label_maps[frame].at<int>(seed_idx4fgRefine[i]) == seed_label4fgRefine){
                cellSegment.cell_label_maps[frame].at<int>(seed_idx4fgRefine[i]) = 0;
            }
        }
    }
    return detected_missing;
}
/**
 * @brief redetectCellinTrackingwithSeed: contains 3 functions in matlab, which are redetectCellinTrackingwithSeed +
 * testMissingCellGotFromOneSeed + testMissingCellGotFromMultiSeeds
 * @param cellSegment
 * @param seed_idx4fgRefine
 * @param seed_label_in_map
 * @param frame
 * @param valid_seeds_loc_idx
 * @param parentKid_idx
 * @param missing_type
 * @param newCells
 * @return
 */
bool cellTrackingMain::redetectCellinTrackingwithSeed(cellSegmentMain &cellSegment, vector<size_t> seed_idx4fgRefine,
                                    int seed_label_in_map, int frame, vector<vector<size_t>> valid_seeds_loc_idx,
                                    size_t parentKid_idx[2], int missing_type, vector<newFoundCellInfo> &newCells,
                                    vector<size_t> &uptCell_idxs){
    if(seed_idx4fgRefine.size() < cellSegment.p4segVol.min_seed_size){
        return false;
    }
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];
    unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*frame; // sub-matrix pointer
    Mat data_grayim3d(3, cellSegment.normalized_data4d.size, CV_8U, ind);
    singleCellSeed seed;
    cellSegment.cropSeed(seed_label_in_map, seed_idx4fgRefine, &data_grayim3d, nullptr, &cellSegment.cell_label_maps[frame],
            frame, seed, cellSegment.p4segVol);
    int touchBnd = cellSegment.refineSeed2Region(seed, cellSegment.p4odStats, cellSegment.p4segVol);

    if(touchBnd != NO_TOUCH){ // region too small
        if(touchBnd == XY_TOUCH || touchBnd == XYZ_TOUCH){
            cellSegment.p4segVol.shift_yxz[0] *=2;
            cellSegment.p4segVol.shift_yxz[1] *=2;
        }
        if(touchBnd == Z_TOUCH || touchBnd == XYZ_TOUCH){
            cellSegment.p4segVol.shift_yxz[2] *=2;
        }
        Range new_range[3];
        getRange(seed.y, cellSegment.p4segVol.shift_yxz[0], cellSegment.cell_label_maps[frame].size[0], new_range[0]);
        getRange(seed.x, cellSegment.p4segVol.shift_yxz[1], cellSegment.cell_label_maps[frame].size[1], new_range[1]);
        getRange(seed.z, cellSegment.p4segVol.shift_yxz[2], cellSegment.cell_label_maps[frame].size[2], new_range[2]);
        if(new_range[0].size() > seed.validSearchAreaMap.size[0] ||
                new_range[1].size() > seed.validSearchAreaMap.size[1] ||
                new_range[2].size() > seed.validSearchAreaMap.size[2]){

            cellSegment.cropSeed(seed_label_in_map, seed_idx4fgRefine, &data_grayim3d, nullptr,
                                 &cellSegment.cell_label_maps[frame], frame, seed, cellSegment.p4segVol);
            cellSegment.refineSeed2Region(seed, cellSegment.p4odStats, cellSegment.p4segVol);
        }
        cellSegment.reset_shift();
    }

//    Mat1b tmp0=seed.otherIdMap > 0, tmp1 = seed.outputIdMap > 0, tmp2;
//    bitwise_and(tmp0, tmp1, tmp2);
//    if(!isempty(tmp2, CV_8U)){
//        ccShowSlice3Dmat(tmp2, CV_8U);
//        ccShowSliceLabelMat(seed.outputIdMap);
//        ccShowSliceLabelMat(seed.idMap);
//        qFatal("region grow is problematic");
//    }

    if(seed.bestFgThreshold <= 0 || isempty(seed.outputIdMap, CV_32S)){
        return false;
    }
    //ccShowSliceLabelMat(seed.outputIdMap);
    //// testMissingCellGotFromMultiSeeds: if we detected multiple cells with
    ///  several seeds, we only test if they are compatible with its own
    ///  parents/kids, if so, keep them, otherwise remove those bad ones.
    if(valid_seeds_loc_idx.size() > 1){
        /** step 1. split the fg into different cells */
        Mat1i label_map = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_32S);
        for(int i=0; i<valid_seeds_loc_idx.size(); i++){
            for(size_t tmp_idx : valid_seeds_loc_idx[i]){
                label_map.at<int>(tmp_idx) = i+1;
            }
        }
        Mat1b fgMap = seed.outputIdMap>0;
        Mat1i grownSeedMap2d, grownSeedMap3d;
        bool bg2sink = true;
        //ccShowSliceLabelMat(seedMap);
        regionGrow(&label_map, 2, grownSeedMap2d, &seed.score2d, &fgMap,
                   cellSegment.p4segVol.growConnectInTest, cellSegment.p4segVol.graph_cost_design, bg2sink);
        //ccShowSliceLabelMat(grownSeedMap2d);
        bg2sink = false;
        regionGrow(&grownSeedMap2d, 2, grownSeedMap3d, &seed.scoreMap, &fgMap,
                   cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, bg2sink);
        //ccShowSliceLabelMat(grownSeedMap3d);
        //grownSeedMap3d.copyTo(seed.idMap);
        /** step 2. extract idx and save them into newCells*/
        vector<vector<size_t>> newCells_idx_small(valid_seeds_loc_idx.size());
        extractVoxIdxList(&grownSeedMap3d, newCells_idx_small, (int)valid_seeds_loc_idx.size());
//        void coordinateTransfer(vector<size_t> &in_idx, int org_sz_yxz[3],
//                    vector<size_t> &out_idx, int crop_start_yxz[3], int crop_sz_yxz[3]);
//        void coordinateTransfer(vector<size_t> &in_idx, MatSize org_sz_yxz,
//                    vector<size_t> &out_idx, int crop_start_yxz[3], MatSize crop_sz_yxz);
        int crop_start_yxz[] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start,
                                -seed.crop_range_yxz[2].start};
        bool indeed_add_cell = false;
        for(int i=0; i<newCells_idx_small.size(); i++){
            newFoundCellInfo tmp;
            tmp.frame = frame;
            tmp.threshold = seed.bestFgThreshold;
            coordinateTransfer(newCells_idx_small[i], grownSeedMap3d.size, tmp.voxIdx,
                               crop_start_yxz, cellSegment.cell_label_maps[frame].size);
            if(parentOrKidValidLinkTest(tmp.voxIdx, frame, parentKid_idx, missing_type, cellSegment.cell_label_maps[frame].size)){
                newCells.push_back(tmp);
                indeed_add_cell = true;
            }
        }
        return indeed_add_cell;
    }else{
        // check the neighboring cells, should we combine them?
        /** step 1. find neighboring cells */
        vector<size_t> adj_cell_ids;
        Mat1b fgMap = seed.outputIdMap>0, fgMap_dilate;
        int dilate_sz[] = {1,1,0};
        volumeDilate(&fgMap, fgMap_dilate, dilate_sz, MORPH_CROSS);
        fgMap_dilate = fgMap_dilate - fgMap;
        vector<size_t> adj_cells_idx;
        vector<int> adj_cells_labels = extractValsGivenMask_type<int>(&seed.outputIdMap, CV_32S, &fgMap_dilate, 0);
        unordered_map<int, size_t> adj_cells_adj_sz = frequecy_cnt(adj_cells_labels);
        vector<size_t> adj_cell_id (1);
        if(adj_cells_adj_sz.size() > 0){
            size_t tmp_max_size = 0;
            int adj_cell_label = -1;
            for(auto it : adj_cells_adj_sz){
                if (it.first!=0 && it.second > tmp_max_size){
                    tmp_max_size = it.second;
                    adj_cell_label = it.first;
                }
            }
            if(adj_cell_label > 0){
                adj_cell_id[0] = adj_cell_label;
                adj_cell_id[0] -= 1; // label 1 indicates node id 0 in the vector
                adj_cell_id[0] += frame == 0 ? 0:cumulative_cell_nums[frame-1];
            }else{
                adj_cell_id.resize(0);
            }
        }else{
            adj_cell_id.resize(0);
        }
        vector<vector<size_t>> new_cells_loc_idx;
        extractVoxIdxList(&seed.outputIdMap, new_cells_loc_idx, seed.outCell_num);
        vector<size_t> new_cell_loc_idx, tmp = vec_merge(new_cells_loc_idx);
        int crop_start_yxz[] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start,
                                -seed.crop_range_yxz[2].start};
        coordinateTransfer(tmp, seed.outputIdMap.size, new_cell_loc_idx,
                           crop_start_yxz, cellSegment.cell_label_maps[frame].size);
        vector<vector<size_t>> extra_new_cells_loc_idx;
        pair<int, int> res = newlyAddedCellValidTest(cellSegment, seed, new_cell_loc_idx, frame,
                                                     adj_cell_id, parentKid_idx, missing_type,
                                                     extra_new_cells_loc_idx);
        if(res.second == 0){
            return false;
        }else{
            if(res.first == 2 || res.first == 5){
                newFoundCellInfo tmp;
                tmp.frame = frame;
                tmp.threshold = seed.bestFgThreshold;
                tmp.voxIdx = vec_merge(new_cell_loc_idx, movieInfo.voxIdx[adj_cell_id[0]]);
                newCells.push_back(tmp);
//                if(find(uptCell_idxs.begin(), uptCell_idxs.end(), adj_cell_id[0]) != uptCell_idxs.end()){
//                    qDebug("duplicated update found!");
//                }
                uptCell_idxs.push_back(adj_cell_id[0]);
                movieInfo.nodes[adj_cell_id[0]].node_id = -1;
            }else if(res.first == 3){
                assert(extra_new_cells_loc_idx.size() == 2);
                newFoundCellInfo tmp;
                tmp.frame = frame;
                tmp.threshold = cellSegment.threshold_maps[frame].at<unsigned char>(movieInfo.voxIdx[adj_cell_id[0]][0]);
                tmp.voxIdx = extra_new_cells_loc_idx[0];
                newCells.push_back(tmp);
                tmp.voxIdx = extra_new_cells_loc_idx[1];
                tmp.threshold = seed.bestFgThreshold;
                newCells.push_back(tmp);
//                if(find(uptCell_idxs.begin(), uptCell_idxs.end(), adj_cell_id[0]) != uptCell_idxs.end()){
//                    qDebug("duplicated update found!");
//                }
                uptCell_idxs.push_back(adj_cell_id[0]);
                movieInfo.nodes[adj_cell_id[0]].node_id = -1;
            }else if(res.first == 6){
                assert(extra_new_cells_loc_idx.size() == 1);
                newFoundCellInfo tmp;
                tmp.frame = frame;
                tmp.voxIdx = extra_new_cells_loc_idx[0];
                tmp.threshold = seed.bestFgThreshold;
                newCells.push_back(tmp);
                //uptCell_idxs.push_back(adj_cell_id[0]);
            }else{ // case 1 and 4, we can directly add the new cell
                newFoundCellInfo tmp;
                tmp.frame = frame;
                tmp.voxIdx = new_cell_loc_idx;
                tmp.threshold = seed.bestFgThreshold;
                newCells.push_back(tmp);
                //uptCell_idxs.push_back(adj_cell_id[0]);
            }
            return true;
        }
    }
}
//% there are multiple cases that are valid for an newly found cell
//% 1. a1->miss->a3
//% 2. a1->miss,b->a3
//% 3. a1,b1->miss,b2->a3,b3
//% 4. a1,b1->miss->a3,b3
//% 5. a1+b1->miss,b2->a3+b3
//% 6. a1->miss+noise_bright_region->a3 (try split)
//% 7. a1+b1->a2+b2_half,half_missing->a3,b3
//% PS: '+' means merged region, ',' means independent region
//% case 1 and 4 has not region to merge
pair<int, int> cellTrackingMain::newlyAddedCellValidTest(cellSegmentMain &cellSegment, singleCellSeed &seed, vector<size_t> &new_cell_idx,
                                               int new_cell_frame, vector<size_t> ajd_cell_idx, size_t parentKid_idx[2],
                                                int missing_type, vector<vector<size_t>> &extra_new_cells_loc_idx){
    float baseCost1, baseCost2;
    pair<int, int> res = make_pair(0,0);
    MatSize sz = cellSegment.cell_label_maps[new_cell_frame].size;
    bool parentKidLink_valid = parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx, missing_type, sz, baseCost1, baseCost2);
    if(ajd_cell_idx.size() == 0){
        if(parentKidLink_valid){
            res = make_pair(1,1);
        }else{/*% first check case 6: the region can be split by principal
                                    % curvauture. The processing steps are exactly the same as the
                                    % function 'segmentCurrentRegion.m'.*/
            Mat1i label_map = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_32S);
            Mat1b fgMap = seed.idMap>0;
            Mat1b seedMap = fgMap - seed.gap3dMap;
            int numCC = connectedComponents3d(&seedMap, label_map, cellSegment.p4segVol.neiMap);
            vector<vector<size_t>> seeds_idxes(numCC);
            extractVoxIdxList(&label_map, seeds_idxes, numCC);
            bool overlap_with_seed = false, not_overlap_with_seed = false;
            Mat1i label_map_new = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_32S);
            for(auto &it : seeds_idxes){
                if(it.size() > cellSegment.p4segVol.min_seed_size){
                    bool overlapped = false;
                    for(size_t tmp_idx : it){
                        if(seed.seedMap.at<unsigned char>(tmp_idx) > 0){
                            overlapped = true;
                            break;
                        }
                    }
                    if(overlapped){
                        overlap_with_seed = true;
                        setValMat(label_map_new, CV_32S, it, 1);
                    }else{
                        setValMat(label_map_new, CV_32S, it, 2);
                        not_overlap_with_seed = true;
                    }
                }
            }
            vector<size_t> new_cell_with_append_idx;
            if(not_overlap_with_seed && overlap_with_seed){
                bool bg2sink = false;
                label_map = Mat::zeros(seed.idMap.dims, seed.idMap.size, CV_32S);
                regionGrow(&label_map_new, 2, label_map, &seed.scoreMap, &fgMap,
                           cellSegment.p4segVol.growConnectInRefine, cellSegment.p4segVol.graph_cost_design, bg2sink);
                vector<size_t> tmp_idx;
                extractVoxIdxGivenId(&label_map, tmp_idx, 1);
                int crop_start_yxz[] = {-seed.crop_range_yxz[0].start, -seed.crop_range_yxz[1].start,
                                        -seed.crop_range_yxz[2].start};
                coordinateTransfer(tmp_idx, fgMap.size, new_cell_with_append_idx,
                                   crop_start_yxz, cellSegment.cell_label_maps[new_cell_frame].size);
                parentKidLink_valid = parentOrKidValidLinkTest(new_cell_with_append_idx, new_cell_frame, parentKid_idx, missing_type, sz);
            }
            if(parentKidLink_valid){
                res = make_pair(6,1);
                extra_new_cells_loc_idx.clear();
                extra_new_cells_loc_idx.push_back(new_cell_with_append_idx);
            }else if(missing_type == MISS_AS_JUMP){
                //case 4: TOO RARE, so we only check it loosely with a strict criterion
                vector<size_t> adj_pk_vec(2);
                if(multiNeis_check(cellSegment, parentKid_idx[0], parentKid_idx[1], new_cell_idx, adj_pk_vec)){
                    vector<vector<size_t>> double_p_k_idx (2);
                    double_p_k_idx[0] = {parentKid_idx[0], adj_pk_vec[0]};
                    double_p_k_idx[1] = {parentKid_idx[1], adj_pk_vec[1]};
                    if(parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, double_p_k_idx, missing_type, sz)){
                        res = make_pair(4,1);
                    }else{
                        res = make_pair(4,0);
                    }
                }else{
                    res = make_pair(4,0);
                }

            }else{
                res = make_pair(1,0);
            }
        }
    }else{ // there indeed an adjacent cell, case 2,3,5
        vector<size_t> append_idxPlus = movieInfo.voxIdx[ajd_cell_idx[0]];
        append_idxPlus.reserve(append_idxPlus.size() + new_cell_idx.size());
        append_idxPlus.insert(append_idxPlus.end(), new_cell_idx.begin(), new_cell_idx.end());
        if(movieInfo.nodes[ajd_cell_idx[0]].parent_num > 0 &&
                movieInfo.nodes[ajd_cell_idx[0]].kid_num > 0){
            float mergeCost1, mergeCost2;
            bool merge_link_pk_valid = parentOrKidValidLinkTest(append_idxPlus, new_cell_frame,
                                                                parentKid_idx, missing_type, sz, mergeCost1, mergeCost2);
            if(merge_link_pk_valid){
                if(mergeCost1 <= baseCost1 && mergeCost2 <= baseCost2){
                    res = make_pair(2,1); // merge is better than split
                }else if(parentKidLink_valid){
                    res = make_pair(1,1); // merge not as good as split
                }else{
                    res = make_pair(2,1); // merge is the only choice
                }
            }else if(parentKidLink_valid){
                res = make_pair(1,1); // not merge is the only choice
            }else{
                res = make_pair(2,0);
            }
        }else{//case 3 and 5: 3. a1,b1->miss,b2->a3,b3; 5. a1+b1->miss,b2->a3+b3
            vector<size_t> both_parents; // to see if a1 and b1 both exist
            vector<size_t> both_kids; // to see if a3 and b3 both exist
            if(movieInfo.nodes[ajd_cell_idx[0]].parent_num > 0){
                assert(movieInfo.nodes[ajd_cell_idx[0]].parent_num == 1); // should be exactly one
                both_parents.push_back(movieInfo.nodes[ajd_cell_idx[0]].parents[0]);
            }
            if(missing_type != MISS_AT_TRACK_START){
                both_parents.push_back(parentKid_idx[0]);
            }
            if(movieInfo.nodes[ajd_cell_idx[0]].kid_num > 0){
                assert(movieInfo.nodes[ajd_cell_idx[0]].kid_num == 1); // should be exactly one
                both_kids.push_back(movieInfo.nodes[ajd_cell_idx[0]].kids[0]);
            }
            if(missing_type != MISS_AT_TRACK_END){
                both_kids.push_back(parentKid_idx[1]);
            }
            if(both_parents.size()==1 && both_kids.size()==1){
                // scenario 5: both_parent can jump to both_kid; merged region can link to both of them.
                size_t pk[] = {both_parents[0], both_kids[0]};
                float p2cur_Cost, cur2k_Cost;
                bool p2miss2k_valid = parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, pk,
                                                               missing_type, sz, p2cur_Cost, cur2k_Cost);
                // scenario 7: there is "another region" have not been considered.
                // If we add that region, both parents->merged region->both kids+the "another region". or
                // both parents+the "another region"->merged region->both kids
                bool p2miss2k_and_one_more_valid = false;
                if (p2miss2k_valid || multiParentsKidsValidLinkTest(cellSegment, both_parents[0], both_parents[1], new_cell_idx, new_cell_frame,
                                                                    p2cur_Cost, cur2k_Cost)){
                    res = make_pair(5,1);
                }else if(parentKidLink_valid){
                    res = make_pair(1,1);
                }else{
                    res = make_pair(5,0);
                }
            }else{
                vector<vector<size_t>> double_p_k_idx(2);
                double_p_k_idx[0] = both_parents;
                double_p_k_idx[1] = both_kids;
                bool p_valid = false, k_valid = false;
                if(both_parents.size() == 2){
                    p_valid = mergeSplitBothValid(cellSegment, true, both_parents, new_cell_idx,
                                                  new_cell_frame, extra_new_cells_loc_idx);
                    if(p_valid){
                        if(extra_new_cells_loc_idx[0].size() == 0 || extra_new_cells_loc_idx[1].size() == 0 ){
                            extra_new_cells_loc_idx.clear();
                            res = make_pair(1, 1);
                        }else{
                            res = make_pair(3, 1);
                        }
                    }
                }
                if(!p_valid && both_kids.size() == 2){
                    k_valid = mergeSplitBothValid(cellSegment, false, both_kids, new_cell_idx,
                                                  new_cell_frame, extra_new_cells_loc_idx);

                    if(k_valid){
                        if(extra_new_cells_loc_idx[0].size() == 0 || extra_new_cells_loc_idx[1].size() == 0 ){
                            extra_new_cells_loc_idx.clear();
                            res = make_pair(1, 1);
                        }else{
                            res = make_pair(3, 1);
                        }
                    }
                }
                if(res.first == 0){
                    if(parentKidLink_valid){
                        res = make_pair(1, 1);
                    }else{
                        res = make_pair(3, 0);
                    }
                }
            }
        }

    }
    return res;
}
//check if the newly detected cell indeed corresponds to two cells in previous or latter frame.
// we only check the case of two cells
bool cellTrackingMain::multiNeis_check(cellSegmentMain &cellSegment, size_t exist_parent_idx,
                                       size_t exist_kid_idx, vector<size_t> &new_cell_idx, vector<size_t> &new_p_k_pair){
    int f = movieInfo.frames[exist_parent_idx];
    vector<size_t> labels = extractValsGivenIdx_type<size_t>(&cellSegment.cell_label_maps[f], new_cell_idx, CV_32S);
    labels = vec_Add(labels, cumulative_cell_nums[f-1]-1);
    unordered_map<size_t, size_t> freq = frequecy_cnt(labels);
    if (freq.size() < 2) return false;
    size_t max_2nd, max_2nd_sz = 0;
    for (auto p : freq){
        if (p.first > 0 && p.first != exist_parent_idx && p.second > max_2nd_sz){
            max_2nd = p.first;
        }
    }
    if (movieInfo.nodes[max_2nd].kid_num!=1 ||
            movieInfo.frames[movieInfo.nodes[max_2nd].kids[0]]!=movieInfo.frames[exist_kid_idx]) return false;

    new_p_k_pair.resize(2);
    new_p_k_pair[0] = max_2nd;
    new_p_k_pair[1] = movieInfo.nodes[max_2nd].kids[0];
}

/**
 * @brief multiParentsKidsValidLinkTest: for scenario 7: a1->a2+missed->b3+c3 or a1+c1->a2+missed->b3.
 * @param new_cell_idx
 * @param new_cell_frame
 * @param node_idx
 * @param sz
 * @param cost
 */
bool cellTrackingMain::multiParentsKidsValidLinkTest(cellSegmentMain &cellSegment, size_t exist_parent_idx,
                                                     size_t exist_kid_idx, vector<size_t> &new_cell_idx, int new_cell_frame,
                                                     float cost1_lb, float cost2_lb){
    if(movieInfo.nodes[exist_parent_idx].kid_num + movieInfo.nodes[exist_kid_idx].parent_num != 1){
        // at most one of node has valid node in the new_cell_frame
        return false;
    }
    float cost1 = INFINITY;
    float cost2 = INFINITY;
    if(movieInfo.nodes[exist_parent_idx].kid_num == 1){
        size_t curr_node = movieInfo.nodes[exist_parent_idx].kids[0];
        assert(movieInfo.frames[curr_node] == new_cell_frame);
        if (movieInfo.nodes[exist_parent_idx].kid_cost[0] <= cost2_lb) return false;

        size_t best_id, best_ov_sz = 0;
        for(nodeRelation nr : movieInfo.nodes[curr_node].neighbors){
            if (nr.node_id != exist_kid_idx &&
                    movieInfo.frames[nr.node_id] == movieInfo.frames[exist_kid_idx]){
                if (nr.overlap_size > best_ov_sz){
                    best_ov_sz = nr.overlap_size;
                    best_id = nr.node_id;
                }
            }
        }
        if(best_ov_sz > 0){
            cost1 = 0;
            vector<size_t> node_idx = {best_id, exist_kid_idx};
            parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, node_idx,
                                     cellSegment.cell_label_maps[0].size, cost2);
        }
    }else if (movieInfo.nodes[exist_kid_idx].parent_num == 1){
        size_t curr_node = movieInfo.nodes[exist_kid_idx].parents[0];
        assert(movieInfo.frames[curr_node] == new_cell_frame);
        if (movieInfo.nodes[exist_kid_idx].parent_cost[0] <= cost2_lb) return false;
        size_t best_id, best_ov_sz = 0;
        for(nodeRelation nr : movieInfo.nodes[curr_node].preNeighbors){
            if (nr.node_id != exist_parent_idx &&
                    movieInfo.frames[nr.node_id] == movieInfo.frames[exist_parent_idx]){
                if (nr.overlap_size > best_ov_sz){
                    best_ov_sz = nr.overlap_size;
                    best_id = nr.node_id;
                }
            }
        }
        if(best_ov_sz > 0){
            cost2 = 0;
            vector<size_t> node_idx = {best_id, exist_parent_idx};
            parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, node_idx,
                                     cellSegment.cell_label_maps[0].size, cost1);
        }
    }else{
        qFatal("number of valid node is wrong");
    }
    if(MAX(cost1, cost2) < abs(p4tracking.observationCost)){
        return true;
    }else return false;
}
// Both parents can link to both kids and if merged, they also can link each other
bool cellTrackingMain::mergeSplitBothValid(cellSegmentMain &cellSegment, bool parent_flag, vector<size_t> givenSeedCells,
                                           vector<size_t> valid_loc_idx, int valid_loc_frame,
                                           vector<vector<size_t>> &splitted_reg_loc_idx){
    assert(givenSeedCells.size() == 2);
    if (movieInfo.frames[givenSeedCells[0]] != movieInfo.frames[givenSeedCells[1]]) return false;

    splitted_reg_loc_idx.resize(2);

    // merge test
    vector<vector<size_t>> reg4seeds(2);
    reg4seeds[0] = movieInfo.voxIdx[givenSeedCells[0]];
    reg4seeds[1] = movieInfo.voxIdx[givenSeedCells[1]];
    vector<size_t> both_seeds = reg4seeds[0];
    both_seeds.reserve(both_seeds.size() + reg4seeds[1].size());
    both_seeds.insert(both_seeds.end(), reg4seeds[1].begin(), reg4seeds[1].end());
    int reg4seeds_frame = movieInfo.frames[givenSeedCells[0]];
    float dummy_cost;
    if(parentOrKidValidLinkTest(valid_loc_idx, valid_loc_frame, givenSeedCells,
                                cellSegment.cell_label_maps[0].size, dummy_cost)){
        bool gapBased = true, usePriorGapMap = false;
        float dummy_cost[2];
        if(bisectValidTest(cellSegment, movieInfo.labelInMap.size(), valid_loc_idx, valid_loc_frame,
                           reg4seeds, reg4seeds_frame, gapBased, usePriorGapMap,
                            splitted_reg_loc_idx, dummy_cost) ||
                bisectValidTest(cellSegment, movieInfo.labelInMap.size(), valid_loc_idx, valid_loc_frame,
                                           reg4seeds, reg4seeds_frame, !gapBased, usePriorGapMap,
                                            splitted_reg_loc_idx, dummy_cost)){
            long long oppo_side_nei[2] = {-1, -1};
            if(parent_flag){
                if (movieInfo.nodes[givenSeedCells[0]].kid_num == 1){
                    oppo_side_nei[0] = movieInfo.nodes[givenSeedCells[0]].kids[0];
                }
                if (movieInfo.nodes[givenSeedCells[1]].kid_num == 1){
                    oppo_side_nei[1] = movieInfo.nodes[givenSeedCells[1]].kids[0];
                }
            }else{
                if (movieInfo.nodes[givenSeedCells[0]].parent_num == 1){
                    oppo_side_nei[0] = movieInfo.nodes[givenSeedCells[0]].parents[0];
                }
                if (movieInfo.nodes[givenSeedCells[1]].kid_num == 1){
                    oppo_side_nei[1] = movieInfo.nodes[givenSeedCells[1]].parents[0];
                }
            }
            float oppo_seed1_valid = true, oppo_seed2_valid = true;

            if(oppo_side_nei[0] == -1 && oppo_side_nei[1] == -1){
                return true;
            }else {
                float dummy_cost;
                if(oppo_side_nei[0] != -1){
                    oppo_seed1_valid = parentOrKidValidLinkTest(splitted_reg_loc_idx[0], valid_loc_frame, oppo_side_nei[0],
                            cellSegment.cell_label_maps[0].size, dummy_cost);
                }
                if(oppo_side_nei[1] != -1){
                    oppo_seed2_valid = parentOrKidValidLinkTest(splitted_reg_loc_idx[1], valid_loc_frame, oppo_side_nei[1],
                            cellSegment.cell_label_maps[0].size, dummy_cost);
                }
                if(oppo_seed1_valid && oppo_seed2_valid){
                    return true;
                }else{
                    return false;
                }
            }
        }else{
            return false;
        }
    }else{
        return false;
    }
}
/**
 * @brief checkSeedCoveredByExistingCell: when retrieve missing cell, the seed we get may already be covered by an
 * existing cell. It is worthy to check if the cell is real or an
 * over-merged one, because it clearly cannot be linked to the cell represented by the seed.
 * @param cellSegment
 * @param missing_type
 * @param parent_idx
 * @param kid_idx
 * @param missing_frame
 * @param in_seed_loc_idx
 * @param out_seed_loc_idx
 */
bool cellTrackingMain::checkSeedCoveredByExistingCell(cellSegmentMain &cellSegment, int missing_type,
                                    size_t parent_idx, size_t kid_idx, int missing_frame, size_t min_seed_sz,
                                    vector<size_t> &in_seed_loc_idx, vector<int> &seed_loc_existing_labels,
                                    unordered_set<size_t> &cell_idx_can_be_removed,
                                    vector<pair<size_t, int>> &seeds_missing_type,
                                    vector<vector<size_t>> &seeds_loc_idx){
    vector<size_t> out_seed_loc_idx (0);
    // conditon 1. the existing cell overlapped with seed is an orphan.
    vector<int>overlapped_labels = vec_atrange(seed_loc_existing_labels, (int)INFINITY, 0, true);
    unordered_map<int, size_t> freq_map = frequecy_cnt(overlapped_labels);
    unordered_set<int> labels_can_be_removed;
    //labels_can_be_removed.insert(0);
    for(auto it : freq_map){
        size_t curr_exist_idx = it.first;
        curr_exist_idx -= 1;
        curr_exist_idx += missing_frame == 0? 0:cumulative_cell_nums[missing_frame-1];
        if(movieInfo.nodes[curr_exist_idx].parent_num == 0 && movieInfo.nodes[curr_exist_idx].kid_num == 0){
            labels_can_be_removed.insert(it.first);
        }
    }
    if(labels_can_be_removed.size() > 0){
        for(size_t i=0; i<in_seed_loc_idx.size(); i++){
            int label = seed_loc_existing_labels[i];
            if (label==0 || labels_can_be_removed.find(label) != labels_can_be_removed.end()){
                out_seed_loc_idx.push_back(in_seed_loc_idx[i]);
            }
        }
        if(out_seed_loc_idx.size() >= min_seed_sz){
            if(missing_type == MISS_AT_TRACK_START){
                seeds_missing_type.push_back(make_pair(kid_idx, missing_type));
            }else{
                seeds_missing_type.push_back(make_pair(parent_idx, missing_type));
            }
            seeds_loc_idx.push_back(out_seed_loc_idx);
            size_t append_idx = missing_frame == 0? 0:cumulative_cell_nums[missing_frame-1];
            for(size_t it : labels_can_be_removed){
                cell_idx_can_be_removed.insert(it-1+append_idx);
            }
            return true;
        }
    }
    // condition 2, test if the largest existing overlapped cell is indeed a merged region.
    size_t max_overlapped_sz = 0;
    int max_overlapped_cell_label;
    for(unordered_map<int, size_t>::iterator it = freq_map.begin(); it != freq_map.end(); it ++){
        if(it->second > max_overlapped_sz){
            max_overlapped_cell_label = it->first;
            max_overlapped_sz = it->second;
        }
    }
    if(max_overlapped_sz + out_seed_loc_idx.size() < min_seed_sz){
        return false;
    }
    size_t candidate_cell_idx = max_overlapped_cell_label;
    candidate_cell_idx -= 1;
    candidate_cell_idx += missing_frame == 0? 0:cumulative_cell_nums[missing_frame-1];

//    vector<pair<size_t, int>> seeds_missing_type;
//    vector<vector<size_t>> seeds_loc_idx;
    bool seeds_found = false;
    if(missing_type != MISS_AT_TRACK_START){
        bool parent_flag = true;
        seeds_found = extractSeedInGivenCell(cellSegment, parent_flag, candidate_cell_idx, movieInfo.frames[parent_idx],
                               seeds_missing_type, seeds_loc_idx);
    }
    if(!seeds_found && missing_type != MISS_AT_TRACK_END){
        bool parent_flag = false;
        seeds_found = extractSeedInGivenCell(cellSegment, parent_flag, candidate_cell_idx, movieInfo.frames[kid_idx],
                               seeds_missing_type, seeds_loc_idx);
    }
    if(!seeds_found){
        return false;
    }
    cell_idx_can_be_removed.insert(candidate_cell_idx);
    return true;
}

/**
 * @brief extractSeedInGivenCell: extract all the seeds that contained in current cell
 * seedBase indicate which frame should we consider
 * @param cellSegment
 * @param missing_type
 * @param parent_idx
 * @param kid_idx
 * @param seed_loc_idx
 * @return
 */
bool cellTrackingMain::extractSeedInGivenCell(cellSegmentMain &cellSegment, bool parent_flag, size_t givenCell,
                                              int cell4seed_frame, vector<pair<size_t, int>> seeds_missing_type,
                                              vector<vector<size_t>> &seeds_loc_idx){
    int frame4seed = cell4seed_frame;
    int frame4givenCell = movieInfo.frames[givenCell];
    if(parent_flag){
        if(movieInfo.nodes[givenCell].parent_num > 0 &&
                frame4seed == movieInfo.frames[movieInfo.nodes[givenCell].parents[0]]){
            // If the given cell has a proper parent/kid in the frame for seed detection,
            // we should not remove the given cell. Thus no need to extract seed anymore.
            return false;
        }
    }else if(!parent_flag){
        if(movieInfo.nodes[givenCell].kid_num > 0 &&
                frame4seed == movieInfo.frames[movieInfo.nodes[givenCell].kids[0]]){
            // If the given cell has a proper parent/kid in the frame for seed detection,
            // we should not remove the given cell. Thus no need to extract seed anymore.
            return false;
        }
    }
    // Now the given cell has no valid linking in frame frame4seed, we can test if more than two seeds are there.
//    [possible_cells, neis_PorK_frame] = ...
//            extractParentsOrKidsGivenCell(...
//            movieInfo, cell_id, test_frame, frame4detect_cell_family, parent_flag);
    vector<size_t> givenCell_Neis;
    if(parent_flag){
        for(nodeRelation nr : movieInfo.nodes[givenCell].neighbors){
            if(movieInfo.frames[nr.node_id] == frame4seed){
                givenCell_Neis.push_back(nr.node_id);
            }
        }
    }else{
        for(nodeRelation nr : movieInfo.nodes[givenCell].preNeighbors){
            if(movieInfo.frames[nr.node_id] == frame4seed){
                givenCell_Neis.push_back(nr.node_id);
            }
        }
    }
    vector<size_t> candidate_included_cells(0);
    for(size_t nei : givenCell_Neis){
        size_t cur_nei_best_buddy;
        float cur_nei_best_buddy_cost;
        if(movieInfo.nodes[nei].nodeId2trackId>=0 &&
                findBestNeighbor(nei, cur_nei_best_buddy, cur_nei_best_buddy_cost, frame4givenCell) &&
                cur_nei_best_buddy == givenCell){
            // one more thing: nei should not have a valid parent/kid in frame4givenCell
            if(parent_flag){
                if(movieInfo.nodes[nei].kid_num == 0 ||
                        movieInfo.frames[movieInfo.nodes[nei].kids[0]] != frame4givenCell){
                    candidate_included_cells.push_back(nei);
                }
            }else{
                if(movieInfo.nodes[nei].parent_num == 0 ||
                        movieInfo.frames[movieInfo.nodes[nei].parents[0]] != frame4givenCell){
                    candidate_included_cells.push_back(nei);
                }
            }
        }
    }
    //if only 2 cell, then no need to process here; if more than 2 cell, split/merge test cannot handle
    if(candidate_included_cells.size() < 3){
        return false;
    }
    seeds_loc_idx.resize(0);
    seeds_missing_type.resize(0);
    int givenCell_label = movieInfo.labelInMap[givenCell];
    for(size_t candidate : candidate_included_cells){
        pair<size_t, int> curr_seeds_missing_type;
        if(parent_flag){
            if (movieInfo.nodes[candidate].kid_num > 0){
                curr_seeds_missing_type = make_pair(candidate, MISS_AS_JUMP);
            }else{
                curr_seeds_missing_type = make_pair(candidate, MISS_AT_TRACK_END);
            }
        }else{
            if (movieInfo.nodes[candidate].parent_num > 0){
                curr_seeds_missing_type = make_pair(movieInfo.nodes[candidate].parents[0], MISS_AS_JUMP);
            }else{
                curr_seeds_missing_type = make_pair(candidate, MISS_AT_TRACK_START);
            }
        }

        int missing_type = curr_seeds_missing_type.second;
        vector<int> dummry_missing_frames;
        vector<size_t> cur_seed_loc_idx;
        size_t parent_idx, kid_idx;
        if (missing_type == MISS_AS_JUMP){
            parent_idx = curr_seeds_missing_type.first;
            kid_idx = movieInfo.nodes[parent_idx].parents[0];
        }else if(missing_type == MISS_AT_TRACK_START){
            parent_idx = -1;
            kid_idx = curr_seeds_missing_type.first;
        }else{
            parent_idx = curr_seeds_missing_type.first;
            kid_idx = -1;
        }
        if(extractSeedFromGivenCell(cellSegment, missing_type, parent_idx,
                                    kid_idx, dummry_missing_frames, cur_seed_loc_idx)){
            vector<size_t> cur_seed_loc_idx_valid(0);
            for(size_t one_seed_idx : cur_seed_loc_idx){
                if(cellSegment.cell_label_maps[frame4givenCell].at<int>(one_seed_idx) == 0 ||
                        cellSegment.cell_label_maps[frame4givenCell].at<int>(one_seed_idx) == givenCell_label){
                    cur_seed_loc_idx_valid.push_back(one_seed_idx);
                }
            }
            if(cur_seed_loc_idx_valid.size() *2 >= MAX(cur_seed_loc_idx.size(), cellSegment.p4segVol.min_seed_size*2)){
                seeds_loc_idx.push_back(cur_seed_loc_idx_valid);
                seeds_missing_type.push_back(curr_seeds_missing_type);
            }
        }
    }
    if(seeds_loc_idx.size() < 3){
        seeds_loc_idx.clear();
        seeds_missing_type.clear();
        return false;
    }
    return true;
}
/**
 * @brief extractSeedFromGivenCell: given a cell, extract the seed region (for its parent/kid cells)
 * in adjacent frames based on it.
 * @param cellSegment
 * @param missing_type
 * @param parent_idx
 * @param kid_idx
 * @param missing_frames
 * @param seed_loc_idx
 * @return
 */
bool cellTrackingMain::extractSeedFromGivenCell(cellSegmentMain &cellSegment, int missing_type,
                                                 size_t parent_idx, size_t kid_idx, vector<int> &missing_frames,
                                                 vector<size_t> &seed_loc_idx){
    long sz_single_frame = cellSegment.data_rows_cols_slices[0]*
            cellSegment.data_rows_cols_slices[1]*cellSegment.data_rows_cols_slices[2];
    missing_frames.resize(0); // there can be two jumped frames
    unsigned char threshold;
    seed_loc_idx.resize(0);
    if (missing_type == MISS_AS_JUMP){
        for(int f = movieInfo.frames[parent_idx] + 1; f < movieInfo.frames[kid_idx]; f++){
            missing_frames.push_back(f);
        }
        seed_loc_idx = intersection(movieInfo.voxIdx[parent_idx], movieInfo.voxIdx[kid_idx]);
    }else if(missing_type == MISS_AT_TRACK_START){
        if(movieInfo.frames[kid_idx] == 0){
            return false;
        }
        missing_frames = {movieInfo.frames[kid_idx] - 1};
        threshold =
                cellSegment.threshold_maps[movieInfo.frames[kid_idx]].at<unsigned char>(movieInfo.voxIdx[kid_idx][0]);
        if(threshold <= 0){
            qDebug("check point");
        }
        assert(threshold > 0);
        unsigned char *ind = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*missing_frames[0]; // sub-matrix pointer
        Mat single_frame(3, cellSegment.normalized_data4d.size, CV_8U, ind);
        for(size_t i : movieInfo.voxIdx[kid_idx]){
            if(single_frame.at<unsigned char>(i) > threshold){
                seed_loc_idx.push_back(i);
            }
        }
    }else if(missing_type == MISS_AT_TRACK_END){
        if(movieInfo.frames[parent_idx] >= cellSegment.threshold_maps.size()-1){
            return false;
        }
        missing_frames = {movieInfo.frames[parent_idx] + 1};
        threshold =
                cellSegment.threshold_maps[movieInfo.frames[parent_idx]].at<unsigned char>(movieInfo.voxIdx[parent_idx][0]);
        if(threshold <= 0){
            qDebug("check point");
        }
        assert(threshold > 0);
        unsigned char *ind2 = (unsigned char*)cellSegment.normalized_data4d.data + sz_single_frame*missing_frames[0]; // sub-matrix pointer
        Mat single_frame2(3, cellSegment.normalized_data4d.size, CV_8U, ind2);
        for(size_t i : movieInfo.voxIdx[parent_idx]){
            if(single_frame2.at<unsigned char>(i) > threshold){
                seed_loc_idx.push_back(i);
            }
        }
    }else{
        qFatal("unKNOWN missing_type");
    }
    if (seed_loc_idx.size() < cellSegment.p4segVol.min_seed_size){
        seed_loc_idx.clear();// equals to seed_loc_idx.resize(0)
        return false;
    }else{
        // pick the largest one
        Mat1b tightBw;
        int start_yxz[3];
        idx2tightBwMap(seed_loc_idx, cellSegment.normalized_data4d.size, tightBw, start_yxz);
        Mat1i label_map;
        int n = connectedComponents3d(&tightBw, label_map, 6);
        if (n > 1){
            int kept = largestRegionIdExtract(&label_map, n);
            vector<size_t> kept_idx;
            extractVoxIdxGivenId(&label_map, kept_idx, kept);
            seed_loc_idx.clear();
            start_yxz[0] = -start_yxz[0];
            start_yxz[1] = -start_yxz[1];
            start_yxz[2] = -start_yxz[2];
            coordinateTransfer(kept_idx, tightBw.size,
                        seed_loc_idx, start_yxz, cellSegment.normalized_data4d.size);
        }
        if(seed_loc_idx.size() < cellSegment.p4segVol.min_seed_size){
            seed_loc_idx.clear();
            return false;
        }else{
            return true;
        }
    }
}

/**
 * @brief parentOrKidValidLinkTest: given a newly detected cell, see if it could be linked to existing cells
 * @param new_cell_idx
 * @param new_cell_frame
 * @param parent_idx
 * @param kid_idx
 * @return
 */
bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                size_t parentKid_idx[2], int missing_type, MatSize sz){
    float dummy_cost;
    if(missing_type == MISS_AS_JUMP){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, dummy_cost) &&
                parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, dummy_cost);
    }else if(missing_type == MISS_AT_TRACK_START){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, dummy_cost);
    }else{
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz,dummy_cost);
    }
}

bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                size_t parentKid_idx[2], int missing_type, MatSize sz,
                                                float &cost1, float &cost2){
    cost1 = INFINITY;
    cost2 = INFINITY;
    if(missing_type == MISS_AS_JUMP){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, cost1) &&
                parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, cost2);
    }else if(missing_type == MISS_AT_TRACK_START){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, cost1);
    }else{
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, cost2);
    }
}

bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                size_t node_idx, MatSize sz, float &cost){
    float dummy_c2n, dummy_n2c;
    float dist = voxelwise_avg_distance(new_cell_idx, new_cell_frame,
                                        movieInfo.voxIdx[node_idx], movieInfo.frames[node_idx], sz, dummy_c2n, dummy_c2n);
    cost = distance2cost(dist);
    return (cost<abs(p4tracking.observationCost));
}


// if input parents or kids are more than one
bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                vector<vector<size_t>> parentKid_idx, int missing_type, MatSize sz){
    float dummy_cost;
    if(missing_type == MISS_AS_JUMP){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, dummy_cost) &&
                parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, dummy_cost);
    }else if(missing_type == MISS_AT_TRACK_START){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, dummy_cost);
    }else{
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz,dummy_cost);
    }
}

bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                vector<vector<size_t>> parentKid_idx, int missing_type, MatSize sz,
                                                float &cost1, float &cost2){
    cost1 = INFINITY;
    cost2 = INFINITY;
    if(missing_type == MISS_AS_JUMP){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, cost1) &&
                parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, cost2);
    }else if(missing_type == MISS_AT_TRACK_START){
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[1], sz, cost1);
    }else{
        return parentOrKidValidLinkTest(new_cell_idx, new_cell_frame, parentKid_idx[0], sz, cost2);
    }
}

bool cellTrackingMain::parentOrKidValidLinkTest(vector<size_t> &new_cell_idx, int new_cell_frame,
                                                vector<size_t> node_idx, MatSize sz, float &cost){
    float dummy_c2n, dummy_n2c;
    vector<size_t> node_loc_idx;
    FOREACH_i(node_idx){
        node_loc_idx.reserve(node_loc_idx.size() + movieInfo.voxIdx[node_idx[i]].size());
        node_loc_idx.insert(node_loc_idx.end(), movieInfo.voxIdx[node_idx[i]].begin(), movieInfo.voxIdx[node_idx[i]].end());
    }
    float dist = voxelwise_avg_distance(new_cell_idx, new_cell_frame,
                                        node_loc_idx, movieInfo.frames[node_idx[0]], sz, dummy_c2n, dummy_c2n);
    cost = distance2cost(dist);
    return (cost<abs(p4tracking.observationCost));
}

///--------------------------------------combine results from batch processing--------------------------//
/// The key principal of cell fusion is to determine if two cells from overlapped regions are the same or not.
/// We used a IoU=0.5 as principal, if ov > IoU, they are the same. Otherwise, we keep the one that further to
/// the spatial or temporal boundary.
cellTrackingMain::cellTrackingMain(const vector<int> &fixed_crop_sz, const vector<int> &overlap_sz, const QString &dataFolderName, const QString &resFolderName){
    //vector<int> fixed_crop_sz = {493, 366, 259};
    //vector<int> overlap_sz = {50, 50, 24, 5}; // y, x, z, t
    init_parameter();
    batchResultsFusion(dataFolderName, resFolderName, fixed_crop_sz, overlap_sz);
    tracking_sucess =false;
    bool get_res_from_txt = true, get_jumpCost_only = false;
    mccTracker_one2one(get_jumpCost_only, get_res_from_txt);
    tracking_sucess = true;
}
bool cellTrackingMain::batchResultsFusion(const QString &dataFolderName, const QString &resFolderName,
                                          const vector<int> &fixed_crop_sz, const vector<int> &overlap_sz){
    overlapped_frames.reserve(p4tracking.k);
    // step 1. fuse data in one batch; get the mapping function: batch_id, frame, label ==> newIdx
    QDir root_directory(dataFolderName);
    QStringList data_batches = root_directory.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    for(int batch_id=0; batch_id<data_batches.count(); batch_id++){
        oneBatchResultsFusion(batch_id, dataFolderName+"/"+data_batches[batch_id], fixed_crop_sz, overlap_sz);
    }
    // step 2. load all the cell information (movieInfo.txt), build movieInfo
    movieInfo.nodes.resize(fuse_batch_processed_cell_cnt);
    movieInfo.frames.resize(fuse_batch_processed_cell_cnt);
    movieInfo.xCoord.resize(fuse_batch_processed_cell_cnt);
    movieInfo.yCoord.resize(fuse_batch_processed_cell_cnt);
    movieInfo.zCoord.resize(fuse_batch_processed_cell_cnt);
    FOREACH_i(movieInfo.nodes){
        movieInfo.nodes[i].node_id = i;
    }
    for(auto &nodeId_yxz : fusedFrameLabel2centerYXZ){
        movieInfo.xCoord[nodeId_yxz.first] = nodeId_yxz.second[1];
        movieInfo.yCoord[nodeId_yxz.first] = nodeId_yxz.second[0];
        movieInfo.zCoord[nodeId_yxz.first] = nodeId_yxz.second[2];
    }
    movieInfo.overall_neighbor_num = 0;
    for(int batch_id=0; batch_id<data_batches.count(); batch_id++){
        oneBatchMovieInfoParse(batch_id, dataFolderName+"/"+data_batches[batch_id]);
    }
}

void cellTrackingMain::oneBatchMovieInfoParse(int batch_id, const QString &subfolderName){
    vector<QString> crop_names = {"frontleft", "frontright", "backleft", "back_right"};
    QRegExp rx("(\\d+)");
    if(rx.indexIn(subfolderName) == -1) qFatal("Wrong file name");
    int start_frame = rx.cap(1).toInt();
    for(int i=0; i<crop_names.size(); i++){
        QString txt_file_name = subfolderName + "/" + crop_names[i] + "/" + "movieInfo.txt";
        QFile txt_file(txt_file_name);
        if (!txt_file.open(QIODevice::ReadOnly)){
            QMessageBox::information(0, "error", txt_file.errorString());
        }
        QTextStream in(&txt_file);
        long cell_num;
        sscanf(in.readLine().toStdString().c_str(), "p cell num: %ld", &cell_num);

        for(int i=0; i<2; i++){
            // skip the first 3 lines
            //            p cell num: 65164
            //            p gamma: 1.50283 8.1355
            //            p jumpRatio: 0.691873 0.162243 0.145884
            in.readLine();
        }
        double en_cost, ex_cost, boz_cost;
        sscanf(in.readLine().toStdString().c_str(), "p enexObz: %lf %lf %lf", &en_cost, &ex_cost, &boz_cost);
        unordered_map<int, size_t> cell_localId2Key;
        while(!in.atEnd()) {
            string line = in.readLine().toStdString();
            // start parsing each line
            switch (line[0]) {
            case 'c':                  /* skip lines with comments */
            case '\n':                 /* skip empty lines   */
            case 'n':{
                int id, frame, labelInMap;
                sscanf(line.c_str(), "%*c %d %d %d", &id, &frame, &labelInMap);
                if(labelInMap != 0){
                    cell_localId2Key[id] = key(batch_id, start_frame+frame, i, labelInMap);
                    auto it = oldinfo2newIdx.find(cell_localId2Key[id]);
                    if(it != oldinfo2newIdx.end()){
                        movieInfo.frames[it->second] = frame;
                        movieInfo.nodes[it->second].in_cost = en_cost;
                        movieInfo.nodes[it->second].out_cost = ex_cost;
                    }
                }
            }
            case '\0':                 /* skip empty lines at the end of file */
                break;
            case 'p':
            case 'a':{
                int tail = 0;
                int head = 0;
                double distance, link_cost;
                sscanf(line.c_str(), "%*c %d %d %lf %lf", &tail, &head, &distance, &link_cost);
                size_t real_tail, real_tailKey, real_head, real_headKey;
                real_tailKey = cell_localId2Key[tail];
                real_headKey = cell_localId2Key[head];
                auto it_tail = oldinfo2newIdx.find(real_tailKey);
                auto it_head = oldinfo2newIdx.find(real_headKey);
                if(it_tail!=oldinfo2newIdx.end() && it_head!=oldinfo2newIdx.end()){
                    real_tail = it_tail->second;
                    real_head = it_head->second;
                    bool append_flag = true;
                    for(auto &nn : movieInfo.nodes[real_tail].neighbors){
                        if(nn.node_id == real_head){
                            append_flag = false;
                            break;
                        }
                    }
                    if(append_flag){
                        nodeRelation tmp;
                        tmp.node_id = real_head;
                        tmp.dist_c2n = distance;
                        tmp.link_cost = link_cost;
                        movieInfo.nodes[real_tail].neighbors.push_back(tmp);
                        movieInfo.overall_neighbor_num ++;
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }
}

/**
 * @brief spaceFusion: fuse data in vertical and horizontal directions
 * @param subfolderName
 */
void cellTrackingMain::oneBatchResultsFusion(int batch_id, const QString &subfolderName, const vector<int> &fixed_crop_sz, const vector<int> &overlap_sz){
    vector<QString> crop_names = {"frontleft", "frontright", "backleft", "back_right"};
//    vector<int> fixed_crop_sz = {493, 366, 259};
//    vector<int> overlap_sz = {50, 50, 24};
    QDir root_directory(subfolderName+"/"+crop_names[0]);
    QStringList images = root_directory.entryList(QStringList() << "*.bin" ,QDir::Files);// <<"*.tif" << "*.JPG" if we want two type images

    //// Get the frame that temporally overlapped and use them to update overlapped_frames
    /// The principal is get the last p4tracking.k frames
    QRegExp rx("(\\d+)");
    if(rx.indexIn(images[images.count()-1]) == -1) qFatal("Wrong file name");
    int frame2save = 0;
    for (int i=0; i<images.count(); i++){
        rx.indexIn(images[i]);
        frame2save = MAX(frame2save, rx.cap(1).toInt());
    }
    vector<int> data_process_order (images.count());
    for(int i=0; i<images.count(); i++){
        rx.indexIn(images[i]);
        //qDebug("%d", rx.cap(1).toInt()-1+(images.count()-frame2save));
        data_process_order[rx.cap(1).toInt()-1+(images.count()-frame2save)] = i;
    }
    frame2save -= (p4tracking.k - 1);
    //if(overlapped_frames.size() == 0)
    int frame_saved_cnt = 0;
    for(int ii = 0; ii<data_process_order.size(); ii++) { //QT's version of for_each
        QString filename = images[data_process_order[ii]];
        QRegExp rx("(\\d+)");
        if(rx.indexIn(filename) == -1) qFatal("Wrong file name");
        int frame = rx.cap(1).toInt();
        if(overlapped_frames.size()>0 && frame < overlapped_frames.begin()->first){
            //// if temperally processed && out of the scope of jump
            continue;
        }
        vector<Mat1i> mat_crops(4); //mat_fl, mat_fr, mat_bl, mat_br;
        for(int i=0; i<crop_names.size(); i++){
            QString label_file_name = subfolderName + "/" + crop_names[i] + "/" + filename;
            QFile label_file(label_file_name);
            if (!label_file.open(QIODevice::ReadOnly)){
                qFatal("lost files!");
            }
            //QByteArray tmp = label_file.readAll();
            // tmp is local variable, which will be released soon, so we need copyTo
            Mat(3, fixed_crop_sz.data(), CV_32S, label_file.readAll().data()).copyTo(mat_crops[i]);
            label_file.close();
        }
        // do data fusion at spatial domain
        Mat lr_fuse_mat_u, lr_fuse_mat_d, ud_fuse_mat;
        vector<vector<int>> u_label_map_lr, d_label_map_lr, label_map_ud;
        spaceFusion_leftRight(mat_crops[0], mat_crops[1], lr_fuse_mat_u, overlap_sz[1], u_label_map_lr);
        //ccShowSliceLabelMat(lr_fuse_mat_u);
        spaceFusion_leftRight(mat_crops[2], mat_crops[3], lr_fuse_mat_d, overlap_sz[1], d_label_map_lr);
        spaceFusion_upDown(lr_fuse_mat_u, lr_fuse_mat_d, ud_fuse_mat, overlap_sz[0], label_map_ud);
        //ccShowSliceLabelMat(ud_fuse_mat);
        //// if temperally never processed
        if(overlapped_frames.size()==0 || frame > overlapped_frames.rbegin()->first){
            int max_label = 0;
            size_t new_idx;
            FOREACH_i(u_label_map_lr){
                FOREACH_j(u_label_map_lr[i]){
                    if(u_label_map_lr[i][j]==0 || label_map_ud[0][u_label_map_lr[i][j]]==0){
                        continue;
                    }
                    new_idx = fuse_batch_processed_cell_cnt + label_map_ud[0][u_label_map_lr[i][j]] - 1;
                    oldinfo2newIdx[key(batch_id, frame, i, j)] = new_idx;
                    newIdx2newinfo[new_idx] = newkey(frame, label_map_ud[0][u_label_map_lr[i][j]]);
                    newinfo2newIdx[newkey(frame, label_map_ud[0][u_label_map_lr[i][j]])] = new_idx;
                    max_label = MAX(max_label, label_map_ud[0][u_label_map_lr[i][j]]);
                }
            }
            FOREACH_i(d_label_map_lr){
                FOREACH_j(d_label_map_lr[i]){
                    if(d_label_map_lr[i][j]==0 || label_map_ud[1][d_label_map_lr[i][j]]==0){
                        continue;
                    }
                    new_idx = fuse_batch_processed_cell_cnt + label_map_ud[1][d_label_map_lr[i][j]] - 1;
                    oldinfo2newIdx[key(batch_id, frame, i+2, j)] = new_idx;
                    newIdx2newinfo[new_idx] = newkey(frame, label_map_ud[1][d_label_map_lr[i][j]]);
                    newinfo2newIdx[newkey(frame, label_map_ud[1][d_label_map_lr[i][j]])] = new_idx;
                    max_label = MAX(max_label, label_map_ud[1][d_label_map_lr[i][j]]);
                }
            }
            fuse_batch_processed_cell_cnt += max_label; // count start from 1 and idx starts from 0;

            //// update the fusedFrameLabel2centerYXZ
            vector<vector<size_t>> cell_voxIdx;
            extractVoxIdxList(&ud_fuse_mat, cell_voxIdx, max_label);
            FOREACH_j(cell_voxIdx){
                if(cell_voxIdx[j].size() == 0) continue;
                vector<int> y, x, z;
                vec_ind2sub(cell_voxIdx[j], y, x, z, ud_fuse_mat.size);
                fusedFrameLabel2centerYXZ[newinfo2newIdx[newkey(frame, j+1)]] = {vec_mean(y),
                        vec_mean(x), vec_mean(z)};
            }
            //// save the 16bit unsigned results as label map
            //Mat ud_fuse_mat_u16;
            //ud_fuse_mat.convertTo(ud_fuse_mat_u16, CV_16U);
            QString full_im_name = subfolderName + "/" + filename.left(12)+"fused.bin";
            //imwrite(full_im_name.toStdString(), ud_fuse_mat_u16);
            QFileInfo check_file(full_im_name);
            if(!(check_file.exists() && check_file.isFile())){
                ofstream fused_file(full_im_name.toStdString(), ios::binary);
                if (fused_file.is_open()){
                    fused_file.write((const char*)(ud_fuse_mat.data), ud_fuse_mat.elemSize() * ud_fuse_mat.total());
                    fused_file.close();
                }
            }
        }else{ //// if temperally processed,
            int ov_label_map_id = frame - overlapped_frames.begin()->first;
            temporalFusion(overlapped_frames[ov_label_map_id].second, ud_fuse_mat, batch_id, frame,
                           u_label_map_lr, d_label_map_lr, label_map_ud);
        }
        //// start to save the overlapped results for next round
        if(frame == frame2save){
            if(overlapped_frames.size() == p4tracking.k){
                overlapped_frames[frame_saved_cnt++] = make_pair(frame, ud_fuse_mat);
            }else{
                overlapped_frames.emplace_back(make_pair(frame, ud_fuse_mat));
            }
            frame2save++;
        }
    }
}

/**
 * @brief spaceFusion_leftRight
 * @param left
 * @param right
 * @param fusedMat
 * @param ov_sz : on the x-axis
 * @param oldLabel2newLabel
 */
void cellTrackingMain::spaceFusion_leftRight(Mat &left, Mat &right, Mat &fusedMat, int ov_sz, vector<vector<int>> &oldLabel2newLabel){
    double l_max_id, r_max_id;
    minMaxIdx(left, nullptr, &l_max_id);
    minMaxIdx(right, nullptr, &r_max_id);
    vector<vector<size_t>> l_cell_voxIdx, r_cell_voxIdx;
    extractVoxIdxList(&left, l_cell_voxIdx, (int)l_max_id);
    extractVoxIdxList(&right, r_cell_voxIdx, (int)r_max_id);

    oldLabel2newLabel.resize(2);
    oldLabel2newLabel[0].resize(l_max_id+1);
    FOREACH_i(oldLabel2newLabel[0]){ // use the left one as baseline
        oldLabel2newLabel[0][i] = i;
    }
    oldLabel2newLabel[1].resize(r_max_id+1);
    fill(oldLabel2newLabel[1].begin(), oldLabel2newLabel[1].end(), 0);
    // check the right one
    FOREACH_i(r_cell_voxIdx){
        vector<int> y, x, z;
        if(r_cell_voxIdx[i].size()==0) continue;
        vec_ind2sub(r_cell_voxIdx[i], y, x, z, right.size);
        if (vec_min(x) >= ov_sz){ // located in non-overlapping area
            oldLabel2newLabel[1][i+1] = i + 1 + l_max_id;
        }else{
            vector<size_t> to_left_idx;
            to_left_idx.reserve(x.size());
            size_t tmp;
            for(int j=0; j<x.size(); j++){
                if(x[j] < ov_sz){
                    vol_sub2ind(tmp, y[j], x[j] + left.size[1]-ov_sz, z[j], left.size);
                    to_left_idx.emplace_back(tmp);
                }
            }
            vector<float> left_vals = extractValsGivenIdx(&left, to_left_idx, CV_32S);
            unordered_map<float, size_t> freqs = frequecy_cnt(left_vals);
            int best_id = 0, best_ov_sz = 0;
            for(auto ele:freqs){
                if(ele.second>best_ov_sz){
                    best_id = (int) ele.first;
                    best_ov_sz = ele.second;
                }
            }
            //first test if there is a cell in left with IoU>0.5
            if(best_id != 0 && best_ov_sz*2>=(r_cell_voxIdx[i].size() + l_cell_voxIdx[best_id-1].size())){
                oldLabel2newLabel[1][i+1] = oldLabel2newLabel[0][best_id];
            }else{ // if there is no cells in left with IOU>0.5
                if(vec_mean(x) > (left.size[1]-ov_sz/2)){
                    //case 1: the cell is the right half panel=> set cells in left half panel as 0
                    oldLabel2newLabel[1][i+1] = i + 1 + l_max_id;
                    for(auto ele:freqs){
                        oldLabel2newLabel[0][ele.first] = 0;
                    }
                }else{
                    //case 2: the cell is the left half panel=> set this cell as 0
                    oldLabel2newLabel[1][i+1] = 0;
                }
            }
        }
    }
    int fuse_sz[3] = {left.size[0], left.size[1], left.size[2]};
    fuse_sz[1] += (right.size[1] - ov_sz);
    fusedMat = Mat::zeros(3, fuse_sz, CV_32S);
    unordered_set<int> used;
    for(int i=0; i<l_cell_voxIdx.size(); i++){
        if (oldLabel2newLabel[0][i+1] == 0){ // removed regions
            continue;
        }
        vector<int> y, x, z;
        vec_ind2sub(l_cell_voxIdx[i], y, x, z, left.size);
        vector<size_t> idx_fuse_mat;
        vec_sub2ind(idx_fuse_mat, y, x, z, fusedMat.size);
        used.insert(oldLabel2newLabel[0][i+1]);
        setValMat(fusedMat, CV_32S, idx_fuse_mat, (float)oldLabel2newLabel[0][i+1]);
    }

    for(int i=0; i<r_cell_voxIdx.size(); i++){
        if (oldLabel2newLabel[1][i+1] == 0 ||
                used.find(oldLabel2newLabel[1][i+1]) != used.end() ){ // removed regions
            continue;
        }
        vector<int> y, x, z;
        vec_ind2sub(r_cell_voxIdx[i], y, x, z, right.size);
        vector<size_t> idx_fuse_mat;
        vec_sub2ind(idx_fuse_mat, y, vec_Add(x, left.size[1]-ov_sz), z, fusedMat.size);
        setValMat(fusedMat, CV_32S, idx_fuse_mat, oldLabel2newLabel[1][i+1]);
    }
}
/**
 * @brief spaceFusion_upDown: the last step of spatial fusion
 * @param up
 * @param down
 * @param fusedMat
 * @param ov_sz: along y-axis
 * @param oldLabel2newLabel
 */
void cellTrackingMain::spaceFusion_upDown(Mat &up, Mat &down, Mat &fusedMat, int ov_sz, vector<vector<int>> &oldLabel2newLabel){
    double u_max_id, d_max_id;
    minMaxIdx(up, nullptr, &u_max_id);
    minMaxIdx(down, nullptr, &d_max_id);
    vector<vector<size_t>> u_cell_voxIdx, d_cell_voxIdx;
    extractVoxIdxList(&up, u_cell_voxIdx, (int)u_max_id);
    extractVoxIdxList(&down, d_cell_voxIdx, (int)d_max_id);

    oldLabel2newLabel.resize(2);
    oldLabel2newLabel[0].resize(u_max_id+1);
    FOREACH_i(oldLabel2newLabel[0]){ // use the left one as baseline
        oldLabel2newLabel[0][i] = i;
    }
    oldLabel2newLabel[1].resize(d_max_id+1);
    fill(oldLabel2newLabel[1].begin(), oldLabel2newLabel[1].end(), 0);
    // check the right one
    FOREACH_i(d_cell_voxIdx){
        vector<int> y, x, z;
        if(d_cell_voxIdx[i].size()==0) continue;
        vec_ind2sub(d_cell_voxIdx[i], y, x, z, down.size);
        if (vec_min(y) >= ov_sz){ // located in non-overlapping area
            oldLabel2newLabel[1][i+1] = i + 1 + u_max_id;
        }else{
            vector<size_t> to_up_idx;
            to_up_idx.reserve(x.size());
            size_t tmp;
            for(int j=0; j<y.size(); j++){
                if(y[j] < ov_sz){
                    vol_sub2ind(tmp, y[j]+up.size[0]-ov_sz, x[j], z[j], up.size);
                    to_up_idx.emplace_back(tmp);
                }
            }
            vector<float> up_vals = extractValsGivenIdx(&up, to_up_idx, CV_32S);
            unordered_map<float, size_t> freqs = frequecy_cnt(up_vals);
            int best_id = 0, best_ov_sz = 0;
            for(auto ele:freqs){
                if(ele.second>best_ov_sz){
                    best_id = (int) ele.first;
                    best_ov_sz = ele.first;
                }
            }
            //first test if there is a cell in left with IoU>0.5
            if(best_id != 0 && best_ov_sz*2>=(d_cell_voxIdx[i].size() + u_cell_voxIdx[best_id-1].size())){
                oldLabel2newLabel[1][i+1] = oldLabel2newLabel[0][best_id];
            }else{ // if there is no cells in left with IOU>0.5
                if(vec_mean(y) > (up.size[0]-ov_sz/2)){
                    //case 1: the cell is the right half panel=> set cells in left half panel as 0
                    oldLabel2newLabel[1][i+1] = i + 1 + u_max_id;
                    for(auto ele:freqs){
                        oldLabel2newLabel[0][ele.first] = 0;
                    }
                }else{
                    //case 2: the cell is the left half panel=> set this cell as 0
                    oldLabel2newLabel[1][i+1] = 0;
                }
            }
        }
    }
    int fuse_sz[3] = {up.size[0], up.size[1], up.size[2]};
    fuse_sz[0] += (down.size[0] - ov_sz);
    fusedMat = Mat::zeros(3, fuse_sz, CV_32S);
    unordered_set<int> used;
    for(int i=0; i<u_cell_voxIdx.size(); i++){
        if (oldLabel2newLabel[0][i+1] == 0){ // removed regions
            continue;
        }
        vector<int> y, x, z;
        vec_ind2sub(u_cell_voxIdx[i], y, x, z, up.size);
        vector<size_t> idx_fuse_mat;
        vec_sub2ind(idx_fuse_mat, y, x, z, fusedMat.size);
        used.insert(oldLabel2newLabel[0][i+1]);
        setValMat(fusedMat, CV_32S, idx_fuse_mat, (float)oldLabel2newLabel[0][i+1]);
    }

    for(int i=0; i<d_cell_voxIdx.size(); i++){
        if (oldLabel2newLabel[1][i+1] == 0 ||
                used.find(oldLabel2newLabel[1][i+1]) != used.end() ){ // removed regions
            continue;
        }
        vector<int> y, x, z;
        vec_ind2sub(d_cell_voxIdx[i], y, x, z, down.size);
        vector<size_t> idx_fuse_mat;
        vec_sub2ind(idx_fuse_mat, vec_Add(y, up.size[0]-ov_sz), x, z, fusedMat.size);
        setValMat(fusedMat, CV_32S, idx_fuse_mat, oldLabel2newLabel[1][i+1]);
    }
}
/**
 * @brief temporalFusion: merge two overlapped frames
 * @param kept: as reference; cells in this frame are kept
 * @param mov: the merged label maps (fl, fr, bl, br)
 * @param mov_batch_id
 * @param frame
 * @param u_label_map_lr
 * @param d_label_map_lr
 * @param label_map_ud
 */
void cellTrackingMain::temporalFusion(Mat &kept, Mat &mov, int mov_batch_id, int frame,
                                      vector<vector<int>> &u_label_map_lr, vector<vector<int>> &d_label_map_lr,
                                      vector<vector<int>> &label_map_ud){
    double kept_max_id, mov_max_id;
    minMaxIdx(kept, nullptr, &kept_max_id);
    minMaxIdx(mov, nullptr, &mov_max_id);
    unordered_map<int, vector<vector<int>>> mov_label2crop;
    FOREACH_i(u_label_map_lr){
        FOREACH_j(u_label_map_lr[i]){
            if(u_label_map_lr[i][j]==0 || label_map_ud[0][u_label_map_lr[i][j]]==0){
                continue;
            }
            mov_label2crop[label_map_ud[0][u_label_map_lr[i][j]]].push_back({mov_batch_id, frame, (int)i, (int)j});
        }
    }
    FOREACH_i(d_label_map_lr){
        FOREACH_j(d_label_map_lr[i]){
            if(d_label_map_lr[i][j]==0 || label_map_ud[1][d_label_map_lr[i][j]]==0){
                continue;
            }
            if(mov_label2crop.find(label_map_ud[1][d_label_map_lr[i][j]]) == mov_label2crop.end()){
                mov_label2crop[label_map_ud[1][d_label_map_lr[i][j]]].push_back({mov_batch_id, frame, (int)i+2, (int)j});
            }
        }
    }
    vector<vector<size_t>> kept_cell_voxIdx, mov_cell_voxIdx;
    extractVoxIdxList(&kept, kept_cell_voxIdx, (int)kept_max_id);
    extractVoxIdxList(&mov, mov_cell_voxIdx, (int)mov_max_id);

    // keep those cells in mov that IoU>0.5 with cells in kept
    FOREACH_i(mov_cell_voxIdx){
        int labelInMap = i+1;
        if(mov_label2crop.find(labelInMap) == mov_label2crop.end()){
            continue;
        }
        vector<float> kept_vals = extractValsGivenIdx(&kept, mov_cell_voxIdx[i], CV_32S);
        unordered_map<float, size_t> freqs = frequecy_cnt(kept_vals);
        int best_id = 0, best_ov_sz = 0;
        for(auto ele:freqs){
            if(ele.second>best_ov_sz){
                best_id = (int) ele.first;
                best_ov_sz = ele.first;
            }
        }
        //first test if there is a cell in left with IoU>0.5
        if(best_id != 0 && best_ov_sz*2>=(mov_cell_voxIdx[i].size() + kept_cell_voxIdx[best_id-1].size())){
            size_t kept_key = newkey(frame, best_id);
            if(newinfo2newIdx.find(kept_key) == newinfo2newIdx.end()){
                qFatal("refer to a non-exist cell in reference frame");
            }
            for(auto &key_eles : mov_label2crop[labelInMap]){
                size_t mov_key = key(key_eles[0], key_eles[1], key_eles[2], key_eles[3]);
                oldinfo2newIdx[mov_key] = newinfo2newIdx[kept_key];
            }
        }

    }
}
