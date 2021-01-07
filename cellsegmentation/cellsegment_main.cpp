#include "cellsegment_main.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>

#include "synquant_simple.h"

enum boundary_touched { NO_TOUCH = 0, XY_TOUCH = 1, Z_TOUCH = 2, XYZ_TOUCH = 3};
cellSegmentMain::cellSegmentMain(void *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
{
    init_parameter();
    data_rows_cols_slices = new int[3];

    if (bufSize[0] > INT_MAX || bufSize[1] > INT_MAX || bufSize[2] > INT_MAX){
        qDebug("Data is too large to process: %d!", INT_MAX);
    }
    data_rows_cols_slices[0] = bufSize[0];
    data_rows_cols_slices[1] = bufSize[1];
    data_rows_cols_slices[2] = bufSize[2];
    if (bufSize[3] != 1){
        //qFatal("Input is not a gray image\n");
        Q_ASSERT(bufSize[3] != 1);
    }
    time_points = bufSize[4];
    data_type = _data_type;
    cell_label_maps.resize(time_points);
    threshold_maps.resize(time_points);
    principalCurv2d.resize(time_points);
    principalCurv3d.resize(time_points);
    varMaps.resize(time_points);
    stblizedVarMaps.resize(time_points);
    varTrends.resize(time_points);
    stblizedVarTrends.resize(time_points);
    variances.resize(time_points);
    number_cells.resize(time_points);
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];

    Mat *data4d;
    int data_sz[4] = {data_rows_cols_slices[0], data_rows_cols_slices[1],
                      data_rows_cols_slices[2], (int)time_points};
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, data_sz, CV_16U, data_grayim4d);
//        int test_time = 1;
//        int sz_single_slice = data4d->size[1] * data4d->size[0];
//        for(int i = 0; i < data_sz[2]; i++){
//            unsigned short *ind = (unsigned short *)data4d->data
//                   + sz_single_slice * i + sz_single_frame*test_time; // sub-matrix pointer
//            Mat *single_slice = new Mat(2, data4d->size, CV_16U, ind);
//            imshow(to_string(i), *single_slice);
//            waitKey(0);
//        }
    }
    else if(data_type == V3D_UINT8){
        data4d = new Mat(4, data_sz, CV_8U, data_grayim4d);
    }else{
        qFatal("Unsupported data type\n");
    }
    //normalize(*data4d, *data4d, 0, 255, NORM_MINMAX, CV_8U);

//    Mat tmp;
//    data4d->copyTo(tmp);
//    qInfo("%d - %d - %d - %d \n", data4d->size[0], data4d->size[1],
//            data4d->size[2], data4d->size[3]);
//    normalize(tmp, normalized_data4d, 255, 0, NORM_MINMAX, CV_8U);
    Mat normalized_data4d;
    normalize(*data4d, normalized_data4d, 255, 0, NORM_MINMAX, CV_8U);
    assert(normalized_data4d.type() == CV_8U);
    for (int i = 0; i < time_points; i++){
        curr_time_point = i;
        unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*i; // sub-matrix pointer
        Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);
        cell_label_maps[i] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label
        threshold_maps[i] = Mat::zeros(3, normalized_data4d.size, CV_8U);
        cellSegmentSingleFrame(single_frame , i);
    }
}
void cellSegmentMain::init_parameter(){
    debug_folder = "/home/ccw/Desktop/embryo_res_folder/";
    default_name = debug_folder + "test.tiff";

    p4segVol.min_intensity = 0.0;
    p4segVol.fdr = .05;
    p4segVol.min_cell_sz =100;
    p4segVol.max_cell_sz = 3000;
    p4segVol.min_fill = 0.0001;
    p4segVol.max_WHRatio = 100;
    p4segVol.min_seed_size = 10;
    p4segVol.graph_cost_design[0] = ARITHMETIC_AVERAGE; //default 1, GEOMETRIC_AVERAGE = 2;
    p4segVol.graph_cost_design[1] = 2;
    p4segVol.growConnectInTest = 4;
    p4segVol.growConnectInRefine = 6;
    p4segVol.edgeConnect = 48;
    p4segVol.neiMap = 26;
    p4segVol.connect4fgGapRemoval = 26;
    p4segVol.shift_yxz[0] = 20;
    p4segVol.shift_yxz[1] = 20;
    p4segVol.shift_yxz[2] = 4;
    p4segVol.shrink_flag = true;
    p4segVol.shrink_scale_yxz[0] = 4;
    p4segVol.shrink_scale_yxz[1] = 4;
    p4segVol.shrink_scale_yxz[2] = 4;
    p4segVol.fgBoundaryHandle = LEAVEALONEFIRST;
    p4segVol.gapTestMinMaxRadius[0] = 2;
    p4segVol.gapTestMinMaxRadius[1] = 4;
    p4segVol.growSeedInTracking = false;

    p4odStats.gap4varTrendEst = 2;
    p4odStats.gap4fgbgCompare = 0;
    p4odStats.roundNum4fgbgCompare = 3;
    p4odStats.varAtRatio = 0.95;
    p4odStats.fgSignificanceTestWay = KSEC;
    p4odStats.minGapWithOtherCell_yxz[0] = 3;
    p4odStats.minGapWithOtherCell_yxz[1] = 3;
    p4odStats.minGapWithOtherCell_yxz[2] = 1;
    p4odStats.connectInSeedRefine = 6;
    p4odStats.gapTestMethod = GAP_LOCALORDERSTATS;
    p4odStats.gapTestSkippedBandWidth = 2;
    p4odStats.gapTestThreshold = 0.01;
}
void cellSegmentMain::reset_shift(){
    p4segVol.shift_yxz[0] = 20;
    p4segVol.shift_yxz[1] = 20;
    p4segVol.shift_yxz[2] = 4;
}
/**
 * @brief cellSegmentMain::cellSegmentSingleFrame
 * @param data_grayim3d: The data is mirrored by the y-direction, this is because
 * OpenGL texture are loaded left to right, bottom to top. Most image loaders
 * however will store the image in memory left to right, top to bottom.
 */
void cellSegmentMain::cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame)
{
    //data_grayim3d is uint8 0-255 datatype
    Mat *dataVolFloat = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    data_grayim3d->convertTo(*dataVolFloat, CV_32F);
    //ccShowSlice3Dmat(dataVolFloat, CV_32F);
    /******** start to do cell segmentation *******/
    float sigma2d[3] = {3.0, 3.0, 0.0};
    principalCv2d(dataVolFloat, principalCurv2d[curr_frame], sigma2d, p4segVol.min_intensity);
    float sigma3d[3] = {5.0, 5.0, 1.0};
    principalCv3d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, p4segVol.min_intensity);

    variances[curr_frame] = calVarianceStablization(dataVolFloat, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);

    Mat *stblizedVol = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    float stb_term = 3/8;
    FOREACH_i_ptrMAT(stblizedVol){
        stblizedVol->at<float>(i) = sqrt(dataVolFloat->at<float>(i) + stb_term);
    }
    calVarianceStablization(stblizedVol, stblizedVarMaps[curr_frame], stblizedVarTrends[curr_frame],
                                                       p4odStats.varAtRatio, p4odStats.gap4varTrendEst);
    // first use synQuant to get 1-tier seed regions
    synQuantSimple seeds_from_synQuant(data_grayim3d, variances[curr_frame], p4segVol, p4odStats);
    // second refine the seed regions
    vector<int> test_ids(0);
    regionWiseAnalysis4d(data_grayim3d, dataVolFloat, stblizedVol, seeds_from_synQuant.idMap/*int*/,
                         seeds_from_synQuant.cell_num, &principalCurv2d[curr_frame],
                         &principalCurv3d[curr_frame], &varMaps[curr_frame], &stblizedVarMaps[curr_frame],
                         test_ids);
    // third get 2-tier seed regions
    //seed_from_synQuant.retrieve_seeds();
    // fourth refine 2-tier seed regions

    // Finally, return
}

void boundaryTouchedTest(Mat *label_map, Mat *fgMap, bool &xy_touched, bool &z_touched){
    int n_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    int n_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    //n_z = {  0,  0,  0,  0, 0, 0,  0, 0 };
    int z_id, y, x;
    size_t remain;
    size_t page_sz = label_map->size[0] * label_map->size[1];
    FOREACH_i_ptrMAT(label_map){
        if(label_map->at<int>(i) > 0){
            z_id = i / page_sz;
            if(z_id == 0 || z_id == (label_map->size[2]-1)){
                z_touched = true;
                break;
            }
        }
    }
    int im_sz[] = {label_map->size[0], label_map->size[1]};
    FOREACH_i_ptrMAT(label_map){
        if(label_map->at<int>(i) > 0){
            z_id = i / page_sz;
            remain = i - z_id * page_sz;
            x = remain / label_map->size[0];
            y = remain - x * label_map->size[0];
            for(int i=0; i < 8; i++){
                if(inField(y + n_y[i], x + n_x[i], im_sz)
                        && fgMap->at<unsigned char>(y + n_y[i], x + n_x[i], z_id) == 0){
                    xy_touched = true;
                    return;
                }
            }
        }
    }
}
/**
 * @brief cellSegmentMain::regionWiseAnalysis4d, refine the region based on the results from synQuant
 * @param idMap
 * @param eigMap
 * @param vid
 * @param varMap
 * @param test_ids
 */
void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
                                           Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, vector<int> test_ids){
    //1. sort the seeds based on intensity levels
    vector<float> seed_intensity(seed_num);
    regionAvgIntensity(dataVolFloat, idMap, seed_intensity);
    vector<size_t> seed_intensity_order;
    seed_intensity_order = sort_indexes(seed_intensity, false); // false->descending
    //2. for each seed, refine it region
    vector<vector<size_t>> voxIdxList(seed_num);
    extractVoxIdxList(idMap, voxIdxList, seed_num);

    int cell_cnt = 0;
    FOREACH_i(seed_intensity_order){
        int seed_id = seed_intensity_order[i];
        singleCellSeed seed;
        cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, volStblizedFloat, idMap, eigMap2d,
                        eigMap3d, varMap, stblizedVarMap, seed, p4segVol);
        int touchBnd = refineSeed2Region(seed, p4odStats, p4segVol);
        if(touchBnd != NO_TOUCH){ // region too small
            if(touchBnd == XY_TOUCH || touchBnd == XYZ_TOUCH){
                p4segVol.shift_yxz[0] *=2;
                p4segVol.shift_yxz[1] *=2;
            }
            if(touchBnd == Z_TOUCH || touchBnd == XYZ_TOUCH){
                p4segVol.shift_yxz[2] *=2;
            }
            Range new_range[3];
            getRange(seed.y, p4segVol.shift_yxz[0], idMap->size[0], new_range[0]);
            getRange(seed.x, p4segVol.shift_yxz[1], idMap->size[1], new_range[1]);
            getRange(seed.z, p4segVol.shift_yxz[2], idMap->size[2], new_range[2]);
            if(new_range[0].size() > seed.fgMap.size[0] ||
                    new_range[1].size() > seed.fgMap.size[1] ||
                    new_range[2].size() > seed.fgMap.size[2]){

                cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, volStblizedFloat, idMap, eigMap2d,
                                eigMap3d, varMap, stblizedVarMap, seed, p4segVol);
                refineSeed2Region(seed, p4odStats, p4segVol);
            }
            reset_shift();
        }
        // assign cells to outmap
        Mat cropLabelMap = cell_label_maps[curr_time_point](seed.crop_range_yxz);
        FOREACH_i_MAT(seed.outputIdMap){
            if(seed.outputIdMap.at<int>(i) > 0){
                assert(cropLabelMap.at<int>(i)==0);
                cropLabelMap.at<int>(i) = cell_cnt + seed.outputIdMap.at<int>(i);
            }
        }
        cell_cnt += seed.outCell_num;
        Mat cropThresholdMap = threshold_maps[curr_time_point](seed.crop_range_yxz);

        FOREACH_i_MAT(seed.outputIdMap){
            if(seed.outputIdMap.at<int>(i) > 0){
                cropThresholdMap.at<int>(i) = seed.bestFgThreshold;
            }
        }
    }
    number_cells[curr_time_point] = cell_cnt;
}

void cellSegmentMain::cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d, Mat *data_stbized, Mat *idMap, Mat *eigMap2d,
                               Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, singleCellSeed &seed, segParameter p4segVol){
    seed.id = seed_id;
    seed.idx_yxz = idx_yxz; // shallow copy
    seed.y.resize(idx_yxz.size());
    seed.x.resize(idx_yxz.size());
    seed.z.resize(idx_yxz.size());
    vec_ind2sub(idx_yxz, seed.y, seed.x, seed.z, idMap->size);// MatSize itself is a vector

    getRange(seed.y, p4segVol.shift_yxz[0], idMap->size[0], seed.crop_range_yxz[0]);
    getRange(seed.x, p4segVol.shift_yxz[1], idMap->size[1], seed.crop_range_yxz[1]);
    getRange(seed.z, p4segVol.shift_yxz[2], idMap->size[2], seed.crop_range_yxz[2]);

    seed.seedMap = seed.idMap == seed.id;
    seed.eigMap2d = (*eigMap2d)(seed.crop_range_yxz); // all shallow copy
    seed.eigMap3d = (*eigMap3d)(seed.crop_range_yxz); // all shallow copy
    seed.varMap = (*varMap)(seed.crop_range_yxz); // all shallow copy
    seed.stblizedVarMap = (*stblizedVarMap)(seed.crop_range_yxz); // all shallow copy
    seed.gap2dMap = seed.eigMap2d > 0;
    seed.gap3dMap = seed.eigMap3d > 0;
    seed.volUint8 = (*data_grayim3d)(seed.crop_range_yxz); // all shallow copy
    seed.volStblizedFloat = (*data_stbized)(seed.crop_range_yxz); // all shallow copy
    seed.idMap = (*idMap)(seed.crop_range_yxz); // all shallow copy

    seed.outputIdMap = Mat(seed.idMap.dims, seed.idMap.size, CV_32S);

    vec_sub2ind(seed.idx_yxz_cropped, vec_Minus(seed.y, seed.crop_range_yxz[0].start),
            vec_Minus(seed.x, seed.crop_range_yxz[1].start),
            vec_Minus(seed.z, seed.crop_range_yxz[2].start), seed.idMap.size);// MatSize itself is a vector

    //seed.idMap.copyTo(seed.otherIdMap);
    seed.otherIdMap = seed.idMap > 0;
    FOREACH_i(seed.idx_yxz_cropped){
        seed.otherIdMap.at<int>(seed.idx_yxz_cropped[i]) = 0;
    }
    // generate the valid fg map
    Mat regMap = seed.idMap == seed_id;
    volumeDilate(&regMap/*cv_8u*/, seed.fgMap, p4segVol.shift_yxz, MORPH_ELLIPSE);
}

int cellSegmentMain::refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol){
    bool bnd_check = false;
    if (seed.bestFgThreshold < 0) bnd_check = true; // this is the first time we did seed refine
    // foreground Detection reusing synQuant
    synQuantSimple synQuant_refine_seed(seed, p4segVol, p4odStats);

    if(bnd_check){
        // if the fg from refinement touching the fg or not, if so, enlarge the fg
        bool xy_touched = false, z_touched = false;
        boundaryTouchedTest(synQuant_refine_seed.idMap, &seed.fgMap, xy_touched, z_touched);
        if(z_touched && xy_touched){
            return XYZ_TOUCH;
        }else if(z_touched){
            return Z_TOUCH;
        }else if(xy_touched){
            return XY_TOUCH;
        }else
            return NO_TOUCH;
    }else{
        return NO_TOUCH;
    }
}
