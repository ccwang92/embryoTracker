#include "cellsegment_main.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <chrono> // time elapsed
using namespace cv;
using namespace std;
//using namespace volproc;

cellSegmentMain::cellSegmentMain(void *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
{
    init_parameter();
    data_rows_cols_slices.resize(3);

    if (bufSize[0] > INT_MAX || bufSize[1] > INT_MAX || bufSize[2] > INT_MAX){
        qDebug("Data is too large to process: %d!", INT_MAX);
    }
    data_rows_cols_slices[0] = bufSize[1]; // !!! V3D transpose the image
    data_rows_cols_slices[1] = bufSize[0];
    data_rows_cols_slices[2] = bufSize[2];
    if (bufSize[3] != 1){
        //qFatal("Input is not a gray image\n");
        Q_ASSERT(bufSize[3] != 1);
    }
    time_points = bufSize[4];
    time_points_processed.resize(time_points);
    fill(time_points_processed.begin(), time_points_processed.end(), false);
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

    Mat *data4d;
    int data_sz[4] = {data_rows_cols_slices[0], data_rows_cols_slices[1],
                      data_rows_cols_slices[2], (int)time_points};
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, data_sz, CV_16U, data_grayim4d);
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
    //Mat normalized_data4d;
    normalize(*data4d, normalized_data4d, 255, 0, NORM_MINMAX, CV_8U);
    //ccShowSlice3Dmat(data4d, CV_16U);
    assert(normalized_data4d.type() == CV_8U);
}

void cellSegmentMain::processSingleFrameAndReturn(RayCastCanvas *glWidget, QString fileName){
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    curr_time_point = glWidget->curr_timePoint_in_canvas;
    ///
    /// \brief We fist test if there is saved binary file for use
    ///
    if(!fileName.isEmpty()){
        QString fileNameNoExt = fileName.left(fileName.lastIndexOf('.'));
        // read label data
        QString label_file_name = fileNameNoExt + "_label_map_int32.bin";
        QFile label_file(label_file_name);
        if (!label_file.open(QIODevice::ReadOnly)) return;
        //QByteArray tmp = label_file.readAll();
        // tmp is local variable, which will be released soon, so we need copyTo
        Mat(3, normalized_data4d.size, CV_32S, label_file.readAll().data()).copyTo(cell_label_maps[curr_time_point]);
        //cell_label_maps[curr_time_point] = Mat(3, normalized_data4d.size, CV_32S, label_file.readAll().data());
        //ccShowSliceLabelMat(cell_label_maps[curr_time_point]);
        label_file.close();

        // read threshold data
        QString threshod_file_name = fileNameNoExt + "_threshold_map_uint8.bin";
        QFile threshold_file(threshod_file_name);
        if (!threshold_file.open(QIODevice::ReadOnly)) return;
        Mat(3, normalized_data4d.size, CV_8U, threshold_file.readAll().data()).copyTo(threshold_maps[curr_time_point]);
        threshold_file.close();

        // read 2d principal map
        QString p2d_file_name = fileNameNoExt + "_principal2d_map_single.bin";
        QFile p2d_file(p2d_file_name);
        if (!p2d_file.open(QIODevice::ReadOnly)) return;
        Mat(3, normalized_data4d.size, CV_32F, p2d_file.readAll().data()).copyTo(principalCurv2d[curr_time_point]);
        p2d_file.close();

        // read 3d principal map
        QString p3d_file_name = fileNameNoExt + "_principal3d_map_single.bin";
        QFile p3d_file = QFile(p3d_file_name);
        if (!p3d_file.open(QIODevice::ReadOnly)) return;
        Mat(3, normalized_data4d.size, CV_32F, p3d_file.readAll().data()).copyTo(principalCurv3d[curr_time_point]);
        p3d_file.close();

        // read variance map
        QString varmap_file_name = fileNameNoExt + "_var_map_single.bin";
        QFile varmap_file = QFile(varmap_file_name);
        if (!varmap_file.open(QIODevice::ReadOnly)) return;
        Mat(3, normalized_data4d.size, CV_32F, varmap_file.readAll().data()).copyTo(varMaps[curr_time_point]);
        varmap_file.close();

        // read stablized variance map
        QString stbVarmap_file_name = fileNameNoExt + "_stb_var_map_single.bin";
        QFile stbVarmap_file = QFile(stbVarmap_file_name);
        if (!stbVarmap_file.open(QIODevice::ReadOnly)) return;
        Mat(3, normalized_data4d.size, CV_32F, stbVarmap_file.readAll().data()).copyTo(stblizedVarMaps[curr_time_point]);
        stbVarmap_file.close();

        // read variance trend
        QString vartrend_file_name = fileNameNoExt + "_var_trend_single.bin";
        QFile vartrend_file = QFile(vartrend_file_name);
        if (!vartrend_file.open(QIODevice::ReadOnly)) return;
        //QByteArray arr = vartrend_file.read(4);
        QByteArray tmp_in = vartrend_file.read(4);
//        QDataStream stream(tmp_in.data()); // #0 element
//        stream.setFloatingPointPrecision(QDataStream::SinglePrecision); // for float, default is double
//        stream >> variances[curr_time_point];
        memcpy(&variances[curr_time_point], tmp_in.data(), tmp_in.size());
        QByteArray arr = vartrend_file.readAll();// #1-200 elements
        varTrends[curr_time_point].resize(arr.size() / 4);
        memcpy(varTrends[curr_time_point].data(), arr.data(), arr.size());

        vartrend_file.close();

        // read stablized variance trend
        QString stbVartrend_file_name = fileNameNoExt + "_stb_var_trend_single.bin";
        QFile stbVartrend_file = QFile(stbVartrend_file_name);
        if (!stbVartrend_file.open(QIODevice::ReadOnly)) return;
        //QByteArray arr = vartrend_file.read(4);
        QByteArray arr2 = stbVartrend_file.readAll();// #1-200 elements
        stblizedVarTrends[curr_time_point].resize(arr2.size() / 4 - 1);
        memcpy(stblizedVarTrends[curr_time_point].data(), arr2.data()+4, arr2.size()); // skip the first element
        stbVartrend_file.close();

        time_points_processed[curr_time_point] = true;
    }


    ///
    /// \brief If no saved results, start from scratch to detect cells
    ///
    unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*curr_time_point; // sub-matrix pointer
    Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);
    if(!time_points_processed[curr_time_point]){
        cell_label_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label
        threshold_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_8U);
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        cellSegmentSingleFrame(single_frame, curr_time_point);
        time_points_processed[curr_time_point] = true;
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);
    }
    //ccShowSliceLabelMat(cell_label_maps[curr_time_point]);
    //// display volume in canvas
    Mat4b rgb_mat4display;
    label2rgb3d(cell_label_maps[curr_time_point], *single_frame, rgb_mat4display);
    glWidget->setMode("Alpha blending rgba");
    glWidget->getRenderer()->transfer_volume((unsigned char *)rgb_mat4display.data, 0, 255, data_rows_cols_slices[1],
            data_rows_cols_slices[0], data_rows_cols_slices[2], 4);
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
    //ccShowSlice3Dmat(data_grayim3d, CV_8U);
    /******** start to do cell segmentation *******/
    float sigma2d[3] = {3.0, 3.0, 0.0};
    principalCv2d(dataVolFloat, principalCurv2d[curr_frame], sigma2d, p4segVol.min_intensity);
    //ccShowSlice3Dmat(&principalCurv2d[curr_frame], CV_32F, 3);
    float sigma3d[3] = {5.0, 5.0, 1.0};
    principalCv3d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, p4segVol.min_intensity);
    //ccShowSlice3Dmat(&principalCurv3d[curr_frame], CV_32F);
    //Mat tmp = principalCurv3d[curr_frame]>0.01;
    //ccShowSlice3Dmat(tmp, CV_8U);
    double tmp_min, tmp_max;
    minMaxIdx(principalCurv2d[curr_frame], &tmp_min, &tmp_max);
    //ccShowSlice3Dmat(&principalCurv3d[curr_frame], CV_32F);
    variances[curr_frame] = calVarianceStablization(dataVolFloat, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);
    //ccShowSlice3Dmat(dataVolFloat, CV_32F, 3);
    //ccShowSlice3Dmat(&varMaps[curr_frame], CV_32F, 3);
    Mat *stblizedVol = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    float stb_term = 3/8;
    FOREACH_i_ptrMAT(stblizedVol){
        stblizedVol->at<float>(i) = sqrt(dataVolFloat->at<float>(i) + stb_term);
    }
    float stb_var = calVarianceStablization(stblizedVol, stblizedVarMaps[curr_frame], stblizedVarTrends[curr_frame],
                                                       p4odStats.varAtRatio, p4odStats.gap4varTrendEst);
    //ccShowSlice3Dmat(&stblizedVarMaps[curr_frame], CV_32F, 3);

    //////////////////////////////////////////////////////////////////
    //           1. use synQuant to get 1-tier seed regions         //
    //////////////////////////////////////////////////////////////////
    Mat smoothed_VolFloat, smoothed_grayim3d;
    dataVolFloat->copyTo(smoothed_VolFloat);
    float smooth4seedExtract[3] = {1.0, 1.0, 1.0};
    gaussianSmooth3Ddata(smoothed_VolFloat, smooth4seedExtract);
    smoothed_VolFloat.convertTo(smoothed_grayim3d, CV_8U);
    //ccShowSlice3Dmat(smoothed_grayim3d, CV_8U);
    synQuantSimple seeds_from_synQuant(&smoothed_grayim3d, variances[curr_frame], p4segVol, p4odStats);
    //ccShowSliceLabelMat(seeds_from_synQuant.idMap);
    //if directly use non-smoothed data
    //synQuantSimple seeds_from_synQuant(data_grayim3d, variances[curr_frame], p4segVol, p4odStats);

    //////////////////////////////////////////////////////////////////
    //                 2. refine the seed regions                   //
    //////////////////////////////////////////////////////////////////
    bool RoundOne = true;
    regionWiseAnalysis4d(data_grayim3d, dataVolFloat, stblizedVol, seeds_from_synQuant.idMap/*int*/,
                         seeds_from_synQuant.cell_num, curr_frame, RoundOne);
    //ccShowSliceLabelMat(cell_label_maps[curr_frame]);
    //////////////////////////////////////////////////////////////////
    //               3. get 2-tier seed regions                  //
    //////////////////////////////////////////////////////////////////
    Mat cellGapMap = principalCurv3d[curr_frame] > 0.001;
    Mat label_map_2ndRound;
    int seed_num_2ndRound;
    retrieve_seeds(&smoothed_VolFloat, &cell_label_maps[curr_frame], number_cells[curr_frame],
                   &cellGapMap, label_map_2ndRound, seed_num_2ndRound);

    //////////////////////////////////////////////////////////////////
    //              4. refine 2-tier seed regions                   //
    //////////////////////////////////////////////////////////////////
    if(seed_num_2ndRound > 0){
        RoundOne = false;
        regionWiseAnalysis4d(data_grayim3d, dataVolFloat, stblizedVol, &label_map_2ndRound/*int*/,
                             seed_num_2ndRound, curr_frame, RoundOne);
    }

}
/**
 * @brief retrieve_seeds: retrieve the seeds of cells lost in the first round of synQuant
 * @param smoothed_VolFloat: 0-255 but in float format
 * @param label_map_1stRound
 * @param cell_num_1stRound
 * @param test_ids
 */
void cellSegmentMain::retrieve_seeds(Mat *smoothed_VolFloat, Mat *label_map_1stRound, size_t cell_num_1stRound,
                                     Mat *cellGapMap, Mat &idMap_2ndRound, int &seed_num_2ndRound){
    // 1. extract the cell size and intensity
    vector<vector<size_t>> voxList_exist_cells;
    extractVoxIdxList(label_map_1stRound, voxList_exist_cells, cell_num_1stRound);
    vector<float> cell_mean_intensity;
    regionAvgIntensity(smoothed_VolFloat, voxList_exist_cells, cell_mean_intensity);
    vector<size_t> exist_cell_size (voxList_exist_cells.size());
    FOREACH_i(voxList_exist_cells) exist_cell_size[i] = voxList_exist_cells[i].size();

    size_t min_exist_cell_size = *min_element(exist_cell_size.begin(), exist_cell_size.end());
    size_t max_exist_cell_size = *max_element(exist_cell_size.begin(), exist_cell_size.end());
    float min_cell_intensity = *min_element(cell_mean_intensity.begin(), cell_mean_intensity.end());

    // 2. extract the new seed map if it satisifies the constraints from existing cells
    Mat seed_map;
    Mat exist_cell_territory, exist_cell_fg;
    exist_cell_fg = *label_map_1stRound > 0;
    int dilate_radius[] = {2, 2, 1};
    volumeDilate(&exist_cell_fg, exist_cell_territory, dilate_radius, MORPH_ELLIPSE);
    //setValMat(*smoothed_VolFloat, CV_32F, &exist_cell_territory, 0);
    //ccShowSlice3Dmat(smoothed_VolFloat, CV_32F);
    //volumeWrite(smoothed_VolFloat, "/home/ccw/Desktop/remain.tif");
    //ccShowSlice3Dmat(cellGapMap, CV_8U);
    //volumeWrite(&exist_cell_territory, "/home/ccw/Desktop/cell_territory.tif");
    bitwise_or(*cellGapMap, exist_cell_territory, exist_cell_territory);
    //ccShowSlice3Dmat(exist_cell_territory, CV_8U);
    //setValMat(*smoothed_VolFloat, CV_32F, &exist_cell_territory, 0);
    //volumeWrite(smoothed_VolFloat, "/home/ccw/Desktop/remain_further.tif");
    //ccShowSlice3Dmat(smoothed_VolFloat, CV_32F);
    //volumeWrite(cellGapMap, "/home/ccw/Desktop/gap.tif");
    //volumeWrite(&exist_cell_territory, "/home/ccw/Desktop/overlap.tif");
    seed_map = exist_cell_territory == 0;
    dilate_radius[0] = 1; dilate_radius[1] = 1; dilate_radius[2] = 0;
    Mat seed_map_dilate;
    volumeDilate(&seed_map, seed_map_dilate, dilate_radius, MORPH_ELLIPSE);
    Mat label_map;
    int numCC;
    numCC = connectedComponents3d(&seed_map_dilate, label_map, 26);
    setValMat(label_map, CV_32S, &exist_cell_territory, 0);
    //ccShowSliceLabelMat(label_map);
    vector<vector<size_t>> voxList_seeds;
    extractVoxIdxList(&label_map, voxList_seeds, numCC);
    vector<float> seeds_mean_intensity;
    regionAvgIntensity(smoothed_VolFloat, voxList_seeds, seeds_mean_intensity);
    idMap_2ndRound = Mat::zeros(label_map.dims, label_map.size, CV_32S);
    seed_num_2ndRound = 0;
    FOREACH_i(voxList_seeds){
        if(voxList_seeds[i].size() >= min_exist_cell_size &&
                voxList_seeds[i].size() <= max_exist_cell_size &&
                seeds_mean_intensity[i] > min_cell_intensity){
            seed_num_2ndRound ++;
            setValMat(idMap_2ndRound, CV_32S, voxList_seeds[i], (float)seed_num_2ndRound);
        }
    }
    //ccShowSliceLabelMat(idMap_2ndRound);
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
            y = remain / label_map->size[1];
            x = remain - y * label_map->size[1];
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
//void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
//                                           Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, vector<int> test_ids){
void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat,
                                           Mat * volStblizedFloat, Mat *idMapIn /*int*/, int seed_num,
                                           int curr_frame, bool roundOne){
    Mat idMap;
    idMapIn->copyTo(idMap);
    if(roundOne){
        removeSmallCC(idMap, seed_num, p4segVol.min_seed_size, true);
    }
    //1. sort the seeds based on intensity levels
    vector<float> seed_intensity;
    regionAvgIntensity(dataVolFloat, &idMap, seed_num, seed_intensity);
    //ccShowSlice3Dmat(dataVolFloat, CV_32F);
    vector<size_t> seed_intensity_order;
    seed_intensity_order = sort_indexes(seed_intensity, false); // false->descending
    //2. for each seed, refine it region
    vector<vector<size_t>> voxIdxList(seed_num);
    extractVoxIdxList(&idMap, voxIdxList, seed_num);

    int cell_cnt = 0;
    if(!roundOne) { // for the second tround
        //1. cell number start from the second round
        cell_cnt += number_cells[curr_time_point];
        //2. idMap needs add other existing cells (label start from seed_num+1)
        FOREACH_i_MAT(cell_label_maps[curr_time_point]){
            if(cell_label_maps[curr_time_point].at<int>(i) > 0){
                idMap.at<int>(i) = seed_num + cell_label_maps[curr_time_point].at<int>(i);
            }
        }
    }
    //int debug_cell_id = 2;
    FOREACH_i(seed_intensity_order){
        int seed_id = seed_intensity_order[i] + 1;
        qInfo("------------------#%ld, start %d seed process, size: %ld---------------------", i, seed_id,
              voxIdxList[seed_id-1].size());

        singleCellSeed seed;
        cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, volStblizedFloat, &idMap,
                curr_frame, seed, p4segVol);
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
            getRange(seed.y, p4segVol.shift_yxz[0], idMap.size[0], new_range[0]);
            getRange(seed.x, p4segVol.shift_yxz[1], idMap.size[1], new_range[1]);
            getRange(seed.z, p4segVol.shift_yxz[2], idMap.size[2], new_range[2]);
            if(new_range[0].size() > seed.validSearchAreaMap.size[0] ||
                    new_range[1].size() > seed.validSearchAreaMap.size[1] ||
                    new_range[2].size() > seed.validSearchAreaMap.size[2]){

                cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, volStblizedFloat, &idMap,
                        curr_frame, seed, p4segVol);
//                if(debug_cell_id == i){
//                    qInfo("may not correct");
//                    ccShowSliceLabelMat(seed.idMap);
//                }
                refineSeed2Region(seed, p4odStats, p4segVol);
//                if(debug_cell_id == i){
//                    qInfo("may not correct");
//                    ccShowSliceLabelMat(seed.outputIdMap);
//                }
            }
            reset_shift();
        }
        qInfo("----------------%d cell output from this seed-------------------", seed.outCell_num);
        //ccShowSliceLabelMat(seed.outputIdMap);
        // assign cells to outmap
//        Mat cropLabelMap = cell_label_maps[curr_time_point](seed.crop_range_yxz);
//        FOREACH_i_MAT(seed.outputIdMap){
//            if(seed.outputIdMap.at<int>(i) > 0){
//                assert(cropLabelMap.at<int>(i)==0);
//                cropLabelMap.at<int>(i) = cell_cnt + seed.outputIdMap.at<int>(i);
//            }
//        }
//        if(debug_cell_id == i){
//            ccShowSliceLabelMat(seed.outputIdMap);
//        }
        subVolReplace(cell_label_maps[curr_time_point], CV_32S, seed.outputIdMap, seed.crop_range_yxz, cell_cnt);
        cell_cnt += seed.outCell_num;

        Mat valid_threshold_map = seed.outputIdMap > 0;
        subVolReplace(threshold_maps[curr_time_point], CV_8U, valid_threshold_map, (float)seed.bestFgThreshold, seed.crop_range_yxz);
//        Mat cropThresholdMap = threshold_maps[curr_time_point](seed.crop_range_yxz);
//        FOREACH_i_MAT(seed.outputIdMap){
//            if(seed.outputIdMap.at<int>(i) > 0){
//                cropThresholdMap.at<int>(i) = seed.bestFgThreshold;
//            }
//        }
        //break; // !!!! just for debug
    }
    //ccShowSliceLabelMat(cell_label_maps[curr_time_point]);
    //ccShowSlice3Dmat(threshold_maps[curr_time_point], CV_8U);
    number_cells[curr_time_point] = cell_cnt;
}

void cellSegmentMain::cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d,
                               Mat *data_stbized, Mat *idMap, int curr_frame, singleCellSeed &seed,
                               segParameter p4segVol){
    bool seed_extract_in_tracking =false;
    if(data_stbized == nullptr){
        seed_extract_in_tracking = true;
    }
    seed.id = seed_id;
    seed.idx_yxz = idx_yxz; // deep copy
    seed.y.resize(idx_yxz.size());
    seed.x.resize(idx_yxz.size());
    seed.z.resize(idx_yxz.size());
    vec_ind2sub(idx_yxz, seed.y, seed.x, seed.z, idMap->size);// MatSize itself is a vector
//    ccShowSliceLabelMat(idMap);
    getRange(seed.y, p4segVol.shift_yxz[0], idMap->size[0], seed.crop_range_yxz[0]);
    getRange(seed.x, p4segVol.shift_yxz[1], idMap->size[1], seed.crop_range_yxz[1]);
    getRange(seed.z, p4segVol.shift_yxz[2], idMap->size[2], seed.crop_range_yxz[2]);

//// debug of label2rgb function
    //ccShowSliceLabelMat(idMap);
//    Mat tmp = idMap->clone();
//    tmp.convertTo(tmp, CV_8U);
//    Range xyz[3];
//    xyz[0] = Range(0, 125);
//    xyz[1] = Range(0, 250);
//    xyz[2] = Range(0, 3);

//    Mat subMat = tmp(xyz).clone();
//    int *ind = (int*) tmp.data;idMap
//    Mat single_s = Mat(2, tmp.size, CV_8U, ind);
//    imshow("test_", single_s>0);
//    waitKey(0);
//    ccShowSliceLabelMat(&subMat);

//    Range xy[2];
//    xy[0] = Range(0, 100);
//    xy[1] = Range(0, 125);
//    int *ind = (int*) idMap->data;
//    Mat single_s = Mat(2, idMap->size, CV_32S, ind);
//    //ccShowSliceLabelMat(&single_s);

////    cv::Rect roi(0, 0, 125, 100);
////    Mat1i subMat2 = single_s;
//    //single_s.convertTo(single_s, CV_32S);
//    Mat subMat3 = single_s(xy).clone();
//    //imshow("test", subMat3 > 0);
//    //waitKey(0);
    //ccShowSliceLabelMat(idMap);
//    Range xyz[3];
//    xyz[0] = Range(0, 100);
//    xyz[1] = Range(0, 250);
//    xyz[2] = Range(0, 1);
//    int size[3] = {xyz[0].end - xyz[0].start,
//                  xyz[1].end - xyz[1].start,
//                  xyz[2].end - xyz[2].start};
//    seed.idMap.create(3, size, CV_32S);
    subVolExtract(idMap, CV_32S, seed.idMap, seed.crop_range_yxz);
    //ccShowSliceLabelMat(&seed.idMap);
    //seed.idMap = (*idMap)(seed.crop_range_yxz).clone(); // all shallow copy
    //ccShowSliceLabelMat(&seed.idMap);
    seed.seedMap = seed.idMap == seed.id;

//    seed.eigMap2d = principalCurv2d[curr_frame](seed.crop_range_yxz); // all shallow copy
//    seed.eigMap3d = principalCurv3d[curr_frame](seed.crop_range_yxz); // all shallow copy
//    seed.varMap = varMaps[curr_frame](seed.crop_range_yxz); // all shallow copy
//    seed.stblizedVarMap = stblizedVarMaps[curr_frame](seed.crop_range_yxz); // all shallow copy
    subVolExtract(&principalCurv2d[curr_frame], CV_32F, seed.eigMap2d, seed.crop_range_yxz);// deep copy
    subVolExtract(&principalCurv3d[curr_frame], CV_32F, seed.eigMap3d, seed.crop_range_yxz);// deep copy
    subVolExtract(&varMaps[curr_frame], CV_32F, seed.varMap, seed.crop_range_yxz);// deep copy
    if(!seed_extract_in_tracking){
        subVolExtract(&stblizedVarMaps[curr_frame], CV_32F, seed.stblizedVarMap, seed.crop_range_yxz);// deep copy
    }
    //ccShowSlice3Dmat(&seed.stblizedVarMap, CV_32F);
    seed.gap2dMap = seed.eigMap2d > 0;
    seed.gap3dMap = seed.eigMap3d > 0;
//    seed.volUint8 = (*data_grayim3d)(seed.crop_range_yxz); // all shallow copy
//    seed.volStblizedFloat = (*data_stbized)(seed.crop_range_yxz); // all shallow copy
    subVolExtract(data_grayim3d, CV_8U, seed.volUint8, seed.crop_range_yxz);// deep copy
    if(!seed_extract_in_tracking){
        subVolExtract(data_stbized, CV_32F, seed.volStblizedFloat, seed.crop_range_yxz);// deep copy
    }
//    ccShowSliceLabelMat(&seed.idMap, 6);
//    ccShowSlice3Dmat(&seed.seedMap, CV_8U, 6);
//    ccShowSlice3Dmat(&seed.volUint8, CV_8U, 6);
    seed.outputIdMap = Mat(seed.idMap.dims, seed.idMap.size, CV_32S);
    //ccShowSliceLabelMat(idMap);
    vec_sub2ind(seed.idx_yxz_cropped, vec_Minus(seed.y, seed.crop_range_yxz[0].start),
            vec_Minus(seed.x, seed.crop_range_yxz[1].start),
            vec_Minus(seed.z, seed.crop_range_yxz[2].start), seed.idMap.size);// MatSize itself is a vector
    //ccShowSliceLabelMat(idMap);
    //seed.idMap.copyTo(seed.otherIdMap);
    seed.otherIdMap = seed.idMap > 0;
    FOREACH_i(seed.idx_yxz_cropped){
        seed.otherIdMap.at<unsigned char>(seed.idx_yxz_cropped[i]) = 0;
    }
    //ccShowSliceLabelMat(idMap);
    // generate the valid fg map
    //Mat regMap = (seed.idMap == seed_id);

    //ccShowSliceLabelMat(idMap);
    //p4segVol.shift_yxz[2] = 0;
    //// NOTE: if there is a segmentation error, check first if the Mat format is wrongly
    /// assigned. This is a problem of using opencv, where data format is predefined everywhere.
    //seed.validSearchAreaMap.create(seed.seedMap.dims, seed.seedMap.size, CV_8U);
    volumeDilate(&seed.seedMap/*cv_8u*/, seed.validSearchAreaMap, p4segVol.shift_yxz, MORPH_ELLIPSE);
    //ccShowSlice3Dmat(&seed.idMap, CV_32S, 7);
//    ccShowSlice3Dmat(data_grayim3d, CV_8U);
//    ccShowSliceLabelMat(idMap);
//    ccShowSliceLabelMat(&seed.idMap);
//    ccShowSlice3Dmat(&seed.seedMap, CV_8U, 7);
//    ccShowSlice3Dmat(&seed.validSearchAreaMap, CV_8U);
}

int cellSegmentMain::refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol){
    bool bnd_check = false;
    if (seed.bestFgThreshold < 0) bnd_check = true; // this is the first time we did seed refine
    // foreground Detection reusing synQuant
    synQuantSimple cellSegFromSynQuant(seed);
    //// 1. find best threshold to get fgMap and apply some basic morphlogical operation (ordstat4fg.m)
    // PS: the region from this step will has no overlap with other cell's terriotory
    // PPS: the cell territory may contain more than one cell, we need verify it using gaps.
    if(seed.bestFgThreshold <= 0){
        // we have no idea what is the best intensity threshold to get the region
        // we use this funcion to get this.fgMap
        cellSegFromSynQuant.cellTerritoryExtractFromSeed(seed, p4odStats);
    } else{
        // we have alreay have some idea what is the best intensity level
        // we use thresholding to get this.fgMap
        bitwise_and(seed.validSearchAreaMap, seed.volUint8 >= seed.bestFgThreshold, cellSegFromSynQuant.fgMap);
        //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    }
    qInfo("after refine, fgmap is empty?: %d", isempty(cellSegFromSynQuant.fgMap, CV_8U));

    refineCellTerritoryWithSeedRegion(cellSegFromSynQuant, seed, p4segVol);
    qInfo("after refine, fgmap is empty?: %d", isempty(cellSegFromSynQuant.fgMap, CV_8U));

    //// 2. update seed's score map based on idMap (fg)
    normalize(seed.eigMap2d, seed.score2d, 0.001, 1, NORM_MINMAX, CV_32F, cellSegFromSynQuant.fgMap);
    normalize(seed.eigMap3d, seed.score3d, 0.001, 1, NORM_MINMAX, CV_32F, cellSegFromSynQuant.fgMap);
    seed.scoreMap = seed.score2d + seed.score3d;
    //// 2.1 append step for simple thresholding step
    if(seed.bestFgThreshold > 0){
        // fgMap from simple thresholding may cover >1 cell seeds.
        //qInfo("after change shift scale, fgmap valid?: %d", isempty(cellSegFromSynQuant.fgMap, CV_8U));
        //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
        removeOtherSeedsInfgMap(cellSegFromSynQuant, seed, p4segVol);
        qInfo("after remove other seeds, fgmap empty?: %d", isempty(cellSegFromSynQuant.fgMap, CV_8U));
    }
    //// 3. segment fgMap into idMap based on gaps from principal curvature
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    //ccShowSlice3Dmat(seed.seedMap, CV_8U);
    gapBasedRegionSegment(cellSegFromSynQuant, seed, p4segVol, p4odStats);

    //// 4. refine the idMap from step 3 based on size and other prior knowledge
    if(cellSegFromSynQuant.cell_num > 1){
        bool link_bg2sink = false; // this can force seeds to grow as much as they can
//        Mat input_seeds;
//        cellSegFromSynQuant.idMap->copyTo(input_seeds);
        //Mat scoreMap = seed.score2d + seed.score3d;
        regionGrow(cellSegFromSynQuant.idMap, cellSegFromSynQuant.cell_num,
                   seed.outputIdMap, &seed.scoreMap, &cellSegFromSynQuant.fgMap, //idMap is fully included in fgMap
                   p4segVol.growConnectInRefine, p4segVol.graph_cost_design,
                   link_bg2sink);
        //TODO: re-assign extra voxels missed in idMap, but contained in this.fgMap
    }else if(cellSegFromSynQuant.cell_num == 1){
        FOREACH_i_MAT(seed.outputIdMap){
            if(cellSegFromSynQuant.fgMap.at<unsigned char>(i) > 0){
                seed.outputIdMap.at<int>(i) = 1;
            }else{
                seed.outputIdMap.at<int>(i) = 0;
            }
        }
    }else{
        qFatal("We did not find any cell from this seed.");
    }

    if(cellSegFromSynQuant.cell_num > 0){ // shrink test
        cellShrinkTest(cellSegFromSynQuant, seed, p4segVol);
    }
    removeSmallCC(seed.outputIdMap, cellSegFromSynQuant.cell_num, p4segVol.min_cell_sz, true);

    seed.outCell_num = cellSegFromSynQuant.cell_num;
    //cellSegFromSynQuant.idMap->copyTo(seed.outputIdMap); // output
    if(seed.bestFgThreshold <= 0){ // the first run which have not got best threshold
        seed.bestFgThreshold = (int) cellSegFromSynQuant.maxZ_intensity_level;
    }



    if(bnd_check){
        // if the fg from refinement touching the fg or not, if so, enlarge the fg
        bool xy_touched = false, z_touched = false;
        boundaryTouchedTest(cellSegFromSynQuant.idMap, &seed.validSearchAreaMap, xy_touched, z_touched);
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



/**
 * @brief removeOtherSeedsInfgMap: if we use thresholding, the fgMap may contains more than one seed. This will
 * not happen if we consider othercellterritory as in cellTerritoryExtractFromSeed();
 * @param seed
 * @param p4segVol
 */
void cellSegmentMain::removeOtherSeedsInfgMap(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    //// 1. find if there is other cells in fgMap
    vector<size_t> fg_idx = fgMapIdx(&cellSegFromSynQuant.fgMap, CV_8U, 0);
    bool other_cell_exist = false;
    FOREACH_i(fg_idx){
        if(seed.idMap.at<int>(fg_idx[i]) > 0 &&
                seed.idMap.at<int>(fg_idx[i]) != seed.id){
            other_cell_exist = true;
            break;
        }
    }
    if(!other_cell_exist){
        return;
    }
    Mat seedMap = Mat::zeros(cellSegFromSynQuant.fgMap.dims, cellSegFromSynQuant.fgMap.size, CV_32S);
    FOREACH_i(fg_idx){
        if(seed.idMap.at<int>(fg_idx[i]) > 0){
            if(seed.idMap.at<int>(fg_idx[i]) == seed.id){
                seedMap.at<int>(fg_idx[i]) = 1;
            }else{
                seedMap.at<int>(fg_idx[i]) = 2;
            }
        }else if(isOnBoundary2d(&seed.validSearchAreaMap, fg_idx[i])){
            // if the boundary of init foreground is also contained, label as sink
            // this is only considered when we double seed.shift_yxz, since the area are
            // too large now. The boudnary of init foreground should not be part of target
            // cell. Before doubling seed.shift_yxz, it is possible.
            seedMap.at<int>(fg_idx[i]) = 2;
        }
    }
    Mat grownSeedMap2d, grownSeedMap3d;
    bool bg2sink = true;
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    //ccShowSliceLabelMat(seedMap);
    regionGrow(&seedMap, 2, grownSeedMap2d, &seed.score2d, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInTest, p4segVol.graph_cost_design, bg2sink);
    //ccShowSliceLabelMat(grownSeedMap2d);
    bg2sink = false;
    regionGrow(&grownSeedMap2d, 2, grownSeedMap3d, &seed.scoreMap, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInRefine, p4segVol.graph_cost_design, bg2sink);
    //ccShowSliceLabelMat(grownSeedMap3d);
    cellSegFromSynQuant.fgMap = grownSeedMap3d == 1;
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
}
void cellSegmentMain::cellShrinkTest(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    int extra_cell = 0;
    for(int i = 1; i <= cellSegFromSynQuant.cell_num; i++){
        Mat cur_cell = seed.outputIdMap == i;
        Mat shrinked_cell, label_shrinked_cell;
        volumeErode(&cur_cell, shrinked_cell, p4segVol.shrink_scale_yxz, MORPH_ELLIPSE);
//        if(seed.crop_range_yxz[0].end == 94 && seed.crop_range_yxz[1].end == 123){
//            ccShowSlice3Dmat(&shrinked_cell, CV_8U);
//        }
        int n = connectedComponents3d(&shrinked_cell, label_shrinked_cell, p4segVol.growConnectInRefine);

        if(n > 1){
            removeSmallCC(label_shrinked_cell, n, p4segVol.min_seed_size, true);
            if(n > 1){
                Mat grown_shrinked_cells;
                bool link_bg2sink = false;
                regionGrow(&label_shrinked_cell, n, grown_shrinked_cells, &seed.scoreMap, &cur_cell,
                           p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);
                setValMat(seed.outputIdMap, CV_32S, &cur_cell, 0.0);
                Mat sub_cell_mask = grown_shrinked_cells == 1;
                setValMat(seed.outputIdMap, CV_32S, &sub_cell_mask, (float)i);

                for(int j=2; j<=n; j++){
                    sub_cell_mask = grown_shrinked_cells == j;
                    extra_cell ++;
                    setValMat(seed.outputIdMap, CV_32S, &sub_cell_mask, (float)(extra_cell + cellSegFromSynQuant.cell_num));
                }
            }
        }
    }
    cellSegFromSynQuant.cell_num += extra_cell;
}

/**
 * @brief refineFgWithSeedRegion, extract the largest region if fg and remove those
 * voxels that has no z-direction neighbors
 * @param seed
 * @param p4segVol
 */
void cellSegmentMain::refineCellTerritoryWithSeedRegion(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    Mat labelMap;
    int n = connectedComponents3d(&cellSegFromSynQuant.fgMap, labelMap, 6);
    //ccShowSlice3Dmat(seed.seedMap, CV_8U);

    if (n > 1){
        int id = largestRegionIdExtract(&labelMap, n, &seed.seedMap);
        cellSegFromSynQuant.fgMap = labelMap == id;
        //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    }
    // remove pixels that happen only at one slice
    Mat tmp_map;
    cellSegFromSynQuant.fgMap.copyTo(tmp_map);
    size_t page_sz = cellSegFromSynQuant.fgMap.size[1] * cellSegFromSynQuant.fgMap.size[0];
    int width = cellSegFromSynQuant.fgMap.size[1];
    FOREACH_ijk_ptrMAT(cellSegFromSynQuant.idMap){
        size_t idx = vol_sub2ind(i,j,k, width, page_sz);
        if(tmp_map.at<unsigned char>(idx) > 0){
            if((k-1>=0 && tmp_map.at<unsigned char>(idx - page_sz) == 0) &&
                (k+1<tmp_map.size[2] && tmp_map.at<unsigned char>(idx + page_sz) == 0)){
                cellSegFromSynQuant.fgMap.at<unsigned char>(idx) = 0;
            }
        }
    }
    int radius[] = {1,1,0};
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    volumeDilate(&cellSegFromSynQuant.fgMap, tmp_map, radius, MORPH_ELLIPSE);
    n = connectedComponents3d(&tmp_map, labelMap, 26);
    removeSmallCC(labelMap, n, p4segVol.min_cell_sz, false);
    cellSegFromSynQuant.fgMap = labelMap > 0;
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
}

/**
 * @brief fgGapRemoval: remove the gap defined by principal curvature. We jointly consider the 2d and 3d
 * principal curvature, since 3d may be too liberal to remove too much areas
 * @param seed
 * @param p4segVol
 */
void cellSegmentMain::fgGapRemoval(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    bitwise_and(cellSegFromSynQuant.fgMap, seed.gap3dMap == 0, cellSegFromSynQuant.fgMapGapRemoved);

    if (!isempty(&cellSegFromSynQuant.fgMapGapRemoved, CV_8U)){
        Mat seedMapFrom2dMap, mapUnion, mapUnion_label;
        bitwise_and(cellSegFromSynQuant.fgMap, seed.gap2dMap == 0, seedMapFrom2dMap);// or use -
        bitwise_or(seedMapFrom2dMap, cellSegFromSynQuant.fgMapGapRemoved, mapUnion); // or use +
        int numCC = connectedComponents3d(&mapUnion, mapUnion_label, p4segVol.connect4fgGapRemoval);
        Mat newSeedMap; // CV_8U
        bool found = findUnrelatedCC(&mapUnion_label, numCC, &cellSegFromSynQuant.fgMapGapRemoved, newSeedMap);
        if(found){
            cellSegFromSynQuant.fgMapGapRemoved = cellSegFromSynQuant.fgMapGapRemoved + newSeedMap; //CV_8U
        }
    }
}
/**
 * @brief gapBasedRegionSegment: use principal curvature to test if current fg contains > 1 cells.
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void cellSegmentMain::gapBasedRegionSegment(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats){
    // 1. remove gaps defined by principal curvature
    fgGapRemoval(cellSegFromSynQuant, seed, p4segVol);
    // 2. check if gaps divided the region into multiple valid seeds
    Mat label_map;
    int n = connectedComponents3d(&cellSegFromSynQuant.fgMapGapRemoved, label_map, p4segVol.neiMap);
    removeSmallCC(label_map, n, p4segVol.min_seed_size, REARRANGE_IDS);
    if(n <= 1){
        cellSegFromSynQuant.cell_num = 1;
        return;
    }
    // 3. gap test: test if the gap is true, if so, split the fg with respect to the multiple seeds
    gapTest2SplitCellTerritory(cellSegFromSynQuant, &label_map, n, seed, p4segVol, p4odStats);
}
/**
 * @brief gapTest2SplitCellTerritory: Given a seeds_map, test if the gap between any two seeds are true or false.
 * If true, keep them. Otherwise, merge the two seeds divided by the gap. If the gap is indeed dimmer than the
 * two seeds it split, it is true. Otherwise false.
 * @param seeds_Map:cv_32s
 * @param n
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void cellSegmentMain::gapTest2SplitCellTerritory(synQuantSimple &cellSegFromSynQuant, Mat* seeds_Map /*CV_32S*/, int n,
                                                singleCellSeed &seed, segParameter &p4segVol,
                                                odStatsParameter &p4odStats){
    double maxId;
    minMaxIdx(*seeds_Map, nullptr, &maxId);
    // 1. first grow all seeds until they touch with each other
    Mat grown_seedMap;
    //Mat scoreMap = seed.score2d + seed.score3d;
    bool link_bg2sink = false; // this can force seeds to grow as much as they can
    regionGrow(seeds_Map, n, grown_seedMap, &seed.scoreMap, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);
    //ccShowSliceLabelMat(seeds_Map);
    //ccShowSliceLabelMat(&grown_seedMap);


    // 2. for the connected seeds, test if the gap among them are true
    // the gap are choosen by dilating the seed regions
    vector<int> gap_tested_true(n*n, 0); // 0 not tested, -1 tested but merge, 1 tested split
    float p_treshold = p4odStats.gapTestThreshold / (p4segVol.gapTestMinMaxRadius[1] - p4segVol.gapTestMinMaxRadius[0]);
    vector<size_t> real_gap_idx;
    for(int r = p4segVol.gapTestMinMaxRadius[0]; r <= p4segVol.gapTestMinMaxRadius[1]; r++){
        vector<vector<size_t>> gap_idx_list;
        extractGapVoxel(&grown_seedMap, &cellSegFromSynQuant.fgMap, n, r, gap_idx_list, gap_tested_true);
        FOREACH_i(gap_idx_list){
            if(gap_tested_true[i] != 1 && gap_idx_list[i].size() > 10){
                int id0 = i / n + 1;
                int id1 = i % n + 1;
                vector<vector<size_t>> nei_list0, nei_list1;
                Mat cur_fg = grown_seedMap == id0;
                neighbor_idx_2d(gap_idx_list[i], &cur_fg, nei_list0, 2*r+2);
                cur_fg = grown_seedMap == id1;
                neighbor_idx_2d(gap_idx_list[i], &cur_fg, nei_list1, 2*r+2);

                double var_sum = 0;
                size_t var_idx_num = 0;
                vector<float> r0, r1, skipped, gap;

                gap = extractValsGivenIdx(&seed.volStblizedFloat, gap_idx_list[i], CV_32F);
                //ccShowSlice3Dmat(&seed.volStblizedFloat, CV_32F);
                var_sum += extractSumGivenIdx(&seed.stblizedVarMap, gap_idx_list[i], CV_32F);
                var_idx_num += gap_idx_list[i].size();
                for(size_t j = p4odStats.gapTestSkippedBandWidth; j < nei_list0.size(); j++){
                    vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list0[j], CV_32F);
                    r0.insert(r0.end(), curr_vals.begin(), curr_vals.end());
                    var_sum += extractSumGivenIdx(&seed.stblizedVarMap, nei_list0[j], CV_32F);
                    var_idx_num += nei_list0[j].size();
                }
                for(size_t j = p4odStats.gapTestSkippedBandWidth; j < nei_list1.size(); j++){
                    vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list1[j], CV_32F);
                    r1.insert(r1.end(), curr_vals.begin(), curr_vals.end());
                    var_sum += extractSumGivenIdx(&seed.stblizedVarMap, nei_list1[j], CV_32F);
                    var_idx_num += nei_list1[j].size();
                }
                float cur_std = (float)sqrt(var_sum / var_idx_num);
                float mu, sigma;
                float p0, p1;
                if(p4odStats.gapTestMethod == GAP_LOCALORDERSTATS){
                    float sum_stats0 = vec_mean(r0) - vec_mean(gap);
                    float sum_stats1 = vec_mean(r1) - vec_mean(gap);
                    //gap.resize(0);
                    for(size_t j = 0; j < p4odStats.gapTestSkippedBandWidth; j++){
                        vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list0[j], CV_32F);
                        skipped.insert(skipped.end(), curr_vals.begin(), curr_vals.end());
                    }
                    for(size_t j = 0; j < p4odStats.gapTestSkippedBandWidth; j++){
                        vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list1[j], CV_32F);
                        skipped.insert(skipped.end(), curr_vals.begin(), curr_vals.end());
                    }
                    vector<float> bg(0);
                    orderStatsKSection(gap, bg, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p0 = zscore2pvalue((sum_stats0 - mu) / sigma);
                    p1 = zscore2pvalue((sum_stats1 - mu) / sigma);
                    //qInfo("p0:%.2f and p1:%.2f", p0, p1);
                }else if(p4odStats.gapTestMethod == GAP_ORDERSTATS){
                    float sum_stats0 = vec_mean(r0) - vec_mean(gap);
                    float sum_stats1 = vec_mean(r1) - vec_mean(gap);

                    orderStatsKSection(r0, gap, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p0 = zscore2pvalue((sum_stats0 - mu) / sigma);

                    orderStatsKSection(r1, gap, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p1 = zscore2pvalue((sum_stats1 - mu) / sigma);
                }else{
                    qFatal("non-defined method!");
                    p0 = 0;
                    p1 = 0;
                }

                if (p0 <= p_treshold && p1 <= p_treshold){ // gap is true
                    qInfo("!!!--->We find a real Gap!!!");
                    gap_tested_true[i] = 1;
                    real_gap_idx.insert(real_gap_idx.end(), gap_idx_list[i].begin(), gap_idx_list[i].end());
                }else{
                    gap_tested_true[i] = -1;
                }
            }
        }
        gap_idx_list.clear();
    }
    if (real_gap_idx.size() > 0){
        vector<vector<int>> groups;
        FOREACH_i(gap_tested_true){
            if(gap_tested_true[i] == -1){
                int target0 = i / n + 1;
                int target1 = i % n + 1;
                groups.push_back({target0, target1});
            }
        }
        mergeIntersectGroups(groups);
        vector<int> label_map (n);
        FOREACH_i(label_map){
            label_map[i] = i + 1;
        }
        FOREACH_i(groups){
            //vec_unique(groups[i]); //sorted also
            for(size_t j = 1; j < groups[i].size(); j++){
                assert(label_map[groups[i][j] - 1] == groups[i][j] || label_map[groups[i][j] - 1] == groups[i][0]);// "Duplicated Regions in two groups");
                label_map[groups[i][j] - 1] = groups[i][0];
            }
        }

        FOREACH_i_MAT(grown_seedMap){
            if(grown_seedMap.at<int>(i) > 0){
                grown_seedMap.at<int>(i) = label_map[grown_seedMap.at<int>(i) - 1];
            }
        }
        //ccShowSliceLabelMat(&grown_seedMap);
        vector<size_t> non_used;
        cellSegFromSynQuant.cell_num = rearrangeIdMap(&grown_seedMap, *cellSegFromSynQuant.idMap, non_used);

        FOREACH_i(real_gap_idx){ // remove the real gaps
            cellSegFromSynQuant.idMap->at<int>(real_gap_idx[i]) = 0;
        }
    }else{
        if(cellSegFromSynQuant.idMap->empty()){
            cellSegFromSynQuant.idMap->create(cellSegFromSynQuant.fgMap.dims,
                                              cellSegFromSynQuant.fgMap.size,
                                              CV_32S);
        }
        FOREACH_i_MAT(cellSegFromSynQuant.fgMap){
            if (cellSegFromSynQuant.fgMap.at<unsigned char>(i) > 0){
                cellSegFromSynQuant.idMap->at<int>(i) = 1;
            }
        }

        cellSegFromSynQuant.cell_num = 1;
    }
}
