#include "cellsegment_main.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>

#include "synquant_simple.h"
cellSegmentMain::cellSegmentMain(unsigned char *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
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
    principalCurv2d.resize(time_points);
    principalCurv3d.resize(time_points);
    varMaps.resize(time_points);
    varTrends.resize(time_points);
    variances.resize(time_points);
    number_cells.resize(time_points);
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];

    Mat *data4d;
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, (int *)bufSize, CV_8UC1, data_grayim4d);
    }
    else if(data_type == V3D_UINT8){
        data4d = new Mat(4, (int *)bufSize, CV_16UC1, data_grayim4d);
    }else{
        qFatal("Unsupported data type\n");
    }
    normalize(*data4d, *data4d, 0, 255, NORM_MINMAX, CV_8U);

    for (int i = 0; i < time_points; i++){
        unsigned char *ind = (unsigned char*)data4d->data + sz_single_frame*i; // sub-matrix pointer
        Mat *single_frame = new Mat(3, data4d->size, CV_8U, ind);
        cell_label_maps[i].create(3, data4d->size, CV_32S); // int label
        cellSegmentSingleFrame(single_frame , i);
    }
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
    /******** start to do cell segmentation *******/
    float sigma2d[3] = {3.0, 3.0, 0.0};
    principalCv2d(dataVolFloat, principalCurv2d[curr_frame], sigma2d, p4segVol.min_intensity);
    float sigma3d[3] = {5.0, 5.0, 1.0};
    principalCv2d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, p4segVol.min_intensity);

    variances[curr_frame] = calVarianceStablization(data_grayim3d, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);

    // first use synQuant to get 1-tier seed regions
    synQuantSimple seed_from_synQuant(data_grayim3d, variances[curr_frame], p4segVol, p4odStats);
    // second refine the seed regions

    // third get 2-tier seed regions
    //seed_from_synQuant.retrieve_seeds();
    // fourth refine 2-tier seed regions

    // Finally, return
}
/**
 * @brief cellSegmentMain::regionWiseAnalysis4d, refine the region based on the results from synQuant
 * @param idMap
 * @param eigMap
 * @param vid
 * @param varMap
 * @param test_ids
 */
void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
                                           Mat *eigMap3d, Mat *varMap, vector<int> test_ids){
    //1. sort the seeds based on intensity levels
    vector<float> seed_intensity(seed_num);
    regionAvgIntensity(dataVolFloat, idMap, seed_intensity);
    vector<size_t> seed_intensity_order;
    seed_intensity_order = sort_indexes(seed_intensity, false); // false->descending
    //2. for each seed, refine it region
    vector<vector<size_t>> voxIdxList(seed_num);
    extractVoxIdxList(idMap, voxIdxList, seed_num);
    FOREACH_i(seed_intensity_order){
        int seed_id = seed_intensity_order[i];
        singleCellSeed seed;
        cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, idMap, eigMap2d,
                        eigMap3d, varMap, seed, p4segVol);
        refineSeed2Region(seed, p4odStats, p4segVol);

    }
}

void cellSegmentMain::cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d, Mat *idMap, Mat *eigMap2d,
                               Mat *eigMap3d, Mat *varMap, singleCellSeed &seed, segParameter p4segVol){
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

    seed.gap2dMap = seed.eigMap2d > 0;
    seed.gap3dMap = seed.eigMap3d > 0;
    seed.volUint8 = (*data_grayim3d)(seed.crop_range_yxz); // all shallow copy
    seed.idMap = (*idMap)(seed.crop_range_yxz); // all shallow copy

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

void cellSegmentMain::refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol){

    // 1st: foreground Detection reusing synQuant
    synQuantSimple synQuant_refine_seed(seed, p4segVol, p4odStats);
}
