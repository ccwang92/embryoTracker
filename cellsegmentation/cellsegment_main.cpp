#include "cellsegment_main.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include "img_basic_proc_declare.h"

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
                                                   p4odStats.varAtRatio, p4odStats.gap4fgbgCompare);


}
