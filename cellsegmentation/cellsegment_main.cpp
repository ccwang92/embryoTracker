#include "cellsegment_main.h"
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
using namespace cv;
cellSegmentMain::cellSegmentMain(unsigned char *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
{
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

    long sz_single_frame_char = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];

    if (data_type == V3D_UINT16) sz_single_frame_char *= 2;
    else if(data_type != V3D_UINT8) qFatal("Unsupported data type\n");

    for (int i = 0; i < time_points; i++){
        cellSegmentSingleFrame(data_grayim4d + sz_single_frame_char*i);
    }
}
void cellSegmentMain::cellSegmentSingleFrame(unsigned char *data_grayim3d)
{
    // The data is mirrored by the y-direction, this is because OpenGL
    // texture are loaded left to right, bottom to top. Most image loaders
    // however will store the image in memory left to right, top to bottom.
    Mat *input_3dmat;
    if (data_type == V3D_UINT8){
        input_3dmat = new Mat(3, data_rows_cols_slices, CV_8UC1, data_grayim3d);
//        for (int s = 0; s < data_rows_cols_slices[2]; s ++)
//        {
//            long offset = data_rows_cols_slices[0]*data_rows_cols_slices[1]*s;
//            input_3dmat = new Mat(data_rows_cols_slices[0], data_rows_cols_slices[1], CV_8UC1, (void *)(data_grayim3d+offset));
//            imshow(std::to_string(s), *input_3dmat);
//            waitKey(0);
//            destroyWindow(std::to_string(s));
//        }
    }else if(data_type == V3D_UINT16){
        input_3dmat = new Mat(3, data_rows_cols_slices, CV_16UC1, data_grayim3d);
//        for (int s = 0; s < data_rows_cols_slices[2]; s ++)
//        {
//            long offset = data_rows_cols_slices[0]*data_rows_cols_slices[1]*s * 2; // 16 bit = 2 * char
//            input_3dmat = new Mat(data_rows_cols_slices[0], data_rows_cols_slices[1], CV_16UC1, (void *)(data_grayim3d+offset));
//            *input_3dmat *= 10;
//            imshow(std::to_string(s), *input_3dmat);
//            waitKey(0);
//            destroyWindow(std::to_string(s));
//        }
    }else{
        //qDebug("unsupported data type!\n");
        MESSAGE_ASSERT("unsupported data type!\n");
        return;
    }
    Mat inputVolFloat;
    input_3dmat->convertTo(inputVolFloat, CV_32F); // data we will work on in the future

}
