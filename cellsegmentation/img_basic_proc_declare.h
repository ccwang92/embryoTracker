#ifndef IMG_BASIC_PROC_DECLARE_H
#define IMG_BASIC_PROC_DECLARE_H

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <opencv2/imgproc.hpp> // image process
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>

void principalCv2d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity = 0);
void smooth3Ddata(Mat &data4smooth, const float *sigma);
#endif // IMG_BASIC_PROC_DECLARE_H
