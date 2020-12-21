#ifndef IMG_BASIC_PROC_DECLARE_H
#define IMG_BASIC_PROC_DECLARE_H

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <opencv2/imgproc.hpp> // image process
//#include <itkImage.h>
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <vector>

using namespace cv;
using namespace std;
void principalCv2d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity);
void principalCv3d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity);
void gaussianSmooth3Ddata(Mat &data4smooth, const float *sigma);
void filterZdirection(Mat* src3d, Mat &dst3d, Mat kernal_z);
int ConnectedComponents3d(Mat* src3d, Mat &dst3d, int connect);

template <typename T> vector<T> vec_smallerthan(vector<T> values, T threshold, bool strict);
template <typename T> vector<T> vec_largerthan(vector<T> values, T threshold, bool strict);
template <typename T> vector<T> vec_atrange(vector<T> values, T ub, T lb, bool strict);
template <typename T> T normalCDF(T x, T m = 0, T s = 1);
template <typename T> T normalPDF(T x, T m = 0, T s = 1);
template <typename T> T vec_stddev(vector<T> const & func);
template <typename T> T vec_variance(vector<T> const & func);
template <typename T> T vec_mean(vector<T> const & func);
#endif // IMG_BASIC_PROC_DECLARE_H
