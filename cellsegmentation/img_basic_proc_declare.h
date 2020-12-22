#ifndef IMG_BASIC_PROC_DECLARE_H
#define IMG_BASIC_PROC_DECLARE_H

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <opencv2/imgproc.hpp> // image process
//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//#include "itkFlatStructuringElement.h"
//#include "itkBinaryDilateImageFilter.h"
//#include "itkBinaryBallStructuringElement.h"
//#include "itkOpenCVImageBridge.h"

#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace cv;
using namespace std;
void principalCv2d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void principalCv3d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void gaussianSmooth3Ddata(Mat &data4smooth, const float sigma[]);
void filterVolume(const Mat* src3d, Mat &dst3d, Mat kernel, unsigned direction);
void filterZdirection(const Mat* src3d, Mat &dst3d, Mat kernal_z);
float truncatedGauss(float mu, float sigma, float lower, float upper, float &t_mu, float &t_sigma);
float varByTruncate(vector<float> vals4var, int numSigma, int numIter);
float calVarianceStablization(const Mat* src3d, Mat & varMap, vector<float> &varTrend, float validRatio, int gap);

template <typename T> vector<T> vec_smallerthan(vector<T> values, T threshold, bool strict = true);
template <typename T> vector<T> vec_largerthan(vector<T> values, T threshold, bool strict = true);
template <typename T> vector<T> vec_atrange(vector<T> values, T ub, T lb, bool strict = true);
template <typename T> T normalCDF(T x, T m = 0, T s = 1);
template <typename T> T zscore2Pvalue(T z);
template <typename T> T normalPDF(T x, T m = 0, T s = 1);
template <typename T> T vec_stddev(vector<T> const & func);
template <typename T> T vec_variance(vector<T> const & func);
template <typename T> T vec_mean(vector<T> const & func);

int connectedComponents3d(const Mat* src3d, Mat &dst3d, int connect);
void extractVoxIdxList(const Mat *label3d, vector<vector<size_t>> &voxList, int numCC);
void removeSmallCC(Mat &label3d, int &numCC, size_t min_size, bool relabel);
void volumeDilate(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);
void volumeErode(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);

void volumeWrite(Mat *src3d, string filename);
//template <class T> const T& max (const T& a, const T& b) {
//  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
//}

//template <class T> const T& min (const T& a, const T& b) {
//  return (a<b)?a:b;     // or: return comp(a,b)?b:a; for version (2)
//}
#endif // IMG_BASIC_PROC_DECLARE_H
