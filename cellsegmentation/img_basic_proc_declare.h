#ifndef IMG_BASIC_PROC_H
#define IMG_BASIC_PROC_H

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

enum filterDirection{DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z};
enum debiasMethods{TTEST2 = 0, TTEST2_VAR_KNOWN, NON_OV_TRUNCATED_NORMAL, OV_TRUNCATED_NORMAL, KSEC, APPROX3SEC};

#define FOREACH_i(x) for(size_t i = 0; i<x.size(); i++)

#define FOREACH_i_MAT(x) for(size_t i=0; i<x.total(); i++)
#define FOREACH_i_ptrMAT(x) for(size_t i=0; i<x->total(); i++)
#define FOREACH_ijk_ptrMAT(x) for(int i=0; i<x->size[0]; i++) \
    for(int j=0; j<x->size[1]; j++) \
    for(int k=0; k<x->size[2]; k++)

#define FOREACH_ijk_MAT(x) for(int i=0; i<x.size[0]; i++) \
    for(int j=0; j<x.size[1]; j++) \
    for(int k=0; k<x.size[2]; k++)


void principalCv2d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void principalCv3d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void gaussianSmooth3Ddata(Mat &data4smooth, const float sigma[]);
void filterVolume(const Mat* src3d, Mat &dst3d, Mat kernel, unsigned direction);
void filterZdirection(const Mat* src3d, Mat &dst3d, Mat kernal_z);
float truncatedGauss(float mu, float sigma, float lower, float upper, float &t_mu, float &t_sigma);
float varByTruncate(vector<float> vals4var, int numSigma, int numIter);
float calVarianceStablization(const Mat *src3d, Mat & varMap, vector<float> &varTrend, float validRatio, int gap);
int connectedComponents3d(const Mat* src3d, Mat &dst3d, int connect);
int floatMap2idMap(Mat* src3d, Mat &dst3d, int connect);
int rearrangeIdMap(Mat* src3d, Mat &dst3d, vector<size_t> &idMap);
void getRange(vector<int> idx_sub, int shift, int bound, Range &out_range);
void regionAvgIntensity(Mat* src3dFloatData, Mat* src3dIdMap, vector<float> &avgIntensities);

void extractVoxIdxList(const Mat *label3d, vector<vector<size_t>> &voxList, int numCC);
void removeSmallCC(Mat &label3d, int &numCC, size_t min_size, bool relabel);
void volumeDilate(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);
void volumeErode(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);

void volumeWrite(Mat *src3d, string filename);

void singleRegionCheck(Mat &vol3d, Mat *binary_mask, int connect);
int largestRegionIdExtract(Mat *label_map, int numCC, Mat *mask);
size_t fgMapSize(Mat *src3d, int datatype, float threshold_in = 0);
bool isempty(Mat *src3d, int datatype, float threshold_in = 0);
vector<size_t> fgMapIdx(Mat *src3d, int datatype, float threshold_in);
vector<float> fgMapVals(Mat *vol3d, Mat *src3d, int datatype, float threshold_in);
bool findUnrelatedCC(Mat *src3d4testing, int numCC, Mat *src3d4reference, Mat &dst3d);
// Function to find t-test of
// two set of statistical data.

//template <class T> const T& max (const T& a, const T& b) {
//  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
//}

//template <class T> const T& min (const T& a, const T& b) {
//  return (a<b)?a:b;     // or: return comp(a,b)?b:a; for version (2)
//}

template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, MatSize size);
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, MatSize size);
template <typename T> void vec_sub2ind(vector<T> &idx, vector<int> y, vector<int> x, vector<int> z, MatSize size);
template <typename T> void vec_ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, vector<int> &z, MatSize size);


template <typename T> vector<T> vec_cumsum(vector<T> v1);
template <typename T> vector<T> vec_pointMultiply(vector<T> v1, vector<T> v2);
template <typename T> vector<T> vec_Minus(vector<T> v1, vector<T> v2);
template <typename T> vector<T> vec_Minus(vector<T> v1, T s2);
template <typename T> vector<T> vec_Add(vector<T> v1, vector<T> v2);
template <typename T> vector<T> vec_pointDivide(vector<T> v1, vector<T> v2);
template <typename T> vector<T> vec_smallerthan(vector<T> values, T threshold, bool strict = true);
template <typename T> vector<T> vec_largerthan(vector<T> values, T threshold, bool strict = true);
template <typename T> vector<T> vec_atrange(vector<T> values, T ub, T lb, bool strict = true);
template <typename T> T normalCDF(T x, T m = 0, T s = 1);
template <typename T> T normalPDF(T x, T m = 0, T s = 1);
template <typename T> T zscore2palue(T z);
template <typename T> T pvalue2zscore(T p);
template <typename T> T normInv(T p, T mu = 0.0, T sigma = 1.0);
template <typename T> T vec_stddev(vector<T> const & func);
template <typename T> T vec_variance(vector<T> const & func);
template <typename T> T vec_mean(vector<T> const & func);
template <typename T> T vec_max(vector<T> const &func, size_t &max_val_idx);
template <typename T> T vec_min(vector<T> const &func, size_t &min_val_idx);
template <typename T> float mat_mean(Mat *src3d, int datatype, vector<T> idx);
template <typename T> vector<size_t> sort_indexes(const vector<T> &v, bool ascending = true, size_t start_id = 0);
template <typename T> T ttest2(vector<T> arr1, vector<T> arr2);

template <typename T> T ttest2_var_known(vector<T> arr1, vector<T> arr2, T var_known);

template <typename T> void nonOV_truncatedGauss(size_t M, size_t N, T &mu, T &sigma);
template <typename T> void OV_truncatedGauss(vector<T> fg, vector<T> bg, T &mu, T &sigma);
template <typename T> void orderStatsKSection(vector<T> fg, vector<T> bg, vector<T> otherVals, float &mu, float &sigma);
template <typename T> size_t overlap_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in = 0);
template <typename T> bool isempty_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in = 0);

template <typename T> void scale_vol(Mat *src3d, int datatype, Mat *dst, float il, float ih, float vl=1.0, float vh=0.0);
#endif // IMG_BASIC_PROC_H


