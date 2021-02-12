#ifndef IMG_BASIC_PROC_H
#define IMG_BASIC_PROC_H

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <opencv2/imgproc.hpp> // image process
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/digamma.hpp>
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
#include <iomanip>      // std::setprecision
#include <set>
#include <iterator>
#include "maxflow_bk/graph.h" //max-flow BK algorithm

#include "cc3d.hpp" //connected component
#include "../celltracking/dt3d.hpp"

using namespace cv;
using namespace std;

//enum volDataTypes{UINT8C1 = 0, UINT8RGB = 1, UINT16C1 = 2, UINT16RGB = 3, UINT8RGBA = 4,  UINT16RGBA = 5};
enum filterDirection{DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z};
enum debiasMethods{TTEST2 = 0, TTEST2_VAR_KNOWN, NON_OV_TRUNCATED_NORMAL, OV_TRUNCATED_NORMAL, KSEC, APPROX3SEC};
enum cost_design_method{ARITHMETIC_AVERAGE = 1, GEOMETRIC_AVERAGE};
enum fgBoundaryHandleMethod {LEAVEALONEFIRST = 0, COMPETE, REPEAT};
enum gapTestMethods {GAP_TTEST = 0, GAP_ORDERSTATS, GAP_LOCALORDERSTATS};

enum trackingVarEstMethod {tve_INDEPENDENT = 0, tve_MEDIAN};
enum costCalMethod {ccm_CHI1SQUARE = 0, ccm_FISHER, ccm_ZSCORE};
enum priorDistrType {pri_GAUSS = 0, pri_WEIGHTEDAVG, pri_GAMMA};
enum slitMergeHandleMethods {NOJUMPALL = 0, NOJUMP, NOT_DECIDE};

enum segmentStableStatus {NOT_STABLE = -1, STABLE_TRACK_HEAD, STABLE_TRACK_MID, STABLE_TRACK_END};
enum parentsKidsConsistentFlag {CONSISTENCY_NOT_SURE = 0, PARENTS_KIDS_CONSISTENT, PARENTS_KIDS_NOT_CONSISTENT};
enum merge_split_decisions {SPLIT_BY_KIDS = 0, SPLIT_BY_PARENTS, MERGE_PARENTS, MERGE_KIDS, MERGE_BOTH_PARENTS_KIDS, MERGE_SPLIT_NOT_SURE};

#define REARRANGE_IDS true
#define NOT_REARRANGE_IDS false
#define FOREACH_i(x) for(size_t i = 0; i<x.size(); i++)
#define FOREACH_j(x) for(size_t j = 0; j<x.size(); j++)

#define FOREACH_i_MAT(x) for(size_t i=0; i<x.total(); i++)
#define FOREACH_j_MAT(x) for(size_t j=0; j<x.total(); j++)
#define FOREACH_i_ptrMAT(x) for(size_t i=0; i<x->total(); i++)
#define FOREACH_j_ptrMAT(x) for(size_t j=0; j<x->total(); j++)
#define FOREACH_ijk_ptrMAT(x) for(int i=0; i<x->size[0]; i++) \
    for(int j=0; j<x->size[1]; j++) \
    for(int k=0; k<x->size[2]; k++)

#define FOREACH_ijk_MAT(x) for(int i=0; i<x.size[0]; i++) \
    for(int j=0; j<x.size[1]; j++) \
    for(int k=0; k<x.size[2]; k++)

#define EPSILON_0 0.00001
const int n8_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
const int n8_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions

void principalCv2d(Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void principalCv3d(Mat* src3d, Mat &dst3d, float sigma[], int minIntensity);
void gaussianSmooth3Ddata(Mat &data4smooth, const float sigma[]);
void filterVolume(Mat* src3d, Mat &dst3d, Mat kernel, unsigned direction);
void filterZdirection(Mat* src3d, Mat &dst3d, Mat kernal_z);
float truncatedGauss(float mu, float sigma, float lower, float upper, float &t_mu, float &t_sigma);
float varByTruncate(vector<float> vals4var, int numSigma, int numIter);
float calVarianceStablization(Mat *src3d, Mat & varMap, vector<float> &varTrend, float validRatio = 0.95, int gap=2);
int connectedComponents3d(Mat* src3d, Mat &dst3d, int connect);
int floatMap2idMap(Mat* src3d, Mat &dst3d, int connect);
int rearrangeIdMap(Mat* src3d, Mat &dst3d, vector<size_t> &idMap);
void getRange(vector<int> idx_sub, int shift, int bound, Range &out_range);
void regionAvgIntensity(Mat* src3dFloatData, Mat* src3dIdMap, int numCC, vector<float> &avgIntensities);
void regionAvgIntensity(Mat* src3dFloatData, vector<vector<size_t>> &voxList, vector<float> &avgIntensities);

void extractVoxIdxList(Mat *label3d, vector<vector<size_t>> &voxList, int numCC, bool bk_extract = false);
void extractVoxIdxList(Mat *label3d, vector<vector<int>> &voxList, int numCC, bool bk_extract = false);
void extractVoxIdxList(vector<Mat> &label3d, vector<vector<size_t>> &voxList, vector<size_t> numCC);
void extractVoxIdxGivenId(Mat *label3d, vector<size_t> &voxList, int id);
void extractVolume(Mat *label3d, vector<size_t> &voxSzList, int numCC);
bool removeSmallCC(Mat &label3d, int &numCC, size_t min_size, bool relabel = true);
void volumeDilate(Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);
void volumeErode(Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type);

void volumeWrite(Mat *src3d, string filename);

void validSingleRegionExtract(Mat &vol3d, Mat *binary_mask, int connect);
int largestRegionIdExtract(Mat *label_map, int numCC, Mat *mask);
size_t fgMapSize(Mat *src3d, int datatype, float threshold_in = 0);
bool isempty(Mat src3d, int datatype, float threshold_in = 0);
bool isempty(Mat *src3d, int datatype, float threshold_in = 0);
vector<size_t> fgMapIdx(Mat *src3d, int datatype, float threshold_in);
vector<float> extractValsGivenMask(Mat *vol3d, int datatype, Mat *src3d, float threshold_in);
vector<float> extractValsGivenIdx(Mat *vol3d, vector<size_t> idx, int datatype);
double extractSumGivenIdx(Mat *vol3d, vector<size_t> idx, int datatype);
bool findUnrelatedCC(Mat *src3d4testing, int numCC, Mat *src3d4reference, Mat &dst3d);
bool findRelatedCC(Mat *src3d4testing, int numCC, Mat *src3d4reference, Mat &dst3d);
void neighbor_idx(vector<size_t> idx, vector<size_t> &center_idx, vector<size_t> &nei_idx, int sz[], int connect);
void regionGrow(Mat *label_map, int numCC, Mat &outLabelMap, Mat *scoreMap,
                Mat *fgMap, int connect, int cost_design[], bool bg2sink = true);
void setValMat(Mat &src, int datatype, vector<size_t> idx, float v);
float getMaxValMat(Mat &src, int datatype, vector<size_t> idx);
void scale_vol(Mat *src3d, int datatype, Mat *dst, float il, float ih, float vl=1.0, float vh=0.0);
void setValMat(Mat &src, int datatype, Mat *mask, float v);
void extractGapVoxel(Mat *label_map, Mat *fgMap, int numCC, int gap_radius, vector<vector<size_t>> &gap_voxIdx, vector<int> tested_flag);
void neighbor_idx_2d(vector<size_t> idx, Mat *fgMap, vector<vector<size_t>> &neighbor_idx_list, int radius);

bool inField( int r, int c, int z, int *sz );
bool inField( int r, int c, int *sz );
bool isOnBoundary2d(Mat *fgMap, size_t idx);
bool isOnBoundary2d(Mat *fgMap, int r, int c, int z);
void subVolExtract(Mat *src, int datatype, Mat &subVol, Range yxz_range[3]);
void subVolReplace(Mat &src, int datatype, Mat &subVol, Range yxz_range[3], int start = 0);
void subVolReplace(Mat &src, int datatype, Mat &subVol, float val, Range yxz_range[3]);

float distanceTransRegion2Region(bool *bw_ref_cell, vector<int> ref_range_xyz,
                                                       bool *bw_mov_cell, vector<int> mov_range_xyz,
                                                       vector<double> shift_xyz, vector<float> dist);
void adjacentRegions(Mat &src, vector<size_t> curr_label_idx, int curr_label, unordered_set<int> &adj_labels, int connect = 10);
bool adjacentRegions(Mat &src, vector<size_t> curr_label_idx, int test_adj_label, int connect = 10);
void coordinateTransfer(vector<size_t> &in_idx, int org_sz_yxz[3],
            vector<size_t> &out_idx, int crop_start_yxz[3], int crop_sz_yxz[3]);
void coordinateTransfer(vector<size_t> &in_idx, MatSize org_sz_yxz,
            vector<size_t> &out_idx, int crop_start_yxz[3], MatSize crop_sz_yxz);
void bwareaopenMat(Mat1b &fgMap, Mat1b &fgMap_new, size_t minSz, int connect);






// fucntions for display
void colorMapGen(Mat *src, Mat3b &colormap, String colorType = "HSV");
void ccShowSlice3Dmat(Mat *src3d, int datatype, int slice = 0 /*default 2d*/, bool binary = false);
void ccShowSliceLabelMat(Mat *src3d, int slice = 0 /*default 2d*/);
void ccShowSlice3Dmat(Mat src3d, int datatype, int slice = 0 /*default 2d*/, bool binary = false);
void ccShowSliceLabelMat(Mat src3d, int slice = 0 /*default 2d*/);
void label2rgb3d(Mat &src_label, Mat &src_intensity, Mat4b &dst);
void label2rgb2d(Mat1i &src, Mat3b &dst);


// Function to find t-test of
// two set of statistical data.

//template <class T> const T& max (const T& a, const T& b) {
//  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
//}

//template <class T> const T& min (const T& a, const T& b) {
//  return (a<b)?a:b;     // or: return comp(a,b)?b:a; for version (2)
//}

//template <typename T> void ind2sub(vector<T> idx, vector<int> &z, vector<int> &y, vector<int> &x, int *sz);
//template <typename T> void ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, int *sz);
//template <typename T> void ind2sub(T idx, int &z, int &y, int &x, int *sz);
//template <typename T> void ind2sub(T idx, int &y, int &x, int *sz);

//template <typename T> void sub2ind(vector<T> &idx, vector<int> z, vector<int> y, vector<int> x, int *sz);
//template <typename T> void sub2ind(vector<T> &idx, vector<int> y, vector<int> x, int *sz);
//template <typename T> void sub2ind(T &idx, int z, int y, int x, int *sz);
//template <typename T> void sub2ind(T &idx, int y, int x, int *sz);

template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, MatSize size);
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, MatSize size);
template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, int *size);
template <typename T> T vol_sub2ind(int y, int x, int z, int col_num, T page_sz);
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, int *size);
template <typename T> void vec_sub2ind(vector<T> &idx, vector<int> y, vector<int> x, vector<int> z, MatSize size);
template <typename T> void vec_ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, vector<int> &z, MatSize size);
template <typename T> void vec_sub2ind(vector<T> &idx, vector<int> y, vector<int> x, vector<int> z, int size[3]);
template <typename T> void vec_ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, vector<int> &z, int size[3]);

template <typename T> vector<double> vec_cumsum(vector<T> v1);
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
template <typename T> T zscore2pvalue(T z);
template <typename T> T pvalue2zscore(T p);
template <typename T> T normInv(T p, T mu = 0.0, T sigma = 1.0);
template <typename T> float binocdf(T num_cls0, T num_all, float p);

template <typename T> double chi2inv(T p, int df = 1); // inverse chi-square distribution
template <typename T> double gammacdf(T x, T a, T b, bool upper = true); // gamma cdf function
template <typename T> void gammafit(vector<T> data, T &a, T &b);// gamma fitting use MLE
template <typename T> void truncatedGammafit(vector<T> data, T &a, T &b); // not done
template <typename T> vector<T> vec_log(vector<T> data);

template <typename T> T vec_stddev(vector<T> const & func);
template <typename T> T vec_variance(vector<T> const & func);
template <typename T> float vec_mean(vector<T> const & func);
template <typename T> float vec_median(vector<T> const & func);
template <typename T> T vec_max(vector<T> const &func, size_t &max_val_idx);
template <typename T> T vec_min(vector<T> const &func, size_t &min_val_idx);
template <typename T> T vec_max(vector<T> const &func);
template <typename T> T vec_min(vector<T> const &func);
template <typename T> void vec_unique(vector<T> & v);
template <typename T> bool vec_find(vector<T> const &func, T target, size_t &idx);
template <typename T> bool set_exist(unordered_set<T> const &func, T target);
template <typename T> bool set_exist(unordered_set<T> const &set_sub, unordered_set<T> const &set_large);
template <typename T> bool set_equal(unordered_set<T> const &func0, unordered_set<T> const &func1);
template <typename T> float mat_mean(Mat *src3d, int datatype, vector<T> idx);
template <typename T> vector<size_t> sort_indexes(const vector<T> &v, bool ascending = true, size_t start_id = 0);
template <typename T> T ttest2(vector<T> arr1, vector<T> arr2);

template <typename T> T ttest2_var_known(vector<T> arr1, vector<T> arr2, T var_known);

template <typename T> void nonOV_truncatedGauss(size_t M, size_t N, T &mu, T &sigma);
template <typename T> void OV_truncatedGauss(vector<T> fg, vector<T> bg, T &mu, T &sigma);
template <typename T> void orderStatsKSection(vector<T> fg, vector<T> bg, vector<T> otherVals, float &mu, float &sigma);
template <typename T> size_t overlap_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in = 0);
template <typename T> bool isempty_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in = 0);

template <typename T> void vec_unique(vector<T> & v);
template <typename T> void vec_ele_wise_abs_diff(vector<T> & v1, vector<T> & v2);

template <typename T> void mergeIntersectGroups(vector<vector<T>> &groups);


template <class Iter> typename std::iterator_traits<Iter>::value_type Mode(Iter first, Iter last);
#endif // IMG_BASIC_PROC_H


