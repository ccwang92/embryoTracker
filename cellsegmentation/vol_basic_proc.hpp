#ifndef VOL_BASIC_PROC_HPP
#define VOL_BASIC_PROC_HPP
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <opencv2/imgproc.hpp> // image process
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iomanip>      // std::setprecision

#include "maxflow_bk/graph.h" //max-flow BK algorithm

#include "cc3d.hpp" //connected component

using namespace cv;
using namespace std;

enum filterDirection{DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z};
enum debiasMethods{TTEST2 = 0, TTEST2_VAR_KNOWN, NON_OV_TRUNCATED_NORMAL, OV_TRUNCATED_NORMAL, KSEC, APPROX3SEC};
enum cost_design_method{ARITHMETIC_AVERAGE = 1, GEOMETRIC_AVERAGE};
enum fgBoundaryHandleMethod {LEAVEALONEFIRST = 0, COMPETE, REPEAT};
enum gapTestMethods {GAP_TTEST = 0, GAP_ORDERSTATS, GAP_LOCALORDERSTATS};

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
namespace volproc {
//double normalCDF(double value) //standard
//{
//   return 0.5 * erfc(-value * M_SQRT1_2);
//}
/** template functions needs to be implemented in header files. */
template <typename T> vector<T> vec_cumsum(vector<T> v1){
    assert(v1.size() > 0);
    vector<T> out(v1.size());
    out[0] = v1[0];
    for (size_t i=1 ;i<v1.size(); i++){
        out[i] = v1[i] + out[i-1];
    }
    return out;
}
template <typename T> vector<T> vec_pointMultiply(vector<T> v1, vector<T> v2){
    assert(v1.size() == v2.size());
    vector<T> out(v1.size());
    for (size_t i=0 ;i<v1.size(); i++){
        out[i] = v1[i] * v2[i];
    }
    return out;
}
template <typename T> vector<T> vec_Minus(vector<T> v1, vector<T> v2){
    assert(v1.size() == v2.size());
    vector<T> out(v1.size());
    for (size_t i=0 ;i<v1.size(); i++){
        out[i] = v1[i] - v2[i];
    }
    return out;
}
template <typename T> vector<T> vec_Minus(vector<T> v1, T s2){
    vector<T> out(v1.size());
    for (size_t i=0 ;i<v1.size(); i++){
        out[i] = v1[i] - s2;
    }
    return out;
}
template <typename T> vector<T> vec_Add(vector<T> v1, vector<T> v2){
    assert(v1.size() == v2.size());
    vector<T> out(v1.size());
    for (size_t i=0 ;i<v1.size(); i++){
        out[i] = v1[i] + v2[i];
    }
    return out;
}
template <typename T> vector<T> vec_pointDivide(vector<T> v1, vector<T> v2){
    assert(v1.size() == v2.size());
    vector<T> out(v1.size());
    for (size_t i=0 ;i<v1.size(); i++){
        out[i] = v1[i] / v2[i];
    }
    return out;
}
template <typename T> vector<T> vec_smallerthan(vector<T> values, T threshold, bool strict){
    vector<T> out;
    if (strict){
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] < threshold)
                out.push_back(values[i]);
        }
    }else{
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] <= threshold)
                out.push_back(values[i]);
        }
    }
    return out;
}
template <typename T> vector<T> vec_largerthan(vector<T> values, T threshold, bool strict){
    vector<T> out;
    if (strict){
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] > threshold)
                out.push_back(values[i]);
        }
    }else{
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] >= threshold)
                out.push_back(values[i]);
        }
    }
    return out;
}
template <typename T> vector<T> vec_atrange(vector<T> values, T ub, T lb, bool strict){
    assert(ub>lb);
    vector<T> out;
    if (strict){
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] > lb && values[i] < ub)
                out.push_back(values[i]);
        }
    }else{
        for (size_t i = 0; i < values.size(); i++){
            if(values[i] >= lb && values[i] <= ub)
                out.push_back(values[i]);
        }
    }
    return out;
}
template <typename T> T normalCDF(T x, T m = 0, T s = 1)
{
    T a = (x - m) / s;
    return 0.5 * erfc(-a * M_SQRT1_2); //erfc(x) = 1-erf(x), erf(-x) = -erf(x)
}
template <typename T> T normalPDF(T x, T m = 0, T s = 1)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-T(0.5) * a * a);
}
template <typename T> T zscore2pvalue(T z){
    return 1-normalCDF(z);
}
template <typename T> T pvalue2zscore(T p){
    const double a1 = -39.6968302866538;
    const double a2 = 220.946098424521;
    const double a3 = -275.928510446969;
    const double a4 = 138.357751867269;
    const double a5 = -30.6647980661472;
    const double a6 = 2.50662827745924;

    const double b1 = -54.4760987982241;
    const double b2 = 161.585836858041;
    const double b3 = -155.698979859887;
    const double b4 = 66.8013118877197;
    const double b5 = -13.2806815528857;

    const double c1 = -7.78489400243029E-03;
    const double c2 = -0.322396458041136;
    const double c3 = -2.40075827716184;
    const double c4 = -2.54973253934373;
    const double c5 = 4.37466414146497;
    const double c6 = 2.93816398269878;

    const double d1 = 7.78469570904146E-03;
    const double d2 = 0.32246712907004;
    const double d3 = 2.445134137143;
    const double d4 = 3.75440866190742;

    //Define break-points
    // using Epsilon is wrong; see link above for reference to 0.02425 value
    //const double pLow = double.Epsilon;
    const double pLow = 0.02425;

    const double pHigh = 1 - pLow;

    //Define work variables
    double q;
    double result = 0;

    // if argument out of bounds.
    // set it to a value within desired precision.
    if (p <= 0)
        p = pLow;

    if (p >= 1)
        p = pHigh;

    if (p < pLow)
    {
        //Rational approximation for lower region
        q = sqrt(-2 * log(p));
        result = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    }
    else if (p <= pHigh)
    {
        //Rational approximation for lower region
        q = p - 0.5;
        double r = q * q;
        result = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    }
    else if (p < 1)
    {
        //Rational approximation for upper region
        q = sqrt(-2 * log(1 - p));
        result = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    }

    return result;
}
template <typename T> T normInv(T p, T mu = 0, T sigma = 1)
{
    // original from https://gist.github.com/kmpm/1211922/6b7fcd0155b23c3dc71e6f4969f2c48785371292
    assert(p >= 0 && p <= 1);
    assert(sigma >= 0);

    if (p == 0)
    {
        return -MAXFLOAT;
    }
    if (p == 1)
    {
        return MAXFLOAT;
    }
    if (sigma == 0)
    {
        return mu;
    }

    T q, r, val;

    q = p - 0.5;

    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
            Produces the normal deviate Z corresponding to a given lower
            tail area of P; Z is accurate to about 1 part in 10**16.
    */
    if (abs(q) <= .425)
    {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
               q * (((((((r * 2509.0809287301226727 +
                          33430.575583588128105) * r + 67265.770927008700853) * r +
                        45921.953931549871457) * r + 13731.693765509461125) * r +
                      1971.5909503065514427) * r + 133.14166789178437745) * r +
                    3.387132872796366608)
               / (((((((r * 5226.495278852854561 +
                        28729.085735721942674) * r + 39307.89580009271061) * r +
                      21213.794301586595867) * r + 5394.1960214247511077) * r +
                    687.1870074920579083) * r + 42.313330701600911252) * r + 1);
    }
    else
    { /* closer than 0.075 from {0,1} boundary */

        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = 1 - p;
        else
            r = p;

        r = sqrt(-log(r));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5)
        { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1);
        }
        else
        { /* very close to  0 or 1 */
            r += -5;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1);
        }

        if (q < 0.0)
        {
            val = -val;
        }
    }

    return mu + sigma * val;
}
template <typename T> T vec_mean(vector<T> const & func)
{
    return accumulate(func.begin(), func.end(), 0.0) / func.size();
}
template <typename T> T vec_variance(vector<T> const & func)
{
    T mean = vec_mean(func);
    T sq_sum = inner_product(func.begin(), func.end(), func.begin(), 0.0,
        [](T const & x, T const & y) { return x + y; },
        [mean](T const & x, T const & y) { return (x - mean)*(y - mean); });
    return sq_sum / ( func.size() - 1 );
}

template <typename T> T vec_stddev(vector<T> const & func)
{
    return sqrt(vec_variance(func));
}


template <typename T> T vec_max(vector<T> const &func, size_t &max_val_idx){
    auto max_val_it = max_element(func.begin(), func.end());
    max_val_idx = distance(func.begin(), max_val_it);
    return *max_val_it;
}
template <typename T> T vec_min(vector<T> const &func, size_t &min_val_idx){
    auto min_val_it = min_element(func.begin(), func.end());
    min_val_idx = distance(func.begin(), min_val_it);
    return *min_val_it;
}

template <typename T> float mat_mean(Mat *src3d, int datatype, vector<T> idx){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    double sum = 0.0;
    if (datatype == CV_8U){
        FOREACH_i(idx){
            sum += src3d->at<unsigned char>(idx[i]);
        }
        return (float) sum/idx.size();
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            sum += src3d->at<float>(idx[i]);
        }
        return (float) sum/idx.size();
    }else{ //CV_32S is the same as float
        FOREACH_i(idx){
            sum += src3d->at<int>(idx[i]); //int_32_t as default
        }
        return (float) sum/idx.size();
    }
}
template <typename T> vector<size_t> sort_indexes(const vector<T> &v, bool ascending = true, size_t start_id = 0) {
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), start_id);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    if (ascending)
        stable_sort(idx.begin(), idx.end(),
                    [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    else
        stable_sort(idx.begin(), idx.end(),
                    [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}

template <typename T> T ttest2(vector<T> arr1, vector<T> arr2)
{
    size_t n = arr1.size();
    size_t m = arr2.size();
    T mean1 = vec_mean(arr1);
    T mean2 = vec_mean(arr2);
    T var1 = vec_variance(arr1);
    T var2 = vec_variance(arr2);
    // !! t-stats not z-score
    return (mean1 - mean2) / sqrt((1.0/n+1.0/m) * (((n-1)*var1 + (m-1)*var2)/(n+m-2)));
}

template <typename T> T ttest2_var_known(vector<T> fg, vector<T> bg, T var_known){
    T sum_st = vec_mean(fg) - vec_mean(bg);
    double n = (double)fg.size();
    double m = (double)bg.size();
    T sigma = sqrt(var_known*(n+m)/(n*m));
    return  sum_st / sigma;
}
float truncatedGauss(float mu, float sigma, float lower, float upper, float &t_mu, float &t_sigma){
    float alpha = (lower-mu) / sigma;
    float beta = (upper-mu) / sigma;

    float Z = normalCDF(beta) - normalCDF(alpha);

    float a_pdf = normalPDF(alpha);
    float b_pdf = normalPDF(beta);

    t_mu = mu + (a_pdf-b_pdf)*sigma/Z;
    if (abs(alpha) > 1e7) alpha = -1e7;
    if (abs(beta) > 1e7) beta = 1e7;

    float t1 = (alpha*a_pdf-beta*b_pdf) / Z;
    float t2 = pow(((a_pdf - b_pdf) / Z), 2);
    float varCorrectionTerm = 1+t1-t2;
    t_sigma = sigma*sqrt(varCorrectionTerm);

    return varCorrectionTerm;
}

template <typename T> void nonOV_truncatedGauss(size_t M, size_t N, T &mu, T &sigma){
        T lower_fg = normInv(1-(double)M/(M+N));
        T f_mu, f_sigma;
        truncatedGauss(0, 1, lower_fg, INFINITY, f_mu, f_sigma);
        f_sigma = f_sigma/sqrt((double)M);

        T upper_nei = lower_fg;
        T n_mu, n_sigma;
        truncatedGauss(0, 1, -INFINITY, upper_nei, n_mu, n_sigma);
        n_sigma = n_sigma/sqrt((double)M);

        mu = f_mu - n_mu;
        sigma = sqrt(f_sigma*f_sigma + n_sigma*n_sigma);

//        T sum_st = mean_fg - mean_bg;
//        sigma = sigma * sqrt(var_known);
//        T z_debias = mu*sqrt(var_known); // mu
//        return (sum_st - z_debias) / sigma; //zscore
}

template <typename T> void OV_truncatedGauss(vector<T> fg, vector<T> bg, T &mu, T &sigma){
    double M = (double)fg.size(), N = (double)bg.size();

    fg.insert(fg.end(), bg.begin(), bg.end());
    vector<size_t> sorted_id = sort_indexes(fg, false, 1);

    T fg_ratio = 0, bg_ratio;
    for (size_t i = 0; i<sorted_id.size(); i++){
        if(sorted_id[i] <= M){
            fg_ratio += i;
        }else{
            bg_ratio += i;
        }
    }
    fg_ratio /= ((M+N)*M);
    bg_ratio /= ((M+N)*N);


    T lower_fg = -norminv(2*fg_ratio);
    T f_mu, f_sigma;
    truncatedGauss(0, 1, lower_fg, INFINITY, f_mu, f_sigma);
    f_sigma = f_sigma/sqrt(M);

    T nei_ratio = 2*(1-bg_ratio);
    T upper_nei = norminv(nei_ratio);
    T n_mu, n_sigma;
    truncatedGauss(0, 1, -INFINITY, upper_nei, n_mu, n_sigma);
    n_sigma = n_sigma/sqrt(N);

    mu = f_mu - n_mu;
    sigma = sqrt(f_sigma^2 + n_sigma^2);

//    T sum_st = vec_mean(fg) - vec_mean(bg);
//    sigma = sigma * sqrt(var_known);
//    T z_debias = mu*sqrt(var_known); // mu

//    return (sum_st - z_debias) / sigma; //zscore
}

/**
 * @brief orderStatsKSection
 * @param midVals: values not used but together with fg&bg to form a intact Gaussian
 */
template <typename T> void orderStatsKSection_f1(T x, T &y, T &xnormcdf, T &xnormpdf){
    xnormcdf = normalCDF(x);
    xnormpdf = normalPDF(x);
    y = x*xnormcdf+xnormpdf;
}
template <typename T> T orderStatsKSection_f2(T x, T xnormcdf, T xnormpdf){
    return 0.5*(xnormcdf*x*x - xnormcdf + xnormpdf*x);
}
template <typename T> void orderStatsKSection(vector<T> fg, vector<T> bg, vector<T> otherVals, float &mu, float &sigma){
    double M = (double)fg.size();
    double N = (double)bg.size();
    double mid_sz = (double)otherVals.size();
    double n = M+N+mid_sz;

    bg.insert(bg.end(), fg.begin(), fg.end());
    bg.insert(bg.end(), otherVals.begin(), otherVals.end());// fg: 1-M, bg, M+1-M+N, mid: M+N+1:n
    vector<size_t> sorted_id = sort_indexes(bg, false, 1);
    vector<byte> sorted_class_id((long)n);
    for (size_t i = 0; i<sorted_id.size(); i++){
        if(sorted_id[i] <= M){
            sorted_class_id[i] = (byte)-1;
        }else if(sorted_id[i] <= M+N){
            sorted_class_id[i] = (byte)1;
        }
        else{
            sorted_class_id[i] = (byte)0;
        }
    }
    vector<size_t> bkpts;
    for (size_t i = 1; i < n ; i ++){
        if (sorted_class_id[i] != sorted_class_id[i-1]){
            bkpts.push_back(i-1);
        }
    }
    //bkpts = find(labels(2:end)-labels(1:end-1));
    vector<float> ai(bkpts.size() + 1);
    //ai = cat(1, labels(bkpts), labels(end));
    for (size_t i=0; i<bkpts.size(); i++){
        if (sorted_class_id[bkpts[i]] == (byte)-1){
            ai[i] = -1.0*n/N;
        }else if(sorted_class_id[bkpts[i]] == (byte)1){
            ai[i] = 1.0*n/M;
        }else{
            ai[i] = 0.0;
        }
    }
    if (sorted_class_id[sorted_class_id.size()-1] == (byte)-1){
        ai[ai.size()-1] = -1.0*n/N;
    }else if(sorted_class_id[sorted_class_id.size()-1] == (byte)1){
        ai[ai.size()-1] = 1.0*n/M;
    }else{
        ai[ai.size()-1] = 0.0;
    }
    float delta = 1.0/n;
    // bi is start, ti is end of the i-th section
    vector<float> bi (ai.size()); bi[0] = 0.0;
    vector<float> ti (ai.size()); ti[ti.size()-1] = 1.0;
    for(size_t i = 0; i< bi.size()-1; i++){
        bi[i+1] = (bkpts[i] + 1) * delta;
        ti[i] = bi[i+1];
    }

    vector<float> Finvbi(ai.size()), Finvti(ai.size());
    for(size_t i=1;i<ai.size();i++){
        Finvbi[i] = normInv(bi[i]);
        Finvti[i-1] = normInv(ti[i-1]);
    }
    Finvbi[0] = -100000.0; // set 0.0 -> -inf
    Finvti[ai.size() - 1] = 100000.0; // set 1.0 -> inf

    mu = 0.0;
    FOREACH_i(ai){
        mu += ai[i] * (normalPDF(Finvbi[i]) - normalPDF(Finvti[i]));
    }


    // mat implementaion
    vector<float> f1Finvti(Finvti.size()), FinvtiNcdf(Finvti.size()), FinvtiNpdf(Finvti.size());
    FOREACH_i(Finvti){
        orderStatsKSection_f1(Finvti[i], f1Finvti[i], FinvtiNcdf[i], FinvtiNpdf[i]);
    }
    vector<float> f1Finvbi(Finvbi.size()), FinvbiNcdf(Finvbi.size()), FinvbiNpdf(Finvbi.size());
    FOREACH_i(Finvti){
        orderStatsKSection_f1(Finvbi[i], f1Finvbi[i], FinvbiNcdf[i], FinvbiNpdf[i]);
    }
    vector<float> f1Finvti_f1Finvbi = vec_Minus(f1Finvti, f1Finvbi);
    vector<float> aixFinvtj_Finvbj = vec_pointMultiply(ai, vec_Minus(Finvti, Finvbi));
    vector<float> cumsum_aixFinvtj_Finvbj = vec_cumsum(aixFinvtj_Finvbj);
    float all_sum = cumsum_aixFinvtj_Finvbj[cumsum_aixFinvtj_Finvbj.size() - 1];
    FOREACH_i(cumsum_aixFinvtj_Finvbj){
        cumsum_aixFinvtj_Finvbj[i] = all_sum - cumsum_aixFinvtj_Finvbj[i];
    }

    float t1 = 0.0, t2=0.0, t3=0.0, B = 0.0; //vector<float> t1_all (ai.size());
    FOREACH_i(ai){
        t1 += ai[i] * cumsum_aixFinvtj_Finvbj[i] * f1Finvti_f1Finvbi[i];
        t2 += ai[i] * ai[i] * Finvti[i] * f1Finvti_f1Finvbi[i];

        t3 += ai[i] * ai[i] * (orderStatsKSection_f2(Finvti[i], FinvtiNcdf[i], FinvtiNpdf[i])-
             orderStatsKSection_f2(Finvbi[i], FinvbiNcdf[i], FinvbiNpdf[i]));

        B += ai[i] * f1Finvti_f1Finvbi[i];
    }
    //vec_pointMultiply(vec_pointMultiply(ai, cumsum_aixFinvtj_Finvbj), f1Finvti_f1Finvbi);
    //float t1 = accumulate(t1_all.begin(), t1_all.end(), 0.0);

    float A = 2*(t1+t2-t3);
    B = B*B;

    sigma = float(sqrt(A-B)/sqrt(n));
}
// !!!NOTE: in openCV the index is stored row by row, not like Matlab which is column by column
template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, MatSize size){
    // assert(
    idx = z*size[0]*size[1] + x + y*size[1];
}
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, MatSize size){
    z = idx / (size[0]*size[1]);
    T rmder = idx % (size[0]*size[1]);
    y = rmder/size[1];
    x = rmder-y*size[0];

}

template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, int *size){
    idx = z*size[0]*size[1] + x + y*size[1];
}

template <typename T> T vol_sub2ind(int y, int x, int z, int col_num, T page_sz){
    return z*page_sz + x + y*col_num;
}
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, int *size){
    z = idx / (size[0]*size[1]);
    T rmder = idx - z*(size[0]*size[1]);
    y = rmder/size[1];
    x = rmder-y*size[1];
}

template <typename T> void vec_sub2ind(vector<T> &idx, vector<int> y, vector<int> x, vector<int> z, MatSize size){
    size_t nPixels_slice = size[0]*size[1];
    FOREACH_i(y){
        idx[i] = z[i]*nPixels_slice + x[i] + y[i]*size[1];
    }
}
template <typename T> void vec_ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, vector<int> &z, MatSize size){
    size_t nPixels_slice = size[0]*size[1];
    T rmder;
    FOREACH_i(idx){
        z[i] = idx[i] / nPixels_slice;
        rmder = idx[i]  - z[i] * nPixels_slice;
        y[i] = rmder/size[1];
        x[i] = rmder-y[i]*size[1];
    }
}


template <typename T> size_t overlap_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    size_t fg_sz = 0;
    if (datatype == CV_8U){
        FOREACH_i(vec_idx){
            if(src3d->at<unsigned char>(vec_idx[i]) > threshold_in){
                fg_sz ++;
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i(vec_idx){
            if(src3d->at<float>(vec_idx[i]) > threshold_in){
                fg_sz ++;
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i(vec_idx){
            if(src3d->at<int>(vec_idx[i]) > threshold_in){
                fg_sz ++;
            }
        }
    }
    return fg_sz;
}

template <typename T> bool isempty_mat_vec(Mat *src3d, int datatype, vector<T> vec_idx, float threshold_in){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    //size_t fg_sz = 0;
    if (datatype == CV_8U){
        FOREACH_i(vec_idx){
            if(src3d->at<unsigned char>(vec_idx[i]) > threshold_in){
                return false;
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i(vec_idx){
            if(src3d->at<float>(vec_idx[i]) > threshold_in){
                return false;
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i(vec_idx){
            if(src3d->at<int>(vec_idx[i]) > threshold_in){
                return false;
            }
        }
    }
    return true;
}


template <typename T> void scale_vol(Mat *src3d, int datatype, Mat *dst, float il, float ih, float vl, float vh){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    if (vl > vh){
        normalize(*src3d, *dst, il, ih, NORM_MINMAX, CV_32FC1);
    }else{
        dst->create(src3d->dims, src3d->size, CV_32F);
        float term = (vh-vl) / (ih - il);

        if (datatype == CV_8U){
            FOREACH_i_ptrMAT(src3d){
                dst->at<unsigned char>(i) = (src3d->at<unsigned char>(i) - il) * term + vl;
            }
        }else if (datatype == CV_32F){
            FOREACH_i_ptrMAT(src3d){
                dst->at<float>(i) = (src3d->at<unsigned char>(i) - il) * term + vl;
            }
        }else if (datatype == CV_32S){
            FOREACH_i_ptrMAT(src3d){
                dst->at<int>(i) = (src3d->at<unsigned char>(i) - il) * term + vl;
            }
        }
    }
}

template <typename T> void vec_unique(vector<T> & v){
    sort(v.begin(), v.end());
    typename vector<T>::iterator it;
    it = unique(v.begin(), v.end());
    v.resize(distance(v.begin(),it));
}

/**
 * @brief filterZdirection
 * @param src3d : float
 * @param dst3d : float
 * @param kernel_z : float (border is replicated)
 */
void filterZdirection(Mat* src3d, Mat &dst3d, Mat kernel_z){
    assert(src3d->dims == 3);
    if (dst3d.empty()) dst3d.create(src3d->dims, src3d->size, CV_32F);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    // Filter Z dimension
    float* kernData = (float *)kernel_z.data;
    unsigned kernSize = kernel_z.total(); // should be odd
    int kernMargin = (kernSize - 1)/2;
    float* lineBuffer = new float[z_size + 2*kernMargin];
    for (int y = 0; y < y_size; y++)
    {
        for (int x = 0; x < x_size; x++)
        {
            // Copy along Z dimension into a line buffer
            float* z_ptr = (float*)src3d->data + y * x_size + x;//same as data4smooth.ptr<float>(0, y, x)
            for (int z = 0; z < z_size; z++, z_ptr += xy_size)
            {
                lineBuffer[z + kernMargin] = *z_ptr;
            }

            // Replicate borders
            for (int m = 0; m < kernMargin; m++)
            {
                lineBuffer[m] = lineBuffer[kernMargin];// replicate left side
                lineBuffer[z_size + 2*kernMargin - 1 - m] = lineBuffer[kernMargin + z_size - 1];//replicate right side
            }

            // Filter line buffer 1D - convolution
            z_ptr = (float*)dst3d.data + y * x_size + x;
            for (int z = 0; z < z_size; z++, z_ptr += xy_size)
            {
                *z_ptr = 0.0f;
                for (unsigned k = 0; k < kernSize; k++)
                {
                    *z_ptr += lineBuffer[z+k]*kernData[k];
                }
            }
        }
    }
    delete [] lineBuffer;
}
/**
 * @brief gaussianSmooth3Ddata
 * @param data4smooth
 * @param sigma
 */
void gaussianSmooth3Ddata(Mat &data4smooth, const float sigma[])
{
    assert(data4smooth.dims == 3);
    int x_size  = data4smooth.size[0];
    int y_size  = data4smooth.size[1];
    int z_size  = data4smooth.size[2];
    size_t xy_size = x_size*y_size;

    int kernDimension = ceil(2*sigma[0])+1;
    Mat kernel_row = getGaussianKernel(kernDimension, sigma[0], CV_32F);
    kernDimension = ceil(2*sigma[1])+1;
    Mat kernel_col = getGaussianKernel(kernDimension, sigma[1], CV_32F);
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)data4smooth.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, data4smooth.size, CV_32F, ind);
        sepFilter2D(subMatrix, subMatrix, CV_32F, kernel_row.t(), kernel_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
    if (sigma[2] > 0){
        kernDimension = ceil(2*sigma[2])+1;
        Mat kernel_z = getGaussianKernel(kernDimension, sigma[2], CV_32F);
        // Filter Z dimension
        Mat copyMat;
        data4smooth.copyTo(copyMat);
        filterZdirection(&copyMat, data4smooth, kernel_z);
//        float* kernGauss = (float *)kernel_z.data;
//        unsigned kernSize = kernel_z.total();
//        int kernMargin = (kernSize - 1)/2;
//        float* lineBuffer = new float[z_size + 2*kernMargin];
//        for (int y = 0; y < y_size; y++)
//        {
//            for (int x = 0; x < x_size; x++)
//            {
//                // Copy along Z dimension into a line buffer
//                float* z_ptr = (float*)data4smooth.data + y * x_size + x;//same as data4smooth.ptr<float>(0, y, x)
//                for (int z = 0; z < z_size; z++, z_ptr += xy_size)
//                {
//                    lineBuffer[z + kernMargin] = *z_ptr;
//                }

//                // Replicate borders
//                for (int m = 0; m < kernMargin; m++)
//                {
//                    lineBuffer[m] = lineBuffer[kernMargin];// replicate left side
//                    lineBuffer[z_size + 2*kernMargin - 1 - m] = lineBuffer[kernMargin + z_size - 1];//replicate right side
//                }

//                // Filter line buffer 1D - convolution
//                z_ptr = (float*)data4smooth.data + y * x_size + x;
//                for (int z = 0; z < z_size; z++, z_ptr += xy_size)
//                {
//                    *z_ptr = 0.0f;
//                    for (unsigned k = 0; k < kernSize; k++)
//                    {
//                        *z_ptr += lineBuffer[z+k]*kernGauss[k];
//                    }
//                }
//            }
//        }
//        delete [] lineBuffer;
    }
}
void filterVolume(Mat* src3d, Mat &dst3d, Mat kernel, unsigned direction){
    if (dst3d.empty()){
        dst3d.create(src3d->dims, src3d->size, CV_32F);
        src3d->copyTo(dst3d);
    }
    int x_size  = dst3d.size[0];
    int y_size  = dst3d.size[1];
    int z_size  = dst3d.size[2];
    size_t xy_size = x_size*y_size;

    if (direction == DIRECTION_X){
        assert(kernel.size[0]==1);
        for (int z = 0; z < z_size; z++)
        {
            float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
            Mat subMatrix(2, dst3d.size, CV_32F, ind);
            filter2D(subMatrix, subMatrix, -1, kernel);
        }
    }else if(direction == DIRECTION_Y){
        assert(kernel.size[1]==1);
        for (int z = 0; z < z_size; z++)
        {
            float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
            Mat subMatrix(2, dst3d.size, CV_32F, ind);
            filter2D(subMatrix, subMatrix, -1, kernel);
        }
    }else if(direction == DIRECTION_Z){
        filterZdirection(src3d, dst3d, kernel);
    }else{
        return;
    }
}
/**
 * @brief principalCv2d
 * @param src3d: float mat
 * @param dst3d: float mat
 */
void principalCv2d(Mat* src3d, Mat &dst3d, float sigma[], int minIntensity = 0){
    src3d->copyTo(dst3d);
    gaussianSmooth3Ddata(dst3d, sigma);
    //ccShowSlice3Dmat(&dst3d, CV_32F);
    int x_size  = dst3d.size[0];
    int y_size  = dst3d.size[1];
    int z_size  = dst3d.size[2];
    size_t xy_size = x_size*y_size;

    // Filter XY dimensions for every Z
    Mat kernelx = (Mat_<float>(1,3)<<-0.5, 0, 0.5);
    Mat kernely = (Mat_<float>(3,1)<<-0.5, 0, 0.5);

    int mat_sizes[] = {2,2};
    Mat mat2x2(2, mat_sizes, CV_32F, Scalar(0));
    Mat eigen_values;
    //float min_intensity = (float) minIntensity;
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_32F, ind);
        Mat lx, ly, lxx, lxy, lyy;//lyx,

        //ccShowSlice3Dmat(&subMatrix, CV_32F);
//        subMatrix.copyTo(lx);
//        ccShowSlice3Dmat(&lx, CV_32F);
        filter2D(subMatrix, lx, -1, kernelx);
        filter2D(subMatrix, ly, -1, kernely);
        filter2D(lx, lxy, -1, kernely);
        filter2D(lx, lxx, -1, kernelx);
        //filter2D(ly, lyx, -1, kernelx); // the same as lxy
        filter2D(ly, lyy, -1, kernely);
//        ccShowSlice3Dmat(&lxy, CV_32F);
//        ccShowSlice3Dmat(&lxx, CV_32F);
//        ccShowSlice3Dmat(&lyy, CV_32F);
//        Mat tmp = subMatrix > 0;
//        ccShowSlice3Dmat(&tmp, CV_8U);

//        float *ind2 = (float*)src3d->data + z * xy_size; // sub-matrix pointer
//        Mat subMatrix2(2, dst3d.size, CV_32F, ind2);
//        for (int c = 0; c < y_size; c ++){
//                qInfo("%f", src3d->at<float>(c, 484, y_size));
//        }
//        ccShowSlice3Dmat(&subMatrix2, CV_32F);
        for (int r = 0; r < x_size; r ++){
            size_t increment = r*y_size;
            for (int c = 0; c < y_size; c ++){
                size_t idx = z*xy_size + c + increment;
                if (src3d->at<float>(idx) <= minIntensity){
                    subMatrix.at<float>(r,c) = 0;
                    continue;
                }
                mat2x2.at<float>(0,0) = lxx.at<float>(r,c);
                mat2x2.at<float>(0,1) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,0) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,1) = lyy.at<float>(r,c);
                eigen(mat2x2, eigen_values);
                //float eig_max = MAX(eigen_values.at<float>(0), eigen_values.at<float>(1));
                //subMatrix.at<float>(r,c) = eigen_values.at<float>(0); // the largest eigen value
                dst3d.at<float>(idx) = eigen_values.at<float>(0);
//                if (eig_max > EPSILON_0){
//                    if (new_r){
//                        new_r = false;
//                        qInfo("%d, %d %.4f\n", r, c, eig_max);
//                    }
//                    subMatrix.at<float>(r,c) = eig_max; // the largest eigen value
//                }
            }
        }

        //ccShowSlice3Dmat(&dst3d, CV_32F);
    }
}
/**
 * @brief principalCv3d
 * @param src3d: float mat
 * @param dst3d: float mat
 * @param sigma
 * @param minIntensity
 */
void principalCv3d(Mat* src3d, Mat &dst3d, float sigma[], int minIntensity = 0){
    src3d->copyTo(dst3d);
    gaussianSmooth3Ddata(dst3d, sigma);

    int x_size  = dst3d.size[0];
    int y_size  = dst3d.size[1];
    int z_size  = dst3d.size[2];
    size_t xy_size = x_size*y_size;

    // Filter XY dimensions for every Z
    Mat kernelx = (Mat_<float>(1,3)<<-0.5, 0, 0.5);//(x0+x2)/2
    Mat kernely = (Mat_<float>(3,1)<<-0.5, 0, 0.5);//(y0+y2)/2
    Mat kernelz = (Mat_<float>(1,3)<<-0.5, 0, 0.5);//(z0+z2)/2;
    Mat Dz, Dzz;
    Dz.create(dst3d.dims, dst3d.size, CV_32F);
    Dzz.create(dst3d.dims, dst3d.size, CV_32F);
    filterZdirection(&dst3d, Dz, kernelz);
    //ccShowSlice3Dmat(&Dz, CV_32F, 5);
    filterZdirection(&Dz, Dzz, kernelz);
    //ccShowSlice3Dmat(&Dzz, CV_32F, 5);
    int mat_sizes[] = {3,3};
    Mat mat3x3(2, mat_sizes, CV_32F, Scalar(0));
    Mat eigen_values;
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_32F, ind);
        Mat lx, ly, lxx, lyy, lxy, lzx, lzy;
        filter2D(subMatrix, lx, -1, kernelx);
        filter2D(subMatrix, ly, -1, kernely);
        filter2D(lx, lxx, -1, kernelx);
        filter2D(ly, lyy, -1, kernelx);
        filter2D(lx, lxy, -1, kernely);

        Mat lz(2, dst3d.size, CV_32F, (float*)Dz.data + z * xy_size);
        filter2D(lz, lzx, -1, kernelx);
        filter2D(lz, lzy, -1, kernely);
        //ccShowSlice3Dmat(&Dzz, CV_32F, 5);
        Mat lzz(2, dst3d.size, CV_32F, (float*)Dzz.data + z * xy_size);
        for (int r = 0; r < x_size; r ++){
            size_t increment = r * y_size;
            for (int c = 0; c < y_size; c ++){
                size_t idx = z*xy_size + c + increment;
                if (src3d->at<float>(idx) <= minIntensity){
                    dst3d.at<float>(idx) = 0;
                    continue;
                }
                mat3x3.at<float>(0,0) = lxx.at<float>(r,c);
                mat3x3.at<float>(1,1) = lyy.at<float>(r,c);
                mat3x3.at<float>(2,2) = lzz.at<float>(r,c);

                mat3x3.at<float>(0,1) = lxy.at<float>(r,c);
                mat3x3.at<float>(1,0) = lxy.at<float>(r,c);

                mat3x3.at<float>(0,2) = lzx.at<float>(r,c);
                mat3x3.at<float>(2,0) = lzx.at<float>(r,c);

                mat3x3.at<float>(1,2) = lzy.at<float>(r,c);
                mat3x3.at<float>(2,1) = lzy.at<float>(r,c);


                eigen(mat3x3, eigen_values);
                dst3d.at<float>(idx) = eigen_values.at<float>(0); // the largest eigen value
            }
        }
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernel_row.t(), kernel_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
}




/**
 * @brief varByTruncate: estimate variance using truncated gaussian
 * @param vals4var
 * @param numSigma
 * @param numIter
 * @return
 */
float varByTruncate(vector<float> vals4var, int numSigma, int numIter){
    float sigma = vec_stddev(vals4var);
    float mu = vec_mean(vals4var);
    float t_mu, t_sigma;
    for(int i = 0 ; i< numIter; i++){
        float ub = mu+numSigma*sigma;
        float lb = mu-numSigma*sigma;
        vector<float> tmpVals = vec_atrange(vals4var, ub, lb, true);
        float sigmaTr = vec_stddev(tmpVals);
        //sigmaTr = sqrt(mean(vals(vals<ub & vals>lb).^2));
        float varCorrectionTerm = truncatedGauss(mu, sigma, lb, ub, t_mu, t_sigma); // only lb and ub is usefull if we only want the third output
        sigma = sigmaTr/sqrt(varCorrectionTerm);
        mu = vec_mean(tmpVals);
    }
    return sigma*sigma;
}

/**
 * @brief calVarianceStablization
 * @param src3d: float mat 3d
 * @param validRatio: consider the validRatio pixels in the data (remove pixels with saturated intensity)
 * @param gap
 */
float calVarianceStablization(Mat *src3d, Mat & varMap, vector<float> &varTrend, float validRatio, int gap){
    Mat meanVal; //(src3d->dims, src3d->size, CV_32F, Scalar(0));
    src3d->copyTo(meanVal);
    Mat kernelx, kernely, kernelz, kernel;
    float denominator = 0;
    if(gap == 0){
        kernelx = (Mat_<float>(1,3)<<1.0, 1.0, 1.0);
        kernely = (Mat_<float>(3,1)<<1.0, 1.0, 1.0);
        denominator = 9;
        //kernelz = (Mat_<float>(1,3)<<1.0, 1.0, 1.0);
    }else if(gap == 1){
        kernelx = (Mat_<float>(1,5)<<1.0, 0.0, 0.0, 0.0, 1.0);
        kernely = (Mat_<float>(5,1)<<1.0, 0.0, 0.0, 0.0, 1.0);
        denominator = 16;
        //kernelz = (Mat_<float>(1,5)<<1.0, 0.0, 0.0, 0.0, 1.0);
    }else if(gap == 2){
        kernelx = (Mat_<float>(1,7)<<1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        kernely = (Mat_<float>(7,1)<<1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        kernel = Mat::zeros(7,7,CV_32F);
        kernel.at<float>(3) = 1;
        kernel.at<float>(8) = 1; kernel.at<float>(9) = 1; kernel.at<float>(11) = 1; kernel.at<float>(12) = 1;
        kernel.at<float>(15) = 1; kernel.at<float>(19) = 1; kernel.at<float>(21) = 1; kernel.at<float>(27) = 1;
        kernel.at<float>(29) = 1; kernel.at<float>(33) = 1; kernel.at<float>(36) = 1; kernel.at<float>(37) = 1;
        kernel.at<float>(39) = 1; kernel.at<float>(40) = 1; kernel.at<float>(45) = 1;
        denominator = 16;
        kernel /= denominator;
        //kernelz = (Mat_<float>(1,7)<<1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)meanVal.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, meanVal.size, CV_32F, ind);
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernely, kernelx, Point(-1,-1), 0.0, BORDER_REPLICATE);
        filter2D(subMatrix, subMatrix, -1, kernel);
    }
    //sepFilter2D(*src3d, meanVal, CV_32F, kernelx, kernely, Point(-1,-1), 0.0, BORDER_REPLICATE);
    double min_intensity, max_intensity;
//    minMaxIdx(meanVal, &min_intensity, &max_intensity);
//    meanVal /= denominator;
    Mat diff = *src3d - meanVal;
    int levels = 200;

    minMaxIdx(meanVal, &min_intensity, &max_intensity);

    double unit_intensity = (max_intensity - min_intensity) / (levels-1);
    varTrend.resize(levels);
    //for (int i = 0; i < levels; i++) xx[0] = 0;
    fill(varTrend.begin(), varTrend.end(), 0);


    vector<vector<long>> numElements (levels);
    long testedElements = 0, validElements = 0;
    int cur_level = 0;
    for (int i = gap+1; i < x_size-gap; i++){
        for (int j = gap+1; j< y_size-gap; j++){
            for (int k = 0; k < z_size; k++){
                size_t idx = j + i * y_size + k * xy_size;
                cur_level = int((meanVal.at<float>(idx) - min_intensity) / unit_intensity);
                if(cur_level > 0){
                    validElements ++;
                    numElements[cur_level].push_back(idx);
                }
            }
        }
    }
    int target_level = levels;
    cur_level = 1;
    int nei_vec_step = 1;

    while(cur_level < levels){
        nei_vec_step = 1;
        vector<long> cur_level_idx = numElements[cur_level];
        while (cur_level_idx.size() < 100 && nei_vec_step < 5){
            if (cur_level - nei_vec_step > 0){
                cur_level_idx.insert(cur_level_idx.end(),
                                              numElements[cur_level-nei_vec_step].begin(),
                                            numElements[cur_level-nei_vec_step].end());
            }
            if (cur_level + nei_vec_step < levels){
                cur_level_idx.insert(cur_level_idx.end(),
                                              numElements[cur_level+nei_vec_step].begin(),
                                            numElements[cur_level+nei_vec_step].end());
            }
            nei_vec_step ++;
        }
        vector<float> vals4var (cur_level_idx.size());
        for (size_t j = 0; j < cur_level_idx.size(); j++){
            vals4var[j] = diff.at<float>(cur_level_idx[j]);
        }
        float varCur = varByTruncate(vals4var, 2, 3);
        varTrend[cur_level] = (denominator/(denominator+1))*varCur;
        testedElements += cur_level_idx.size();
        if (cur_level < target_level &&
                ((float)testedElements/validElements > validRatio)){
                //if we only want to use a fixed point to reprsent all the variance
                target_level = cur_level;
        }
        cur_level++;
    }
    if (max_intensity>20){
        for (int i=0;i<3;i++) varTrend[i] = 0;
    }
    else{
        for (int i=0;i<15;i++) varTrend[i] = 0;
    }
    if (varMap.empty()){
        varMap = Mat::zeros(meanVal.dims, meanVal.size, CV_32F);
    }
    float min_valid_intensity = min_intensity+unit_intensity;
    for(size_t i = 0; i < xy_size*z_size; i++){
        if (meanVal.at<float>(i) > min_valid_intensity){
            cur_level = (int) ((meanVal.at<float>(i) - min_intensity)/unit_intensity);
            if (cur_level >= levels) cur_level = levels - 1;
            if (cur_level < 0) cur_level = 0;
            varMap.at<float>(i) = varTrend[cur_level];
        }
    }

    //varTrend = *xx;
    return varTrend[target_level];
}

/**
 * @brief connectedComponents3d
 * @param src3d: boolean but represented by CV_8U
 * @param dst3d: int represented by CV_32S
 * @param connect: 4, 8 for 2d and 6,10,26 for 3d
 */
int connectedComponents3d(Mat* src3d, Mat &dst3d, int connect){
    //assert(src3d->type()==bool);
    assert(src3d->dims == 3);
    assert(connect == 4 || connect == 8 || connect == 6 || connect == 10 || connect == 26);
    size_t numCC;
    if (connect == 4 || connect == 8){
        int x_size  = src3d->size[0];
        int y_size  = src3d->size[1];
        int z_size  = src3d->size[2];
        size_t xy_size = x_size*y_size;
        numCC = 0;
        //dst3d = Mat(src3d->dims, src3d->size, CV_32S);
        for (int z = 0; z < z_size; z++)
        {
            int *ind = (int*)dst3d.data + z * xy_size; // sub-matrix pointer
            Mat subMatrix_dst(2, dst3d.size, CV_32S, ind);
            unsigned char *ind2 = (unsigned char*)src3d->data + z * xy_size; // sub-matrix pointer
            Mat subMatrix_src(2, src3d->size, CV_8U, ind2);
            numCC += connectedComponents(subMatrix_src, subMatrix_dst, connect);
            //double Max_label, min_label;
            //minMaxIdx(subMatrix_dst, &min_label, &Max_label);
            //ccShowSliceLabelMat(&subMatrix_dst);
        }
    }else{
        //dst3d.create(src3d->dims, src3d->size, CV_32S);
        int* out = cc3d::connected_components3d<unsigned char, int>(src3d->data, src3d->size[0], src3d->size[1],
                src3d->size[2], connect);
        FOREACH_i_MAT(dst3d){
            dst3d.at<int>(i) = out[i];
        }
        numCC = *max_element(out, out + dst3d.total());
        delete[] out;
        //int current_id = 0;
//        if(connect == 6){
//            for (int i = 0; i < x_size; i++){
//                for (int j = 0; j< y_size; j++){
//                    for (int k = 0; k < z_size; k++){
//                        size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                        if (src3d->at<unsigned char>(idx) == 0){
//                            continue;
//                        }
//                        if (dst3d.at<int>(idx) == 0){
//                            numCC ++;
//                            dst3d.at<int>(idx) = numCC;
//                        }
//                        current_id = dst3d.at<int>(idx);
//                        if ((i + 1) < x_size && src3d->at<unsigned char>(idx + y_size)>0){
//                            dst3d.at<int>(idx + y_size) = current_id;
//                        }
//                        if ((j + 1) < y_size && src3d->at<unsigned char>(idx + 1)>0){
//                            dst3d.at<int>(idx + 1) = current_id;
//                        }
//                        if ((k + 1) < z_size && src3d->at<unsigned char>(idx + xy_size)>0){
//                            dst3d.at<int>(idx + xy_size) = current_id;
//                        }
//                    }
//                }
//            }
//        }

//        if(connect == 10){
//            for (int i = 0; i < x_size; i++){
//                for (int j = 0; j< y_size; j++){
//                    for (int k = 0; k < z_size; k++){
//                        size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                        if (src3d->at<unsigned char>(idx) == 0){
//                            continue;
//                        }
//                        if (dst3d.at<int>(idx) == 0){
//                            numCC ++;
//                            dst3d.at<int>(idx) = numCC;
//                        }
//                        current_id = dst3d.at<int>(idx);
//                        if ((i + 1) < x_size && src3d->at<unsigned char>(idx + y_size)>0){
//                            dst3d.at<int>(idx + y_size) = current_id;
//                        }
//                        if ((j + 1) < y_size && src3d->at<unsigned char>(idx + 1)>0){
//                            dst3d.at<int>(idx + 1) = current_id;
//                        }
//                        if ((k + 1) < z_size && src3d->at<unsigned char>(idx + xy_size)>0){
//                            dst3d.at<int>(idx + xy_size) = current_id;
//                        }

//                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned char>(idx + y_size + 1)>0){
//                            dst3d.at<int>(idx + y_size + 1) = current_id;
//                        }
//                    }
//                }
//            }
//        }

//        if(connect == 26){
//            for (int i = 0; i < x_size; i++){
//                for (int j = 0; j< y_size; j++){
//                    for (int k = 0; k < z_size; k++){
//                        size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                        if (src3d->at<unsigned char>(idx) == 0){
//                            continue;
//                        }
//                        if (dst3d.at<int>(idx) == 0){
//                            numCC ++;
//                            dst3d.at<int>(idx) = numCC;
//                        }
//                        current_id = dst3d.at<int>(idx);
//                        if ((i + 1) < x_size && src3d->at<unsigned char>(idx + y_size)>0){
//                            dst3d.at<int>(idx + y_size) = current_id;
//                        }
//                        if ((j + 1) < y_size && src3d->at<unsigned char>(idx + 1)>0){
//                            dst3d.at<int>(idx + 1) = current_id;
//                        }
//                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned char>(idx + y_size + 1)>0){
//                            dst3d.at<int>(idx + y_size + 1) = current_id;
//                        }

//                        if ((k + 1) < z_size){
//                            if (src3d->at<unsigned char>(idx + xy_size)>0){
//                                dst3d.at<int>(idx + xy_size) = current_id;
//                            }
//                            if ((j + 1) < y_size && src3d->at<unsigned char>(idx + xy_size + 1)>0){
//                                dst3d.at<int>(idx + xy_size + 1) = current_id;
//                            }
//                            if ((i + 1) < x_size && src3d->at<unsigned char>(idx + xy_size + y_size)>0){
//                                dst3d.at<int>(idx + xy_size + y_size) = current_id;
//                            }
//                            if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned char>(idx + xy_size + y_size + 1)>0){
//                                dst3d.at<int>(idx + xy_size + y_size + 1) = current_id;
//                            }
//                        }
//                    }
//                }
//            }
//        }

    }

    return numCC;
}
int floatMap2idMap(Mat* src3d, Mat &dst3d, int connect){
    assert(src3d->dims == 3);
    assert(connect == 6 || connect == 10 || connect == 26);
//    int x_size  = src3d->size[0];
//    int y_size  = src3d->size[1];
//    int z_size  = src3d->size[2];
    int numCC = 0;
    //cc3d::connected_components3d<float, int>(src3d->data, src3d->size[0], src3d->size[1], src3d->size[2], connect);
    float *ind = (float *) src3d->data;
    int* out = cc3d::connected_components3d<float, int>(ind, src3d->size[0], src3d->size[1],
            src3d->size[2], connect);
    //dst3d = Mat(src3d->dims, src3d->size, CV_32S);
    FOREACH_i_MAT(dst3d){
        dst3d.at<int>(i) = out[i];
    }
    numCC = *max_element(out, out + dst3d.total());
//    size_t xy_size = x_size*y_size;
//    size_t numCC = 0;

//    if (dst3d.empty()) dst3d = Mat::zeros(src3d->dims, src3d->size, CV_32S);
//    int current_id = 0;
//    if(connect == 6){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
//                for (int k = 0; k < z_size; k++){
//                    size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                    float curr_val = src3d->at<float>(idx);
//                    if (curr_val < 1e-5){
//                        continue;
//                    }
//                    if (dst3d.at<int>(idx) == 0){
//                        numCC ++;
//                        dst3d.at<int>(idx) = numCC;
//                    }
//                    current_id = dst3d.at<int>(idx);
//                    if ((i + 1) < x_size && (src3d->at<float>(idx + y_size) - curr_val) < 1e-5){
//                        dst3d.at<int>(idx + y_size) = current_id;
//                    }
//                    if ((j + 1) < y_size && (src3d->at<float>(idx + 1)- curr_val) < 1e-5){
//                        dst3d.at<int>(idx + 1) = current_id;
//                    }
//                    if ((k + 1) < z_size && (src3d->at<float>(idx + xy_size)- curr_val) < 1e-5){
//                        dst3d.at<int>(idx + xy_size) = current_id;
//                    }
//                }
//            }
//        }
//    }

//    if(connect == 10){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
//                for (int k = 0; k < z_size; k++){
//                    size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                    float curr_val = src3d->at<float>(idx);
//                    if (curr_val < EPSILON_0){
//                        continue;
//                    }
//                    if (dst3d.at<int>(idx) == 0){
//                        numCC ++;
//                        dst3d.at<int>(idx) = numCC;
//                    }
//                    current_id = dst3d.at<int>(idx);
//                    if ((i + 1) < x_size && (src3d->at<float>(idx + y_size) - curr_val) < EPSILON_0){
//                        dst3d.at<int>(idx + y_size) = current_id;
//                    }
//                    if ((j + 1) < y_size && (src3d->at<float>(idx + 1) - curr_val) < EPSILON_0){
//                        dst3d.at<int>(idx + 1) = current_id;
//                    }
//                    if ((k + 1) < z_size && (src3d->at<float>(idx + xy_size) - curr_val) < EPSILON_0){
//                        dst3d.at<int>(idx + xy_size) = current_id;
//                    }

//                    if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(idx + y_size + 1) - curr_val) < EPSILON_0){
//                        dst3d.at<int>(idx + y_size + 1) = current_id;
//                    }
//                }
//            }
//        }
//    }

//    if(connect == 26){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
//                for (int k = 0; k < z_size; k++){
//                    size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
//                    float curr_val = src3d->at<float>(idx);
//                    if (curr_val < 1e-5){
//                        continue;
//                    }
//                    if (dst3d.at<int>(idx) == 0){
//                        numCC ++;
//                        dst3d.at<int>(idx) = numCC;
//                    }
//                    current_id = dst3d.at<int>(idx);
//                    if ((i + 1) < x_size && (src3d->at<float>(idx + y_size) - curr_val) < 1e-5){
//                        dst3d.at<int>(idx + y_size) = current_id;
//                    }
//                    if ((j + 1) < y_size && (src3d->at<float>(idx + 1) - curr_val) < 1e-5){
//                        dst3d.at<int>(idx + 1) = current_id;
//                    }
//                    if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(idx + y_size + 1) - curr_val) < 1e-5){
//                        dst3d.at<int>(idx + y_size + 1) = current_id;
//                    }

//                    if ((k + 1) < z_size){
//                        if ((src3d->at<float>(idx + xy_size) - curr_val) < 1e-5){
//                            dst3d.at<int>(idx + xy_size) = current_id;
//                        }
//                        if ((j + 1) < y_size && (src3d->at<float>(idx + xy_size + 1) - curr_val) < 1e-5){
//                            dst3d.at<int>(idx + xy_size + 1) = current_id;
//                        }
//                        if ((i + 1) < x_size && (src3d->at<float>(idx + xy_size + y_size) - curr_val) < 1e-5){
//                            dst3d.at<int>(idx + xy_size + y_size) = current_id;
//                        }
//                        if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(idx + xy_size + y_size + 1) - curr_val) < 1e-5){
//                            dst3d.at<int>(idx + y_size + 1+1) = current_id;
//                        }
//                    }
//                }
//            }
//        }
//    }

    return numCC;
}
int rearrangeIdMap(Mat* src3d, Mat &dst3d, vector<size_t> &idMap){
    double min, max;
    minMaxIdx(*src3d, &min, &max);
    idMap.resize(round(max) + 1);
    fill(idMap.begin(), idMap.end(), 0);
    for(size_t i = 0; i < src3d->total(); i++){
        if (src3d->at<int>(i) > 0)
            idMap[src3d->at<int>(i)] = 1;
    }
    int dst_region_num = 0;
    FOREACH_i(idMap){
        if(idMap[i] > 0){
            dst_region_num++;
            idMap[i] = dst_region_num;
        }
    }
    for(size_t i = 0; i < src3d->total(); i++){
        if (src3d->at<int>(i) > 0)
            dst3d.at<int>(i) = idMap[src3d->at<int>(i)];
    }
    return dst_region_num;
}
void getRange(vector<int> idx_sub, int shift, int bound, Range &out_range){
    int min_v = max(0, *min_element(idx_sub.begin(), idx_sub.end()) - shift);
    int max_v = min(bound, *max_element(idx_sub.begin(), idx_sub.end()) + shift);
    out_range = Range(min_v, max_v);
}
void regionAvgIntensity(Mat* src3dFloatData, Mat* src3dIdMap, vector<float> &avgIntensities){
    vector<size_t> region_sz(avgIntensities.size());
    fill(region_sz.begin(), region_sz.end(), 0);
    for(size_t i=0; i<src3dFloatData->total(); i++){
        if (src3dIdMap->at<int>(i)>0){
            avgIntensities[src3dIdMap->at<int>(i)-1] += src3dFloatData->at<float>(i);
            region_sz[src3dIdMap->at<int>(i) - 1] ++;
        }
    }
    FOREACH_i(region_sz){
        avgIntensities[i] /= region_sz[i];
    }
}
/**
 * @brief extractVoxList: return the voxel idx of all connected component
 * @param label3d
 * @param voxList
 * @param numCC
 * @param bk_extract: background also count as a region
 */
void extractVoxIdxList(Mat *label3d, vector<vector<size_t>> &voxList, int numCC, bool bk_extract = false){
    if (bk_extract){
        voxList.resize(numCC+1);
        for (size_t i = 0; i < label3d->total(); i++){
             voxList[label3d->at<int>(i)].push_back(i);
        }
    }else{
        voxList.resize(numCC);
        for (size_t i = 0; i < label3d->total(); i++){
            if(label3d->at<int>(i) > 0){
                voxList[label3d->at<int>(i)-1].push_back(i);
            }
        }
    }
}

/**
 * @brief extractVolume: return the size of all connected component
 * @param label3d: cv_32s
 * @param voxSzList
 * @param numCC
 */
void extractVolume(Mat *label3d, vector<size_t> &voxSzList, int numCC){
    voxSzList.resize(numCC);
    fill(voxSzList.begin(), voxSzList.end(), 0);
    for (size_t i = 0; i < label3d->total(); i++){
        if(label3d->at<int>(i) > 0){
            voxSzList[label3d->at<int>(i)-1]++;
        }
    }
}
/**
 * @brief removeSmallCC: remove the connected component that smaller than min_size
 * @param label3d
 * @param numCC
 * @param min_size
 * @param relabel
 */
bool removeSmallCC(Mat &label3d, int &numCC, size_t min_size, bool relabel = true){
    vector<size_t> cc_size(numCC);
    fill(cc_size.begin(), cc_size.end(), 0);
    for (size_t i = 0; i < label3d.total(); i++){
        if(label3d.at<int>(i) > 0){
            cc_size[label3d.at<int>(i) - 1]++;
        }
    }
    bool cc_removed = false;
    vector<size_t> newlabel(numCC);
    fill(newlabel.begin(), newlabel.end(), 0);
    int rm_cc_cnt = 0;
    for (int i = 0; i < numCC; i ++){
        if (cc_size[i] >= min_size){
            newlabel[i] = rm_cc_cnt;
            rm_cc_cnt ++ ;
        }
    }
    if(numCC > rm_cc_cnt){
        cc_removed = true;
        if(relabel){
            numCC = rm_cc_cnt;
            for (size_t i = 0; i < label3d.total(); i++){
                if(label3d.at<int>(i) > 0){
                    label3d.at<int>(i) = newlabel[label3d.at<int>(i) - 1];
                }
            }
        }
        else{
            // no change to numCC, only set invalid CC to 0
            for (size_t i = 0; i < label3d.total(); i++){
                if(label3d.at<int>(i) > 0 && cc_size[label3d.at<int>(i) - 1] < min_size){
                    label3d.at<int>(i) = 0;
                }
            }
        }
    }
    return cc_removed;
}

/**
 * @brief volumeDilate: x,y direction is formal, while z direction is simplified
 * @param src3d: boolean (but represented by CV_8U)
 * @param dst3d
 * @param kernel
 */
void volumeDilate(Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type){
//    int dilation_type = 0;
//    if( dilation_elem == 0 ){ dilation_type = MORPH_RECT; }
//    else if( dilation_elem == 1 ){ dilation_type = MORPH_CROSS; }
//    else if( dilation_elem == 2) { dilation_type = MORPH_ELLIPSE; }
    assert(src3d->dims == 3);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;

    src3d->copyTo(dst3d);
    Mat element = getStructuringElement( dilation_type,
                           Size( 2*radiusValues[0] + 1, 2*radiusValues[1]+1 ),
                           Point( radiusValues[0], radiusValues[1] ) );
    int min_valid_z = z_size - 1, max_valid_z = 0;
    bool valid_z;
    for (int z = 0; z < z_size; z++)
    {
        valid_z = false;
        unsigned char *ind = (unsigned char*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_8U, ind);
        for (size_t i = 0; i < subMatrix.total(); i++){
            if(subMatrix.at<int>(i) > 0){
                if(min_valid_z > z ){
                    min_valid_z = z;
                }else if(max_valid_z < z){
                    max_valid_z = z;
                }
                valid_z = true;
                break;
            }
        }
        if(valid_z){
            dilate(subMatrix, subMatrix, element);
        }
    }
    //unsigned char max_val;
    if (radiusValues[2] > 0 && max_valid_z >= min_valid_z){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
        unsigned char *ind = (unsigned char*)dst3d.data + min_valid_z * xy_size; // sub-matrix pointer
        Mat subMatrix_up(2, dst3d.size, CV_8U, ind);
        for (int k = min_valid_z - 1; k >= max(0, min_valid_z - radiusValues[2]); k--){
            ind = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind);
            subMatrix_up.copyTo(subMatrix2);
        }

        ind = (unsigned char*)dst3d.data + max_valid_z * xy_size; // sub-matrix pointer
        Mat subMatrix_down(2, dst3d.size, CV_8U, ind);
        for (int k = max_valid_z + 1; k <= max(z_size - 1, max_valid_z + radiusValues[2]); k++){
            ind = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind);
            subMatrix_down.copyTo(subMatrix2);
        }

    }
}

/**
 * @brief volumeErode: x,y direction is formal, while z direction is fixed to morph_corss
 * @param src3d: boolean (but represented by CV_8U)
 * @param dst3d
 * @param kernel
 */
void volumeErode(Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type){
    // there is a question here: if dst3d and src3d point to the same memory address?? what will happen
//    int dilation_type = 0;
//    if( dilation_elem == 0 ){ dilation_type = MORPH_RECT; }
//    else if( dilation_elem == 1 ){ dilation_type = MORPH_CROSS; }
//    else if( dilation_elem == 2) { dilation_type = MORPH_ELLIPSE; }
    //assert(src3d != dst3d);
    assert(src3d->dims == 3);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    //Mat org_src3d;
    src3d->copyTo(dst3d);
    Mat element = getStructuringElement( dilation_type,
                           Size( 2*radiusValues[0] + 1, 2*radiusValues[1]+1 ),
                           Point( radiusValues[0], radiusValues[1] ) );
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_8U, ind);
        erode(subMatrix, subMatrix, element);
    }
    unsigned char min_val;
    if (radiusValues[2] > 0){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    size_t idx = vol_sub2ind(i, j, k, y_size, xy_size);
                    min_val = src3d->at<unsigned char>(idx);
                    for (int kk = 1; kk < radiusValues[2]; kk++){
                        if (k+kk < z_size){
                            min_val = min(src3d->at<unsigned char>(idx + kk * xy_size), min_val);
                        }
                        if (k-kk >= 0){
                            min_val = min(src3d->at<unsigned char>(idx - kk * xy_size), min_val);
                        }
                    }
                    dst3d.at<int>(idx) = min_val;
                }
            }
        }
    }
}

void volumeWrite(Mat *src3d, string filename){
    vector<Mat> imgs;
    split(*src3d, imgs);
    imwrite(filename, imgs);
    printf("Multiple files saved in test.tiff\n");
}
/**
 * @brief findUnrelatedCC: return the Mat containing regions related to reference
 * @param src3d4testing : CV_32S
 * @param src3d4reference : CV_8U
 * @param dst3d : CV_8U
 */
bool findRelatedCC(Mat *src3d4testing, int numCC, Mat *src3d4reference, Mat &dst3d){
    vector<bool> existInReference(numCC);
    fill(existInReference.begin(), existInReference.end(), false);

    FOREACH_i_ptrMAT(src3d4testing){
        if(src3d4testing->at<int>(i) > 0 && src3d4reference->at<unsigned char>(i) > 0){
            existInReference[src3d4testing->at<int>(i)-1] = true;
        }
    }

    dst3d = Mat::zeros(src3d4reference->dims, src3d4reference->size, CV_8U);
    bool found = false;
    FOREACH_i_MAT(dst3d){
        if(src3d4testing->at<int>(i) > 0 && existInReference[src3d4testing->at<int>(i)-1]){
            dst3d.at<unsigned char>(i) = 255;
            found = true;
        }
    }
    return found;
}
void extractVoxIdxList(Mat *label3d, vector<vector<int>> &voxList, int numCC, bool bk_extract = false){
    if (bk_extract){
        voxList.resize(numCC+1);
        for (size_t i = 0; i < label3d->total(); i++){
             voxList[label3d->at<int>(i)].push_back(i);
        }
    }else{
        voxList.resize(numCC);
        for (size_t i = 0; i < label3d->total(); i++){
            if(label3d->at<int>(i) > 0){
                voxList[label3d->at<int>(i)-1].push_back(i);
            }
        }
    }
}
/**
 * @brief validRegionExtract: keep the regions that related to the mask, if no mask is provided
 * keep the region with largest mask
 * @param binary_3d
 * @param binary_mask
 * @param connect
 */
void validSingleRegionExtract(Mat &binary_3d, Mat *binary_mask, int connect){
    Mat label_map;
    int n = connectedComponents3d(&binary_3d, label_map, connect);
    if(n<=1) return;

    if(binary_mask == nullptr){ // only keep the largest one
        vector<vector<size_t>> voxIdxList;
        extractVoxIdxList(&label_map, voxIdxList, n);
        size_t max_id = 0, max_sz = 0;
        FOREACH_i(voxIdxList){
            if (voxIdxList[i].size() > max_sz){
                max_id = i + 1;
                max_sz = voxIdxList[i].size();
            }
        }
        for(size_t i = 0; i<binary_3d.total(); i++){
            if (binary_3d.at<unsigned char>(i) != max_id){
                binary_3d.at<unsigned char>(i) = 0;
            }
        }
    }else{ // keep the region related the mask
//        vector<vector<size_t>> voxIdxList;
//        extractVoxIdxList(&label_map, voxIdxList, n);

//        FOREACH_i(voxIdxList){
//            bool mask_covered = false;
//            for(size_t j = 0; j<voxIdxList[i].size(); j++){
//                if (binary_mask->at<unsigned char>(voxIdxList[i][j]) > 0){
//                    mask_covered = true;
//                    break;
//                }
//            }
//            if (!mask_covered){
//                for(size_t j = 0; j<voxIdxList[i].size(); j++){
//                    binary_3d.at<unsigned char>(voxIdxList[i][j]) = 0;
//                }
//            }
//        }
        findRelatedCC(&label_map, n, binary_mask, binary_3d);
    }

}
/**
 * @brief largestRegionIdExtract
 * @param label_map: CV_32S (int) label map
 * @param mask
 * @param connect
 * @return
 */
int largestRegionIdExtract(Mat *label_map, int numCC, Mat *mask){
    vector<size_t> cc_sz(numCC);
    fill(cc_sz.begin(), cc_sz.end(), 0);
    if (mask == nullptr){
        FOREACH_i_ptrMAT(label_map){
            if(label_map->at<int>(i) > 0){
                cc_sz[label_map->at<int>(i) - 1] ++;
            }
        }
    }else{
        FOREACH_i_ptrMAT(label_map){
            if(label_map->at<int>(i) > 0 && mask->at<int>(i) > 0){
                cc_sz[label_map->at<int>(i) - 1] ++;
            }
        }
    }
    size_t largest_id;
    vec_max(cc_sz, largest_id);

    return (int)largest_id;
}
size_t fgMapSize(Mat *src3d, int datatype, float threshold_in = 0){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    size_t fg_sz = 0;
    if (datatype == CV_8U){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<unsigned char>(i) > threshold_in){
                fg_sz ++;
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<float>(i) > threshold_in){
                fg_sz ++;
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<int>(i) > threshold_in){
                fg_sz ++;
            }
        }
    }
    return fg_sz;
}
bool isempty(Mat *src3d, int datatype, float threshold_in = 0){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    //size_t fg_sz = 0;
    if (datatype == CV_8U){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<unsigned char>(i) > threshold_in){
                return false;
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<float>(i) > threshold_in){
                return false;
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<int>(i) > threshold_in){
                return false;
            }
        }
    }
    return true;
}

vector<size_t> fgMapIdx(Mat *src3d, int datatype, float threshold_in){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    vector<size_t> fg_Idx;
    if (datatype == CV_8U){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<unsigned char>(i) > threshold_in){
                fg_Idx.push_back(i);
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<float>(i) > threshold_in){
                fg_Idx.push_back(i);
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<int>(i) > threshold_in){
                fg_Idx.push_back(i);
            }
        }
    }
    return fg_Idx;
}
/**
 * @brief extractValsGivenMask: return the values of foreground
 * @param val3d
 * @param src3d: mask for value extraction
 * @param datatype
 * @param threshold_in
 * @return
 */
vector<float> extractValsGivenMask(Mat *val3d, Mat *src3d, int datatype, float threshold_in){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    vector<float> fg_vals;
    if (datatype == CV_8U){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<unsigned char>(i) > threshold_in){
                fg_vals.push_back((float)val3d->at<unsigned char>(i));
            }
        }
    }else if (datatype == CV_32F){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<float>(i) > threshold_in){
                fg_vals.push_back(val3d->at<float>(i));
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<int>(i) > threshold_in){
                fg_vals.push_back((float)val3d->at<int>(i));
            }
        }
    }
    return fg_vals;
}

vector<float> extractValsGivenIdx(Mat *vol3d, vector<size_t> idx, int datatype){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    vector<float> fg_vals;
    if (datatype == CV_8U){
        FOREACH_i(idx){
            fg_vals.push_back((float)vol3d->at<unsigned char>(i));
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            fg_vals.push_back(vol3d->at<float>(i));
        }
    }else if (datatype == CV_32S){
        FOREACH_i(idx){
            fg_vals.push_back((float)vol3d->at<int>(i));
        }
    }
    return fg_vals;
}

double extractSumGivenIdx(Mat *vol3d, vector<size_t> idx, int datatype){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    double sum = 0;
    if (datatype == CV_8U){
        FOREACH_i(idx){
            sum += (double)vol3d->at<unsigned char>(i);
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            sum += (double)vol3d->at<float>(i);
        }
    }else if (datatype == CV_32S){
        FOREACH_i(idx){
            sum += (double)(float)vol3d->at<int>(i);
        }
    }
    return sum;
}
/**
 * @brief findUnrelatedCC: return the Mat containing regions unrelated to reference
 * @param src3d4testing : CV_32S
 * @param src3d4reference : CV_8U
 * @param dst3d : CV_8U
 */
bool findUnrelatedCC(Mat *src3d4testing, int numCC, Mat *src3d4reference, Mat &dst3d){
    vector<bool> existInReference(numCC);
    fill(existInReference.begin(), existInReference.end(), false);

    FOREACH_i_ptrMAT(src3d4testing){
        if(src3d4testing->at<int>(i) > 0 && src3d4reference->at<unsigned char>(i) > 0){
            existInReference[src3d4testing->at<int>(i)-1] = true;
        }
    }

    dst3d = Mat::zeros(src3d4reference->dims, src3d4reference->size, CV_8U);
    bool found = false;
    FOREACH_i_MAT(dst3d){
        if(src3d4testing->at<int>(i) > 0 && !existInReference[src3d4testing->at<int>(i)-1]){
            dst3d.at<unsigned char>(i) = 255;
            found = true;
        }
    }
    return found;
}


bool inField( int r, int c, int z, int *sz )
{
  if( r < 0 || r >= sz[0] ) return false;
  if( c < 0 || c >= sz[1] ) return false;
  if( z < 0 || z >= sz[2] ) return false;
    return true;
}
bool inField( int r, int c, int *sz )
{
  if( r < 0 || r >= sz[0] ) return false;
  if( c < 0 || c >= sz[1] ) return false;
    return true;
}
bool isOnBoundary2d(Mat *fgMap, int y, int x, int z){
    int im_sz[] = {fgMap->size[0], fgMap->size[1]};
    int page_sz = im_sz[0] * im_sz[1];
    for(int i=0; i < 8; i++){
        if(inField(y + n8_y[i], x + n8_x[i], im_sz)
                && fgMap->at<unsigned char>(vol_sub2ind(y + n8_y[i], x + n8_x[i], z, im_sz[1], page_sz)) == 0){
            return true;
        }
    }
    return false;
}

bool isOnBoundary2d(Mat *fgMap, size_t idx){
    int y,x,z;
    vol_ind2sub(idx, y, x, z, fgMap->size);
    return isOnBoundary2d(fgMap, y, x, z);
}


void neighbor_idx(vector<size_t> idx, vector<size_t> &center_idx, vector<size_t> &nei_idx, int sz[], int connect){
    nei_idx.resize(idx.size() * connect);
    center_idx.resize(idx.size() * connect);

    vector<int> n_y(connect), n_x(connect), n_z(connect);
    if(connect == 8){
        n_y = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
        n_x = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
        n_z = {  0,  0,  0,  0, 0, 0,  0, 0 };
    }else if(connect == 4){
        n_y = { -1,  1,  0, 0 };// 8 shifts to neighbors
        n_x = {  0,  0, -1, 1 };// used in functions
        n_z = {  0,  0,  0, 0 };
    }else if (connect == 6){
        n_y = { -1, 1,  0, 0, 0,  0 };
        n_x = { 0,  0, -1, 1, 0,  0 };
        n_z = { 0,  0,  0, 0, 1, -1 };
    }else if(connect == 26){
        n_y = { -1, -1, -1,  1, 1, 1,  0, 0, -1, -1, -1,  1,  1,  1,  0,  0, 0,  -1, -1, -1,  1, 1, 1,  0, 0, 0 };// 8 shifts to neighbors
        n_x = { -1,  0,  1, -1, 0, 1, -1, 1, -1,  0,  1, -1,  0,  1, -1,  1, 0,  -1,  0,  1, -1, 0, 1, -1, 1, 0 };// used in functions
        n_z = {  0,  0,  0,  0, 0, 0,  0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1, 1, 1,  1, 1, 1 };
    }else if(connect == 10){
        n_y = { -1, -1, -1,  1, 1, 1,  0, 0, 0,  0 };
        n_x = { -1,  0,  1, -1, 0, 1, -1, 1, 0,  0 };
        n_z = {  0,  0,  0,  0, 0, 0,  0, 0, 1, -1 };
    }else if(connect == 18){
        n_y = { -1, -1, -1,  1, 1, 1,  0, 0, -1,  1,  0,  0, 0, -1, 1,  0, 0, 0 };// 8 shifts to neighbors
        n_x = { -1,  0,  1, -1, 0, 1, -1, 1,  0,  0, -1,  1, 0,  0, 0, -1, 1, 0 };// used in functions
        n_z = {  0,  0,  0,  0, 0, 0,  0, 0, -1, -1, -1, -1, -1, 1, 1,  1, 1, 1 };
    }
    int x, y, z, remain;
    int page_sz = sz[0] * sz[1];
    int cur_nei_idx;
    size_t linkage_cnt = 0;
    vector<bool> checked_label(sz[2] * page_sz, false);
    FOREACH_i(idx){
        z = idx[i] / page_sz;
        remain = idx[i] - (z*page_sz);
        x = remain / sz[0];
        y = remain - x * sz[0];
        for(int j = 0; j < n_y.size(); j++){
            if(inField(y + n_y[j], x + n_x[j], z + n_z[j], sz)){
                cur_nei_idx = x + n_x[j] + (y + n_y[j])*sz[1] + (z + n_z[j])*page_sz;
                if(!checked_label[cur_nei_idx]){
                    center_idx[linkage_cnt] = idx[i];
                    nei_idx[linkage_cnt] = cur_nei_idx;
                    linkage_cnt++;
                }
                center_idx[cur_nei_idx] = true;
            }
        }
    }
    if (linkage_cnt < center_idx.size()){ // trim the unused memory
        center_idx.resize(linkage_cnt);
        nei_idx.resize(linkage_cnt);
    }
}
//int maxflow(vector<size_t> head, vector<size_t> tail, vector<float> capacity, vector<vector<size_t>> node_list,
//             vector<size_t> src_ids, vector<size_t> sink_ids, int numNodes, vector<size_t> &out_src_node_ids){

//    typedef Graph<size_t,size_t,float> GraphType;
//    GraphType *g = new GraphType(numNodes, head.size());

//    g -> add_node(numNodes);
//    FOREACH_i(src_ids){
//        g -> add_tweights( src_ids[i],  INFINITY, 0);
//    }
//    FOREACH_i(sink_ids){
//        g -> add_tweights( sink_ids[i],  0, INFINITY);
//    }
//    FOREACH_i(head){
//        g -> add_edge( head[i], tail[i],   capacity[i],   capacity[i] );
//    }

//    int flow = g -> maxflow();

//    for(int i = 0; i < numNodes; i++){
//        if (g->what_segment(i) == GraphType::SOURCE){
//            out_src_node_ids.push_back(i);
//        }
//    }
////    printf("Flow = %d\n", flow);
////    printf("Minimum cut:\n");
////    if (g->what_segment(0) == GraphType::SOURCE)
////        printf("node0 is in the SOURCE set\n");
////    else
////        printf("node0 is in the SINK set\n");
////    if (g->what_segment(1) == GraphType::SOURCE)
////        printf("node1 is in the SOURCE set\n");
////    else
////        printf("node1 is in the SINK set\n");

//    delete g;

//    return flow;
//}
/**
 * @brief regionGrow: grow the region using max-flow (combined idea of watershed and graph-cut), only voxels in fg will
 * be considered for regionGrow.
 * @param label_map: cv_32s
 * @param numCC
 * @param outLabelMap
 * @param scoreMap
 * @param fgMap
 * @param connect
 * @param cost_design
 * @param bg2sink
 */
void regionGrow(Mat *label_map, int numCC, Mat &outLabelMap, Mat *scoreMap,
                Mat *fgMap, int connect, int cost_design[], bool bg2sink){
    outLabelMap = Mat::zeros(label_map->dims, label_map->size, CV_32S);
    // for each connected component, run max-flow

    vector<size_t> arc_head_in_fg, arc_tail;
    int mat_sz[] = {label_map->size[0], label_map->size[1], label_map->size[2]};
    vector<size_t> valid_fg_idx = fgMapIdx(fgMap, CV_8U, 0);
    neighbor_idx(valid_fg_idx, arc_head_in_fg, arc_tail, mat_sz, connect);
    vector<int> arc_capacity (arc_head_in_fg.size()); //BK algorithm accept only integer capacity
    float p1, p2;
    FOREACH_i(arc_head_in_fg){
        p1 = scoreMap->at<float>(arc_head_in_fg[i]);
        p2 = scoreMap->at<float>(arc_tail[i]);
        if (cost_design[0]==ARITHMETIC_AVERAGE){
            arc_capacity[i] = (int)(pow((2/(p1+p2)), cost_design[1]) * 10000);// we keep 1/1000 accuracy
        }else if (cost_design[0]==GEOMETRIC_AVERAGE){
            if (p2==0) p2 = p1; // there should be no non-negative score, so this line should never be used
            arc_capacity[i] = (int)(pow((1/sqrt(p1*p2)), cost_design[1]) * 10000); // we keep 1/1000 accuracy
        }
    }
    vector<vector<int>> label_voxIdx;
    extractVoxIdxList(label_map, label_voxIdx, numCC, false);
    vector<size_t> sink_ids (arc_tail.size()); //valid sink nodes will not be larger than tail size
    size_t sink_ids_cnt = 0;
    typedef Graph<int,int,int> GraphType;
    GraphType *g;
    for(int i = 1; i<numCC; i++){
        // build sink_ids
        if (bg2sink){
            FOREACH_j(arc_tail){
                if (label_map->at<int>(arc_tail[j]) != i){
                    sink_ids[sink_ids_cnt] = arc_tail[j];
                    sink_ids_cnt++;
                }
            }
        }else{
            FOREACH_j(arc_tail){
                if (label_map->at<int>(arc_tail[j])>0 && label_map->at<int>(arc_tail[j]) != i){
                    sink_ids[sink_ids_cnt] = arc_tail[j];
                    sink_ids_cnt++;
                }
            }
        }
        sink_ids.resize(sink_ids_cnt);
        assert (sink_ids_cnt > 0 && label_voxIdx[i-1].size() > 0);
        // max-flow the get the grown region
        g = new GraphType(label_map->total(), (int)arc_head_in_fg.size());
        g -> add_node(label_map->total());
        FOREACH_j(label_voxIdx[i-1]){
            g -> add_tweights( label_voxIdx[i-1][j],  INFINITY, 0);
        }
        FOREACH_i(sink_ids){
            g -> add_tweights( sink_ids[i],  0, INFINITY);
        }
        FOREACH_i(arc_head_in_fg){
            g -> add_edge( arc_head_in_fg[i], arc_tail[i],   arc_capacity[i],   arc_capacity[i] );
        }

        int flow = g -> maxflow();
        FOREACH_j(valid_fg_idx){
            if (g->what_segment(valid_fg_idx[j]) == GraphType::SOURCE){
                outLabelMap.at<int>(valid_fg_idx[j]) = i;
            }
        }
    }
    delete g;
}

void setValMat(Mat &vol3d, int datatype, vector<size_t> idx, float v){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);

    if (datatype == CV_8U){
        unsigned char v0 = (unsigned char)v;
        FOREACH_i(idx){
            vol3d.at<unsigned char>(i) = v0;
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            vol3d.at<float>(i) = v;
        }
    }else if (datatype == CV_32S){
        int v0 = (int)v;
        FOREACH_i(idx){
            vol3d.at<int>(i) = v0;
        }
    }
}

void setValMat(Mat &src, int datatype, Mat *mask, float v){
    vector<size_t> idx = fgMapIdx(mask, CV_8U, 0);
    setValMat(src, datatype, idx, v);
}
/**
 * @brief gapRefine: remove voxels that not directly between two regions (tail area of the gap)
 * @param label_map
 * @param target_label0
 * @param target_label1
 * @param gap_idx
 */
void gapRefine(Mat *label_map, int target_label0, int target_label1, vector<size_t> &gap_idx){
    // 1. find the gap_idx that touched one of the region
    vector<size_t> sub_idx0;
    vector<size_t> sub_idx1;
    int im_sz[] = {label_map->size[0], label_map->size[1]};
    int x, y, z, remain;
    int page_sz = im_sz[0] * im_sz[1];

    int n_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    int n_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    FOREACH_i(gap_idx){
        if (label_map->at<int>(gap_idx[i]) == target_label0){
            z = gap_idx[i] / page_sz;
            remain = gap_idx[i] - (z*page_sz);
            x = remain / im_sz[0];
            y = remain - x * im_sz[0];

            for(int j = 0; j < 8; j++){
                if(inField(y+n_y[j], x+n_x[j], im_sz) &&
                        label_map->at<int>( gap_idx[i] + n_y[j] * im_sz[1] + n_x[j])==0){
                    sub_idx0.push_back(gap_idx[i]);
                    sub_idx0.push_back(y);
                    sub_idx0.push_back(x);
                    break;
                }
            }
        }else if(label_map->at<int>(gap_idx[i]) == target_label1){
            z = gap_idx[i] / page_sz;
            remain = gap_idx[i] - (z*page_sz);
            x = remain / im_sz[0];
            y = remain - x * im_sz[0];
            for(int j = 0; j < 8; j++){
                if(inField(y+n_y[j], x+n_x[j], im_sz) &&
                        label_map->at<int>(gap_idx[i] + n_y[j] * im_sz[1] + n_x[j])==0){
                    sub_idx1.push_back(gap_idx[i]);
                    sub_idx0.push_back(y);
                    sub_idx0.push_back(x);
                    break;
                }
            }
        }
    }
    assert(sub_idx0.size()>=6 && sub_idx1.size()>=6);
    // find the two nodes far away
    size_t dy, dx, d;
    if(sub_idx0.size() > 6){
        size_t max_dist = 0, target_i;
        for(size_t i = 3; i < sub_idx0.size(); i+=3){
            dy = sub_idx0[i+1]-sub_idx0[1];
            dx = sub_idx0[i+2]-sub_idx0[2];
            d = dy*dy + dx*dx;
            if(d > max_dist){
                max_dist = d;
                target_i = i;
            }
        }
        sub_idx0[3] = sub_idx0[target_i];
        sub_idx0[4] = sub_idx0[target_i+1];
        sub_idx0[5] = sub_idx0[target_i+2];
        sub_idx0.resize(6);
    }
    if(sub_idx1.size() > 6){
        size_t max_dist = 0, target_i;
        for(size_t i = 3; i < sub_idx1.size(); i+=3){
            dy = sub_idx1[i+1]-sub_idx1[1];
            dx = sub_idx1[i+2]-sub_idx1[2];
            d = dy*dy + dx*dx;
            if(d > max_dist){
                max_dist = d;
                target_i = i;
            }
        }
        sub_idx1[3] = sub_idx1[target_i];
        sub_idx1[4] = sub_idx1[target_i+1];
        sub_idx1[5] = sub_idx1[target_i+2];
        sub_idx1.resize(6);
    }
    // trim points that outside the quadrilateral
    // 1. get the two lines
    float k_b0[2], k_b1[2];
    if(sub_idx1[5] == sub_idx0[5]) k_b0[0] = (float)(sub_idx1[4] - sub_idx0[4]);
    else k_b0[0] = ((float)sub_idx1[4] - sub_idx0[4]) / (sub_idx1[5] - sub_idx0[5]);
    k_b0[1] = sub_idx1[4] -  k_b0[0] * sub_idx1[5];

    if(sub_idx1[2] == sub_idx0[2]) k_b1[0] = (float)(sub_idx1[1] - sub_idx0[1]);
    else k_b1[0] = ((float)sub_idx1[1] - sub_idx0[1]) / (sub_idx1[2] - sub_idx0[2]);
    k_b1[1] = sub_idx1[1] -  k_b1[0] * sub_idx1[2];

    float intersect_x;
    if(k_b0[0] == k_b1[0]) intersect_x = INFINITY;
    else intersect_x = ((k_b1[1] - k_b0[1]))/(k_b0[0] - k_b1[0]);
    //float intersect_y = k_b1[1] * intersect_x + k_b1[0];
    if(intersect_x <= MAX(sub_idx0[5], sub_idx1[5]) && intersect_x >= MIN(sub_idx0[5], sub_idx1[5])){
            //intersect_y <= MAX(sub_idx0[4], sub_idx1[4]) && intersect_x >= MIN(sub_idx0[4], sub_idx1[4]) ){
        size_t tmp = sub_idx0[0];
        sub_idx0[0] = sub_idx0[3];
        sub_idx0[3] = tmp;
        tmp = sub_idx0[1];
        sub_idx0[1] = sub_idx0[4];
        sub_idx0[4] = tmp;
        tmp = sub_idx0[2];
        sub_idx0[2] = sub_idx0[5];
        sub_idx0[5] = tmp;

//        if(sub_idx1[5] == sub_idx0[5]) k_b0[0] = (sub_idx1[4] - sub_idx0[4]);
//        else k_b0[0] = (sub_idx1[4] - sub_idx0[4]) / (sub_idx1[5] - sub_idx0[5]);
//        k_b0[1] = sub_idx1[4] -  k_b0[0] * sub_idx1[5];

//        if(sub_idx1[2] == sub_idx0[2]) k_b1[0] = (sub_idx1[1] - sub_idx0[1]);
//        else k_b1[0] = (sub_idx1[1] - sub_idx0[1]) / (sub_idx1[2] - sub_idx0[2]);
//        k_b1[1] = sub_idx1[1] -  k_b1[0] * sub_idx1[2];
    }
    // remove pionts outside the quadrilateral
    vector<Point> contour(4);
    contour[0].y = sub_idx0[1];
    contour[0].x = sub_idx0[2];
    contour[1].y = sub_idx0[4];
    contour[1].x = sub_idx0[5];
    contour[2].y = sub_idx1[4];
    contour[2].x = sub_idx1[5];
    contour[3].y = sub_idx1[1];
    contour[3].x = sub_idx1[2];
    size_t valid_idx_cnt = 0;

    FOREACH_i(gap_idx){
        if (label_map->at<int>(gap_idx[i]) == target_label0 ||
                label_map->at<int>(gap_idx[i]) == target_label1){
            gap_idx[valid_idx_cnt] = gap_idx[i];
            valid_idx_cnt ++;
        }else{
            z = gap_idx[i] / page_sz;
            remain = gap_idx[i] - (z*page_sz);
            x = remain / im_sz[0];
            y = remain - x * im_sz[0];
            if(pointPolygonTest(contour, Point2f((float)y, (float)x), false)>=0){
                gap_idx[valid_idx_cnt] = gap_idx[i];
                valid_idx_cnt ++;
            }
        }
    }
    gap_idx.resize(valid_idx_cnt);
}
/**
 * @brief extractGapVoxel: 2d gap test
 * @param label_map
 * @param numCC
 * @param gap_radius
 * @param gap_voxIdx
 */
void extractGapVoxel(Mat *label_map, Mat *fgMap, int numCC, int gap_radius,
                     vector<vector<size_t>> &gap_voxIdx, vector<bool> tested_flag){
    gap_voxIdx.resize(numCC*numCC);
    vector<size_t> valid_fg_idx = fgMapIdx(fgMap, CV_8U, 0);
    int im_sz[] = {label_map->size[0], label_map->size[1]};
    int x, y, z, remain;
    int page_sz = im_sz[0] * im_sz[1];
    int nei_seed_ids[2] = {0,0}, nei_seed_cnt = 0;
    int cur_label;
    FOREACH_i(valid_fg_idx){
        z = valid_fg_idx[i] / page_sz;
        remain = valid_fg_idx[i] - (z*page_sz);
        x = remain / im_sz[0];
        y = remain - x * im_sz[0];
        nei_seed_cnt = 0;
        for(int n_x = -gap_radius; n_x <= gap_radius; n_x++){
            for(int n_y = -gap_radius; n_y <= gap_radius; n_y++){
                if(n_x != 0 || n_y != 0){
                    if(inField(y + n_y, x + n_x, im_sz)){
                        cur_label = label_map->at<int>(valid_fg_idx[i] + n_y * im_sz[1] + n_x);
                        if(cur_label >0){
                            if (nei_seed_cnt == 0 || (nei_seed_cnt == 1 && cur_label != nei_seed_ids[0])){
                                nei_seed_ids[nei_seed_cnt++] = cur_label;
                            }else if(cur_label != nei_seed_ids[0] && cur_label != nei_seed_ids[1]){
                                nei_seed_cnt ++;
                                break;
                            }
                        }
                    }
                }
            }
            if(nei_seed_cnt > 2) break; // do not consider gap pixels among 3 or more regions
        }
        if (nei_seed_cnt == 2){
            size_t tmp_idx = (MIN(nei_seed_ids[0], nei_seed_ids[1])-1) * numCC + (MAX(nei_seed_ids[0], nei_seed_ids[1])-1);
            //if(tested_flag[tmp_idx] == false){
            gap_voxIdx[tmp_idx].push_back(valid_fg_idx[i]);
            //}
        }
    }
    FOREACH_i(gap_voxIdx){
        if(gap_voxIdx[i].size() > 0){
            if (tested_flag[i]){ // this gap does not need to test
                gap_voxIdx[i].resize(0);
            }else{
                gapRefine(label_map, i / numCC + 1, i % numCC + 1, gap_voxIdx[i]);
            }
        }
    }
}
/**
 * @brief neighbor_idx_2d: return the neighbor index of idx circle by circle
 * @param idx
 * @param fgMap
 * @param neighbor_idx_list
 * @param radius: how many circles we would like to dilate
 */
void neighbor_idx_2d(vector<size_t> idx, Mat *fgMap, vector<vector<size_t>> &neighbor_idx_list, int radius){
    Mat curr_fg;
    fgMap->copyTo(curr_fg);
    neighbor_idx_list.resize(radius);
    FOREACH_i(idx){
        curr_fg.at<unsigned char>(idx[i]) = 0;
    }
    int im_sz[] = {fgMap->size[0], fgMap->size[1]};
    int x, y, z, remain;
    int page_sz = im_sz[0] * im_sz[1];
    int n_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    int n_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    FOREACH_i(neighbor_idx_list){
        vector<size_t> curr_idx;
        size_t tmp_idx;
        if(i == 0){
            curr_idx = idx;
        }else{
            curr_idx = neighbor_idx_list[i-1];
        }
        FOREACH_j(curr_idx){
            z = curr_idx[i] / page_sz;
            remain = curr_idx[i] - (z*page_sz);
            x = remain / im_sz[0];
            y = remain - x * im_sz[0];
            for(int k = 0; k < 8; k++){
                if(inField(y+n_y[k], x+n_x[k], im_sz)){
                    vol_sub2ind(tmp_idx, y+n_y[k], x+n_x[k], z, fgMap->size);
                    if(fgMap->at<unsigned char>(tmp_idx) > 0){
                        neighbor_idx_list[i].push_back(tmp_idx);
                        fgMap->at<unsigned char>(tmp_idx) = 0;
                    }
                }
            }
        }
    }
}




/**************************Functions for debug*******************************/

void ccShowSlice3Dmat(Mat *src3d, int datatype, int slice, bool binary){
    int sz_single_frame = src3d->size[1] * src3d->size[0];
    Mat *single_slice;
    double min_v, max_v;
    if (binary){
        Mat bi_mat = *src3d > EPSILON_0;
        unsigned char *ind = (unsigned char*)bi_mat.data + sz_single_frame*slice; // sub-matrix pointer
        single_slice = new Mat(2, src3d->size, CV_8U, ind); // note, this pointer will be null out of this if
        double min_v, max_v;
        minMaxIdx(*single_slice, &min_v, &max_v);
        stringstream s_min, s_max;
        s_min << fixed << setprecision(2) << min_v;
        s_max << fixed << setprecision(2) << max_v;
        string title = "Slice:"+to_string(slice) + ", min:" + s_min.str() + ", max:" + s_max.str();
        imshow(title, *single_slice);
        waitKey(0);
        destroyWindow(title);
    }else{
        if(datatype == CV_8U){
            unsigned char *ind = (unsigned char*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
            single_slice = new Mat(2, src3d->size, CV_8U, ind);
            minMaxIdx(*single_slice, &min_v, &max_v);
        }else if(datatype == CV_16U){
            unsigned short *ind = (unsigned short*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
            single_slice = new Mat(2, src3d->size, CV_16U, ind);
            minMaxIdx(*single_slice, &min_v, &max_v);
        }else if(datatype == CV_32F){
            float *ind = (float*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
            single_slice = new Mat(2, src3d->size, CV_32F, ind);
            minMaxIdx(*single_slice, &min_v, &max_v);
            normalize(*single_slice, *single_slice, 1, 0, NORM_MINMAX);
        }else{
            single_slice = nullptr;
            qFatal("Unsupported data type!");
        }

        stringstream s_min, s_max;
        s_min << fixed << setprecision(2) << min_v;
        s_max << fixed << setprecision(2) << max_v;
        string title = "Slice:"+to_string(slice) + ", min:" + s_min.str() + ", max:" + s_max.str();
        imshow(title, *single_slice);
        waitKey(0);
        //destroyWindow(title);
        //destroyAllWindows();
    }
    delete single_slice;
}

/**
 * @brief label2rgb3d
 * @param src: CV_32S
 * @param dst: CV_8UC3, RGB data
 */
void label2rgb3d(Mat &src, Mat &dst)
{
    // Create JET colormap
    double m;
    minMaxLoc(src, nullptr, &m);
    m++;

    int n = ceil(m / 4);
    Mat1d u(n*3-1, 1, double(1.0));

    for (int i = 1; i <= n; ++i) {
        u(i-1) = double(i) / n;
        u((n*3-1) - i) = double(i) / n;
    }

    vector<double> g(n * 3 - 1, 1);
    vector<double> r(n * 3 - 1, 1);
    vector<double> b(n * 3 - 1, 1);
    for (int i = 0; i < g.size(); ++i)
    {
        g[i] = ceil(double(n) / 2) - (int(m)%4 == 1 ? 1 : 0) + i + 1;
        r[i] = g[i] + n;
        b[i] = g[i] - n;
    }

    g.erase(remove_if(g.begin(), g.end(), [m](double v){ return v > m;}), g.end());
    r.erase(remove_if(r.begin(), r.end(), [m](double v){ return v > m; }), r.end());
    b.erase(remove_if(b.begin(), b.end(), [](double v){ return v < 1.0; }), b.end());

    Mat1d cmap(m, 3, double(0.0));
    for (int i = 0; i < r.size(); ++i) { cmap(int(r[i])-1, 2) = u(i); }
    for (int i = 0; i < g.size(); ++i) { cmap(int(g[i])-1, 1) = u(i); }
    for (int i = 0; i < b.size(); ++i) { cmap(int(b[i])-1, 0) = u(u.rows - b.size() + i); }

    Mat3d cmap3 = cmap.reshape(3);

    Mat3b colormap;
    cmap3.convertTo(colormap, CV_8U, 255.0);


    // Apply color mapping
//    dst = Mat3b(src.rows, src.cols, Vec3b(0,0,0));
//    for (int r = 0; r < src.rows; ++r)
//    {
//        for (int c = 0; c < src.cols; ++c)
//        {
//            dst(r, c) = colormap(src(r,c));
//        }
//    }
}

/**
 * @brief label2rgb3d
 * @param src: CV_32S
 * @param dst: CV_8UC3, RGB data
 */
void label2rgb2d(Mat1i &src, Mat3b &dst)
{
    // Create JET colormap
    double m;
    minMaxLoc(src, nullptr, &m);
    m++;

    int n = ceil(m / 4);
    Mat1d u(n*3-1, 1, double(1.0));

    for (int i = 1; i <= n; ++i) {
        u(i-1) = double(i) / n;
        u((n*3-1) - i) = double(i) / n;
    }

    vector<double> g(n * 3 - 1, 1);
    vector<double> r(n * 3 - 1, 1);
    vector<double> b(n * 3 - 1, 1);
    for (int i = 0; i < g.size(); ++i)
    {
        g[i] = ceil(double(n) / 2) - (int(m)%4 == 1 ? 1 : 0) + i + 1;
        r[i] = g[i] + n;
        b[i] = g[i] - n;
    }

    g.erase(remove_if(g.begin(), g.end(), [m](double v){ return v > m;}), g.end());
    r.erase(remove_if(r.begin(), r.end(), [m](double v){ return v > m; }), r.end());
    b.erase(remove_if(b.begin(), b.end(), [](double v){ return v < 1.0; }), b.end());

    Mat1d cmap(m, 3, double(0.0));
    for (int i = 0; i < r.size(); ++i) { cmap(int(r[i])-1, 2) = u(i); }
    for (int i = 0; i < g.size(); ++i) { cmap(int(g[i])-1, 1) = u(i); }
    for (int i = 0; i < b.size(); ++i) { cmap(int(b[i])-1, 0) = u(u.rows - b.size() + i); }

    Mat3d cmap3 = cmap.reshape(3);

    Mat3b colormap;
    cmap3.convertTo(colormap, CV_8U, 255.0);


    // Apply color mapping
    dst = Mat3b(src.rows, src.cols, Vec3b(0,0,0));
    for (int r = 0; r < src.rows; ++r)
    {
        for (int c = 0; c < src.cols; ++c)
        {
            dst(r, c) = colormap(src(r,c));
        }
    }
}

void ccShowSliceLabelMat(Mat *src3d, int slice){
    int sz_single_frame = src3d->size[1] * src3d->size[0];
    int *ind = (int*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
    Mat single_slice(2, src3d->size, CV_32S, ind);
    Mat1i mat2show = single_slice;
    double min_v, max_v;

    minMaxIdx(mat2show, &min_v, &max_v);
    stringstream s_min, s_max;
    s_min << fixed << setprecision(2) << min_v;
    s_max << fixed << setprecision(2) << max_v;
    string title = "Slice:"+to_string(slice) + ", min:" + s_min.str() + ", max:" + s_max.str();
    Mat3b mat2display;
    label2rgb2d(mat2show, mat2display);
    imshow(title, mat2display);
    waitKey(0);
}

};
#endif // VOL_BASIC_PROC_HPP
