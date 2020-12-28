/** the basic function used in volume image processing */

#include "img_basic_proc_declare.h"

// all the template functions

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
template <typename T> T normalCDF(T x, T m, T s)
{
    T a = (x - m) / s;
    return 0.5 * erfc(-a * M_SQRT1_2); //erfc(x) = 1-erf(x), erf(-x) = -erf(x)
}
template <typename T> T normalPDF(T x, T m, T s)
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
template <typename T> T normInv(T p, T mu, T sigma)
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
template <typename T> T vec_stddev(vector<T> const & func)
{
    return sqrt(vec_variance(func));
}
template <typename T> T vec_variance(vector<T> const & func)
{
    T mean = vec_mean(func);
    T sq_sum = inner_product(func.begin(), func.end(), func.begin(), 0.0,
        [](T const & x, T const & y) { return x + y; },
        [mean](T const & x, T const & y) { return (x - mean)*(y - mean); });
    return sq_sum / ( func.size() - 1 );
}
template <typename T> T vec_mean(vector<T> const & func)
{
    return accumulate(func.begin(), func.end(), 0.0) / func.size();
}
template <typename T> vector<size_t> sort_indexes(const vector<T> &v, bool ascending, size_t start_id) {
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
    size_t n = length(fg);
    size_t m = length(bg);
    T sigma = sqrt(var_known*(n+m)/(n*m));
    return  sum_st / sigma;
}
template <typename T> void nonOV_truncatedGauss(size_t M, size_t N, T &mu, T &sigma){
        T lower_fg = normInv(1-M/(M+N));
        T f_mu, f_sigma;
        truncatedGauss(0, 1, lower_fg, INFINITY, f_mu, f_sigma);
        f_sigma = f_sigma/sqrt(M);

        T upper_nei = lower_fg;
        T n_mu, n_sigma;
        truncatedGauss(0, 1, -INFINITY, upper_nei, n_mu, n_sigma);
        n_sigma = n_sigma/sqrt(N);

        mu = f_mu - n_mu;
        sigma = sqrt(f_sigma*f_sigma + n_sigma*n_sigma);

//        T sum_st = mean_fg - mean_bg;
//        sigma = sigma * sqrt(var_known);
//        T z_debias = mu*sqrt(var_known); // mu
//        return (sum_st - z_debias) / sigma; //zscore
}

template <typename T> void OV_truncatedGauss(vector<T> fg, vector<T> bg, T &mu, T &sigma){
    size_t M = fg.size(), N = bg.size();

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
template <typename T> T orderStatsKSection(vector<T> fg, vector<T> bg, vector<T> midVals, T &mu, T &sigma){
    size_t M = fg.size();
    size_t N = bg.size();
    size_t mid_sz = midVals.size();
    size_t n = M+N+mid_sz;

    bg.insert(bg.end(), fg.begin(), fg.end());
    bg.insert(bg.end(), midVals.begin(), midVals.end());// fg: 1-M, bg, M+1-M+N, mid: M+N+1:n
    vector<size_t> sorted_id = sort_indexes(bg, false, 1);
    vector<byte> sorted_class_id(n);
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

    sigma = double(sqrt(A-B)/sqrt(n));
}

template <typename T> void vol_sub2ind(T &idx, int y, int x, int z, MatSize size){
    // assert(
    idx = z*size[0]*size[1] + x*size[0] + y;
}
template <typename T> void vol_ind2sub(T idx, int &y, int &x, int &z, MatSize size){
    z = idx / (size[0]*size[1]);
    T rmder = idx % (size[0]*size[1]);
    y = rmder/size[1];
    x = rmder-y*size[0];

}

template <typename T> void vec_sub2ind(vector<T> &idx, vector<int> y, vector<int> x, vector<int> z, MatSize size){
    size_t nPixels_slice = size[0]*size[1];
    FOREACH_i(y){
        idx[i] = z[i]*nPixels_slice + x[i]*size[0] + y[i];
    }
}
template <typename T> void vec_ind2sub(vector<T> idx, vector<int> &y, vector<int> &x, vector<int> &z, MatSize size){
    size_t nPixels_slice = size[0]*size[1];
    T rmder;
    FOREACH_i(idx){
        z[i] = idx / nPixels_slice;
        rmder = idx % nPixels_slice;
        y[i] = rmder/size[1];
        x[i] = rmder-y[i]*size[0];
    }
}
