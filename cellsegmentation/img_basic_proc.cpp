/** the basic function used in volume image processing */

#include "img_basic_proc.h"

/**
 * @brief principalCv2d
 * @param src3d: float mat
 * @param dst3d: float mat
 */
void principalCv2d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity = 0){
    src3d->copyTo(dst3d);
    gaussianSmooth3Ddata(dst3d, sigma);

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
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_32F, ind);
        Mat lx, ly, lxx, lxy, lyy;

        filter2D(subMatrix, lx, -1, kernelx);
        filter2D(subMatrix, ly, -1, kernely);
        filter2D(lx, lxy, -1, kernely);
        filter2D(lx, lxx, -1, kernelx);
        //filter2D(ly, lyx, -1, kernelx); // the same as lxy
        filter2D(ly, lyy, -1, kernelx);
        for (int r = 0; r < x_size; r ++){
            for (int c = 0; c < y_size; c ++){
                if (src3d->at<float>(r,c,z) <= minIntensity){
                    dst3d.at<float>(r,c,z) = 0;
                    continue;
                }
                mat2x2.at<float>(0,0) = lxx.at<float>(r,c);
                mat2x2.at<float>(0,1) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,0) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,1) = lyy.at<float>(r,c);
                eigen(mat2x2, eigen_values);
                dst3d.at<float>(r,c,z) = eigen_values.at<float>(0); // the largest eigen value
            }
        }
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernel_row.t(), kernel_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
}
/**
 * @brief principalCv3d
 * @param src3d: float mat
 * @param dst3d: float mat
 * @param sigma
 * @param minIntensity
 */
void principalCv3d(const Mat* src3d, Mat &dst3d, float sigma[], int minIntensity = 0){
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
    filterZdirection(&Dz, Dzz, kernelz);
    int mat_sizes[] = {2,2};
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
        Mat lzz(2, dst3d.size, CV_32F, (float*)Dzz.data + z * xy_size);
        for (int r = 0; r < x_size; r ++){
            for (int c = 0; c < y_size; c ++){
                if (src3d->at<float>(r,c,z) <= minIntensity){
                    dst3d.at<float>(r,c,z) = 0;
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
                dst3d.at<float>(r,c,z) = eigen_values.at<float>(0); // the largest eigen value
            }
        }
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernel_row.t(), kernel_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
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

    int kernDimension = 2*ceil(2*sigma[0])+1;
    Mat kernel_row = getGaussianKernel(kernDimension, sigma[0], CV_32F);
    kernDimension = 2*ceil(2*sigma[1])+1;
    Mat kernel_col = getGaussianKernel(kernDimension, sigma[1], CV_32F);
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)data4smooth.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, data4smooth.size, CV_32F, ind);
        sepFilter2D(subMatrix, subMatrix, CV_32F, kernel_row.t(), kernel_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
    if (sigma[2] > 0){
        kernDimension = 2*ceil(2*sigma[2])+1;
        Mat kernel_z = getGaussianKernel(kernDimension, sigma[2], CV_32F);
        // Filter Z dimension
        filterZdirection(&data4smooth, data4smooth, kernel_z);
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
void filterVolume(const Mat* src3d, Mat &dst3d, Mat kernel, unsigned direction){
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
 * @brief filterZdirection
 * @param src3d : float
 * @param dst3d : float
 * @param kernel_z : float (border is replicated)
 */
void filterZdirection(const Mat* src3d, Mat &dst3d, Mat kernel_z){
    assert(src3d->dims == 3);
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
    return sigma;
}

/**
 * @brief calVarianceStablization
 * @param src3d: float mat 3d
 * @param validRatio: consider the validRatio pixels in the data (remove pixels with saturated intensity)
 * @param gap
 */
float calVarianceStablization(const Mat* src3d, Mat & varMap, vector<float> &varTrend, float validRatio = 0.95, int gap=2){
    Mat meanVal(src3d->dims, src3d->size, CV_32F, Scalar(0));
    Mat kernelx, kernely, kernelz;
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
        denominator = 24;
        //kernelz = (Mat_<float>(1,7)<<1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }
    sepFilter2D(*src3d, meanVal, CV_32F, kernelx, kernely, Point(-1,-1), 0.0, BORDER_REPLICATE);
    meanVal /= denominator;
    Mat diff = *src3d - meanVal;
    int levels = 200;
    double min_intensity, max_intensity;
    minMaxLoc(meanVal, &min_intensity, &max_intensity);

    double unit_intensity = (max_intensity - min_intensity) / (levels-1);
    varTrend.resize(levels);
    //for (int i = 0; i < levels; i++) xx[0] = 0;
    fill(varTrend.begin(), varTrend.end(), 0);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;

    vector<vector<long>> numElements (levels);
    long testedElements = 0, validElements = 0;
    int cur_level = 0;
    for (int i = gap+1; i < x_size-gap; i++){
        for (int j = gap+1; j< y_size-gap; j++){
            for (int k = 0; k < z_size; k++){
                cur_level = int((meanVal.at<float>(i,j,k) - min_intensity) / unit_intensity);
                if(cur_level > 0){
                    validElements ++;
                    numElements[cur_level].push_back(j * x_size + i + k * xy_size);
                }
            }
        }
    }
    int target_level = 0;
    cur_level = 1;
    int nei_vec_step = 1;

    while(cur_level < levels){
        nei_vec_step = 1;
        while (numElements[cur_level].size() < 100){
            if (cur_level - nei_vec_step > 0){
                numElements[cur_level].insert(numElements[cur_level].end(),
                                              numElements[cur_level-1].begin(),
                                            numElements[cur_level-1].end());
            }
            if (cur_level + nei_vec_step < levels){
                numElements[cur_level].insert(numElements[cur_level].end(),
                                              numElements[cur_level+1].begin(),
                                            numElements[cur_level+1].end());
            }
            nei_vec_step ++;
        }
        vector<float> vals4var (numElements[cur_level].size());
        for (size_t j = 0; j < numElements[cur_level].size(); j++){
            vals4var[j] = diff.at<float>(numElements[cur_level][j]);
        }
        float varCur = varByTruncate(vals4var, 2, 3);
        varTrend[cur_level] = (denominator/(denominator+1))*varCur;
        testedElements += numElements[cur_level].size();
        if (cur_level < target_level &&
                (testedElements/validElements > validRatio)){
                //if we only want to use a fixed point to reprsent all the variance
                target_level = cur_level;
        }
    }
    if (max_intensity>20){
        for (int i=0;i<3;i++) varTrend[i] = 0;
    }
    else{
        for (int i=0;i<15;i++) varTrend[i] = 0;
    }
    for(size_t i = 0; i < xy_size*z_size; i++){
        cur_level = (int) ((meanVal.at<float>(i) - min_intensity)/unit_intensity);
        if (cur_level >= levels) cur_level = levels - 1;
        if (cur_level < 0) cur_level = 0;
        varMap.at<float>(i) = varTrend[cur_level];
    }

    //varTrend = *xx;
    return varTrend[target_level];
}

//double normalCDF(double value) //standard
//{
//   return 0.5 * erfc(-value * M_SQRT1_2);
//}
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
    assert (p >= 0 && p <= 1)
    assert (sigma >= 0);

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
        sigma = sqrt(f_sigma^2 + n_sigma^2);

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
/**
 * @brief connectedComponents3d
 * @param src3d: boolean but represented by CV_8U
 * @param dst3d
 * @param connect: 4, 8 for 2d and 6,10,26 for 3d
 */
int connectedComponents3d(const Mat* src3d, Mat &dst3d, int connect){
    //assert(src3d->type()==bool);
    assert(src3d->dims == 3);
    assert(connect == 4 || connect == 8 || connect == 6 || connect == 10 || connect == 26);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    size_t numCC = 0;
    if (connect == 4 || connect == 8){
        src3d->copyTo(dst3d);
        for (int z = 0; z < z_size; z++)
        {
            float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
            Mat subMatrix(2, dst3d.size, CV_32F, ind);
            numCC += connectedComponents(subMatrix, subMatrix, connect);
        }
    }else{
        dst3d = Mat::zeros(src3d->dims, src3d->size, CV_32S);
        int current_id = 0;
        if(connect == 6){
            for (int i = 0; i < x_size; i++){
                for (int j = 0; j< y_size; j++){
                    for (int k = 0; k < z_size; k++){
                        if (src3d->at<unsigned short>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<unsigned short>(i+1,j,k)>0){
                            dst3d.at<unsigned short>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<unsigned short>(i,j+1,k)>0){
                            dst3d.at<unsigned short>(i,j+1,k) = current_id;
                        }
                        if ((k + 1) < z_size && src3d->at<unsigned short>(i,j,k+1)>0){
                            dst3d.at<unsigned short>(i,j,k+1) = current_id;
                        }
                    }
                }
            }
        }

        if(connect == 10){
            for (int i = 0; i < x_size; i++){
                for (int j = 0; j< y_size; j++){
                    for (int k = 0; k < z_size; k++){
                        if (src3d->at<unsigned short>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<unsigned short>(i+1,j,k)>0){
                            dst3d.at<unsigned short>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<unsigned short>(i,j+1,k)>0){
                            dst3d.at<unsigned short>(i,j+1,k) = current_id;
                        }
                        if ((k + 1) < z_size && src3d->at<unsigned short>(i,j,k+1)>0){
                            dst3d.at<unsigned short>(i,j,k+1) = current_id;
                        }

                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned short>(i+1,j+1,k)>0){
                            dst3d.at<unsigned short>(i+1,j+1,k) = current_id;
                        }
                    }
                }
            }
        }

        if(connect == 26){
            for (int i = 0; i < x_size; i++){
                for (int j = 0; j< y_size; j++){
                    for (int k = 0; k < z_size; k++){
                        if (src3d->at<unsigned short>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<unsigned short>(i+1,j,k)>0){
                            dst3d.at<unsigned short>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<unsigned short>(i,j+1,k)>0){
                            dst3d.at<unsigned short>(i,j+1,k) = current_id;
                        }
                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned short>(i+1,j+1,k)>0){
                            dst3d.at<unsigned short>(i+1,j+1,k) = current_id;
                        }

                        if ((k + 1) < z_size){
                            if (src3d->at<unsigned short>(i,j,k+1)>0){
                                dst3d.at<unsigned short>(i,j,k+1) = current_id;
                            }
                            if ((j + 1) < y_size && src3d->at<unsigned short>(i,j+1,k+1)>0){
                                dst3d.at<unsigned short>(i,j+1,k+1) = current_id;
                            }
                            if ((i + 1) < x_size && src3d->at<unsigned short>(i+1,j,k+1)>0){
                                dst3d.at<unsigned short>(i+1,j,k+1) = current_id;
                            }
                            if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<unsigned short>(i+1,j+1,k+1)>0){
                                dst3d.at<unsigned short>(i+1,j+1,k+1) = current_id;
                            }
                        }
                    }
                }
            }
        }
    }

    return numCC;
}
int floatMap2idMap(Mat* src3d, Mat &dst3d, int connect){
    assert(src3d->dims == 3);
    assert(connect == 6 || connect == 10 || connect == 26);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    //size_t xy_size = x_size*y_size;
    size_t numCC = 0;

    dst3d = Mat::zeros(src3d->dims, src3d->size, CV_32S);
    int current_id = 0;
    if(connect == 6){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    float curr_val = src3d->at<float>(i,j,k);
                    if (curr_val < 1e-5){
                        continue;
                    }
                    if (dst3d.at<int>(i,j,k) == 0){
                        numCC ++;
                        dst3d.at<int>(i,j,k) = numCC;
                    }
                    current_id = dst3d.at<int>(i,j,k);
                    if ((i + 1) < x_size && (src3d->at<float>(i+1,j,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i+1,j,k) = current_id;
                    }
                    if ((j + 1) < y_size && (src3d->at<float>(i,j+1,k)- curr_val) < 1e-5){
                        dst3d.at<int>(i,j+1,k) = current_id;
                    }
                    if ((k + 1) < z_size && (src3d->at<float>(i,j,k+1)- curr_val) < 1e-5){
                        dst3d.at<int>(i,j,k+1) = current_id;
                    }
                }
            }
        }
    }

    if(connect == 10){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    float curr_val = src3d->at<float>(i,j,k);
                    if (curr_val < 1e-5){
                        continue;
                    }
                    if (dst3d.at<int>(i,j,k) == 0){
                        numCC ++;
                        dst3d.at<int>(i,j,k) = numCC;
                    }
                    current_id = dst3d.at<int>(i,j,k);
                    if ((i + 1) < x_size && (src3d->at<float>(i+1,j,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i+1,j,k) = current_id;
                    }
                    if ((j + 1) < y_size && (src3d->at<float>(i,j+1,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i,j+1,k) = current_id;
                    }
                    if ((k + 1) < z_size && (src3d->at<float>(i,j,k+1) - curr_val) < 1e-5){
                        dst3d.at<int>(i,j,k+1) = current_id;
                    }

                    if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(i+1,j+1,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i+1,j+1,k) = current_id;
                    }
                }
            }
        }
    }

    if(connect == 26){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    float curr_val = src3d->at<float>(i,j,k);
                    if (curr_val < 1e-5){
                        continue;
                    }
                    if (dst3d.at<int>(i,j,k) == 0){
                        numCC ++;
                        dst3d.at<int>(i,j,k) = numCC;
                    }
                    current_id = dst3d.at<int>(i,j,k);
                    if ((i + 1) < x_size && (src3d->at<float>(i+1,j,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i+1,j,k) = current_id;
                    }
                    if ((j + 1) < y_size && (src3d->at<float>(i,j+1,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i,j+1,k) = current_id;
                    }
                    if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(i+1,j+1,k) - curr_val) < 1e-5){
                        dst3d.at<int>(i+1,j+1,k) = current_id;
                    }

                    if ((k + 1) < z_size){
                        if ((src3d->at<float>(i,j,k+1) - curr_val) < 1e-5){
                            dst3d.at<int>(i,j,k+1) = current_id;
                        }
                        if ((j + 1) < y_size && (src3d->at<float>(i,j+1,k+1) - curr_val) < 1e-5){
                            dst3d.at<int>(i,j+1,k+1) = current_id;
                        }
                        if ((i + 1) < x_size && (src3d->at<float>(i+1,j,k+1) - curr_val) < 1e-5){
                            dst3d.at<int>(i+1,j,k+1) = current_id;
                        }
                        if ((i + 1) < x_size && (j + 1) < y_size && (src3d->at<float>(i+1,j+1,k+1) - curr_val) < 1e-5){
                            dst3d.at<int>(i+1,j+1,k+1) = current_id;
                        }
                    }
                }
            }
        }
    }

    return numCC;
}
/**
 * @brief extractVoxList
 * @param label3d
 * @param voxList
 * @param numCC
 */
void extractVoxIdxList(const Mat *label3d, vector<vector<size_t>> &voxList, int numCC){
    voxList.resize(numCC);
    for (size_t i = 0; i < label3d->total(); i++){
        if(label3d->at<int>(i) > 0){
            voxList[label3d->at<int>(i)].push_back(i);
        }
    }
}

/**
 * @brief removeSmallCC
 * @param label3d
 * @param numCC
 * @param min_size
 * @param relabel
 */
void removeSmallCC(Mat &label3d, int &numCC, size_t min_size, bool relabel = true){
    vector<size_t> cc_size(numCC);
    fill(cc_size.begin(), cc_size.end(), 0);
    for (size_t i = 0; i < label3d.total(); i++){
        if(label3d.at<int>(i) > 0){
            cc_size[label3d.at<int>(i) - 1]++;
        }
    }
    if(relabel){
        vector<size_t> newlabel(numCC);
        fill(newlabel.begin(), newlabel.end(), 0);
        int rm_cc_cnt = 0;
        for (int i = 0; i < numCC; i ++){
            if (cc_size[i] >= min_size){
                newlabel[i] = rm_cc_cnt;
                rm_cc_cnt ++ ;
            }
        }
        numCC = rm_cc_cnt;
        for (size_t i = 0; i < label3d.total(); i++){
            if(label3d.at<int>(i) > 0){
                label3d.at<int>(i) = newlabel[label3d.at<int>(i) - 1];
            }
        }
    }else{
        // no change to numCC, only set invalid CC to 0
        for (size_t i = 0; i < label3d.total(); i++){
            if(label3d.at<int>(i) > 0 && cc_size[label3d.at<int>(i) - 1] < min_size){
                label3d.at<int>(i) = 0;
            }
        }
    }
}

/**
 * @brief volumeDilate: x,y direction is formal, while z direction is simplified
 * @param src3d: boolean (but represented by CV_8U)
 * @param dst3d
 * @param kernel
 */
void volumeDilate(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type){
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
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_8U, ind);
        dilate(subMatrix, subMatrix, element);
    }
    unsigned short max_val;
    if (radiusValues[2] > 0){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    max_val = src3d->at<unsigned short>(i,j,k);
                    for (int kk = 1; kk < radiusValues[2]; kk++){
                        if (k+kk < z_size){
                            max_val = max(src3d->at<unsigned short>(i,j,k+kk), max_val);
                        }
                        if (k-kk >= 0){
                            max_val = max(src3d->at<unsigned short>(i,j,k-kk), max_val);
                        }
                    }
                    dst3d.at<unsigned short>(i,j,k) = max_val;
                }
            }
        }
    }
}

/**
 * @brief volumeErode: x,y direction is formal, while z direction is fixed to morph_corss
 * @param src3d: boolean (but represented by CV_8U)
 * @param dst3d
 * @param kernel
 */
void volumeErode(const Mat *src3d, Mat &dst3d, int *radiusValues, int dilation_type){
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
    unsigned short min_val;
    if (radiusValues[2] > 0){
        for (int i = 0; i < x_size; i++){
            for (int j = 0; j< y_size; j++){
                for (int k = 0; k < z_size; k++){
                    min_val = src3d->at<unsigned short>(i,j,k);
                    for (int kk = 1; kk < radiusValues[2]; kk++){
                        if (k+kk < z_size){
                            min_val = min(src3d->at<unsigned short>(i,j,k+kk), min_val);
                        }
                        if (k-kk >= 0){
                            min_val = min(src3d->at<unsigned short>(i,j,k-kk), min_val);
                        }
                    }
                    dst3d.at<unsigned short>(i,j,k) = min_val;
                }
            }
        }
    }
}

void volumeWrite(Mat *src3d, string filename){
    vector<Mat> imgs;
    split(*src3d, imgs);
    imwrite(filename, imgs);
}



