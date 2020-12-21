#ifndef IMG_BASIC_PROC_H
#define IMG_BASIC_PROC_H

/** the basic function used in volume image processing */

#include "img_basic_proc_declare.h"

enum filterDirection{DIRECTION_X = 0, DIRECTION_Y, DIRECTION_Z};
/**
 * @brief principalCv2d
 * @param src3d: float mat
 * @param dst3d: float mat
 */
void principalCv2d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity = 0){
    src3d->copyTo(dst3d);
    gaussianSmooth3Ddata(dst3d, sigma);

    int x_size  = dst3d.size[0];
    int y_size  = dst3d.size[1];
    int z_size  = dst3d.size[2];
    size_t xy_size = x_size*y_size;

    // Filter XY dimensions for every Z
    Mat kernelx = (Mat_<float>(1,3)<<-0.5, 0, 0.5);
    Mat kernely = (Mat_<float>(3,1)<<-0.5, 0, 0.5);

    //int mat_sizes[] = {2,2};
    Mat mat2x2(2, (int []){2,2}, CV_32F, Scalar(0));
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
void principalCv3d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity = 0){
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
    //int mat_sizes[] = {2,2};
    Mat mat3x3(2, (int []){3,3}, CV_32F, Scalar(0));
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
void gaussianSmooth3Ddata(Mat &data4smooth, const float *sigma)
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
//                    }data4smooth
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
 * @brief filterZdirection
 * @param src3d : float
 * @param dst3d : float
 * @param kernel_z : float (border is replicated)
 */
void filterZdirection(Mat* src3d, Mat &dst3d, Mat kernel_z){
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
    if (abs(alpha) > 1e7) beta = 1e7;

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
        vector<float> tmpVals = vec_atrange(vals4var, ub, lb);
        float sigmaTr = vec_stddev(tmpVals);
        //sigmaTr = sqrt(mean(vals(vals<ub & vals>lb).^2));
        float varCorrectionTerm = truncatedGauss(mu, sigma, lb, ub, t_mu, t_sigma); // only lb and ub is usefull if we only want the third output
        sigma = sigmaTr/sqrt(varCorrectionTerm);
        mu = vec_mean(tmpVals);
    }
}

/**
 * @brief calVarianceStablization
 * @param src3d: float mat 3d
 * @param validRatio: consider the validRatio pixels in the data (remove pixels with saturated intensity)
 * @param gap
 */
float calVarianceStablization(Mat* src3d, Mat & varMap, float &varTrend, float validRatio = 0.95, int gap=2){
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
    float *xx = new float[levels];
    for (int i = 0; i < levels; i++) xx[0] = 0;
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
        xx[cur_level] = (denominator/(denominator+1))*varCur;
        testedElements += numElements[cur_level].size();
        if (cur_level < target_level &&
                (testedElements/validElements > validRatio)){
                //if we only want to use a fixed point to reprsent all the variance
                target_level = cur_level;
        }
    }
    if (max_intensity>20){
        for (int i=0;i<3;i++) xx[i] = 0;
    }
    else{
        for (int i=0;i<15;i++) xx[i] = 0;
    }
    for(size_t i = 0; i < xy_size*z_size; i++){
        cur_level = (int) ((meanVal.at<float>(i) - min_intensity)/unit_intensity);
        if (cur_level >= levels) cur_level = levels - 1;
        if (cur_level < 0) cur_level = 0;
        varMap.at<float>(i) = xx[cur_level];
    }

    varTrend = *xx;
    return xx[target_level];
}

//double normalCDF(double value) //standard
//{
//   return 0.5 * erfc(-value * M_SQRT1_2);
//}
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
    return 0.5 * erfc(-a * M_SQRT1_2); //erfc(x) = 1-erf(x)
}
template <typename T> T normalPDF(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-T(0.5) * a * a);
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
/**
 * @brief ConnectedComponents3d
 * @param src3d
 * @param dst3d
 * @param connect: 4, 8 for 2d and 6,10,26 for 3d
 */
int ConnectedComponents3d(Mat* src3d, Mat &dst3d, int connect){
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
                        if (src3d->at<bool>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<bool>(i+1,j,k)>0){
                            dst3d.at<bool>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<bool>(i,j+1,k)>0){
                            dst3d.at<bool>(i,j+1,k) = current_id;
                        }
                        if ((k + 1) < z_size && src3d->at<bool>(i,j,k+1)>0){
                            dst3d.at<bool>(i,j,k+1) = current_id;
                        }
                    }
                }
            }
        }

        if(connect == 10){
            for (int i = 0; i < x_size; i++){
                for (int j = 0; j< y_size; j++){
                    for (int k = 0; k < z_size; k++){
                        if (src3d->at<bool>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<bool>(i+1,j,k)>0){
                            dst3d.at<bool>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<bool>(i,j+1,k)>0){
                            dst3d.at<bool>(i,j+1,k) = current_id;
                        }
                        if ((k + 1) < z_size && src3d->at<bool>(i,j,k+1)>0){
                            dst3d.at<bool>(i,j,k+1) = current_id;
                        }

                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<bool>(i+1,j+1,k)>0){
                            dst3d.at<bool>(i+1,j+1,k) = current_id;
                        }
                    }
                }
            }
        }

        if(connect == 26){
            for (int i = 0; i < x_size; i++){
                for (int j = 0; j< y_size; j++){
                    for (int k = 0; k < z_size; k++){
                        if (src3d->at<bool>(i,j,k) == 0){
                            continue;
                        }
                        if (dst3d.at<int>(i,j,k) == 0){
                            numCC ++;
                            dst3d.at<int>(i,j,k) = numCC;
                        }
                        current_id = dst3d.at<int>(i,j,k);
                        if ((i + 1) < x_size && src3d->at<bool>(i+1,j,k)>0){
                            dst3d.at<bool>(i+1,j,k) = current_id;
                        }
                        if ((j + 1) < y_size && src3d->at<bool>(i,j+1,k)>0){
                            dst3d.at<bool>(i,j+1,k) = current_id;
                        }
                        if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<bool>(i+1,j+1,k)>0){
                            dst3d.at<bool>(i+1,j+1,k) = current_id;
                        }

                        if ((k + 1) < z_size){
                            if (src3d->at<bool>(i,j,k+1)>0){
                                dst3d.at<bool>(i,j,k+1) = current_id;
                            }
                            if ((j + 1) < y_size && src3d->at<bool>(i,j+1,k+1)>0){
                                dst3d.at<bool>(i,j+1,k+1) = current_id;
                            }
                            if ((i + 1) < x_size && src3d->at<bool>(i+1,j,k+1)>0){
                                dst3d.at<bool>(i+1,j,k+1) = current_id;
                            }
                            if ((i + 1) < x_size && (j + 1) < y_size && src3d->at<bool>(i+1,j+1,k+1)>0){
                                dst3d.at<bool>(i+1,j+1,k+1) = current_id;
                            }
                        }
                    }
                }
            }
        }
    }

    return numCC;
}
#endif // IMG_BASIC_PROC_H


