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
    dst3d.create(src3d->dims, src3d->size, CV_32F);
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
float calVarianceStablization(const Mat *src3d, Mat & varMap, vector<float> &varTrend, float validRatio, int gap){
    Mat meanVal; //(src3d->dims, src3d->size, CV_32F, Scalar(0));
    src3d->copyTo(meanVal);
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
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)meanVal.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, meanVal.size, CV_32F, ind);
        sepFilter2D(subMatrix, subMatrix, CV_32F, kernely, kernelx, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
    //sepFilter2D(*src3d, meanVal, CV_32F, kernelx, kernely, Point(-1,-1), 0.0, BORDER_REPLICATE);
    double min_intensity, max_intensity;
//    minMaxIdx(meanVal, &min_intensity, &max_intensity);
    meanVal /= denominator;
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
        cur_level++;
    }
    if (max_intensity<20){
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
        dst3d = Mat(src3d->dims, src3d->size, CV_32S);
        //src3d->copyTo(dst3d);
        for (int z = 0; z < z_size; z++)
        {
            int *ind = (int*)dst3d.data + z * xy_size; // sub-matrix pointer
            Mat subMatrix_dst(2, dst3d.size, CV_32S, ind);
            unsigned char *ind2 = (unsigned char*)src3d->data + z * xy_size; // sub-matrix pointer
            Mat subMatrix_src(2, src3d->size, CV_8U, ind2);
            numCC += connectedComponents(subMatrix_src, subMatrix_dst, connect);
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
void extractVoxIdxList(const Mat *label3d, vector<vector<size_t>> &voxList, int numCC, bool bk_extract){
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
void extractVoxIdxList(const Mat *label3d, vector<vector<int>> &voxList, int numCC, bool bk_extract){
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
void extractVolume(const Mat *label3d, vector<size_t> &voxSzList, int numCC){
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
    //unsigned short max_val;
    if (radiusValues[2] > 0 && max_valid_z >= min_valid_z){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
        unsigned char *ind = (unsigned char*)dst3d.data + min_valid_z * xy_size; // sub-matrix pointer
        Mat *subMatrix = new Mat(2, dst3d.size, CV_8U, ind);
        for (int k = min_valid_z - 1; k >= max(0, min_valid_z - radiusValues[2]); k--){
            ind = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind);
            subMatrix->copyTo(subMatrix2);
        }

        ind = (unsigned char*)dst3d.data + max_valid_z * xy_size; // sub-matrix pointer
        subMatrix = new Mat(2, dst3d.size, CV_8U, ind);
        for (int k = max_valid_z + 1; k <= max(z_size - 1, max_valid_z + radiusValues[2]); k++){
            ind = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind);
            subMatrix->copyTo(subMatrix2);
        }
//            }
//        }
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
    printf("Multiple files saved in test.tiff\n");
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
size_t fgMapSize(Mat *src3d, int datatype, float threshold_in){
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
bool isempty(Mat *src3d, int datatype, float threshold_in){
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
bool isOnBoundary2d(Mat *fgMap, size_t idx){
    int y,x,z;
    vol_ind2sub(idx, y, x, z, fgMap->size);
    return isOnBoundary2d(fgMap, y, x, z);
}
bool isOnBoundary2d(Mat *fgMap, int y, int x, int z){
    int im_sz[] = {fgMap->size[0], fgMap->size[1]};
    for(int i=0; i < 8; i++){
        if(inField(y + n8_y[i], x + n8_x[i], im_sz)
                && fgMap->at<unsigned char>(y + n8_y[i], x + n8_x[i], z) == 0){
            return true;
        }
    }
    return false;
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
                cur_nei_idx = y + n_y[j] + (x + n_x[j])*sz[0] + (z + n_z[j])*page_sz;
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
                        label_map->at<int>(y+n_y[j], x+n_x[j], z)==0){
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
                        label_map->at<int>(y+n_y[j], x+n_x[j], z)==0){
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
                        cur_label = label_map->at<int>(y+n_y, x+n_x, z);
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

void ccShowSlice3Dmat(Mat *src3d, int datatype, int slice){
    int sz_single_frame = src3d->size[1] * src3d->size[0];
    Mat *single_slice;
    if(datatype == CV_8U){
        unsigned char *ind = (unsigned char*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
        single_slice = new Mat(2, src3d->size, CV_8U, ind);
    }else if(datatype == CV_16U){
        unsigned short *ind = (unsigned short*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
        single_slice = new Mat(2, src3d->size, CV_16U, ind);
    }else if(datatype == CV_32F){
        float *ind = (float*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
        single_slice = new Mat(2, src3d->size, CV_32F, ind);
        normalize(*single_slice, *single_slice, 1, 0, NORM_MINMAX);
    }else{
        single_slice = nullptr;
        qFatal("Unsupported data type!");
    }
    imshow("Slice:"+to_string(slice), *single_slice);
    waitKey(0);
    destroyWindow("Slice:"+to_string(slice));
    //destroyAllWindows();
}
