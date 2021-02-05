/** the basic function used in volume image processing */

#include "img_basic_proc.h"

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
    //ccShowSlice3Dmat(dst3d, CV_32F);
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
        filter2D(ly, lyy, -1, kernely);
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
    //volumeWrite(&data4smooth, "/home/ccw/Desktop/1.tiff");
    //ccShowSlice3Dmat(data4smooth, CV_32F);
    if (sigma[2] > 0){
        kernDimension = ceil(2*sigma[2])+1;
        Mat kernel_z = getGaussianKernel(kernDimension, sigma[2], CV_32F);
        // Filter Z dimension
        Mat copyMat;
        data4smooth.copyTo(copyMat);
        filterZdirection(&copyMat, data4smooth, kernel_z);

        //ccShowSlice3Dmat(data4smooth, CV_32F);
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
            float* z_ptr = (float*)src3d->data + x * y_size + y;//same as data4smooth.ptr<float>(0, y, x)
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
            z_ptr = (float*)dst3d.data + x * y_size + y;
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
        if (lb < ub){
            vector<float> tmpVals = vec_atrange(vals4var, ub, lb, true);
            float sigmaTr = vec_stddev(tmpVals);
            //sigmaTr = sqrt(mean(vals(vals<ub & vals>lb).^2));
            float varCorrectionTerm = truncatedGauss(mu, sigma, lb, ub, t_mu, t_sigma); // only lb and ub is usefull if we only want the third output
            sigma = sigmaTr/sqrt(varCorrectionTerm);
            mu = vec_mean(tmpVals);
        }else{
            return 0;
        }
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
    //ccShowSlice3Dmat(&kernel, CV_32F);
    int x_size  = src3d->size[0];
    int y_size  = src3d->size[1];
    int z_size  = src3d->size[2];
    size_t xy_size = x_size*y_size;
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)meanVal.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, meanVal.size, CV_32F, ind);
        //ccShowSlice3Dmat(src3d, CV_32F);
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernely, kernelx, Point(-1,-1), 0.0, BORDER_REPLICATE);
        filter2D(subMatrix, subMatrix, -1, kernel);
        //ccShowSlice3Dmat(&subMatrix, CV_32F);
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
        if (cur_level_idx.size() < 100){ // valid pionts are limited
            if (cur_level > 0) varTrend[cur_level] = varTrend[cur_level-1];
        }else{
            vector<float> vals4var (cur_level_idx.size());
            for (size_t j = 0; j < cur_level_idx.size(); j++){
                vals4var[j] = diff.at<float>(cur_level_idx[j]);
            }
            float varCur = varByTruncate(vals4var, 2, 3);
            varTrend[cur_level] = (denominator/(denominator+1))*varCur;
        }
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
    //ccShowSlice3Dmat(&meanVal, CV_32F, 3);
    for(size_t i = 0; i < xy_size*z_size; i++){
        if (meanVal.at<float>(i) > min_valid_intensity){
            cur_level = (int) ((meanVal.at<float>(i) - min_intensity)/unit_intensity);
            if (cur_level >= levels) cur_level = levels - 1;
            if (cur_level < 0) cur_level = 0;
            varMap.at<float>(i) = varTrend[cur_level];
        }
    }
    //ccShowSlice3Dmat(&varMap, CV_32F, 3);
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
    if(dst3d.empty()){
        dst3d = Mat(src3d->dims, src3d->size, CV_32S);
    }
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
    if (dst3d.empty()) dst3d.create(src3d->dims, src3d->size, CV_32S);
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
    int max_v = min(bound, *max_element(idx_sub.begin(), idx_sub.end()) + shift+1);
    out_range = Range(min_v, max_v);
}
void regionAvgIntensity(Mat* src3dFloatData, Mat* src3dIdMap, int numCC, vector<float> &avgIntensities){
    vector<size_t> region_sz(numCC);
    avgIntensities.resize(numCC);
    fill(region_sz.begin(), region_sz.end(), 0);
    fill(avgIntensities.begin(), avgIntensities.end(), 0);
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
void regionAvgIntensity(Mat* src3dFloatData, vector<vector<size_t>> &voxList, vector<float> &avgIntensities){
    avgIntensities.resize(voxList.size());
    fill(avgIntensities.begin(), avgIntensities.end(), 0);
    for(size_t i = 0; i<voxList.size(); i++){
        for(size_t j = 0; j<voxList[i].size(); j++){
            avgIntensities[i] += src3dFloatData->at<float>(voxList[i][j]);
        }
    }
    FOREACH_i(voxList){
        avgIntensities[i] /= voxList[i].size();
    }
}
/**
 * @brief extractVoxList: return the voxel idx of all connected component
 * @param label3d
 * @param voxList
 * @param numCC
 * @param bk_extract: background also count as a region
 */
void extractVoxIdxList(Mat *label3d, vector<vector<size_t>> &voxList, int numCC, bool bk_extract){
    voxList.resize(numCC);
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
void extractVoxIdxList(Mat *label3d, vector<vector<int>> &voxList, int numCC, bool bk_extract){
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
void extractVoxIdxList(vector<Mat> &label3d, vector<vector<size_t>> &voxList, vector<size_t> numCC){
    long total_CC = accumulate(numCC.begin(), numCC.end(), 0);
    voxList.resize(total_CC);
    long increment = 0;
    FOREACH_i(label3d){
        increment += i > 0 ? numCC[i-1]:0;
        FOREACH_j_MAT(label3d[i]){
            if(label3d[i].at<int>(j) > 0){
                voxList[label3d[i].at<int>(j)-1 + increment].push_back(j);
            }
        }
    }
}
void extractVoxIdxGivenId(Mat *label3d, vector<size_t> &voxList, int id){
    for (size_t i = 0; i < label3d->total(); i++){
        if(label3d->at<int>(i) == id){
            voxList.push_back(i);
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
            newlabel[i] = rm_cc_cnt + 1;
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
            if(subMatrix.at<unsigned char>(i) > 0){
                if(min_valid_z > z ){
                    min_valid_z = z;
                }
                if(max_valid_z < z){
                    max_valid_z = z;
                }
                valid_z = true;
                break;
            }
        }
        if(valid_z){
            //ccShowSlice3Dmat(&subMatrix, CV_8U);
            dilate(subMatrix, subMatrix, element);
            //ccShowSlice3Dmat(&subMatrix, CV_8U);
//            unsigned char *ind_up = (unsigned char*)dst3d.data; // sub-matrix pointer
//            Mat subMatrix_up(2, dst3d.size, CV_8U, ind_up);
//            ccShowSlice3Dmat(&subMatrix_up, CV_8U);
//            subMatrix.copyTo(subMatrix_up);
//            ccShowSlice3Dmat(&subMatrix_up, CV_8U);
        }

    }
    //unsigned char max_val;
    if (radiusValues[2] > 0 && max_valid_z >= min_valid_z){
//        for (int i = 0; i < x_size; i++){
//            for (int j = 0; j< y_size; j++){
        unsigned char *ind_up = (unsigned char*)dst3d.data + min_valid_z * xy_size; // sub-matrix pointer
        Mat subMatrix_up(2, dst3d.size, CV_8U, ind_up);
        //ccShowSlice3Dmat(&subMatrix_up, CV_8U);
        for (int k = min_valid_z - 1; k >= max(0, min_valid_z - radiusValues[2]); k--){
            ind_up = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind_up);
            subMatrix_up.copyTo(subMatrix2);
        }
        //ccShowSlice3Dmat(&subMatrix_up, CV_8U);
        unsigned char *ind_down = (unsigned char*)dst3d.data + max_valid_z * xy_size; // sub-matrix pointer
        Mat subMatrix_down(2, dst3d.size, CV_8U, ind_down);
        //ccShowSlice3Dmat(&subMatrix_down, CV_8U);
        for (int k = max_valid_z + 1; k <= min(z_size - 1, max_valid_z + radiusValues[2]); k++){
            ind_down = (unsigned char*)dst3d.data + k * xy_size; // sub-matrix pointer
            Mat subMatrix2(2, dst3d.size, CV_8U, ind_down);
            subMatrix_down.copyTo(subMatrix2);
        }
        //ccShowSlice3Dmat(&subMatrix_down, CV_8U);
        //ccShowSlice3Dmat(&dst3d, CV_8U);
    }

}

//// This volume dilation and erosion algorithm should be implemented using distance transform with
/// O(n) complexity. (Though now is already O(n), the constant part is large.)
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
    //if(dst3d.empty()) dst3d = Mat(src3d->dims, src3d->size, CV_8U);
    src3d->copyTo(dst3d);
    //ccShowSlice3Dmat(&dst3d, CV_8U);
    Mat element = getStructuringElement( dilation_type,
                           Size( 2*radiusValues[0] + 1, 2*radiusValues[1]+1 ),
                           Point( radiusValues[0], radiusValues[1] ) );
    for (int z = 0; z < z_size; z++)
    {
        unsigned char *ind = (unsigned char*)dst3d.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, dst3d.size, CV_8U, ind);
        erode(subMatrix, subMatrix, element);
    }
    //ccShowSlice3Dmat(&dst3d, CV_8U);
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
                    dst3d.at<unsigned char>(idx) = min_val;
                }
            }
        }
    }
    //ccShowSlice3Dmat(&dst3d, CV_8U);
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
            if(label_map->at<int>(i) > 0 && mask->at<unsigned char>(i) > 0){
                cc_sz[label_map->at<int>(i) - 1] ++;
            }
        }
    }
    size_t largest_id;
    vec_max(cc_sz, largest_id);

    return (int)largest_id + 1;
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
bool isempty(Mat src3d, int datatype, float threshold_in){
    return isempty(&src3d, datatype, threshold_in);
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
 * @param val3d: CV_8U, CV_32F, CV_32S
 * @param src3d: mask for value extraction, CV_8U
 * @param datatype
 * @param threshold_in
 * @return
 */
vector<float> extractValsGivenMask(Mat *val3d, int datatype, Mat *src3d, float threshold_in){
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
            if(src3d->at<unsigned char>(i) > threshold_in){
                fg_vals.push_back(val3d->at<float>(i));
            }
        }
    }else if (datatype == CV_32S){
        FOREACH_i_ptrMAT(src3d){
            if(src3d->at<unsigned char>(i) > threshold_in){
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
            fg_vals.push_back((float)vol3d->at<unsigned char>(idx[i]));
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            fg_vals.push_back(vol3d->at<float>(idx[i]));
        }
    }else if (datatype == CV_32S){
        FOREACH_i(idx){
            fg_vals.push_back((float)vol3d->at<int>(idx[i]));
        }
    }
    return fg_vals;
}

double extractSumGivenIdx(Mat *vol3d, vector<size_t> idx, int datatype){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    double sum = 0;
    if (datatype == CV_8U){
        FOREACH_i(idx){
            sum += (double)vol3d->at<unsigned char>(idx[i]);
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            sum += (double)vol3d->at<float>(idx[i]);
        }
    }else if (datatype == CV_32S){
        FOREACH_i(idx){
            sum += (double)vol3d->at<int>(idx[i]);
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
    int page_sz = im_sz[0] * im_sz[1];
    vector<int> n_y(8), n_x(8);
    n_y = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    n_x = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    for(int i=0; i < 8; i++){
        if(inField(y + n_y[i], x + n_x[i], im_sz)
                && fgMap->at<unsigned char>(vol_sub2ind(y + n_y[i], x + n_x[i], z, im_sz[1], page_sz)) == 0){
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
        y = remain / sz[1];
        x = remain - y * sz[1];
        for(int j = 0; j < n_y.size(); j++){
            if(inField(y + n_y[j], x + n_x[j], z + n_z[j], sz)){
                cur_nei_idx = x + n_x[j] + (y + n_y[j])*sz[1] + (z + n_z[j])*page_sz;
                if(!checked_label[cur_nei_idx]){
                    center_idx[linkage_cnt] = idx[i];
                    nei_idx[linkage_cnt] = cur_nei_idx;
                    linkage_cnt++;
                }
            }
        }
        checked_label[idx[i]] = true;
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

    vector<size_t> arc_tail_in_fg, arc_head;
    int mat_sz[] = {label_map->size[0], label_map->size[1], label_map->size[2]};
    vector<size_t> valid_fg_idx = fgMapIdx(fgMap, CV_8U, 0);
    //// undirected neighbor extraction, for each neighboring pair, they only appear once in
    /// the output.
    neighbor_idx(valid_fg_idx, arc_tail_in_fg, arc_head, mat_sz, connect);
    vector<double> arc_capacity (arc_tail_in_fg.size()); //BK algorithm accept only integer capacity
    float p1, p2;
    double cost;
    FOREACH_i(arc_tail_in_fg){
        p1 = scoreMap->at<float>(arc_tail_in_fg[i]);
        p2 = scoreMap->at<float>(arc_head[i]);
        if(p1+p2 == 0 || p1==0){
            arc_capacity[i] = numeric_limits<int>::max();
        }else{
            if (cost_design[0]==ARITHMETIC_AVERAGE){
                cost = pow((2/(p1+p2)), cost_design[1]);

            }else if (cost_design[0]==GEOMETRIC_AVERAGE){
                if (p2==0) p2 = p1; // there should be no non-negative score, so this line should never be used
                cost = pow((1/sqrt(p1*p2)), cost_design[1]);
            }
            if(cost < 0){
                qFatal("Wrong cap value!");
            }{
                if (cost > INFINITY) arc_capacity[i] = INFINITY; //float infinity is enough
                else arc_capacity[i] = cost;//
            }
        }
    }
    vector<vector<int>> label_voxIdx;
    extractVoxIdxList(label_map, label_voxIdx, numCC, false /*no_bk*/);
    vector<size_t> label_voxIdx_sz(label_voxIdx.size());
    FOREACH_i(label_voxIdx) label_voxIdx_sz[i] = label_voxIdx[i].size();
    size_t overall_seed_vox_num = accumulate(label_voxIdx_sz.begin(), label_voxIdx_sz.end(), 0);

    typedef Graph<double,double,double> GraphType; // the three types are forced to be the same: kind of wired
    set<int> bgsinkIds;
    if(bg2sink){
        FOREACH_j(arc_head){
            if (fgMap->at<unsigned char>(arc_head[j])==0 &&
                    label_map->at<int>(arc_head[j]) == 0){
                bgsinkIds.insert(arc_head[j]);
            }
        }
    }
    GraphType *g = new GraphType(label_map->total(),
                      (arc_tail_in_fg.size() + overall_seed_vox_num + bgsinkIds.size())*2);
    for(int region_id = 1; region_id <= numCC; region_id++){
        //// build sink_ids: valid sink nodes will not be larger than node size
        /// but may duplicate in the vector of sink_ids, which should be OK
//        vector<size_t> sink_ids (2*arc_head.size());
        //fill(sink_ids.begin(), sink_ids.end(), -1);
//        size_t sink_ids_cnt = 0;

//        if(bg2sink){//(true) {//
//            FOREACH_j(arc_head){
//                if (fgMap->at<unsigned char>(arc_head[j])==0 ||
//                        (label_map->at<int>(arc_head[j])>0 &&
//                        label_map->at<int>(arc_head[j]) != region_id)){
//                    sink_ids[sink_ids_cnt] = arc_head[j];
//                    sink_ids_cnt++;
//                }
//                if (fgMap->at<unsigned char>(arc_tail_in_fg[j])==0 ||
//                        (label_map->at<int>(arc_tail_in_fg[j])>0 &&
//                        label_map->at<int>(arc_tail_in_fg[j]) != region_id)){
//                    sink_ids[sink_ids_cnt] = arc_tail_in_fg[j];
//                    sink_ids_cnt++;
//                }
//            }
//        }
//        else{
//            FOREACH_j(arc_head){
//                if (fgMap->at<unsigned char>(arc_head[j])>0 &&
//                        label_map->at<int>(arc_head[j])>0 &&
//                        label_map->at<int>(arc_head[j]) != region_id){
//                    sink_ids[sink_ids_cnt] = arc_head[j];
//                    sink_ids_cnt++;
//                }
//                if (fgMap->at<unsigned char>(arc_tail_in_fg[j])>0 &&
//                        label_map->at<int>(arc_tail_in_fg[j])>0 &&
//                        label_map->at<int>(arc_tail_in_fg[j]) != region_id){
//                    sink_ids[sink_ids_cnt] = arc_tail_in_fg[j];
//                    sink_ids_cnt++;
//                }
//            }
//        }
//        sink_ids.resize(sink_ids_cnt);
//        if(sink_ids_cnt == 0 || label_voxIdx[i-1].size() == 0){
//            if(isempty(*label_map == 1, CV_8U)){
//                qInfo("Fatal error: no label 1.");
//            }
//            if(isempty(*label_map == 2, CV_8U)){
//                qInfo("Fatal error: no label 2.");
//            }
//            if(isempty(*label_map == 3, CV_8U)){
//                qInfo("Fatal error: no label 3.");
//            }
//            ccShowSlice3Dmat(*label_map == 1, CV_8U);
//            ccShowSlice3Dmat(*label_map == 2, CV_8U);
//            ccShowSlice3Dmat(*label_map == 3, CV_8U);
//            ccShowSliceLabelMat(label_map);
//            qInfo("Fatal error: no sink node.");
//        }
//        assert (sink_ids_cnt > 0 && label_voxIdx[region_id-1].size() > 0);

        /////////////////////////////////////////////
        //     max-flow the get the grown region   //
        /////////////////////////////////////////////
        g -> add_node(label_map->total());
        FOREACH_i(label_voxIdx){
            if((int)i == (region_id - 1)){
                FOREACH_j(label_voxIdx[i]){
                    g -> add_tweights( label_voxIdx[i][j],  INFINITY, 0);
                }
            }else{
                FOREACH_j(label_voxIdx[i]){
                    g -> add_tweights( label_voxIdx[i][j], 0,  INFINITY);
                }
            }
        }
        if(bgsinkIds.size() > 0){
            for(auto it : bgsinkIds){
                g -> add_tweights( it, 0,  INFINITY);
            }
        }
        FOREACH_j(arc_tail_in_fg){
            g -> add_edge( arc_tail_in_fg[j], arc_head[j],   arc_capacity[j],   arc_capacity[j] );
        }
//        if (bg2sink){
//            FOREACH_j(label_voxIdx[region_id-1]){
//                g -> add_tweights( label_voxIdx[region_id-1][j],  numeric_limits<int>::max(), 0);
//            }
//            FOREACH_j(sink_ids){
//                g -> add_tweights( sink_ids[j],  0, numeric_limits<int>::max());
//            }
//            FOREACH_j(arc_tail_in_fg){
//                g -> add_edge( arc_tail_in_fg[j], arc_head[j],   arc_capacity[j],   arc_capacity[j] );
//            }
//        }else{
//            FOREACH_j(label_voxIdx[region_id-1]){
//                g -> add_tweights( label_voxIdx[region_id-1][j],  numeric_limits<int>::max(), 0);
//            }
//            g -> add_edge( label_voxIdx[region_id-1][0], sink_ids[0],   10000,   10000);
//            FOREACH_j(sink_ids){
//                g -> add_tweights( sink_ids[j],  0, numeric_limits<int>::max());
//                int f = g -> maxflow();
//                qInfo("Node problem? %ld, %d", j, f);
//            }
////            FOREACH_j(arc_tail_in_fg){
////                g -> add_edge( arc_tail_in_fg[j], arc_head[j],   arc_capacity[j],   arc_capacity[j] );
////            }
//            //g -> add_tweights( label_voxIdx[region_id-1][0],  numeric_limits<int>::max(), 0);
//            //g -> add_tweights( sink_ids[0],  0, numeric_limits<int>::max());

////            FOREACH_j(arc_tail_in_fg){
////                if(fgMap->at<unsigned char>(arc_head[j])>0){
////                    //g -> add_edge( arc_tail_in_fg[j], arc_head[j],   arc_capacity[j],   arc_capacity[j]);
////                    g -> add_edge( arc_tail_in_fg[j], arc_head[j],   0,   0);
////                }
////            }
//            //g -> add_edge( label_voxIdx[region_id-1][0], sink_ids[0],   10000,   10000);
//        }
//        qInfo("sink:%ld, src:%ld, arc:%ld", sink_ids_cnt,
//              label_voxIdx[region_id-1].size(), arc_tail_in_fg.size());
        double flow = g -> maxflow();
        if(flow < 0){
            qFatal("The flow number is wrong!");
        }
        size_t region_sz = 0;
        FOREACH_j(valid_fg_idx){
            if (g->what_segment(valid_fg_idx[j]) == GraphType::SOURCE){
//                if (outLabelMap.at<int>(valid_fg_idx[j]) > 0){
//                    qInfo("error occured");
//                }
                outLabelMap.at<int>(valid_fg_idx[j]) = region_id;
                region_sz ++;
            }
        }
        qInfo("%ld out of %ld voxel in region %d.", region_sz, valid_fg_idx.size(),
              region_id);
        g->reset();
    }
    delete g;
}

void setValMat(Mat &vol3d, int datatype, vector<size_t> idx, float v){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);

    if (datatype == CV_8U){
        unsigned char v0 = (unsigned char)v;
        FOREACH_i(idx){
            vol3d.at<unsigned char>(idx[i]) = v0;
        }
    }else if (datatype == CV_32F){
        FOREACH_i(idx){
            vol3d.at<float>(idx[i]) = v;
        }
    }else if (datatype == CV_32S){
        int v0 = (int)v;
        FOREACH_i(idx){
            vol3d.at<int>(idx[i]) = v0;
        }
    }
}

void setValMat(Mat &src, int datatype, Mat *mask, float v){
    vector<size_t> idx = fgMapIdx(mask, CV_8U, 0);
    setValMat(src, datatype, idx, v);
}
void findGapGatePoints(Mat* label_map, int target_label0, vector<size_t> &gap_idx_y_x){
    float dy, dx;
    vector<size_t> target_voxIdx0;
    vector<int> y, x, z;
    extractVoxIdxGivenId(label_map, target_voxIdx0, target_label0);
    vec_ind2sub(target_voxIdx0, y, x, z, label_map->size);
    float y_center = vec_mean(y);
    float x_center = vec_mean(x);
    float min_dist2center = INFINITY;
    vector<float> dist2center(gap_idx_y_x.size()/3);
    for(size_t i = 0; i <gap_idx_y_x.size(); i+=3){
        dy = gap_idx_y_x[i+1]-y_center;
        dx = gap_idx_y_x[i+2]-x_center;
        dist2center[i/3] = sqrt(dy*dy + dx*dx);
        if(dist2center[i/3] < min_dist2center){
            min_dist2center = dist2center[i/3];
        }
    }
    //update y and x
    vector<float> scaled_y(gap_idx_y_x.size()/3), scaled_x(gap_idx_y_x.size()/3);
    for(size_t i = 0; i < gap_idx_y_x.size(); i+=3){
        dy = gap_idx_y_x[i+1]-y_center;
        dx = gap_idx_y_x[i+2]-x_center;
        scaled_y[i/3] = y_center + dy * (min_dist2center / dist2center[i/3]);
        scaled_x[i/3] = x_center + dx * (min_dist2center / dist2center[i/3]);
    }
    float max_dist = 0, curr_dist;
    size_t target_i, target_j;
    for(int i = 0 ; i < scaled_y.size(); i ++){
        for (int j = i + 1; j < scaled_x.size(); j++){
            curr_dist = pow((scaled_y[i] - scaled_y[j]),2) +
                    pow((scaled_x[i] - scaled_x[j]),2);
            if(curr_dist > max_dist){
                max_dist = curr_dist;
                target_i = i * 3;
                target_j = j * 3;
            }
        }
    }

    if(max_dist == 0){
        max_dist = 0;
        for(int i = 0 ; i < gap_idx_y_x.size(); i += 3){
            for (int j = i + 3; j < gap_idx_y_x.size(); j+= 3){
                curr_dist = pow((gap_idx_y_x[i+1] - gap_idx_y_x[j+1]),2) +
                        pow((gap_idx_y_x[i+2] - gap_idx_y_x[j+2]),2);
                if(curr_dist > max_dist){
                    max_dist = curr_dist;
                    target_i = i * 3;
                    target_j = j * 3;
                }
            }
        }
    }
//    qInfo("%ld element and %ld element is the best pair out of %ld ",
//          target_i, target_j, gap_idx_y_x.size());
    gap_idx_y_x[0] = gap_idx_y_x[target_i];
    gap_idx_y_x[1] = gap_idx_y_x[target_i+1];
    gap_idx_y_x[2] = gap_idx_y_x[target_i+2];
    gap_idx_y_x[3] = gap_idx_y_x[target_j];
    gap_idx_y_x[4] = gap_idx_y_x[target_j+1];
    gap_idx_y_x[5] = gap_idx_y_x[target_j+2];
    gap_idx_y_x.resize(6);
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
    vector<vector<size_t>> sub_idx0(label_map->size[2]);
    vector<vector<size_t>> sub_idx1(label_map->size[2]);
    int im_sz[] = {label_map->size[0], label_map->size[1]};
    int x, y, z, remain;
    int page_sz = im_sz[0] * im_sz[1];

    int n_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    int n_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    vector<vector<size_t>> gap_idx_z_divid(label_map->size[2]);
    FOREACH_i(gap_idx){
        z = gap_idx[i] / page_sz;
        remain = gap_idx[i] - (z*page_sz);
        y = remain / im_sz[1];
        x = remain - y * im_sz[1];
        gap_idx_z_divid[z].push_back(gap_idx[i]);
        if (label_map->at<int>(gap_idx[i]) == target_label0){
            for(int j = 0; j < 8; j++){
                if(inField(y+n_y[j], x+n_x[j], im_sz) &&
                        label_map->at<int>( gap_idx[i] + n_y[j] * im_sz[1] + n_x[j])!=target_label0){
                    sub_idx0[z].push_back(gap_idx[i]);
                    sub_idx0[z].push_back(y);
                    sub_idx0[z].push_back(x);
                    break;
                }
            }
        }else if(label_map->at<int>(gap_idx[i]) == target_label1){
            for(int j = 0; j < 8; j++){
                if(inField(y+n_y[j], x+n_x[j], im_sz) &&
                        label_map->at<int>(gap_idx[i] + n_y[j] * im_sz[1] + n_x[j])!=target_label1){
                    sub_idx1[z].push_back(gap_idx[i]);
                    sub_idx1[z].push_back(y);
                    sub_idx1[z].push_back(x);
                    break;
                }
            }
        }
    }
    vector<size_t> refined_gap_idx(gap_idx.size());
    size_t refined_gap_idx_size = 0;
    for (int z_s = 0; z_s < label_map->size[2]; z_s++){
        if(sub_idx0[z_s].size()<=6 || sub_idx1[z_s].size()<=6){
            for(size_t i = 0; i < sub_idx0[z_s].size(); i+=3){
                refined_gap_idx[refined_gap_idx_size++] = sub_idx0[z_s][i];
            }
            for(size_t i = 0; i < sub_idx1[z_s].size(); i+=3){
                refined_gap_idx[refined_gap_idx_size++] = sub_idx1[z_s][i];
            }
            continue;
        }
        //// find the two nodes far away (with respect to the region)
        findGapGatePoints(label_map, target_label0, sub_idx0[z_s]);
        findGapGatePoints(label_map, target_label1, sub_idx1[z_s]);

        // trim points that outside the quadrilateral
        // 1. get the two lines
        float k_b0[2], k_b1[2];
        if(sub_idx1[z_s][5] == sub_idx0[z_s][5]) k_b0[0] = (float)(sub_idx1[z_s][4] - sub_idx0[z_s][4]);
        else k_b0[0] = ((float)sub_idx1[z_s][4] - sub_idx0[z_s][4]) / (sub_idx1[z_s][5] - sub_idx0[z_s][5]);
        k_b0[1] = sub_idx1[z_s][4] -  k_b0[0] * sub_idx1[z_s][5];

        if(sub_idx1[z_s][2] == sub_idx0[z_s][2]) k_b1[0] = (float)(sub_idx1[z_s][1] - sub_idx0[z_s][1]);
        else k_b1[0] = ((float)sub_idx1[z_s][1] - sub_idx0[z_s][1]) / (sub_idx1[z_s][2] - sub_idx0[z_s][2]);
        k_b1[1] = sub_idx1[z_s][1] -  k_b1[0] * sub_idx1[z_s][2];

        float intersect_x;
        if(k_b0[0] == k_b1[0]) intersect_x = INFINITY;
        else intersect_x = ((k_b1[1] - k_b0[1]))/(k_b0[0] - k_b1[0]);
        //float intersect_y = k_b1[1] * intersect_x + k_b1[0];
        if(intersect_x <= MAX(sub_idx0[z_s][5], sub_idx1[z_s][5]) && intersect_x >= MIN(sub_idx0[z_s][5], sub_idx1[z_s][5])){
            //intersect_y <= MAX(sub_idx0[z_s][4], sub_idx1[z_s][4]) && intersect_x >= MIN(sub_idx0[z_s][4], sub_idx1[z_s][4]) ){
            size_t tmp = sub_idx0[z_s][0];
            sub_idx0[z_s][0] = sub_idx0[z_s][3];
            sub_idx0[z_s][3] = tmp;
            tmp = sub_idx0[z_s][1];
            sub_idx0[z_s][1] = sub_idx0[z_s][4];
            sub_idx0[z_s][4] = tmp;
            tmp = sub_idx0[z_s][2];
            sub_idx0[z_s][2] = sub_idx0[z_s][5];
            sub_idx0[z_s][5] = tmp;

            //        if(sub_idx1[5] == sub_idx0[5]) k_b0[0] = (sub_idx1[4] - sub_idx0[4]);
            //        else k_b0[0] = (sub_idx1[4] - sub_idx0[4]) / (sub_idx1[5] - sub_idx0[5]);
            //        k_b0[1] = sub_idx1[4] -  k_b0[0] * sub_idx1[5];

            //        if(sub_idx1[2] == sub_idx0[2]) k_b1[0] = (sub_idx1[1] - sub_idx0[1]);
            //        else k_b1[0] = (sub_idx1[1] - sub_idx0[1]) / (sub_idx1[2] - sub_idx0[2]);
            //        k_b1[1] = sub_idx1[1] -  k_b1[0] * sub_idx1[2];
        }
        // remove pionts outside the quadrilateral
        vector<Point> contour(4);
        contour[0].y = sub_idx0[z_s][1];
        contour[0].x = sub_idx0[z_s][2];
        contour[1].y = sub_idx0[z_s][4];
        contour[1].x = sub_idx0[z_s][5];
        contour[2].y = sub_idx1[z_s][4];
        contour[2].x = sub_idx1[z_s][5];
        contour[3].y = sub_idx1[z_s][1];
        contour[3].x = sub_idx1[z_s][2];
        //size_t valid_idx_cnt = 0;

        FOREACH_i(gap_idx_z_divid[z_s]){
            if (label_map->at<int>(gap_idx_z_divid[z_s][i]) == target_label0 ||
                    label_map->at<int>(gap_idx_z_divid[z_s][i]) == target_label1){
                refined_gap_idx[refined_gap_idx_size++] = gap_idx_z_divid[z_s][i];
                //gap_idx[valid_idx_cnt] = gap_idx[i];
                //valid_idx_cnt ++;
            }else{
                z = gap_idx_z_divid[z_s][i] / page_sz; // z == z_s
                assert(z == z_s);
                remain = gap_idx_z_divid[z_s][i] - (z*page_sz);
                y = remain / im_sz[1];
                x = remain - y * im_sz[1];
                if(pointPolygonTest(contour, Point2f((float)y, (float)x), false)>=0){
                    //gap_idx[valid_idx_cnt] = gap_idx[i];
                    //valid_idx_cnt ++;
                    refined_gap_idx[refined_gap_idx_size++] = gap_idx_z_divid[z_s][i];
                }
            }
//            if(refined_gap_idx_size == 108){
//                qInfo("Degug!");
//            }
        }
        //gap_idx.resize(valid_idx_cnt);
    }
    refined_gap_idx.resize(refined_gap_idx_size);
    if(refined_gap_idx_size < 1 || refined_gap_idx[refined_gap_idx_size - 1] > label_map->total()){
        qInfo("no valid gap is found");
        gap_idx.resize(0);
    }else{
        gap_idx.resize(0);
        gap_idx.insert(gap_idx.begin(), refined_gap_idx.begin(), refined_gap_idx.end());
    }
}
/**
 * @brief extractGapVoxel: 2d gap test
 * @param label_map
 * @param numCC
 * @param gap_radius
 * @param gap_voxIdx
 */
void extractGapVoxel(Mat *label_map, Mat *fgMap, int numCC, int gap_radius,
                     vector<vector<size_t>> &gap_voxIdx, vector<int> tested_flag){
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
        y = remain / im_sz[1];
        x = remain - y * im_sz[1];
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
            int tmp_idx = (MIN(nei_seed_ids[0]-1, nei_seed_ids[1]-1)) * numCC + (MAX(nei_seed_ids[0]-1, nei_seed_ids[1]-1));
            //if(tested_flag[tmp_idx] == false){
            gap_voxIdx[tmp_idx].push_back(valid_fg_idx[i]);
            //}
        }
    }
    FOREACH_i(gap_voxIdx){
        if(gap_voxIdx[i].size() > 10){
            if (tested_flag[i] > 0){ // this gap is alreay true, does not need to test
                gap_voxIdx[i].resize(0);
            }else{
//                if(gap_voxIdx[i].size() <= 4){
//                    qInfo("too small size of gap");
//                }
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
            z = curr_idx[j] / page_sz;
            remain = curr_idx[j] - (z*page_sz);
            y = remain / im_sz[1];
            x = remain - y * im_sz[1];
            for(int k = 0; k < 8; k++){
                if(inField(y+n_y[k], x+n_x[k], im_sz)){
                    vol_sub2ind(tmp_idx, y+n_y[k], x+n_x[k], z, curr_fg.size);
                    if(curr_fg.at<unsigned char>(tmp_idx) > 0){
                        neighbor_idx_list[i].push_back(tmp_idx);
                        curr_fg.at<unsigned char>(tmp_idx) = 0;
                    }
                }
            }
        }
    }
}
/**
 * @brief subVolExtract: given range, extrat the sub volume from src voluem
 * @param src
 * @param datatype
 * @param subVol
 * @param xyz_range
 */
void subVolExtract(Mat *src, int datatype, Mat &subVol, Range yxz_range[3]){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    int size[3] = {yxz_range[0].end - yxz_range[0].start,
                  yxz_range[1].end - yxz_range[1].start,
                  yxz_range[2].end - yxz_range[2].start};
    subVol = Mat(3, size, datatype);
    size_t src_pgsz = src->size[0] * src->size[1];
    size_t sub_pgsz = size[0] * size[1];
    size_t src_page_shift, sub_page_shift, src_row_shift, sub_row_shift;
    if (datatype == CV_8U){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src->size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i) =
                            src->at<unsigned char>(src_page_shift + src_row_shift + i + yxz_range[1].start);
                }
            }
        }
    }else if(datatype == CV_32F){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src->size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    subVol.at<float>(sub_page_shift + sub_row_shift + i) =
                            src->at<float>(src_page_shift + src_row_shift + i + yxz_range[1].start);
                }
            }
        }
    }else if(datatype == CV_32S){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src->size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    subVol.at<int>(sub_page_shift + sub_row_shift + i) =
                            src->at<int>(src_page_shift + src_row_shift + i + yxz_range[1].start);
                }
            }
        }
    }
}
/**
 * @brief subVolReplace: replace the data in src by the subVol
 * @param src
 * @param datatype
 * @param subVol
 * @param xyz_range
 */
void subVolReplace(Mat &src, int datatype, Mat &subVol, Range yxz_range[3], int start){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    int size[3] = {yxz_range[0].end - yxz_range[0].start,
                  yxz_range[1].end - yxz_range[1].start,
                  yxz_range[2].end - yxz_range[2].start};

    size_t src_pgsz = src.size[0] * src.size[1];
    size_t sub_pgsz = size[0] * size[1];
    size_t src_page_shift, sub_page_shift, src_row_shift, sub_row_shift;
    if (datatype == CV_8U){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if(subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i) > 0){
                        src.at<unsigned char>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                         = subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i) + (unsigned char)start;
                    }
                }
            }
        }
    }else if(datatype == CV_32F){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if(subVol.at<float>(sub_page_shift + sub_row_shift + i) > 0){
                        src.at<float>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                                 = subVol.at<float>(sub_page_shift + sub_row_shift + i) + (float)start;
                    }
                }
            }
        }
    }else if(datatype == CV_32S){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if(subVol.at<int>(sub_page_shift + sub_row_shift + i) > 0){
                        src.at<int>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                                = subVol.at<int>(sub_page_shift + sub_row_shift + i) + start;
                    }
                }
            }
        }
    }
}

/**
 * @brief subVolReplace: replace the data in src by a unique val with mask defined in subVol
 * @param src
 * @param datatype
 * @param subVol
 * @param xyz_range
 */
void subVolReplace(Mat &src, int datatype, Mat &subVol, float val, Range yxz_range[3]){
    assert(datatype == CV_8U || datatype == CV_32F || datatype == CV_32S);
    assert(subVol.type() == CV_8U);
    int size[3] = {yxz_range[0].end - yxz_range[0].start,
                  yxz_range[1].end - yxz_range[1].start,
                  yxz_range[2].end - yxz_range[2].start};

    size_t src_pgsz = src.size[0] * src.size[1];
    size_t sub_pgsz = size[0] * size[1];
    size_t src_page_shift, sub_page_shift, src_row_shift, sub_row_shift;
    if (datatype == CV_8U){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if (subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i) > 0){
                        src.at<unsigned char>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                         = (unsigned char) val;
                    }
                }
            }
        }
    }else if(datatype == CV_32F){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if(subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i) > 0){
                        src.at<float>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                                 = val;
                    }
                }
            }
        }
    }else if(datatype == CV_32S){
        for(int k = 0; k < size[2]; k++){
            src_page_shift = (k + yxz_range[2].start) * src_pgsz;
            sub_page_shift = k * sub_pgsz;
            for(int j = 0; j < size[0]; j++){
                src_row_shift = (j + yxz_range[0].start) * src.size[1];
                sub_row_shift = j * size[1];
                for(int i = 0; i < size[1]; i++){
                    if(subVol.at<unsigned char>(sub_page_shift + sub_row_shift + i)>0){
                        src.at<int>(src_page_shift + src_row_shift + i + yxz_range[1].start)
                                = (int)val;
                    }
                }
            }
        }
    }
}


/** distanceTransRegion2Region: distance transform between two regions
 *
 */
float distanceTransRegion2Region(bool *bw_ref_cell, vector<int> ref_range_xyz,
                                                       bool *bw_mov_cell, vector<int> mov_range_xyz,
                                                       vector<double> shift_xyz, vector<float> dist){
    dist.resize(2);//distance from r1 to r2 and counter direction
    // use the slightly edited mex-version previously designed for matlab
    // for matlab, the matrix is saved in column-by-column
    size_t *ref_sz_yxz = new size_t[3];
    ref_sz_yxz[0] = ref_range_xyz[1];
    ref_sz_yxz[1] = ref_range_xyz[0];
    ref_sz_yxz[2] = ref_range_xyz[2];

    size_t *mov_sz_yxz = new size_t[3];
    mov_sz_yxz[0] = mov_range_xyz[1];
    mov_sz_yxz[1] = mov_range_xyz[0];
    mov_sz_yxz[2] = mov_range_xyz[2];

    float *shift_yxz = new float[3];
    shift_yxz[0] = shift_xyz[1];
    shift_yxz[1] = shift_xyz[0];
    shift_yxz[2] = shift_xyz[2];

    size_t ref_l = ref_range_xyz[0]*ref_range_xyz[1]*ref_range_xyz[2];
    size_t mov_l = mov_range_xyz[0]*mov_range_xyz[1]*mov_range_xyz[2];
    // ref cell to mov cell
    float *dist_voxwise = new float[mov_l];
    dt3d(bw_ref_cell, ref_sz_yxz, mov_sz_yxz, shift_yxz, dist_voxwise);
    double dist_sum = 0;
    size_t valid_n = 0;
    for(size_t i = 0; i < mov_l; i++){
        if(bw_mov_cell[i]){
            dist_sum += dist_voxwise[i];
            valid_n ++;
        }
    }
    dist[0] = dist_sum / valid_n;
    // mov cell to ref cell
    float *dist2_voxwise = new float[ref_l];
    shift_yxz[0] = -shift_xyz[1];
    shift_yxz[1] = -shift_xyz[0];
    shift_yxz[2] = -shift_xyz[2];
    dt3d(bw_mov_cell, mov_sz_yxz, ref_sz_yxz, shift_yxz, dist2_voxwise);
    dist_sum = 0;
    valid_n = 0;
    for(size_t i = 0; i < ref_l; i++){
        if(bw_ref_cell[i]){
            dist_sum += dist_voxwise[i];
            valid_n ++;
        }
    }
    dist[1] = dist_sum / valid_n;

    return MAX(dist[0], dist[1]);
}
/**************************Functions for debug*******************************/
void ccShowSlice3Dmat(Mat src3d, int datatype, int slice, bool binary){
    ccShowSlice3Dmat(&src3d, datatype, slice, binary);
}
void ccShowSlice3Dmat(Mat *src3d, int datatype, int slice, bool binary){
    assert(datatype == CV_8U || datatype == CV_16U || datatype == CV_32F || datatype == CV_32S);
    int sz_single_frame = src3d->size[1] * src3d->size[0];
    Mat *single_slice;
    Mat slice4display;
    double min_v, max_v;
    while (true){
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
            //waitKeyEx(0);
            int key = waitKeyEx(0);
            if(src3d->dims <3) break;
            if(key == 65361) slice = MAX(slice-1, 0); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
            else if(key == 65363) slice = MIN(slice+1, src3d->size[2] - 1); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
            else {
                break;
            }
            //destroyWindow(title);
        }else{
            if(datatype == CV_8U){
                unsigned char *ind = (unsigned char*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
                single_slice = new Mat(2, src3d->size, CV_8U, ind);
                single_slice->copyTo(slice4display);
                minMaxIdx(slice4display, &min_v, &max_v);
            }else if(datatype == CV_16U){
                unsigned short *ind = (unsigned short*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
                single_slice = new Mat(2, src3d->size, CV_16U, ind);
                minMaxIdx(*single_slice, &min_v, &max_v);
                single_slice->copyTo(slice4display);
                slice4display *= 10;
                //normalize(*single_slice, *single_slice, 1, 0, NORM_MINMAX, CV_32F);
            }else if(datatype == CV_32S){
                int *ind = (int*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
                single_slice = new Mat(2, src3d->size, CV_32S, ind);
                minMaxIdx(*single_slice, &min_v, &max_v);
                single_slice->copyTo(slice4display);
                //normalize(slice4display, slice4display, 1, 0, NORM_MINMAX, CV_32F);
            }else if(datatype == CV_32F){
                float *ind = (float*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
                single_slice = new Mat(2, src3d->size, CV_32F, ind);
                minMaxIdx(*single_slice, &min_v, &max_v);
                single_slice->copyTo(slice4display);
                normalize(slice4display, slice4display, 1, 0, NORM_MINMAX);
            }else{
                single_slice = nullptr;
                qFatal("Unsupported data type!");
            }

            stringstream s_min, s_max;
            s_min << fixed << setprecision(2) << min_v;
            s_max << fixed << setprecision(2) << max_v;
            string title = "Slice:"+to_string(slice) + ", min:" + s_min.str() + ", max:" + s_max.str();
            imshow(title, slice4display);
            int key = waitKeyEx(0);
            if(src3d->dims <3) break;
            if(key == 65361) {
                slice = MAX(slice-1, 0); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
                destroyWindow(title);
            }
            else if(key == 65363) {
                slice = MIN(slice+1, src3d->size[2] - 1); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
                destroyWindow(title);
            }
            else {
                break;
                destroyWindow(title);
            }
            //destroyWindow(title);
            //destroyAllWindows();
        }
    }
    //delete single_slice;
}

/**
 * @brief label2rgb3d
 * @param src: CV_32S
 * @param dst: CV_8UC3, RGB data
 */
void label2rgb3d(Mat &src_label, Mat &src_intensity, Mat4b &dst)
{
    Mat3b colormap;
    colorMapGen(&src_label, colormap);
    // Apply color mapping
    dst = Mat4b(src_label.dims, src_label.size, Vec4b(0,0,0));
    Vec3b cur_cl;
    //float a;
    FOREACH_i_MAT(src_label){
        cur_cl = colormap(src_label.at<int>(i));
//        if (src_label.at<int>(i) > 0){
//            dst.at<Vec4b>(i)(0) = 255;
//            dst.at<Vec4b>(i)(1) = 0;
//            dst.at<Vec4b>(i)(2) = 0;
//            dst.at<Vec4b>(i)(3) = 255;
//        }
        if (src_label.at<int>(i) > 0){
            dst.at<Vec4b>(i)(0) = cur_cl(0);
            dst.at<Vec4b>(i)(1) = cur_cl(1);
            dst.at<Vec4b>(i)(2) = cur_cl(2);
            dst.at<Vec4b>(i)(3) = 255;
        }else{
            dst.at<Vec4b>(i)(1) = src_intensity.at<unsigned char>(i); //g
            dst.at<Vec4b>(i)(3) = src_intensity.at<unsigned char>(i)/2;//the lower the more transparent it is
        }
    }
}

//hsv rgb2hsv(rgb in)
//{
//    hsv         out;
//    double      min, max, delta;

//    min = in.r < in.g ? in.r : in.g;
//    min = min  < in.b ? min  : in.b;

//    max = in.r > in.g ? in.r : in.g;
//    max = max  > in.b ? max  : in.b;

//    out.v = max;                                // v
//    delta = max - min;
//    if (delta < 0.00001)
//    {
//        out.s = 0;
//        out.h = 0; // undefined, maybe nan?
//        return out;
//    }
//    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
//        out.s = (delta / max);                  // s
//    } else {
//        // if max is 0, then r = g = b = 0
//        // s = 0, h is undefined
//        out.s = 0.0;
//        out.h = NAN;                            // its now undefined
//        return out;
//    }
//    if( in.r >= max )                           // > is bogus, just keeps compilor happy
//        out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
//    else
//    if( in.g >= max )
//        out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
//    else
//        out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

//    out.h *= 60.0;                              // degrees

//    if( out.h < 0.0 )
//        out.h += 360.0;

//    return out;
//}

void HsvToRgb(double h, double s, double v, double &r, double &g, double &b)
{
    double      hh, p, q, t, ff;
    long        i;

    if(s <= 0.0) {       // < is bogus, just shuts up warnings
        r = v;
        g = v;
        b = v;
        return;
    }
    hh = h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * ff));
    t = v * (1.0 - (s * (1.0 - ff)));

    switch(i) {
    case 0:
        r = v;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = v;
        b = p;
        break;
    case 2:
        r = p;
        g = v;
        b = t;
        break;

    case 3:
        r = p;
        g = q;
        b = v;
        break;
    case 4:
        r = t;
        g = p;
        b = v;
        break;
    case 5:
    default:
        r = v;
        g = p;
        b = q;
        break;
    }
}
/**
 * @brief colorMapGen
 * @param src
 * @param colorMap
 */
void colorMapGen(Mat *src, Mat3b &colormap, String colorType){
    double m;
    minMaxIdx(*src, nullptr, &m); // for 3 or more dims; for 2d, use minMaxLoc
    m++;

    if (colorType.compare("JET") == 0){
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
        srand(time(0));
        random_shuffle(g.begin() + 1, g.end());
        random_shuffle(r.begin() + 1, r.end());
        random_shuffle(b.begin() + 1, b.end());

        g.erase(remove_if(g.begin(), g.end(), [m](double v){ return v > m;}), g.end());
        r.erase(remove_if(r.begin(), r.end(), [m](double v){ return v > m; }), r.end());
        b.erase(remove_if(b.begin(), b.end(), [](double v){ return v < 1.0; }), b.end());

        Mat1d cmap(m, 3, double(0.0));
        int c, l;
        for (int i = 0; i < r.size(); ++i)
        {
            c = int(r[i])-1;
            cmap(c, 2) = u(i);
        }
        for (int i = 0; i < g.size(); ++i)
        {
            cmap(int(g[i])-1, 1) = u(i);
        }
        for (int i = 0; i < b.size(); ++i)
        {
            cmap(int(b[i])-1, 0) = u(u.rows - b.size() + i);
        }

        Mat3d cmap3 = cmap.reshape(3);

        //Mat3b colormap;
        cmap3.convertTo(colormap, CV_8U, 255.0);
    }else{ //"HSV"
        vector<double> g(m, 1);
        vector<double> r(m, 1);
        vector<double> b(m, 1);
        double h, s, v, rand_val;
        //float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); // 0-1
        //float r3 = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO))); // LO-HI
        for (int i = 0; i < m; i++)
        {
            h = i * 360.0/m;
            rand_val = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            s = 90.0 + rand_val * 10;
            rand_val = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            v = 80.0 + rand_val * 10;
            HsvToRgb(h,s/100,v/100, r[i], g[i], b[i]);
        }
        Mat1d cmap(m, 3, double(0.0));
        for (int i = 0; i < m; ++i)
        {
            cmap(i, 2) = b[i];
            cmap(i, 1) = g[i];
            cmap(i, 0) = r[i];
        }
        Mat3d cmap3 = cmap.reshape(3);
        //Mat3b colormap;
        cmap3.convertTo(colormap, CV_8U, 255.0);
    }
}
/**
 * @brief label2rgb3d
 * @param src: CV_32S
 * @param dst: CV_8UC3, RGB data
 */
void label2rgb2d(Mat1i &src, Mat3b &colormap, Mat3b &dst)
{
//    // Create JET colormap
//    double m;
//    minMaxLoc(src, nullptr, &m);
//    m++;

//    int n = ceil(m / 4);
//    Mat1d u(n*3-1, 1, double(1.0));

//    for (int i = 1; i <= n; ++i) {
//        u(i-1) = double(i) / n;
//        u((n*3-1) - i) = double(i) / n;
//    }

//    vector<double> g(n * 3 - 1, 1);
//    vector<double> r(n * 3 - 1, 1);
//    vector<double> b(n * 3 - 1, 1);
//    for (int i = 0; i < g.size(); ++i)
//    {
//        g[i] = ceil(double(n) / 2) - (int(m)%4 == 1 ? 1 : 0) + i + 1;
//        r[i] = g[i] + n;
//        b[i] = g[i] - n;
//    }

//    g.erase(remove_if(g.begin(), g.end(), [m](double v){ return v > m;}), g.end());
//    r.erase(remove_if(r.begin(), r.end(), [m](double v){ return v > m; }), r.end());
//    b.erase(remove_if(b.begin(), b.end(), [](double v){ return v < 1.0; }), b.end());

//    Mat1d cmap(m, 3, double(0.0));
//    for (int i = 0; i < r.size(); ++i) { cmap(int(r[i])-1, 2) = u(i); }
//    for (int i = 0; i < g.size(); ++i) { cmap(int(g[i])-1, 1) = u(i); }
//    for (int i = 0; i < b.size(); ++i) { cmap(int(b[i])-1, 0) = u(u.rows - b.size() + i); }

//    Mat3d cmap3 = cmap.reshape(3);

//    Mat3b colormap;
//    cmap3.convertTo(colormap, CV_8U, 255.0);


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
void ccShowSliceLabelMat(Mat src3d, int slice){
    ccShowSliceLabelMat(&src3d, slice);
};
void ccShowSliceLabelMat(Mat *src3d, int slice){
    Mat3b colormap;
    colorMapGen(src3d, colormap);
    while (true){
        int sz_single_frame = src3d->size[1] * src3d->size[0];
        int *ind = (int*)src3d->data + sz_single_frame*slice; // sub-matrix pointer
        //// The following line is very dangerous if the input
        /// matrix is not saved in continuous memory. E.g. if the input
        /// matrix is just a shallow crop of a intact mat, the index is
        /// not continuous.
        Mat single_slice(2, src3d->size, CV_32S, ind);
//        imshow("test_", single_slice > 0);
//        waitKey(0);
        Mat1i mat2show = single_slice;
        double min_v, max_v;

        minMaxIdx(mat2show, &min_v, &max_v);
        stringstream s_min, s_max;
        s_min << fixed << setprecision(2) << min_v;
        s_max << fixed << setprecision(2) << max_v;
        string title = "Slice:"+to_string(slice) + ", min:" + s_min.str() + ", max:" + s_max.str() + ".";
    //    if (getWindowProperty(title, WND_PROP_AUTOSIZE) >= 0){
    //        title += " -1";
    //    }
        int i = rand();
        title += to_string(i);
        Mat3b mat2display;
        label2rgb2d(mat2show, colormap, mat2display);
        string folder = "/home/ccw/Desktop/embryo_res_folder/crop_embryo_data_500x500x30x40/imgs/";
        //imwrite(folder + title + ".png", mat2show);
//        imshow(title, single_slice>0);
//        waitKey(0);
//        imshow(title, mat2show>0);
//        waitKey(0);
        imshow(title, mat2display);
        int key = waitKeyEx(0);
        if(src3d->dims <3) break;
        if(key == 65361) {
            slice = MAX(slice-1, 0); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
            destroyWindow(title);
        }
        else if(key == 65363) {
            slice = MIN(slice+1, src3d->size[2] - 1); //Left: 2424832 Up: 2490368 Right: 2555904 Down: 2621440
            destroyWindow(title);
        }
        else {
            break;
            destroyWindow(title);
        }
    }
}
