#ifndef IMG_BASIC_PROC_H
#define IMG_BASIC_PROC_H

/** the basic function used in volume image processing */

#include "img_basic_proc_declare.h"
using namespace cv;
/**
 * @brief principalCv2d
 * @param src3d: float mat
 * @param dst3d: float mat
 */
void principalCv2d(Mat* src3d, Mat &dst3d, float *sigma, int minIntensity = 0){
    src3d->copyTo(dst3d);
    smooth3Ddata(dst3d, sigma);

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
        Mat lx, ly, lxx, lyx, lxy, lyy;

        filter2D(subMatrix, lx, -1, kernelx);
        filter2D(subMatrix, ly, -1, kernely);
        filter2D(lx, lxy, -1, kernely);
        filter2D(lx, lxx, -1, kernelx);
        filter2D(ly, lyx, -1, kernelx); // the same as lxy
        filter2D(ly, lyy, -1, kernelx);
        for (int r = 0; r < x_size; r ++){
            for (int c = 0; c < y_size; c ++){
                if (src3d->at<float>(r,c,z) <= minIntensity){
                    dst3d.at<float>(r,c,z) = 0;
                    continue;
                }
                mat2x2.at<float>(0,0) = lxx.at<float>(r,c);
                mat2x2.at<float>(0,1) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,0) = lyx.at<float>(r,c);
                mat2x2.at<float>(1,1) = lyy.at<float>(r,c);
                eigen(mat2x2, eigen_values);
                dst3d.at<float>(r,c,z) = eigen_values.at<float>(0); // the largest eigen value
            }
        }
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernal_row.t(), kernal_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
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
    smooth3Ddata(dst3d, sigma);

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
        Mat lx, ly, lxx, lyx, lxy, lyy;

        filter2D(subMatrix, lx, -1, kernelx);
        filter2D(subMatrix, ly, -1, kernely);
        filter2D(lx, lxy, -1, kernely);
        filter2D(lx, lxx, -1, kernelx);
        filter2D(ly, lyx, -1, kernelx); // the same as lxy
        filter2D(ly, lyy, -1, kernelx);
        for (int r = 0; r < x_size; r ++){
            for (int c = 0; c < y_size; c ++){
                if (src3d->at<float>(r,c,z) <= minIntensity){
                    dst3d.at<float>(r,c,z) = 0;
                    continue;
                }
                mat2x2.at<float>(0,0) = lxx.at<float>(r,c);
                mat2x2.at<float>(0,1) = lxy.at<float>(r,c);
                mat2x2.at<float>(1,0) = lyx.at<float>(r,c);
                mat2x2.at<float>(1,1) = lyy.at<float>(r,c);
                eigen(mat2x2, eigen_values);
                dst3d.at<float>(r,c,z) = eigen_values.at<float>(0); // the largest eigen value
            }
        }
        //sepFilter2D(subMatrix, subMatrix, CV_32F, kernal_row.t(), kernal_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
}
void smooth3Ddata(Mat &data4smooth, const float *sigma)
{
    assert(data4smooth.dims == 3);
    int x_size  = data4smooth.size[0];
    int y_size  = data4smooth.size[1];
    int z_size  = data4smooth.size[2];
    size_t xy_size = x_size*y_size;

    int kernDimension = 2*ceil(2*sigma[0])+1;
    Mat kernal_row = getGaussianKernel(kernDimension, sigma[0], CV_32F);
    kernDimension = 2*ceil(2*sigma[1])+1;
    Mat kernal_col = getGaussianKernel(kernDimension, sigma[1], CV_32F);
    // Filter XY dimensions for every Z
    for (int z = 0; z < z_size; z++)
    {
        float *ind = (float*)data4smooth.data + z * xy_size; // sub-matrix pointer
        Mat subMatrix(2, data4smooth.size, CV_32F, ind);
        sepFilter2D(subMatrix, subMatrix, CV_32F, kernal_row.t(), kernal_col, Point(-1,-1), 0.0, BORDER_REPLICATE);
    }
    if (sigma[2] > 0){
        kernDimension = 2*ceil(2*sigma[2])+1;
        Mat kernal_z = getGaussianKernel(kernDimension, sigma[2], CV_32F);
        // Filter Z dimension
        float* kernGauss = (float *)kernal_z.data;
        unsigned kernSize = kernal_z.total();
        int kernMargin = (kernSize - 1)/2;
        float* lineBuffer = new float[z_size + 2*kernMargin];
        for (int y = 0; y < y_size; y++)
        {
            for (int x = 0; x < x_size; x++)
            {
                // Copy along Z dimension into a line buffer
                float* z_ptr = (float*)data4smooth.data + y * x_size + x;//same as data4smooth.ptr<float>(0, y, x)
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
                z_ptr = (float*)data4smooth.data + y * x_size + x;
                for (int z = 0; z < z_size; z++, z_ptr += xy_size)
                {
                    *z_ptr = 0.0f;
                    for (unsigned k = 0; k < kernSize; k++)
                    {
                        *z_ptr += lineBuffer[z+k]*kernGauss[k];
                    }
                }
            }
        }
        delete [] lineBuffer;
    }
}

#endif // IMG_BASIC_PROC_H


