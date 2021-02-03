#ifndef DT3D_HPP
#define DT3D_HPP
/* This function is a c++ version of edt_3dMex.cpp, which is from the matlab
 * version of the cell tracking framework.
 * Ideas are from Felzenszwalb and Huttenlocher's paper about distance
 * transform. Our distance transfrom follows the z->x->y function.
 * More info can be found in edt_3d_regFixed.m.
 *
 * Author: Congchao Wang
 * Yu's Lab, Virginia Tech
 * Date: Jan, 2021
 *
 */
#include <stdio.h>
//#include <mex.h>
#include <algorithm>
#include <limits>

const float INFINITY = std::numeric_limits<float>::max();

float square(float x) {
    return x*x;
}
/*
 * Set a set of values in f to nan (-1)
 * Parameters:
 * f: values
 * n: the number of nan elements
 * stride: step size
 */
void toNan(float *f,
        const size_t n,
        const size_t stride)
{
    for (int i = 0; i < n; i++) {
        f[i * stride] = -1;
    }
}
/* dt of 1d function using squared distance
 * Parameters:
 * p: location of valid parabolas: sorted
 * f: function of p
 * n: size of p or f
 * stride: step size
 * t: locations needed for distance transform: sorted
 * t_n: size of t
 */
void dt(int *p,
        float *f,
        size_t n,
        float *d,
        const size_t stride, /*for d*/
        float *t,
        const size_t t_n)
{

    int *v = new int[n];
    float *z = new float[n+1];
    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] = +INFINITY;
    float s;
    for (int q = 1; q < n; q++) {
        s  = (f[q] - f[v[k]])/(2*p[q]-2*p[v[k]])
            + (float)(p[q] + p[v[k]])/2;

        while (s <= z[k]) {
            k--;
            s  = (f[q] - f[v[k]])/(2*p[q]-2*p[v[k]])
                + (float)(p[q] + p[v[k]])/2;
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k+1] = +INFINITY;
    }
    k = 0;
    for (int q = 0; q < t_n; q++) {
        while (z[k+1] < t[q])
            k++;
        d[q * stride] = square(t[q] - p[v[k]]) + f[v[k]];
    }

    delete [] v;
    delete [] z;

    //return d;
}
/* dt of 1d function using squared distance with no function on p
 * Parameters:
 * p: location of valid parabolas: sorted
 * n: size of p or f
 * stride: step size
 * t: locations needed for distance transform: sorted
 * t_n: size of t
 */
void dt_noF(int *p,
        size_t n,
        float *d,
        const size_t stride,
        float *t,
        const size_t t_n)
{
    int *v = new int[n];
    float *z = new float[n+1];
    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] = +INFINITY;
    float s;
    for (int q = 1; q < n; q++) {
        s  = (float)(p[q] + p[v[k]])/2;
        while (s <= z[k]) {
            k--;
            s  = (float)(p[q] + p[v[k]])/2;
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k+1] = +INFINITY;
    }
    k = 0;
    for (int q = 0; q < t_n; q++) {
        while (z[k+1] < t[q])
            k++;
        d[q * stride] = square(t[q] - p[v[k]]);
    }

    delete [] v;
    delete [] z;
}

/* dt of 1d faster function using squared distance with no function on p
 * Parameters:
 * p: location of valid parabolas: sorted
 * n: size of p or f
 * stride: step size
 * t: locations needed for distance transform: sorted
 * t_n: size of t
 */
void dt_noF_faster(int *p,
        size_t n,
        float *d,
        const size_t stride,
        float *t,
        const size_t t_n)
{
    int k = 0;
    for (int q = 0; q < t_n; q++) {
        if (t[q] <= p[0])
            d[q * stride] = square(t[q] - p[0]);
        else if (t[q] >= p[n-1])
            d[q * stride] = square(t[q] - p[n-1]);
        else{
            while (p[k+1] < t[q])
                k++;
            if ( (t[q] - p[k]) > (p[k+1] - t[q]) )
                d[q * stride] = square(t[q] - p[k+1]);
            else
                d[q * stride] = square(t[q] - p[k]);
        }
    }
}

/*distance transform of 3d image */
void dt3d(bool *ref_cell, size_t *s1_yxz, size_t *s2_yxz, float *shift, float *d2map_zxy)
{
    /*printf("Reference cell: \n");
    for (int i=0; i<s1_yxz[0]; i++){
        for (int j=0; j<s1_yxz[1]; j++){
            for (int k=0; k<s1_yxz[2]; k++)
                printf("%d ", ref_cell[k*s1_yxz[0]*s1_yxz[1] + j*s1_yxz[1] + i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("Moving cell: \n");
    for (int i=0; i<s2_yxz[0]; i++){
        for (int j=0; j<s2_yxz[1]; j++){
            for (int k=0; k<s2_yxz[2]; k++)
                printf("%d ", mov_cell[k*s2_yxz[0]*s2_yxz[1] + j*s2_yxz[1] + i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("shift: %f %f %f \n", shift[0], shift[1], shift[2]);
    mexErrMsgTxt("debug use.");*/
    /***distance transform along z-direction***/
    float *d2map_z = new float[s1_yxz[0] * s1_yxz[1] * s2_yxz[2]];

    float *target_idx = new float[s2_yxz[0] + s2_yxz[1] + s2_yxz[2]];

    for (int i=0; i<s2_yxz[2]; i++){
        target_idx[i] = i + shift[2];
        //printf("%f ", target_idx[i]);
    }
    //printf("\n");
    //mexErrMsgTxt("debug use.");
    int *parabolas_mu = new int[s1_yxz[0] + s1_yxz[1] + s1_yxz[2]];

    const size_t numEle_s1xy = s1_yxz[0] * s1_yxz[1];
    for (int y = 0; y < s1_yxz[0]; y++){
        for (int x = 0; x < s1_yxz[1]; x++){
            size_t mu_cnt = 0;
            for (int z = 0; z < s1_yxz[2]; z++){
                if (ref_cell[z*numEle_s1xy + x*s1_yxz[0] + y]){
                    parabolas_mu[mu_cnt] = z;
                    mu_cnt++;
                }
            }
            if (mu_cnt > 0){
                //d = edt_1dMex(parabolas_mu,  target_idx);
                //dt_noF(parabolas_mu,
                dt_noF_faster(parabolas_mu,
                            mu_cnt,
                            d2map_z + x*s1_yxz[0] + y,
                            numEle_s1xy,
                            target_idx,
                            s2_yxz[2]);
            }else{
                /*int *q = d2map_z + x*s1_yxz[0] + y;
                for (int tmp_z = 0; tmp_z < s2_yxz[2]; tmp_z++) {
                    q[tmp_z * numEle_s1xy] = -1;
                }*/
                toNan(d2map_z + x*s1_yxz[0] + y,
                        s2_yxz[2],
                        numEle_s1xy);
            }

            /*for (int k=0; k<s2_yxz[2]; k++){
                printf("%f ", d2map_z[x*s1_yxz[0] + y + k*numEle_s1xy]);
            }
            printf("\n");*/
        }
    }
    /*printf("d2map_z: \n");
    for (int k=0; k<s2_yxz[2]; k++){
        for (int i=0; i<s1_yxz[0]; i++){
            for (int j=0; j<s1_yxz[1]; j++)
                printf("%.2f ", d2map_z[k*s1_yxz[0]*s1_yxz[1] + j*s1_yxz[0] + i]);
            printf("\n");
        }
        printf("\n");
    }
    mexErrMsgTxt("debug use.");*/

    /***distance transform along x-direction***/
    float *d2map_zx = new float[s1_yxz[0] * s2_yxz[1] * s2_yxz[2]];

    //float *target_idx = new float[s2_yxz[1]];
    for (int i=0; i<s2_yxz[1]; i++)
        target_idx[i] = i + shift[1];

    //int *parabolas_mu = new int[s1_yxz[1]];
    float *f = new float[s1_yxz[0] + s1_yxz[1]];

    const size_t numEle_s1ys2x = s1_yxz[0] * s2_yxz[1];
    for (int y = 0; y < s1_yxz[0]; y++){
        for (int z = 0; z < s2_yxz[2]; z++){
            size_t mu_cnt = 0;
            for (int x = 0; x < s1_yxz[1]; x++){

                float tmp = d2map_z[z*numEle_s1xy + x*s1_yxz[0] + y];
                if (tmp >= 0){
                    f[mu_cnt] = tmp;
                    parabolas_mu[mu_cnt] = x;
                    mu_cnt++;
                }
            }
            if (mu_cnt > 0){
                dt(parabolas_mu,
                        f,
                        mu_cnt,
                        d2map_zx + z*numEle_s1ys2x + y,
                        s1_yxz[0],
                        target_idx,
                        s2_yxz[1]);
            }else{
                /*int *q = d2map_zx + z*numEle_s1xy + y;
                for (int tmp_z = 0; tmp_z < s1_yxz[0]; tmp_z++) {
                    q[tmp_z * s1_yxz[0]] = -1;
                }*/
                toNan(d2map_zx + z*numEle_s1ys2x + y,
                        s2_yxz[1],
                        s1_yxz[0]);
            }
        }
    }
    /*printf("d2map_zx: \n");
    for (int k=0; k<s2_yxz[2]; k++){
        for (int i=0; i<s1_yxz[0]; i++){
            for (int j=0; j<s2_yxz[1]; j++)
                printf("%.2f ", d2map_zx[k*s1_yxz[0]*s2_yxz[1] + j*s1_yxz[0] + i]);
            printf("\n");
        }
        printf("------------# %d--------------\n", k+1);
    }
    mexErrMsgTxt("debug use.");*/

    /***distance transform along y-direction ==> output***/
    //float *d2map_zxy = (float *)mxGetPr(plhs[0]);

    //float *target_idx = new float[s2_yxz[0]];
    for (int i=0; i<s2_yxz[0]; i++)
        target_idx[i] = i + shift[0];

    size_t mu_cnt = 0;
    //int *parabolas_mu = new int[s1_yxz[0]];
    for (int i=0; i<s1_yxz[0]; i++){
        if (d2map_zx[i] >= 0){
            parabolas_mu[mu_cnt] = i;
            mu_cnt++;
        }
    }
    //float *f = new float[s1_yxz[0]];
    //printf("%d %d %f\n", s2_yxz[0], shift[0]);
    //mexErrMsgTxt("debug use.");
    const size_t numEle_s2yx = s2_yxz[0] * s2_yxz[1];
    //printf("%d %d %d %f\n", s2_yxz[0], s2_yxz[1], s2_yxz[2]);
    for (int x = 0; x < s2_yxz[1]; x++){
        for (int z = 0; z < s2_yxz[2]; z++){
            for (int y = 0; y < mu_cnt; y++){
                f[y] = d2map_zx[z * numEle_s1ys2x
                        + x * s1_yxz[0]
                        + parabolas_mu[y]];
            }
            dt(parabolas_mu,
                    f,
                    mu_cnt,
                    d2map_zxy + z*numEle_s2yx + x * s2_yxz[0],
                    1,
                    target_idx,
                    s2_yxz[0]);
        }
    }

    /*printf("d2map_zxy: \n");
    for (int k=0; k<s2_yxz[2]; k++){
        for (int i=0; i<s2_yxz[0]; i++){
            for (int j=0; j<s2_yxz[1]; j++)
                printf("%.2f ", d2map_zxy[k*s2_yxz[0]*s2_yxz[1] + j*s2_yxz[0] + i]);
            printf("\n");
        }
        printf("\n");
    }
    mexErrMsgTxt("debug use.");*/

    delete [] d2map_z;
    delete [] d2map_zx;
    delete [] target_idx;
    delete [] parabolas_mu;
    delete [] f;

}
#endif // DT3D_HPP
