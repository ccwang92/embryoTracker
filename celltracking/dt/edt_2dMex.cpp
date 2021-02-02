
/* This function is a c++ version of edt_2d_regFixed.m
 * Ideas are from Felzenszwalb and Huttenlocher's paper about distance 
 * transform. Our distance transfrom follows the x->y function. 
 * More info can be found in edt_2d_regFixed.m.
 * 
 * Author: Congchao Wang
 * Yu's Lab, Virginia Tech
 * Date: September, 2020
 * 
 */
#include <stdio.h>
#include <mex.h>
#include <algorithm>
#include <limits>

const int INFINITY = std::numeric_limits<int>::max();

float square(float x) {
    return x*x;
}

void toNan(float *f, 
        const mwSize n, 
        const mwSize stride){
    for (int i = 0; i < n; i++) {
        f[i * stride] = -1;
    }
}
/* dt of 1d function using squared distance 
 * with function f defined.
 * Parameters:
 *
 */
void dt(int *p, 
        float *f, 
        mwSize n, 
        float *d, 
        const mwSize stride, /*for d*/ 
        float *t, 
        const mwSize t_n) {
    
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
/* dt of 1d function using squared distance 
 * with function f undefined (all are 0s).
 * Parameters:
 *
 */
void dt_noF(int *p, 
        mwSize n, 
        float *d, 
        const mwSize stride, 
        float *t, 
        const mwSize t_n) {
    
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
/*distance transform of 2d image */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     if (nlhs > 1 || nrhs != 3)
        mexErrMsgTxt("Syntax:\n\tH = edt_2dMex(p,f,n,t,t_n)");

    if (mxGetClassID(prhs[0]) != mxLOGICAL_CLASS 
             || mxGetClassID(prhs[1]) != mxLOGICAL_CLASS
             || mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
        mexErrMsgTxt("inputs are: 2d_logic, 2d_logic, array_single.");
    
    //int *p->parabolas_mu, float *f, int n, float *t, int t_n
    bool *ref_cell = (bool*)mxGetPr(prhs[0]);
    const mwSize *s1_yx = mxGetDimensions(prhs[0]); // 1x2;

    
    bool *mov_cell = (bool*)mxGetPr(prhs[1]);
    const mwSize *s2_yx = mxGetDimensions(prhs[1]); // 1x2
    
    float *shift = (float*)mxGetPr(prhs[2]); // 1x2

    /***distance transform along x-direction***/
    float *d2map_x = new float[s1_yx[0] * s2_yx[1]];
    
    float *target_idx = new float[s2_yx[0] + s2_yx[1]];
    
    for (int i=0; i<s2_yx[1]; i++){
        target_idx[i] = i + shift[1];
        //printf("%f ", target_idx[i]);
    }
    //printf("\n");
    //mexErrMsgTxt("debug use.");
    int *parabolas_mu = new int[s1_yx[0] + s1_yx[1]];
    
    for (int y = 0; y < s1_yx[0]; y++){
        mwSize mu_cnt = 0;
        for (int x = 0; x < s1_yx[1]; x++){
            if (ref_cell[x*s1_yx[0] + y]){
                parabolas_mu[mu_cnt] = x;
                mu_cnt++;
            }
        }
        if (mu_cnt > 0){
            //d = edt_1dMex(parabolas_mu,  target_idx);
            dt_noF(parabolas_mu,
                    mu_cnt,
                    d2map_x + y,
                    s1_yx[0],
                    target_idx,
                    s2_yx[1]);
        }else{
            /*int *q = d2map_z + x*s1_yx[0] + y;
             * for (int tmp_z = 0; tmp_z < s2_yx[2]; tmp_z++) {
             * q[tmp_z * numEle_s1xy] = -1;
             * }*/
            toNan(d2map_x + y,
                    s2_yx[1],
                    s1_yx[0]);
        }
    }
    
    /*printf("d2map_x: \n");
    for (int i=0; i<s1_yx[0]; i++){
        for (int j=0; j<s2_yx[1]; j++)
            printf("%.2f ", d2map_x[j*s1_yx[0] + i]);
        printf("\n");
    }
    mexErrMsgTxt("debug use.");*/

    /***distance transform along y-direction ==> output***/
    plhs[0] = mxCreateNumericMatrix(1, s2_yx[0] * s2_yx[1],
            mxSINGLE_CLASS,mxREAL);
    
    float *d2map_xy = (float *)mxGetPr(plhs[0]);
    
    //float *target_idx = new float[s2_yx[0]];
    for (int i=0; i<s2_yx[0]; i++)
        target_idx[i] = i + shift[0];
    
    mwSize mu_cnt = 0;
    //int *parabolas_mu = new int[s1_yx[0]];
    for (int i=0; i<s1_yx[0]; i++){
        if (d2map_x[i] >= 0){
            parabolas_mu[mu_cnt] = i;
            mu_cnt++;
        }
    }
    float *f = new float[s1_yx[0]];
    //printf("%d %d %f\n", s2_yx[0], shift[0]);
    //mexErrMsgTxt("debug use.");
    for (int x = 0; x < s2_yx[1]; x++){
        for (int y = 0; y < mu_cnt; y++){
            f[y] = d2map_x[x * s1_yx[0] + parabolas_mu[y]];
        }
        dt(parabolas_mu,
                f,
                mu_cnt,
                d2map_xy + x * s2_yx[0],
                1,
                target_idx,
                s2_yx[0]);
    }
    
    /*printf("d2map_xy: \n");
    for (int i=0; i<s2_yx[0]; i++){
        for (int j=0; j<s2_yx[1]; j++){
            printf("%.2f ", d2map_xy[j*s2_yx[0] + i]);
        }
        printf("\n");
    }
    mexErrMsgTxt("debug use.");*/
    
    delete [] d2map_x;
    delete [] target_idx;
    delete [] parabolas_mu;
    delete [] f;
        
}