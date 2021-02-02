#include <stdio.h>
#include "mex.h"
//#include "matrix.h"
#include <algorithm>
#include <limits>

const int INFINITY = std::numeric_limits<int>::max();

/*distance transform of 1d image */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     if (nlhs > 1 || nrhs != 3)
        mexErrMsgTxt("Syntax:\n\tH = edt_1dMex(p,f,n,t,t_n)");

    if (mxGetClassID(prhs[0]) != mxUINT16_CLASS)
        mexErrMsgTxt("inputs are: unit16, single, int, single, int.");

    
    //int *p->parabolas_mu, float *f, int n, float *t, int t_n
    uint16_T *p = (uint16_T *)mxGetPr(prhs[0]);
    float *f = (float*)mxGetPr(prhs[1]);
    int n = mxGetNumberOfElements(prhs[0]);
    float *t = (float*)mxGetPr(prhs[2]);
    int t_n = mxGetNumberOfElements(prhs[2]);
    
    //for (int q = 0; q < n; q++) {
    //    if (f[q] > 0.00001)
    //       printf("%d %f", p[q], f[q]);
    //}
    //printf("\n");
    
    int *v = new int[n];
    float *z = new float[n+1];
    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] = +INFINITY;
    float s;
    for (int q = 1; q < n; q++) {
        s  = (f[q] - f[v[k]])/(2*p[q]-2*p[v[k]]) + (float)(p[q] + p[v[k]])/2;
        //if (p[q] == 4102)
        //    printf("%d %d %f->", p[q]*p[q] - p[v[k]]*p[v[k]], 2*p[q]-2*p[v[k]], s);
        while (s <= z[k]) {
            k--;
            s  = (f[q] - f[v[k]])/(2*p[q]-2*p[v[k]]) + (float)(p[q] + p[v[k]])/2;
            
            //if (p[q] == 4102)
            // printf("%d %f %f %f ->", p[v[k]], f[q], f[v[k]], s);
        }
        //if (p[q] == 4102)
        //    printf("\n");
        k++;
        v[k] = q;
        z[k] = s;
        z[k+1] = +INFINITY;
    }
    /*for (int q = 0; q < n; q++) {
        printf("%d: %f->%f \n", p[v[q]], z[q], z[q+1]);
        if (z[q+1] == +INFINITY)
            break;
    }*/
    
    plhs[0] = mxCreateNumericMatrix(1,t_n,mxSINGLE_CLASS,mxREAL);

    float *d = (float *)mxGetPr(plhs[0]);
    
    k = 0;
    for (int q = 0; q < t_n; q++) {
        while (z[k+1] < t[q])
            k++;
        d[q] = (t[q] - p[v[k]])*(t[q] - p[v[k]]) + f[v[k]];
    }
    
    delete [] v;
    delete [] z;
}