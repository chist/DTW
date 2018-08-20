/* Chistiakov Ivan, 2018 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"

typedef struct Minstruct
{
    double value;
    int index;
} Minstruct;

/* find minimum of 3 double values and return index (1 to 3) */
Minstruct min3(double a, double b, double c)
{
    Minstruct res;
    if (a <= b && a <= c) {
        res.value = a;
        res.index = 1;
    } else if (b <= a && b <= c) {
        res.value = b;
        res.index = 2;
    } else {
        res.value = c;
        res.index = 3;
    }
    return res;
}

/* create matrix of all pairwise distances between element of x and y */
double *create_d_matrix(double *x, double *y, int xcols, int ycols, int bound)
{
    double *d_matrix = calloc(xcols * ycols, sizeof(double));
    if (d_matrix == NULL) {
        printf("%d %d\n", xcols, ycols);
        mexErrMsgIdAndTxt("DTW:memory", "Calloc error.");
    }
    
    for (int j = 0; j < ycols; j++) {
        for (int i = 0; i < xcols; i++) {
            // skip some cells to reduce execution time
            if (abs(i - j) > bound) {
                continue;
            }
            
            // calculate euclidian distance
            double sq_dist = (x[i] - y[j]) * (x[i] - y[j]) + 
                (x[xcols + i] - y[ycols + j]) * (x[xcols + i] - y[ycols + j]);
            d_matrix[j * xcols + i] = sqrt(sq_dist);
        }
    }
    
    return d_matrix;
}

void fill_matrices(double *t_matrix, int *p_matrix, double *d_matrix, 
        int xcols, int ycols, int bound)
{
    t_matrix[0] = d_matrix[0];
    for (int i = 1; i < xcols; i++) {
        t_matrix[i] = t_matrix[i - 1] + d_matrix[i];
        p_matrix[i] = 1;
    }
    for (int j = 1; j < ycols; j++) {
        t_matrix[j * xcols] = t_matrix[(j - 1) * xcols] + d_matrix[j * xcols];
        p_matrix[j * xcols] = 3;
    }
    for (int j = 1; j < ycols; j++) {
        for (int i = 1; i < xcols; i++) {
            // skip some cells to reduce execution time
            if (abs(i - j) > bound) {
                continue;
            }
            
            // select neighbours to the current cell
            double el1 = t_matrix[j * xcols + i - 1];
            double el2 = t_matrix[(j - 1) * xcols + i - 1];
            double el3 = t_matrix[(j - 1) * xcols + i];
            
            Minstruct m;
            if (abs(i - j) < bound) {
                m = min3(el1, el2, el3);
            } else {
                // case of (abs(i - j) == bound)
                
                if (i >= j) {
                    m = min3(el1, el2, DBL_MAX);
                } else {
                    m = min3(DBL_MAX, el2, el3);
                }
            }
            t_matrix[j * xcols + i] = d_matrix[j * xcols + i] + m.value;
            p_matrix[j * xcols + i] = m.index;
        }
    }
}

int find_optimal_way(int *p_matrix, int xcols, int ycols, int *ix, int *iy)
{
    int k = 0;
    ix[k] = xcols;
    iy[k] = ycols;
    int i = xcols;
    int j = ycols;
    
    while (i + j > 2) {
        switch (p_matrix[(j - 1) * xcols + (i - 1)]) {
            case 1:
                i = i - 1;
                break;
            case 2:
                i = i - 1;
                j = j - 1;
                break;
            case 3:
                j = j - 1;
                break;
            default:
                i = 1;
                j = 1;
        }
        k = k + 1;
        ix[k] = i;
        iy[k] = j;
    }
    return k + 1;
}

double find_dtw(double *x, double *y, int window_size, int xcols, int ycols,
        int *ix, int *iy, int *length)
{
    // set restriction for the DTW algorithm to reduce execution time
    int bound = abs(xcols - ycols) + window_size;
    
    // create matrix of pairwise euclidian distances
    double *d_matrix = create_d_matrix(x, y, xcols, ycols, bound);
    
    // allocate transformation matrix
    double *t_matrix = calloc(xcols * ycols, sizeof(double));
    if (t_matrix == NULL) {
        mexErrMsgIdAndTxt("DTW:memory", "Calloc error");
    }
    
    // allocate matrix of directions
    int *p_matrix = calloc(xcols * ycols, sizeof(int));
    if (p_matrix == NULL) {
        mexErrMsgIdAndTxt("DTW:memory", "Calloc error");
    }
    /* possible directions of the #-cell are the following:
     *     23
     *     1#      */
    
    // fill distance and direction matrices
    fill_matrices(t_matrix, p_matrix, d_matrix, xcols, ycols, bound);
    
    // find optimal path in p_matrix and return path length
    *length = find_optimal_way(p_matrix, xcols, ycols, ix, iy);
    
    // calculate normalized DTW distance
    double dist = t_matrix[ycols * xcols - 1] / *length;
    
    // free allocated memory
    free(d_matrix);
    free(t_matrix);
    free(p_matrix);
    
    return dist;
}

/* remove rows from array of (length)*2 size and save new
 * length with (*length_ptr) pointer */
double *reduce_array(double *ptr, int *length_ptr, int coef)
{
    // get array dimension
    int length = *length_ptr;
    
    // calculate linear size of reduced array
    int new_length = length / coef;
    
    // save reduced array's dimension
    *length_ptr = new_length;
    
    // allocate memory
    double *new_ptr = calloc(2 * new_length, sizeof(double));
    if (new_ptr == NULL) {
        mexErrMsgIdAndTxt("DTW:memory", "Calloc error.");
    }
    
    // copy data
    for (int i = 0; i < new_length; i++) {
        new_ptr[i] = ptr[i * coef];
    }
    for (int i = 0; i < new_length; i++) {
        new_ptr[new_length + i] = ptr[length + i * coef];
    }
    
    return new_ptr;
}

int reduce_data(double **x, double **y, int *xcols, int *ycols)
{   
    long long x_size = *xcols, y_size = *ycols;
    if (x_size * y_size > (long long)INT_MAX) {   
        int coef = (int)sqrt((double)(x_size * y_size / 
                (long long)INT_MAX)) + 1;
        mexWarnMsgIdAndTxt("DTW:shrink", 
                "Arrays were reduced to satisfy memory limits.");
        *x = reduce_array(*x, xcols, coef);
        *y = reduce_array(*y, ycols, coef);
        return coef;
    } else {
        return 1;
    }
}

void check_arguments(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("DTW:nrhs", "Two inputs required.");
    } else if (nrhs > 3) {
        mexErrMsgIdAndTxt("DTW:nrhs", "Too many input arguments.");
    }
    if (nlhs < 1) {
        mexErrMsgIdAndTxt("DTW:nlhs", "Output required.");
    } else if (nlhs > 3) {
        mexErrMsgIdAndTxt("DTW:nlhs", "Too many output arguments.");
    }
    
    if (mxGetN(prhs[0]) != 2 || mxGetN(prhs[1]) != 2) {
        mexErrMsgIdAndTxt("DTW:wrongMatrixSize",
                "Input matrices must have two columns.");
    }
    if (nrhs >= 3 && (mxGetM(prhs[2]) != 1 || mxGetM(prhs[2]) != 1)) {
        mexErrMsgIdAndTxt("DTW:wrongMatrixSize",
                "Window size parameter must be scalar.");
    }
}
    
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    check_arguments(nlhs, plhs, nrhs, prhs);

    /* create pointers to the real data in the input matrices;
     * each matrix columns contain x and y coordinates of curve's points */
    double *x_matrix = mxGetPr(prhs[0]);
    double *y_matrix = mxGetPr(prhs[1]);
    /* create pointers to store the reduced arrays (if necessary) */
    double *x = x_matrix, *y = y_matrix;

    /* get lengths of the input matrices */
    int xcols = (int)mxGetM(prhs[0]);
    int ycols = (int)mxGetM(prhs[1]);
    
    int reduce_coef = reduce_data(&x, &y, &xcols, &ycols);
    
    /* get the value of the window size or set it to the default value  */
    int window_size;
    if (nrhs < 3) {
        window_size = min(xcols, ycols);
    } else {
        window_size = (int)mxGetScalar(prhs[2]);
        if (window_size > min(xcols, ycols)) {
            window_size = min(xcols, ycols);
            mexWarnMsgIdAndTxt("DTW:window", 
                "Window size too large! It was reduced.");
        }
    }
    
    /* allocate fixed-size arrays for the indices */
    int fixed_size = 2 * (xcols + ycols);
    int *ix_long = calloc(fixed_size, sizeof(int));
    int *iy_long = calloc(fixed_size, sizeof(int));
    
    /* allocate memory for the optimal sequence length to pass indices 
     * correctly to the function output */ 
    int length = 0;
    
    double dist;
    dist = find_dtw(x_matrix, y_matrix, window_size, xcols, ycols, 
            ix_long, iy_long, &length);
    
    /* pass results to output */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    double *data = mxGetPr(plhs[0]);
    data[0] = dist;
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericMatrix(1, length, mxINT32_CLASS, mxREAL);
        int32_T *ix = (int32_T *)mxGetData(plhs[1]);
        for (int i = 0; i < length; i++) {
            ix[i] = (int32_T)(ix_long[i] * reduce_coef);
        }
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateNumericMatrix(1, length, mxINT32_CLASS, mxREAL);
        int32_T *iy = (int32_T *)mxGetData(plhs[2]);
        for (int i = 0; i < length; i++) {
            iy[i] = (int32_T)(iy_long[i] * reduce_coef);
        }
    }
    
    /* free allocated memory */
    if (reduce_coef != 1) {
        free(x);
        free(y);
    }
    free(ix_long);
    free(iy_long);
}
