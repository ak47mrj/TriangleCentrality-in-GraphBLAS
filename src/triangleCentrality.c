#include <stdio.h>
#include "GraphBLAS.h"
#include <time.h>

GrB_Info triangleCentrality(GrB_Vector *result, GrB_Matrix A) {
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);

    GrB_Matrix T;
    GrB_Matrix_new(&T, GrB_FP64, n, n);
    GrB_Vector T_y;
    GrB_Vector_new(&T_y,GrB_FP64,n);
    GrB_Vector y;
    GrB_Vector_new(&y, GrB_FP64, n);
    GrB_Vector_new(result, GrB_FP64, n);
    double k;

    GrB_mxm(T, A, NULL, GxB_PLUS_TIMES_FP64, A, A, NULL);
    GrB_reduce(y, NULL, NULL, GrB_PLUS_FP64, T, NULL);
    GrB_reduce(&k, NULL, (GrB_Monoid) GrB_PLUS_MONOID_FP64, y, NULL);

    GrB_mxv(*result,NULL,NULL,GxB_PLUS_TIMES_FP64,A,y,NULL);
    GrB_apply(*result,NULL,NULL,GrB_TIMES_FP64,3,*result,NULL);
    GrB_mxv(T_y,NULL,NULL,GxB_PLUS_SECOND_FP64,T,y,NULL);
    GrB_apply(T_y,NULL,NULL,GrB_TIMES_FP64,2,T_y,NULL);
    GrB_eWiseAdd(*result,NULL,GrB_PLUS_FP64,GrB_MINUS_FP64,y,T_y,NULL);
    GrB_apply(*result, NULL, NULL, GrB_TIMES_FP64, 1 / k, *result, NULL);

    GrB_free(&T);
    GrB_free(&T_y);
    GrB_free(&y);

    return GrB_SUCCESS;

}