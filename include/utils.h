//
// Created by Fuhuan Li on 7/4/21.
//

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void print_vector(GrB_Vector vec, char const *label)
{
    GrB_Index N;
    GrB_Vector_size(&N, vec);

    double val;
    GrB_Info ret_val;

    printf("Vector: %s =\n", label);

    printf("[");
    for (GrB_Index idx = 0; idx < N; ++idx)
    {
        ret_val = GrB_Vector_extractElement_FP64(&val, vec, idx);
        if (GrB_SUCCESS == ret_val)
        {
            if (idx == 0)
            {
                printf("%lf", (double)val);
            }
            else
            {
                printf(", %lf", (double)val);
            }
        }
        else if (GrB_NO_VALUE == ret_val)
        {
            if (idx == 0)
            {
                printf("  -");
            }
            else
            {
                printf(",   -");
            }
        }
        else
        {
            if (idx == 0)
            {
                printf("  ERR");
            }
            else
            {
                printf(", ERR");
            }
        }
    }
    printf("]\n");

}

GrB_Info readMatrix(GrB_Matrix *graph, FILE *f, bool one_based) {
    int64_t len = 256;
    int64_t ntuples = 0;
    void *X1 = malloc (len * sizeof(double)) ;
    double *Xdouble ;
    void *X2 = NULL ;

    GrB_Index *I = (GrB_Index *) malloc(len * sizeof(GrB_Index)), *I2 = NULL;
    GrB_Index *J = (GrB_Index *) malloc(len * sizeof(GrB_Index)), *J2 = NULL;
    Xdouble = (double *) X1;

    double i2, j2;
    while (fscanf(f, "%lg %lg\n", &i2, &j2) != EOF) {
        int64_t i = (int64_t) i2;
        int64_t j = (int64_t) j2;
        if (ntuples >= len) {
            I2 = (GrB_Index *) realloc(I, 2 * len * sizeof(GrB_Index));
            J2 = (GrB_Index *) realloc(J, 2 * len * sizeof(GrB_Index));
            X2 = realloc(X1, 2 * len * sizeof(double));
            I = I2;
            I2 = NULL;
            J = J2;
            J2 = NULL;
            X1 = X2;
            X2 = NULL;
            len = len * 2;
            Xdouble = (double *) X1;
        }
        if (one_based) {
            i--;
            j--;
        }
        I[ntuples] = i;
        J[ntuples] = j;
        Xdouble[ntuples] = 1.0;
        ntuples++;
    }

    int64_t nrows = 0;
    int64_t ncols = 0;

    for (int64_t kk = 0; kk < ntuples; kk++) {
        nrows = MAX(nrows, I[kk]);
        ncols = MAX(ncols, J[kk]);
    }
    nrows = MAX(nrows, ncols);
    nrows += 1;

    GrB_Matrix_new(graph, GrB_FP64, nrows, nrows);
    GrB_Matrix_build(*graph, I, J, Xdouble, ntuples,
                     GrB_LOR);
    GrB_Descriptor dt2;
    GrB_Descriptor_new(&dt2);
    GrB_Descriptor_set(dt2, GrB_INP1, GrB_TRAN);
    GrB_Matrix_eWiseAdd_BinaryOp(*graph, NULL, NULL, GrB_LOR, *graph, *graph, dt2);

    free(I);
    free(J);
    free(X1);
    free(I2);
    free(J2);
    free(X2);
    GrB_free(&dt2);
    return GrB_SUCCESS;
}