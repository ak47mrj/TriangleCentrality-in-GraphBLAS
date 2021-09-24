//
// Created by Fuhuan Li
//

#include <stdio.h>
#include <time.h>
#include "GraphBLAS.h"
#include "../include/utils.h"

GrB_Info triangleCentrality(GrB_Vector *result, GrB_Matrix A);

int main(int argc,char** argv) {

    GrB_init(GrB_NONBLOCKING);

    GrB_Matrix graph;
    FILE *fp;
    fp = fopen("../graphs/testgraph.txt", "r");

    readMatrix(&graph, fp, true);
    fclose(fp);

    GrB_Vector result;

    struct timespec start, end;
    uint64_t duration;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    triangleCentrality(&result, graph);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    duration = (end.tv_sec - start.tv_sec) * 1000000000 + (end.tv_nsec - start.tv_nsec);
    printf("%llu\n", duration);

    GrB_free(&result);
    GrB_free(&graph);
    GrB_finalize();

    return 0;

}