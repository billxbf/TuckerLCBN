#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <mkl.fi>
//#include <blas.f90>


double getMean(double* arr)
{
    double sum = 0;
    int size = sizeof(arr)/sizeof(arr[0]);
    for (int i = 0; i < size; i++){
        sum += arr[i];
    }
    return sum / size;
}

double getNorm(double* arr)
{
    double norm = 0;
    int size = sizeof(arr)/sizeof(arr[0]);
    for(int i = 0; i < size; i++){
        norm += arr[i] * arr[i];
    }
    return sqrt(norm);
}

int main(){

    long rowcol[2];
    long numRow;
    long numCol;

    //Read in metadata
    //rowcol = (long *)malloc(sizeof(long)*2);
    FILE* fp = fopen("metadata.txt","r");
    for (int i=0; i < 2; i++){
        fscanf(fp, "%d", &rowcol[i]);
    }

    printf("%d", rowcol[0]);


//    numRow = rowcol[0];
//    numCol = rowcol[1];
//
//    double* myMatrix;
//    double* result;
//
//    //Read in Matrix
//    myMatrix = (double *)malloc(sizeof(double)*numRow*numCol);
//    result = (double *)malloc(sizeof(double)*numCol*numCol);
//
//    FILE* fp2 = fopen("matrix.bin", "r");
//    fread(myMatrix, sizeof(double), numRow*numCol, fp2);
//    printf("fuckkkkkk\n");
//    printf("############");
//    //Calculate a - mean(a) for each column a
//    for (int j = 0; j < numCol; j++){
//        double* col;
//        col = (double *)malloc(sizeof(double)*numRow);
//        for (int i = 0; i < numRow; i++){
//            col[i] =  myMatrix[numRow*j + i];
//        }
//        // Calculate mean(a)
//        double mean = getMean(col);
//        //  Calculate a = a - mean(a)
//        for(int i = 0; i < numRow; i++){
//            myMatrix[numRow*j + i] = myMatrix[numRow*j + i] - mean;
//        }
//    }
//
//    //Normalization
//    for (int j = 0; j < numCol; j++){
//        double* col;
//        col = (double *)malloc(sizeof(double)*numRow);
//        for (int i = 0; i < numRow; i++){
//            col[i] =  myMatrix[numRow*j + i];
//        }
//        // Calculate norm(a)
//        double norm = getNorm(col);
//        //  Calculate a = a/norm(a)
//        for(int i = 0; i < numRow; i++){
//            myMatrix[numRow*j + i] = myMatrix[numRow*j + i] / norm;
//        }
//    }

    //Calculate Correlation Matrix using CBLAS (X' * X)
    //dsyrk('L','C', numRow, numCol, 1.0, &myMatrix, numRow, 0.0, result, numCol);


}

