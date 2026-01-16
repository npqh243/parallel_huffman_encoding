#include<vector>
#include<iostream>
#include<omp.h>

using namespace std;

// Dot product of two vectors parrallelized
double vector_dot_product(vector<double> vector1, vector<double> vector2, int n) {
    double sum = 0;

    #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < n; i++) {
            sum += vector1[i] * vector2[i];
        }

    return sum;
}


int main() {
    const int vectorSize = 10000000;
    vector<double> vector1, vector2;
    vector1.resize(vectorSize);
    vector2.resize(vectorSize);
    double result;
    
    // initialize vectors for different values
    /*
    #pragma omp parallel for
        for (int i = 0; i < vectorSize; i++) {
            vector1[i] = i * 1.0;
            vector2[i] = (vectorSize - i) * 1.0;
        }
    */
    /*
    #pragma omp parallel for
        for (int i = 0; i < vectorSize; i++) {
            vector1[i] = 1.0;
            vector2[i] = 1.0;
        }
    */

    double time = omp_get_wtime();
    result = vector_dot_product(vector1, vector2, vectorSize);
    time = omp_get_wtime() - time;
    
    printf("Time taken: %.5lf seconds\n", time);
    printf("\nDot product: %.2lf\n", result);

    return 0;
}