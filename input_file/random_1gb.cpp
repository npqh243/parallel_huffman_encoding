// 1024 MB file of random ascii characters
#include <random>
#include <fstream>
#include <iostream>
#include <omp.h>

using namespace std;

int main() {
    FILE* f = fopen("random_1gb.huffinp", "w");
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(32, 126);
    const long long FILE_SIZE = 1024LL * 1024 * 1024;
    #pragma omp parallel for 
    for (long long i = 0; i < FILE_SIZE; i++) {
        fputc(dis(gen), f);
    }
    fclose(f);
    cout << "done";
}