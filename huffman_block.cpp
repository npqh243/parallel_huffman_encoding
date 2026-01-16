#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <omp.h>     
#include <queue>
#include <iomanip>
#include "huffman.cpp"

using namespace std;

int main() {
    cout << "--- MA HOA HUFFMAN (BLOCK) ---\n";
    
    omp_set_num_threads(THREADS_NUM);
    
    FILE* inpFile = fopen(FILE_NAME, "r");
    if (!inpFile) {
        cerr << "Khong the mo file input\n";
        return 1;
    }

    fseek(inpFile, 0, SEEK_END);
    long long fileSize = ftell(inpFile);
    fseek(inpFile, 0, SEEK_SET);

    double time0 = omp_get_wtime();

    // --- BƯỚC 1: ĐẾM TẦN SUẤT  ---
    long long freq[256] = { 0 };
    countFrequencyBlock(inpFile, fileSize, freq);

    double time1 = omp_get_wtime();
    cout << "Thoi gian chay:        " << time1 - time0 << " giay\n";

    // --- BƯỚC 2: DỰNG CÂY HUFFMAN ---
    map<char, string> huffmanCode;
    generateNaiveCode(freq, huffmanCode);

    double time2 = omp_get_wtime();
    cout << "Thoi gian chay:        " << time2 - time1 << " giay\n";

    // --- BƯỚC 3: TÍNH KÍCH THƯỚC SAU NÉN ---
    fseek(inpFile, 0, SEEK_SET);

    long long totalBits = 0;
    encodeBlock(inpFile, fileSize, totalBits, huffmanCode);
    
    double time3 = omp_get_wtime();
    cout << "Thoi gian chay:        " << time3 - time2 << " giay\n";

    fclose(inpFile);

    // --- KẾT QUẢ ---
    cout << "=== KET QUA BLOCK HUFFMAN ===\n";
    cout << "Tong thoi gian chay:   " << time3 - time0 << " giay\n";

    double originalSizeMB = (double)fileSize / 1024 / 1024;    // / KB / MB
    double compressedSizeMB = (double)((totalBits / 8) + 1) / 1024 / 1024;  // / KB / MB
    double blockSizeMB = (double)BLOCK_SIZE / 1024 / 1024;  // / KB / MB

    cout << fixed << setprecision(2);
    cout << "So luong song song:    " << THREADS_NUM << "\n";
    cout << "Kich thuoc block:      " << blockSizeMB << " MB\n";
    cout << "Kich thuoc TRUOC nen:  " << originalSizeMB << " MB\n";
    cout << "Kich thuoc SAU nen:    " << compressedSizeMB << " MB\n";
    cout << "Ti le nen:             " << (compressedSizeMB / originalSizeMB) * 100 << "%\n";
    cout << "Toc do xu ly:          " << originalSizeMB / (time3 - time0) << " MB/s\n";

    return 0;
}