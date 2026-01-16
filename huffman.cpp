#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <algorithm>

using namespace std;

#define BLOCK_SIZE 1024 * 1024
#define FILE_NAME "input_file/random_1gb.huffinp"
#define THREADS_NUM 16

// ========================
// NAIVE
// ========================

struct Node {
    char ch;
    int freq;
    Node* left, * right;

    Node(char character, int frequency, Node* l = nullptr, Node* r = nullptr) {
        ch = character;
        freq = frequency;
        left = l;
        right = r;
    }
};

struct Compare {
    bool operator()(Node* l, Node* r) {
        return l->freq > r->freq;
    }
};

void countFrequencyNaive(FILE *inpFile, long long freq[]) {    
    cout << "--- [PHASE 1] DEM TAN SUAT (NAIVE) ---\n";

    while (!feof(inpFile)) {
            int ch = fgetc(inpFile);
            if (ch != EOF) {
            freq[(unsigned char)ch]++;
        }
    }

    cout << "Hoan thanh dem tan suat\n";
}

void generateCodes(Node* root, string str, map<char, string>& huffmanCode) {
    if (!root) return;
    // Nếu là lá (Leaf node) thì gán mã
    if (!root->left && !root->right) {
        huffmanCode[root->ch] = str;
    }
    generateCodes(root->left, str + "0", huffmanCode);
    generateCodes(root->right, str + "1", huffmanCode);
}

void deleteTree(Node* root) {
    if (!root) return;
    deleteTree(root->left);
    deleteTree(root->right);
    delete root;
}

void generateNaiveCode(long long freq[], map<char, string>& huffmanCode) {
    cout << "--- [PHASE 2] TAO CAU TRUC MA (NAIVE) ---\n";

    priority_queue<Node*, vector<Node*>, Compare> pq;
    for (int i = 0; i < 256; i++) {
        if (freq[i] > 0) {
            pq.push(new Node((char)i, freq[i]));
        }
    }

    while (pq.size() != 1) {
        Node* left = pq.top(); pq.pop();
        Node* right = pq.top(); pq.pop();
        int sum = left->freq + right->freq;
        pq.push(new Node('\0', sum, left, right));
    }
    Node* root = pq.top();

    generateCodes(root, "", huffmanCode);
    deleteTree(root);
    
    cout << "Hoan thanh tao ma huffman\n";
}

void encodeNaive(FILE* inpFile, long long& totalBits, map<char, string>& huffmanCode) {
    cout << "--- [PHASE 3] GIA LAP MA HOA (NAIVE) ---\n";

    while (!feof(inpFile)) {
        int ch = fgetc(inpFile);
        if (ch != EOF) {
            totalBits += huffmanCode[(char)ch].length();
        }
    }

    cout << "Hoan thanh gia lap ma hoa\n";
}

// ========================
// BLOCKING
// ========================

void countFrequencyBlock(FILE* inpFile, long long fileSize, long long freq[]) {
    cout << "--- [PHASE 1] DEM TAN SUAT (BLOCK) ---\n";

    int numBlocks = fileSize / BLOCK_SIZE;
    char* buffer = new char[BLOCK_SIZE];

    while (ftell(inpFile) - fileSize != 0) {
        long long currentBlockSize = (fileSize - ftell(inpFile) >= BLOCK_SIZE) ? BLOCK_SIZE : fileSize - ftell(inpFile);
        fread(buffer, sizeof(unsigned char), currentBlockSize, inpFile);

#pragma omp parallel
        {
            long long local_freq[256] = { 0 };
#pragma omp for nowait
            for (int i = 0; i < currentBlockSize; i++) {
                local_freq[(unsigned char)buffer[i]]++;
            }
#pragma omp critical
            {
                for (int i = 0; i < 256; i++) freq[i] += local_freq[i];
            }
        }
    }
    cout << "Hoan thanh dem tan suat!\n";
    delete[] buffer;
}

void encodeBlock(FILE* inpFile, long fileSize, long long& totalBits, map<char, string>& codeMap) {
    int numBlocks = fileSize / BLOCK_SIZE;
    char* buffer = new char[BLOCK_SIZE];

    cout << "--- [PHASE 3] GIA LAP MA HOA (BLOCK) ---\n";

    while (ftell(inpFile) - fileSize != 0) {
        long long currentBlockSize = (fileSize - ftell(inpFile) >= BLOCK_SIZE) ? BLOCK_SIZE : fileSize - ftell(inpFile);
        fread(buffer, sizeof(unsigned char), currentBlockSize, inpFile);

        long long blockBits = 0;
#pragma omp parallel reduction(+:blockBits)
        {
            long long localBits = 0;
#pragma omp for
            for (int i = 0; i < currentBlockSize; i++) {
                localBits += codeMap[buffer[i]].length();
            }
            blockBits += localBits;
        }
        totalBits += blockBits;
    }

    delete[] buffer;

    cout << "Hoan thanh gia lap ma hoa\n";
}

// ========================
// CANON
// ========================

struct NodeRef { unsigned long long freq; int index; bool isLeaf; };
struct Leaf { unsigned long long freq; int leader; Leaf(unsigned long long f=0): freq(f), leader(-1) {} };
struct INode { unsigned long long freq; int leader; INode(unsigned long long f=0): freq(f), leader(-1) {} };

static string codeToBits_uint64(unsigned long long code, int len) {
    if (len <= 0) return string();
    if (len > 64) return string(len, '0');
    string s(len, '0');
    for (int i = 0; i < len; ++i) s[len - 1 - i] = ((code & (1ULL << i)) ? '1' : '0');
    return s;
}

void generateCanonCode(long long freq[], map<char, string>& huffmanCode) {
    cout << "--- [PHASE 2] TAO CAU TRUC MA (CANON) ---\n";
    
    const int N = 256;

    // collect symbols with positive frequency
    vector<pair<int, unsigned long long>> symbols;
    symbols.reserve(N);
    for (int i = 0; i < N; ++i) if (freq[i] > 0) symbols.emplace_back(i, (unsigned long long)freq[i]);

    if (symbols.empty()) {
        huffmanCode.clear();
        for (int i = 0; i < N; ++i) huffmanCode[(char)i] = "";
        return;
    }
    if (symbols.size() == 1) {
        huffmanCode.clear();
        for (int i = 0; i < N; ++i) huffmanCode[(char)i] = "";
        huffmanCode[(char)symbols[0].first] = "0";
        return;
    }

    // sort by frequency ascending (stable)
    std::stable_sort(symbols.begin(), symbols.end(),
        [](const pair<int,unsigned long long>& a, const pair<int,unsigned long long>& b){
            if (a.second != b.second) return a.second < b.second;
            return a.first < b.first;
        });

    // CLGeneration (parallel-friendly)
    int m = (int)symbols.size();
    vector<Leaf> lNodes(m);
    for (int i = 0; i < m; ++i) { lNodes[i].freq = symbols[i].second; lNodes[i].leader = -1; }
    vector<INode> iNodes; iNodes.reserve(m);
    vector<int> CL(m, 0);

    int lNodesCur = 0, iNodesFront = 0;

    while (true) {
        unsigned long long midFreq[4];
        bool midIsLeaf[4];
        int midIndex[4];
        int midCount = 0;

        if (lNodesCur < m) {
            midFreq[midCount] = lNodes[lNodesCur].freq; midIsLeaf[midCount] = true; midIndex[midCount] = lNodesCur; midCount++;
            if (lNodesCur + 1 < m) { midFreq[midCount] = lNodes[lNodesCur+1].freq; midIsLeaf[midCount] = true; midIndex[midCount] = lNodesCur+1; midCount++; }
        }
        if (iNodesFront < (int)iNodes.size()) {
            midFreq[midCount] = iNodes[iNodesFront].freq; midIsLeaf[midCount] = false; midIndex[midCount] = iNodesFront; midCount++;
            if (iNodesFront + 1 < (int)iNodes.size()) { midFreq[midCount] = iNodes[iNodesFront+1].freq; midIsLeaf[midCount] = false; midIndex[midCount] = iNodesFront+1; midCount++; }
        }
        if (midCount < 2) break;

        int a=-1,b=-1;
        for (int i=0;i<midCount;++i){
            if (a==-1 || midFreq[i] < midFreq[a] || (midFreq[i]==midFreq[a] && midIsLeaf[i] && !midIsLeaf[a])) { b=a; a=i; }
            else if (b==-1 || midFreq[i] < midFreq[b] || (midFreq[i]==midFreq[b] && midIsLeaf[i] && !midIsLeaf[b])) { b=i; }
        }
        unsigned long long MinFreq = midFreq[a] + midFreq[b];

        int lastLeafIdx = lNodesCur;
        while (lastLeafIdx < m && lNodes[lastLeafIdx].freq <= MinFreq) ++lastLeafIdx;
        int CurLeavesNum = lastLeafIdx - lNodesCur;

        // copy selected leaves deterministically
        vector<NodeRef> Copy(CurLeavesNum);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < CurLeavesNum; ++i) {
            int li = lNodesCur + i;
            Copy[i].freq = lNodes[li].freq;
            Copy[i].isLeaf = true;
            Copy[i].index = li;
        }

        int availableINodes = (int)iNodes.size() - iNodesFront;
        int totalParticipants = CurLeavesNum + availableINodes;
        int mergeFront = iNodesFront;
        int mergeRear = (int)iNodes.size() - 1;
        if (totalParticipants % 2 == 1) {
            bool dropLeaf = false;
            if (CurLeavesNum > 0 && availableINodes > 0) {
                unsigned long long lastLeafFreq = lNodes[lastLeafIdx-1].freq;
                unsigned long long lastINodeFreq = iNodes.back().freq;
                if (lastLeafFreq > lastINodeFreq) dropLeaf = true;
                else if (lastLeafFreq < lastINodeFreq) dropLeaf = false;
                else dropLeaf = false;
            } else if (CurLeavesNum > 0) dropLeaf = true;
            else dropLeaf = false;
            if (dropLeaf) { if (!Copy.empty()) { Copy.pop_back(); CurLeavesNum--; lastLeafIdx--; } }
            else { if (availableINodes > 0) mergeRear = (int)iNodes.size() - 2; }
        }

        // merge Copy and iNodes[mergeFront..mergeRear] (both sorted)
        vector<NodeRef> Temp; Temp.reserve(Copy.size() + max(0, mergeRear - mergeFront + 1));
        int p1 = 0, p2 = mergeFront;
        while (p1 < (int)Copy.size() || p2 <= mergeRear) {
            if (p1 < (int)Copy.size() && (p2 > mergeRear || Copy[p1].freq <= iNodes[p2].freq)) Temp.push_back(Copy[p1++]);
            else { NodeRef r; r.freq = iNodes[p2].freq; r.isLeaf = false; r.index = p2; Temp.push_back(r); ++p2; }
        }

        int TempLength = (int)Temp.size();
        int newINodes = TempLength / 2;
        int baseIndex = (int)iNodes.size();
        iNodes.resize(baseIndex + newINodes);

        // Meld: create new internal nodes and update immediate leaders and CL for leaves
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < newINodes; ++i) {
            NodeRef left = Temp[2*i];
            NodeRef right = Temp[2*i + 1];
            unsigned long long fsum = left.freq + right.freq;
            int ind = baseIndex + i;
            iNodes[ind].freq = fsum;
            iNodes[ind].leader = -1;
            if (left.isLeaf) { lNodes[left.index].leader = ind; 
                #pragma omp atomic 
                CL[left.index]++; }
            else iNodes[left.index].leader = ind;
            if (right.isLeaf) { lNodes[right.index].leader = ind; 
                #pragma omp atomic 
                CL[right.index]++; }
            else iNodes[right.index].leader = ind;
        }

        // Propagate leaders: ensure every leaf climbs up to the highest leader created this round
        while (true) {
            int anyFlag = 0;
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < m; ++i) {
                int curLeader = lNodes[i].leader;
                if (curLeader != -1) {
                    int parentOfLeader = iNodes[curLeader].leader; // may be -1
                    if (parentOfLeader != -1) {
                        lNodes[i].leader = parentOfLeader;
                        #pragma omp atomic
                        CL[i]++;
                        #pragma omp critical
                        { anyFlag = 1; }
                    }
                }
            }
            if (!anyFlag) break;
        }

        // advance iterators
        lNodesCur = lastLeafIdx;
        iNodesFront = mergeRear + 1;
        if (iNodesFront >= (int)iNodes.size() - 1 && lNodesCur >= m) break;
    }

    // CWGeneration: canonical codes (lengths ascending)
    vector<pair<int,int>> lenIdx(m);
    for (int i = 0; i < m; ++i) lenIdx[i] = { CL[i], symbols[i].first };
    std::stable_sort(lenIdx.begin(), lenIdx.end(),
        [](const pair<int,int>& a, const pair<int,int>& b){
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });

    vector<int> lengths(m);
    for (int i = 0; i < m; ++i) lengths[i] = lenIdx[i].first;

    vector<unsigned long long> codesNum(m, 0ULL);
    codesNum[0] = 0ULL;
    for (int i = 1; i < m; ++i) {
        int diff = lengths[i] - lengths[i-1];
        unsigned long long prev = codesNum[i-1];
        unsigned long long nextBase = prev + 1ULL;
        if (diff > 0) {
            if (diff >= 64) codesNum[i] = 0ULL;
            else codesNum[i] = (nextBase << diff);
        } else {
            codesNum[i] = nextBase;
        }
    }

    vector<string> codesByChar(256, "");
    for (int i = 0; i < m; ++i) {
        int origChar = lenIdx[i].second;
        int len = lengths[i];
        if (len > 0) codesByChar[origChar] = codeToBits_uint64(codesNum[i], len);
        else codesByChar[origChar].clear();
    }

    // fill map (assign empty string for non-occurring chars)
    huffmanCode.clear();
    for (int c = 0; c < N; ++c) huffmanCode[(char)c] = codesByChar[c];

    cout << "Hoan thanh tao ma huffman\n";
}