// g++ -std=c++20 c_ans.cpp
// g++ コマンドは、デフォルトのままだと古いC++標準（C++14やC++17など）でコンパイルしようとするため、ranges という名前空間を認識できずにエラーを出してしまう
// →よって、rangesを使うには上記のように新しいversionをコンパイル時に明示する必要がある
#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    using namespace std;

    unsigned N;
    cin >> N;

    vector<pair<unsigned, unsigned>> takahashi;
    for (unsigned i = 0; i < N; ++i) {
        unsigned H, L;
        cin >> H >> L;
        while (!empty(takahashi) && takahashi.back().first <= H) // より早く退室するより背が低い高橋くんがいたら
            takahashi.pop_back(); // はじめに退室してもらう
        takahashi.emplace_back(H, L);
    }

    unsigned Q;
    cin >> Q;

    for (unsigned i = 0; i < Q; ++i) {
        unsigned T;
        cin >> T;
        // pair の第二要素で二分探索し、T 分より後に退出する中で最も背が高い高橋くんを見つける
        cout << ranges::upper_bound(takahashi, T, {}, &pair<unsigned, unsigned>::second)->first << endl;
    }
    return 0;
}
