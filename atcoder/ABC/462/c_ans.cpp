#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    // コンテスト向けの入出力高速化
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    // Pythonのリストのタプル xy = [(x, y), ...] に相当する pair のベクター
    vector<pair<int, int>> xy(n);
    for (int i = 0; i < n; ++i) {
        cin >> xy[i].first >> xy[i].second;
    }

    // X座標で昇順にソート (C++のpairは自動的に1つ目の要素(first)を基準にソートされます)
    sort(xy.begin(), xy.end());

    int ans = 0;
    
    // Pythonの float('inf') の代わりに、問題の制約(N <= 300,000)よりも十分に大きな値を設定します
    // 1e9 は 10の9乗 (1,000,000,000) を意味します
    int min_y = 1e9; 

    for (int i = 0; i < n; ++i) {
        // Pythonの x_i, y_i = xy[i][0], xy[i][1] に相当
        int x_i = xy[i].first;
        int y_i = xy[i].second;
        
        if (y_i < min_y) {
            ans++;
            min_y = y_i;
        }
    }

    cout << ans << "\n";

    return 0;
}